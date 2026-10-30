"""Hardy et al. (1987) empirical auroral conductance model.

Tables 6 and 7 and equations 4--5 describe independent Epstein--Fourier
fits to height-integrated Hall and Pedersen conductance.  These fitted
conductances are not obtained by applying Robinson to Hardy's published
full-spectrum electron moments.
"""

#%% Imports and provenance

from functools import lru_cache
from importlib.resources import files

import numpy as np


HARDY_CONDUCTANCE_MODEL_DOI = "10.1029/JA092iA11p12275"
HARDY_CONDUCTANCE_COEFFICIENT_VERSION = "hardy1987-tables6-7-v1"
HARDY_CONDUCTANCE_KP_POLICY = "clip_0_6_interpolate_evaluated_conductance"
HARDY_CONDUCTANCE_LATITUDE_RANGE = (50.0, 90.0)
EQUATORWARD_CONDUCTANCE_BACKGROUND = 0.0
POLEWARD_CONDUCTANCE_BACKGROUND = 0.55
HARDY_EFFECTIVE_MINIMUM_ENERGY_KEV = 0.5
HARDY_EFFECTIVE_METHOD = "inverse_robinson_from_hardy_conductance"


#%% Published coefficient tables

@lru_cache(maxsize=2)
def _load_coefficients(name):
    """Load one immutable (Kp, Fourier term, Epstein parameter) table."""

    resource = files("icphysics").joinpath("data", f"hardy_{name}_coefficients.csv")
    with resource.open("r", encoding="utf-8") as stream:
        table_lines = [line for line in stream if line.count(",") == 5]
    coefficients = np.loadtxt(table_lines, delimiter=",", skiprows=1, usecols=(2, 3, 4, 5))
    coefficients = coefficients.reshape(7, 13, 4)
    coefficients.setflags(write=False)
    return coefficients


def _fourier_parameters(coefficients, kp_level, mlt):
    """Return peak, peak latitude, equatorward slope and poleward slope."""

    parameters = coefficients[kp_level, 0, :].copy()
    for harmonic in range(1, 7):
        phase = np.pi * harmonic * mlt / 12
        parameters += coefficients[kp_level, harmonic, :] * np.cos(phase)[..., None]
        parameters += coefficients[kp_level, harmonic + 6, :] * np.sin(phase)[..., None]
    return parameters


def _epstein_conductance_with_background(parameters, mlat):
    """Evaluate equation 4 and identify use of published backgrounds."""

    peak, peak_latitude, equatorward_slope, poleward_slope = np.moveaxis(parameters, -1, 0)
    valid_slopes = (equatorward_slope > 0) & (poleward_slope < 0)
    slope_ratio = np.where(
        valid_slopes,
        -equatorward_slope / np.where(poleward_slope != 0, poleward_slope, -1),
        1
    )
    latitude_offset = mlat - peak_latitude
    transition = np.logaddexp(0, np.log(slope_ratio) + latitude_offset) - np.log1p(slope_ratio)
    conductance = peak + equatorward_slope * latitude_offset
    conductance += (poleward_slope - equatorward_slope) * transition
    conductance = np.where(valid_slopes, conductance, np.nan)

    equatorward_background = (
        (mlat < peak_latitude)
        & (conductance < EQUATORWARD_CONDUCTANCE_BACKGROUND)
    )
    poleward_background = (
        (mlat > peak_latitude)
        & (conductance < POLEWARD_CONDUCTANCE_BACKGROUND)
    )
    conductance = np.where(
        equatorward_background,
        EQUATORWARD_CONDUCTANCE_BACKGROUND,
        conductance
    )
    conductance = np.where(
        poleward_background,
        POLEWARD_CONDUCTANCE_BACKGROUND,
        conductance
    )
    return conductance, equatorward_background | poleward_background


def _epstein_conductance(parameters, mlat):
    """Evaluate equation 4 stably and apply the published backgrounds."""

    conductance, _ = _epstein_conductance_with_background(parameters, mlat)
    return conductance


def _conductance_at_kp_level(coefficients, kp_level, mlt, mlat):
    parameters = _fourier_parameters(coefficients, kp_level, mlt)
    return _epstein_conductance(parameters, mlat)


def _conductance_and_background_at_kp_level(
    coefficients, kp_level, mlt, mlat
):
    parameters = _fourier_parameters(coefficients, kp_level, mlt)
    return _epstein_conductance_with_background(parameters, mlat)


def _evaluate_hardy_conductance(kp, mlt, mlat):
    """Evaluate conductance and track whether either Kp endpoint used a floor."""

    kp, mlt, mlat = np.broadcast_arrays(*[
        np.asarray(np.ma.filled(np.ma.asarray(value, dtype=float), np.nan))
        for value in (kp, mlt, mlat)
    ])
    valid = np.isfinite(kp) & np.isfinite(mlt) & np.isfinite(mlat)
    valid &= (mlat >= HARDY_CONDUCTANCE_LATITUDE_RANGE[0])
    valid &= (mlat <= HARDY_CONDUCTANCE_LATITUDE_RANGE[1])

    kp_safe = np.clip(np.where(valid, kp, 0), 0, 6)
    lower_kp = np.floor(kp_safe).astype(int)
    upper_kp = np.ceil(kp_safe).astype(int)
    fraction = kp_safe - lower_kp
    mlt_safe = np.mod(np.where(valid, mlt, 0), 24)
    mlat_safe = np.where(valid, mlat, HARDY_CONDUCTANCE_LATITUDE_RANGE[0])

    result = {}
    background = {}
    for name in ("hall", "pedersen"):
        coefficients = _load_coefficients(name)
        lower, lower_background = _conductance_and_background_at_kp_level(
            coefficients, lower_kp, mlt_safe, mlat_safe
        )
        upper, upper_background = _conductance_and_background_at_kp_level(
            coefficients, upper_kp, mlt_safe, mlat_safe
        )
        conductance = lower + fraction * (upper - lower)
        result[name] = np.where(
            valid & np.isfinite(conductance), conductance, np.nan
        )

        # Fractional Kp is valid only when every contributing endpoint comes
        # from the fitted auroral profile rather than a published background.
        endpoint_background = lower_background
        endpoint_background |= (fraction > 0) & upper_background
        background[name] = endpoint_background | ~valid
    return result, background


#%% Public spatial model

def hardy_conductance(kp, mlt, mlat):
    """Return Hardy electron-produced Hall and Pedersen conductance in S.

    Inputs broadcast normally: numerical Kp, MLT in hours, and northern
    corrected geomagnetic latitude in degrees. MLT wraps at 24 hours; finite
    Kp clips to 0--6 and fractional Kp linearly interpolates the evaluated
    integer-Kp conductances. Latitude must be within 50--90 degrees. Masked,
    nonfinite, out-of-domain, or invalid fitted profiles return NaN.

    The result is a dictionary with ``hall`` and ``pedersen`` arrays. These
    are the auroral electron contributions only; solar-EUV conductance is not
    included and no model uncertainty is assigned.
    """

    result, _ = _evaluate_hardy_conductance(kp, mlt, mlat)
    return result


def _inverse_robinson_conductance(hall_conductance, pedersen_conductance):
    """Invert positive finite Hall/Pedersen pairs through Robinson."""

    hall_conductance, pedersen_conductance = np.broadcast_arrays(
        np.asarray(hall_conductance, dtype=float),
        np.asarray(pedersen_conductance, dtype=float),
    )
    invertible = (
        np.isfinite(hall_conductance)
        & np.isfinite(pedersen_conductance)
        & (hall_conductance > 0)
        & (pedersen_conductance > 0)
    )

    effective_energy = np.full(hall_conductance.shape, np.nan)
    effective_energy[invertible] = (
        hall_conductance[invertible]
        / (0.45 * pedersen_conductance[invertible])
    ) ** (1 / 0.85)

    effective_energy_flux = np.full(hall_conductance.shape, np.nan)
    effective_energy_flux[invertible] = (
        pedersen_conductance[invertible]
        * (16 + effective_energy[invertible]**2)
        / (40 * effective_energy[invertible])
    ) ** 2
    return effective_energy, effective_energy_flux, invertible


def hardy_robinson_effective_precipitation(kp, mlt, mlat):
    """Return Robinson-effective precipitation for native Hardy conductance.

    Inputs follow :func:`hardy_conductance`: numerical Kp, MLT in hours and
    northern corrected geomagnetic latitude in degrees broadcast normally.
    ``effective_energy`` is in keV and ``effective_energy_flux`` is in
    mW m-2. Substitution into the Robinson equations reproduces the native
    Hardy Hall and Pedersen conductances wherever ``valid`` is true.

    The inverse values are conductance-equivalent parameters, not Hardy's
    published full-spectrum mean energy/flux or the unpublished moments after
    its approximately 500-eV cutoff. Separate smoothing and functional fits
    of Hall and Pedersen conductance can make those quantities differ.

    Raw inverse values remain available wherever both conductances are finite
    and positive. ``valid`` additionally requires energy >= 0.5 keV and rejects
    points where either conductance used Hardy's equatorward or poleward
    background. No upper-energy or energy-flux threshold is imposed.
    """

    conductance, background = _evaluate_hardy_conductance(kp, mlt, mlat)
    hall_conductance = conductance["hall"]
    pedersen_conductance = conductance["pedersen"]
    effective_energy, effective_energy_flux, invertible = (
        _inverse_robinson_conductance(
            hall_conductance, pedersen_conductance
        )
    )

    valid = invertible
    valid &= ~background["hall"] & ~background["pedersen"]
    valid &= effective_energy >= HARDY_EFFECTIVE_MINIMUM_ENERGY_KEV
    return {
        "effective_energy": effective_energy,
        "effective_energy_flux": effective_energy_flux,
        "valid": valid,
        "hall": hall_conductance,
        "pedersen": pedersen_conductance,
        "hall_background": background["hall"],
        "pedersen_background": background["pedersen"],
    }
