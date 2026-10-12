"""Full spatial Hardy electron flux model, Hardy et al. (1987).

Tables 4/5 and equations 4/5 (printed pp. 12284--12285) fit log10 integral
number/energy flux. Mean energy is their ratio (Hardy et al., 1985, p. 4231).
The moments include the full electron spectrum, measured 50 eV--20 keV and
extrapolated to 100 keV. They differ from the >=661 eV moments used for the
1987 paper's fitted conductivities. Feeding this mean energy to Robinson
therefore does not reproduce Hardy's conductivity fits.

Printed equation 5 omits pi from the harmonic phase although T is MLT hours.
We infer the normalization n*pi*T/12 required for a 24-hour Fourier series.
Fractional-Kp interpolation is an implementation choice, not a prescription
found in this paper: completed log10 fluxes are linearly interpolated after
applying the published backgrounds. The highest group is Kp>=6-.
"""

#%% Imports and provenance

import numpy as np

from .hardy_electron_coefficients import (
    ENERGY_FLUX_COEFFICIENTS, NUMBER_FLUX_COEFFICIENTS,
    HARDY_ELECTRON_MODEL_DOI, HARDY_ELECTRON_COEFFICIENT_VERSION
)

HARDY_ELECTRON_KP_POLICY = "clip_0_6_interpolate_evaluated_log10_flux"
HARDY_ELECTRON_LATITUDE_RANGE = (50.0, 90.0)
HARDY_ELECTRON_UNCERTAINTY_METHOD = "measurement_conditional_on_fixed_hardy_energy"
EQUATORWARD_FLUX_BACKGROUND = 1e6
POLEWARD_FLUX_BACKGROUND = 1e7


#%% Published latitude profile and Fourier parameters

def _fourier_parameters(coefficients, kp_level, mlt):
    """Return r, h0, S1, S2 from constant,cos1..6,sin1..6 rows."""

    # Select one Fourier row at a time: an orbit-sized spatial grid should
    # not materialize all 52 coefficients at every pixel simultaneously.
    parameters = coefficients[kp_level, 0, :].copy()
    for harmonic in range(1, 7):
        phase = np.pi * harmonic * mlt / 12
        parameters += coefficients[kp_level, harmonic, :] * np.cos(phase)[..., None]
        parameters += coefficients[kp_level, harmonic + 6, :] * np.sin(phase)[..., None]
    return parameters


def _epstein_log_flux(parameters, mlat):
    """Equation 4 with a stable natural-log transition and log10-flux output."""

    r, h0, S1, S2 = np.moveaxis(parameters, -1, 0)
    valid_slopes = (S1 > 0) & (S2 < 0)
    # These signs hold across the published integer-Kp/MLT model. Invalid
    # Fourier extrapolations are unavailable rather than complex logarithms.
    q = np.where(valid_slopes, -S1 / np.where(S2 != 0, S2, -1), 1)
    latitude_offset = mlat - h0
    transition = np.logaddexp(0, np.log(q) + latitude_offset) - np.log1p(q)
    log_flux = r + S1 * latitude_offset + (S2 - S1) * transition
    return np.where(valid_slopes, log_flux, np.nan)


def _log_flux_at_kp_level(coefficients, kp_level, mlt, mlat):
    """Reconstruct parameters, latitude profile, then published backgrounds."""

    parameters = _fourier_parameters(coefficients, kp_level, mlt)
    log_flux = _epstein_log_flux(parameters, mlat)
    h0 = parameters[..., 1]
    log_flux = np.where((mlat < h0) & (log_flux < 6), 6, log_flux)
    log_flux = np.where((mlat > h0) & (log_flux < 7), 7, log_flux)
    return log_flux


#%% Public spatial model

def hardy_electron_precipitation(kp, mlt, mlat):
    """Return Hardy electron number flux, energy flux and mean energy.

    Inputs broadcast normally: Kp (numerical), MLT (hours, wrapped 0--24),
    northern corrected geomagnetic latitude (degrees, 50--90 inclusive).
    Out-of-domain latitude, masked or nonfinite inputs return NaN. There is
    no implicit hemisphere reflection. Modified Apex latitude is an explicit
    caller approximation, not the coordinates used to fit the original data.

    Number flux is electrons cm^-2 s^-1 sr^-1; energy flux is keV cm^-2 s^-1
    sr^-1; their ratio is mean energy in keV. Finite Kp clips 0--6; fractional
    values interpolate completed log10 flux profiles. No model uncertainty
    is assigned. Constants expose the source DOI, version and domain policy.
    """

    kp, mlt, mlat = np.broadcast_arrays(*[
        np.asarray(np.ma.filled(np.ma.asarray(value, dtype=float), np.nan))
        for value in (kp, mlt, mlat)
    ])
    valid = np.isfinite(kp) & np.isfinite(mlt) & np.isfinite(mlat)
    valid = valid & (mlat >= 50) & (mlat <= 90)

    kp_safe = np.clip(np.where(valid, kp, 0), 0, 6)
    lower_kp = np.floor(kp_safe).astype(int)
    upper_kp = np.ceil(kp_safe).astype(int)
    fraction = kp_safe - lower_kp
    mlt_safe = np.mod(np.where(valid, mlt, 0), 24)
    mlat_safe = np.where(valid, mlat, 50)

    number_lower = _log_flux_at_kp_level(NUMBER_FLUX_COEFFICIENTS, lower_kp, mlt_safe, mlat_safe)
    number_upper = _log_flux_at_kp_level(NUMBER_FLUX_COEFFICIENTS, upper_kp, mlt_safe, mlat_safe)
    energy_lower = _log_flux_at_kp_level(ENERGY_FLUX_COEFFICIENTS, lower_kp, mlt_safe, mlat_safe)
    energy_upper = _log_flux_at_kp_level(ENERGY_FLUX_COEFFICIENTS, upper_kp, mlt_safe, mlat_safe)

    number_flux = 10 ** (number_lower + fraction * (number_upper - number_lower))
    energy_flux = 10 ** (energy_lower + fraction * (energy_upper - energy_lower))
    mean_energy = energy_flux / number_flux
    valid = valid & np.isfinite(mean_energy) & (mean_energy > 0)
    return {
        "number_flux": np.where(valid, number_flux, np.nan),
        "energy_flux": np.where(valid, energy_flux, np.nan),
        "mean_energy": np.where(valid, mean_energy, np.nan)
    }
