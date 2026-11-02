"""Fang (2025) electron-impact Hall and Pedersen conductance.

The parameterization implements equations 7--10 and Tables 2--3 of Fang
(2025), DOI 10.1029/2025JA033835.  It requires Maxwellian *average* electron
energy, which is twice the characteristic energy under the paper's
definition.  The magnetic-field magnitude must be evaluated at 150 km by the
caller; this module does not evaluate a field model.
"""

#%% Imports

import numpy as np


#%% Published model definition

FANG_2025_MODEL_DOI = "10.1029/2025JA033835"
FANG_2025_ENERGY_RANGE_KEV = (1.0, 20.0)
FANG_2025_REFERENCE_FIELD_MICROTESLA = 50.0

# Equation 9, Table 2: coefficients are ordered g0 through g3.
FANG_2025_PEDERSEN_GAMMA_COEFFICIENTS = (
    2.582618, 1.126099, 0.005476, 1.276113
)
FANG_2025_HALL_GAMMA_COEFFICIENTS = (
    7.075787, 6.946879, -0.014756, 6.201520
)

# Equation 10, Table 3: coefficients are ordered c0 through c6.
FANG_2025_PEDERSEN_ENERGY_COEFFICIENTS = (
    -3.107071, 11.393469, 0.185083, 0.000885,
    0.944358, 0.261833, -0.000190
)
FANG_2025_HALL_ENERGY_COEFFICIENTS = (
    -0.310866, -0.348437, 3.159323, 0.001002,
    0.165759, 0.225693, 0.000918
)


#%% Numerical helpers

def _as_float64(value):
    """Convert numeric or masked input to float64, preserving missingness."""
    return np.ma.asarray(value, dtype=np.float64).filled(np.nan)


def _broadcast_inputs(*values):
    """Broadcast float64 inputs and give shape errors a model-specific message."""
    arrays = [_as_float64(value) for value in values]
    try:
        return np.broadcast_arrays(*arrays)
    except ValueError as error:
        raise ValueError("all Fang 2025 inputs must broadcast together") from error


def _check_low_energy_policy(low_energy_policy):
    if low_energy_policy not in ("clip", "reject"):
        raise ValueError("low_energy_policy must be 'clip' or 'reject'")


def _gamma_and_derivative(energy, coefficients):
    """Evaluate equation 9 and its derivative with respect to energy."""
    g0, g1, g2, g3 = coefficients
    numerator = g0 + g1 * energy + g2 * energy**2
    denominator = 1.0 + g3 * energy
    gamma = numerator / denominator

    numerator_derivative = g1 + 2.0 * g2 * energy
    gamma_derivative = (
        numerator_derivative * denominator - numerator * g3
    ) / denominator**2
    return gamma, gamma_derivative


def _energy_term_and_log_derivative(energy, coefficients):
    """Evaluate equation 10 and ``d(log(S)) / dE``."""
    c0, c1, c2, c3, c4, c5, c6 = coefficients
    numerator = c0 + c1 * energy + c2 * energy**2 + c3 * energy**3
    denominator = 1.0 + c4 * energy + c5 * energy**2 + c6 * energy**3
    energy_term = numerator / denominator

    numerator_derivative = c1 + 2.0 * c2 * energy + 3.0 * c3 * energy**2
    denominator_derivative = c4 + 2.0 * c5 * energy + 3.0 * c6 * energy**2
    log_derivative = (
        numerator_derivative / numerator
        - denominator_derivative / denominator
    )
    return energy_term, log_derivative


def _evaluate_one_conductance(energy, energy_flux, magnetic_field,
                              gamma_coefficients, energy_coefficients):
    """Evaluate one conductance and its smooth-domain Jacobian."""
    gamma, gamma_derivative = _gamma_and_derivative(
        energy, gamma_coefficients
    )
    energy_term, energy_log_derivative = _energy_term_and_log_derivative(
        energy, energy_coefficients
    )

    field_ratio = magnetic_field / FANG_2025_REFERENCE_FIELD_MICROTESLA
    conductance = np.sqrt(energy_flux) * energy_term * field_ratio**(-gamma)

    derivative_energy = conductance * (
        energy_log_derivative - gamma_derivative * np.log(field_ratio)
    )
    derivative_field = -gamma * conductance / magnetic_field
    derivative_flux = np.divide(
        conductance, 2.0 * energy_flux,
        out=np.zeros_like(conductance), where=energy_flux > 0
    )
    return conductance, derivative_energy, derivative_flux, derivative_field


def _nominal_from_broadcast(energy, energy_flux, magnetic_field,
                            low_energy_policy):
    """Evaluate already-broadcast inputs and retain internal Jacobians."""
    minimum_energy, maximum_energy = FANG_2025_ENERGY_RANGE_KEV
    finite_energy = np.isfinite(energy)
    within_model_domain = (
        finite_energy & (energy >= minimum_energy) & (energy <= maximum_energy)
    )
    energy_clipped_low = finite_energy & (energy > 0) & (energy < minimum_energy)
    apply_low_energy_floor = energy_clipped_low & (low_energy_policy == "clip")
    energy_usable = within_model_domain | apply_low_energy_floor

    evaluated_energy = np.where(
        within_model_domain, energy,
        np.where(apply_low_energy_floor, minimum_energy, np.nan)
    )
    input_valid = (
        energy_usable
        & np.isfinite(energy_flux) & (energy_flux >= 0)
        & np.isfinite(magnetic_field) & (magnetic_field > 0)
    )

    # Safe placeholders prevent invalid cells from producing warnings while
    # retaining NaN in every public result at those cells.
    safe_energy = np.where(input_valid, evaluated_energy, minimum_energy)
    safe_flux = np.where(input_valid, energy_flux, 0.0)
    safe_field = np.where(
        input_valid, magnetic_field, FANG_2025_REFERENCE_FIELD_MICROTESLA
    )

    pedersen = _evaluate_one_conductance(
        safe_energy, safe_flux, safe_field,
        FANG_2025_PEDERSEN_GAMMA_COEFFICIENTS,
        FANG_2025_PEDERSEN_ENERGY_COEFFICIENTS
    )
    hall = _evaluate_one_conductance(
        safe_energy, safe_flux, safe_field,
        FANG_2025_HALL_GAMMA_COEFFICIENTS,
        FANG_2025_HALL_ENERGY_COEFFICIENTS
    )

    valid = input_valid & np.isfinite(pedersen[0]) & np.isfinite(hall[0])
    result = {
        "P": np.where(valid, pedersen[0], np.nan),
        "H": np.where(valid, hall[0], np.nan),
        "valid": valid,
        "within_model_domain": within_model_domain,
        "energy_clipped_low": energy_clipped_low,
        "input_energy": np.array(energy, dtype=np.float64, copy=True),
        "evaluated_energy": np.where(energy_usable, evaluated_energy, np.nan)
    }
    return result, pedersen[1:], hall[1:]


#%% Public nominal calculation

def fang_2025_conductance(mean_energy, energy_flux, magnetic_field, *,
                          low_energy_policy="clip"):
    """Return Fang (2025) electron-impact conductance in siemens.

    Parameters
    ----------
    mean_energy : float or array
        Maxwellian average electron energy [keV].  This is twice the
        characteristic energy under Fang's convention.
    energy_flux : float or array
        Incident electron energy flux [mW m-2], numerically equal to
        ``erg cm-2 s-1``.
    magnetic_field : float or array
        Positive magnetic-field magnitude at 150 km [microtesla].
    low_energy_policy : {"clip", "reject"}, optional
        ``"clip"`` evaluates positive energies below 1 keV at 1 keV and
        reports positive sub-1-keV inputs in ``energy_clipped_low``.
        ``"reject"`` retains the same input flag but returns unavailable
        conductance for those cells. Energies above 20 keV are always rejected.

    Returns
    -------
    dict
        ``P`` and ``H`` are Pedersen and Hall conductance [S]. ``valid`` marks
        finite operational estimates. ``within_model_domain`` describes the
        original 1--20 keV input, while ``input_energy`` and
        ``evaluated_energy`` keep the lower-floor operation explicit.

    Notes
    -----
    The result is the electron-precipitation contribution only. It does not
    include proton or solar-EUV conductance. The paper calibrates Maxwellian
    precipitation over 1--20 keV average energy; the lower floor is an
    operational extrapolation selected for IMAGE coverage.
    """
    _check_low_energy_policy(low_energy_policy)
    energy, flux, field = _broadcast_inputs(
        mean_energy, energy_flux, magnetic_field
    )
    result, _, _ = _nominal_from_broadcast(
        energy, flux, field, low_energy_policy
    )
    return result


#%% First-order uncertainty

def fang_2025_conductance_uncertainty(
        mean_energy, energy_flux, magnetic_field, *,
        mean_energy_std, energy_flux_std, magnetic_field_std,
        energy_flux_covariance=0.0, low_energy_policy="clip"):
    """Return Fang conductance with first-order input uncertainty.

    Standard deviations use the same units as their corresponding nominal
    inputs. ``energy_flux_covariance`` is covariance between mean energy and
    energy flux [keV mW m-2]. Other input covariances are assumed zero.
    ``P_std`` and ``H_std`` are in siemens and ``PH_covariance`` is in S2.

    Propagation uses analytical Jacobians. Below the operational 1-keV floor,
    the energy derivative is zero. ``uncertainty_straddles_low_energy_boundary``
    identifies one-standard-deviation intervals crossing that nonlinear
    boundary; sampled propagation should be preferred there. Invalid errors
    leave nominal conductance unchanged but make uncertainty unavailable.

    At exactly zero nominal flux, analytical propagation is valid only when
    flux uncertainty and energy--flux covariance are both zero. This is the
    differentiable fixed-zero-flux case; otherwise the square-root flux
    derivative is singular.
    """
    _check_low_energy_policy(low_energy_policy)
    arrays = _broadcast_inputs(
        mean_energy, energy_flux, magnetic_field,
        mean_energy_std, energy_flux_std, magnetic_field_std,
        energy_flux_covariance
    )
    energy, flux, field, energy_std, flux_std, field_std, covariance = arrays
    result, pedersen_jacobian, hall_jacobian = _nominal_from_broadcast(
        energy, flux, field, low_energy_policy
    )

    error_inputs_valid = (
        np.isfinite(energy_std) & (energy_std >= 0)
        & np.isfinite(flux_std) & (flux_std >= 0)
        & np.isfinite(field_std) & (field_std >= 0)
        & np.isfinite(covariance)
    )
    covariance_limit = energy_std * flux_std
    covariance_valid = np.abs(covariance) <= covariance_limit * (1.0 + 1e-14)
    zero_flux_smooth = (flux > 0) | ((flux_std == 0) & (covariance == 0))
    uncertainty_valid = (
        result["valid"] & error_inputs_valid & covariance_valid
        & zero_flux_smooth
    )

    p_energy, p_flux, p_field = pedersen_jacobian
    h_energy, h_flux, h_field = hall_jacobian

    # The operational floor is constant with respect to energy below 1 keV.
    p_energy = np.where(result["energy_clipped_low"], 0.0, p_energy)
    h_energy = np.where(result["energy_clipped_low"], 0.0, h_energy)

    p_variance = (
        (p_energy * energy_std)**2
        + (p_flux * flux_std)**2
        + (p_field * field_std)**2
        + 2.0 * p_energy * p_flux * covariance
    )
    h_variance = (
        (h_energy * energy_std)**2
        + (h_flux * flux_std)**2
        + (h_field * field_std)**2
        + 2.0 * h_energy * h_flux * covariance
    )
    ph_covariance = (
        p_energy * h_energy * energy_std**2
        + p_flux * h_flux * flux_std**2
        + p_field * h_field * field_std**2
        + (p_energy * h_flux + p_flux * h_energy) * covariance
    )

    finite_propagation = (
        np.isfinite(p_variance) & np.isfinite(h_variance)
        & np.isfinite(ph_covariance)
    )
    uncertainty_valid &= finite_propagation
    p_variance = np.maximum(p_variance, 0.0)
    h_variance = np.maximum(h_variance, 0.0)

    lower = energy - energy_std
    upper = energy + energy_std
    straddles_low_boundary = (
        result["valid"] & error_inputs_valid & (energy_std > 0)
        & (lower < FANG_2025_ENERGY_RANGE_KEV[0])
        & (upper > FANG_2025_ENERGY_RANGE_KEV[0])
    )

    result.update({
        "P_std": np.where(uncertainty_valid, np.sqrt(p_variance), np.nan),
        "H_std": np.where(uncertainty_valid, np.sqrt(h_variance), np.nan),
        "PH_covariance": np.where(
            uncertainty_valid, ph_covariance, np.nan
        ),
        "uncertainty_valid": uncertainty_valid,
        "uncertainty_straddles_low_energy_boundary": straddles_low_boundary
    })
    return result
