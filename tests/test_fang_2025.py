"""Independent checks of the Fang (2025) conductance parameterization."""

import numpy as np

from icphysics import (
    fang_2025_conductance,
    fang_2025_conductance_uncertainty,
)
from icphysics.fang_2025 import (
    FANG_2025_HALL_ENERGY_COEFFICIENTS,
    FANG_2025_HALL_GAMMA_COEFFICIENTS,
    FANG_2025_PEDERSEN_ENERGY_COEFFICIENTS,
    FANG_2025_PEDERSEN_GAMMA_COEFFICIENTS,
)


def literal_fang(mean_energy, energy_flux, magnetic_field):
    """Literal equations 7--10, separate from production helpers."""
    energy = np.asarray(mean_energy, dtype=np.float64)
    flux = np.asarray(energy_flux, dtype=np.float64)
    field = np.asarray(magnetic_field, dtype=np.float64)

    coefficients = {
        "P": {
            "g": (2.582618, 1.126099, 0.005476, 1.276113),
            "c": (-3.107071, 11.393469, 0.185083, 0.000885,
                  0.944358, 0.261833, -0.000190)
        },
        "H": {
            "g": (7.075787, 6.946879, -0.014756, 6.201520),
            "c": (-0.310866, -0.348437, 3.159323, 0.001002,
                  0.165759, 0.225693, 0.000918)
        }
    }

    result = {}
    for name, values in coefficients.items():
        g0, g1, g2, g3 = values["g"]
        c0, c1, c2, c3, c4, c5, c6 = values["c"]
        gamma = (g0 + g1 * energy + g2 * energy**2) / (1 + g3 * energy)
        energy_term = (
            (c0 + c1 * energy + c2 * energy**2 + c3 * energy**3)
            / (1 + c4 * energy + c5 * energy**2 + c6 * energy**3)
        )
        result[name] = (
            np.sqrt(flux) * energy_term * (field / 50.0)**(-gamma)
        )
    return result


def central_derivative(function, value, step):
    return (function(value + step) - function(value - step)) / (2 * step)


def test_published_table_coefficients_are_transcribed_exactly():
    assert FANG_2025_PEDERSEN_GAMMA_COEFFICIENTS == (
        2.582618, 1.126099, 0.005476, 1.276113
    )
    assert FANG_2025_HALL_GAMMA_COEFFICIENTS == (
        7.075787, 6.946879, -0.014756, 6.201520
    )
    assert FANG_2025_PEDERSEN_ENERGY_COEFFICIENTS == (
        -3.107071, 11.393469, 0.185083, 0.000885,
        0.944358, 0.261833, -0.000190
    )
    assert FANG_2025_HALL_ENERGY_COEFFICIENTS == (
        -0.310866, -0.348437, 3.159323, 0.001002,
        0.165759, 0.225693, 0.000918
    )


def test_scalar_checkpoints_from_literal_published_equations():
    checkpoints = (
        (1.0, 1.0, 50.0, 3.840599346963125, 1.7962337597046765),
        (5.0, 4.0, 40.0, 12.320642292438832, 27.120614055460827),
        (20.0, 0.25, 65.0, 0.948401523946965, 4.622332276585024)
    )
    for energy, flux, field, expected_p, expected_h in checkpoints:
        actual = fang_2025_conductance(energy, flux, field)
        np.testing.assert_allclose(actual["P"], expected_p, rtol=2e-15)
        np.testing.assert_allclose(actual["H"], expected_h, rtol=2e-15)
        assert actual["valid"]


def test_direct_equation_parity_over_energy_flux_and_field_grid():
    energy = np.array([1.0, 2.5, 8.0, 20.0])[:, None, None]
    flux = np.array([0.0, 0.2, 3.0])[None, :, None]
    field = np.array([32.0, 50.0, 68.0])[None, None, :]
    expected = literal_fang(energy, flux, field)
    actual = fang_2025_conductance(energy, flux, field)

    assert actual["P"].shape == (4, 3, 3)
    np.testing.assert_allclose(actual["P"], expected["P"], rtol=2e-15)
    np.testing.assert_allclose(actual["H"], expected["H"], rtol=2e-15)
    assert np.all(actual["valid"])


def test_broadcasting_preserves_inputs_and_float64_outputs():
    energy = np.array([[1.0], [5.0]])
    flux = np.array([0.5, 1.0, 2.0])
    field = np.float32(47.0)
    original_energy = energy.copy()
    original_flux = flux.copy()

    actual = fang_2025_conductance(energy, flux, field)

    assert actual["P"].shape == (2, 3)
    assert actual["P"].dtype == np.float64
    assert actual["input_energy"].dtype == np.float64
    np.testing.assert_array_equal(energy, original_energy)
    np.testing.assert_array_equal(flux, original_flux)
    np.testing.assert_array_equal(actual["input_energy"], np.broadcast_to(energy, (2, 3)))


def test_masked_nonfinite_and_nonphysical_inputs_are_unavailable():
    energy = np.ma.array(
        [5.0, 5.0, np.nan, 0.0, -1.0, 21.0, 5.0, 5.0],
        mask=[False, True, False, False, False, False, False, False]
    )
    flux = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -0.1, 1.0])
    field = np.array([50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 0.0])

    actual = fang_2025_conductance(energy, flux, field)

    np.testing.assert_array_equal(
        actual["valid"], [True, False, False, False, False, False, False, False]
    )
    assert np.all(np.isnan(actual["P"][1:]))
    assert np.all(np.isnan(actual["H"][1:]))
    assert np.isnan(actual["input_energy"][1])


def test_zero_flux_is_exactly_zero_when_energy_and_field_are_valid():
    actual = fang_2025_conductance([0.5, 1.0, 20.0], 0.0, 50.0)
    np.testing.assert_array_equal(actual["P"], np.zeros(3))
    np.testing.assert_array_equal(actual["H"], np.zeros(3))
    assert np.all(actual["valid"])


def test_energy_boundaries_floor_and_domain_flags_are_explicit():
    energy = np.array([0.0, 0.4, 0.999, 1.0, 20.0, 20.001])
    actual = fang_2025_conductance(energy, 1.0, 50.0)

    np.testing.assert_array_equal(
        actual["within_model_domain"], [False, False, False, True, True, False]
    )
    np.testing.assert_array_equal(
        actual["energy_clipped_low"], [False, True, True, False, False, False]
    )
    np.testing.assert_array_equal(actual["input_energy"], energy)
    np.testing.assert_allclose(actual["evaluated_energy"][1:4], [1.0, 1.0, 1.0])
    assert np.isnan(actual["evaluated_energy"][0])
    assert np.isnan(actual["evaluated_energy"][-1])
    assert np.isnan(actual["P"][0]) and np.isnan(actual["P"][-1])
    np.testing.assert_array_equal(actual["P"][1:4], np.repeat(actual["P"][3], 3))

    rejected = fang_2025_conductance(0.5, 1.0, 50.0, low_energy_policy="reject")
    assert not rejected["valid"]
    assert rejected["energy_clipped_low"]
    assert np.isnan(rejected["P"])

    with np.testing.assert_raises_regex(ValueError, "low_energy_policy"):
        fang_2025_conductance(5.0, 1.0, 50.0, low_energy_policy="extrapolate")


def test_magnetic_scaling_uses_separate_energy_dependent_exponents():
    energy = 4.0
    lower = fang_2025_conductance(energy, 1.0, 35.0)
    upper = fang_2025_conductance(energy, 1.0, 65.0)
    reference = literal_fang(energy, 1.0, 50.0)

    gp = (2.582618 + 1.126099 * energy + 0.005476 * energy**2) / (
        1 + 1.276113 * energy
    )
    gh = (7.075787 + 6.946879 * energy - 0.014756 * energy**2) / (
        1 + 6.201520 * energy
    )
    np.testing.assert_allclose(upper["P"] / lower["P"], (65 / 35)**(-gp))
    np.testing.assert_allclose(upper["H"] / lower["H"], (65 / 35)**(-gh))
    assert not np.isclose(gp, gh)
    np.testing.assert_allclose(
        fang_2025_conductance(energy, 1.0, 50.0)["P"], reference["P"]
    )


def test_uncertainty_matches_independent_finite_difference_jacobian():
    energy, flux, field = 5.0, 3.0, 46.0
    energy_std, flux_std, field_std = 0.08, 0.12, 0.25
    energy_flux_covariance = 0.003
    actual = fang_2025_conductance_uncertainty(
        energy, flux, field,
        mean_energy_std=energy_std,
        energy_flux_std=flux_std,
        magnetic_field_std=field_std,
        energy_flux_covariance=energy_flux_covariance
    )

    covariance = np.array([
        [energy_std**2, energy_flux_covariance, 0.0],
        [energy_flux_covariance, flux_std**2, 0.0],
        [0.0, 0.0, field_std**2]
    ])
    steps = (1e-5, 1e-5, 1e-4)
    gradients = {}
    for name in ("P", "H"):
        gradients[name] = np.array([
            central_derivative(
                lambda value: literal_fang(value, flux, field)[name],
                energy, steps[0]
            ),
            central_derivative(
                lambda value: literal_fang(energy, value, field)[name],
                flux, steps[1]
            ),
            central_derivative(
                lambda value: literal_fang(energy, flux, value)[name],
                field, steps[2]
            )
        ])

    expected_p_variance = gradients["P"] @ covariance @ gradients["P"]
    expected_h_variance = gradients["H"] @ covariance @ gradients["H"]
    expected_ph_covariance = gradients["P"] @ covariance @ gradients["H"]
    np.testing.assert_allclose(actual["P_std"], np.sqrt(expected_p_variance), rtol=2e-9)
    np.testing.assert_allclose(actual["H_std"], np.sqrt(expected_h_variance), rtol=2e-9)
    np.testing.assert_allclose(actual["PH_covariance"], expected_ph_covariance, rtol=2e-9)
    assert actual["uncertainty_valid"]


def test_small_error_monte_carlo_matches_first_order_covariance():
    rng = np.random.default_rng(2025)
    nominal = np.array([7.0, 4.0])
    energy_std, flux_std, field_std = 0.025, 0.04, 0.08
    energy_flux_covariance = 0.0003
    input_covariance = np.array([
        [energy_std**2, energy_flux_covariance],
        [energy_flux_covariance, flux_std**2]
    ])
    energy_flux_samples = rng.multivariate_normal(
        nominal, input_covariance, size=100_000
    )
    field_samples = rng.normal(52.0, field_std, size=energy_flux_samples.shape[0])
    samples = fang_2025_conductance(
        energy_flux_samples[:, 0], energy_flux_samples[:, 1], field_samples
    )
    analytical = fang_2025_conductance_uncertainty(
        nominal[0], nominal[1], 52.0,
        mean_energy_std=energy_std,
        energy_flux_std=flux_std,
        magnetic_field_std=field_std,
        energy_flux_covariance=energy_flux_covariance
    )
    sampled_covariance = np.cov(np.stack([samples["P"], samples["H"]]), ddof=1)

    np.testing.assert_allclose(
        analytical["P_std"], np.sqrt(sampled_covariance[0, 0]), rtol=0.012
    )
    np.testing.assert_allclose(
        analytical["H_std"], np.sqrt(sampled_covariance[1, 1]), rtol=0.012
    )
    np.testing.assert_allclose(
        analytical["PH_covariance"], sampled_covariance[0, 1], rtol=0.02
    )


def test_floor_derivative_is_zero_and_boundary_crossing_is_flagged():
    below = fang_2025_conductance_uncertainty(
        [0.5, 0.9, 1.0], 2.0, 50.0,
        mean_energy_std=[0.1, 0.2, 0.1],
        energy_flux_std=0.0,
        magnetic_field_std=0.0
    )

    np.testing.assert_array_equal(below["P_std"][:2], [0.0, 0.0])
    np.testing.assert_array_equal(below["H_std"][:2], [0.0, 0.0])
    np.testing.assert_array_equal(
        below["uncertainty_straddles_low_energy_boundary"],
        [False, True, True]
    )
    assert below["P_std"][2] > 0


def test_invalid_uncertainty_does_not_change_nominal_conductance():
    nominal = fang_2025_conductance([5.0, 5.0, 5.0], 2.0, 50.0)
    actual = fang_2025_conductance_uncertainty(
        [5.0, 5.0, 5.0], 2.0, 50.0,
        mean_energy_std=[0.1, -0.1, 0.1],
        energy_flux_std=[0.2, 0.2, 0.2],
        magnetic_field_std=[0.3, 0.3, 0.3],
        energy_flux_covariance=[0.0, 0.0, 0.03]
    )

    np.testing.assert_array_equal(actual["P"], nominal["P"])
    np.testing.assert_array_equal(actual["H"], nominal["H"])
    np.testing.assert_array_equal(actual["uncertainty_valid"], [True, False, False])
    assert np.all(np.isnan(actual["P_std"][1:]))
    assert np.all(np.isnan(actual["H_std"][1:]))


def test_zero_flux_uncertainty_distinguishes_fixed_zero_from_singular_flux():
    fixed = fang_2025_conductance_uncertainty(
        5.0, 0.0, 50.0,
        mean_energy_std=0.2,
        energy_flux_std=0.0,
        magnetic_field_std=0.5
    )
    singular = fang_2025_conductance_uncertainty(
        5.0, 0.0, 50.0,
        mean_energy_std=0.2,
        energy_flux_std=0.1,
        magnetic_field_std=0.5
    )

    assert fixed["uncertainty_valid"]
    assert fixed["P_std"] == 0 and fixed["H_std"] == 0
    assert fixed["PH_covariance"] == 0
    assert not singular["uncertainty_valid"]
    assert np.isnan(singular["P_std"])
