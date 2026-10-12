"""Independent checks of conditional electron-spectrum measurement errors.

Neutral atmosphere, field, channel centers and angular/model assumptions are
fixed. These tests do not validate their physical uncertainty or establish
original MSIS-file/IDL-runtime parity.
"""

import numpy as np
import pytest

from icphysics import conductance_from_spectrum
from test_spectrum import scalar_idl_reference, synthetic_atmosphere


def spectrum_inputs():
    energy = np.array([.03, .1, .4, 1.1, 3., 8., 20.])
    flux = np.array([1e-5, 2e-5, 3e-5, 2e-5, 1e-5, 3e-6, 1e-7])
    return energy, flux, synthetic_atmosphere()


def scalar_finite_difference_covariance(energy, flux, standard_deviation, field, atmosphere, samples=500):
    """Differentiate the independent scalar IDL transcription numerically."""
    derivative = np.zeros((2, energy.size))
    for channel in range(energy.size):
        step = flux[channel] * 1e-4
        high, low = flux.copy(), flux.copy()
        high[channel] += step
        low[channel] -= step
        high_result = scalar_idl_reference(energy, high, field, atmosphere, samples=samples)
        low_result = scalar_idl_reference(energy, low, field, atmosphere, samples=samples)
        for component, name in enumerate(("P", "H")):
            derivative[component, channel] = (high_result[name] - low_result[name]) / (2 * step)
    weighted = derivative * standard_deviation
    return weighted @ weighted.T


def test_analytical_uncertainty_matches_independent_scalar_finite_differences():
    energy, flux, atmosphere = spectrum_inputs()
    std = flux * np.array([.01, .03, .02, .04, .015, .025, .05])
    nominal = conductance_from_spectrum(energy, flux, 49000, atmosphere)
    result = conductance_from_spectrum(energy, flux, 49000, atmosphere, differential_energy_flux_std=std)
    expected = scalar_finite_difference_covariance(energy, flux, std, 49000, atmosphere)
    assert set(nominal) == {"P", "H"}
    for name in ("P", "H"):
        np.testing.assert_array_equal(result[name], nominal[name])
    assert result["uncertainty_valid"]
    np.testing.assert_allclose(result["P_std"]**2, expected[0, 0], rtol=5e-7)
    np.testing.assert_allclose(result["H_std"]**2, expected[1, 1], rtol=5e-7)
    np.testing.assert_allclose(result["PH_covariance"], expected[0, 1], rtol=5e-7)
    assert 0 < result["PH_covariance"] <= result["P_std"] * result["H_std"] * (1 + 1e-14)


def test_small_independent_channel_errors_match_joint_monte_carlo():
    energy, flux, atmosphere = spectrum_inputs()
    std = flux * .005
    result = conductance_from_spectrum(energy, flux, 50000, atmosphere, altitude_samples=101, differential_energy_flux_std=std)
    rng = np.random.default_rng(32017)
    draws = rng.normal(flux, std, size=(12000, energy.size))
    assert np.all(draws > 0)  # no clipping or censored-input approximation
    samples = conductance_from_spectrum(energy, draws, 50000, atmosphere, altitude_samples=101)
    covariance = np.cov(np.stack([samples["P"], samples["H"]]), ddof=1)
    np.testing.assert_allclose(result["P_std"], np.sqrt(covariance[0, 0]), rtol=.025)
    np.testing.assert_allclose(result["H_std"], np.sqrt(covariance[1, 1]), rtol=.025)
    np.testing.assert_allclose(result["PH_covariance"], covariance[0, 1], rtol=.04)


def test_one_uncertain_channel_has_perfect_joint_correlation_and_linear_error_scaling():
    energy, flux, atmosphere = spectrum_inputs()
    std = np.zeros_like(flux)
    std[3] = .02 * flux[3]
    result = conductance_from_spectrum(energy, flux, 50000, atmosphere, differential_energy_flux_std=std)
    doubled = conductance_from_spectrum(energy, flux, 50000, atmosphere, differential_energy_flux_std=2 * std)
    np.testing.assert_allclose(result["PH_covariance"], result["P_std"] * result["H_std"], rtol=1e-14)
    for name in ("P_std", "H_std"):
        np.testing.assert_allclose(doubled[name], 2 * result[name], rtol=1e-14)
    np.testing.assert_allclose(doubled["PH_covariance"], 4 * result["PH_covariance"], rtol=1e-14)


def test_reversal_shared_profile_and_batched_profiles_preserve_uncertainty():
    energy, flux, base = spectrum_inputs()
    atmosphere = {name: values if name == "altitude_km" else np.stack([values, values * 1.2])[:, None, :] for name, values in base.items()}
    spectra = np.array([1., 3., 9.])[None, :, None] * flux
    std = flux * .02
    fields = np.array([44000, 49000, 55000])
    result = conductance_from_spectrum(energy, spectra, fields, atmosphere, differential_energy_flux_std=std)
    reverse = conductance_from_spectrum(energy[::-1], spectra[..., ::-1], fields, atmosphere, differential_energy_flux_std=std[::-1])
    for name in ("P", "H", "P_std", "H_std", "PH_covariance", "uncertainty_valid"):
        assert result[name].shape == (2, 3)
        np.testing.assert_allclose(result[name], reverse[name], rtol=1e-14)
    for row in range(2):
        profile = {name: values if name == "altitude_km" else values[row, 0] for name, values in atmosphere.items()}
        for column in range(3):
            separate = conductance_from_spectrum(energy, spectra[0, column], fields[column], profile, differential_energy_flux_std=std)
            for name in result:
                np.testing.assert_allclose(result[name][row, column], separate[name], rtol=2e-14)


@pytest.mark.parametrize("invalid", [-1., np.inf, -np.inf])
def test_invalid_standard_deviation_preserves_nominal_conductance(invalid):
    energy, flux, atmosphere = spectrum_inputs()
    std = np.stack([flux * .01, flux * .02])
    std[1, 3] = invalid
    spectra = np.stack([flux, flux])
    nominal = conductance_from_spectrum(energy, spectra, 50000, atmosphere)
    result = conductance_from_spectrum(energy, spectra, 50000, atmosphere, differential_energy_flux_std=std)
    for name in ("P", "H"):
        np.testing.assert_array_equal(result[name], nominal[name])
    np.testing.assert_array_equal(result["uncertainty_valid"], [True, False])
    for name in ("P_std", "H_std", "PH_covariance"):
        assert np.isfinite(result[name][0]) and np.isnan(result[name][1])
    np.testing.assert_array_equal(result["uncertainty_channel_count"], [energy.size, energy.size - 1])


def test_masked_uncertainty_is_missing_and_inputs_are_unchanged():
    energy, flux, atmosphere = spectrum_inputs()
    std = np.ma.array(flux * .03, mask=np.arange(flux.size) == 2)
    original_flux, original_std = flux.copy(), std.copy()
    result = conductance_from_spectrum(energy[::-1], flux[::-1], 50000, atmosphere, differential_energy_flux_std=std[::-1])
    omitted = std.filled(0.)
    reference = conductance_from_spectrum(energy, flux, 50000, atmosphere, differential_energy_flux_std=omitted)
    assert np.isfinite(result["P"]) and np.isfinite(result["H"])
    assert result["uncertainty_valid"]
    assert result["uncertainty_channel_count"] == energy.size - 1
    for name in ("P", "H", "P_std", "H_std", "PH_covariance"):
        np.testing.assert_array_equal(result[name], reference[name])
    np.testing.assert_array_equal(flux, original_flux)
    np.testing.assert_array_equal(std.data, original_std.data)
    np.testing.assert_array_equal(std.mask, original_std.mask)


def test_zero_signal_zero_error_is_valid_but_depositing_error_is_singular():
    energy, flux, atmosphere = spectrum_inputs()
    zeros = np.zeros_like(flux)
    zero_error = conductance_from_spectrum(energy, zeros, 50000, atmosphere, differential_energy_flux_std=zeros)
    positive_error = conductance_from_spectrum(energy, zeros, 50000, atmosphere, differential_energy_flux_std=flux * .01)
    assert zero_error["uncertainty_valid"]
    for name in ("P", "H", "P_std", "H_std", "PH_covariance"):
        assert zero_error[name] == 0
    assert not positive_error["uncertainty_valid"]
    assert positive_error["P"] == positive_error["H"] == 0
    assert np.isnan(positive_error["P_std"]) and np.isnan(positive_error["H_std"])


def test_error_only_in_nondepositing_channel_does_not_create_singularity():
    # This supplied flat profile has depth between the 1- and 10-keV ranges.
    altitude = np.array([70., 295.])
    atmosphere = {
        "altitude_km": altitude, "oxygen_cm3": np.full(2, 1e12),
        "nitrogen_cm3": np.full(2, 1e12), "molecular_oxygen_cm3": np.full(2, 1e12),
        "mass_density_g_cm3": np.full(2, 1e-10), "temperature_k": np.full(2, 300.)
    }
    for flux in ([0., 0.], [0., 1e-5]):
        result = conductance_from_spectrum([1., 10.], flux, 50000, atmosphere, differential_energy_flux_std=[1e-6, 0.])
        assert result["uncertainty_valid"]
        for name in ("P_std", "H_std", "PH_covariance"):
            assert result[name] == 0


def test_zero_ionization_in_added_uncertain_deposition_region_is_singular():
    atmosphere = synthetic_atmosphere()
    nominal = conductance_from_spectrum([1., 10.], [1e-5, 0.], 50000, atmosphere)
    result = conductance_from_spectrum([1., 10.], [1e-5, 0.], 50000, atmosphere, differential_energy_flux_std=[1e-6, 1e-6])
    assert nominal["P"] > 0 and nominal["H"] > 0
    np.testing.assert_array_equal(result["P"], nominal["P"])
    assert not result["uncertainty_valid"]
    assert np.isnan(result["P_std"]) and np.isnan(result["PH_covariance"])


def test_invalid_nominal_input_cannot_have_valid_uncertainty():
    energy, flux, atmosphere = spectrum_inputs()
    flux[2] = np.nan
    result = conductance_from_spectrum(energy, flux, 50000, atmosphere, differential_energy_flux_std=np.zeros_like(flux))
    assert not result["uncertainty_valid"]
    for name in ("P", "H", "P_std", "H_std", "PH_covariance"):
        assert np.isnan(result[name])


def test_uncertain_ionization_only_at_omitted_endpoint_does_not_invalidate_integral():
    # With three grid samples only the middle height enters the integral.
    # This synthetic pressure profile allows deposition only at the last
    # (omitted) grid height; uncertainty there cannot affect conductance.
    altitude = np.array([70., 295.])
    neutral = 1e14 * np.exp(-(altitude - 70) / 7.5)
    atmosphere = {
        "altitude_km": altitude, "oxygen_cm3": neutral,
        "nitrogen_cm3": neutral, "molecular_oxygen_cm3": neutral,
        "mass_density_g_cm3": 1e-9 * np.exp(-(altitude - 70) / 7.5),
        "temperature_k": np.full(2, 300.)
    }
    result = conductance_from_spectrum([1., 2.], [0., 0.], 50000, atmosphere, altitude_samples=3, differential_energy_flux_std=[1e-6, 1e-6])
    assert result["uncertainty_valid"]
    for name in ("P", "H", "P_std", "H_std", "PH_covariance"):
        assert result[name] == 0


def test_wrong_uncertainty_channel_count_is_rejected():
    energy, flux, atmosphere = spectrum_inputs()
    with pytest.raises(ValueError, match="(channel|final axis|energy_keV)"):
        conductance_from_spectrum(energy, flux, 50000, atmosphere, differential_energy_flux_std=np.ones(energy.size - 1))


def test_uncertainty_cannot_change_nominal_batch_dimensions():
    energy, flux, atmosphere = spectrum_inputs()
    with pytest.raises(ValueError, match="broadcast"):
        conductance_from_spectrum(energy, flux, 50000, atmosphere, differential_energy_flux_std=np.ones((2, energy.size)))


def test_partial_nan_errors_match_omitted_variance_terms_and_keep_full_nominal_spectrum():
    energy, flux, atmosphere = spectrum_inputs()
    std = flux * np.array([.01, .03, .02, .04, .015, .025, .05])
    std[[1, 4, 6]] = np.nan
    omitted = np.where(np.isnan(std), 0., std)
    nominal = conductance_from_spectrum(energy, flux, 49000, atmosphere)
    reference = conductance_from_spectrum(energy, flux, 49000, atmosphere, differential_energy_flux_std=omitted)
    result = conductance_from_spectrum(energy, flux, 49000, atmosphere, differential_energy_flux_std=std)
    independent = scalar_finite_difference_covariance(energy, flux, omitted, 49000, atmosphere)
    for name in ("P", "H"):
        np.testing.assert_array_equal(result[name], nominal[name])
    for name in ("P_std", "H_std", "PH_covariance"):
        np.testing.assert_array_equal(result[name], reference[name])
    np.testing.assert_allclose(result["P_std"]**2, independent[0, 0], rtol=5e-7)
    np.testing.assert_allclose(result["H_std"]**2, independent[1, 1], rtol=5e-7)
    np.testing.assert_allclose(result["PH_covariance"], independent[0, 1], rtol=5e-7)
    assert result["uncertainty_valid"]
    assert result["uncertainty_channel_count"] == 4
    assert result["uncertainty_channel_count"].dtype == np.int16


@pytest.mark.parametrize("zero_flux", [False, True])
def test_all_missing_errors_unavailable_but_one_reported_zero_can_give_zero_error(zero_flux):
    energy, flux, atmosphere = spectrum_inputs()
    if zero_flux:
        flux = np.zeros_like(flux)
    nominal = conductance_from_spectrum(energy, flux, 50000, atmosphere)
    std = np.full_like(flux, np.nan)
    missing = conductance_from_spectrum(energy, flux, 50000, atmosphere, differential_energy_flux_std=std)
    assert not missing["uncertainty_valid"]
    assert missing["uncertainty_channel_count"] == 0
    for name in ("P_std", "H_std", "PH_covariance"):
        assert np.isnan(missing[name])
    std[2] = 0.
    one_reported = conductance_from_spectrum(energy, flux, 50000, atmosphere, differential_energy_flux_std=std)
    assert one_reported["uncertainty_valid"]
    assert one_reported["uncertainty_channel_count"] == 1
    for name in ("P_std", "H_std", "PH_covariance"):
        assert one_reported[name] == 0
    for name in ("P", "H"):
        np.testing.assert_array_equal(missing[name], nominal[name])
        np.testing.assert_array_equal(one_reported[name], nominal[name])


def test_missing_errors_in_batched_profiles_preserve_reversal_and_channel_counts():
    energy, flux, base = spectrum_inputs()
    atmosphere = {name: values if name == "altitude_km" else np.stack([values, values * 1.2]) for name, values in base.items()}
    spectra = np.stack([flux, flux * 2])
    std = spectra * .02
    std[0, [0, 2]] = np.nan
    std[1, [3, 5, 6]] = np.nan
    result = conductance_from_spectrum(energy, spectra, 50000, atmosphere, differential_energy_flux_std=std)
    reverse = conductance_from_spectrum(energy[::-1], spectra[..., ::-1], 50000, atmosphere, differential_energy_flux_std=std[..., ::-1])
    omitted = conductance_from_spectrum(energy, spectra, 50000, atmosphere, differential_energy_flux_std=np.nan_to_num(std))
    np.testing.assert_array_equal(result["uncertainty_channel_count"], [5, 4])
    assert result["uncertainty_channel_count"].dtype == np.int16
    for name in result:
        np.testing.assert_array_equal(result[name], reverse[name])
        if name != "uncertainty_channel_count":
            np.testing.assert_array_equal(result[name], omitted[name])


def test_negative_or_infinite_error_is_not_hidden_by_nan_omission():
    energy, flux, atmosphere = spectrum_inputs()
    std = np.full((3, energy.size), np.nan)
    std[:, 1] = flux[1] * .02
    std[:, 3] = [-1e-6, np.inf, -np.inf]
    nominal = conductance_from_spectrum(energy, np.stack([flux] * 3), 50000, atmosphere)
    result = conductance_from_spectrum(energy, np.stack([flux] * 3), 50000, atmosphere, differential_energy_flux_std=std)
    np.testing.assert_array_equal(result["uncertainty_channel_count"], [1, 1, 1])
    assert not np.any(result["uncertainty_valid"])
    for name in ("P", "H"):
        np.testing.assert_array_equal(result[name], nominal[name])
    for name in ("P_std", "H_std", "PH_covariance"):
        assert np.isnan(result[name]).all()
