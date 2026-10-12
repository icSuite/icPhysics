"""Conditional parity with IDL equations using a supplied synthetic atmosphere.

The unavailable msis120065.xdr atmosphere and an actual IDL runtime are not
reference inputs here. The scalar transcription isolates float32 Q rounding
when requested; other IDL literal-precision effects are not reproduced.
"""

import math

import numpy as np
import pytest

from icphysics import conductance_from_spectrum
from icphysics.spectrum import _log_interpolate


def synthetic_atmosphere():
    altitude = np.arange(60.0, 311.0, 10.0)
    return {
        "altitude_km": altitude,
        "oxygen_cm3": 5e11 * np.exp(-(altitude - 100) / 28),
        "nitrogen_cm3": 1e13 * np.exp(-(altitude - 100) / 8),
        "molecular_oxygen_cm3": 2e12 * np.exp(-(altitude - 100) / 8),
        "mass_density_g_cm3": 6e-10 * np.exp(-(altitude - 100) / 9),
        "temperature_k": 180 + 2.5 * (altitude - 60)
    }


def scalar_idl_reference(energy, flux, field_nt, atmosphere, samples=500, float32_q=False):
    """Transcribe the source's channel/height loops independently."""
    order = np.argsort(energy)
    energies = [float(energy[i]) for i in order]
    fluxes = [float(flux[i]) for i in order]
    field = field_nt * 1e-9
    altitudes = [10 ** (math.log10(70) + j * (math.log10(295) - math.log10(70)) / samples) for j in range(samples)]
    profile_alt = atmosphere["altitude_km"]
    number_flux = []
    for i, energy_i in enumerate(energies):
        lower_half = (energy_i - energies[i - 1]) / 2 if i else (energies[1] - energy_i) / 2
        upper_half = (energies[i + 1] - energy_i) / 2 if i < len(energies) - 1 else lower_half
        bin_flux = fluxes[i] * math.pi * (lower_half + upper_half) * 1000
        number_flux.append(bin_flux * 6.242e11 / (energy_i * 1000))

    p_profile = []
    h_profile = []
    for altitude in altitudes:
        lower = max(0, int(np.searchsorted(profile_alt, altitude, side="left")) - 1)
        fraction = (altitude - profile_alt[lower]) / (profile_alt[lower + 1] - profile_alt[lower])
        values = {}
        for key in ("oxygen_cm3", "nitrogen_cm3", "molecular_oxygen_cm3", "mass_density_g_cm3", "temperature_k"):
            left = atmosphere[key][lower]
            right = atmosphere[key][lower + 1]
            if fraction == 0:
                values[key] = left
            elif fraction == 1:
                values[key] = right
            elif left == 0 or right == 0:
                values[key] = 0.0
            else:
                values[key] = math.exp(math.log(left) + fraction * (math.log(right) - math.log(left)))
        o1 = values["oxygen_cm3"]
        n2 = values["nitrogen_cm3"]
        o2 = values["molecular_oxygen_cm3"]
        temp = values["temperature_k"]
        pressure = temp * 82.05 * (o1 + o2 + n2) / 6.022e23
        depth = 10130 * pressure / 9.82

        q_total = 0.0
        for i, energy_i in enumerate(energies):
            electron_range = 4.3e-7 + 5.36e-6 * energy_i**1.67
            x = depth / electron_range
            if x <= 1:
                deposition = 0.751031 + 2.526201*x + 0.276534*x**2 - 8.344593*x**3 + 4.796477*x**4
                q = number_flux[i] * energy_i * deposition * values["mass_density_g_cm3"] / (0.035 * electron_range)
                q_total += float(np.float32(q)) if float32_q else q
        alpha = 4.3e-6 * math.exp(-2.42e-2 * altitude) + 8.16e12 * math.exp(-0.524 * altitude)
        density = math.sqrt(q_total / alpha)
        sigma_o2 = (2.2 + 5.1 * math.sqrt(3 * 8.63e-5 * temp)) * 1e-20
        sigma_n2 = (18.3 - 7.3 * (3 * 8.63e-5 * temp)) * math.sqrt(3 * 8.63e-5 * temp) * 1e-20
        nue = 1e6 * (o2 * sigma_o2 + n2 * sigma_n2 + o1 * 1e-19) * math.sqrt(8 * 1.381e-23 * temp / (math.pi * 9.11e-31))
        nui = 3.75e-10 * (o2 + n2 + o1)
        oe = 1.602e-19 * field / 9.11e-31
        oi = 1.602e-19 * field / (30 * 1.672e-27)
        p_profile.append(density * 1e6 * 1.602e-19 * (nue*oe / (oe**2 + nue**2) + nui*oi / (oi**2 + nui**2)) / field)
        h_profile.append(density * 1e6 * 1.602e-19 * (oe**2 / (oe**2 + nue**2) - oi**2 / (oi**2 + nui**2)) / field)

    p = 0.0
    h = 0.0
    for j in range(1, samples - 1):
        height_step = (altitudes[j + 1] - altitudes[j - 1]) * 1000 / 2
        p += p_profile[j] * height_step
        h += h_profile[j] * height_step
    return {"P": p, "H": h}


@pytest.mark.parametrize("samples", [37, 500])
def test_matches_independent_scalar_idl_equations(samples):
    energy = np.array([0.03, 0.1, 0.4, 1.1, 3.0, 8.0, 20.0])
    flux = np.array([1e-6, 2e-5, 3e-5, 4e-5, 2e-5, 3e-6, 1e-7])
    atmosphere = synthetic_atmosphere()
    actual = conductance_from_spectrum(energy, flux, 49000, atmosphere, altitude_samples=samples)
    reference = scalar_idl_reference(energy, flux, 49000, atmosphere, samples=samples)
    for name in ("P", "H"):
        np.testing.assert_allclose(actual[name], reference[name], rtol=2e-14)


def test_float64_difference_from_source_float32_q_is_small():
    energy = np.geomspace(0.03, 30, 19)
    flux = np.geomspace(1e-4, 1e-8, 19)
    atmosphere = synthetic_atmosphere()
    actual = conductance_from_spectrum(energy, flux, 50000, atmosphere)
    rounded = scalar_idl_reference(energy, flux, 50000, atmosphere, float32_q=True)
    for name in ("P", "H"):
        np.testing.assert_allclose(actual[name], rounded[name], rtol=5e-8)


def test_batch_broadcast_zero_flux_sqrt_amplitude_and_channel_order():
    energy = np.array([0.2, 1, 4, 10])
    flux = np.array([1e-5, 3e-5, 2e-5, 1e-6])
    atmosphere = synthetic_atmosphere()
    batch_flux = np.array([0, 1, 9])[:, None, None] * flux
    fields = np.array([40000, 55000])
    actual = conductance_from_spectrum(energy, batch_flux, fields, atmosphere)
    reversed_order = conductance_from_spectrum(energy[::-1], batch_flux[..., ::-1], fields, atmosphere)
    for name in ("P", "H"):
        assert actual[name].shape == (3, 2)
        np.testing.assert_array_equal(actual[name][0], [0, 0])
        np.testing.assert_allclose(actual[name][2], 3 * actual[name][1], rtol=1e-14)
        np.testing.assert_array_equal(actual[name], reversed_order[name])
        for j, field in enumerate(fields):
            scalar = scalar_idl_reference(energy, flux, field, atmosphere)
            np.testing.assert_allclose(actual[name][1, j], scalar[name], rtol=2e-14)


def test_invalid_spectra_and_field_do_not_become_partial_spectra():
    energy = np.array([0.1, 1, 10])
    spectra = np.array([[1e-5, 2e-5, 1e-6], [np.nan, 2e-5, 1e-6], [1e-5, -1, 1e-6], [1e-5, 2e-5, np.inf]])
    actual = conductance_from_spectrum(energy, spectra, 50000, synthetic_atmosphere())
    for name in ("P", "H"):
        assert np.isfinite(actual[name][0])
        assert np.all(np.isnan(actual[name][1:]))
    masked = np.ma.array(spectra[0], mask=[False, True, False])
    actual = conductance_from_spectrum(energy, masked, 50000, synthetic_atmosphere())
    assert np.isnan(actual["P"]) and np.isnan(actual["H"])
    masked_integer = np.ma.array([0, 0, 0], mask=[False, True, False])
    actual = conductance_from_spectrum(energy, masked_integer, 50000, synthetic_atmosphere())
    assert np.isnan(actual["P"]) and np.isnan(actual["H"])
    actual = conductance_from_spectrum(energy, spectra[0], [0, -1, np.nan, np.inf], synthetic_atmosphere())
    assert np.all(np.isnan(actual["P"])) and np.all(np.isnan(actual["H"]))


def test_inputs_are_not_mutated():
    energy = np.array([10.0, 1.0, 0.1])
    flux = np.array([[1e-6, 2e-5, 1e-5]])
    field = np.array([50000.0])
    atmosphere = synthetic_atmosphere()
    originals = [a.copy() for a in (energy, flux, field, *atmosphere.values())]
    conductance_from_spectrum(energy, flux, field, atmosphere)
    for actual, original in zip((energy, flux, field, *atmosphere.values()), originals):
        np.testing.assert_array_equal(actual, original)


@pytest.mark.parametrize("energy", [[1], [0, 1], [1, 1], [1, np.nan], [1, 3, 2], [[1, 2]]])
def test_invalid_energy_grid_is_rejected(energy):
    with pytest.raises(ValueError, match="energy_keV"):
        conductance_from_spectrum(energy, [1e-5, 1e-5], 50000, synthetic_atmosphere())


@pytest.mark.parametrize("samples", [2, 3.5, True])
def test_invalid_altitude_samples_are_rejected(samples):
    with pytest.raises(ValueError, match="altitude_samples"):
        conductance_from_spectrum([1, 2], [1e-5, 1e-5], 50000, synthetic_atmosphere(), altitude_samples=samples)


def test_invalid_atmosphere_and_channel_shape_are_rejected():
    atmosphere = synthetic_atmosphere()
    with pytest.raises(ValueError, match="final axis"):
        conductance_from_spectrum([1, 2], [1e-5], 50000, atmosphere)
    atmosphere["nitrogen_cm3"][0] = -1
    with pytest.raises(ValueError, match="nitrogen_cm3"):
        conductance_from_spectrum([1, 2], [1e-5, 1e-5], 50000, atmosphere)
    atmosphere = synthetic_atmosphere()
    atmosphere["altitude_km"][-1] = 290
    with pytest.raises(ValueError, match="altitude_km"):
        conductance_from_spectrum([1, 2], [1e-5, 1e-5], 50000, atmosphere)


def test_zero_species_use_geometric_limit_with_exact_endpoints():
    source_altitude = np.array([70.0, 80.0, 90.0, 100.0])
    source_density = np.array([0.0, 4.0, 16.0, 0.0])
    altitude = np.array([70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0])
    expected = np.array([0.0, 0.0, 4.0, 8.0, 16.0, 0.0, 0.0])
    np.testing.assert_allclose(_log_interpolate(altitude, source_altitude, source_density), expected)


def test_msis_zero_low_altitude_oxygen_remains_finite_and_matches_scalar_reference():
    energy = np.array([0.1, 1.0, 10.0])
    flux = np.array([1e-5, 3e-5, 2e-6])
    atmosphere = synthetic_atmosphere()
    atmosphere["oxygen_cm3"][atmosphere["altitude_km"] < 72.5] = 0.0
    actual = conductance_from_spectrum(energy, flux, 50000, atmosphere)
    reference = scalar_idl_reference(energy, flux, 50000, atmosphere)
    for name in ("P", "H"):
        assert np.isfinite(actual[name]) and actual[name] > 0
        np.testing.assert_allclose(actual[name], reference[name], rtol=2e-14)


@pytest.mark.parametrize("name", ["mass_density_g_cm3", "temperature_k"])
def test_mass_density_and_temperature_cannot_be_zero(name):
    atmosphere = synthetic_atmosphere()
    atmosphere[name][0] = 0.0
    with pytest.raises(ValueError, match=name):
        conductance_from_spectrum([1, 2], [1e-5, 1e-5], 50000, atmosphere)


def test_different_batched_atmospheres_match_per_profile_and_scalar_calculations():
    energy = np.array([0.05, 0.2, 0.8, 3.0, 12.0])
    base_flux = np.array([1e-5, 2e-5, 4e-5, 1e-5, 1e-6])
    flux = np.array([1.0, 2.5])[:, None, None] * base_flux
    fields = np.array([41000.0, 49000.0, 56000.0])
    base = synthetic_atmosphere()
    atmosphere = {name: values.copy() for name, values in base.items()}
    atmosphere["oxygen_cm3"] = np.array([0.6, 1.8])[:, None, None] * base["oxygen_cm3"]
    atmosphere["oxygen_cm3"][0, 0, base["altitude_km"] < 90] = 0
    atmosphere["nitrogen_cm3"] = np.array([1.4, 0.7])[:, None, None] * base["nitrogen_cm3"]
    atmosphere["mass_density_g_cm3"] = np.array([0.8, 1.3])[:, None, None] * base["mass_density_g_cm3"]
    atmosphere["temperature_k"] = np.array([0.9, 1.0, 1.2])[None, :, None] * base["temperature_k"]
    originals = {name: values.copy() for name, values in atmosphere.items()}
    actual = conductance_from_spectrum(energy, flux, fields, atmosphere)
    for name in ("P", "H"):
        assert actual[name].shape == (2, 3)
        assert np.all(np.isfinite(actual[name]))
    for i in range(2):
        for j in range(3):
            profile = {"altitude_km": base["altitude_km"]}
            for name, values in atmosphere.items():
                if name != "altitude_km":
                    profile[name] = np.broadcast_to(values, (2, 3, base["altitude_km"].size))[i, j]
            shared_profile = conductance_from_spectrum(energy, flux[i, 0], fields[j], profile)
            scalar = scalar_idl_reference(energy, flux[i, 0], fields[j], profile)
            for name in ("P", "H"):
                np.testing.assert_allclose(actual[name][i, j], shared_profile[name], rtol=2e-14)
                np.testing.assert_allclose(actual[name][i, j], scalar[name], rtol=2e-14)
    for name, original in originals.items():
        np.testing.assert_array_equal(atmosphere[name], original)


def test_batched_atmosphere_interpolation_preserves_zero_limits_and_endpoints():
    profile_alt = np.array([70.0, 80.0, 90.0])
    values = np.array([[0.0, 4.0, 0.0], [4.0, 16.0, 64.0]])
    targets = np.array([70.0, 75.0, 80.0, 85.0, 90.0])
    expected = np.array([[0.0, 0.0, 4.0, 0.0, 0.0], [4.0, 8.0, 16.0, 32.0, 64.0]])
    np.testing.assert_allclose(_log_interpolate(targets, profile_alt, values), expected)


def test_invalid_spectrum_masks_remain_invalid_with_batched_atmosphere():
    atmosphere = synthetic_atmosphere()
    for name in atmosphere:
        if name != "altitude_km":
            atmosphere[name] = np.stack([atmosphere[name], atmosphere[name] * 1.1])
    spectra = np.ma.array([[1e-5, 2e-5], [1e-5, 2e-5]], mask=[[False, False], [False, True]])
    actual = conductance_from_spectrum([0.5, 3], spectra, 50000, atmosphere)
    for name in ("P", "H"):
        assert np.isfinite(actual[name][0])
        assert np.isnan(actual[name][1])


def test_masked_atmosphere_profiles_raise_instead_of_becoming_zero():
    atmosphere = synthetic_atmosphere()
    profile = np.stack([atmosphere["oxygen_cm3"], atmosphere["oxygen_cm3"]])
    profile = np.ma.array(profile, mask=np.zeros_like(profile, dtype=bool))
    profile.mask[1, 5] = True
    atmosphere["oxygen_cm3"] = profile
    with pytest.raises(ValueError, match="oxygen_cm3"):
        conductance_from_spectrum([0.5, 3], [1e-5, 2e-5], 50000, atmosphere)


@pytest.mark.parametrize("target", ["spectrum", "atmosphere"])
def test_incompatible_leading_dimensions_raise_descriptive_error(target):
    atmosphere = synthetic_atmosphere()
    flux = np.full((2, 2), 1e-5)
    field = np.full(3, 50000) if target == "spectrum" else 50000
    if target == "atmosphere":
        atmosphere["oxygen_cm3"] = np.stack([atmosphere["oxygen_cm3"]] * 3)
    with pytest.raises(ValueError, match="dimensions must broadcast together"):
        conductance_from_spectrum([0.5, 3], flux, field, atmosphere)
