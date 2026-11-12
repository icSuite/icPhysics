"""Small regression check for the first copied physics routines."""

import numpy as np

from icphysics import (
    hall, halluncertainty, ped, peduncertainty, robinson_conductance,
    robinson_effective_precipitation,
)


def test_robinson_reference_values():
    inputs = (2.5, 4.0, 0.3, 0.5, -0.02)

    np.testing.assert_allclose(ped(*inputs[:2]), 8.98876404494382)
    np.testing.assert_allclose(hall(*inputs[:2]), 8.813765704000613)
    np.testing.assert_allclose(
        peduncertainty(*inputs), 0.6842658428542764
    )
    np.testing.assert_allclose(
        halluncertainty(*inputs), 1.3998674616220634
    )


def test_zero_flux_uses_one_sided_excursion():
    E0 = 4.0
    dFe = 0.25

    assert peduncertainty(E0, 0.0, 1.0, dFe, 0.0) == ped(E0, dFe)
    assert halluncertainty(E0, 0.0, 1.0, dFe, 0.0) == hall(E0, dFe)


def test_array_entry_point_matches_scalar_equations():
    E0 = np.array([2.0, 4.0, 3.0])
    Fe = np.array([1.5, 0.0, np.nan])
    dE0 = np.array([0.2, 0.3, 0.4])
    dFe = np.array([0.1, 0.2, 0.3])
    covariance = np.array([-0.01, 0.0, 0.0])

    result = robinson_conductance(E0, Fe, dE0, dFe, covariance)

    np.testing.assert_allclose(result["P"][:2], ped(E0[:2], Fe[:2]))
    np.testing.assert_allclose(result["H"][:2], hall(E0[:2], Fe[:2]))
    np.testing.assert_allclose(
        result["dP"][0],
        peduncertainty(E0[0], Fe[0], dE0[0], dFe[0], covariance[0]),
    )
    np.testing.assert_allclose(result["dP"][1], ped(E0[1], dFe[1]))
    assert np.isnan(result["dH"][2])


def test_robinson_effective_precipitation_round_trip_grid():
    energy = np.geomspace(.05, 100, 31)[:, None]
    energy_flux = np.geomspace(1e-4, 1e3, 29)[None, :]
    result = robinson_effective_precipitation(
        hall(energy, energy_flux), ped(energy, energy_flux)
    )
    assert result["invertible"].all()
    expected_energy = np.broadcast_to(energy, result["effective_energy"].shape)
    expected_flux = np.broadcast_to(energy_flux, result["effective_energy_flux"].shape)
    np.testing.assert_allclose(result["effective_energy"], expected_energy, rtol=1e-15)
    np.testing.assert_allclose(
        result["effective_energy_flux"], expected_flux, rtol=4e-15
    )


def test_robinson_effective_precipitation_broadcast_masks_and_invalid_inputs():
    hall_input = np.ma.array([[1.0], [2.0], [3.0]], mask=[[False], [True], [False]])
    pedersen_input = np.array([[1.0, 0.0, -1.0, np.nan]])
    result = robinson_effective_precipitation(hall_input, pedersen_input)
    assert result["effective_energy"].shape == (3, 4)
    expected = np.array([
        [True, False, False, False],
        [False, False, False, False],
        [True, False, False, False]
    ])
    np.testing.assert_array_equal(result["invertible"], expected)
    assert np.isfinite(result["effective_energy"][expected]).all()
    assert np.isnan(result["effective_energy"][~expected]).all()
    assert np.isnan(result["effective_energy_flux"][~expected]).all()

    scalar = robinson_effective_precipitation(2.0, 1.0)
    assert scalar["effective_energy"].shape == ()
    assert scalar["invertible"].shape == ()


def test_robinson_effective_precipitation_does_not_clip_extreme_ratio():
    result = robinson_effective_precipitation(1e12, 1e-6)
    expected_energy = (1e12 / (0.45e-6)) ** (1.0 / .85)
    np.testing.assert_allclose(result["effective_energy"], expected_energy)
    assert result["effective_energy"] > 1e20
    assert np.isfinite(result["effective_energy_flux"])
    assert result["invertible"]
