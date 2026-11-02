"""Synthetic checks for the three Zhang--Paxton latitude collapses."""

import numpy as np

from icphysics.robinson import hall, ped
from icphysics.zhang_paxton import (
    _q_weighted_mlt_profile,
    collapsed_zhang_paxton_q_weighted,
    collapse_zhang_paxton,
    collapse_latitude_slice,
    latitude_cell_area_weights,
)
import icphysics.zhang_paxton as zhang_paxton_module


def test_area_and_q_weighted_means_and_structural_spreads():
    energy = np.array([1.0, 3.0, 6.0])
    flux = np.array([1.0, 2.0, 4.0])
    edges = np.array([50.0, 60.0, 70.0, 80.0])
    area = latitude_cell_area_weights(edges)

    result = collapse_latitude_slice(energy, flux, edges, threshold=0.05)

    area_mean = np.average(energy, weights=area)
    area_spread = np.sqrt(np.average((energy - area_mean) ** 2, weights=area))
    q_weights = area * flux
    q_mean = np.average(energy, weights=q_weights)
    q_spread = np.sqrt(np.average((energy - q_mean) ** 2, weights=q_weights))

    np.testing.assert_allclose(result["area_mean_energy"], area_mean)
    np.testing.assert_allclose(result["area_mean_energy_spread"], area_spread)
    np.testing.assert_allclose(result["q_weighted_energy"], q_mean)
    np.testing.assert_allclose(result["q_weighted_energy_spread"], q_spread)

    # Historical names remain exact area-mean aliases.
    assert result["representative_energy"] == result["area_mean_energy"]
    assert result["weighted_spread"] == result["area_mean_energy_spread"]


def test_conductance_equivalent_collapse_round_trips_area_mean_conductance():
    energy = np.array([0.8, 2.5, 7.0])
    flux = np.array([0.2, 1.4, 3.0])
    edges = np.array([50.0, 58.0, 69.0, 82.0])
    area = latitude_cell_area_weights(edges)

    result = collapse_latitude_slice(energy, flux, edges, threshold=0.05)
    expected_hall = np.average(hall(energy, flux), weights=area)
    expected_pedersen = np.average(ped(energy, flux), weights=area)

    np.testing.assert_allclose(
        hall(
            result["conductance_effective_energy"],
            result["conductance_effective_energy_flux"],
        ),
        expected_hall,
    )
    np.testing.assert_allclose(
        ped(
            result["conductance_effective_energy"],
            result["conductance_effective_energy_flux"],
        ),
        expected_pedersen,
    )


def test_conductance_structural_spread_is_centered_on_effective_energy():
    energy = np.array([1.0, 4.0, 8.0])
    flux = np.array([0.25, 1.0, 9.0])
    edges = np.array([50.0, 60.0, 70.0, 80.0])
    area = latitude_cell_area_weights(edges)

    result = collapse_latitude_slice(energy, flux, edges, threshold=0.05)
    effective_energy = result["conductance_effective_energy"]
    weights = area * np.sqrt(flux)
    expected_spread = np.sqrt(
        np.average((energy - effective_energy) ** 2, weights=weights)
    )

    np.testing.assert_allclose(
        result["conductance_effective_energy_spread"], expected_spread
    )
    assert not np.isclose(effective_energy, np.average(energy, weights=weights))


def test_selection_is_shared_by_all_three_methods():
    energy = np.array([50.0, 1.0, 3.0, 40.0])
    flux = np.array([0.01, 1.0, 2.0, 0.01])
    edges = np.array([50.0, 60.0, 70.0, 80.0, 90.0])

    result = collapse_latitude_slice(energy, flux, edges, threshold=0.05)

    np.testing.assert_array_equal(result["selected"], [False, True, True, False])
    assert 1.0 <= result["area_mean_energy"] <= 3.0
    assert 1.0 <= result["q_weighted_energy"] <= 3.0
    assert np.isfinite(result["conductance_effective_energy"])
    assert np.isfinite(result["conductance_effective_energy_flux"])


def test_batched_collapse_exposes_all_method_shapes(monkeypatch):
    def synthetic_model(kp, mlt, mlat):
        shape = np.broadcast_shapes(kp.shape, mlt.shape, mlat.shape)
        latitude = np.broadcast_to(mlat, shape)
        energy = 1 + (latitude - latitude.min(axis=-1, keepdims=True)) / 20
        flux = np.ones(shape)
        return energy, flux

    monkeypatch.setattr(zhang_paxton_module, "zhang_paxton", synthetic_model)
    result = collapse_zhang_paxton(
        kp=np.array([[1.0], [2.0]]),
        mlt=np.array([[0.0, 6.0, 12.0]]),
        lower_mlat=50,
        upper_mlat=90,
        latitude_step=10,
    )

    for name in (
        "area_mean_energy",
        "area_mean_energy_spread",
        "q_weighted_energy",
        "q_weighted_energy_spread",
        "conductance_effective_energy",
        "conductance_effective_energy_spread",
        "conductance_effective_energy_flux",
    ):
        assert result[name].shape == (2, 3)
        assert np.isfinite(result[name]).all()


def test_cached_q_weighted_evaluator_is_periodic_and_rounds_kp(monkeypatch):
    def synthetic_collapse(kp, mlt, **settings):
        return {
            "q_weighted_energy": kp + mlt / 24,
            "q_weighted_energy_spread": np.full(np.shape(mlt), 0.4),
        }

    monkeypatch.setattr(
        zhang_paxton_module, "collapse_zhang_paxton", synthetic_collapse
    )
    _q_weighted_mlt_profile.cache_clear()
    result = collapsed_zhang_paxton_q_weighted(
        [1.519, 1.519, np.nan],
        [-1.0, 23.0, 2.0],
        latitude_step=10,
        mlt_step=1,
    )
    _q_weighted_mlt_profile.cache_clear()

    np.testing.assert_allclose(result["energy"][:2], [1.52 + 23 / 24] * 2)
    np.testing.assert_allclose(result["energy_spread"][:2], 0.4)
    np.testing.assert_allclose(result["kp"][:2], 1.52)
    assert np.isnan(result["energy"][2])
    np.testing.assert_array_equal(result["valid"], [True, True, False])
