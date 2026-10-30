"""Tests for collapsed Hardy Robinson-effective precipitation."""

import numpy as np
import pytest

from icphysics import (
    collapse_hardy_effective_slice,
    collapse_hardy_robinson_effective,
    collapsed_hardy_robinson_effective_energy,
    hall,
    hardy_robinson_effective_precipitation,
    ped,
)
from icphysics.zhang_paxton import regular_latitude_grid


def test_collapsed_effective_precipitation_round_trips_robinson():
    result = collapse_hardy_robinson_effective(
        kp=np.array([0, 2.5, 6]),
        mlt=np.array([23.75, 1.5, 5]),
    )

    assert result["valid"].all()
    np.testing.assert_allclose(
        hall(result["representative_energy"], result["representative_effective_flux"]),
        result["collapsed_hall"],
        rtol=2e-15,
    )
    np.testing.assert_allclose(
        ped(result["representative_energy"], result["representative_effective_flux"]),
        result["collapsed_pedersen"],
        rtol=2e-15,
    )


def test_collapse_inverts_mean_conductance_not_mean_local_energy():
    edges = np.array([60.0, 61.0, 62.0, 63.0])
    local_energy = np.array([3.0, 3.0, 1000.0])
    local_pedersen = np.array([5.0, 5.0, 1e-9])
    local_hall = 0.45 * local_energy**0.85 * local_pedersen

    result = collapse_hardy_effective_slice(
        local_hall, local_pedersen, np.ones(3, dtype=bool), edges
    )

    direct_energy_mean = np.average(
        local_energy, weights=result["area_weights"]
    )
    expected_hall = np.average(local_hall, weights=result["area_weights"])
    expected_pedersen = np.average(
        local_pedersen, weights=result["area_weights"]
    )
    np.testing.assert_allclose(result["collapsed_hall"], expected_hall)
    np.testing.assert_allclose(
        result["collapsed_pedersen"], expected_pedersen
    )
    assert result["representative_energy"] < 3.01
    assert direct_energy_mean > 300
    assert not np.isclose(result["representative_energy"], direct_energy_mean)
    np.testing.assert_allclose(
        hall(result["representative_energy"], result["representative_effective_flux"]),
        result["collapsed_hall"],
    )


def test_disconnected_valid_region_away_from_peak_is_rejected():
    edges = np.arange(60.0, 66.0)
    pedersen = np.array([2.0, 10.0, 8.0, 9.0, 9.0])
    hall_conductance = pedersen.copy()
    local_valid = np.array([True, True, False, True, True])

    result = collapse_hardy_effective_slice(
        hall_conductance, pedersen, local_valid, edges
    )

    np.testing.assert_equal(result["selected"], [True, True, False, False, False])
    assert result["peak_mlat"] == 61.5
    assert result["selected_lower_mlat"] == 60
    assert result["selected_upper_mlat"] == 62


def test_native_backgrounds_are_excluded_from_selected_interval():
    centres, edges = regular_latitude_grid(50, 90, .05)
    local = hardy_robinson_effective_precipitation(0, 0, centres)
    result = collapse_hardy_effective_slice(
        local["hall"], local["pedersen"], local["valid"], edges
    )

    assert local["hall_background"].any()
    assert local["pedersen_background"].any()
    assert not local["hall_background"][result["selected"]].any()
    assert not local["pedersen_background"][result["selected"]].any()
    assert result["valid"]


def test_empty_and_sampling_boundary_diagnostics():
    empty = collapse_hardy_effective_slice(
        [1, 2], [1, 2], [False, False], [60, 61, 62]
    )
    assert empty["empty"]
    assert not empty["valid"]
    assert empty["selected_cell_count"] == 0
    assert np.isnan(empty["representative_energy"])

    bounded = collapse_hardy_robinson_effective(
        2, 1.5, lower_mlat=64, upper_mlat=66, latitude_step=.05
    )
    assert bounded["touches_equatorward_sampling_limit"]
    assert bounded["touches_poleward_sampling_limit"]
    assert not bounded["reaches_physical_pole"]


def test_broadcast_fractional_kp_masks_and_invalid_inputs():
    result = collapse_hardy_robinson_effective(
        np.array([0, 2.5])[:, None],
        np.array([0, 6, 12])[None, :],
        latitude_step=.05,
    )
    assert result["representative_energy"].shape == (2, 3)
    assert result["valid"].all()
    np.testing.assert_equal(result["kp"][:, 0], [0, 2.5])

    masked = collapse_hardy_robinson_effective(
        np.ma.array([2, 2], mask=[False, True]), 0, latitude_step=.05
    )
    assert np.isfinite(masked["representative_energy"][0])
    assert np.isnan(masked["representative_energy"][1])
    assert masked["empty"][1]

    invalid = collapse_hardy_robinson_effective(
        [np.nan, np.inf], 0, latitude_step=.05
    )
    assert np.isnan(invalid["representative_energy"]).all()
    assert invalid["empty"].all()


def test_cached_convenience_evaluator_is_periodic_and_preserves_nan():
    energy = collapsed_hardy_robinson_effective_energy(
        [2, 2, 2, np.nan],
        [0, 24, -24, 0],
        latitude_step=.05,
        mlt_step=6,
    )

    np.testing.assert_allclose(energy[:3], energy[0])
    assert np.isnan(energy[3])


@pytest.mark.parametrize("mlt_step", [0, -1, np.nan, 5])
def test_convenience_evaluator_rejects_invalid_mlt_grid(mlt_step):
    with pytest.raises(ValueError):
        collapsed_hardy_robinson_effective_energy(2, 0, mlt_step=mlt_step)


def test_collapse_rejects_invalid_latitude_grid_and_profile_shapes():
    with pytest.raises(ValueError):
        collapse_hardy_robinson_effective(2, 0, latitude_step=.03)
    with pytest.raises(ValueError):
        collapse_hardy_effective_slice(
            [[1, 2]], [1, 2], [True, True], [60, 61, 62]
        )
    with pytest.raises(ValueError):
        collapse_hardy_effective_slice(
            [1, 2], [1], [True, True], [60, 61, 62]
        )
