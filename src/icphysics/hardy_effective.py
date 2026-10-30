"""Collapse native Hardy conductance to one Robinson-effective energy.

The latitude reduction deliberately averages Hall and Pedersen conductance
before applying the inverse Robinson equations. It never averages local
effective energies. This makes strong-conductance cells control the collapsed
Hall/Pedersen ratio and prevents fitted-boundary energy spikes from dominating.
"""

#%% Imports and public defaults

from functools import lru_cache

import numpy as np

from .hardy_conductance import (
    HARDY_EFFECTIVE_MINIMUM_ENERGY_KEV,
    _inverse_robinson_conductance,
    hardy_robinson_effective_precipitation,
)
from .zhang_paxton import (
    DEFAULT_MLAT_STEP_DEGREES,
    latitude_cell_area_weights,
    regular_latitude_grid,
)


DEFAULT_HARDY_EFFECTIVE_MLAT_STEP_DEGREES = DEFAULT_MLAT_STEP_DEGREES
DEFAULT_HARDY_EFFECTIVE_MLT_STEP_HOURS = 0.05
MODEL_EVALUATION_BATCH_SIZE = 64


#%% Collapse one latitude slice

def _contiguous_interval_around_peak(valid, peak_index):
    """Return the connected valid interval containing the selected peak."""

    selected = np.zeros(valid.size, dtype=bool)
    if not valid[peak_index]:
        return selected

    lower = peak_index
    while lower > 0 and valid[lower - 1]:
        lower -= 1

    upper = peak_index
    while upper + 1 < valid.size and valid[upper + 1]:
        upper += 1

    selected[lower : upper + 1] = True
    return selected


def _empty_slice(area_weights, peak_mlat=np.nan, peak_pedersen=np.nan):
    """Return explicit missing diagnostics for an empty latitude slice."""

    return {
        "representative_energy": np.nan,
        "representative_effective_flux": np.nan,
        "collapsed_hall": np.nan,
        "collapsed_pedersen": np.nan,
        "selected_lower_mlat": np.nan,
        "selected_upper_mlat": np.nan,
        "peak_mlat": float(peak_mlat),
        "peak_pedersen": float(peak_pedersen),
        "selected_area_weight": 0.0,
        "selected_area_fraction": 0.0,
        "selected_cell_count": 0,
        "valid": False,
        "empty": True,
        "touches_equatorward_sampling_limit": False,
        "touches_poleward_sampling_limit": False,
        "reaches_physical_pole": False,
        "selected": np.zeros(area_weights.size, dtype=bool),
        "area_weights": area_weights,
    }


def collapse_hardy_effective_slice(
    hall_conductance,
    pedersen_conductance,
    local_valid,
    latitude_edges,
):
    """Collapse one background-screened Hardy latitude profile.

    Hall and Pedersen conductance are spherical-area averaged over the
    contiguous locally valid interval containing the Pedersen peak. The
    representative effective precipitation is then obtained by inverse
    Robinson. ``local_valid`` is expected to come from
    :func:`hardy_robinson_effective_precipitation`.
    """

    hall_conductance = np.asarray(hall_conductance, dtype=float)
    pedersen_conductance = np.asarray(pedersen_conductance, dtype=float)
    local_valid = np.asarray(local_valid, dtype=bool)
    latitude_edges = np.asarray(latitude_edges, dtype=float)
    area_weights = latitude_cell_area_weights(latitude_edges)

    if (
        hall_conductance.ndim != 1
        or pedersen_conductance.ndim != 1
        or local_valid.ndim != 1
    ):
        raise ValueError("conductance and validity profiles must be one-dimensional")
    if not (
        hall_conductance.shape
        == pedersen_conductance.shape
        == local_valid.shape
        == area_weights.shape
    ):
        raise ValueError("profiles and latitude edges must describe the same cells")

    local_valid = (
        local_valid
        & np.isfinite(hall_conductance)
        & np.isfinite(pedersen_conductance)
        & (hall_conductance > 0)
        & (pedersen_conductance > 0)
    )
    if not local_valid.any():
        return _empty_slice(area_weights)

    latitude_centres = (latitude_edges[:-1] + latitude_edges[1:]) / 2
    peak_index = int(np.argmax(np.where(local_valid, pedersen_conductance, -np.inf)))
    peak_mlat = latitude_centres[peak_index]
    peak_pedersen = pedersen_conductance[peak_index]
    selected = _contiguous_interval_around_peak(local_valid, peak_index)
    if not selected.any():
        return _empty_slice(area_weights, peak_mlat, peak_pedersen)

    indices = np.flatnonzero(selected)
    lower_index = indices[0]
    upper_index = indices[-1]
    selected_area = area_weights[selected]
    selected_area_sum = selected_area.sum()
    collapsed_hall = float(np.average(hall_conductance[selected], weights=selected_area))
    collapsed_pedersen = float(np.average(pedersen_conductance[selected], weights=selected_area))
    energy, flux, invertible = _inverse_robinson_conductance(
        collapsed_hall, collapsed_pedersen
    )
    energy = float(energy)
    flux = float(flux)
    collapsed_valid = bool(
        invertible
        and np.isfinite(energy)
        and energy >= HARDY_EFFECTIVE_MINIMUM_ENERGY_KEV
    )

    return {
        "representative_energy": energy,
        "representative_effective_flux": flux,
        "collapsed_hall": collapsed_hall,
        "collapsed_pedersen": collapsed_pedersen,
        "selected_lower_mlat": float(latitude_edges[lower_index]),
        "selected_upper_mlat": float(latitude_edges[upper_index + 1]),
        "peak_mlat": float(peak_mlat),
        "peak_pedersen": float(peak_pedersen),
        "selected_area_weight": float(selected_area_sum),
        "selected_area_fraction": float(selected_area_sum / area_weights.sum()),
        "selected_cell_count": int(selected.sum()),
        "valid": collapsed_valid,
        "empty": False,
        "touches_equatorward_sampling_limit": lower_index == 0,
        "touches_poleward_sampling_limit": (
            upper_index == hall_conductance.size - 1
            and not np.isclose(latitude_edges[-1], 90)
        ),
        "reaches_physical_pole": (
            upper_index == hall_conductance.size - 1
            and np.isclose(latitude_edges[-1], 90)
        ),
        "selected": selected,
        "area_weights": area_weights,
    }


#%% Evaluate and collapse the spatial model

FLOAT_FIELDS = [
    "representative_energy",
    "representative_effective_flux",
    "collapsed_hall",
    "collapsed_pedersen",
    "selected_lower_mlat",
    "selected_upper_mlat",
    "peak_mlat",
    "peak_pedersen",
    "selected_area_weight",
    "selected_area_fraction",
]

INTEGER_FIELDS = ["selected_cell_count"]

FLAG_FIELDS = [
    "valid",
    "empty",
    "touches_equatorward_sampling_limit",
    "touches_poleward_sampling_limit",
    "reaches_physical_pole",
]


def _as_unmasked_float(value):
    return np.asarray(np.ma.filled(np.ma.asarray(value, dtype=float), np.nan))


def collapse_hardy_robinson_effective(
    kp,
    mlt,
    lower_mlat=50.0,
    upper_mlat=90.0,
    latitude_step=DEFAULT_HARDY_EFFECTIVE_MLAT_STEP_DEGREES,
):
    """Return one conductance-consistent Hardy energy per Kp/MLT pair.

    The model is evaluated on a regular MLAT grid. Among locally valid cells,
    the connected interval containing the Pedersen-conductance peak is kept.
    Hall and Pedersen are spherical-area averaged over that interval before
    inverse Robinson is applied. No upper-energy or flux threshold is used.
    """

    kp_grid, mlt_grid = np.broadcast_arrays(
        _as_unmasked_float(kp),
        _as_unmasked_float(mlt),
    )
    latitude_centres, latitude_edges = regular_latitude_grid(
        lower_mlat, upper_mlat, latitude_step
    )

    result = {
        name: np.full(kp_grid.shape, np.nan, dtype=float)
        for name in FLOAT_FIELDS
    }
    result.update({
        name: np.zeros(kp_grid.shape, dtype=int)
        for name in INTEGER_FIELDS
    })
    result.update({
        name: np.zeros(kp_grid.shape, dtype=bool)
        for name in FLAG_FIELDS
    })

    flat_kp = kp_grid.ravel()
    flat_mlt = mlt_grid.ravel()
    for start in range(0, flat_kp.size, MODEL_EVALUATION_BATCH_SIZE):
        stop = min(start + MODEL_EVALUATION_BATCH_SIZE, flat_kp.size)
        local = hardy_robinson_effective_precipitation(
            flat_kp[start:stop, None],
            flat_mlt[start:stop, None],
            latitude_centres[None, :],
        )

        # Keep the scientific reduction explicit and independently testable.
        for row, flat_index in enumerate(range(start, stop)):
            collapsed = collapse_hardy_effective_slice(
                local["hall"][row],
                local["pedersen"][row],
                local["valid"][row],
                latitude_edges,
            )
            output_index = np.unravel_index(flat_index, kp_grid.shape)
            for name in FLOAT_FIELDS + INTEGER_FIELDS + FLAG_FIELDS:
                result[name][output_index] = collapsed[name]

    finite_kp = np.isfinite(kp_grid)
    finite_mlt = np.isfinite(mlt_grid)
    result["input_kp"] = kp_grid.copy()
    result["kp"] = np.where(finite_kp, np.clip(kp_grid, 0, 6), np.nan)
    result["mlt"] = np.where(finite_mlt, np.mod(mlt_grid, 24), np.nan)
    result["latitude_edges"] = latitude_edges
    result["latitude_step"] = float(latitude_step)
    return result


#%% Cached periodic convenience evaluator

@lru_cache(maxsize=None)
def _collapsed_effective_energy_mlt_profile(
    kp,
    lower_mlat,
    upper_mlat,
    latitude_step,
    mlt_step,
):
    """Calculate one reusable periodic effective-energy MLT profile."""

    mlt = np.arange(0, 24, mlt_step)
    collapsed = collapse_hardy_robinson_effective(
        np.full(mlt.shape, kp),
        mlt,
        lower_mlat=lower_mlat,
        upper_mlat=upper_mlat,
        latitude_step=latitude_step,
    )
    return mlt, collapsed["representative_energy"]


def collapsed_hardy_robinson_effective_energy(
    kp,
    mlt,
    lower_mlat=50.0,
    upper_mlat=90.0,
    latitude_step=DEFAULT_HARDY_EFFECTIVE_MLAT_STEP_DEGREES,
    mlt_step=DEFAULT_HARDY_EFFECTIVE_MLT_STEP_HOURS,
):
    """Evaluate collapsed Hardy effective energy efficiently at arbitrary MLT.

    A regular MLT profile is cached for each effective Kp and linearly
    interpolated with periodic continuity at 0/24 hours. NaN inputs and NaN
    portions of a cached profile remain unavailable.
    """

    if not np.isfinite(mlt_step) or mlt_step <= 0 or mlt_step > 24:
        raise ValueError("mlt_step must be finite and between 0 and 24 hours")
    if not np.isclose(24 / mlt_step, round(24 / mlt_step)):
        raise ValueError("mlt_step must divide 24 hours exactly")

    # Validate the latitude grid even if every requested input is unavailable.
    regular_latitude_grid(lower_mlat, upper_mlat, latitude_step)
    kp_grid, mlt_grid = np.broadcast_arrays(
        _as_unmasked_float(kp),
        _as_unmasked_float(mlt),
    )
    energy = np.full(kp_grid.shape, np.nan, dtype=float)
    finite = np.isfinite(kp_grid) & np.isfinite(mlt_grid)
    effective_kp = np.clip(kp_grid, 0, 6)

    for kp_value in np.unique(effective_kp[finite]):
        profile_mlt, profile_energy = _collapsed_effective_energy_mlt_profile(
            float(kp_value),
            float(lower_mlat),
            float(upper_mlat),
            float(latitude_step),
            float(mlt_step),
        )
        selected = finite & (effective_kp == kp_value)
        periodic_mlt = np.append(profile_mlt, 24.0)
        periodic_energy = np.append(profile_energy, profile_energy[0])
        energy[selected] = np.interp(
            np.mod(mlt_grid[selected], 24),
            periodic_mlt,
            periodic_energy,
        )
    return energy
