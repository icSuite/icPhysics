"""Reduce a Zhang--Paxton latitude profile to one energy per Kp and MLT.

The model gives electron mean energy ``E0`` and energy flux ``Q`` as functions
of Kp, magnetic local time (MLT), and magnetic latitude (MLAT).  For each MLT
slice we:

1. find the largest Q value;
2. keep the contiguous cells around that peak where Q is above a threshold;
3. calculate three complementary latitude collapses over that same interval:
   spherical-area mean, spherical-area times Q mean, and a Robinson-
   conductance-equivalent precipitation pair.

The default ``Q > 0.05 mW m-2`` threshold follows the criterion used for the
global-average energy in Zhang and Paxton (2008), Figure 8.  The 0.01-degree
latitude step resolves movement of the threshold boundary; it does not imply
that the empirical model has 0.01-degree physical resolution.

Functions return ordinary dictionaries so that the calculation and its
outputs remain easy to inspect in a research notebook or script.
"""

from functools import lru_cache

import numpy as np

from zhangpaxton2008 import zhang_paxton

from .robinson import hall, ped


DEFAULT_MLAT_STEP_DEGREES = 0.01
FIGURE8_THRESHOLD = 0.05
BOUNDARY_THRESHOLD = 0.25
MODEL_EVALUATION_BATCH_SIZE = 64
DEFAULT_COLLAPSED_MLT_STEP_HOURS = 0.1


def regular_latitude_grid(
    lower_mlat=50.0,
    upper_mlat=90.0,
    step=DEFAULT_MLAT_STEP_DEGREES,
):
    """Return latitude cell centres and edges in degrees."""

    if not np.isfinite([lower_mlat, upper_mlat, step]).all():
        raise ValueError("latitude limits and step must be finite")
    if not 0 <= lower_mlat < upper_mlat <= 90:
        raise ValueError("require 0 <= lower_mlat < upper_mlat <= 90 degrees")
    if step <= 0:
        raise ValueError("latitude step must be positive")

    width = upper_mlat - lower_mlat
    number_of_cells = int(round(width / step))
    if number_of_cells < 1 or not np.isclose(number_of_cells * step, width):
        raise ValueError("latitude step must divide the requested range exactly")

    edges = np.linspace(lower_mlat, upper_mlat, number_of_cells + 1)
    centres = (edges[:-1] + edges[1:]) / 2
    return centres, edges


def latitude_cell_area_weights(latitude_edges):
    """Return exact relative spherical areas for latitude cells.

    For a fixed MLT width, cell area is proportional to
    ``sin(latitude_upper) - sin(latitude_lower)``.  The common Earth-radius
    and MLT-width factors cancel from every weighted statistic.
    """

    edges = np.asarray(latitude_edges, dtype=float)
    if edges.ndim != 1 or edges.size < 2:
        raise ValueError("latitude_edges must be a one-dimensional edge array")
    if not np.all(np.isfinite(edges)) or not np.all(np.diff(edges) > 0):
        raise ValueError("latitude_edges must be finite and increasing")
    if edges[0] < -90 or edges[-1] > 90:
        raise ValueError("latitude edges must lie between -90 and 90 degrees")

    edges = np.deg2rad(edges)
    return np.sin(edges[1:]) - np.sin(edges[:-1])


def weighted_median(values, weights):
    """Return the value where cumulative positive weight reaches one half."""

    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if values.ndim != 1 or weights.ndim != 1 or values.shape != weights.shape:
        raise ValueError("values and weights must be one-dimensional and aligned")
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("values must be non-empty and finite")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("weights must be finite and non-negative")
    if weights.sum() <= 0:
        raise ValueError("at least one weight must be positive")

    order = np.argsort(values)
    cumulative_weight = np.cumsum(weights[order])
    median_index = np.searchsorted(cumulative_weight, weights.sum() / 2)
    return float(values[order][median_index])


def _contiguous_cells_around_peak(above_threshold, peak_index):
    """Select the connected above-threshold interval containing the Q peak."""

    selected = np.zeros(above_threshold.size, dtype=bool)
    if not above_threshold[peak_index]:
        return selected

    lower = peak_index
    while lower > 0 and above_threshold[lower - 1]:
        lower -= 1

    upper = peak_index
    while upper + 1 < above_threshold.size and above_threshold[upper + 1]:
        upper += 1

    selected[lower : upper + 1] = True
    return selected


def _empty_slice(area_weights, threshold, peak_mlat=np.nan, peak_flux=np.nan):
    """Return explicit missing values when no oval interval is selected."""

    return {
        "representative_energy": np.nan,
        "area_mean_energy": np.nan,
        "area_weighted_median_energy": np.nan,
        "weighted_spread": np.nan,
        "area_mean_energy_spread": np.nan,
        "q_weighted_energy": np.nan,
        "q_weighted_energy_spread": np.nan,
        "conductance_effective_energy": np.nan,
        "conductance_effective_energy_spread": np.nan,
        "conductance_effective_energy_flux": np.nan,
        "area_mean_hall_conductance": np.nan,
        "area_mean_pedersen_conductance": np.nan,
        "selected_lower_mlat": np.nan,
        "selected_upper_mlat": np.nan,
        "selected_area_weight": 0.0,
        "selected_area_fraction": 0.0,
        "threshold_cutoff": float(threshold),
        "peak_mlat": float(peak_mlat),
        "peak_energy_flux": float(peak_flux),
        "empty": True,
        "touches_equatorward_sampling_limit": False,
        "touches_poleward_sampling_limit": False,
        "reaches_physical_pole": False,
        "selected": np.zeros(area_weights.size, dtype=bool),
        "area_weights": area_weights,
    }


def collapse_latitude_slice(
    mean_energy,
    energy_flux,
    latitude_edges,
    threshold=FIGURE8_THRESHOLD,
):
    """Collapse one latitude profile inside its principal auroral interval.

    All methods use the same contiguous interval around the Q peak where
    ``Q > threshold``:

    - ``area_mean_energy`` uses spherical cell area and its structural spread
      is the area-weighted RMS deviation about that mean.
    - ``q_weighted_energy`` uses spherical area times Q and its structural
      spread is the matching weighted RMS deviation about that mean.
    - ``conductance_effective_energy`` area-averages local Robinson Hall and
      Pedersen conductances, then inverts that pair to effective E0 and Q.
      Its structural spread is the area-times-sqrt(Q)-weighted RMS deviation
      of local E0 about the effective E0. Robinson conductance amplitude is
      proportional to sqrt(Q), so this weights cells by their common ability
      to contribute conductance without preferring Hall or Pedersen.

    These spreads quantify MLAT structure discarded by each collapse. They
    are not coefficient, measurement, or predictive uncertainties. The old
    ``representative_energy`` and ``weighted_spread`` names remain aliases for
    the area-mean method.
    """

    energy = np.asarray(mean_energy, dtype=float)
    flux = np.asarray(energy_flux, dtype=float)
    edges = np.asarray(latitude_edges, dtype=float)
    area_weights = latitude_cell_area_weights(edges)

    if energy.ndim != 1 or flux.ndim != 1:
        raise ValueError("mean_energy and energy_flux must be one-dimensional")
    if energy.shape != flux.shape or energy.size != area_weights.size:
        raise ValueError("energy, flux, and latitude edges must describe the same cells")
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("threshold must be finite and non-negative")

    latitude_centres = (edges[:-1] + edges[1:]) / 2
    finite_flux = np.isfinite(flux)
    if not finite_flux.any():
        # Match the historical diagnostic: no finite Q means the cutoff was
        # never evaluated, so peak and cutoff diagnostics are all missing.
        return _empty_slice(area_weights, np.nan)

    # Non-finite values are unavailable and can never enter the selected oval.
    peak_index = int(np.argmax(np.where(finite_flux, flux, -np.inf)))
    peak_mlat = latitude_centres[peak_index]
    peak_flux = flux[peak_index]
    above_threshold = finite_flux & np.isfinite(energy) & (flux > threshold)
    selected = _contiguous_cells_around_peak(above_threshold, peak_index)
    if not selected.any():
        return _empty_slice(area_weights, threshold, peak_mlat, peak_flux)

    indices = np.flatnonzero(selected)
    lower_index = indices[0]
    upper_index = indices[-1]
    selected_energy = energy[selected]
    selected_area = area_weights[selected]
    selected_area_sum = selected_area.sum()

    mean = float(np.average(selected_energy, weights=selected_area))
    median = weighted_median(selected_energy, selected_area)
    variance = np.average((selected_energy - mean) ** 2, weights=selected_area)

    # This statistic emphasizes locations carrying most modeled energy flux.
    # It is not the primary spatial mean or a particle-number-weighted energy.
    power_weights = selected_area * np.clip(flux[selected], 0, None)
    if power_weights.sum() > 0:
        q_weighted_mean = float(np.average(selected_energy, weights=power_weights))
        q_weighted_variance = np.average(
            (selected_energy - q_weighted_mean) ** 2,
            weights=power_weights,
        )
        q_weighted_spread = float(np.sqrt(max(q_weighted_variance, 0)))
    else:
        q_weighted_mean = np.nan
        q_weighted_spread = np.nan

    selected_flux = flux[selected]
    local_hall = hall(selected_energy, selected_flux)
    local_pedersen = ped(selected_energy, selected_flux)
    mean_hall = float(np.average(local_hall, weights=selected_area))
    mean_pedersen = float(np.average(local_pedersen, weights=selected_area))

    conductance_effective_energy = (
        mean_hall / (0.45 * mean_pedersen)
    ) ** (1 / 0.85)
    conductance_effective_flux = (
        mean_pedersen
        * (16 + conductance_effective_energy**2)
        / (40 * conductance_effective_energy)
    ) ** 2

    conductance_weights = selected_area * np.sqrt(selected_flux)
    conductance_variance = np.average(
        (selected_energy - conductance_effective_energy) ** 2,
        weights=conductance_weights,
    )
    conductance_spread = float(np.sqrt(max(conductance_variance, 0)))

    return {
        "representative_energy": mean,
        "area_mean_energy": mean,
        "area_weighted_median_energy": median,
        "weighted_spread": float(np.sqrt(max(variance, 0))),
        "area_mean_energy_spread": float(np.sqrt(max(variance, 0))),
        "q_weighted_energy": q_weighted_mean,
        "q_weighted_energy_spread": q_weighted_spread,
        "conductance_effective_energy": float(conductance_effective_energy),
        "conductance_effective_energy_spread": conductance_spread,
        "conductance_effective_energy_flux": float(conductance_effective_flux),
        "area_mean_hall_conductance": mean_hall,
        "area_mean_pedersen_conductance": mean_pedersen,
        "selected_lower_mlat": float(edges[lower_index]),
        "selected_upper_mlat": float(edges[upper_index + 1]),
        "selected_area_weight": float(selected_area_sum),
        "selected_area_fraction": float(selected_area_sum / area_weights.sum()),
        "threshold_cutoff": float(threshold),
        "peak_mlat": float(peak_mlat),
        "peak_energy_flux": float(peak_flux),
        "empty": False,
        "touches_equatorward_sampling_limit": lower_index == 0,
        "touches_poleward_sampling_limit": (
            upper_index == energy.size - 1 and not np.isclose(edges[-1], 90)
        ),
        "reaches_physical_pole": (
            upper_index == energy.size - 1 and np.isclose(edges[-1], 90)
        ),
        "selected": selected,
        "area_weights": area_weights,
    }


def collapse_zhang_paxton(
    kp,
    mlt,
    threshold=FIGURE8_THRESHOLD,
    lower_mlat=50.0,
    upper_mlat=90.0,
    latitude_step=DEFAULT_MLAT_STEP_DEGREES,
):
    """Evaluate and collapse all broadcast pairs of Kp and MLT.

    Array fields in the returned dictionary have the broadcast shape of
    ``kp`` and ``mlt``.  Scalar input therefore gives zero-dimensional arrays.
    """

    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("threshold must be finite and non-negative")

    kp_grid, mlt_grid = np.broadcast_arrays(
        np.asarray(kp, dtype=float),
        np.asarray(mlt, dtype=float),
    )
    latitude_centres, latitude_edges = regular_latitude_grid(
        lower_mlat,
        upper_mlat,
        latitude_step,
    )

    float_fields = [
        "representative_energy",
        "area_mean_energy",
        "area_weighted_median_energy",
        "weighted_spread",
        "area_mean_energy_spread",
        "q_weighted_energy",
        "q_weighted_energy_spread",
        "conductance_effective_energy",
        "conductance_effective_energy_spread",
        "conductance_effective_energy_flux",
        "area_mean_hall_conductance",
        "area_mean_pedersen_conductance",
        "selected_lower_mlat",
        "selected_upper_mlat",
        "selected_area_weight",
        "selected_area_fraction",
        "threshold_cutoff",
        "peak_mlat",
        "peak_energy_flux",
    ]
    flag_fields = [
        "empty",
        "touches_equatorward_sampling_limit",
        "touches_poleward_sampling_limit",
        "reaches_physical_pole",
    ]
    result = {
        name: np.full(kp_grid.shape, np.nan, dtype=float)
        for name in float_fields
    }
    result.update({
        name: np.zeros(kp_grid.shape, dtype=bool)
        for name in flag_fields
    })

    flat_kp = kp_grid.ravel()
    flat_mlt = mlt_grid.ravel()
    for start in range(0, flat_kp.size, MODEL_EVALUATION_BATCH_SIZE):
        stop = min(start + MODEL_EVALUATION_BATCH_SIZE, flat_kp.size)
        energy, flux = zhang_paxton(
            flat_kp[start:stop, None],
            flat_mlt[start:stop, None],
            latitude_centres[None, :],
        )

        # Keep this slice loop explicit: it is the scientific reduction and
        # costs little compared with evaluating the empirical model.
        for row, flat_index in enumerate(range(start, stop)):
            collapsed = collapse_latitude_slice(
                energy[row],
                flux[row],
                latitude_edges,
                threshold,
            )
            output_index = np.unravel_index(flat_index, kp_grid.shape)
            for name in float_fields + flag_fields:
                result[name][output_index] = collapsed[name]

    result["kp"] = kp_grid.copy()
    result["mlt"] = np.mod(mlt_grid, 24)
    result["latitude_edges"] = latitude_edges
    result["threshold"] = float(threshold)
    return result


#%% Cached flux-weighted collapsed evaluator

@lru_cache(maxsize=None)
def _q_weighted_mlt_profile(
    kp,
    lower_mlat,
    upper_mlat,
    latitude_step,
    mlt_step,
):
    """Calculate one reusable periodic Q-weighted MLT profile."""

    mlt = np.arange(0, 24, mlt_step)
    collapsed = collapse_zhang_paxton(
        np.full(mlt.shape, kp),
        mlt,
        lower_mlat=lower_mlat,
        upper_mlat=upper_mlat,
        latitude_step=latitude_step,
    )
    return (
        mlt,
        collapsed["q_weighted_energy"],
        collapsed["q_weighted_energy_spread"],
    )


def collapsed_zhang_paxton_q_weighted(
    kp,
    mlt,
    lower_mlat=50.0,
    upper_mlat=90.0,
    latitude_step=DEFAULT_MLAT_STEP_DEGREES,
    mlt_step=DEFAULT_COLLAPSED_MLT_STEP_HOURS,
):
    """Evaluate flux-weighted collapsed Zhang--Paxton energy at any MLT.

    The returned ``energy_spread`` is the Q-weighted RMS latitude variation
    discarded by the collapse. It is structural spread, not predictive model
    uncertainty. Kp is rounded to the nearest hundredth to match the lookup
    tables and periodic MLT profiles are cached for efficient orbit use.
    """

    if not np.isfinite(mlt_step) or mlt_step <= 0 or mlt_step > 24:
        raise ValueError("mlt_step must be finite and between 0 and 24 hours")
    if not np.isclose(24 / mlt_step, round(24 / mlt_step)):
        raise ValueError("mlt_step must divide 24 hours exactly")

    regular_latitude_grid(lower_mlat, upper_mlat, latitude_step)
    kp_grid, mlt_grid = np.broadcast_arrays(
        np.asarray(np.ma.filled(np.ma.asarray(kp, dtype=float), np.nan)),
        np.asarray(np.ma.filled(np.ma.asarray(mlt, dtype=float), np.nan)),
    )
    energy = np.full(kp_grid.shape, np.nan, dtype=float)
    spread = np.full(kp_grid.shape, np.nan, dtype=float)
    finite = (
        np.isfinite(kp_grid) & (kp_grid >= 0) & (kp_grid <= 9)
        & np.isfinite(mlt_grid)
    )
    effective_kp = np.floor(kp_grid * 100 + 0.5 + 1e-12) / 100

    for kp_value in np.unique(effective_kp[finite]):
        profile_mlt, profile_energy, profile_spread = (
            _q_weighted_mlt_profile(
                float(kp_value),
                float(lower_mlat),
                float(upper_mlat),
                float(latitude_step),
                float(mlt_step),
            )
        )
        selected = finite & (effective_kp == kp_value)
        periodic_mlt = np.append(profile_mlt, 24.0)
        periodic_energy = np.append(profile_energy, profile_energy[0])
        periodic_spread = np.append(profile_spread, profile_spread[0])
        target_mlt = np.mod(mlt_grid[selected], 24)
        energy[selected] = np.interp(
            target_mlt, periodic_mlt, periodic_energy
        )
        spread[selected] = np.interp(
            target_mlt, periodic_mlt, periodic_spread
        )

    return {
        "energy": energy,
        "energy_spread": spread,
        "kp": np.where(finite, effective_kp, np.nan),
        "mlt": np.where(finite, np.mod(mlt_grid, 24), np.nan),
        "valid": np.isfinite(energy) & np.isfinite(spread),
    }
