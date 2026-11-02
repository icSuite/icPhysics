"""Read legacy and detector-grid Zhang--Paxton electron-energy lookups."""

#%% Imports and public constants

from hashlib import sha256
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


KP_VALUES = np.arange(901) / 100

LEGACY_GRID_ID = "image_wic_36x36_legacy_v1"
DETECTOR_GRID_ID = "image_apex_130km_46x46_v1"

DEFAULT_LOOKUP_PATH = Path(__file__).parent / "data" / "zhang_paxton_e0_lookup.nc"
DETECTOR_LOOKUP_PATH = (
    Path(__file__).parent / "data" / "zhang_paxton_e0_image_apex_130km_46x46_v1.nc"
)
DETECTOR_COORDINATE_TEMPLATE_PATH = (
    Path(__file__).parent / "data" / "image_apex_130km_46x46_v1_coordinates.nc"
)

# Filled after exporting the canonical coordinates from icBuilder. Keeping the
# expected value in code makes coordinated changes to a file and its own
# checksum insufficient to hide accidental grid drift.
DETECTOR_COORDINATE_CHECKSUM = (
    "96076278e1d9125b48555f09d8ca4ece1f93c4afd131fd697504121b34a4308f"
)

ENERGY_VARIANTS = (
    "full",
    "area_mean",
    "q_weighted",
    "conductance_effective",
    "collapsed",
)
EXPECTED_DETECTOR_SHAPE = (46, 46)

VARIANT_FIELDS = {
    "area_mean": ("E0_area_mean", "dE0_area_mean", "E0_area_median"),
    "q_weighted": ("E0_q_weighted", "dE0_q_weighted"),
    "conductance_effective": (
        "E0_conductance_effective",
        "dE0_conductance_effective",
        "Q_conductance_effective",
    ),
}


#%% Shared validation

def kp_to_index(kp):
    """Return the nearest-hundredth Kp index."""

    kp = np.asarray(kp, dtype=float)
    if np.any(~np.isfinite(kp)) or np.any((kp < 0) | (kp > 9)):
        raise ValueError("Kp must be finite and between 0 and 9")
    return np.floor(kp * 100 + 0.5 + 1e-12).astype(int)


def coordinate_checksum(grid_id, xi, eta, mlat, mlt):
    """Return a stable SHA-256 checksum for one fixed coordinate grid."""

    digest = sha256()
    digest.update(str(grid_id).encode("utf-8"))
    for name, values in (
        ("xi", xi), ("eta", eta), ("mlat", mlat), ("mlt", mlt)
    ):
        values = np.asarray(values, dtype="<f8")
        digest.update(name.encode("ascii"))
        digest.update(np.asarray(values.shape, dtype="<i8").tobytes())
        digest.update(np.ascontiguousarray(values).tobytes())
    return digest.hexdigest()


def _read_kp_coordinates(nc, expected_shape):
    stored_kp = np.asarray(nc["kp"][:])
    if stored_kp.shape != (901,) or not np.allclose(stored_kp, KP_VALUES):
        raise ValueError("lookup must contain Kp 0.00 to 9.00 in 0.01 steps")

    coordinates = {
        name: np.asarray(nc[name][:], dtype=float)
        for name in ("xi", "eta", "mlt")
    }
    if "mlat" in nc.variables:
        coordinates["mlat"] = np.asarray(nc["mlat"][:], dtype=float)

    for name, values in coordinates.items():
        if values.shape != expected_shape or np.any(~np.isfinite(values)):
            raise ValueError(
                f"lookup coordinate {name} must be finite with shape {expected_shape}"
            )
    return stored_kp, coordinates


def _read_selected_layers(nc, names, index, expected_shape):
    unique_index, inverse = np.unique(index.reshape(-1), return_inverse=True)
    output_shape = index.shape + expected_shape
    fields = {}
    for name in names:
        if name not in nc.variables:
            raise ValueError(f"lookup is missing required field {name}")
        values = np.asarray(nc[name][unique_index, :, :])[inverse]
        fields[name] = values.reshape(output_shape)
    return fields


def _read_area_mean_fields(nc, index, expected_shape, legacy_names):
    """Read explicit area fields, falling back to the frozen legacy aliases."""

    explicit_names = VARIANT_FIELDS["area_mean"]
    has_explicit = all(name in nc.variables for name in explicit_names)
    has_legacy = all(name in nc.variables for name in legacy_names)
    if not has_explicit and not has_legacy:
        raise ValueError("lookup is missing area-mean collapse fields")

    if has_explicit:
        fields = _read_selected_layers(nc, explicit_names, index, expected_shape)
        if has_legacy:
            aliases = _read_selected_layers(nc, legacy_names, index, expected_shape)
            for explicit_name, legacy_name in zip(explicit_names, legacy_names):
                if not np.array_equal(fields[explicit_name], aliases[legacy_name]):
                    raise ValueError(
                        f"lookup has ambiguous area-mean fields: {explicit_name} "
                        f"does not match {legacy_name}"
                    )
        return fields

    aliases = _read_selected_layers(nc, legacy_names, index, expected_shape)
    return {
        explicit_name: aliases[legacy_name]
        for explicit_name, legacy_name in zip(explicit_names, legacy_names)
    }


#%% Legacy 36-by-36 lookup

def _load_legacy_lookup(kp, path, energy_variant=None):
    index = kp_to_index(kp)
    with Dataset(path) as nc:
        _, coordinates = _read_kp_coordinates(nc, (36, 36))
        if energy_variant is None or energy_variant == "collapsed":
            fields = _read_selected_layers(
                nc, ("E0", "dE0", "E0_median"), index, (36, 36)
            )
        elif energy_variant == "area_mean":
            fields = _read_area_mean_fields(
                nc, index, (36, 36), ("E0", "dE0", "E0_median")
            )
        elif energy_variant in ("q_weighted", "conductance_effective"):
            fields = _read_selected_layers(
                nc, VARIANT_FIELDS[energy_variant], index, (36, 36)
            )
        else:
            raise ValueError(
                "legacy grid energy_variant must be area_mean, q_weighted, "
                "conductance_effective, or collapsed"
            )
        provenance = {
            "description": nc.description,
            "threshold_energy_flux_mW_m2": float(nc.threshold_energy_flux_mW_m2),
            "lower_mlat_degrees": float(nc.lower_mlat_degrees),
            "upper_mlat_degrees": float(nc.upper_mlat_degrees),
            "latitude_step_degrees": float(nc.latitude_step_degrees),
            "zhang_paxton_package_version": nc.zhang_paxton_package_version,
        }

    for name, values in fields.items():
        if np.any(~np.isfinite(values)):
            raise ValueError(f"legacy lookup field {name} contains non-finite values")

    result = {"kp": index / 100, **fields, **coordinates, "provenance": provenance}
    if energy_variant is not None:
        result["grid_id"] = LEGACY_GRID_ID
        result["energy_variant"] = energy_variant
    return result


#%% Detector 46-by-46 lookup

def _validate_detector_provenance(nc, coordinates):
    if getattr(nc, "grid_id", None) != DETECTOR_GRID_ID:
        raise ValueError(f"lookup grid_id must be {DETECTOR_GRID_ID!r}")

    stored_checksum = getattr(nc, "coordinate_template_sha256", None)
    calculated_checksum = coordinate_checksum(
        DETECTOR_GRID_ID,
        coordinates["xi"], coordinates["eta"],
        coordinates["mlat"], coordinates["mlt"],
    )
    if stored_checksum != calculated_checksum:
        raise ValueError("lookup coordinate checksum does not match its coordinates")
    if stored_checksum != DETECTOR_COORDINATE_CHECKSUM:
        raise ValueError("lookup coordinates do not match the canonical detector grid")

    expected_attributes = {
        "lower_mlat_degrees": 50.0,
        "upper_mlat_degrees": 90.0,
        "threshold_energy_flux_mW_m2": 0.05,
        "latitude_step_degrees": 0.01,
        "kp_selection_policy": "nearest_hundredth_no_interpolation",
    }
    for name, expected in expected_attributes.items():
        actual = getattr(nc, name, None)
        if isinstance(expected, float):
            if actual is None or not np.isclose(actual, expected):
                raise ValueError(f"lookup provenance attribute {name} is invalid")
        elif actual != expected:
            raise ValueError(f"lookup provenance attribute {name} is invalid")


def _load_detector_lookup(kp, path, energy_variant):
    if energy_variant not in ENERGY_VARIANTS:
        raise ValueError(f"energy_variant must be one of {ENERGY_VARIANTS}")

    index = kp_to_index(kp)
    with Dataset(path) as nc:
        _, coordinates = _read_kp_coordinates(nc, EXPECTED_DETECTOR_SHAPE)
        if "mlat" not in coordinates:
            raise ValueError("detector lookup must contain mlat coordinates")
        _validate_detector_provenance(nc, coordinates)

        if energy_variant == "full":
            names = ("E0_full", "Q_full", "full_valid")
            fields = _read_selected_layers(
                nc, names, index, EXPECTED_DETECTOR_SHAPE
            )
        elif energy_variant == "collapsed":
            names = (
                "E0_collapsed", "dE0_collapsed", "E0_collapsed_median",
                "collapsed_touches_equatorward_limit",
            )
            fields = _read_selected_layers(
                nc, names, index, EXPECTED_DETECTOR_SHAPE
            )
        elif energy_variant == "area_mean":
            fields = _read_area_mean_fields(
                nc, index, EXPECTED_DETECTOR_SHAPE,
                ("E0_collapsed", "dE0_collapsed", "E0_collapsed_median"),
            )
            fields.update(_read_selected_layers(
                nc, ("collapsed_touches_equatorward_limit",),
                index, EXPECTED_DETECTOR_SHAPE,
            ))
        else:
            fields = _read_selected_layers(
                nc,
                VARIANT_FIELDS[energy_variant]
                + ("collapsed_touches_equatorward_limit",),
                index,
                EXPECTED_DETECTOR_SHAPE,
            )
        provenance = {name: getattr(nc, name) for name in nc.ncattrs()}

    if energy_variant == "full":
        fields["full_valid"] = fields["full_valid"].astype(bool)
        expected_valid = (
            (coordinates["mlat"] >= 50) & (coordinates["mlat"] <= 90)
        )
        expected_valid = np.broadcast_to(
            expected_valid, index.shape + EXPECTED_DETECTOR_SHAPE
        )
        if not np.array_equal(fields["full_valid"], expected_valid):
            raise ValueError("full_valid does not match the 50--90 degree domain")
        for name in ("E0_full", "Q_full"):
            finite = np.isfinite(fields[name])
            if not np.array_equal(finite, fields["full_valid"]):
                raise ValueError(
                    f"{name} must be finite exactly where full_valid is true"
                )
    else:
        fields["collapsed_touches_equatorward_limit"] = fields[
            "collapsed_touches_equatorward_limit"
        ].astype(bool)
        for name, values in fields.items():
            if name == "collapsed_touches_equatorward_limit":
                continue
            if np.any(~np.isfinite(values)):
                raise ValueError(f"detector lookup field {name} contains non-finite values")

    return {
        "kp": index / 100,
        "grid_id": DETECTOR_GRID_ID,
        "energy_variant": energy_variant,
        **fields,
        **coordinates,
        "provenance": provenance,
    }


#%% Public loader

def load_zhang_paxton_lookup(
    kp, path=None, *, grid_id=None, energy_variant=None
):
    """Load one explicitly selected Zhang--Paxton energy lookup.

    Calling with only ``kp`` preserves the legacy 36-by-36 collapsed lookup
    API. The detector-grid lookup requires explicit ``grid_id`` and
    ``energy_variant`` selections so full and collapsed energy cannot be
    confused. Kp is rounded to the nearest hundredth without interpolation.
    """

    legacy_call = grid_id is None and energy_variant is None
    if legacy_call:
        selected_path = DEFAULT_LOOKUP_PATH if path is None else Path(path)
        return _load_legacy_lookup(kp, selected_path)

    if grid_id is None or energy_variant is None:
        raise ValueError("grid_id and energy_variant must be selected together")

    if grid_id == LEGACY_GRID_ID:
        if energy_variant == "full":
            raise ValueError("the legacy grid does not contain full spatial energy")
        selected_path = DEFAULT_LOOKUP_PATH if path is None else Path(path)
        return _load_legacy_lookup(kp, selected_path, energy_variant)

    if grid_id == DETECTOR_GRID_ID:
        selected_path = DETECTOR_LOOKUP_PATH if path is None else Path(path)
        return _load_detector_lookup(kp, selected_path, energy_variant)

    raise ValueError(f"unsupported Zhang--Paxton lookup grid_id {grid_id!r}")
