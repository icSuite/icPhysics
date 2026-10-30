"""Shared precipitation and ionospheric-conductance calculations."""

from .robinson import (
    hall, halluncertainty, ped, peduncertainty, robinson_conductance,
)
from .image import (
    PROTON_CORRECTION_UNCERTAINTY_MODES,
    PROTON_RESPONSE_ENERGY_RANGE,
    precipitation_from_ratio,
    precipitation_from_hardy,
    precipitation_from_hardy_effective,
    precipitation_from_collapsed_hardy_effective,
    ELECTRON_RESPONSE_ENERGY_RANGE,
    precipitation_from_zhang_paxton,
    proton_correct_images,
)
from .zhang_paxton import collapse_zhang_paxton
from .zhang_paxton_lookup import load_zhang_paxton_lookup
from .hardy import hardy_ion_precipitation
from .hardy_electron import (
    hardy_electron_precipitation, HARDY_ELECTRON_MODEL_DOI,
    HARDY_ELECTRON_COEFFICIENT_VERSION, HARDY_ELECTRON_KP_POLICY,
    HARDY_ELECTRON_LATITUDE_RANGE, HARDY_ELECTRON_UNCERTAINTY_METHOD
)
from .hardy_conductance import (
    hardy_conductance, hardy_robinson_effective_precipitation,
    HARDY_CONDUCTANCE_MODEL_DOI,
    HARDY_CONDUCTANCE_COEFFICIENT_VERSION, HARDY_CONDUCTANCE_KP_POLICY,
    HARDY_CONDUCTANCE_LATITUDE_RANGE, HARDY_EFFECTIVE_MINIMUM_ENERGY_KEV,
    HARDY_EFFECTIVE_METHOD
)
from .hardy_effective import (
    collapse_hardy_effective_slice,
    collapse_hardy_robinson_effective,
    collapsed_hardy_robinson_effective_energy,
    DEFAULT_HARDY_EFFECTIVE_MLAT_STEP_DEGREES,
    DEFAULT_HARDY_EFFECTIVE_MLT_STEP_HOURS,
)
from .spectrum import conductance_from_spectrum

__all__ = [
    "ped", "hall", "peduncertainty", "halluncertainty",
    "robinson_conductance",
    "precipitation_from_ratio", "precipitation_from_zhang_paxton",
    "proton_correct_images", "PROTON_RESPONSE_ENERGY_RANGE",
    "PROTON_CORRECTION_UNCERTAINTY_MODES",
    "collapse_zhang_paxton", "load_zhang_paxton_lookup",
    "hardy_ion_precipitation", "hardy_electron_precipitation",
    "precipitation_from_hardy", "precipitation_from_hardy_effective",
    "precipitation_from_collapsed_hardy_effective",
    "ELECTRON_RESPONSE_ENERGY_RANGE",
    "HARDY_ELECTRON_MODEL_DOI", "HARDY_ELECTRON_COEFFICIENT_VERSION",
    "HARDY_ELECTRON_KP_POLICY", "HARDY_ELECTRON_LATITUDE_RANGE",
    "HARDY_ELECTRON_UNCERTAINTY_METHOD",
    "hardy_conductance", "hardy_robinson_effective_precipitation",
    "HARDY_CONDUCTANCE_MODEL_DOI",
    "HARDY_CONDUCTANCE_COEFFICIENT_VERSION", "HARDY_CONDUCTANCE_KP_POLICY",
    "HARDY_CONDUCTANCE_LATITUDE_RANGE", "HARDY_EFFECTIVE_MINIMUM_ENERGY_KEV",
    "HARDY_EFFECTIVE_METHOD",
    "collapse_hardy_effective_slice", "collapse_hardy_robinson_effective",
    "collapsed_hardy_robinson_effective_energy",
    "DEFAULT_HARDY_EFFECTIVE_MLAT_STEP_DEGREES",
    "DEFAULT_HARDY_EFFECTIVE_MLT_STEP_HOURS",
    "conductance_from_spectrum",
]
