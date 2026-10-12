"""Electron-spectrum conductance using the supplied legacy Rees calculation."""

import numpy as np


def _log_interpolate(altitude, profile_alt, values):
    """Interpolate final-axis profiles, including zero species-density limits."""
    lower = np.clip(np.searchsorted(profile_alt, altitude, side="right") - 1, 0, profile_alt.size - 2)
    upper = lower + 1
    fraction = (altitude - profile_alt[lower]) / (profile_alt[upper] - profile_alt[lower])
    left = values[..., lower]
    right = values[..., upper]
    # Temporary ones only permit log evaluation; zero endpoints explicitly
    # yield zero in the open interval, without introducing a density floor.
    log_left = np.log(np.where(left > 0, left, 1.0))
    log_right = np.log(np.where(right > 0, right, 1.0))
    interpolated = np.exp(log_left + fraction * (log_right - log_left))
    interpolated = np.where((left == 0) | (right == 0), 0.0, interpolated)
    interpolated = np.where(fraction == 0, left, interpolated)
    return np.where(fraction == 1, right, interpolated)


def conductance_from_spectrum(energy_keV, differential_energy_flux, magnetic_field_nt, atmosphere, *, altitude_samples=500, differential_energy_flux_std=None):
    """Infer Pedersen ``P`` and Hall ``H`` conductance in siemens.

    ``energy_keV`` is a positive, strictly monotonic 1D channel-center array
    with at least two channels. ``differential_energy_flux`` has final axis
    matching those channels, in erg cm^-2 s^-1 eV^-1 sr^-1. Leading spectrum
    dimensions, atmosphere-profile dimensions, and ``magnetic_field_nt``
    (positive field magnitude, in nT) broadcast together. Descending channel
    order is accepted. Units are never
    inferred. Any missing, nonfinite, or negative channel invalidates that
    complete spectrum; invalid field magnitudes also produce NaN results.

    ``atmosphere`` supplies nonnegative species profiles on ``altitude_km``:
    ``oxygen_cm3`` (atomic O), ``nitrogen_cm3`` (molecular N2),
    ``molecular_oxygen_cm3`` (O2), ``mass_density_g_cm3``, and
    ``temperature_k``. Mass density and temperature must be strictly positive.
    Altitudes must be a shared 1D increasing array spanning 70 through 295 km.
    Profiles have shape ``(..., n_altitudes)``; 1D profiles describe a shared
    atmosphere, while leading dimensions permit a different atmosphere per
    spectrum or ensemble member. Invalid or masked atmosphere values raise
    ValueError rather than being replaced. No atmosphere or field model is
    evaluated here.

    Optional ``differential_energy_flux_std`` supplies absolute one-standard-
    deviation errors in the same flux units, with a matching channel axis and
    broadcastable leading dimensions. Errors are propagated analytically to
    ``P_std``, ``H_std`` (S) and ``PH_covariance`` (S squared), with an
    ``uncertainty_valid`` flag and int16 ``uncertainty_channel_count`` counting
    finite nonnegative supplied errors (including zero). Energy-channel errors
    are assumed independent;
    atmosphere, B and channel energies are held fixed. Any uncertainty included
    in the supplied channel errors is propagated, but common calibration or
    physical-model errors and inter-channel covariance are not added.
    P/H retain their original values even when errors are invalid. NaN and
    masked errors are omitted from variance and covariance sums; uncertainty
    therefore describes only channels with reported errors, without imputing
    an error for missing channels. All-missing errors remain unavailable.
    Any supplied negative or infinite error invalidates the entire uncertainty
    result. A zero-ionization height affected by a nonzero reported error has
    a singular square-root derivative and also invalidates uncertainty;
    reported errors that are all exactly zero give zero uncertainty.

    This ports energy_to_conductance.idl's scientific equations: log-spaced
    70..295 km grid excluding its upper endpoint, log-interpolated atmosphere,
    midpoint bin widths with extrapolated endpoints, a pi solid-angle factor,
    the unidirectional deposition polynomial and revised Rees range,
    altitude-only recombination, and 30-amu ions. Channel ionization is summed
    before taking the electron-density square root. Conductivity integration
    omits both grid endpoints, as in the source. The original combination of
    isotropic angular integration and unidirectional deposition is preserved;
    it is an assumption requiring scientific review. Arithmetic uses float64
    instead of the source's float32 ionization array. Missing channels remain
    invalid rather than being silently replaced by zero as in the IDL source.
    Zero species densities (such as MSIS atomic O below 72.5 km) are supported
    by the geometric-interpolation limit: an interval touching zero is zero
    in its interior, while exact source endpoints retain their values. This
    explicitly resolves the source's undefined log(0) endpoint arithmetic.
    """
    # Input shapes, units, and valid support.
    energy = np.ma.asarray(energy_keV, dtype=np.float64).filled(np.nan)
    flux = np.ma.asarray(differential_energy_flux, dtype=np.float64).filled(np.nan)
    flux_std = None if differential_energy_flux_std is None else np.ma.asarray(differential_energy_flux_std, dtype=np.float64).filled(np.nan)
    field_nt = np.ma.asarray(magnetic_field_nt, dtype=np.float64).filled(np.nan)
    if energy.ndim != 1 or energy.size < 2 or not np.all(np.isfinite(energy) & (energy > 0)):
        raise ValueError("energy_keV must be a positive finite 1D array with at least two channels")
    if flux.ndim < 1 or flux.shape[-1] != energy.size:
        raise ValueError("differential_energy_flux must have a final axis matching energy_keV")
    if flux_std is not None and (flux_std.ndim < 1 or flux_std.shape[-1] != energy.size):
        raise ValueError("differential_energy_flux_std must have a final axis matching energy_keV")
    if np.all(np.diff(energy) < 0):
        energy = energy[::-1]
        flux = flux[..., ::-1]
        if flux_std is not None:
            flux_std = flux_std[..., ::-1]
    elif not np.all(np.diff(energy) > 0):
        raise ValueError("energy_keV must be strictly increasing or decreasing")
    if isinstance(altitude_samples, (bool, np.bool_)) or not isinstance(altitude_samples, (int, np.integer)) or altitude_samples < 3:
        raise ValueError("altitude_samples must be an integer of at least three")

    # Interpolate neutral atmosphere exponentially onto the legacy grid.
    profile_alt = np.ma.asarray(atmosphere["altitude_km"], dtype=np.float64).filled(np.nan)
    if profile_alt.ndim != 1 or profile_alt.size < 2 or not np.all(np.isfinite(profile_alt)) or not np.all(np.diff(profile_alt) > 0):
        raise ValueError("atmosphere altitude_km must be finite and strictly increasing")
    if profile_alt[0] > 70 or profile_alt[-1] < 295:
        raise ValueError("atmosphere altitude_km must span 70 through 295 km")
    altitude = 10 ** (np.log10(70.0) + np.arange(altitude_samples) * (np.log10(295.0) - np.log10(70.0)) / altitude_samples)
    source_profiles = {}
    for name in ("oxygen_cm3", "nitrogen_cm3", "molecular_oxygen_cm3", "mass_density_g_cm3", "temperature_k"):
        values = np.ma.asarray(atmosphere[name], dtype=np.float64).filled(np.nan)
        species = name in ("oxygen_cm3", "nitrogen_cm3", "molecular_oxygen_cm3")
        valid_values = (values >= 0) if species else (values > 0)
        if values.ndim < 1 or values.shape[-1] != profile_alt.size or not np.all(np.isfinite(values) & valid_values):
            bound = "nonnegative" if species else "positive"
            raise ValueError(f"atmosphere {name} must be {bound} finite with a final axis matching altitude_km")
        source_profiles[name] = values
    try:
        atmosphere_shape = np.broadcast_shapes(*(values.shape[:-1] for values in source_profiles.values()))
        batch_shape = np.broadcast_shapes(flux.shape[:-1], field_nt.shape, atmosphere_shape)
    except ValueError as error:
        raise ValueError("leading spectrum, magnetic-field, and atmosphere-profile dimensions must broadcast together") from error
    profiles = {}
    for name, values in source_profiles.items():
        values = np.broadcast_to(values, atmosphere_shape + (profile_alt.size,))
        profiles[name] = _log_interpolate(altitude, profile_alt, values)

    flux = np.broadcast_to(flux, batch_shape + (energy.size,))
    if flux_std is not None:
        try:
            flux_std = np.broadcast_to(flux_std, flux.shape)
        except ValueError as error:
            raise ValueError("differential_energy_flux_std must broadcast to the spectrum/atmosphere/field batch") from error
    field_nt = np.broadcast_to(field_nt, batch_shape)
    valid = np.all(np.isfinite(flux) & (flux >= 0), axis=-1) & np.isfinite(field_nt) & (field_nt > 0)
    safe_flux = np.where(valid[..., None], flux, 0.0)
    field = np.where(valid, field_nt, 50000.0) * 1e-9
    oxygen = profiles["oxygen_cm3"]
    nitrogen = profiles["nitrogen_cm3"]
    molecular_oxygen = profiles["molecular_oxygen_cm3"]
    mass_density = profiles["mass_density_g_cm3"]
    temperature = profiles["temperature_k"]
    pressure = temperature * 82.05 * (oxygen + molecular_oxygen + nitrogen) / 6.022e23
    depth = 10130.0 * pressure / 9.82

    # Convert differential energy flux to per-channel electron number flux.
    energy_steps = np.diff(energy)
    bin_widths = np.empty_like(energy)
    bin_widths[0] = energy_steps[0]
    bin_widths[-1] = energy_steps[-1]
    bin_widths[1:-1] = (energy_steps[:-1] + energy_steps[1:]) / 2
    energy_flux = safe_flux * np.pi * bin_widths * 1000.0
    number_flux = energy_flux * 6.242e11 / (energy * 1000.0)

    # Sum monoenergetic ionization profiles before calculating electron density.
    electron_range = 4.3e-7 + 5.36e-6 * energy ** 1.67
    relative_depth = depth[..., None, :] / electron_range[:, None]
    depositing = relative_depth <= 1
    x = np.where(depositing, relative_depth, 0.0)
    deposition = 0.751031 + 2.526201 * x + 0.276534 * x**2 - 8.344593 * x**3 + 4.796477 * x**4
    ionization_kernel = np.where(depositing, energy[:, None] * deposition * mass_density[..., None, :] / (0.035 * electron_range[:, None]), 0.0)
    if not atmosphere_shape:
        # A shared profile has one fixed kernel for all spectra.
        total_ionization = number_flux @ ionization_kernel
    else:
        # Each atmosphere has its own energy-by-altitude deposition kernel.
        total_ionization = np.sum(number_flux[..., :, None] * ionization_kernel, axis=-2)
    recombination = 4.3e-6 * np.exp(-2.42e-2 * altitude) + 8.16e12 * np.exp(-0.524 * altitude)
    electron_density = np.sqrt(total_ionization / recombination)

    # Electron and ion collision and gyro frequencies.
    thermal_energy = 3 * 8.63e-5 * temperature
    sigma_o2 = (2.2 + 5.1 * np.sqrt(thermal_energy)) * 1e-20
    sigma_n2 = (18.3 - 7.3 * thermal_energy) * np.sqrt(thermal_energy) * 1e-20
    electron_collisions = 1e6 * (molecular_oxygen * sigma_o2 + nitrogen * sigma_n2 + oxygen * 1e-19) * np.sqrt(8 * 1.381e-23 * temperature / (np.pi * 9.11e-31))
    ion_collisions = 3.75e-10 * (molecular_oxygen + nitrogen + oxygen)
    electron_gyro = 1.602e-19 * field[..., None] / 9.11e-31
    ion_gyro = 1.602e-19 * field[..., None] / (30 * 1.672e-27)

    # Conductivity and the source's interior-point height integration.
    density_factor = electron_density * 1e6 * 1.602e-19 / field[..., None]
    pedersen = density_factor * (electron_collisions * electron_gyro / (electron_gyro**2 + electron_collisions**2) + ion_collisions * ion_gyro / (ion_gyro**2 + ion_collisions**2))
    hall = density_factor * (electron_gyro**2 / (electron_gyro**2 + electron_collisions**2) - ion_gyro**2 / (ion_gyro**2 + ion_collisions**2))
    height_weights = (altitude[2:] - altitude[:-2]) * 1000.0 / 2
    p = np.sum(pedersen[..., 1:-1] * height_weights, axis=-1)
    h = np.sum(hall[..., 1:-1] * height_weights, axis=-1)
    result = {"P": np.where(valid, p, np.nan), "H": np.where(valid, h, np.nan)}
    if flux_std is None:
        return result

    # Derivatives reuse the ionization and conductivity already calculated.
    reported = np.isfinite(flux_std) & (flux_std >= 0)
    channel_count = np.sum(reported, axis=-1, dtype=np.int16)
    invalid_error = np.any((flux_std < 0) | np.isinf(flux_std), axis=-1)
    uncertainty_valid = valid & (channel_count > 0) & ~invalid_error
    safe_std = np.where(uncertainty_valid[..., None] & reported, flux_std, 0.)
    number_std = safe_std * np.pi * bin_widths * 1000.0 * 6.242e11 / (energy * 1000.0)
    kernel = ionization_kernel[..., :, 1:-1]
    ionization = total_ionization[..., 1:-1]
    supported = np.any(kernel != 0, axis=-2)
    zero_supported = (ionization == 0) & supported
    if np.any(zero_supported & uncertainty_valid[..., None]):
        if not atmosphere_shape:
            ionization_variance = number_std**2 @ (kernel**2)
        else:
            ionization_variance = np.sum(number_std[..., :, None]**2 * kernel**2, axis=-2)
        uncertainty_valid &= ~np.any(zero_supported & (ionization_variance > 0), axis=-1)

    propagated = []
    for conductivity in (pedersen, hall):
        sensitivity = np.divide(conductivity[..., 1:-1], 2 * ionization,
                                out=np.zeros_like(ionization), where=ionization > 0)
        weighted = sensitivity * height_weights
        if not atmosphere_shape:
            gradient = weighted @ kernel.T
        else:
            gradient = np.sum(weighted[..., None, :] * kernel, axis=-1)
        propagated.append(gradient * number_std)
    p_variance = np.sum(propagated[0]**2, axis=-1)
    h_variance = np.sum(propagated[1]**2, axis=-1)
    covariance = np.sum(propagated[0] * propagated[1], axis=-1)
    uncertainty_valid &= np.isfinite(p_variance) & np.isfinite(h_variance) & np.isfinite(covariance)
    result.update({
        "P_std": np.where(uncertainty_valid, np.sqrt(p_variance), np.nan),
        "H_std": np.where(uncertainty_valid, np.sqrt(h_variance), np.nan),
        "PH_covariance": np.where(uncertainty_valid, covariance, np.nan),
        "uncertainty_valid": uncertainty_valid,
        "uncertainty_channel_count": channel_count
    })
    return result
