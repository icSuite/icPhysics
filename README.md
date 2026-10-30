# icPhysics

Lightweight precipitation and ionospheric-conductance physics shared by
`icBuilder` and `icAnalyzer`.

The first copied routines implement the Robinson et al. (1987) Hall and
Pedersen conductance equations and their existing uncertainty propagation.
The remaining extraction is tracked in [`TODO.md`](TODO.md).

```python
from icphysics import hall, ped

Sigma_P = ped(E0, Fe)
Sigma_H = hall(E0, Fe)
```

The crude precipitation API is also available:

```python
from icphysics import precipitation_from_ratio
from icphysics import precipitation_from_zhang_paxton
```

`ZhangPaxton2008` supplies the published empirical model. icPhysics contains
the IMAGE-specific latitude collapse and bundled Kp lookup table.

The Hardy et al. (1991) statistical ion-precipitation model is available as a
separate pure-array calculation:

```python
from icphysics import hardy_ion_precipitation

ions = hardy_ion_precipitation(kp=2, mlt=0, mlat=67)
proton_energy = ions["mean_energy"]  # keV
```

The Hardy coordinates are corrected geomagnetic latitude and MLT. Modified
Apex latitude, as used by the IMAGE pipeline, is a close but not identical
approximation.

Electron spectra can be converted to conductance with a supplied neutral
atmosphere and reference field magnitude:

```python
from icphysics import conductance_from_spectrum

conductance = conductance_from_spectrum(energy_keV, differential_energy_flux, magnetic_field_nt, atmosphere)
Sigma_P = conductance["P"]  # siemens
Sigma_H = conductance["H"]
```

Optional analytical measurement errors reuse the same calculation:

```python
conductance = conductance_from_spectrum(
    energy_keV, differential_energy_flux, magnetic_field_nt, atmosphere,
    differential_energy_flux_std=absolute_flux_std
)
Sigma_P_std = conductance["P_std"]
Sigma_H_std = conductance["H_std"]
covariance_PH = conductance["PH_covariance"]  # S squared
uncertainty_valid = conductance["uncertainty_valid"]
channel_count = conductance["uncertainty_channel_count"]  # int16
```

`absolute_flux_std` uses the same units/channel axis as flux and broadcasts
to the existing nominal batch. Propagation uses first-order derivatives and
assumes independent channel errors. It retains covariance between P/H from
the shared spectrum. Source errors may include calibration contributions;
no additional common calibration, atmosphere, field, energy or model error
is added. Following the [published DMSP convention](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JA023339),
NaN/masked error terms are omitted from variance/covariance sums; nominal P/H
still use the complete spectrum. `uncertainty_channel_count` counts finite
nonnegative supplied errors, including zero. At least one is required.
All-missing errors, supplied negative/infinite errors, or invalid nominal
conductances return NaN errors and a false uncertainty flag. Zero signal with zero errors has
zero uncertainty; nonzero errors affecting a zero-ionization height make the
square-root derivative singular and invalidate analytical uncertainty.
Small-error Monte Carlo and independent scalar finite differences verify the
propagation; large errors and low-count censoring require separate care.

The differential energy flux must have units erg cm^-2 s^-1 eV^-1 sr^-1,
with energy channels on its final axis. Leading spectrum dimensions, field
magnitudes, and atmosphere dimensions broadcast. The atmosphere is a dictionary
of profiles with shape `(..., n_altitudes)`:
`altitude_km`, `oxygen_cm3` (O), `nitrogen_cm3` (N2),
`molecular_oxygen_cm3` (O2), `mass_density_g_cm3`, and `temperature_k`.
Its shared 1D altitude grid must span 70 through 295 km. Use 1D profiles for a
shared atmosphere, or leading dimensions for one atmosphere per measurement
or ensemble member. Missing or negative channels
invalidate the complete spectrum; they are not silently zero-filled.
Species densities may be zero; temperature and mass density must be positive.
Missing or otherwise invalid atmosphere values raise `ValueError`.
Zero species are interpolated using the geometric-mean limit: interiors of
intervals touching zero are zero, with exact endpoint values preserved. This
supports MSIS atomic O below 72.5 km and explicitly resolves undefined
`log(0)` endpoint arithmetic in the original routine, without a density floor.

This translates the supplied `energy_to_conductance.idl` equations using
float64, including its pi angular factor, unidirectional Rees deposition,
altitude-dependent recombination, 30-amu ions, and interior-point integration.
It sums channel ionization before calculating electron density. It infers
conductance from spectra under those assumptions; it does not return a direct
conductance measurement. Its optional uncertainty covers supplied spectral
measurement errors under the stated assumptions. The angular assumptions
need scientific review. Tests reproduce a separate scalar transcription on a
synthetic atmosphere; the original `msis120065.xdr` atmosphere was unavailable,
so original-atmosphere or IDL-runtime parity has not been established.

## Full spatial Hardy electron option

`hardy_electron_precipitation(kp, mlt, mlat)` returns Hardy 1987 integral
`number_flux` (electrons cm^-2 s^-1 sr^-1), `energy_flux` (keV cm^-2 s^-1
sr^-1), and their ratio `mean_energy` (keV). Inputs broadcast; MLT wraps
24 hours, Kp clips 0--6, and masked/nonfinite inputs or latitudes outside the
northern 50--90 degree domain return NaN. There is no hemisphere reflection.
Tables 4/5 contain the four-parameter electron profiles, independently audited
against the supplied journal PDF; the ion model is unchanged.

The model uses corrected geomagnetic latitude. Callers supplying Modified
Apex coordinates must record that approximation. Fourier phase uses the
24-hour normalization `n*pi*MLT/12`; printed equation 5 omits pi. Fractional
Kp uses linear interpolation of evaluated log10 flux after applying the
published 1e6 equatorward/1e7 poleward backgrounds. This interpolation is a
local convention, not a prescription found in Hardy 1987. Source, coefficient
version, latitude range and Kp policy are exported as `HARDY_ELECTRON_*`
constants.

`precipitation_from_hardy(wic_corrected, dwic_corrected, kp, mlt, mlat)`
returns `E0`, `dE0`, `Fe`, `dFe`, `varE0Fe`, and `E0_response_range_flag`.
Hardy sets the raw energy; positive proton-corrected WIC sets event energy
flux in mW m^-2. SI13 is not required. The WIC response retains its existing
endpoint behavior outside 0.2--25 keV; flags are int8: 0 within, 1 below, 2 above,
-1 unavailable. `ELECTRON_RESPONSE_ENERGY_RANGE` exposes that table domain.

Uncertainties are conditional on fixed Hardy energy: `dE0=0` and
`varE0Fe=0` on nominal-valid cells exclude statistical model uncertainty.
`dFe=dWIC/W_response` propagates corrected-WIC measurement uncertainty.
Unavailable, negative or nonfinite measurement STD makes only `dFe` missing;
nominal E0/Fe remain available. No statistical Hardy predictive error is
invented. Hardy flux moments include the full spectrum (measured 50 eV--20 keV,
extrapolated 100 keV), whereas its fitted conductivities use energies >=661 eV;
using its mean energy with Robinson does not reproduce those conductance
fits. DOI: https://doi.org/10.1029/JA092iA11p12275.

## Hardy Hall and Pedersen conductance

`hardy_conductance(kp, mlt, mlat)` evaluates the independent Hardy et al.
(1987) Tables 6--7 Epstein--Fourier fits and returns auroral electron-produced
`hall` and `pedersen` conductance in siemens. It does not apply Robinson to
Hardy's full-spectrum mean energy, and it does not include solar-EUV
conductance.

```python
from icphysics import hardy_conductance

conductance = hardy_conductance(kp=2, mlt=1.5, mlat=66)
Sigma_H = conductance["hall"]
Sigma_P = conductance["pedersen"]
```

Inputs broadcast with the same northern corrected-geomagnetic 50--90 degree
domain and wrapped MLT convention as the Hardy electron model. Finite Kp clips
to 0--6. Fractional Kp linearly interpolates the evaluated conductances; this
is an explicit package convention because the paper publishes integer-Kp
states. The published 0 S equatorward and 0.55 S poleward backgrounds are
applied separately around each fitted peak. No hemisphere reflection, solar
background, or model uncertainty is supplied. Provenance and domain-policy
constants are exported as `HARDY_CONDUCTANCE_*`.

`hardy_robinson_effective_precipitation(kp, mlt, mlat)` inverts those native
Hardy conductances through the Robinson equations and returns
`effective_energy` in keV, `effective_energy_flux` in mW m^-2, and a `valid`
mask. Substituting the effective values back into Robinson exactly recovers
the returned native `hall` and `pedersen` conductances on valid points.

```python
from icphysics import hardy_robinson_effective_precipitation

effective = hardy_robinson_effective_precipitation(kp=2, mlt=1.5, mlat=66)
E_effective = effective["effective_energy"]
Q_effective = effective["effective_energy_flux"]
```

These are conductance-equivalent parameters, not the published full-spectrum
Hardy moments or a reconstruction of its unpublished moments above the
approximately 500-eV cutoff. Raw inverse values are retained diagnostically,
but `valid` rejects energy below 0.5 keV and any location where either native
conductance used the published equatorward or poleward background. For
fractional Kp, every contributing integer-Kp endpoint must come from the fitted
profile. No upper-energy or energy-flux threshold is imposed.

`collapse_hardy_robinson_effective(kp, mlt)` reduces the valid auroral
latitude profile to one conductance-consistent energy for each Kp/MLT pair.
It selects the contiguous valid interval containing the Pedersen-conductance
peak, calculates spherical-area-weighted Hall and Pedersen conductance there,
and only then applies inverse Robinson. It does not average local effective
energies and does not add a flux threshold.

```python
from icphysics import collapse_hardy_robinson_effective

collapsed = collapse_hardy_robinson_effective(kp=2, mlt=1.5)
E_collapsed = collapsed["representative_energy"]
```

The result includes effective energy/flux, collapsed Hall/Pedersen,
selected-latitude and area diagnostics, peak latitude, and explicit
valid/empty/boundary flags. The default 0.01-degree MLAT grid resolves fitted
validity boundaries rather than representing the physical model resolution.
`collapsed_hardy_robinson_effective_energy()` is a cached convenience
evaluator for many samples. It uses a periodic regular-MLT profile with a
configurable 0.05-hour default step and preserves unavailable values.

Licensed under the MIT License.
