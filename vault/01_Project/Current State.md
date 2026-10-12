# Current State

Last reviewed: 2026-10-11 (Hardy electron implementation; spectrum checkpoint retained)
Repository snapshot: `feat/spectrum-conductance` at `eecc263`
Worktree state: Hardy electron and spectrum physics/export/tests/README and memory updates are
uncommitted; pre-existing memory changes preserved.

## Hardy electron model implemented

Confirmed 2026-10-11: the user selected full spatial MLAT/MLT/Kp evaluation
and requested implementation using the supplied 1987 functional paper.
`hardy_electron_precipitation()` now returns the published number/energy flux
and mean-energy ratio; `precipitation_from_hardy()` combines that raw energy
with corrected WIC. [[Hardy Electron Model]] owns sources, equations, domain,
Fourier normalization, Kp convention, response flags and conditional errors.
All 728 coefficients independently audited; independent scalar checks agree
within 4.33e-14 relative. Full package suite 72 tests passes. Existing ion and
spectrum work is preserved. Builder's real 122-frame orbit passed downstream;
no default model, corpus, commit or branch changed.

Next: review empirical spatial/particle validation and coordinate sensitivity
before production selection; predictive Hardy error remains unquantified.

## Spectrum-inferred electron conductance

**Confirmed by the user (2026-10-06):** omit NaN channel-error terms following
the published DMSP convention. The previous all-channel-errors requirement
and pending alternatives are superseded.

Optional differential_energy_flux_std now returns P_std/H_std, PH_covariance,
uncertainty_valid and int16 uncertainty_channel_count. Finite nonnegative
errors contribute (including defined zeros); NaN/masked terms are omitted.
All-missing errors remain unavailable. Supplied negative/infinite errors,
invalid nominal spectra and singular zero-ionization derivatives invalidate
uncertainty. Channel count reports finite nonnegative supplied error terms,
even if nominal/error propagation is invalid. No-STD calls retain the P/H-only
API and nominal P/H stay bit-for-bit unchanged. Shared/batched atmospheres
remain supported. Independent channels, fixed atmosphere/B/energies and no
additional common calibration/model covariance are explicit assumptions.

63 physics tests pass; partial errors exactly match omitted-term references
including covariance, masking, batching and channel reversal. Independent
scalar/Monte Carlo derivative checks remain passing. Builder's135 focused tests
and real daily/crossing/analysis gates pass: Jan1 uncertainty available in
19,696/20,316 nominal-valid spectra; all870 selected real CS test pairs eligible.
The620 allmissing daily spectra remain unavailable. Detailed source/numerical
and scratch evidence belongs to icBuilder's DMSP Spectrum Conductance Validation.

Confirmed workflow: generate yearly DMSP files locally with updated packages,
transfer them to the server and rebuild crossings there, then rerun CS analysis. No production writes/full-year runs or commits made here.
Variance of omitted channels is not modelled; this is the selected published
missing-error convention. All source/tests/docs/memory remain uncommitted.

The user approved a Python port of `energy_to_conductance.idl` and its use in
icBuilder yearly DMSP processing. `conductance_from_spectrum` is exported from
`src/icphysics/spectrum.py`: explicit atmosphere profiles, channel energies in
keV, differential energy flux in erg cm-2 s-1 eV-1 sr-1, and magnetic-field
magnitude in nT give P/H in S. Atmosphere fields can be shared 1D or batched
`(..., altitude)` on one altitude grid; leading atmosphere/spectrum/B dimensions
broadcast together. The shared-profile fast path remains unchanged.
The port preserves the inherited Rees/unidirectional deposition, pi angular
factor, equilibrium recombination, 30-amu ions, and integration endpoints.
It uses float64 and invalidates missing/negative spectra rather than silently
turning absent channel ionization into zero. Zero neutral species use the
geometric interpolation limit; density/temperature remain positive.

All 63 tests pass. Independent scalar transcription agrees below 1.56e-15
relative on supplied synthetic profiles, including broadcast mixed dimensions;
batched and separate-profile calls agree exactly. Isolated source float32-Q
rounding changes P/H by about +4.2e-10/-3.3e-10 relative. The original XDR
atmosphere and actual IDL-runtime parity remain unavailable. Optional spectral measurement errors are now returned as described above.
README documents the API and physical assumptions.
Atmosphere generation, IGRF and DMSP I/O remain in icBuilder. Its 50 focused
index/integration tests and full daily CDF scratch conversion passed. The
default yearly pathway now uses local historical indices and user-selected
2° latitude × 5° longitude × 1-hour shared MSIS profiles, retaining an exact
per-measurement comparison option. No file/time/model acquisition was added
to this numerical package. The shared-1D-profile numerical fast path is reused.

## Measurement-uncertainty contract

`proton_correct_images()` now accepts an explicit `uncertainty_mode`.
`legacy_add_poisson`, the backward-compatible default, preserves the extracted
historical equations that add image signal as a counting-variance term.
`measurement` treats the supplied WIC/SI12/SI13 standard deviations as full
detector measurement uncertainty and does not add counting noise again. Proton
response and flux uncertainty are still propagated into corrected WIC/SI13.

icBuilder detector Product 2 uses `measurement` because fuvpy Product-1
variance already contains Poisson or compound-Poisson detector noise. A
negative-background test verifies that this mode remains finite and matches
the direct uncertainty equations; the frozen legacy comparison is unchanged.
All 17 icPhysics tests pass. On full icBuilder orbit 0085, central Fp/E0/Fe
and P/H are unchanged while propagated uncertainties become finite on every
method-valid cell.

## Current position

The crude shared-physics implementation is in place. Alongside Robinson, the
package now contains the copied IMAGE response tables, SI12 proton correction,
WIC/SI13 ratio method, WIC energy-flux conversion, uncertainty and E0--Fe
covariance, Zhang--Paxton collapse, bundled lookup table, loader, and lookup
generator.

The full Hardy et al. (1991) Table 1 is transcribed in
`hardy_coefficients.py`, and `hardy.py` now implements the published Fourier,
generalized-Epstein, background-limit, and Kp-interpolation equations. The
public `hardy_ion_precipitation()` function returns ion number flux, energy
flux, and mean energy for broadcast-compatible Kp, MLT, and latitude arrays.
The 1092 coefficients were visually checked against the rendered paper table.
The evaluator reproduces the morphology of the paper's Kp 0, 2, and 4 model
maps and agrees within 12% at four independent mean-energy checkpoints from
the original 1989 statistical tables.

Hardy is now connected to icBuilder's modular Product-2 proton correction.
icPhysics publicly exposes the 0.47--46.7 keV domain of the Frey proton camera-
response tables so the caller can preserve raw model energy while clipping the
response input explicitly. The response functions accept spatial Ep/dEp
arrays; SI12 still supplies the event-specific proton flux.

`proton_correct_images()` now exposes SI12-based proton correction separately.
`precipitation_from_ratio()` accepts corrected WIC/SI13 arrays, while
`precipitation_from_zhang_paxton()` accepts corrected WIC plus Kp and MLT.
icBuilder performs correction once before selecting either precipitation
method. A frozen scalar comparison matches the legacy combined ratio routine.

`robinson_conductance()` is the array entry point for Product 3. It returns
Hall/Pedersen conductance and propagated uncertainty, including the established
one-sided uncertainty at zero energy flux.

Zhang--Paxton ownership is now settled: `ZhangPaxton2008` retains the published
base model, while icPhysics will own the latitude collapse, uncertainty and
sampling, lookup generation, bundled lookup file, and loader. icBuilder will
eventually call icPhysics and contain no Zhang--Paxton implementation.

## Confirmed boundary

- icAnalyzer must not depend on the full icBuilder distribution because that
  would propagate pipeline dependencies such as fuvpy.
- icPhysics sits below both projects and depends on neither.
- `ZhangPaxton2008` remains the authoritative implementation of the published
  Zhang--Paxton equations.
- icBuilder owns IMAGE preprocessing, coordinate handling, Kp assignment,
  orchestration, and files.
- icAnalyzer owns VAE and statistical-analysis workflows.
- Shared numerical equations should be pure, array-based, lightweight, and
  tested against the existing implementation before replacement.

## Not yet decided

- final distribution and import name;
- initial module layout and public function names;
- minimal declared dependency set;
- versioning and first alpha release criteria.

## Next action

Review/install the new spectrum function in the DMSP processing environment,
then regenerate selected yearly products after the reference atmosphere is
chosen. Original-atmosphere recovery remains useful for legacy parity.
The earlier measurement-uncertainty branch merge/server installation work is
not resolved by this checkpoint and remains a separate deployment action.

The pre-existing icAnalyzer batched-array import check, cleanup, and removal of
old icBuilder copies remain deferred.
