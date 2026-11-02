# Current State

Last reviewed: 2026-11-02 (Fang 2025 conductance implemented)
Repository snapshot: `feat/hardy-conductance` from `496e841`
Worktree state: verified implementation and project-memory changes are uncommitted.

## Flux-weighted collapsed Zhang--Paxton orbit evaluator added

`collapsed_zhang_paxton_q_weighted()` now evaluates the 50--90-degree,
Q-weighted latitude collapse efficiently at arbitrary detector MLT. It rounds
Kp to the lookup convention, caches periodic 0.1-hour profiles, and returns
both energy and the Q-weighted within-slice RMS spread.

`precipitation_from_collapsed_zhang_paxton()` combines this energy with
corrected WIC counts, propagating WIC measurement error and the structural
spread through the camera response, including E0--Fe covariance. The spread
is explicitly not predictive model uncertainty. The icBuilder detector
Product-2 route now exposes this calculation as
`--precipitation-method zhang_paxton`; model provenance and response-range
flags pass through Product 3 and detector-CS products. The icPhysics suite
passes 123 tests, and 48 focused icBuilder pipeline tests pass.

Next: run the explicitly selected cZP Product-2 corpus. Do not reinterpret its
stored `dE0` as a calibrated uncertainty or change the production default.

## Zhang--Paxton Robinson/Fang conductance comparison plot added

`scripts/plot_ZP_conductance.py` writes one 3-by-3 `polplot` figure for every
integer Kp 0--9. All displayed panels start from the same full spatial Zhang--
Paxton average-energy and energy-flux fields. Hall and Pedersen rows compare
Robinson, Fang and their Fang-minus-Robinson difference. The ratio row shows
Robinson and Fang Hall/Pedersen without a difference panel. Each model pair
shares a color range, while conductance differences use zero-centred diverging
scales.

Fang uses a uniform 50-microtesla field at 150 km by default, matching its
published reference field; `--magnetic-field` makes this diagnostic assumption
explicit and configurable. The comparison uses strict Fang-domain evaluation:
energies below 1 or above 20 keV remain unavailable rather than being clipped.
The script does not calculate IGRF or claim that the fixed field represents a
historical IMAGE map. Coarse Kp 0 and 9 figures rendered successfully in the
`tConductance` environment.

## Fang 2025 conductance implemented

`fang_2025_conductance()` now evaluates Fang's equations 7--10 for Maxwellian
average energy [keV], energy flux [mW m-2] and caller-supplied 150-km magnetic-
field magnitude [microtesla]. It broadcasts arrays, preserves masks and source
energy diagnostics, applies the selected positive sub-1-keV operational floor,
and rejects nonpositive or above-20-keV energy. A strict `reject` policy is
available for the calibrated-domain cohort. Robinson and consumer defaults
are unchanged.

`fang_2025_conductance_uncertainty()` propagates explicit energy, flux and
field standard deviations and energy--flux covariance analytically, returns
P/H standard deviations and covariance, and flags uncertainty crossing the
nonlinear 1-keV floor. Independent literal equations, finite differences and
100,000-draw Monte Carlo sampling verify the calculation. The complete suite
passes 121 tests.

Next: icBuilder must supply historical IGRF magnitude at 150 km and compare
the operational and strict Fang cohorts alongside Robinson against matched
DMSP spectrum conductance before any production-model selection.

## Zhang--Paxton model plotting script added

`scripts/plot_ZP_models.py` writes one 12-panel `polplot` figure for each
integer Kp 0--9 directly from the scientific model APIs. Panels show full
spatial mean energy and energy flux, plus area-mean, area-times-Q-weighted,
conductance-effective and selected-interval maximum collapsed energies. The
maximum is the highest local energy inside the same contiguous Q > 0.05
mW m-2 latitude interval and remains a plotting-only diagnostic.

Full-model Robinson Hall/Pedersen and their ratio are joined by strict-domain
Fang Hall/Pedersen and Fang Hall/Pedersen ratio. The two ratio panels share a
color range. Fang uses a configurable uniform 150-km field, defaulting to 50
microtesla, and remains unavailable outside 1--20 keV without clipping.
Collapsed MLT profiles are repeated over latitude only for polar display. All
energy panels share one range and all four conductance panels share another.

Kp 0 and 9 12-panel smoke figures rendered successfully with the standard
0.01-degree collapse grid. No lookup table was read, generated or modified.

## Three Zhang--Paxton collapse methods implemented

The 50--90-degree collapse now returns area-mean, area-times-Q-weighted, and
Robinson-conductance-equivalent E0 definitions over one shared contiguous
Q > 0.05 mW m-2 interval. Matching structural MLAT spreads are supplied for
all three. The conductance-equivalent method also returns effective Q and uses
area-times-sqrt(Q) for its RMS local-energy deviation about effective E0.
These spreads are not predictive/model uncertainties.

Both generator modes and the detector schema can write all three methods.
Loaders expose explicit `area_mean`, `q_weighted`, and
`conductance_effective` variants, reject mismatched/ambiguous fields, and
retain the old legacy and `collapsed` area-mean APIs. Synthetic and loader
fixture tests plus the complete suite pass: 108 tests.

No lookup-generation function or CLI was run and no NetCDF lookup/template
artifact was changed during this checkpoint. Existing live lookup artifacts
remain separate concurrent/recent work; consumer selection and scientific
comparison of the three methods are still pending.

## Detector-grid Zhang--Paxton lookup tooling implemented

The 46-by-46 plan is implemented in code: canonical checksummed coordinates,
separate legacy/detector generator modes, full and collapsed detector fields,
atomic detector writing, explicit grid/model loader selection and focused
validation tests. The full model evaluates only the 1,747 cells inside
50--90 degrees and supplies no invented uncertainty. The collapsed model uses
a verified periodic 0.1-hour MLT profile; maximum checked errors are 0.001230
keV for mean, 0.000589 keV for spread and 0.004835 keV for median.

The legacy 36-by-36 API and frozen values remain passing. The complete suite
passes 102 tests. No 901-layer detector lookup was generated, and no consumer
or production default changed. A detector-CS orbit smoke test therefore remains
pending. See [[Zhang Paxton 46x46 Lookup Plan]] for the schema and evidence.

## Zhang--Paxton latitude boundary retained

Confirmed: 50 degrees MLAT remains the hard equatorward limit for the
Zhang--Paxton collapse. A bounded diagnostic found that extending to 40
degrees did not close every thresholded interval and that closure at 20
degrees would require unsupported low-latitude extrapolation. The partial
20-degree source/test changes were removed and lookup generation was stopped;
the existing 50--90-degree, fixed 36-by-36 lookup was not modified.

The implementation requirements for a separate 46-by-46 detector-grid lookup
are recorded in [[Zhang Paxton 46x46 Lookup Plan]]. It will carry both local
full-model energy and an explicitly named collapsed energy/spread/median while
preserving the legacy 36-by-36 artifact.

## Hardy model plotting script added

Confirmed: isolated local Robinson-effective energies reach thousands to
millions of keV near a zero of the independent Pedersen fit while Hall stays
finite, dominating automatic color limits. The inversion round-trips correctly;
the existing `valid` mask does not ensure physical energy plausibility. See
[[Hardy Electron Model]] for reproduced values. Per user choice, the local
Robinson-effective mean-energy panel now caps displayed values at 25 keV,
with the cap stated in its title. A Kp 4 smoke plot rendered successfully.
Underlying model arrays, other panels and scientific validity are unchanged.

`scripts/plot_hardy_models.py` now writes one ten-panel `polplot` figure for
each integer Kp level 0--6. The panels show native Hardy Hall and Pedersen
conductance, published energy and number flux, their derived mean-energy ratio,
local Robinson-effective mean energy and flux, and collapsed Robinson-effective
mean energy. Invalid local effective values are masked; the latitude-independent
collapsed result is repeated over MLAT only for polar display. The additional
panels show full spatial Zhang--Paxton mean energy from `zhangpaxton2008` and
the existing `collapse_zhang_paxton` representative energy, repeated over MLAT.
The latter uses its standard 0.05 mW m-2 threshold and 0.01-degree collapse
grid. The figure layout is 3 by 4, with two unused axes hidden. Coarse-grid
smoke figures for Kp 0 and 6 rendered successfully with all ten panels.

The script uses explicit cell edges so the plotted 50--90-degree domain covers
the full polar disk without a midnight seam. A default-grid smoke run produced
all seven PNGs, and the full package suite passes 97 tests. `polplot` and
Matplotlib remain script-only dependencies; the lightweight scientific package
dependency boundary is unchanged.

## Collapsed Hardy Robinson-effective model implemented

`collapse_hardy_robinson_effective(kp, mlt)` now reduces each native Hardy
latitude profile to one conductance-consistent representative energy. It finds
the Pedersen peak among locally valid cells, selects only its contiguous valid
interval, spherical-area averages Hall and Pedersen conductance, and then
applies inverse Robinson. It deliberately does not average local effective
energies, extrapolate across invalid gaps, or add an upper-energy/flux cutoff.

The full diagnostic result includes representative energy/flux, collapsed
Hall/Pedersen, selected MLAT bounds, peak MLAT/Pedersen, selected area and cell
count, and valid/empty/sampling-boundary flags. The default MLAT step is 0.01
degrees. `collapsed_hardy_robinson_effective_energy()` caches a periodic
0.05-hour MLT profile for efficient crossing evaluation; all grid choices are
explicit and NaNs remain unavailable.

All seven integer Kp levels on a 0.5-hour MLT grid produced valid solutions;
representative energy ranged 1.99--20.51 keV. The full suite passes 95 tests,
including independent inverse-of-mean behavior, resistance to extreme local
energy spikes, disconnected-region rejection, native-background exclusion,
exact Robinson round trips, periodic interpolation, broadcasting, fractional
Kp, masks, empty slices and grid validation. Existing Hardy APIs remain
numerically unchanged.

Next: use the cached collapsed effective energy in the bounded icBuilder Paper
1 comparison and retain the current full-spectrum collapse as a diagnostic.

## Hardy Robinson-effective precipitation implemented

The public array API
`hardy_robinson_effective_precipitation(kp, mlt, mlat)` now inverts the native
Hardy Hall and Pedersen fits through the exact Robinson equations. It returns
raw conductance-equivalent energy [keV] and energy flux [mW m-2], the source
conductances, per-conductance background flags and a validity mask. Valid
values round-trip to native Hardy Hall/Pedersen conductance to numerical
precision.

Validity requires positive finite conductances, energy >=0.5 keV and genuine
fitted-profile Hall and Pedersen values. Published 0-S equatorward and 0.55-S
poleward backgrounds are detected from the clipping decisions rather than
output equality. Every contributing integer-Kp endpoint must be profile-based
for fractional Kp. Raw inverse values remain available diagnostically; no
upper-energy or energy-flux threshold is imposed. These are Robinson-effective
parameters, not the published full-spectrum moments or reconstructed
unpublished cutoff moments.

The existing native conductance and full-spectrum electron APIs are unchanged.
Focused tests include exact round trips, backgrounds,
fractional-Kp policy, the unclipped 0.5-keV validity boundary, broadcasting,
masked/nonfinite inputs and existing frozen Hardy conductance checkpoints.

Next: integrate the effective-energy option in the selected consumer and
quantify its empirical behavior before adopting it as a production default.

## Hardy Hall and Pedersen conductance implemented

The public array API `hardy_conductance(kp, mlt, mlat)` now evaluates the
Hardy et al. (1987) Tables 6--7 empirical Hall/Pedersen conductance fits
directly. This is not Robinson applied to Hardy's full-spectrum electron mean
energy. It returns the auroral electron contribution in siemens, excluding
solar EUV and unquantified predictive uncertainty.

All 728 coefficients were cross-checked against the supplied journal PDF and
an established independent local transcription. Independent scalar equations,
published maxima, boundaries, backgrounds, broadcasting, MLT wrapping,
fractional/clipped Kp and invalid inputs are tested. A wheel build succeeded
and includes both coefficient CSVs. No consumer,
production default, server installation, commit, merge or push was changed.

Next: review and commit the feature branch, then add bounded icBuilder and
icAnalyzer integration only when a consumer is selected. Predictive model
uncertainty and corrected-geomagnetic/Modified-Apex sensitivity remain open.

## Main integration verified

User approved integrating the spectrum branch into main. Fetched GitHub refs
and confirmed main was an ancestor; fast-forwarded c3484fd to 8db326f without
conflicts. This includes eecc263 measurement-uncertainty handling, spectrum
conductance/uncertainty and the optional Hardy electron model. The full package
suite passed 72 tests in tConductance on 2026-10-12 before the unchanged code
was fast-forwarded. Existing feature branches are retained.

Next: update the server from origin/main and reinstall the shared package in
its processing environment. No server installation, corpus job or default
electron-model change was performed. Broader Hardy validation remains open.

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
transfer them to the server and rebuild crossings there, then rerun CS analysis.
No production writes/full-year runs were performed during implementation.
Variance of omitted channels is not modelled; this is the selected published
missing-error convention. This implementation is committed in main.

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
