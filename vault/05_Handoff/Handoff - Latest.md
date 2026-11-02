# Handoff - Latest

Last updated: 2026-11-02
Repository snapshot: `feat/hardy-conductance` from `496e841`
Worktree state: verified implementation and memory changes are uncommitted.

## Latest checkpoint: cZP detector Product-2 route is ready

The public `collapsed_zhang_paxton_q_weighted(kp, mlt)` API evaluates the
flux-weighted 50--90-degree Zhang--Paxton collapse at arbitrary MLT using
cached periodic 0.1-hour profiles. Kp follows the lookup convention by
rounding to the nearest hundredth. It returns the collapsed energy and the
Q-weighted RMS latitude variation discarded by the collapse.

`precipitation_from_collapsed_zhang_paxton()` combines those values with
corrected WIC and its measurement uncertainty, returning E0, dE0, Fe, dFe,
their covariance and the camera-response range flag. Here dE0 remains a
structural sensitivity scale, not predictive Zhang--Paxton uncertainty.

icBuilder exposes the route as `--precipitation-method zhang_paxton` with the
derived label `CollapsedZhangPaxtonQWeightedElectron_HardyProton...` when
Hardy proton energy is retained. Scientific metadata and response flags pass
through detector precipitation, detector conductance and both CS products.
The icPhysics suite passes 123 tests; 48 focused icBuilder Product-2/Product-3/
CS tests pass. A real four-point calculation at Kp 2 and MLT 0/6/12/18 hours
returned finite E0, dE0 and WIC-derived flux.

Next: run the selected cZP Product-2 corpus. This enables a bounded run; it
does not select cZP as the production default or resolve predictive energy
uncertainty.

## Portfolio impact

- Central update needed: Yes
- Changes: flux-weighted collapsed Zhang--Paxton is now available as an
  explicit icBuilder detector Product-2 configuration.
- Sync summary: a cZP Product-2 run can proceed; defaults are unchanged and
  stored dE0 remains structural rather than calibrated uncertainty.

## Latest checkpoint: Zhang--Paxton Robinson/Fang conductance plot added

`scripts/plot_ZP_conductance.py` produces one 3-by-3 polar figure for every
integer Kp 0--9. Each figure evaluates one shared full Zhang--Paxton average-
energy and energy-flux map. The top row shows Robinson Hall, Fang Hall and
Fang-minus-Robinson Hall; the middle row shows the corresponding Pedersen
fields. The bottom row shows Robinson Hall/Pedersen, Fang Hall/Pedersen and the
third axis is unused. Each model pair shares its own per-Kp color range, and
the Hall and Pedersen differences use symmetric zero-centred diverging ranges.
Explicit cell edges preserve the complete 50--90-degree MLAT disk and periodic
MLT seam.

Fang uses its 50-microtesla 150-km reference field by default. The uniform
field is exposed as `--magnetic-field`; this diagnostic does not calculate or
approximate a historical spatial IGRF map. The plot uses strict Fang-domain
evaluation: energy below 1 or above 20 keV is unavailable, with no operational
floor. Coarse Kp 0 and 9 smoke figures rendered successfully; syntax and
whitespace checks pass. No repository figure, lookup artifact, model API or
default changed.

## Portfolio impact

- Central update needed: No
- Changes: None
- Sync summary: added a diagnostic Zhang--Paxton Robinson/Fang comparison
  plot; research priority, consumer state and scientific defaults are
  unchanged.

## Latest checkpoint: Fang 2025 conductance implemented

`fang_2025_conductance()` implements equations 7--10 and Tables 2--3 from the
supplied Fang (2025) paper. Inputs are Maxwellian average energy [keV], energy
flux [mW m-2] and caller-prepared magnetic-field magnitude at 150 km
[microtesla]. NumPy inputs broadcast, masked/nonfinite and nonphysical values
remain unavailable, valid zero flux returns exact zero, and returned masks
distinguish the published 1--20-keV domain from the operational positive
sub-1-keV floor. `low_energy_policy="reject"` supports the strict cohort;
above-20-keV energy is never clipped or extrapolated.

`fang_2025_conductance_uncertainty()` separately propagates energy, flux and
field standard deviations plus energy--flux covariance with analytical
Jacobians. It returns P/H standard deviations and covariance without changing
nominal results when errors are invalid, sets zero energy sensitivity below
the floor, flags one-standard-deviation intervals crossing 1 keV, and rejects
the singular case of nonzero flux uncertainty at zero nominal flux.

Every coefficient was checked against the local PDF. Independent literal
equations and scalar checkpoints, domain and mask cases, multidimensional
broadcasting, separate Hall/Pedersen field exponents, finite differences and
a 100,000-draw Monte Carlo comparison pass. The complete suite passes 121
tests. README and public exports are updated. No consumer, IGRF calculation,
lookup generation, production output or default was changed.

Next: integrate historical IGRF at 150 km in icBuilder's bounded Paper 1
comparison and report both operational and strict Fang-domain cohorts against
matched DMSP spectrum conductance before considering a default change.

## Portfolio impact

- Central update needed: Yes
- Changes: the optional Fang 2025 numerical and uncertainty APIs are verified
  and ready for bounded icBuilder comparison; project priority and defaults
  are unchanged.
- Sync summary: consumer IGRF integration and empirical Fang-versus-Robinson
  validation are the next actions; no production selection has been made.

## Latest checkpoint: Zhang--Paxton model plotting script added

`scripts/plot_ZP_models.py` produces one 12-panel polar figure for every
integer Kp 0--9. It evaluates the full Zhang--Paxton mean-energy and energy-
flux fields directly, then displays area-mean, Q-weighted, inverse-Robinson
conductance-effective and selected-interval maximum collapsed energies. The
maximum is evaluated over the same contiguous Q > 0.05 mW m-2 latitude
interval and is not added to the scientific collapse API or lookup schema.

The script applies Robinson to full E/Q for Hall, Pedersen and H/P, and adds
strict-domain Fang Hall, Pedersen and H/P using a configurable uniform 150-km
field (50 microtesla by default). The Robinson and Fang ratio panels share a
color range. Fang values outside 1--20 keV remain unavailable; no lower-energy
floor is applied. Collapsed values are repeated over MLAT only for display.
Energy panels share a range and all four conductance panels share another.

Coarse display-grid 12-panel smoke figures for Kp 0 and 9 rendered successfully
while retaining the standard 0.01-degree scientific collapse grid. Syntax and
`git diff --check` pass. No lookup artifact was read or generated.

## Portfolio impact

- Central update needed: No
- Changes: None
- Sync summary: Added and smoke-tested a diagnostic Zhang--Paxton comparison
  plot; model APIs, lookup artifacts and production choices are unchanged.

## Latest checkpoint: three Zhang--Paxton collapse methods implemented

`collapse_latitude_slice()` and `collapse_zhang_paxton()` now provide three
50--90-degree reductions over the same contiguous Q > 0.05 mW m-2 interval:
area-weighted arithmetic mean, area-times-Q-weighted mean, and a
conductance-equivalent method that area-averages local Robinson Hall/Pedersen
and inverse-Robinson recovers effective E0/Q. Their structural dE0 values are,
respectively, the matching weighted RMS spreads and an area-times-sqrt(Q)-
weighted RMS local-E0 deviation centered on effective E0. The sqrt(Q) factor
matches the common Robinson conductance amplitude and does not prefer H or P.

Both legacy and detector generator paths now carry explicit fields for all
three energies/spreads and effective Q. Historical `E0`/`dE0` and detector
`E0_collapsed`/`dE0_collapsed` remain exact area-mean aliases. Loaders expose
`area_mean`, `q_weighted`, and `conductance_effective`, retain `collapsed`, and
reject grid mismatch, missing variants, or inconsistent explicit/alias area
fields.

Focused synthetic/fixture tests cover formulas, weighting centers, exact
Robinson conductance round trips, selection, shape, schema, all variants,
legacy values and ambiguity rejection. The complete suite passes 108 tests.
No generator function, CLI, multiprocessing calculation, or smoke generation
was run.

Lookup/template SHA-256 values recorded before and after this task are
unchanged:

- coordinate template: `4a5f2cc233334bc8c269e4fa30c7ee6c8b31e0f9f9e6b3cbe418d25a6b162e7c`;
- detector lookup: `7aa440a052adf0e248e5283d4faab9641792dc21d497c9b1d6585155084887bc`;
- legacy lookup: `ba9db0f96c80312a832761fae023157dbb8984c709a155316d356e8259bf41c8`.

Scientific caveat: every dE0 is within-slice structural variability discarded
by that definition. None represents Zhang--Paxton coefficient uncertainty,
model bias, Kp uncertainty, or full predictive uncertainty. The preferred
collapse must be selected from the downstream DMSP comparison; this
implementation does not choose it.

Next: generate artifacts only when explicitly authorized, then compare all
three variants in the Paper 1 consumer before selecting a product default.

## Portfolio impact

- Central update needed: No
- Changes: None
- Sync summary: three explicit collapsed Zhang--Paxton variants are ready for
  later artifact generation and comparison; no priority, deadline, or
  production default changed.

## Latest checkpoint: 46x46 Zhang--Paxton lookup tooling implemented

The agreed detector-grid lookup design is implemented without running the
901-layer generator. A 75-KB coordinate template freezes icBuilder grid
`image_apex_130km_46x46_v1`, its provenance and SHA-256. The generator now has
separate legacy/detector modes, skips full-model evaluation below 50 degrees,
produces explicit full/collapsed fields, records verified interpolation error,
validates all arrays and provenance, and atomically installs only a complete
detector artifact.

The loader preserves the legacy no-argument 36-by-36 result and requires
explicit `grid_id` and `energy_variant` for the detector table. It rejects grid,
shape, checksum, domain, provenance and nonfinite-data mismatches. Full output
has `E0_full`, diagnostic `Q_full` and `full_valid`, with no `dE0`; collapsed
output has its mean, MLAT spread, median and equatorward-boundary flag.

The 0.1-hour collapsed profile was compared with direct collapse at 240 actual
grid-MLT/Kp samples. Maximum errors are 0.001230 keV (mean), 0.000589 keV
(spread) and 0.004835 keV (median); boundary flags had zero mismatches. The
complete suite passes 102 tests and `git diff --check` passes.

No detector lookup table was created. The user explicitly clarified that the
implementation request did not authorize running the generator. Generation,
artifact-level random-cell checks and the one-orbit detector-CS consumer smoke
test remain pending.

## Portfolio impact

- Central update needed: No
- Changes: None
- Sync summary: 46x46 Zhang--Paxton lookup tooling and canonical coordinates
  are implemented and tested; artifact generation and consumer integration are
  pending, with no production selection or priority/deadline change.

## Latest decision: retain the Zhang--Paxton 50-degree boundary

The user confirmed 50 degrees MLAT as the hard equatorward limit. Do not
extrapolate Zhang--Paxton below its validated high-latitude domain to force the
0.05 mW m-2 selection interval to close. The interrupted 20-degree changes
were removed, lookup generation stopped, and the existing 50--90-degree
36-by-36 lookup remains byte-for-byte unchanged.

The agreed follow-on design is documented in
[[Zhang Paxton 46x46 Lookup Plan]]. It adapts the existing generator for the
frozen detector-CS grid and stores distinct full and collapsed Zhang--Paxton
energy fields without inventing full-model uncertainty or replacing WIC flux.

## Latest checkpoint: Hardy model plotting script added

Boundary diagnosis 2026-11-02: all Kp maps contain isolated divergent local
effective energies where the Pedersen fit approaches zero with finite Hall.
Median valid energies are 3.39--4.78 keV but maxima reach 5.10e3--1.20e6 keV
on the default plotting grid. Robinson round trips pass; `valid` screens
backgrounds and the lower energy bound but does not ensure physical energy
plausibility. [[Hardy Electron Model]] owns the numerical evidence. The
collapsed model check remains 1.99--20.51 keV. No masking policy or code
change was made during this diagnosis.

User-selected display update: the Robinson-effective mean-energy panel clips
values above 25 keV to 25 keV and states the cap in its title. Clipping occurs
only during figure construction; calculated model arrays and other panels are
unchanged. A Kp 4 smoke plot rendered successfully.

`scripts/plot_hardy_models.py` writes one PNG for each integer Kp level 0--6.
Each figure contains ten `polplot` panels: Hardy Hall conductance, Pedersen
conductance, published energy flux, published number flux, Hardy mean energy
derived as energy flux divided by number flux, local Robinson-effective mean
energy, local Robinson-effective flux, collapsed Robinson-effective mean
energy, full spatial Zhang--Paxton mean energy, and collapsed Zhang--Paxton
mean energy. Both collapsed MLT profiles are repeated over latitude for
display; local effective values are masked using the model's validity flag.

The Zhang--Paxton full model calls the maintained `zhangpaxton2008` package.
Its collapsed panel uses `collapse_zhang_paxton` with the default contiguous
peak interval above 0.05 mW m-2 and 0.01-degree collapse grid. A 3-by-4 layout
keeps the original eight panels in place and hides the two unused axes.
Coarse-grid Kp 0 and 6 smoke figures rendered successfully with all ten panels.

The plotting grid uses explicit cell edges with `plotimg(..., corr=False)`, so
the full 50--90-degree latitude domain and periodic MLT seam are rendered.
The default 0.1-hour by 0.25-degree run produced all seven figures. The full
package suite passes 97 tests. Matplotlib and `polplot` remain script-only
dependencies, with no change to the package's scientific dependency boundary.

Next: review figures with the chosen 25-keV display cap. Whether consumers need
an additional physical-validity criterion remains a separate open question.

## Portfolio impact

- Central update needed: No
- Changes: None
- Sync summary: Added a verified diagnostic plotting script; scientific APIs,
  project priority, deadlines and production defaults are unchanged. Diagnosed
  local effective-energy divergence at near-zero Pedersen boundaries.

## Latest checkpoint: collapsed Hardy effective energy implemented

`collapse_hardy_robinson_effective(kp, mlt)` now evaluates native Hardy
conductance over a regular 50--90-degree MLAT grid, finds the Pedersen peak
among background-screened valid cells, selects its contiguous valid interval,
and spherical-area averages Hall and Pedersen before inverse Robinson. This
produces one conductance-consistent energy/flux pair without averaging local
effective energies or allowing fitted-boundary spikes and disconnected valid
regions to dominate.

The full result includes representative energy/flux, collapsed Hall/Pedersen,
selected MLAT limits, peak and area diagnostics, selected cell count, and
valid/empty/sampling-boundary flags. The default MLAT step is 0.01 degrees;
there is no upper-energy or flux threshold. Empty profiles remain NaN and no
poleward extrapolation occurs. The cached convenience API
`collapsed_hardy_robinson_effective_energy()` uses a configurable periodic
MLT profile with a 0.05-hour default step.

The complete suite passes 95 tests. Independent synthetic tests prove that the
result is inverse(mean Hall, mean Pedersen), not mean(local energy), remains
stable beside a 1000-keV local boundary spike, and rejects disconnected cells.
Native Hardy backgrounds, round trips, fractional Kp, scalar/broadcast/masked
inputs, periodic MLT, empty slices and grids are covered. Existing Hardy APIs
and frozen values remain passing. A 7-Kp x 48-MLT sanity grid was entirely
valid with representative energy 1.99--20.51 keV. Default cached-profile
construction took about 3.0 s for one Kp locally. `git diff --check`, syntax
parsing and the full suite pass. No consumer, commit, merge or default changed.

Next: integrate the optional cached energy in the bounded icBuilder Paper 1
comparison, while retaining the old full-spectrum collapse as a diagnostic.

## Portfolio impact

- Central update needed: Yes
- Changes: verified collapsed, conductance-consistent Hardy effective-energy
  API is ready for consumer integration; priorities and deadlines unchanged.
- Sync summary: icBuilder integration and empirical comparison remain.

## Latest checkpoint: Hardy Robinson-effective precipitation implemented

`hardy_robinson_effective_precipitation(kp, mlt, mlat)` now wraps the native
Hardy Hall/Pedersen conductance model and exactly inverts the Robinson
equations to return effective energy [keV] and energy flux [mW m-2]. It also
returns source conductances, explicit Hall/Pedersen background flags and a
validity mask. Existing `hardy_conductance()` and full-spectrum Hardy electron
outputs are unchanged.

Raw inverse values are retained wherever positive finite conductances permit
inversion. Validity additionally requires energy >=0.5 keV and rejects any
point where either conductance used Hardy's published equatorward or poleward
background. Background provenance comes from the Epstein evaluator's clipping
decision. At fractional Kp, every integer endpoint with nonzero interpolation
weight must be profile-based. No upper-energy or flux threshold was added.
The quantities are conductance-equivalent, not Hardy's published full-spectrum
moments or a reconstruction of its unpublished cutoff moments.

The full suite passes 83 tests. Focused checks cover exact Robinson round trips,
both backgrounds, conservative fractional-Kp flags, the raw-but-invalid
sub-0.5-keV boundary, broadcasting, masked/nonfinite inputs and frozen native
Hardy values. Syntax and `git diff --check` pass. No consumer, production
default, commit, merge, push or server state changed.

Next: integrate this optional effective-energy model in icBuilder and compare
it empirically before selecting a production energy source.

## Portfolio impact

- Central update needed: Yes
- Changes: optional conductance-consistent Hardy effective energy/flux API is
  implemented and verified; no portfolio priority or deadline changed.
- Sync summary: icBuilder integration and empirical model comparison remain.

## Latest checkpoint: direct Hardy conductance model implemented

`hardy_conductance(kp, mlt, mlat)` now provides the Hardy et al. (1987)
Tables 6--7 Hall/Pedersen Epstein--Fourier fits as a broadcast array API.
It returns electron-produced conductance in siemens and deliberately excludes
solar EUV, hemisphere reflection and unsupported predictive uncertainty.
It is independent of Robinson and Hardy's full-spectrum precipitation moments.

The implementation packages both 91-row coefficient tables, applies the
published 0 S equatorward/0.55 S poleward backgrounds, wraps MLT, clips finite
Kp 0--6 and linearly interpolates evaluated conductance for fractional Kp.
Latitude support is northern corrected geomagnetic 50--90 degrees. All 728
coefficients were checked against the supplied PDF and an independent local
transcription. Independent scalar-grid and published-maxima tests pass.
The full suite passes 78 tests; a wheel build succeeds and contains the CSVs.

Next: review/commit the feature branch. Consumer integration, default changes,
server installation, predictive uncertainty and coordinate sensitivity were
not part of this checkpoint.

## Portfolio impact

- Central update needed: Yes
- Changes: optional direct Hardy Hall/Pedersen conductance capability is
  implemented and verified; no priority, deadline or production default changed.
- Sync summary: feature branch is uncommitted and not integrated into main.

## Latest checkpoint: spectrum and uncertainty branches integrated into main

User authorized merging feat/spectrum-conductance into main and publishing the
result. GitHub fetch succeeded; main c3484fd was an ancestor of 8db326f, which
already contains fix/measurement-uncertainty eecc263. Main fast-forwarded to
8db326f without conflicts. No code edits or branch deletions were needed.
The full suite passed 72 tests in tConductance on 2026-10-12. Publication is
part of the authorized closeout; verify the final main/origin/main hashes in
live Git. Current State and Decision Log record the integration.

Next: server fetch/switch main/pull --ff-only and reinstall icPhysics in the
processing environment. Retain the independently updated icBuilder/icReader
branches. No server login/install, production run or default-model change.
Hardy predictive uncertainty and empirical/coordinate validation remain open.

## Portfolio impact

- Central update needed: Yes
- Changes: measurement-uncertainty, spectrum and optional Hardy electron
  implementations are integrated on main; server installation remains.
- Sync summary: 72 package tests pass. The previous feature-branch merge
  prerequisite is resolved; no priority, deadline or model-policy change.

## Previous checkpoint: Hardy electron model implemented

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

## Portfolio impact

- Central update needed: Yes
- Changes: Optional full spatial Hardy electron capability implemented; no
  priority, deadline or default production-model change.
- Sync summary: Numerical/source and integration gates passed. Broader
  empirical/coordinate validation and predictive uncertainty remain open.
  Existing unrelated checkpoint evidence remains as previously recorded.

## Previous checkpoint: published NaN-error omission verified

User explicitly chose the published DMSP convention: omit NaN channel-error
terms. This supersedes the all-19-errors rule and proposed alternatives.
conductance_from_spectrum retains unchanged P/H and optional analytical STD/
covariance. New uncertainty_channel_count int16 counts finite nonnegative
supplied errors, including0, even if nominal/propagation invalid. NaN/masked
terms are omitted; allmissing remains unavailable. Negative/Inf supplied
errors and singular zero-ionization derivatives remain invalid. Independent
channels and fixed atmosphere/B/energies remain explicit; no additional
common calibration/model covariance is added.

63 physics tests pass, including independent finite-difference/MC checks and
partial-NaN covariance, allmissing/definedzero, masked/batched/reversed cases.
Builder's135 focused tests pass. Real Jan1 output now has19,696 uncertainty-
valid of20,316 nominal-valid; all870 selected real CS test pairs resample.
Nominal outputs exactly match older scratch results. Detailed source and
execution evidence is owned by icBuilder's DMSP Spectrum Conductance Validation.

Builder carries channel counts, STD/covariance and separate uncertainty flags
through yearly and both crossing files. Confirmed workflow: regenerate yearly
files locally with updated packages, transfer them to the server and rebuild
crossings there, then rerun CS analysis. Annotations reusable.
No full jobs/production writes/commit/tag/merge performed. Preserve existing
worktree changes; server jobs are not to be monitored by the agent.

## Portfolio impact

- Central update needed: No
- Changes: None
- Sync summary: User-selected NaN omission is implemented and verified across
  icPhysics and icBuilder; data regeneration remains a server task.

## Previous checkpoint: spectrum conductance supports batched neutral profiles

`conductance_from_spectrum` now returns electron P/H from full spectra,
explicit atmosphere and reference B, with no I/O or new package dependencies.
The scientific equations preserve the user-provided IDL calculation. Inputs
have explicit units and broadcast shapes; missing/negative spectra are
invalid. Zero neutral species use a geometric-limit interpolation rather than
undefined log(0), enabling public NRLMSISE-00 profiles down to 70 km.
All 43 tests pass, including independent scalar parity on supplied synthetic
atmospheres, amplitude scaling, reverse ordering, missing values and zero flux.
Atmosphere fields accept shared 1D or batched (..., altitude) arrays on a shared
altitude grid; their leading dimensions broadcast with spectra and B. Shared
profiles preserve the fixed-kernel fast path. Batched results match individual
profile calls exactly and independent scalar references within 1.56e-15.
Original-XDR/IDL-runtime parity remains unavailable, not claimed.

icBuilder's yearly DMSP processor calls this function with user-selected
2° latitude × 5° longitude × 1-hour shared MSIS profiles, driven by local
historical daily indices. Per-measurement MSIS remains a comparison option.
It retains an explicit
fixed-profile alternative and evaluates IGRF at mapped 130-km footprints by
default. Its 50 index/integration tests and full daily CDF scratch conversion
passed. No production yearly corpus was changed.

Next: review/install this checkout before selected yearly regeneration with
local historical indices. User may investigate original MSIS table provenance.
The existing measurement-uncertainty merge/deployment remains unresolved;
this feature branch was created from its `eecc263` checkpoint.

## Portfolio impact

- Central update needed: No
- Changes: None
- Sync summary: Added user-approved spectrum conductance for DMSP validation.
  The spectrum calculation now also propagates measurement errors through
  icBuilder yearly/crossing/analysis products. This does not change the IMAGE
  production method, priority or deadline.

## Previous checkpoint: explicit detector measurement uncertainty

`proton_correct_images()` now exposes two explicit contracts. The default
`legacy_add_poisson` preserves existing callers and the frozen legacy result.
The new `measurement` mode assumes supplied standard deviations already
contain detector counting noise and therefore propagates them without adding
WIC/SI12/SI13 signal as a second Poisson term. It still combines WIC/SI13
measurement uncertainty with SI12-derived proton-correction uncertainty.

This fixes icBuilder detector Product 2, whose fuvpy Product-1 variance is
already full conditional detector measurement variance. The negative-signal
regression and all 17 icPhysics tests pass. A full orbit-0085 integration
comparison leaves Fp/E0/Fe and P/H unchanged and restores finite propagated
uncertainty on every method-valid cell.

## Project state

`icPhysics` now contains the crude shared calculation chain: Robinson,
camera-response tables, SI12 proton correction, the R method, WIC-derived
electron energy flux and uncertainty, E0--Fe covariance, Zhang--Paxton
collapse, the bundled 8.4-MB lookup table, loader, and generator. Two public
functions match the `PrecipitationImage` input/output boundary.

SI12 proton correction is now the public `proton_correct_images()` step.
`precipitation_from_ratio()` consumes corrected WIC/SI13 arrays and
`precipitation_from_zhang_paxton()` consumes corrected WIC plus Kp and MLT.
icBuilder calls correction once before method selection. The separated ratio
path matches the legacy combined calculation in a frozen reference case.

The public array-based `robinson_conductance()` entry point now supports the
modular icBuilder Product-3 boundary and preserves the scalar equations and
zero-flux uncertainty behavior in focused tests.

Hardy is now the default proton-energy model in icBuilder Product 2. The
public `PROTON_RESPONSE_ENERGY_RANGE` exposes the Frey table domain
(0.47--46.7 keV); icBuilder retains raw Hardy energy and clips only the array
passed to the camera responses. SI12 remains the proton-flux source. The
full icPhysics test suite passes (15 tests), and an
isolated real orbit-0364 Product-2/Product-3 round trip passed.

The copied table has the same SHA-256 as icBuilder. A small ratio-method
comparison matched the original scalar routine exactly. The new and old
Zhang--Paxton paths also matched exactly for every saved Product-2 field over
complete example orbit 0085. The two icPhysics Robinson tests and nine focused
icBuilder Product-2 tests pass.

Hardy et al. (1991) Table 1 is transcribed in
`src/icphysics/hardy_coefficients.py`. It contains both integral number-flux
and energy-flux coefficients with shape `(7, 13, 6)`, indexed by Kp level,
Fourier term, and generalized-Epstein parameter. All 1092 scalar values were
visually compared with the rendered journal table; no ambiguous cells remain.

`src/icphysics/hardy.py` implements equations 5--11 and exposes
`hardy_ion_precipitation(kp, mlt, mlat)`. Fractional Kp is handled by
interpolating evaluated log-flux profiles, never coefficients. The model
returns native integral number flux, integral energy flux, and their ratio as
mean ion energy in keV. Its Kp 0, 2, and 4 map morphology agrees visually with
the published model panels, and four independent checkpoints agree within 12%
with mean energies tabulated from the original 1989 statistical data. The full
icPhysics suite passes 14 tests.

The dependency boundary is confirmed: icAnalyzer must not depend on the full
icBuilder package and inherit fuvpy or other production-pipeline dependencies.
icPhysics depends on neither consumer. The separate `ZhangPaxton2008` package
remains authoritative for the published model equations.

## Next action

Main integration is complete. Update the UiB server installation before
regenerating detector Products 2 and 3. The deferred icAnalyzer batched-array
import check remains separate.

Do not start the staged file-product redesign or VAE integration in this
repository. Those consuming-project changes follow only after the shared
calculation is stable.

Repeat the complete-corpus image-ratio diagnostics with Hardy proton energy and
retain the Modified-Apex coordinate approximation in the interpretation.

## Risks and open decisions

- The existing functions may encode undocumented clipping, singular-boundary,
  mask, or uncertainty behavior that must be preserved during extraction.
- Zhang--Paxton mean-energy versus characteristic-energy compatibility remains
  a scientific gate outside simple code movement.
- The probabilistic E0 interpretation remains undecided, but the lookup and
  collapse belong in icPhysics.
- Adding a generalized model interface before a second forward model exists
  would create unnecessary framework code.
- The Hardy evaluator is connected to the IMAGE proton-correction pipeline,
  but the complete corpus has not yet been regenerated with it.
- Hardy is defined in corrected geomagnetic latitude and MLT, whereas the
  IMAGE pipeline currently supplies Modified Apex coordinates at a 130-km
  reference height. Coumans et al. (2004) made the same comparison and cited a
  maximum latitude difference of 0.17 degrees near 68 degrees MLAT for Apex at
  110 km versus ground-level corrected geomagnetic coordinates, below the
  IMAGE-FUV resolution. Keep the coordinate convention explicit and perform a
  small grid-wide sensitivity check before pipeline integration.

## Portfolio impact

- Central update needed: Yes
- Changes: the shared proton-correction API now distinguishes legacy residual
  uncertainty from full detector measurement uncertainty.
- Sync summary: icBuilder detector Products 2 and 3 must be regenerated after
  the icPhysics update; Product 1 remains valid.

## Entry points

- `AGENTS.md`
- `README.md`
- `vault/01_Project/Project Brief.md`
- `vault/01_Project/Current State.md`
- `vault/02_Algorithm/Dependency and Model Boundary.md`
- `vault/03_Decisions/Decision Log.md`
