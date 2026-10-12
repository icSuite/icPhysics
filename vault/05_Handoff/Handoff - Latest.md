# Handoff - Latest

Last updated: 2026-10-12
Repository snapshot: `main`; implementation checkpoint `8db326f`
Worktree state: implementation committed; integration closeout updates project memory.

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
