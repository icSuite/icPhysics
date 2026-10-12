# Decision Log

Last reviewed: 2026-10-12

## 2026-10-12 — Integrate spectrum and measurement uncertainty into main

**Confirmed by the user:** merge feat/spectrum-conductance into main. Its
ancestry includes fix/measurement-uncertainty, so one fast-forward integrates
the measurement contract, spectrum conductance/uncertainty and optional Hardy
electron model. Fetched remote refs, clean worktree and 72 passing package
tests verified the integration. Retain existing feature branches and deploy
main on the server; no production-model selection or corpus rerun is implied.

## 2026-10-11 — Optional full spatial Hardy electron energy

**Confirmed by the user:** implement the optional Hardy electron mean-energy
model using the supplied functional paper and each pixel's MLAT, MLT and Kp.
Electron/proton choices remain independent; no default replacement requested.

**Implementation choices:** correct the Fourier phase to 24 h periodicity;
interpolate evaluated log10 flux for fractional Kp; expose the Northern 50–90°
source domain and Modified Apex approximation. Retain raw E0 and flag endpoint
WIC response evaluation. Propagate measurement errors conditional on fixed
Hardy energy, explicitly excluding prediction error. These conventions are
documented implementation choices, not newly confirmed scientific policies.
Broader particle/coordinate validation precedes production model selection.


## 2026-10-06 — Follow published DMSP NaN-error omission

**Confirmed by the user:** follow the published method and ignore NaN channel-
error terms. This supersedes our strict complete-channel-error policy and the
unselected reversion/low-count/moment alternatives. Retain the full nominal
spectrum, propagate available finite nonnegative errors with P/H covariance,
and save channel availability count plus separate uncertainty validity.
All-missing errors remain unavailable; supplied invalid errors and singular
propagation remain invalid. Missing terms are omitted, not a recovered noise
model. Source/method evidence and end-to-end checks are owned by icBuilder's
DMSP Spectrum Conductance Validation note.

## 2026-10-06 — Propagate spectral measurement errors analytically

**Confirmed by the user:** Implement first-order spectrum-error propagation
and verify against independent finite differences and small-error Monte Carlo.

**Implementation:** Extend conductance_from_spectrum with optional absolute
differential_energy_flux_std, retaining the original P/H-only result without
errors and unchanged nominal values when supplied. Return P_std, H_std,
PH_covariance and uncertainty_valid. Independent channel errors are explicit;
atmosphere/B/energies are fixed and no additional common calibration covariance
is added beyond components already present in supplied channel STD. Missing/
negative/nonfinite errors or a singular zero-ionization derivative invalidate
errors only. Preserve exact zero-error and nondepositing-channel limits.

The source DMSP method leaves low-count errors undefined; this numerical API
does not impute missing errors or implement a censored-count model. Production
source acquisition/quality/provenance and cohort selection remain in icBuilder.

## 2026-10-05 — Own the spectrum-to-conductance numerical calculation

**Confirmed by the user:** Translate the supplied IDL spectrum conductance
routine into icPhysics and use it in icBuilder's yearly DMSP processing.

**Boundary:** icPhysics accepts explicit energy/flux arrays, atmosphere
profiles and magnetic-field magnitude, returning electron P/H. Public-model
atmosphere generation, IGRF evaluation, coordinate mapping, particle-data
reading and yearly output remain outside this numerical package.

**Batched-profile extension:** Accept shared or broadcast batched neutral
profiles on a common altitude axis so per-measurement MSIS can be evaluated
by icBuilder without changing the numerical equations. Preserve the shared
profile fast path and validate batched results against individual profiles
and independent scalar calculations.

**Verification limit:** Synthetic supplied-atmosphere scalar parity is tested;
the original atmosphere is unavailable. Inherited angular and equilibrium
assumptions require physical validation before these estimates serve as
absolute reference conductances.

## 2026-09-20 — Make proton-correction uncertainty semantics explicit

**Decision:** Keep `legacy_add_poisson` as the backward-compatible default for
`proton_correct_images()` and add an explicit `measurement` mode for callers
whose supplied variance already contains detector counting noise.

**Rationale:** fuvpy Product-1 variance is full conditional measurement
variance. Adding image signal again double-counts noise and can create invalid
square roots for negative background-subtracted counts, while silently
changing the default would alter legacy callers.

## 2026-08-11 — Keep the Zhang--Paxton application layer in icPhysics

**Decision:** `ZhangPaxton2008` retains the published base equations.
icPhysics owns the collapse, uncertainty/sampling, lookup generation, bundled
lookup table, and loader. icBuilder will consume those interfaces instead of
retaining Zhang--Paxton code or data.

**Rationale:** This keeps the validated base model separate while avoiding an
application implementation split across ZhangPaxton2008, icPhysics, and
icBuilder.

## 2026-08-10 — Create a lightweight shared physics package

**Decision:** Use the dedicated `icPhysics` repository for numerical
precipitation and conductance routines shared by icBuilder and icAnalyzer.
Neither consumer depends on the other.

**Rationale:** icBuilder requires pipeline packages such as fuvpy that are not
needed in icAnalyzer. A direct icAnalyzer dependency on icBuilder would couple
the VAE environment to the complete data-production stack. Copying equations
would instead create scientific drift.

## 2026-08-10 — Keep project-specific orchestration outside icPhysics

**Decision:** Keep image processing, coordinate conversion, file products,
orbit execution, Kp acquisition and frame matching, and VAE workflows in their
own projects. icPhysics accepts prepared arrays and numerical model inputs.

**Rationale:** The shared unit is the scientific calculation, not the complete
workflow. A narrow boundary keeps dependencies small and makes numerical
equivalence testable.

## 2026-08-10 — Retain ZhangPaxton2008 as the model authority

**Decision:** Do not copy the published Zhang--Paxton implementation into
icPhysics. Depend on `ZhangPaxton2008` and add only shared higher-level collapse,
lookup, or sampling behavior where it belongs.

**Rationale:** The focused repository already owns transcription, validation,
figures, and publication provenance for that model.
