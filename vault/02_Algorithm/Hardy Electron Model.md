# Hardy Electron Model

Last reviewed: 2026-10-11
Status: Optional numerical API implemented and independently audited.

## Sources and energy definition

Hardy, Gussenhoven and Holeman (1985), DOI 10.1029/JA090iA05p04229,
printed p. 4231: mean energy is integral energy flux / integral number flux.
Hardy et al. (1987), DOI 10.1029/JA092iA11p12275, printed pp. 12284–12285:
equations 4–5 and tables 4–5 give the electron functional representation.
User-supplied PDFs are archived in sibling icBuilder `literature/`:
`hardy_et_al_1985_A_statistical_model_of_auroral_electron_precipitation.pdf`
and `hardy_et_al_1987_Statistical_and_functional_representations.pdf`.

Number flux units: electrons cm^-2 s^-1 sr^-1. Energy flux units:
keV cm^-2 s^-1 sr^-1. Their ratio is mean energy in keV. Measurements span
50eV–20keV, but the published full-spectrum integrals extrapolate to 100 keV.
The 1987 conductivity fits instead use moments above 661 eV; applying Robinson
to this full-spectrum mean energy does not reproduce Hardy's conductivity fits.

## Equations and explicit conventions

Each flux has364 coefficients:7 Kp states x13 Fourier terms x4 Epstein
parameters(r,h0,S1,S2). Coefficient axes are Kp, constant/cos1..6/sin1..6,
parameter. This differs from the six-parameter ion model.

For x=latitude-h0 and q=-S1/S2, equation 4 is evaluated as:

```text
log10(flux) = r + S1*x + (S2-S1)*[ln(1+q*exp(x)) - ln(1+q)]
```

`logaddexp` avoids overflow in the transition. Valid slopes require S1>0,
S2<0. Below-profile backgrounds use log10 floors 6 equatorward and 7 poleward
of h0 for both number and energy flux, as stated in the paper.

The harmonic phase is n*pi*MLT/12. Printed equation 5 omits pi; the correction
is inferred from the required 24 h Fourier periodicity and is documented in
code and provenance. The electron paper does not specify fractional-Kp treatment:
we interpolate completed log10 fluxes after floors, clipping finite Kp 0–6.
The highest source group is Kp>=6-. This is an implementation convention.

## API and package boundary

`hardy_electron_precipitation(kp, mlt, mlat)` broadcasts inputs and returns
number_flux, energy_flux, mean_energy. MLT wraps 24 h. Northern corrected
geomagnetic latitude 50–90° is supported; outside-domain, nonfinite or masked
inputs return NaN. There is no hemisphere reflection or coordinate conversion.
icBuilder's Modified Apex 130 km approximation is recorded by that consumer.

`precipitation_from_hardy(wic_corrected, dwic_corrected, kp, mlt, mlat)`
retains raw Hardy E0 and uses corrected WIC / WIC response(E0) for event Fe.
Nonpositive corrected WIC is unavailable. Response endpoints hold outside
0.2–25keV; signed int8 flag codes are -1 unavailable, 0 within, 1 below, 2 above.
Missing/invalid measurement STD preserves nominal outputs but makes dFe
unavailable. dE0=0 and covariance=0 express fixed-energy conditional
measurement propagation, excluding Hardy prediction error.

No I/O, Kp acquisition, coordinate dependencies, IMAGE processing or particle
calibration were added. The existing ion evaluator remains separate.

## Verification and limits

All 728 table scalars were visually checked against independent PDF renders.
One energy coefficient at zero-based [2, 4, 2] was corrected to +0.0582433 during
the audit. Independent scalar reconstruction over10,175 points agrees within
2.66e-14 (number), 3.06e-14 (energy), 4.33e-14 (mean) relative. Memory optimization
retains only 4 parameters per pixel; an independent 10,000-point recheck
agrees within 8.44e-15 relative. A 1024 × 1024 numerical grid took~1.69 s locally.

Full package suite 72 passes, including boundaries, periodicity, fractional Kp,
broadcast/masked inputs, reference values and WIC response/error semantics.
Sibling builder completed a real 122-frame orbit through detector and CS
products and added an uncalibrated Hardy particle-comparison option.

These are numerical/integration gates. Spatial empirical validation,
corrected-geomagnetic versus Modified Apex sensitivity, out-of-response-domain
behavior and predictive uncertainty remain open. No default model changed.
