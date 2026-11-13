"""Reference check for the separated IMAGE precipitation steps."""

import numpy as np

import icphysics.image as image
from icphysics.image import E0_eflux_propagated, proton_response
from icphysics import (
    PROTON_RESPONSE_ENERGY_RANGE,
    precipitation_from_ratio,
    proton_correct_images,
)


def legacy_scalar_ratio(wic, dwic, si13, dsi13):
    """Frozen reference for the former per-pixel implementation."""

    arrays = [
        np.asarray(values, dtype=float)
        for values in (wic, dwic, si13, dsi13)
    ]
    shape = arrays[0].shape
    E0, dE0, Fe, dFe, R, dR = [
        np.full(shape, np.nan) for _ in range(6)
    ]

    for index in np.ndindex(shape):
        Wprime, dWprime, Sprime, dSprime = (
            array[index] for array in arrays
        )
        if not np.all(np.isfinite([Wprime, dWprime, Sprime, dSprime])):
            continue

        if Sprime != 0 and Wprime != 0:
            R[index] = image.fR(Wprime, Sprime)
            dR[index] = image.fdR(
                Wprime, Sprime, dWprime, dSprime
            )

        if Sprime < 3 and Wprime < 50:
            E0[index], dE0[index] = 0.2, 0.08
        elif Sprime < 3:
            E0[index], dE0[index] = 1.0, 25.0
        else:
            E0[index] = np.clip(image.fE0m(R[index]), 0, 35)
            dE0[index] = np.abs(
                image.fdE0m_dR(R[index]) * dR[index]
            )

        Wm = image.fWm(E0[index])
        dWm = np.abs(image.fdWm_dE0(E0[index]) * dE0[index])
        Fe[index] = image.fFe(Wprime, Wm)
        dFe[index] = image.fdFe(
            Wprime, Wm, dWprime, dWm
        )

    valid = np.isfinite(E0) & np.isfinite(Fe)
    varE0Fe = np.full(shape, np.nan)
    varE0Fe[valid] = 0.0
    return {
        "E0": E0, "dE0": dE0, "Fe": Fe, "dFe": dFe,
        "varE0Fe": varE0Fe, "R": R, "dR": dR,
    }


def test_proton_response_energy_range_matches_frey_table():
    assert PROTON_RESPONSE_ENERGY_RANGE == (0.47, 46.7)


def test_separate_proton_correction_matches_legacy_ratio_calculation():
    counts = [500.0, 50.0, 20.0]
    uncertainties = [5.0, 2.0, 2.0]
    proton_energy = 2.0
    dproton_energy = 0.1

    corrected = proton_correct_images(
        wic=np.array([counts[0]]),
        dwic=np.array([uncertainties[0]]),
        si12=np.array([counts[1]]),
        dsi12=np.array([uncertainties[1]]),
        si13=np.array([counts[2]]),
        dsi13=np.array([uncertainties[2]]),
        proton_energy=proton_energy,
        proton_energy_uncertainty=dproton_energy,
    )
    separated = precipitation_from_ratio(
        corrected["wic_corrected"],
        corrected["dwic_corrected"],
        corrected["si13_corrected"],
        corrected["dsi13_corrected"],
    )

    legacy = E0_eflux_propagated(
        counts,
        [0.0, 0.0, 0.0],
        uncertainties,
        proton_energy,
        dproton_energy,
    )
    names = ["E0", "Fe", "dE0", "dFe", "R", "dR"]

    for name, expected in zip(names, legacy):
        np.testing.assert_allclose(separated[name][0], expected)


def test_vectorized_ratio_matches_former_scalar_calculation():
    rng = np.random.default_rng(42)
    shape = (2, 4, 5)
    wic = rng.uniform(-20, 600, shape)
    si13 = rng.uniform(-2, 30, shape)
    dwic = rng.uniform(0, 20, shape)
    dsi13 = rng.uniform(0, 5, shape)

    cases = [
        (500.0, 20.0, 5.0, 2.0),
        (49.0, 2.9, 3.0, 1.0),
        (50.0, 2.9, 3.0, 1.0),
        (0.0, 5.0, 3.0, 1.0),
        (50.0, 0.0, 3.0, 1.0),
        (0.0, 0.0, 0.0, 0.0),
        (-1.0, 4.0, 3.0, 1.0),
        (50.0, 3.0, 3.0, 1.0),
        (np.nan, 10.0, 3.0, 1.0),
        (10.0, np.inf, 3.0, 1.0),
    ]
    for flat_index, (Wprime, Sprime, dWprime, dSprime) in enumerate(cases):
        wic.flat[flat_index] = Wprime
        si13.flat[flat_index] = Sprime
        dwic.flat[flat_index] = dWprime
        dsi13.flat[flat_index] = dSprime

    expected = legacy_scalar_ratio(wic, dwic, si13, dsi13)
    actual = precipitation_from_ratio(wic, dwic, si13, dsi13)

    for name in expected:
        np.testing.assert_allclose(
            actual[name], expected[name], rtol=0, atol=0, equal_nan=True
        )


def test_measurement_uncertainty_does_not_add_counting_noise_twice():
    wic = np.array([-10.0])
    si12 = np.array([-5.0])
    si13 = np.array([-4.0])
    dwic = np.array([3.0])
    dsi12 = np.array([2.0])
    dsi13 = np.array([1.0])
    proton_energy = np.array([2.0])
    dproton_energy = np.array([0.0])

    corrected = proton_correct_images(
        wic=wic,
        dwic=dwic,
        si12=si12,
        dsi12=dsi12,
        si13=si13,
        dsi13=dsi13,
        proton_energy=proton_energy,
        proton_energy_uncertainty=dproton_energy,
        uncertainty_mode="measurement",
    )

    response = proton_response(proton_energy, dproton_energy)
    expected_dfp = dsi12 / response["Tmodel"]
    expected_dwic = np.sqrt(
        dwic**2 + (response["Cpwm"] * expected_dfp)**2
    )
    expected_dsi13 = np.sqrt(
        dsi13**2 + (response["Cp13m"] * expected_dfp)**2
    )

    np.testing.assert_allclose(corrected["dFp"], expected_dfp)
    np.testing.assert_allclose(corrected["dwic_corrected"], expected_dwic)
    np.testing.assert_allclose(
        corrected["dsi13_corrected"], expected_dsi13
    )
    assert np.isfinite(corrected["dwic_corrected"]).all()
    assert np.isfinite(corrected["dsi13_corrected"]).all()


def test_proton_correction_rejects_unknown_uncertainty_mode():
    with np.testing.assert_raises_regex(ValueError, "uncertainty_mode"):
        proton_correct_images(
            wic=np.array([1.0]),
            dwic=np.array([1.0]),
            si12=np.array([1.0]),
            dsi12=np.array([1.0]),
            proton_energy=np.array([2.0]),
            proton_energy_uncertainty=np.array([0.0]),
            uncertainty_mode="unknown",
        )
