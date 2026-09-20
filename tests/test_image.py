"""Reference check for the separated IMAGE precipitation steps."""

import numpy as np

from icphysics.image import E0_eflux_propagated, proton_response
from icphysics import (
    PROTON_RESPONSE_ENERGY_RANGE,
    precipitation_from_ratio,
    proton_correct_images,
)


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
