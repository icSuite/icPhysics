"""Tests for collapsed Zhang--Paxton precipitation with observed WIC."""

import numpy as np

from icphysics import precipitation_from_collapsed_zhang_paxton
from icphysics.image import e0_fe_covariance, fdFe, fdWm_dE0, fWm
import icphysics.zhang_paxton as zhang_paxton_module


def test_collapsed_zhang_paxton_propagates_structural_spread(monkeypatch):
    energy = np.array([2.0, 4.0, 3.0])
    spread = np.array([0.3, 0.5, 0.4])

    def fixed_model(kp, mlt):
        return {
            "energy": energy,
            "energy_spread": spread,
            "valid": np.array([True, True, False]),
        }

    monkeypatch.setattr(
        zhang_paxton_module,
        "collapsed_zhang_paxton_q_weighted",
        fixed_model,
    )
    wic = np.array([400.0, 800.0, 600.0])
    dwic = np.array([40.0, 80.0, 60.0])
    result = precipitation_from_collapsed_zhang_paxton(
        wic, dwic, 2.0, np.array([0.0, 1.0, 2.0])
    )

    response = fWm(energy[:2])
    response_std = np.abs(fdWm_dE0(energy[:2]) * spread[:2])
    flux = wic[:2] / response
    np.testing.assert_allclose(result["E0"][:2], energy[:2])
    np.testing.assert_allclose(result["dE0"][:2], spread[:2])
    np.testing.assert_allclose(result["Fe"][:2], flux)
    np.testing.assert_allclose(
        result["dFe"][:2],
        fdFe(wic[:2], response, dwic[:2], response_std),
    )
    np.testing.assert_allclose(
        result["varE0Fe"][:2],
        e0_fe_covariance(energy[:2], flux, spread[:2]),
    )
    assert np.isnan(result["E0"][2])
    assert result["E0_response_range_flag"][2] == -1
