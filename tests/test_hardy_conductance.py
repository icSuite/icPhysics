"""Published-source and independent-equation checks for Hardy conductance."""

import math

import numpy as np

from icphysics import (
    hall, hardy_conductance, hardy_robinson_effective_precipitation, ped,
)
from icphysics.hardy_conductance import _load_coefficients


def scalar_conductance(coefficients, kp, mlt, latitude):
    """Separate scalar transcription of Hardy equations 4--5 and floors."""

    parameters = []
    for column in range(4):
        value = float(coefficients[kp, 0, column])
        for harmonic in range(1, 7):
            phase = math.pi * harmonic * mlt / 12
            value += coefficients[kp, harmonic, column] * math.cos(phase)
            value += coefficients[kp, harmonic + 6, column] * math.sin(phase)
        parameters.append(value)

    peak, peak_latitude, equatorward_slope, poleward_slope = parameters
    offset = latitude - peak_latitude
    transition = math.log(
        (1 - equatorward_slope / (poleward_slope * math.exp(-offset)))
        / (1 - equatorward_slope / poleward_slope)
    )
    conductance = peak + equatorward_slope * offset
    conductance += (poleward_slope - equatorward_slope) * transition
    floor = 0 if latitude < peak_latitude else .55
    return max(conductance, floor)


def test_coefficient_tables_match_published_cells():
    hall = _load_coefficients("hall")
    pedersen = _load_coefficients("pedersen")
    for coefficients in (hall, pedersen):
        assert coefficients.shape == (7, 13, 4)
        assert coefficients.dtype == np.float64
        assert not coefficients.flags.writeable

    np.testing.assert_equal(hall[0, 0], [2.7583333, 70.9270833, .5430208, -.6710083])
    np.testing.assert_equal(hall[3, 0], [10.0956250, 68.5687500, 2.3899479, -1.8148187])
    np.testing.assert_equal(hall[6, 12], [.3858664, -.4469092, .1524355, -.2234162])
    np.testing.assert_equal(pedersen[0, 0], [1.9752083, 71.7708333, .5717104, -.2543417])
    np.testing.assert_equal(pedersen[6, 12], [.2929164, -.0907934, .0712223, -.1096790])


def test_frozen_source_checkpoints():
    locations = [(0, 23.75, 68), (2, 1.75, 66), (3, 1.25, 65), (4, .25, 66), (6, 1.25, 62)]
    expected = np.array([
        [3.4898474222977893, 2.7319522849276825],
        [11.886787916468622, 6.843438416968734],
        [16.509579999431963, 8.711797458210759],
        [16.10873700135969, 8.470996790000129],
        [23.33593670752146, 11.944518985837703]
    ])
    actual = [[hardy_conductance(*point)[name] for name in ("hall", "pedersen")] for point in locations]
    np.testing.assert_allclose(actual, expected, rtol=2e-14)


def test_literal_equation_reference_across_spatial_grid():
    for name in ("hall", "pedersen"):
        coefficients = _load_coefficients(name)
        for kp in range(7):
            for mlt in (0, 3.75, 6, 11.25, 17, 23.5):
                for latitude in (50, 62.5, 70, 77.25, 90):
                    actual = hardy_conductance(kp, mlt, latitude)[name]
                    expected = scalar_conductance(coefficients, kp, mlt, latitude)
                    np.testing.assert_allclose(actual, expected, rtol=5e-14, atol=1e-14)


def test_broadcast_fractional_kp_mlt_wrap_and_domain():
    kp = np.array([0, 2.25, 6])[:, None, None]
    mlt = np.array([0, 6, 12, 18])[None, :, None]
    latitude = np.array([50, 70, 90])[None, None, :]
    result = hardy_conductance(kp, mlt, latitude)
    assert result["hall"].shape == (3, 4, 3)

    for name in result:
        np.testing.assert_equal(hardy_conductance(-4, -6, 70)[name], hardy_conductance(0, 18, 70)[name])
        np.testing.assert_equal(hardy_conductance(9, 48, 70)[name], hardy_conductance(6, 0, 70)[name])
        lower = hardy_conductance(2, 6, 70)[name]
        upper = hardy_conductance(3, 6, 70)[name]
        middle = hardy_conductance(2.25, 6, 70)[name]
        np.testing.assert_allclose(middle, .75 * lower + .25 * upper)

    invalid = hardy_conductance([2, np.nan, np.inf, 2, 2, 2], 0, [49, 70, 70, 90.01, -70, np.nan])
    assert all(np.isnan(value).all() for value in invalid.values())
    masked = hardy_conductance(np.ma.array([2, 2], mask=[False, True]), 0, 70)
    assert np.isfinite(masked["hall"][0]) and np.isnan(masked["hall"][1])


def test_published_background_levels():
    for name in ("hall", "pedersen"):
        assert hardy_conductance(0, 0, 50)[name] == 0
        assert hardy_conductance(0, 0, 90)[name] == .55


def test_paper_tabulated_maxima_are_reproduced_at_bin_centres():
    # Tables 1--2 quote the smoothed source-map maxima, not the functional-fit
    # maxima, so agreement is approximate at the listed half-hour/bin centres.
    midnight = hardy_conductance(0, 23.75, 68)
    noon = hardy_conductance(0, 12, 73)
    np.testing.assert_allclose([midnight["pedersen"], midnight["hall"]], [3.09, 4.05], rtol=.16)
    np.testing.assert_allclose([noon["pedersen"], noon["hall"]], [2.15, 3.80], rtol=.03)


def test_effective_precipitation_round_trips_native_conductance():
    kp = np.array([0, 2.5, 6])
    mlt = np.array([23.75, 1.5, 5])
    mlat = np.array([68, 66, 64])
    result = hardy_robinson_effective_precipitation(kp, mlt, mlat)
    native = hardy_conductance(kp, mlt, mlat)

    assert result["valid"].all()
    np.testing.assert_array_equal(result["hall"], native["hall"])
    np.testing.assert_array_equal(result["pedersen"], native["pedersen"])
    np.testing.assert_allclose(
        hall(result["effective_energy"], result["effective_energy_flux"]),
        result["hall"],
        rtol=2e-15,
    )
    np.testing.assert_allclose(
        ped(result["effective_energy"], result["effective_energy_flux"]),
        result["pedersen"],
        rtol=2e-15,
    )


def test_effective_precipitation_rejects_published_backgrounds():
    result = hardy_robinson_effective_precipitation(0, 0, [50, 90])

    np.testing.assert_equal(result["valid"], [False, False])
    np.testing.assert_equal(result["hall_background"], [True, True])
    np.testing.assert_equal(result["pedersen_background"], [True, True])
    assert np.isnan(result["effective_energy"][0])
    assert np.isfinite(result["effective_energy"][1])
    np.testing.assert_allclose(result["effective_energy"][1], 2.5584991552262224)
    np.testing.assert_allclose(result["effective_energy_flux"][1], 0.014681508604020935)


def test_fractional_kp_rejects_background_from_either_endpoint():
    lower = hardy_robinson_effective_precipitation(0, 0, 61.525)
    upper = hardy_robinson_effective_precipitation(1, 0, 61.525)
    middle = hardy_robinson_effective_precipitation(.25, 0, 61.525)

    assert lower["pedersen_background"]
    assert not upper["pedersen_background"]
    assert middle["pedersen_background"]
    assert np.isfinite(middle["effective_energy"])
    assert not middle["valid"]


def test_effective_energy_floor_invalidates_without_clipping():
    result = hardy_robinson_effective_precipitation(0, 5.25, 63.9)

    assert not result["hall_background"]
    assert not result["pedersen_background"]
    assert 0 < result["effective_energy"] < .5
    assert np.isfinite(result["effective_energy_flux"])
    assert not result["valid"]


def test_effective_precipitation_broadcast_mask_and_invalid_inputs():
    result = hardy_robinson_effective_precipitation(
        np.array([0, 2])[:, None],
        np.array([0, 6, 12])[None, :],
        70,
    )
    assert result["effective_energy"].shape == (2, 3)
    assert result["valid"].dtype == np.bool_

    masked = hardy_robinson_effective_precipitation(
        np.ma.array([2, 2], mask=[False, True]), 0, 70
    )
    assert np.isfinite(masked["effective_energy"][0])
    assert np.isnan(masked["effective_energy"][1])
    assert not masked["valid"][1]

    invalid = hardy_robinson_effective_precipitation(
        [2, np.nan, np.inf, 2], 0, [49, 70, 70, 90.01]
    )
    assert np.isnan(invalid["effective_energy"]).all()
    assert not invalid["valid"].any()

    scalar = hardy_robinson_effective_precipitation(2, 1.5, 66)
    assert scalar["effective_energy"].shape == ()
    assert scalar["valid"].shape == ()
