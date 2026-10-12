"""Electron model source checkpoints and independent literal-equation checks."""

import math

import numpy as np

from icphysics import hardy_electron_precipitation, precipitation_from_hardy
from icphysics.hardy_electron import _epstein_log_flux, _fourier_parameters
from icphysics.hardy_electron_coefficients import ENERGY_FLUX_COEFFICIENTS, NUMBER_FLUX_COEFFICIENTS


def scalar_flux(coefficients, kp, mlt, latitude):
    """Separate scalar transcription of printed equations 4/5 and backgrounds."""

    parameters = []
    for column in range(4):
        value = float(coefficients[kp, 0, column])
        for harmonic in range(1, 7):
            phase = math.pi * harmonic * mlt / 12
            value += coefficients[kp, harmonic, column] * math.cos(phase)
            value += coefficients[kp, harmonic + 6, column] * math.sin(phase)
        parameters.append(value)
    r, h0, s1, s2 = parameters
    offset = latitude - h0
    transition = math.log((1 - s1 / (s2 * math.exp(-offset))) / (1 - s1 / s2))
    profile = r + s1 * offset + (s2 - s1) * transition
    floor = 6 if latitude < h0 else 7
    return 10 ** max(profile, floor)


def test_coefficients_source_shape_and_visually_resolved_cells():
    for coefficients in (ENERGY_FLUX_COEFFICIENTS, NUMBER_FLUX_COEFFICIENTS):
        assert coefficients.shape == (7, 13, 4)
        assert coefficients.dtype == np.float64
        assert not coefficients.flags.writeable
    np.testing.assert_equal(ENERGY_FLUX_COEFFICIENTS[0, 0], [7.7896976, 73.0708333, .2191803, -.0859725])
    assert ENERGY_FLUX_COEFFICIENTS[2, 4, 2] == .0582433
    assert ENERGY_FLUX_COEFFICIENTS[2, 5, 2] == .0166424
    assert NUMBER_FLUX_COEFFICIENTS[3, 1, 2] == .005109
    np.testing.assert_equal(NUMBER_FLUX_COEFFICIENTS[6, 12], [.0203793, -.7063537, -.0000407, .0016903])


def test_frozen_independent_source_checkpoints():
    # Independently transcribed scalar equations following the full visual audit.
    locations = [(0,23.25,64), (1,9.75,70), (2,10.75,70), (3,10.25,69), (4,10.25,67), (6,9.75,65)]
    expected = np.array([
        [8440780.572585989,13575456.962490907,1.6083177196409235],
        [19264008.409364805,62003318.371242434,3.2186093908213196],
        [21021644.079028167,124989341.89885204,5.945745319869877],
        [26005353.676444996,150182550.4475617,5.775062793458308],
        [37281202.83530745,89807598.37708998,2.408924378696199],
        [105127661.49080186,93858021.8760111,.892800434681249]
    ])
    actual = [list(hardy_electron_precipitation(*point).values()) for point in locations]
    np.testing.assert_allclose(actual, expected, rtol=2e-14)


def test_literal_equation_reference_across_spatial_grid():
    for kp in range(7):
        for mlt in (0,3.75,6,11.25,17,23.5):
            for latitude in (50,62.5,70,77.25,90):
                actual = hardy_electron_precipitation(kp,mlt,latitude)
                for name, table in (('number_flux',NUMBER_FLUX_COEFFICIENTS),('energy_flux',ENERGY_FLUX_COEFFICIENTS)):
                    np.testing.assert_allclose(actual[name],scalar_flux(table,kp,mlt,latitude),rtol=5e-14)


def test_peak_and_asymptotic_slopes_are_stable():
    parameters = _fourier_parameters(ENERGY_FLUX_COEFFICIENTS,np.asarray(2),np.asarray(12))
    r,h0,s1,s2 = parameters
    np.testing.assert_allclose(_epstein_log_flux(parameters,h0),r,atol=1e-14)
    # Stable evaluator remains finite where the literal exponential would overflow.
    flux = _epstein_log_flux(parameters,np.array([-1e4,1e4]))
    assert np.all(np.isfinite(flux))
    low_slope = _epstein_log_flux(parameters,-1e4+1) - flux[0]
    high_slope = _epstein_log_flux(parameters,1e4+1) - flux[1]
    np.testing.assert_allclose([low_slope,high_slope],[s1,s2],atol=1e-12)


def test_broadcast_fractional_kp_mlt_wrap_and_domain():
    kp = np.array([0,2.25,6])[:,None,None]
    mlt = np.array([0,6,12,18])[None,:,None]
    latitude = np.array([50,70,90])[None,None,:]
    result = hardy_electron_precipitation(kp,mlt,latitude)
    assert result['mean_energy'].shape == (3,4,3)
    np.testing.assert_allclose(result['mean_energy'],result['energy_flux']/result['number_flux'])
    for name in result:
        np.testing.assert_equal(hardy_electron_precipitation(-4,-6,70)[name],hardy_electron_precipitation(0,18,70)[name])
        np.testing.assert_equal(hardy_electron_precipitation(9,48,70)[name],hardy_electron_precipitation(6,0,70)[name])
        lower = hardy_electron_precipitation(2,6,70)[name]
        upper = hardy_electron_precipitation(3,6,70)[name]
        middle = hardy_electron_precipitation(2.25,6,70)[name]
        np.testing.assert_allclose(np.log10(middle),.75*np.log10(lower)+.25*np.log10(upper))
    invalid = hardy_electron_precipitation([2,np.nan,np.inf,2,2,2],0,[49,70,70,90.01,-70,np.nan])
    assert all(np.isnan(value).all() for value in invalid.values())
    masked = hardy_electron_precipitation(np.ma.array([2,2],mask=[False,True]),0,70)
    assert np.isfinite(masked['mean_energy'][0]) and np.isnan(masked['mean_energy'][1])


def test_published_background_levels():
    for name in ('energy_flux','number_flux'):
        assert hardy_electron_precipitation(0,0,50)[name] == 1e6
        assert hardy_electron_precipitation(0,0,90)[name] == 1e7


def test_observed_wic_sets_flux_and_measurement_uncertainty():
    result = precipitation_from_hardy(np.array([1000,2000]),np.array([10,20]),2,10,70)
    e0 = hardy_electron_precipitation(2,10,70)['mean_energy']
    response = np.interp(e0,[.2,.5,1,5,10,25],[446,470,511,377,223,101])
    np.testing.assert_allclose(result['E0'],e0)
    np.testing.assert_allclose(result['Fe'],np.array([1000,2000])/response)
    np.testing.assert_allclose(result['dFe'],np.array([10,20])/response)
    np.testing.assert_equal(result['dE0'],0)
    np.testing.assert_equal(result['varE0Fe'],0)
    # Nominal values and fixed-energy conditional terms survive unavailable STD.
    missing = precipitation_from_hardy(1000,[np.nan,-1,np.inf,0],2,10,70)
    assert np.all(np.isfinite(missing['E0'])) and np.all(np.isfinite(missing['Fe']))
    assert np.isnan(missing['dFe'][:3]).all() and missing['dFe'][3] == 0
    np.testing.assert_equal(missing['dE0'],0)


def test_response_flags_preserve_raw_energy_and_endpoint_response(monkeypatch):
    import icphysics.hardy_electron as hardy
    monkeypatch.setattr(hardy,'hardy_electron_precipitation',lambda kp,mlt,mlat: {'mean_energy': np.asarray(kp)})
    result = precipitation_from_hardy(1000,20,[.1,.2,25,30,np.nan],0,70)
    np.testing.assert_equal(result['E0_response_range_flag'],[1,0,0,2,-1])
    assert result['E0_response_range_flag'].dtype == np.int8
    np.testing.assert_equal(result['E0'][:4],[.1,.2,25,30])
    np.testing.assert_allclose(result['Fe'][:4],1000/np.array([446,446,101,101]))
    np.testing.assert_allclose(result['dFe'][:4],20/np.array([446,446,101,101]))


def test_nonpositive_or_missing_wic_and_model_support_invalid():
    result = precipitation_from_hardy([0,-1,np.nan,1000,1000],10,2,0,[70,70,70,49,70])
    assert np.isnan(result['E0'][:4]).all()
    assert np.isnan(result['Fe'][:4]).all()
    assert np.isnan(result['dFe'][:4]).all()
    np.testing.assert_equal(result['E0_response_range_flag'],[-1,-1,-1,-1,0])
