#%% Import 
import numpy as np
from scipy.interpolate import interp1d,UnivariateSpline

#%% Predefined things

PE = np.array([0.47,2.00,8.00,25.0,46.7])  # Proton characteristic energy (keV)
P1 = np.array([145,319,554,601,562])       # WIC counts in response to 1 mW/m² proton energy flux, Table V
P2 = np.array([23.7,35.6,30.2,17.0,11.7])  # SI-12 proton counts per 1 mW/m², nominal exposure of 5 s, Table III
P3 = np.array([2.14,4.39,7.09,7.20,6.51])  # SI-13 counts in response to 1 mW/m² proton energy flux, Table VIII

# The camera-response tables contain no information outside this energy range.
# Calling pipelines can preserve a model's raw energy while clipping the value
# supplied to these response curves.
PROTON_RESPONSE_ENERGY_RANGE = (float(PE[0]), float(PE[-1]))

fP1 = interp1d(PE,P1,kind='linear')
fP2 = interp1d(PE,P2,kind='linear')
fP3 = interp1d(PE,P3,kind='linear')

# Electron conversion values from Frey et al. (2003)
EE = np.array([0.2,0.5,1.0,5.0,10.0,25.0])     # Electron characteristic energy (keV)
E1 = np.array([446,470,511,377,223,101])       #WIC counts in response to 1 mW/m² e- energy flux, Table IV
E3 = np.array([12.8,11.3,8.75,4.26,2.11,0.74]) #SI-13 counts in response to 1 mW/m² e- energy flux, Table VII
wic_to_s13 = E1/E3
EE_ratio = np.copy(EE)

# fE1 = interp1d(EE,E1,fill_value=(446,101),bounds_error=False)
# fE3 = interp1d(EE,E3,fill_value=(12.8,0.74),bounds_error=False)
# fE1divE3 = interp1d(wic_to_s13,EE,fill_value=(0.2,25),bounds_error=False)

# Translation between Notability doc and this python code
# fWm = fE1                       # Model of WIC counts in response to 1mW/m² e- eflux
# E0m = fE1divE3                  # Model of mean electron energy in keV
# Cpwm = P1
# Cp13m = P3
# Tm = P2

extrapolate_and_scale_wic_to_s13_fraction = None # set to None to not use this. The number is a scaling of the ratio between wic to s13 from the Frey table
use_interp1d = True
use_spline = False

if not extrapolate_and_scale_wic_to_s13_fraction is None:
    x_, y_ = np.log10(EE_ratio), np.log10(wic_to_s13 * extrapolate_and_scale_wic_to_s13_fraction)
    p = np.polyfit(x_, y_, 1)

    # make new arrays for EE and wic_to_s13:
    EE_ratio = np.linspace(1e-5, 1e3, 1000)
    wic_to_s13 = 10**np.polyval(p, np.log10(EE_ratio))





if use_spline:
    splinekw = dict(k=3)

    if splinekw['k'] == 1:
        maxE0 = 37.4
    elif splinekw['k'] == 2:
        maxE0 = 35
    elif (splinekw['k'] == 3) or (splinekw['k'] is None):
        maxE0 = 26.19

## Model of WIC electron counts for 1mW/m² energy flux as function of electron char energy in keV
#fWm = UnivariateSpline(EE,E1,**splinekw)
# fdWm_dE0 = fWm.derivative(1)

# Limit range of electron char energy to [0,26.16] such that fWm > 0
if use_spline:
    fWm = lambda x: UnivariateSpline(EE,E1,**splinekw)(np.clip(x,0,maxE0))
    fdWm_dE0 = lambda x: UnivariateSpline(EE,E1,**splinekw).derivative(1)(np.clip(x,0,maxE0))

if use_interp1d:
    fWm = interp1d(EE,E1,fill_value=(446,101),bounds_error=False)
    fdWm_dE0 = lambda x: (fWm(x+0.01)-fWm(x-0.01))/0.02


## energy flux model
# wprime: W' = W - dg_w - Cpw
# w     : WIC counts
# dg_w  : dayglow model for WIC
# Cpw   : proton counts predicted from SI-12
fFe = lambda wprime, wm: wprime/wm

#%% Support functions

def fdFe(wprime, wm, dwprime, dwm):
    return np.sqrt( (dwprime/wm)**2 + (wprime*dwm/wm**2)**2 )


def e0_fe_covariance(E0, Fe, dE0):
    """Covariance induced because the WIC-derived Fe depends on E0.

    At fixed proton-corrected WIC counts, Fe = Wprime / Wm(E0). The
    first-order covariance is therefore

        -Wprime * Wm'(E0) / Wm(E0)**2 * dE0**2.
    """

    Wm = fWm(E0)
    Wprime = Fe * Wm
    return -Wprime * fdWm_dE0(E0) / Wm**2 * dE0**2


## electron characteristic energy model
if use_spline:
    #fE0m = UnivariateSpline(wic_to_s13*2,EE_ratio,**splinekw)
    fE0m = UnivariateSpline(wic_to_s13, EE_ratio, **splinekw)
    fdE0m_dR = fE0m.derivative(1)

if use_interp1d:
    #fE0m = interp1d(wic_to_s13*2,EE_ratio,fill_value=(0.2,25),bounds_error=False)
    fE0m = interp1d(wic_to_s13, EE_ratio, fill_value=(0.2,25), bounds_error=False)
    fdE0m_dR = lambda x: (fE0m(x+1)-fE0m(x-1))/2


##############################
## MODELS OF PROTON COUNTS IN EACH CAMERA FOR 1 mW/m² input

## model of proton counts in WIC for 1 mW/m² input
if use_spline:
    fCpwm = UnivariateSpline(PE,P1,**splinekw)
    fdCpwm_dEp = fCpwm.derivative(1)

if use_interp1d:
    # fCpwm = interp1d(PE,P1,kind='linear')
    # Construct each response interpolator once.  These response tables never
    # change during an orbit, so rebuilding a SciPy object for every image
    # pixel only adds overhead to the scientific calculation.
    _fCpwm_interp = interp1d(
        PE, P1, kind='linear', fill_value='extrapolate'
    )
    fCpwm = lambda x: np.clip(_fCpwm_interp(x), 0, None)
    fdCpwm_dEp = lambda x: (fCpwm(x+0.05)-fCpwm(x-0.05))/0.1


def fdCpwm(Ep,dEp):
    """
    Uncertainty in model of proton counts in WIC for 1 mW/m² input

    Ep  : Proton characteristic energy [keV]
    dEp : Uncertainty in Ep            [keV]
    """
    return np.abs( fdCpwm_dEp(Ep) * dEp  )

## model of proton counts in SI-13 for 1 mW/m² input
if use_spline:
    fCp13m = UnivariateSpline(PE,P3,**splinekw)
    fdCp13m_dEp = fCp13m.derivative(1)

if use_interp1d:
    # fCp13m = interp1d(PE,P3,kind='linear')
    _fCp13m_interp = interp1d(
        PE, P3, kind='linear', fill_value='extrapolate'
    )
    fCp13m = lambda x: np.clip(_fCp13m_interp(x), 0, None)
    fdCp13m_dEp = lambda x: (fCp13m(x+0.05)-fCp13m(x-0.05))/0.1

def fdCp13m(Ep,dEp):
    """
    Uncertainty in model of proton counts in SI-13 for 1 mW/m² input

    Ep  : Proton characteristic energy [keV]
    dEp : Uncertainty in Ep            [keV]
    """
    return np.abs( fdCp13m_dEp(Ep) * dEp  )

## model of proton counts in SI-12 for 1 mW/m² input
if use_spline:
    fTm = UnivariateSpline(PE,P2,**splinekw)
    fdTm_dEp = fTm.derivative(1)

if use_interp1d:
    #fTm = interp1d(PE,P2,kind='linear')
    _fTm_interp = interp1d(
        PE, P2, kind='linear', fill_value='extrapolate'
    )
    fTm = lambda x: np.clip(_fTm_interp(x), 0, None)

    fdTm_dEp = lambda x: (fTm(x+0.05)-fTm(x-0.05))/0.1

def fdTm(Ep,dEp):
    """
    Uncertainty in model of proton counts in SI-12 for 1 mW/m² input

    Ep  : Proton characteristic energy [keV]
    dEp : Uncertainty in Ep            [keV]
    """
    return np.abs( fdTm_dEp(Ep) * dEp  )


def proton_response(Ep, dEp):
    """Evaluate the six orbit-invariant proton response quantities.

    ``Ep`` and ``dEp`` may be scalars or broadcast-compatible arrays. The
    ordinary dictionary keeps the names and physical roles visible where the
    values are used by the scalar and vector calculations.
    """

    return {
        'Tmodel': fTm(Ep),
        'dTmodel': fdTm(Ep, dEp),
        'Cpwm': fCpwm(Ep),
        'dCpwm': fdCpwm(Ep, dEp),
        'Cp13m': fCp13m(Ep),
        'dCp13m': fdCp13m(Ep, dEp),
    }

##############################
## TOTAL PREDICTED PROTON COUNTS IN EACH CAMERA GIVEN PROTON ENERGY FLUX

## Proton energy flux and uncertainty
def fFp(Tprime, Tmodel):
    return Tprime/Tmodel

def fdFp(Tprime,Tmodel,dTprime,dTmodel=0):
    return np.sqrt( (dTprime/Tmodel)**2 + (Tprime*dTmodel/Tmodel**2)**2 )

##Proton counts in WIC and uncertainty
def fproton_wic(Ep,Fp):
    """
    Ep : Proton characteristic energy [keV]
    Fp : proton energy flux           [mW/m²]
    """
    return fCpwm(Ep)*Fp

def fdproton_wic(Ep,Fp,dEp,dFp):
    return np.sqrt( (Fp*fdCpwm(Ep,dEp))**2 + (fCpwm(Ep) * dFp )**2  ) 


##Proton counts in SI-13 and uncertainty
def fproton_si13(Ep,Fp):
    """
    Ep : Proton characteristic energy [keV]
    Fp : proton energy flux           [mW/m²]
    """
    return fCp13m(Ep)*Fp

def fdproton_si13(Ep,Fp,dEp,dFp):
    return np.sqrt( (Fp*fdCp13m(Ep,dEp))**2 + (fCp13m(Ep) * dFp )**2  ) 


##############################
## CORRECTED COUNTS IN EACH CAMERA (counts - dayglow model - modeled proton counts)

## W' and uncertainty
def fWprime(W,dayglow_wic, proton_wic,clip=True):
    Wprime = W-dayglow_wic-proton_wic
    if clip:
        Wprime = np.clip(Wprime,0.,None)
    return Wprime

def fdWprime(W,dproton_wic,
             ddayglow_wic=0):
    """
    W            : WIC counts (uncorrected)
    dproton_wic  : Uncertainty in proton counts in WIC
    ddayglow_wic : Uncertainty in dayglow model of WIC counts
    """
    return np.sqrt(W + ddayglow_wic**2 + dproton_wic**2)


## S' and uncertainty
def fSprime(S,dayglow_si13, proton_si13,clip=True):
    Sprime = S-dayglow_si13-proton_si13
    if clip:
        Sprime = np.clip(Sprime,0.,None)
    return Sprime

def fdSprime(S,dproton_si13,
             ddayglow_si13=0):
    """
    S             : SI-13 counts (uncorrected)
    dproton_si13  : Uncertainty in SI-13 proton counts 
    ddayglow_si13 : Uncertainty in dayglow model of SI-13 counts
    """
    return np.sqrt(S + ddayglow_si13**2 + dproton_si13**2)


##T' and uncertainty
def fTprime(T,dayglow_T):
    return T-dayglow_T
def fdTprime(T,ddayglow_T=0):
    return np.sqrt(T + ddayglow_T**2)


## R=W'/S' uncertainty
def fR(Wprime,Sprime):
    return Wprime/Sprime

def fdR(Wprime,Sprime,dWprime,dSprime):
    """
    Wprime : Corrected WIC counts (W' = W - dayglow_WIC - proton_WIC)
    Sprime : Correct SI-13 counts (S' = S - dayglow_SI13 - proton_SI13)
    """

    return np.sqrt( (dWprime / Sprime)**2 + (Wprime*dSprime/Sprime**2)**2 )


#%% Final form

def E0_eflux_propagated(counts_list,
                        dayglowcounts_list,
                        dayglowcounts_unc,
                        Ep,dEp,
                        clip=True,
                        E0=None,
                        dE0=None,
                        proton_response_values=None):
    """
    E0, Fe, dE0, dFe = E0_eflux_propagated(counts_list,
                                           dayglowcounts_list,
                                           dayglowcounts_unc,
                                           Ep,dEp)

    Get predicted e- characteristic energy (E0), e- energy flux (Fe),
    and their uncertainties given counts from WIC, SI-12, and SI-13,
    their dayglow models, and proton characteristic energy and uncertainty

    Inputs
    ======
    counts_list        : List of counts in WIC, SI-12, and SI-13
    dayglowcounts_list : List of dayglow model counts for WIC, SI-12, and SI-13
    dayglowcounts_unc  : List of dayglow model count uncertainties for WIC, SI-12, and SI-13
    Ep                 : Proton characteristic energy [keV]
    dEp                : Uncertainty in proton characteristic energy [keV]
    E0, dE0            : Optional paired electron energy override [keV].
                         When supplied, the WIC/SI13 ratio remains diagnostic
                         but is not used to estimate electron energy.
    proton_response_values : dict, optional
                         Six response quantities returned by
                         ``proton_response``. ``ConductanceImage`` supplies
                         this cache because Ep and dEp are orbit-wide scalars.

    Outputs
    ======
    E0, Fe, dE0, dFe
    """

    if (E0 is None) != (dE0 is None):
        raise ValueError("E0 and dE0 must be supplied together")
    energy_override = E0 is not None

    W, T, S = counts_list
    dayglow_wic, dayglow_T, dayglow_si13 = dayglowcounts_list
    ddayglow_wic, ddayglow_T, ddayglow_si13 = dayglowcounts_unc

    ####################
    # 1. To isolate e- counts in WIC and SI-13, we need to estimate proton counts in WIC and SI-13 so that we can subtract them

    # Get corrected counts for SI-12
    Tprime = T-dayglow_T
    if clip:
        Tprime = np.clip(Tprime,0,None)
    dTprime = np.sqrt(T + ddayglow_T**2)

    # Get the camera responses to a 1 mW/m² proton energy flux.  Standalone
    # callers may omit the cache; ConductanceImage calculates it once because
    # Ep and dEp are constant across the whole orbit.
    if proton_response_values is None:
        proton_response_values = proton_response(Ep, dEp)
    Tmodel = proton_response_values['Tmodel']
    dTmodel = proton_response_values['dTmodel']
    Cpwm = proton_response_values['Cpwm']
    dCpwm = proton_response_values['dCpwm']
    Cp13m = proton_response_values['Cp13m']
    dCp13m = proton_response_values['dCp13m']

    # Estimate proton energy flux
    Fp = Tprime/Tmodel
    dFp = fdFp(Tprime,Tmodel,dTprime,dTmodel=dTmodel)

    # Get predicted proton counts for WIC and SI-13 given 1 mW/m² proton energy flux and estimated proton energy flux Fp

    proton_wic = Cpwm * Fp    # Denoted "Cpw" in Notability
    dproton_wic = np.sqrt((Fp*dCpwm)**2 + (Cpwm*dFp)**2)

    proton_si13 = Cp13m * Fp  # Denoted "Cp13" in Notability
    dproton_si13 = np.sqrt((Fp*dCp13m)**2 + (Cp13m*dFp)**2)

    ####################
    # 2. Isolate e- counts: Subtract dayglow and proton counts from WIC and SI-13

    # Get corrected counts for WIC, SI-13, and SI-12
    Wprime = fWprime(W,dayglow_wic, proton_wic, clip=clip)
    dWprime = fdWprime(W,dproton_wic,ddayglow_wic=ddayglow_wic)

    Sprime = fSprime(S,dayglow_si13, proton_si13,clip=clip)
    dSprime = fdSprime(S,dproton_si13,ddayglow_si13=ddayglow_si13)

    # Ratio of corrected WIC and SI-13 counts, and uncertainty in ratio
    if Sprime == 0 or Wprime == 0:
        R, dR = np.nan, np.nan
    else:
        R = Wprime/Sprime
        #if R > 150: # Cap to max ratio in Figure 7 of Frey et al (2003)
            #R = 150
            #Wprime = R * Sprime # test
        
        dR = fdR(Wprime,Sprime,dWprime,dSprime)
    
    ####################
    # 3. Estimate electron energy unless it was supplied by another model.
    if not energy_override:
        if Sprime < 3 and Wprime < 50:
            E0, dE0 = .2, .08
        elif Sprime < 3:
            E0, dE0 = 1, 25
        else:
            E0 = np.clip(fE0m(R),0,35)
            dE0 = np.abs(fdE0m_dR(R) * dR)

    # Get predicted electron counts for WIC given E0 and assuming 1 mW/m² e- energy flux
    Wm = fWm(E0)
    dWm = np.abs( fdWm_dE0(E0) * dE0 )

    # Estimate e- energy flux
    Fe = Wprime/Wm
    dFe = fdFe(Wprime, Wm, dWprime, dWm)
        
    return E0, Fe, dE0, dFe, R, dR


#%% Array entry points used by icBuilder

PROTON_CORRECTION_UNCERTAINTY_MODES = (
    "legacy_add_poisson",
    "measurement",
)


def proton_correct_images(
    wic,
    dwic,
    si12,
    dsi12,
    proton_energy,
    proton_energy_uncertainty,
    si13=None,
    dsi13=None,
    *,
    uncertainty_mode="legacy_add_poisson",
):
    """Apply SI12-based proton correction to image arrays.

    ``legacy_add_poisson`` preserves the extracted historical calculation:
    the supplied standard deviations exclude counting noise, so image signal
    is added as a Poisson variance term. ``measurement`` means the supplied
    standard deviations already describe full detector measurement variance;
    no additional image-signal variance is then added.
    """

    if uncertainty_mode not in PROTON_CORRECTION_UNCERTAINTY_MODES:
        raise ValueError(
            "uncertainty_mode must be 'legacy_add_poisson' or 'measurement'"
        )

    wic = np.asarray(wic, dtype=float)
    dwic = np.asarray(dwic, dtype=float)
    si12 = np.asarray(si12, dtype=float)
    dsi12 = np.asarray(dsi12, dtype=float)
    if not (wic.shape == dwic.shape == si12.shape == dsi12.shape):
        raise ValueError("WIC and SI12 values and uncertainties must have one shape")

    response = proton_response(proton_energy, proton_energy_uncertainty)

    si12_corrected = np.clip(si12, 0, None)
    if uncertainty_mode == "legacy_add_poisson":
        dsi12_corrected = np.sqrt(si12 + dsi12**2)
    else:
        dsi12_corrected = dsi12.copy()
    Fp = si12_corrected / response["Tmodel"]
    dFp = fdFp(
        si12_corrected,
        response["Tmodel"],
        dsi12_corrected,
        response["dTmodel"],
    )

    proton_wic = response["Cpwm"] * Fp
    dproton_wic = np.sqrt(
        (Fp * response["dCpwm"])**2
        + (response["Cpwm"] * dFp)**2
    )
    wic_corrected = np.clip(wic - proton_wic, 0, None)
    corrected_wic_variance = dwic**2 + dproton_wic**2
    if uncertainty_mode == "legacy_add_poisson":
        corrected_wic_variance = corrected_wic_variance + wic
    dwic_corrected = np.sqrt(corrected_wic_variance)

    result = {
        "wic_corrected": wic_corrected,
        "dwic_corrected": dwic_corrected,
        "Fp": Fp,
        "dFp": dFp,
    }

    if si13 is not None:
        si13 = np.asarray(si13, dtype=float)
        dsi13 = np.asarray(dsi13, dtype=float)
        if si13.shape != wic.shape or dsi13.shape != wic.shape:
            raise ValueError("SI13 values and uncertainties must match WIC")

        proton_si13 = response["Cp13m"] * Fp
        dproton_si13 = np.sqrt(
            (Fp * response["dCp13m"])**2
            + (response["Cp13m"] * dFp)**2
        )
        result["si13_corrected"] = np.clip(si13 - proton_si13, 0, None)
        corrected_si13_variance = dsi13**2 + dproton_si13**2
        if uncertainty_mode == "legacy_add_poisson":
            corrected_si13_variance = corrected_si13_variance + si13
        result["dsi13_corrected"] = np.sqrt(corrected_si13_variance)

    return result


def precipitation_from_ratio(
    wic_corrected,
    dwic_corrected,
    si13_corrected,
    dsi13_corrected,
):
    """Estimate precipitation from proton-corrected WIC and SI13 arrays."""

    arrays = [
        np.asarray(values, dtype=float)
        for values in (
            wic_corrected, dwic_corrected,
            si13_corrected, dsi13_corrected,
        )
    ]
    shape = arrays[0].shape
    if any(values.shape != shape for values in arrays):
        raise ValueError("all image arrays must have one shape")

    E0 = np.full(shape, np.nan)
    dE0 = np.full(shape, np.nan)
    Fe = np.full(shape, np.nan)
    dFe = np.full(shape, np.nan)
    R = np.full(shape, np.nan)
    dR = np.full(shape, np.nan)

    for index in np.ndindex(shape):
        Wprime, dWprime, Sprime, dSprime = (
            array[index] for array in arrays
        )
        if not np.all(np.isfinite([Wprime, dWprime, Sprime, dSprime])):
            continue

        if Sprime != 0 and Wprime != 0:
            R[index] = fR(Wprime, Sprime)
            dR[index] = fdR(Wprime, Sprime, dWprime, dSprime)

        if Sprime < 3 and Wprime < 50:
            E0[index], dE0[index] = 0.2, 0.08
        elif Sprime < 3:
            E0[index], dE0[index] = 1.0, 25.0
        else:
            E0[index] = np.clip(fE0m(R[index]), 0, 35)
            dE0[index] = np.abs(fdE0m_dR(R[index]) * dR[index])

        Wm = fWm(E0[index])
        dWm = np.abs(fdWm_dE0(E0[index]) * dE0[index])
        Fe[index] = fFe(Wprime, Wm)
        dFe[index] = fdFe(Wprime, Wm, dWprime, dWm)

    valid = np.isfinite(E0) & np.isfinite(Fe)
    varE0Fe = np.full(shape, np.nan)
    varE0Fe[valid] = 0.0

    return {
        "E0": E0,
        "dE0": dE0,
        "Fe": Fe,
        "dFe": dFe,
        "varE0Fe": varE0Fe,
        "R": R,
        "dR": dR,
    }


def precipitation_from_zhang_paxton(
    wic_corrected,
    dwic_corrected,
    kp,
    mlt,
):
    """Return Zhang--Paxton E0 and corrected-WIC electron energy flux."""

    from .zhang_paxton_lookup import load_zhang_paxton_lookup

    lookup = load_zhang_paxton_lookup(kp)

    mlt = np.asarray(mlt, dtype=float)
    mlt_difference = np.mod(mlt - lookup["mlt"] + 12, 24) - 12
    if mlt.shape != lookup["mlt"].shape or not np.allclose(mlt_difference, 0):
        raise ValueError("MLT does not match the bundled IMAGE lookup")

    E0 = np.asarray(lookup["E0"], dtype=float)
    dE0 = np.asarray(lookup["dE0"], dtype=float)
    if E0.shape != np.asarray(wic_corrected).shape:
        raise ValueError("Kp and image dimensions do not match")

    Wm = fWm(E0)
    dWm = np.abs(fdWm_dE0(E0) * dE0)
    Wprime = np.asarray(wic_corrected, dtype=float)
    dWprime = np.asarray(dwic_corrected, dtype=float)
    Fe = Wprime / Wm
    dFe = fdFe(Wprime, Wm, dWprime, dWm)
    varE0Fe = e0_fe_covariance(E0, Fe, dE0)

    valid = (
        np.isfinite(Wprime) & np.isfinite(dWprime)
    )
    for values in (E0, dE0, Fe, dFe, varE0Fe):
        values[~valid] = np.nan

    return {
        "E0": E0,
        "dE0": dE0,
        "Fe": Fe,
        "dFe": dFe,
        "varE0Fe": varE0Fe,
    }


#%% Full spatial Hardy energy with observed WIC flux

ELECTRON_RESPONSE_ENERGY_RANGE = (float(EE[0]), float(EE[-1]))


def precipitation_from_hardy(wic_corrected, dwic_corrected, kp, mlt, mlat):
    """Return Hardy mean energy and corrected-WIC event energy flux.

    All inputs broadcast. Corrected WIC counts set Fe [mW m^-2], while the
    full spatial Hardy model sets E0 [keV]. SI13 is not used. Raw E0 is
    retained; the existing WIC response holds endpoint values outside its
    0.2--25 keV table. E0_response_range_flag is int8: 0 within, 1 below,
    2 above, -1 unavailable. A missing/negative/nonfinite measurement STD
    makes dFe unavailable without discarding nominal E0/Fe. dE0 and
    varE0Fe are 0 on nominal-valid cells, conditional on fixed Hardy energy:
    they exclude statistical model error rather than assert exact prediction.
    Nonpositive corrected WIC is invalid for this precipitation route.
    """

    from .hardy_electron import hardy_electron_precipitation

    Wprime, dWprime, kp, mlt, mlat = np.broadcast_arrays(*[
        np.asarray(np.ma.filled(np.ma.asarray(value, dtype=float), np.nan))
        for value in (wic_corrected, dwic_corrected, kp, mlt, mlat)
    ])
    E0 = hardy_electron_precipitation(kp, mlt, mlat)["mean_energy"]
    Wm = fWm(E0)
    valid = np.isfinite(E0) & np.isfinite(Wprime) & (Wprime > 0)
    Fe = np.where(valid, Wprime / Wm, np.nan)
    E0 = np.where(valid, E0, np.nan)

    uncertainty_valid = valid & np.isfinite(dWprime) & (dWprime >= 0)
    dFe = np.where(uncertainty_valid, dWprime / Wm, np.nan)
    dE0 = np.where(valid, 0.0, np.nan)
    varE0Fe = np.where(valid, 0.0, np.nan)

    flag = np.full(E0.shape, -1, dtype=np.int8)
    flag = np.where(valid, 0, flag)
    flag = np.where(valid & (E0 < EE[0]), 1, flag)
    flag = np.where(valid & (E0 > EE[-1]), 2, flag)
    return {
        "E0": E0, "dE0": dE0, "Fe": Fe, "dFe": dFe,
        "varE0Fe": varE0Fe, "E0_response_range_flag": flag.astype(np.int8)
    }
