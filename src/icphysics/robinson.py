"""Robinson et al. (1987) Hall and Pedersen conductance equations.

This is an intentionally direct extraction from icBuilder. Scientific cleanup
and API consolidation can follow after the consuming projects reproduce their
current results.
"""

#%% Imports

import numpy as np


#%% Conductance

def ped(E0, Fe):
    """Return Pedersen conductance in siemens.

    Parameters
    ----------
    E0 : float or array
        Electron average energy [keV].
    Fe : float or array
        Electron energy flux [mW/m2, equivalent to erg/cm2/s].
    """

    return 40 * E0 / (16 + E0**2) * np.sqrt(Fe)


def hall(E0, Fe):
    """Return Hall conductance in siemens."""

    return 18 * E0**1.85 / (16 + E0**2) * np.sqrt(Fe)


def robinson_effective_precipitation(hall_conductance, pedersen_conductance):
    """Invert positive finite Hall/Pedersen pairs through Robinson.

    Inputs broadcast normally. Masked, nonfinite, zero, or negative
    conductances are not invertible and return NaN effective energy [keV] and
    energy flux [mW m-2]. No clipping or model-specific validity criterion is
    applied.
    """

    hall_conductance, pedersen_conductance = np.broadcast_arrays(
        np.ma.asarray(hall_conductance, dtype=float).filled(np.nan),
        np.ma.asarray(pedersen_conductance, dtype=float).filled(np.nan)
    )
    invertible = (
        np.isfinite(hall_conductance)
        & np.isfinite(pedersen_conductance)
        & (hall_conductance > 0)
        & (pedersen_conductance > 0)
    )

    effective_energy = np.full(hall_conductance.shape, np.nan)
    effective_energy[invertible] = (
        hall_conductance[invertible]
        / (0.45 * pedersen_conductance[invertible])
    ) ** (1.0 / 0.85)

    effective_energy_flux = np.full(hall_conductance.shape, np.nan)
    effective_energy_flux[invertible] = (
        pedersen_conductance[invertible]
        * (16.0 + effective_energy[invertible]**2)
        / (40.0 * effective_energy[invertible])
    ) ** 2
    return {
        "effective_energy": effective_energy,
        "effective_energy_flux": effective_energy_flux,
        "invertible": invertible
    }


#%% Uncertainty

def peduncertainty(E0, Fe, dE0, dFe, varE0Fe):
    """Propagate E0 and Fe uncertainty into Pedersen conductance."""

    if Fe == 0:
        # The derivative is singular at zero flux. Use the same one-sided
        # excursion from Fe=0 to Fe=dFe as the current icBuilder code.
        dP = ped(E0, dFe)
    else:
        denominator = 16 + E0**2
        dP_dE0 = (
            40 / denominator - 80 * (E0 / denominator)**2
        ) * np.sqrt(Fe)
        dP_dFe = 40 * E0 / denominator / (2 * np.sqrt(Fe))
        dP = np.sqrt(
            dP_dE0**2 * dE0**2
            + dP_dFe**2 * dFe**2
            + 2 * dP_dE0 * dP_dFe * varE0Fe
        )

    return dP


def halluncertainty(E0, Fe, dE0, dFe, varE0Fe):
    """Propagate E0 and Fe uncertainty into Hall conductance."""

    if Fe == 0:
        dH = hall(E0, dFe)
    else:
        denominator = 16 + E0**2
        dH_dE0 = (
            18 * E0**0.85 / denominator
            * (1.85 - 2 * E0**2 / denominator)
            * np.sqrt(Fe)
        )
        dH_dFe = 9 * E0**1.85 / denominator / np.sqrt(Fe)
        dH = np.sqrt(
            dH_dE0**2 * dE0**2
            + dH_dFe**2 * dFe**2
            + 2 * dH_dE0 * dH_dFe * varE0Fe
        )

    return dH


def robinson_conductance(E0, Fe, dE0, dFe, varE0Fe):
    """Calculate Hall/Pedersen conductance and propagated uncertainty."""

    E0 = np.asarray(E0, dtype=float)
    Fe = np.asarray(Fe, dtype=float)
    dE0 = np.asarray(dE0, dtype=float)
    dFe = np.asarray(dFe, dtype=float)
    varE0Fe = np.asarray(varE0Fe, dtype=float)

    if not (E0.shape == Fe.shape == dE0.shape == dFe.shape == varE0Fe.shape):
        raise ValueError("all precipitation arrays must have one shape")

    P = ped(E0, Fe)
    H = hall(E0, Fe)
    dP = np.full(E0.shape, np.nan)
    dH = np.full(E0.shape, np.nan)

    valid = np.isfinite(E0 + Fe + dE0 + dFe + varE0Fe)
    zero_flux = valid & (Fe == 0)
    dP[zero_flux] = ped(E0[zero_flux], dFe[zero_flux])
    dH[zero_flux] = hall(E0[zero_flux], dFe[zero_flux])

    nonzero = valid & (Fe > 0)
    denominator = 16 + E0[nonzero]**2

    dP_dE0 = (
        40 / denominator - 80 * (E0[nonzero] / denominator)**2
    ) * np.sqrt(Fe[nonzero])
    dP_dFe = 40 * E0[nonzero] / denominator / (2 * np.sqrt(Fe[nonzero]))
    dP[nonzero] = np.sqrt(
        dP_dE0**2 * dE0[nonzero]**2
        + dP_dFe**2 * dFe[nonzero]**2
        + 2 * dP_dE0 * dP_dFe * varE0Fe[nonzero]
    )

    dH_dE0 = (
        18 * E0[nonzero]**0.85 / denominator
        * (1.85 - 2 * E0[nonzero]**2 / denominator)
        * np.sqrt(Fe[nonzero])
    )
    dH_dFe = 9 * E0[nonzero]**1.85 / denominator / np.sqrt(Fe[nonzero])
    dH[nonzero] = np.sqrt(
        dH_dE0**2 * dE0[nonzero]**2
        + dH_dFe**2 * dFe[nonzero]**2
        + 2 * dH_dE0 * dH_dFe * varE0Fe[nonzero]
    )

    return {"P": P, "H": H, "dP": dP, "dH": dH}
