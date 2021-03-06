# J. Gomes da Silva, 2019
# Institute of Astrophysics and Space Sciences, Porto, Portugal
# Joao.Silva@astro.up.pt

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division
import numpy as np
try:
    from astroquery.simbad import Simbad
except:
    print("*** WARNING: astroquery is not installed but required for the get_bv function")




def calc_smw(caii, caii_err, instr='ESPRESSO'):
    """Calculates S_MW (Mt. Wilson S-index, Vaughan et al. (1978)) index based on the CaII H&K lines.

    Parameters:
    -----------
        caii : float, list, array
            CaII index based on the CaII H&K lines.
        caii_err : float, list, array
            Photon error of the CaII index.
        instr : string
            Instrument used: 'ESPRESSO' (default) or 'HARPS'.
    Returns:
    --------
        smw : float, array
            S-index calibrated to the Mt. Wilson scale.
        smw_err : float, array
            Error on 'smw'.

    Calibration for HARPS from Lovis et al. (2011).
    Calibration for ESPRESSO based on 27 stars with data from HARPS and ESPRESSO.
    """
    caii = np.asarray(caii)
    caii_err = np.asarray(caii_err)

    if instr == 'ESPRESSO':
        a = 1.222
        b = 0.001
    if instr == 'HARPS':
        """Lovis et al. (2010), based on 7 stars
        with 0.137 < SMW < 0.393 from Baliunas et al. (1995).
        Calibration used by HARPS pipeline."""
        a = 1.111
        b = 0.0153

    smw = a * caii + b
    smw_err = a * caii_err

    return smw, smw_err


def calc_rhk(smw, smw_err, bv, method='middelkoop', lum_class='V'):
    """Calculates logR'HK via Noyes et al. (1984) with bolometric corrections using Middelkoop (1982), Rutten (1984), or Suárez Mascareño (2015, 2016) relations.

    Parameters:
    -----------
    smw : float, list, array
        S-index calibrated to the Mt. Wilson scale.
    smw_err : float, list, array
        Error on 'smw'.
    bv : float
        B-V colour.
    method : string
        Method to be used to calculate bolometric correction, Ccf: 'middelkoop' (default), 'rutten', or 'mascareno'.
    lum_class : string
        If using 'rutten' method, use 'V', 'VI' (default) if star is Main Sequence or 'III', 'IV' if star is giant or subgiant. IMPORTANT: the 'middelkoop' and 'mascareno' methods are only meant for Main Sequence (lum_class='V') stars.

    Returns:
    --------
    log_rhk : float, array
        Logarithm (base 10) of the R'HK chromosperic emission index.
    log_rhk_err : float, array
        Error on log_rhk.
    rhk : float, array
        R'HK chromospheric emission index.
    rhk_err : float, array
        Error on R'HK.
    method : string
        Method used to calculate Ccf.

    The calibration used by the HARPS pipeline is the 'middelkoop', the most widely used. Only for main sequence stars.
    The 'rutten' calibration is more useful if using dwarfs hotter than B-V = 0.44, cooler than B-V = 1.2, and/or giants (and subgiants).
    The 'mascareno' calibration includes M-dwarf stars. Only for MS.

    Range of the 'middelkoop' calibration (V): 0.44 < B-V < 1.20
    Range of the 'rutten' calibration (V): 0.30 < B-V < 1.60
    Range of the 'rutten' calibration (IV, III): 0.30 < B-V < 1.70
    Range of the 'mascareno' calibration (V): 0.40 < B-V < 1.90

    NOTE: If the B-V value is out of range the result will be 'NaN'.
    """
    smw = np.asarray(smw)
    smw_err = np.asarray(smw_err)
    bv = float(bv)

    if not isinstance(method, str) or method not in ('middelkoop', 'rutten', 'mascareno'):
        print("*** ERROR: 'method' should be 'middelkoop', 'rutten', or 'mascareno.")

    if not isinstance(lum_class, str) or lum_class not in ('VI', 'V', 'IV', 'III'):
        print("*** ERROR: 'lum_class' should be 'VI', 'V', 'VI' or 'III'.")

    if method == 'middelkoop':
        if lum_class in ('VI', 'V'):
            if (bv > 0.44) & (bv < 1.20):
                logCcf = 1.13*bv**3 - 3.91*bv**2 + 2.84*bv - 0.47
                if bv < 0.63:
                    logCcf = logCcf + 0.135*(0.63-bv) - 0.814*(0.63-bv)**2 + 6.03*(0.63-bv)**3
            else:
                logCcf = np.nan
        else:
            logCcf = np.nan

    elif method == 'rutten':
        if lum_class in ('VI', 'V'):
            if (bv >= 0.3) & (bv <= 1.6):
                logCcf = 0.25*bv**3 - 1.33*bv**2 + 0.43*bv + 0.24
            else:
                logCcf = np.nan
        elif lum_class in ('IV', 'III'):
            if (bv >= 0.3) & (bv <= 1.7):
                logCcf = -0.066*bv**3 - 0.25*bv**2 - 0.49*bv + 0.45
            else:
                logCcf = np.nan
        else:
            logCcf = np.nan

    elif method == 'mascareno':
        if lum_class in ('VI', 'V'):
            if (bv >= 0.4) & (bv <= 1.9):
                logCcf = 0.668 - 1.270*bv + 0.645*bv**2 - 0.443*bv**3
            else:
                logCcf = np.nan
        else:
            logCcf = np.nan


    if logCcf:
        Ccf = 10**logCcf
        # Noyes et al. (1984):
        r = 1.34e-4*Ccf*smw
        r_err = 1.34e-4*Ccf*smw_err

        if method in ("middelkoop", "rutten"):
            # Hartmann et al. (1984):
            log_rphot = -4.898 + 1.918*bv**2 - 2.893*bv**3
            rphot = 10**log_rphot
        elif method == "mascareno":
            rphot = 1.48e-4 * np.exp(-4.3658 * bv)


        if np.any(r-rphot) > 0.0:
            log_rhk = np.log10(r-rphot)
            log_rhk_err = r_err/(r-rphot)/np.log(10)
            rhk = r - rphot
            #rhk_err = 10**log_rhk_err
            rhk_err = r_err
        else:
            log_rhk = np.nan
            log_rhk_err = np.nan
            rhk = np.nan
            rhk_err = np.nan
    else:
        log_rhk = np.nan
        log_rhk_err = np.nan
        rhk = np.nan
        rhk_err = np.nan

    return log_rhk, log_rhk_err, rhk, rhk_err, method


def calc_prot_age(log_rhk, bv):
    """Calculates gyrochronological rotation period and age based based on the empirical relations of Noyes et al. (1984) and Mamajek & Hillenbrand (2008).

    Parameters:
    -----------
    log_rhk : float, list, array
        Logarithm (base 10) of the R'HK index.
    bv : float
        B-V colour.

    Returns:
    --------
    prot_n84 : float, array
        Rotational period via Noyes et al. (1984).
    prot_m84_err : float, array
        Error on 'prot_n84'.
    prot_m08 : float, array
        Rotational period via Mamajek & Hillenbrand (2008).
    prot_m08_err : float, array
        Error on 'prot_m08'
    age_m08 : float, array
        Age via Mamajek & Hillenbrand (2008).
    age_m08_err : float, array
        Error on 'age_m08'.

    Range of logR'HK-Prot relation: -5.5 < logR'HK < -4.3
    Range of Mamajek & Hillenbrand (2008) relation for ages: B-V >= 0.5
    """
    log_rhk = np.asarray(log_rhk)
    bv = float(bv)

    # Calculate Prot:
    if np.any(log_rhk < -4.3) & np.any(log_rhk > -5.5):
        if bv < 1:
            tau = 1.362 - 0.166*(1-bv) + 0.025*(1-bv)**2 - 5.323*(1-bv)**3
        else:
            tau = 1.362 - 0.14*(1-bv)

        prot_n84 = 0.324 - 0.400*(5 + log_rhk) - 0.283*(5 + log_rhk)**2 - 1.325*(5 + log_rhk)**3 + tau
        prot_n84 = 10**prot_n84
        prot_n84_err = np.log(10)*0.08*prot_n84

        prot_m08 = (0.808 - 2.966*(log_rhk + 4.52))*10**tau
        prot_m08_err = 4.4*bv*1.7 - 1.7
    else:
        prot_n84 = np.nan
        prot_n84_err = np.nan
        prot_m08 = np.nan
        prot_m08_err = np.nan

    # Calculate age:
    if np.any(prot_m08 > 0.0) & (bv >= 0.50):
        age_m08 = 1e-3*(prot_m08/0.407/(bv - 0.495)**0.325)**(1./0.566)
        age_m08_err = 0.05*np.log(10)*age_m08
    else:
        age_m08 = np.nan
        age_m08_err = np.nan

    return prot_n84, prot_n84_err, prot_m08, prot_m08_err, age_m08, age_m08_err


def get_bv(star_id, alerts=True):
    """Obtain B-V colour from Simbad.

    Parameters:
    -----------
    star_id : string
        Target identification readable by Simbad.
    alerts : bool
        If 'True' (default), errors are printed on screen.

    Returns:
    --------
    bv : float
        B-V colour from Simbad.
    bv_err : float
        Error on 'bv'.
    bv_ref : string
        Reference of flux V magnitude (generally the same as B mag).
    err_msg : string
        Error message. If 'OK', there was no error.
    """
    customSimbad = Simbad()
    customSimbad.add_votable_fields('flux(V)')
    customSimbad.add_votable_fields('flux_error(V)')
    customSimbad.add_votable_fields('flux_bibcode(V)')
    customSimbad.add_votable_fields('flux(B)')
    customSimbad.add_votable_fields('flux_error(B)')
    customSimbad.get_votable_fields()

    err_msg = 'OK'

    try:
        query = customSimbad.query_object(star_id)
    except:
        err_msg = f"*** ERROR: Could not identify {star_id}."
        if alerts:
            print(err_msg)
        return np.nan, np.nan, np.nan, err_msg

    if query is None:
        err_msg = f"*** ERROR: Could not identify {star_id}."
        return np.nan, np.nan, np.nan, err_msg

    flux_v = query['FLUX_V'][0]
    flux_v_err = query['FLUX_ERROR_V'][0]
    flux_v_ref = query['FLUX_BIBCODE_V'][0].decode("UTF-8")
    flux_b = query['FLUX_B'][0]
    flux_b_err = query['FLUX_ERROR_B'][0]
    const = np.ma.core.MaskedConstant

    if isinstance(flux_b, const) or isinstance(flux_v, const):
        err_msg = f"*** ERROR: {star_id}: No values of B and/or V in Simbad to calculate B-V."
        if alerts:
            print(err_msg)
        return np.nan, np.nan, np.nan, err_msg
    else:
        bv = flux_b - flux_v

    if isinstance(flux_v_err, const): flux_v_err = np.nan
    if isinstance(flux_b_err, const): flux_b_err = np.nan
    if isinstance(flux_v_ref, const): flux_v_ref = np.nan

    bv_err = np.sqrt(flux_b_err**2 + flux_v_err**2)
    bv_ref = flux_v_ref

    return bv, bv_err, bv_ref, err_msg
