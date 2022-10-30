"""
J. Gomes da Silva, 2019
Instituto de Astrofísica e Ciências do Espaço (IA)
Centro de Astrofísica da Universidade do Porto (CAUP)
Joao.Silva@astro.up.pt
"""
import numpy as np
try:
    from astroquery.simbad import Simbad
except:
    print("*** WARNING: astroquery is not installed but required for the get_bv function")




def calc_smw(caii, caii_err, instr='HARPS_GDS21'):
    """Calibrates the S-index to the Mt. Wilson scale (Vaughan et al. 1978).

    Parameters:
    -----------
        caii : float, array
            CaII index based on the CaII H&K lines.
        caii_err : float, array
            Photon error of the CaII index.
        instr : string
            Instrument options:
            'HARPS_GDS21':
                (default) Gomes da Silva et al. (2021), S-index calculated with ACTIN (https://github.com/gomesdasilva/ACTIN), based on 43 stars with 0.105 < SMW < 0.496 from Duncan et al. (1991) and Baliunas et al. (1995).
            'HARPS_L11':
                Lovis et al. (2011), based on 7 stars with 0.137 < SMW < 0.393 from Baliunas et al. (1995). Calibration used by HARPS DRS pipeline.
            'ESPRESSO':
                Preliminary calibration based on 27 stars in common between HARPS and ESPRESSO.
    Returns:
    --------
        smw : float, array
            S-index calibrated to the Mt. Wilson scale.
        smw_err : float, array
            Error on 'smw'.
    """
    caii = np.asarray(caii)
    caii_err = np.asarray(caii_err)

    if instr == 'ESPRESSO':
        a = 1.222
        b = 0.001
    if instr == 'HARPS_L11':
        a = 1.111
        b = 0.0153
    if instr == 'HARPS_GDS21':
        a = 1.195
        b = 0.008

    smw = a * caii + b
    smw_err = a * caii_err

    return smw, smw_err


def calc_rhk(smw, smw_err, bv, method='middelkoop', evstage='MS'):
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
        Method used to calculate bolometric correction, Ccf: 'middelkoop' (default), 'rutten', or 'mascareno'.
    evstage : string,
        Evolutionary stage. If using 'rutten' method, use 'MS' if star is in the Main Sequence or 'evol' if star is evolved (giant or subgiant). IMPORTANT: the 'middelkoop' and 'mascareno' methods are only meant for Main Sequence (evstage='MS') stars (default).

    Returns:
    --------
    log_rhk : float, array
        Logarithm (base 10) of the R'HK chromosperic emission ratio.
    log_rhk_err : float, array
        Error on log_rhk.
    rhk : float, array
        R'HK chromospheric emission ratio.
    rhk_err : float, array
        Error on R'HK.

    The calibration used by the HARPS pipeline is the 'middelkoop', the most widely used. Only for main sequence stars.
    The 'rutten' calibration is more useful if using evolved stars (giants and subgiants).
    The 'mascareno' calibration includes cooler M-dwarfs. Only for MS.

    Range of the 'middelkoop' calibration (MS): 0.44 < B-V < 1.20
    Range of the 'rutten' calibration (MS): 0.30 < B-V < 1.60
    Range of the 'rutten' calibration (evol): 0.30 < B-V < 1.70
    Range of the 'mascareno' calibration (MS): 0.40 < B-V < 1.90

    NOTE: If the B-V value is out of range the result will be 'np.nan'.
    """
    smw = np.asarray(smw)
    smw_err = np.asarray(smw_err)

    if not isinstance(method, str) or method not in ('middelkoop', 'rutten', 'mascareno'):
        print("*** ERROR: 'method' should be 'middelkoop', 'rutten', or 'mascareno.")
        return

    if not isinstance(evstage, str) or evstage not in ('MS', 'evol'):
        print("*** ERROR: 'evstage' should be 'MS' or 'evol'.")
        return

    if method == 'middelkoop':
        if evstage in ('MS'):
            if (bv > 0.44) & (bv < 1.20):
                logCcf = 1.13*bv**3 - 3.91*bv**2 + 2.84*bv - 0.47
                if bv < 0.63:
                    logCcf = logCcf + 0.135*(0.63-bv) - 0.814*(0.63-bv)**2 + 6.03*(0.63-bv)**3
            else:
                logCcf = np.nan
        else:
            logCcf = np.nan

    elif method == 'rutten':
        if evstage in ('MS'):
            if (bv >= 0.3) & (bv <= 1.6):
                logCcf = 0.25*bv**3 - 1.33*bv**2 + 0.43*bv + 0.24
            else:
                logCcf = np.nan
        elif evstage in ('evol'):
            if (bv >= 0.3) & (bv <= 1.7):
                logCcf = -0.066*bv**3 - 0.25*bv**2 - 0.49*bv + 0.45
            else:
                logCcf = np.nan
        else:
            logCcf = np.nan

    elif method == 'mascareno':
        if evstage in ('MS'):
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
            #log_rphot_err = (2*1.918*bv_err) - 3 * 2.893 * bv_err**2
            rphot = 10**log_rphot
        elif method == "mascareno":
            rphot = 1.48e-4 * np.exp(-4.3658 * bv)

        if np.any(r-rphot > 0.0):
            log_rhk = np.log10(r-rphot)
            #log_rhk_err = np.sqrt((r_err/(r-rphot)/np.log(10))**2 + log_rphot_err**2)
            log_rhk_err = r_err/(r-rphot)/np.log(10)
            rhk = r - rphot
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

    return log_rhk, log_rhk_err, rhk, rhk_err


def calc_prot_age(log_rhk, bv):
    """Calculates rotation period and age from activity level, based on the empirical relations of Noyes et al. (1984) and Mamajek & Hillenbrand (2008).

    Parameters:
    -----------
    log_rhk : float, list, array
        Logarithm (base 10) of the R'HK index.
    bv : float
        B-V colour.

    Returns:
    --------
    prot_n84 : float, array
        Chromospheric rotational period via Noyes et al. (1984).
    prot_m84_err : float, array
        Error on 'prot_n84'.
    prot_m08 : float, array
        Chromospheric rotational period via Mamajek & Hillenbrand (2008).
    prot_m08_err : float, array
        Error on 'prot_m08'
    age_m08 : float, array
        Gyrochronology age via Mamajek & Hillenbrand (2008).
    age_m08_err : float, array
        Error on 'age_m08'.

    Range of logR'HK-Prot relation: -5.5 < logR'HK < -4.3
    Range of Mamajek & Hillenbrand (2008) relation for ages: 0.5 < B-V < 0.9
    """
    log_rhk = np.asarray(log_rhk)
    bv = float(bv)

    # Calculate chromospheric Prot:
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

    # Calculate gyrochronology age:
    if np.any(prot_m08 > 0.0) & (bv > 0.50) & (bv < 0.9):
        age_m08 = 1e-3*(prot_m08/0.407/(bv - 0.495)**0.325)**(1./0.566)
        #age_m08_err = 0.05*np.log(10)*age_m08
        age_m08_err  = 0.2 * age_m08 * np.log(10) # using 0.2 dex typical error from paper
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
        return np.nan, np.nan, np.nan

    if query is None:
        err_msg = f"*** ERROR: Could not identify {star_id}."
        return np.nan, np.nan, np.nan, err_msg

    flux_v = query['FLUX_V'][0]
    flux_v_err = query['FLUX_ERROR_V'][0]
    flux_v_ref = query['FLUX_BIBCODE_V'][0]#.decode("UTF-8")
    flux_b = query['FLUX_B'][0]
    flux_b_err = query['FLUX_ERROR_B'][0]
    const = np.ma.core.MaskedConstant

    if isinstance(flux_b, const) or isinstance(flux_v, const):
        err_msg = f"*** ERROR: {star_id}: No values of B and/or V in Simbad to calculate B-V."
        if alerts:
            print(err_msg)
        return np.nan, np.nan, np.nan
    else:
        bv = flux_b - flux_v

    if isinstance(flux_v_err, const): flux_v_err = np.nan
    if isinstance(flux_b_err, const): flux_b_err = np.nan
    if isinstance(flux_v_ref, const): flux_v_ref = np.nan

    bv_err = np.sqrt(flux_b_err**2 + flux_v_err**2)
    bv_ref = flux_v_ref

    return bv, bv_err, bv_ref
