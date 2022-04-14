# pyrhk
Python functions to calculate logR'HK for HARPS and ESPRESSO and estimate rotation periods and ages.

What pyrhk can do:
- Calibrate S-index values to the Mt. Wilson scale (SMW) for HARPS (using Gomes da Silva et al. 2021 or Lovis et al. 2011 calibrations) or ESPRESSO (using the calibration shown below) spectrographs.
- Calculate logR'HK via Noyes et al. (1984) using two bolometric corrections:
    - Middelkoop (1982): for 0.44 < B-V < 1.20 (MS stars)
    - Rutten (1984):  for 0.3 < B-V < 1.6 (MS stars) and 0.3 < B-V < 1.7 (evolved stars)
- Calculate chromospheric rotation period and gyrochronology age using activity relations from Noyes et al. (1984) and Mamajek & Hillenbrand (2008).
- Obtain B-V colour from Simbad (requires `astroquery` module installed).


### Calibration of SMW for ESPRESSO:

Using 27 stars with data from HARPS and ESPRESSO.


![SMW calibration for ESPRESSO](smw_espresso_cal_all.png "SMW ESPRESSO calibration")
