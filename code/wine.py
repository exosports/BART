# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

import numpy as np
import kurucz_inten      as ki
import scipy.constants   as sc
import scipy.interpolate as si

"""
WINE: Waveband INtegrated Emission module

This set of routines calculate the integrated emission spectrum of
a signal over specified filter wavebands.
"""

def readfilter(filt):
  """
  Load a filter bandpass from file.

  Notes:
  ------
  - The file can contains empty lines and comments (with '#' character)
    before the data.  No comments or empty lines after the data.
  - The data must come in two columns.  The first column must contain
    the wavelength in microns, the second the filter response, other
    columns will be ignored.

  Parameters:
  -----------
  filt: String
      Filter file name.
  Return:
  -------
  wavenumber: 1D ndarray
     The filter pass band wavenumber in cm^-1.
  transmission: 1D ndarray 
     The filter spectral response. No specific units.
  
  Modification History:
  ---------------------
  2013-01-23  patricio  Initial implementation.   pcubillos@fulbrightmail.org
  2014-03-26  patricio  Changed input to the file name. 
  """
  # Open and read the filter file:
  data = open(filt, "r")
  lines = data.readlines()
  data.close()

  # Remove header comments and empty lines:
  while lines[0].startswith("#") or not lines[0].strip():
    comment = lines.pop(0)

  # Allocate arrays for the wavelength and response:
  nlines = len(lines)
  wavel  = np.zeros(nlines, np.double) # filter's wavelengths  (in microns)
  transm = np.zeros(nlines, np.double) # filter's pass bands

  # Read the data and store in reverse order:
  for i in np.arange(nlines):
    wavel[nlines-1-i], transm[nlines-1-i] = lines[i].strip().split()[0:2]

  m2cm  = 1e-4 # Microns to cm conversion factor
  # Get wavenumber in cm-1:
  waven = 1.0 / (wavel*m2cm)
  # Return statement:
  return waven, transm


def readkurucz(kfile, temperature, logg):
  """
  Load a the Kurucz stellar spectrum with parameters closest to requested
  temperature and log(g).

  Parameters:
  -----------
  kfile: String
     Path to the kurucz file.
  temperature: Scalar
     Surface temperature in K.
  logg: Scalar
     log10 of surface gravity (g in cgs units).

  Returns:
  --------
  starfl: 1D ndarray
     Kurucz stellar flux in ergs s-1 cm-2 cm.
  starwn: 1D ndarray
     Array with wavenumber values in cm^-1.
  tmodel: Scalar
     Surface temperature of the model in K.
  gmodel: Scalar
     log10 of the surface gravity for the model (g in cgs units).

  Modification History:
  ---------------------
  2013-01-23  patricio  Initial implementation.   pcubillos@fulbrightmail.org
  """

  inten, freq, grav, temp, nainten, head = ki.read(kfile, freq=True)

  # Wavenumber in cm^-1
  starwn = freq / sc.c * 1e-2  
  
  # Find the model index with the nearest temp and log(g):
  # Nearest sampled temperature:
  tmodel = temp[np.argmin(np.abs(temp-temperature))]
  # Nearest sampled log(g):
  gmodel = grav[np.argmin(np.abs(grav-logg))]
  imodel = np.where((temp == tmodel) & (grav >= gmodel))[0][0]

  # Get the stellar flux:
  starfl = inten[imodel]  # W m^-2 sr^-1 Hz^-1

  # Convert F_freq to F_wavenumber (Hz-1 --> m):
  #   multiply by c.
  # Convert units MKS to cgs:
  #   W m-2 = 1e3 ergs s-1 cm-2
  # Convert intensity (astrophysical flux) to flux:
  #   sr-1 = pi

  # Flux per wavenumber:  ergs s-1 cm-2 cm
  starfl = starfl * 1e3 * np.pi * (1e2 * sc.c)

  return starfl, starwn, tmodel, gmodel


def resample(specwn, filterwn, filtertr, starwn, starfl):
  """
  Resample the filtertr curve from the filterwn sampling into specwn

  Parameters:
  -----------
  specwn: 1D ndarray
     A wavenumber sampling array (in cm-1).
  filterwn: 1D ndarray
     Filter wavenumber sampling array (in cm-1).
  filtertr: 1D ndarray
     Filter transmission curve sampled as filterwn.
  starwn: 1D ndarray
     Stellar model wavenumber sampling array (in cm-1).
  starfl: 1D ndarray
     Stellar flux.

  Returns:
  --------
  ifilter: 1D ndarray
     The interpolated filter transmission curve.
  istarfl: 1D ndarray
     The interpolated stellar flux.
  indices: 1D ndarray
     The indices of specwn where the filter was interpolated into.

  Modification History:
  ---------------------
  2013-01-23  patricio  Initial implementation. 
  2014-03-26  patricio  Adapted for output converter as resample. 
  """
  # Indices in the spectrum wavenumber array included in the band
  # wavenumber range:
  wnindices = np.where((specwn < filterwn[-1]) & (filterwn[0] < specwn))

  # Make function to spline-interpolate the filter and stellar flux:
  finterp = si.interp1d(filterwn, filtertr)
  sinterp = si.interp1d(starwn,   starfl)

  # Evaluate the stellar spectrum on specwn:
  istarfl = sinterp(specwn[wnindices])
  # Evaluate over the spectrum wavenumber array:
  ifilter = finterp(specwn[wnindices])
  # Normalize to integrate to 1.0:
  nifilter = ifilter/np.trapz(ifilter, specwn[wnindices])

  # Return the normalized interpolated filter and the indices:
  return nifilter, istarfl, wnindices


def bandintegrate(spectrum, specwn, nifilter, wnindices):
  """
  Integrate a spectrum over the band transmission.

  Parameters:
  -----------
  spectrum: 1D ndarray
     Spectral signal to be integrated
  specwn: 1D ndarray
     Wavenumber of spectrum in cm^-1
  nifilter: 1D ndarray
     The normalized interpolated filter transmission curve.
  wnindices: 1D ndarray
     Indices of specwn where bandtr is evaluated.

  Modification History:
  ---------------------
  2014-03-26  patricio  Initial implementation.
  """
  # Flux ratio:
  # fratio = Fplanet / Fstar * rprs**2.0

  return np.trapz(spectrum*nifilter, specwn[wnindices])
