# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

import numpy as np
import scipy.interpolate as si
import scipy.constants   as sc

def interp(inten, grav, temp, wgrav, wtemp, log=False, kf=False):
  """
 NAME:
  INTERP

 PURPOSE:
  This function finds the brightness spectrum of a star by
  interpolating temperature and gravity in a grid of Kurucz
  models.  Wavelengths/frequencies are not interpolated.

 CATEGORY:
  Astronomy

 CALLING SEQUENCE:

  Result = INTERP(Inten, Grav, Temp, Wgrav, Wtemp)

 INPUTS:
  Inten:  [nwavl, nmod] array of model brightnesses
    (erg cm-2 s-1 sr-1 Hz-1).

  Grav:  [nmod] array of the base-10 log of stellar surface
    gravities for the models (cm s-2).

  Temp:  [nmod] array of temperatures for the models (K).

  Wgrav:  Wanted log(gravity) (cm s-2).

  Wtemp:  Wanted temperature (K).

 KEYWORD PARAMETERS:
  LOG:  If set, interpolate in the log.
  KF:  [nwavl, ngrav, ntemp] array of Kurucz models (returned)
    (returned, erg cm-2 s-1 sr-1 Hz-1).

 OUTPUTS:
  This function returns an [nwavl] array of brightnesses
  (erg cm-2 s-1 sr-1 Hz-1).

 PROCEDURE:
  Uses cubic convolution in 2D, done explicitly at each
  wavelength.  Note that /log makes very little difference (see
  example).

  Kurucz's modeled grav and temp do not cover a rectangular
  grid, since he does not calculate low gravities for large
  stars that never have them.  The code currently just takes the
  rectangular region filled by his cooler models and makes a
  grid out of them, ignoring the rest.  If your parameters are
  too close to the edge of this grid, the code warns and exits.
  It will not be hard to recode to accomodate more models.

 EXAMPLE:
  ;Example is still in IDL:
  inten = kurucz_inten_read('fp00ak2odfnew.pck', $
        wave, grav, temp, nainten, head)
  ; find the Sun's spectrum
  wgrav = 4.44d ; log(cm s-2)
  wtemp = 5770d ; Kelvins
  kinten = kurucz_inten_interp(inten, grav, temp, wgrav, wtemp)
  c     = 2.99792458d8 ; speed of light in m/s
  freq  = reverse(c / wave)
  fsun  = reverse(kinten)
  plot, freq, fsun, /xlog, $
    xtitle='Frequency (Hz)', ytitle='I (W m-2 sr-1 Hz-1)'

  ; find the Sun's luminosity
  ; extra factor of !dpi converts brightness to flux
  print, int_tabulated(freq, fsun, /d) * !dpi * 4d * !dpi * 6.96d8^2
  ;   3.8392897e+26

  ; Compare to blackbody...
  pf = planckfreq(freq, wtemp) / !dpi * 1d-7 * 1d4 ; flux -> inten., MKS
  plot, freq, fsun, /xlog, $
    xtitle='Frequency (Hz)', ytitle='I (W m-2 sr-1 Hz-1)'
  oplot, freq, pf
  print, int_tabulated(freq, pf, /d) * !dpi * 4d * !dpi * 6.96d8^2
  ;   3.8260749e+26

  ; check the ratio
  irat = fsun / pf
  plot, freq, irat, /xlog, $
    xtitle='Frequency (Hz)', ytitle='Ikurucz / Iplanck'

  ; how much difference does log interpolation make?
  kinten2 = kurucz_inten_interp(inten, grav, temp, wgrav, wtemp, /log)
  afd = abs( (kinten-kinten2) / ((kinten+kinten2) / 2) )
  plot, wave, afd[where(kinten gt max(kinten) * 1e-2)], /xlog, /ylog

  For the sun, the fractional difference is 2e-3 or better for
  brightnesses greater than 1% of the peak, except for one spike
  to 5e-3.  Differences longward of 0.1 um are below 5e-5 and
  fall rapidly with increasing wavelength.

 MODIFICATION HISTORY:
   Written by:  Joseph Harrington, Cornell.  2006-02-02
      jh@oobleck.astro.cornell.edu
  2006-02-03 jh  Swapped wgrav and wtemp in calling order to
      match grav and temp.  Updated example for
      change of kurucz_flux_read to MKS.
  2006-02-10 jh  Fixed example, added kf keyword, changed name
      to kurucz_inten_interp.
  2006-02-11 jh  Header tweak.
  2006-02-12 jh  Converted variable names flux -> inten.
      Changed messages to printed warnings, since
      they make Monte Carlo trials bomb.
  2008-08-11 kevin
      Converted to Python
  """

  dim    = inten.shape
  nwav   = dim[1]
  iinten = np.zeros(nwav)
  
  # FINDME: HARDWIRED!
  ilograv = 0
  ihigrav = 10 # the full range of gravities
  ilotemp = 0
  ihitemp = 16 # the temp where Kurucz starts skipping gravities
  ngrav = ihigrav - ilograv + 1
  ntemp = ihitemp - ilotemp + 1
  
  kg = grav[ilograv:ihigrav+1]
  it = np.arange(ntemp) * ngrav + ilotemp
  kt = temp[it]
  
  # FINDME: this only works if the low indices are both 0 and we're not
  # at temps where Kurucz starts skipping gravities...
  if kg[ihigrav] < wgrav:   
    print('kurucz_inten.interp: wgrav ' + str(wgrav) +
          ' higher than ' + str(kg[ihigrav]) +
          ', see FINDME in code to fix.')
  if kt[ihitemp] < wtemp:   
    print('kurucz_inten.interp: wtemp ' + str(wtemp) +
          ' higher than ' + str(kt[ihitemp]) +
          ', see FINDME in code to fix.')
  kf = np.reshape(inten[0:ngrav*ntemp, :], (ntemp, ngrav, nwav))

  # Calculate the indices of the parameters we want:
  itemp, igrav = np.mgrid[kt[0]:kt[-1]+1:250, kg[0]:kg[-1]+0.1:0.5]

  if log:
    kf = np.log(kf)

  for i in np.arange(nwav):
    tck = si.bisplrep(itemp, igrav, kf[:,:,i], kx=3, ky=3)
    iinten[i] = si.bisplev(wtemp, wgrav, tck)

  if log:   
    iinten = exp(iinten)

  return iinten


def read(filename, freq=False):
  """
  This function reads a file of stellar spectral intensity models
  from Bob Kurucz (Harvard).

  Parameters:
  -----------
  filename: String
     Name of model file.  These come from http://kurucz.harvard.edu/grids.html

  freq: Boolean
     If True, reverse first dimension of model grid and return frequencies 
     instead of wavelengths in Wave.

  Returns:
  --------
  inten: 2D ndarray
     Array of shape (nmod, nwavl) with the models brightnesses, with
     nwavl the number of wavelength samples, and nmod the number of models.
     These brightnesses in the file have been multiplied by 4, since they
     are Eddington fluxes (W m-2 sr-1 Hz-1).
  wave: 1D ndarray
     Array of size nwavl of the wavelengths (in meters) of inten. 
     If freq==True, wave contains the frequencies (in Hz).
  grav: 1D ndarray
     Array of size nmod with the log10 of the stellar ssurface
     gravities (g in cm s^-2) for the models.
  temp: 1D ndarray
     Array of size nmod with the temperature (in K) of the models.
  nainten: 2D ndarray
     Array of shape (nmod, nwavl) of the models brightnesses without
     line absorption.  Same units as inten.
  head: List of strings
     List of size nmod of the one-line header strings for the models.

  Example:
  --------
  >>> import kurucz_inten as ki
  >>> import time
  >>> import numpy as np
  >>> import matplotlib.pyplot as plt

  >>> # Read the model file:
  >>> kfile = '/home/esp01/ancil/kurucz/fp00ak2odfnew.pck'
  >>> inten, wave, grav, temp, nainten, head = ki.read(kfile)

  >>> # Plot the intensities vs. frequency in Hz:
  >>> nmod  = len(head)
  >>> wait = 0.05
  >>> plt.figure(0)
  >>> for i in np.arange(nmod):
  >>>   plt.clf()
  >>>   a = plt.loglog(wave, inten[i],   "b")
  >>>   a = plt.loglog(wave, nainten[i], "r")
  >>>   a = plt.xlim(5e-9, 5e-4)
  >>>   a = plt.ylim(1e-15, 5e-5)
  >>>   a = plt.title( "model %d:  T=%d  log10(g)=%.1f"%(i, temp[i], grav[i]))
  >>>   plt.draw()
  >>>   plt.show()
  >>>   time.sleep(wait)

  >>> inten, wave, grav, temp, nainten, head = ki.read(kfile, freq=True)
  >>> # Estimate the luminosity of the sun (~4e26 W):
  >>> Tsun = 5770.0  # Sun's surface temperature
  >>> gsun = 4.44    # Sun's surface gravity
  >>> # Find the appropriate model according to Tsun and gsun:
  >>> isun = np.where((temp > Tsun) & (grav > gsun))[0][0]
  >>> fsun = inten[isun]
  >>> # The first pi converts from brightness to flux.
  >>> Lum = np.trapz(fsun, wave) * np.pi * 4 * np.pi * 6.96e8**2
  >>> print("Luminosity  L = %.4e  W/(m^2 sr Hz)"%(Lum))
  >>> Luminosity  L = 4.4718e+26  W/(m^2 sr Hz)
  >>> # The model is for a star with t=6000 and log g=4.5, so expect
  >>> # more than 4e26 W.
 
  Modification History:
  ---------------------
  2006-02-01  jh        Written by Joseph Harrington, Cornell.
                        Based on code provided by Drake Deming.
                                               jh@oobleck.astro.cornell.edu
  2006-02-02  jh        Added solar luminosity example, fixed header,
                        swapped order of grav and temp to match order
                        in grid of models.
  2006 Feb 3  jh        Convert to MKS, update example.
  2006-02-10  jh        Rename to kurucz_inten_read, update header.
  2008-08-11  kevin     Converted to Python by Kevin Stevenson, UCF
                                               kbstvenson@gmail.com
  2014-03-26  patricio  Adapted for Output converter, fixed format and
                        updated docstring by Patricio Cubillos, UCF
                                               pcubillos@fulbrightmail.org
  """
  # Read file into memory:
  f = open(filename, 'r')
  text = f.read()
  f.close()
  text = text.replace('\r','\n')
  filetxt = text.split('\n')

  # Get, parse, and count header lines:
  # Record effective temperatures and gravities:
  head      = []
  temp      = np.zeros(len(filetxt))
  grav      = np.zeros(len(filetxt))-1
  header    = np.zeros(len(filetxt), int)
  startwave = 0
  for i in np.arange(len(filetxt)):
    if filetxt[i].startswith("TEFF"):
      head.append(filetxt[i])
      temp[i]   = float(filetxt[i][ 5:12])
      grav[i]   = float(filetxt[i][22:29])
      header[i] = i
    elif filetxt[i].endswith("END"):
      startwave = i + 1

  temp   = temp  [np.where(temp   !=  0)]
  grav   = grav  [np.where(grav   != -1)]
  header = header[np.where(header !=  0)]
  nmod   = header.size
  nline  = (header[2] - header[1] - 1) // 2  # Omit the header line itself

  # Read and count wavelengths:
  wave = np.zeros(header[0]*len(filetxt[startwave])//10)
  k = 0
  string = ''.join(filetxt[startwave:header[0]])
  for j in np.arange(0, len(string), 10):
    wave[k] = float(string[j:j+10])
    k += 1

  wave = wave[np.where(wave != 0)] * 1e-9  # Convert nm to meters
  nwavl = wave.size

  # Allocate memory for models:
  inten   = np.zeros((nmod, nwavl))
  nainten = np.zeros((nmod, nwavl))

  #LOOP OVER MODELS
  for i in range(0, nmod):
    k = 0
    string1 = ''.join(filetxt[header[i]+1      :header[i]+1+nline  ])
    string2 = ''.join(filetxt[header[i]+1+nline:header[i]+1+nline*2])
    for j in range(0,len(string1),10):
      inten[i,k]   = float(string1[j:j+10])
      nainten[i,k] = float(string2[j:j+10])
      k += 1

  # Convert Eddington fluxes to brightnesses, CGS (erg cm-2) -> MKS (J m-2)
  inten   *= 4.0 * 1e-3
  nainten *= 4.0 * 1e-3

  # Convert to frequency if requested:
  if freq:   
    wave    = np.flipud(sc.c / wave)
    inten   = np.fliplr(inten)
    nainten = np.fliplr(nainten)
  
  return inten, wave, grav, temp, nainten, head
