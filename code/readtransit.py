# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

# This wonderful piece of code reads the output modulation spectrum
# of a transit run and plots it.

import numpy as np
import matplotlib.pyplot as plt

def readplot(tfile, wn=True, fid=0):
  """
  Read transit's output and plot it.
  """
  wave, spectrum = readspectrum(tfile, wn)
  plt.figure(fid)
  plt.clf()
  plt.plot(wave, spectrum)
  plt.title(tfile)

  return wave, spectrum


def readspectrum(tfile, wn=True):
  """
  Read transit's output spectrum.

  Parameters:
  -----------
  tfile: String
     Path to output Transit spectrum file to read.
  wn: Boolean
     If True convert wavelength to wavenumber.

  Modification History:
  ---------------------
  2014-12-16  patricio  Initial version.
  """
  f = open(tfile, "r")

  # Count number of lines in file:
  f.seek(0)
  l = f.readline()  # Ignore first line of comments:
  ndata = 0
  for line in f:
    ndata += 1

  # Initialize arrays:
  wave     = np.zeros(ndata, np.double)
  spectrum = np.zeros(ndata, np.double)
  
  # Return to begining of file:
  f.seek(0)
  f.readline()
  i = 0
  for i in np.arange(ndata):
    l = f.readline().strip().split()
    wave    [i] = np.double(l[ 0])
    spectrum[i] = np.double(l[-1])
  
  # Convert wavelength (micron) to wavenumber (cm-1):
  if wn:
    wave = 1e4/wave
  
  return wave, spectrum
