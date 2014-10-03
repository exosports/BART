#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import PT as pt

def initialPT(date_dir, tepfile, press_file, a1, a2, p1, p3, T3_fac):
  """
  This function generates a non-inverted temperature profile using the
  parametrized model described in Madhusudhan & Seager (2009). It plots
  the profile to screen and save the figures in the output directory.  
  The generated PT profile is a semi-adiabatic profile with a bottom-layer
  temperature corresponding (within certain range) to the planetary
  effective temperature.

  Parameters
  ----------
  date_dir: String
     Directory where to save the output plots.
  tepfile: String
     Name of ASCII tep file with planetary system data.
  press_file: String
     Name of ASCII file with pressure array data.
  a1: Float
     Model exponential factor in Layer 1.
  a2: Float
     Model exponential factor in Layer 2.
  p1: Float
     Pressure boundary between Layers 1 and 2 (in bars).
  p3: Float
     Pressure boundary between Layers 2 and 3 (in bars).
  T3_fac: Float
     Multiplicative factor to set T3 (T3 = Teff * T3_fac).
     Empirically determined to be between (1, 1.5) to account
     for the possible spectral features.

  Notes
  -----
  See model details in Madhusudhan & Seager (2009):
  http://adsabs.harvard.edu/abs/2009ApJ...707...24M

  Returns
  -------
  T_smooth: 1D float ndarray
      Array of temperatures.

  Developers
  ----------
  Jasmina Blecic     jasmina@physics.ucf.edu
  Patricio Cubillos  pcubillos@fulbrightmail.org

  Revisions
  ---------
  2014-04-08  Jasmina   Written by
  2014-07-23  Jasmina   Added date_dir and PT profile arguments.
  2014-08-15  Patricio  Replaced call to PT_Initial by PT_NoInversion.
  2014-09-24  Jasmina   Updated documentation.
  """

  # Calculate the planetary effective temperature from the TEP file
  Teff = pt.planet_Teff(tepfile)

  # Calculate T3 temperature based on Teff
  T3 = float(T3_fac) * Teff

  # Read pressures from file
  p = pt.read_press_file(press_file)

  # Generate initial PT profile
  PT, T_smooth = pt.PT_NoInversion(p, a1, a2, p1, p3, T3)

  # Take temperatures from PT generator
  T, T0, T1, T3 = PT[5], PT[7], PT[8], PT[9]

  # Plot raw PT profile
  plt.figure(1)
  plt.clf()
  plt.semilogy(PT[0], PT[1], '.', color = 'r'     )
  plt.semilogy(PT[2], PT[3], '.', color = 'b'     )
  plt.semilogy(PT[4], PT[5], '.', color = 'orange')
  plt.title('Initial PT', fontsize=14)
  plt.xlabel('T [K]', fontsize=14)
  plt.ylabel('logP [bar]', fontsize=14)
  plt.xlim(0.9*T0, 1.1*T3)
  plt.ylim(max(p), min(p))

  # Save plot to current directory
  plt.savefig(date_dir + '/InitialPT.png') 

  # Plot Smoothed PT profile
  plt.figure(2)
  plt.clf()
  plt.semilogy(T_smooth, p, '-', color = 'b', linewidth=1)
  plt.title('Initial PT Smoothed', fontsize=14)
  plt.xlabel('T [K]'     , fontsize=14)
  plt.ylabel('logP [bar]', fontsize=14)
  plt.ylim(max(p), min(p))
  plt.xlim(0.9*T0, 1.1*T3)

  # Save plot to output directory
  plt.savefig(date_dir + '/InitialPTSmoothed.png') 

  return T_smooth
