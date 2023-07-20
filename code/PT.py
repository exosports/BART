# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.special   as sp
from scipy.ndimage import gaussian_filter1d
import reader as rd


"""
  This code serves as an input generator for BART. It generates
  parametrized PT profile done in a similar fashion as in Madhusudhan and
  Seager 2009, http://adsabs.harvard.edu/abs/2009ApJ...707...24M, but has a
  capability to be extended for other PT profiles.
  

  Functions
  ---------
  read_press_file:
     Read a pressure file and extract a list of pressures.
  planet_Teff:
     Calculate planet effective temperate to constrain the T3 parameter.
  PT_Inversion:
     Generate an inverted PT profile using the method of 
     Madhusudhan & Seager 2009.
  PT_NoInversion:
     Generate a non-inverted PT profile using the method of 
     Madhusudhan & Seager 2009.
  PT_line:
     Generate a PT profile using the method of Line et al. 2013
  PT_iso:
     Generate an isothermal PT profile.
  PT_generator:
     Wrapper that calls either inverted or non-inverted generator.
  plot_PT:
     Plot the PT profile.

  Developers:
  -----------
  Jasmina Blecic     UCF  jasmina@physics.ucf.edu
  Patricio Cubillos  UCF  pcubillos@fulbrightmail.org

  Revisions
  ---------
  2013-11-17  Jasmina  Written by.
  2014-04-05  Jasmina  Added new function planet_Teff, changed free
                       parameter T0 to T3 and equations and functions
                       accordingly
  2014-04-17  Jasmina  Added initialPT profile and free parameter generator
  2014-06-25  Jasmina  Instead from atm file the functions now take the
                       pressure from the pressure file. 
  2014-07-24  Jasmina  Adapted and reordered functions inside the module for
                       BART use.
  2014-08-15  Patricio Cleaned up.
  2014-09-24  Jasmina  Updated documentation.
  2018-05-24  Ryan     Added PT_iso
  2019-02-13  mhimes   Refactored PT_generator
"""


# extracts a pressure array from an pressure file provided
def read_press_file(press_file):
     '''
     Reads a pressure file. The function takes the column with pressures
     and converts it to floats.

     Parameters
     ----------
     press_file: String
        Name ASCII file that contains a pressure array data.

     Returns
     -------
     pressure: 1D float ndarray
        Array of pressures from press_file.

     Revisions
     ---------
     2014-06-19  Jasmina   Written by.
     2014-08-15  Patricio  Cleaned up.
     '''
  
     # Open and read the pressure array file
     f = open(press_file, 'r')
     lines = f.readlines()
     f.close()

     # Store the values in pressure array
     pressure = np.zeros(len(lines)-1, np.double)
     for i in np.arange(1, len(lines)):
          pressure[i-1] = lines[i].split()[1]

     return pressure


# reads the tep file and calculates planet's effective temperature
def planet_Teff(tepfile):
     '''
     Calculates planetary effective temperature. Calls tep reader 
     and gets data needed for effective temperature calculation.
     The effective temperature is calculated assuming zero albedo, and 
     zero redistribution to the night side, i.e uniform dayside 
     redistribution.

     Parameters
     ----------
     tepfile: string
           Name of the tep ASCII file.
 
     Returns
     -------
     Teff: float
           Effective temperature of the planet. 

     Revisions
     ---------
     2014-04-05   Jasmina    Written by.
     2014-08-15   Patricio   Added source for the radius of the sun.
     '''

     # Opens tepfile to read and get data
     tep = rd.File(tepfile)

     # Get stellar temperature in K
     stellarT = tep.getvalue('Ts')
     Tstar    = np.float(stellarT[0])

     # Get stellar radius in units of Rsun
     stellarR = tep.getvalue('Rs')
     Rstar    = np.float(stellarR[0])

     # Get semimajor axis in AU
     semimajor = tep.getvalue('a')
     a         = np.float(semimajor[0])

     # Sun radius in meters:
     # Source: http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
     Rsun = 696000000.0   # m

     # Radius fo the star and semimajor axis in m
     Rstar = Rstar * Rsun # m
     a     = a * sc.au    # m

     # Effective temperature of the planet 
     # Teff^4 = Teff*^4 * f * (Rstar/a)^2 * (1-A)
     # Zero albedo, no energy redistribution to the night side A=0, f=1/2 
     Teff = Tstar * (Rstar/a)**0.5 * (1./2.)**0.25

     return Teff


# generates PT profile for inverted atmosphere
def PT_Inversion(p, a1, a2, p1, p2, p3, T3, verb=False):
     '''
     Calculates PT profile for inversion case based on Equation (2) from
     Madhusudhan & Seager 2009.
     It takes a pressure array (e.g., extracted from a pressure file), and 6
     free parameters for inversion case and generates inverted PT profile. 
     The profile is then smoothed using 1D Gaussian filter. The pressure 
     array needs to be equally spaced in log space.

     Parameters
     ----------
     p:  1D array of floats
         Pressure array needs to be equally spaced in log space from bottom to top 
         of the atmosphere.
     a1: Float
         Model exponential factor in Layer 1, empirically determined to be within
         range (0.2, 0.6).
     a2: Float
         Model exponential factor in Layer 2, empirically determined to be within
         range (0.04, 0.5) 
     p1: Float
     p2: Float
         Pressure boundary between Layer  1 and 2 (in bars).
     p3: Float
         Pressure boundary between Layers 2 and 3 (in bars).
     T3: float
         Temperature in the Layer 3.
      
     Returns
     -------
     PT_Inver:  tupple of arrays that includes:
              - temperature and pressure arrays of every layer of the atmosphere 
                (PT profile)
              - concatenated array of temperatures, 
              - temperatures at point 1, 2 and 3 (see Figure 1, Madhusudhan & 
                Seager 2009)
          T_out:    1D array of floats, temperatures concatenated for all levels
          T_l1:     1D array of floats, temperatures for layer 1
          T_l2_pos: 1D array of floats, temperatures for layer 2 inversion part 
                    (pos-increase in temperatures)
          T_l2_neg: 1D array of floats, temperatures for layer 2 negative part
                    (neg-decrease in temperatures)
          T_l3:     1D array of floats, temperatures for layer 3 
                    (isothermal part)   
          p_l1:     1D array of floats, pressures for layer 1   
          p_l2_pos: 1D array of floats, pressures for layer 2 inversion part 
                    (pos-increase in temperatures)   
          p_l2_neg: 1D array of floats, pressures for layer 2 negative part 
                    (neg-decrease in temperatures)    
          p_l3:     1D array of floats, pressures for layer 3 (isothermal part)   
          T1:       float, temperature at point 1  
          T2:       float, temperature at point 2 
          T3:       float, temperature at point 3  
     T_smooth:  1D array of floats, Gaussian smoothed temperatures, 
                no kinks on Layer boundaries 

     Notes
     -----
     See model details in Madhusudhan & Seager (2009):
     http://adsabs.harvard.edu/abs/2009ApJ...707...24M
 
     Example
     -------    
     # array of pressures, equally spaced in log space 
     p = np.array([  1.00000000e-05,   1.17680000e-05,   1.38480000e-05,
                     1.62970000e-05,   1.91790000e-05,   2.25700000e-05,
                     2.65600000e-05,   3.12570000e-05,   3.67830000e-05,
                     4.32870000e-05,   5.09410000e-05,   5.99400000e-05,
                     7.05480000e-05,   8.30280000e-05,   9.77000000e-05,
                     1.14970000e-04,   1.35300000e-04,   1.59220000e-04,
                     1.87380000e-04,   2.20510000e-04,   2.59500000e-04,
                     3.05380000e-04,   3.59380000e-04,   4.22920000e-04,
                     4.97700000e-04,   5.85700000e-04,   6.89260000e-04,
                     8.11130000e-04,   9.54540000e-04,   1.12330000e-03,
                     1.32190000e-03,   1.55560000e-03,   1.83070000e-03,
                     2.15440000e-03,   2.53530000e-03,   2.98360000e-03,
                     3.51110000e-03,   4.13200000e-03,   4.86260000e-03,
                     5.72230000e-03,   6.73410000e-03,   7.92480000e-03,
                     9.32600000e-03,   1.09740000e-02,   1.29150000e-02,
                     1.51990000e-02,   1.78860000e-02,   2.10490000e-02,
                     2.47700000e-02,   2.91500000e-02,   3.43040000e-02,
                     4.03700000e-02,   4.75080000e-02,   5.59080000e-02,
                     6.57930000e-02,   7.74260000e-02,   9.11160000e-02,
                     1.07220000e-01,   1.26180000e-01,   1.48490000e-01,
                     1.74750000e-01,   2.05650000e-01,   2.42010000e-01,
                     2.84800000e-01,   3.35160000e-01,   3.94420000e-01,
                     4.64150000e-01,   5.46220000e-01,   6.42800000e-01,
                     7.56460000e-01,   8.90210000e-01,   1.04760000e+00,
                     1.23280000e+00,   1.45080000e+00,   1.70730000e+00,
                     2.00920000e+00,   2.36440000e+00,   2.78250000e+00,
                     3.27450000e+00,   3.85350000e+00,   4.53480000e+00,
                     5.33660000e+00,   6.28020000e+00,   7.39070000e+00,
                     8.69740000e+00,   1.02350000e+01,   1.20450000e+01,
                     1.41740000e+01,   1.66810000e+01,   1.96300000e+01,
                     2.31010000e+01,   2.71850000e+01,   3.19920000e+01,
                     3.76490000e+01,   4.43060000e+01,   5.21400000e+01,
                     6.13590000e+01,   7.22080000e+01,   8.49750000e+01,
                     1.00000000e+02])

     # random values imitate DEMC
     a1 = np.random.uniform(0.2  , 0.6 )
     a2 = np.random.uniform(0.04 , 0.5 )
     p3 = np.random.uniform(0.5  , 10  )
     p2 = np.random.uniform(0.01 , 1   )
     p1 = np.random.uniform(0.001, 0.01)
     T3 = np.random.uniform(1500 , 1700)

     # generates raw and smoothed PT profile
     PT_Inv, T_smooth = PT_Inversion(p, a1, a2, p1, p2, p3, T3)

     # returns full temperature array and temperatures at every point
     T, T0, T1, T2, T3 = PT_Inv[8], PT_Inv[9], PT_Inv[10], PT_Inv[11], PT_Inv[12]

     # sets plots in the middle 
     minT= min(T0, T2)*0.75
     maxT= max(T1, T3)*1.25

     # plots raw PT profile with equally spaced points in log space
     plt.figure(1)
     plt.clf()
     plt.semilogy(PT_Inv[0], PT_Inv[1], '.', color = 'r'     )
     plt.semilogy(PT_Inv[2], PT_Inv[3], '.', color = 'b'     )
     plt.semilogy(PT_Inv[4], PT_Inv[5], '.', color = 'orange')
     plt.semilogy(PT_Inv[6], PT_Inv[7], '.', color = 'g'     )
     plt.title('Thermal Inversion Raw', fontsize=14)
     plt.xlabel('T [K]'               , fontsize=14)
     plt.ylabel('logP [bar]'          , fontsize=14)
     plt.xlim(minT  , maxT)
     plt.ylim(max(p), min(p))
     #plt.savefig('ThermInverRaw.png', format='png')
     #plt.savefig('ThermInverRaw.ps' , format='ps' )

     # plots smoothed PT profile
     plt.figure(2)
     plt.clf()
     plt.semilogy(T       , p, color = 'r')
     plt.semilogy(T_smooth, p, color = 'k')
     plt.title('Thermal Inversion Smoothed', fontsize=14)
     plt.xlabel('T [K]'                    , fontsize=14)
     plt.ylabel('logP [bar]'               , fontsize=14)
     plt.xlim(minT  , maxT)
     plt.ylim(max(p), min(p) )
     #plt.savefig('ThermInverSmoothed.png', format='png')
     #plt.savefig('ThermInverSmoothed.ps' , format='ps' )

     Revisions
     ---------
     2013-11-14  Jasmina Blecic, jasmina@physics.ucf.edu   Written by.
     2014-04-05  Jasmina Blecic, jasmina@physics.ucf.edu   Revision
                     added T3 as free parameter instead of T0
                     changed boundary condition equations accordingly
     2014-09-24  Jasmina  Updated documentation.
     2019-02-13  mhimes   Added verb argument.
     '''

     # The following set of equations derived using Equation 2
     # Madhusudhan and Seager 2009

     # Set top of the atmosphere to p0 to have easy understandable equations:
     p0 = np.amin(p)
     if verb:
          print(p0)

     # Temperature at point 2
     # Calculated from boundary condition between layer 2 and 3
     T2 = T3 - (np.log(p3/p2) / a2)**2

     # Temperature at the top of the atmosphere
     # Calculated from boundary condition between layer 1 and 2
     T0 = T2 + (np.log(p1/p2) / -a2)**2 - (np.log(p1/p0) / a1)**2 

     # Temperature at point 1
     T1 = T0 + (np.log(p1/p0) / a1)**2

     # Error message when temperatures ar point 1, 2 or 3 are < 0
     if T0<0 or T1<0 or T2<0 or T3<0:
          if verb:
               print('T0, T1, T2 and T3 temperatures are: ', T0, T1, T2, T3)
          raise ValueError('Input parameters give non-physical profile. Try again.')

     # Defining arrays of pressures for every part of the PT profile
     p_l1     = p[(np.where((p >= min(p)) & (p < p1)))]
     p_l2_pos = p[(np.where((p >= p1)  & (p < p2)))]
     p_l2_neg = p[(np.where((p >= p2)  & (p < p3)))]
     p_l3     = p[(np.where((p >= p3)  & (p <= max(p))))]

     # Sanity check for total number of levels
     check = len(p_l1) + len(p_l2_pos) + len(p_l2_neg) + len(p_l3)
     if verb:
          print('Total number of levels in p: ', len(p))
          print('\nLevels per levels in inversion case (l1, l2_pos, l2_neg, l3) are respectively: ', len(p_l1), len(p_l2_pos), len(p_l2_neg), len(p_l3))
          print('Checking total number of levels in inversion case: ', check)

     # The following set of equations derived using Equation 2
     # Madhusudhan and Seager 2009

     # Layer 1 temperatures
     T_l1     = (np.log(p_l1/p0) / a1)**2 + T0  
 
     # Layer 2 temperatures (inversion part)
     T_l2_pos = (np.log(p_l2_pos/p2) / -a2)**2 + T2
 
     # Layer 2 temperatures (decreasing part)
     T_l2_neg = (np.log(p_l2_neg/p2) / a2)**2 + T2

     # Layer 3 temperatures
     T_l3     = np.linspace(T3, T3, len(p_l3))

     # The output temperatures
     T_out = np.zeros(len(p))
     T_out[np.where((p >= p0) & (p < p1))]          = T_l1
     T_out[np.where((p >= p1) & (p < p2))]          = T_l2_pos
     T_out[np.where((p >= p2) & (p < p3))]          = T_l2_neg
     T_out[np.where((p >= p3) & (p <= np.amax(p)))] = T_l3

     # PT profile
     PT_Inver = (T_l1, p_l1, T_l2_pos, p_l2_pos, T_l2_neg, p_l2_neg, 
                 T_l3, p_l3, T_out, T0, T1, T2, T3)

     # Smoothing with Gaussian_filter1d
     sigma    = 4
     T_smooth = gaussian_filter1d(T_out, sigma, mode='nearest')
    
     return PT_Inver, T_smooth


# generated PT profile for non-inverted atmopshere
def PT_NoInversion(p, a1, a2, p1, p3, T3, verb=False):
     '''
     Calculates PT profile for non-inversion case based on Equation (2) from
     Madhusudhan & Seager 2009.
     It takes a pressure array (e.g., extracted from a pressure file), and 5
     free parameters for non-inversion case and generates non-inverted PT  
     profile. The profile is then smoothed using 1D Gaussian filter. The
     pressure array needs to be equally spaced in log space.

     Parameters
     ----------
     p:  1D array of floats
         Pressure array needs to be equally spaced in log space from bottom to top 
         of the atmosphere.
     a1: Float
         Model exponential factor in Layer 1, empirically determined to be within
         range (0.2, 0.6).
     a2: Float
         Model exponential factor in Layer 2, empirically determined to be within
         range (0.04, 0.5) 
     p1: Float
     p2: Float
         Pressure boundary between Layer  1 and 2 (in bars).
     p3: Float
         Pressure boundary between Layers 2 and 3 (in bars).
     T3: float
         Temperature in the Layer 3.
     verb: Boolean
         If True, print some info to screen.
   
     Returns
     -------
     PT_NoInver:  tupple of arrays that includes:
           - temperature and pressure arrays of every layer of the atmosphere 
             (PT profile)
           - concatenated array of temperatures, 
           - temperatures at point 1 and 3 (see Figure 1, Madhusudhan & 
             Seager 2009)
       T_out:    1D array of floats, temperatures concatenated for all levels
       T_l1:     1D array of floats, temperatures for layer 1
       T_l2_neg: 1D array of floats, temperatures for layer 2
       T_l3:     1D array of floats, temperatures for layer 3 (isothermal part)  
       p_l1:     1D array of floats, pressures for layer 1   
       p_l2_neg: 1D array of floats, pressures for layer 2     
       p_l3:     1D array of floats, pressures for layer 3 (isothermal part)     
       T1:       float, temperature at point 1  
       T3:       float, temperature at point 3  
     T_smooth:  1D array of floats, Gaussian smoothed temperatures, 
             no kinks on layer boundaries 

     Notes
     -----
     The code uses just one equation for layer 2, assuming that decrease 
     in temperature in layer 2 is the same from point 3 to point 2 
     as is from point 2 to point 1.

     Example
     -------
     # array of pressures, equally spaced in log space 
     p = np.array([  1.00000000e-05,   1.17680000e-05,   1.38480000e-05,
                  1.62970000e-05,   1.91790000e-05,   2.25700000e-05,
                  2.65600000e-05,   3.12570000e-05,   3.67830000e-05,
                  4.32870000e-05,   5.09410000e-05,   5.99400000e-05,
                  7.05480000e-05,   8.30280000e-05,   9.77000000e-05,
                  1.14970000e-04,   1.35300000e-04,   1.59220000e-04,
                  1.87380000e-04,   2.20510000e-04,   2.59500000e-04,
                  3.05380000e-04,   3.59380000e-04,   4.22920000e-04,
                  4.97700000e-04,   5.85700000e-04,   6.89260000e-04,
                  8.11130000e-04,   9.54540000e-04,   1.12330000e-03,
                  1.32190000e-03,   1.55560000e-03,   1.83070000e-03,
                  2.15440000e-03,   2.53530000e-03,   2.98360000e-03,
                  3.51110000e-03,   4.13200000e-03,   4.86260000e-03,
                  5.72230000e-03,   6.73410000e-03,   7.92480000e-03,
                  9.32600000e-03,   1.09740000e-02,   1.29150000e-02,
                  1.51990000e-02,   1.78860000e-02,   2.10490000e-02,
                  2.47700000e-02,   2.91500000e-02,   3.43040000e-02,
                  4.03700000e-02,   4.75080000e-02,   5.59080000e-02,
                  6.57930000e-02,   7.74260000e-02,   9.11160000e-02,
                  1.07220000e-01,   1.26180000e-01,   1.48490000e-01,
                  1.74750000e-01,   2.05650000e-01,   2.42010000e-01,
                  2.84800000e-01,   3.35160000e-01,   3.94420000e-01,
                  4.64150000e-01,   5.46220000e-01,   6.42800000e-01,
                  7.56460000e-01,   8.90210000e-01,   1.04760000e+00,
                  1.23280000e+00,   1.45080000e+00,   1.70730000e+00,
                  2.00920000e+00,   2.36440000e+00,   2.78250000e+00,
                  3.27450000e+00,   3.85350000e+00,   4.53480000e+00,
                  5.33660000e+00,   6.28020000e+00,   7.39070000e+00,
                  8.69740000e+00,   1.02350000e+01,   1.20450000e+01,
                  1.41740000e+01,   1.66810000e+01,   1.96300000e+01,
                  2.31010000e+01,   2.71850000e+01,   3.19920000e+01,
                  3.76490000e+01,   4.43060000e+01,   5.21400000e+01,
                  6.13590000e+01,   7.22080000e+01,   8.49750000e+01,
                  1.00000000e+02])

     # random values imitate DEMC
     a1 = np.random.uniform(0.2  , 0.6 )
     a2 = np.random.uniform(0.04 , 0.5 )
     p3 = np.random.uniform(0.5  , 10  )
     p1 = np.random.uniform(0.001, 0.01)
     T3 = np.random.uniform(1500 , 1700)

     # generates raw and smoothed PT profile
     PT_NoInv, T_smooth = PT_NoInversion(p, a1, a2, p1, p3, T3)

     # returns full temperature array and temperatures at every point
     T, T0, T1, T3 = PT_NoInv[6], PT_NoInv[7], PT_NoInv[8], PT_NoInv[9]

     # sets plots in the middle 
     minT= T0*0.75
     maxT= max(T1, T3)*1.25

     # plots raw PT profile with equally spaced points in log space
     plt.figure(3)
     plt.clf()
     plt.semilogy(PT_NoInv[0], PT_NoInv[1], '.', color = 'r'     )
     plt.semilogy(PT_NoInv[2], PT_NoInv[3], '.', color = 'b'     )
     plt.semilogy(PT_NoInv[4], PT_NoInv[5], '.', color = 'orange')
     plt.title('No Thermal Inversion Raw', fontsize=14)
     plt.xlabel('T [K]'                  , fontsize=14)
     plt.ylabel('logP [bar]'             , fontsize=14)
     plt.xlim(minT  , maxT)
     plt.ylim(max(p), min(p))
     #plt.savefig('NoThermInverRaw.png', format='png')
     #plt.savefig('NoThermInverRaw.ps' , format='ps' )

     # plots smoothed PT profile
     plt.figure(4)
     plt.clf()
     plt.semilogy(T       , p, color = 'r')
     plt.semilogy(T_smooth, p, color = 'k')
     plt.title('No Thermal Inversion Smoothed', fontsize=14)
     plt.xlabel('T [K]'                       , fontsize=14)
     plt.ylabel('logP [bar]'                  , fontsize=14)
     plt.xlim(minT  , maxT)
     plt.ylim(max(p), min(p))
     #plt.savefig('NoThermInverSmoothed.png', format='png')
     #plt.savefig('NoThermInverSmoothed.ps' , format='ps' )

     Revisions
     ---------
     2013-11-16  Jasmina   Written by.
     2014-04-05  Jasmina   Added T3 as free parameter instead of T0
                           Changed boundary condition equations accordingly
     2014-08-15  Patricio  Cleaned-up the code. Added verb argument.
     2014-09-24  Jasmina   Updated documentation.
     '''
     if verb:
       print("Pressure range: {} -- {} bar\n"
             "PT params: {} {} {} {} {}\n".format(p[0], p[-1],
                                                  a1, a2, p1, p3, T3))

     # The following set of equations derived using Equation 2
     # Madhusudhan and Seager 2009

     # Set p0 (top of the atmosphere):
     p0 = np.amin(p)

     # Calculate temperature at layer boundaries:
     T1 = T3 - (np.log(p3/p1) / a2)**2.0
     T0 = T1 - (np.log(p1/p0) / a1)**2.0

     # Error message for negative Temperatures:
     if T0 < 0 or T1 < 0 or T3 < 0:
          raise ValueError("Input parameters give non-physical profile:\n"
                     "  T0={:.1f},  T1={:.1f},  T3={:.1f}".format(T0, T1, T3))

     # Defining arrays for every part of the PT profile:
     p_l1     = p[np.where((p >= p0) & (p < p1))]
     p_l2_neg = p[np.where((p >= p1) & (p < p3))]
     p_l3     = p[np.where((p >= p3) & (p <= np.amax(p)))]

     # sanity check for total number of levels:
     check = len(p_l1) + len(p_l2_neg) + len(p_l3)
     if verb:
         print('Total number of layers: {:d}'.format(len(p)))
         print('Number of levels per Layer: Nl1={:d},  Nl2={:d}, Nl3={:d}'.format(
             len(p_l1), len(p_l2_neg), len(p_l3)))
         print('Sum of levels per layer: {:d}'.format(check))

     # Layer 1 temperatures 
     T_l1     = (np.log(p_l1/p0) / a1)**2 + T0

     # Layer 2 temperatures decreasing part
     T_l2_neg = (np.log(p_l2_neg/p1) / a2)**2 + T1

     # Layer 3 temperatures
     T_l3     = np.linspace(T3, T3, len(p_l3))

     # The output temperatures
     T_out = np.zeros(len(p))
     T_out[np.where((p >= p0) & (p < p1))]          = T_l1
     T_out[np.where((p >= p1) & (p < p3))]          = T_l2_neg
     T_out[np.where((p >= p3) & (p <= np.amax(p)))] = T_l3

     # PT profile info:
     PT_NoInver = (T_l1, p_l1, T_l2_neg, p_l2_neg, T_l3, p_l3, 
                   T_out, T0, T1, T3)

     # Smoothed PT profile:
     sigma    = 4
     T_smooth = gaussian_filter1d(T_out, sigma, mode='nearest')

     return PT_NoInver, T_smooth


def PT_line(pressure, kappa,  gamma1, gamma2, alpha, beta, 
            R_star,   T_star, T_int,  sma,    grav,  T_int_type):
  '''
  Generates a PT profile based on input free parameters and pressure array.
  If no inputs are provided, it will run in demo mode, using free
  parameters given by the Line 2013 paper and some dummy pressure
  parameters.

  Inputs
  ------
  pressure: 1D float ndarray
     Array of pressure values in bars.
  kappa : float, in log10. Planck thermal IR opacity in units cm^2/gr
  gamma1: float, in log10. Visible-to-thermal stream Planck mean opacity ratio.
  gamma2: float, in log10. Visible-to-thermal stream Planck mean opacity ratio.
  alpha : float.           Visible-stream partition (0.0--1.0).
  beta  : float.           A 'catch-all' for albedo, emissivity, and day-night
                           redistribution (on the order of unity)
  R_star: Float
     Stellar radius (in meters).
  T_star: Float
     Stellar effective temperature (in Kelvin degrees).
  T_int:  Float
     Planetary internal heat flux (in Kelvin degrees).
  sma:    Float
     Semi-major axis (in meters).
  grav:   Float
     Planetary surface gravity (at 1 bar) in cm/second^2.
  T_int_type: string.
     Method for determining `T_int`: 'const' (for a supplied constant value)
                                     'thorngren' (to use Thorngren et al. 2019)

  Returns
  -------
  T: temperature array

  Example:
  --------
  >>> import PT as pt
  >>> import scipy.constants as sc
  >>> import matplotlib.pyplot as plt
  >>> import numpy as np

  >>> Rsun = 6.995e8 # Sun radius in meters

  >>> # Pressure array (bars):
  >>> p = np.logspace(2, -5, 100)

  >>> # Physical (fixed for each planet) parameters:
  >>> Ts = 5040.0        # K
  >>> Ti =  100.0        # K
  >>> a  = 0.031 * sc.au # m
  >>> Rs = 0.756 * Rsun  # m
  >>> g  = 2192.8        # cm s-2

  >>> # Fitting parameters:
  >>> kappa  = -1.5   # log10(3e-2)
  >>> gamma1 = -0.8   # log10(0.158)
  >>> gamma2 = -0.8   # log10(0.158)
  >>> alpha  = 0.5
  >>> beta   = 1.0
  >>> T0 = pt.PT(p, kappa, gamma1, gamma2, alpha, beta, Rs, Ts, Ti, a, g)

  >>> plt.figure(1)
  >>> plt.clf()
  >>> plt.semilogy(T0, p, lw=2, color="b")
  >>> plt.ylim(p[0], p[-1])
  >>> plt.xlim(800, 2000)
  >>> plt.xlabel("Temperature  (K)")
  >>> plt.ylabel("Pressure  (bars)")

  Developers:
  -----------
  Madison Stemm      astromaddie@gmail.com
  Patricio Cubillos  pcubillos@fulbrightmail.org

  Modification History:
  ---------------------
  2014-09-12  Madison   Initial version, adapted from equations (13)-(16)
                        in Line et al. (2013), Apj, 775, 137.
  2014-12-10  patricio  Reviewed and updated code.
  2015-01-22  patricio  Receive log10 of free parameters now.
  2019-02-13  mhimes    Replaced `params` arg with each parameter for 
                        consistency with other PT models
  2019-09-10  mhimes    Added T_int calculation from Thorngren et al. (2019)
  '''
  # Convert kappa, gamma1, gamma2 from log10
  kappa  = 10**(kappa )
  gamma1 = 10**(gamma1)
  gamma2 = 10**(gamma2)

  if T_int_type == 'thorngren':
      # Planetary internal temperature (Thorngren et al. 2019)
      # Hard-coded values are fitted parameters!
      T_eq  = (R_star/(2.0*sma))**0.5 * T_star
      F     = 4.0 * sc.Stefan_Boltzmann * T_eq**4
      T_int = 1.24 * T_eq * np.exp(-(np.log(F) - 0.14)**2 / 2.96)

  # Stellar input temperature (at top of atmosphere):
  T_irr = beta * (R_star / (2.0*sma))**0.5 * T_star

  # Gray IR optical depth:
  tau = kappa * (pressure*1e6) / grav # Convert bars to barye (CGS)

  xi1 = xi(gamma1, tau)
  xi2 = xi(gamma2, tau)

  # Temperature profile (Eq. 13 of Line et al. 2013):
  temperature = (0.75 * (T_int**4 * (2.0/3.0 + tau) +
                         T_irr**4 * (1-alpha) * xi1 +
                         T_irr**4 * alpha     * xi2 ) )**0.25

  return temperature


def PT_iso(p, T):
  """
  Returns an isothermal atmosphere at given pressures.
  Parameters
  ----------
  pressure: 1D float ndarray
     Array of pressure values in bars.
  Returns
  -------
  1D temperature array

  Revisions
  ---------
  2018-05-24    Ryan        Initial implementation
  """
  return np.ones(len(p)) * T


def xi(gamma, tau):
  """
  Calculate Equation (14) of Line et al. (2013) Apj 775, 137

  Parameters:
  -----------
  gamma: Float
     Visible-to-thermal stream Planck mean opacity ratio.
  tau: 1D float ndarray
     Gray IR optical depth.

  Modification History:
  ---------------------
  2014-12-10  patricio  Initial implemetation.
  """
  return (2.0/3) * \
         (1 + (1./gamma) * (1 + (0.5*gamma*tau-1)*np.exp(-gamma*tau)) +
          gamma*(1 - 0.5*tau**2) * sp.expn(2, gamma*tau)              )

def PT_adiabatic(p, T0, gamma, logp0):
     '''
     Calculates an adiabatic temperature profile. Currently
     makes naive assumptions about temperature, molecular mass,
     and gravity variations with altitude.
     '''
     p0 = 10**logp0
     T = T0 / (1 + (gamma - 1) / gamma * np.log(p0 / p))
     return T


def PT_piette(p, T0, dTbot_32, dT32_10, dT10_0, dT0_1, dT1_01, dT01_001, dT001_top):
    '''
    Non-inverted atmospheric profile based on the "SPT" model of 
    Piette & Madhusudhan (2020), modified to work on arbitrary pressure grids.

    Inputs
    ------
    p: 1D float ndarray
        Atmospheric pressure profile (in bar).
    T0: float
        Temperature at pressure layer closest to 3.2 bar. Reference temperature 
        for the other params.
    dTbot_32: float
        Temperature difference between the bottom of the atmosphere and the 
        pressure layer closest to 32 bar
    dT32_10: float
        Temperature difference between the pressure layer closest to 32 bar and 
        the pressure layer closest to 10 bar
    dT10_0: float
        Temperature difference between the pressure layer closest to 10 bar and 
        the pressure layer closest to 3.2 bar
    dT0_1: float
        Temperature difference between the pressure layer closest to 3.2 bar and
        the pressure layer closest to 1 bar
    dT1_01: float
        Temperature difference between the pressure layer closest to 1 bar and 
        the pressure layer closest to 0.1 bar
    dT01_001: float
        Temperature difference between the pressure layer closest to 100 mbar 
        and the pressure layer closest to 10 mbar
    dT001_top: float
        Temperature difference between the layer closest to 10 mbar and the 
        top of the atmosphere
    '''
    T = np.zeros(p.shape)
    # Indices of pressure layers
    ilaytop = np.argmin(p)
    ilay001 = np.argmin(np.abs(p - 0.01))
    ilay01  = np.argmin(np.abs(p - 0.1))
    ilay1   = np.argmin(np.abs(p - 1))
    ilay0   = np.argmin(np.abs(p - 3.2))
    ilay10  = np.argmin(np.abs(p - 10))
    ilay32  = np.argmin(np.abs(p - 32))
    ilaybot = np.argmax(p)
    # Temperatures at those layers
    T[ilay0]   = T0
    T[ilay10]  = T0         + dT10_0
    T[ilay32]  = T[ilay10]  + dT32_10
    T[ilaybot] = T[ilay32]  + dTbot_32
    T[ilay1]   = T0         - dT0_1
    T[ilay01]  = T[ilay1]   - dT1_01
    T[ilay001] = T[ilay01]  - dT01_001
    T[ilaytop] = T[ilay001] - dT001_top
    # Linear spline interpolate
    import scipy.interpolate as si
    ilays = np.array([ilaytop, ilay001, ilay01, ilay1, ilay0, ilay10, ilay32, ilaybot])
    rep = si.splrep(np.log10(p[ilays]), T[ilays], k=1)
    T   = si.splev(np.log10(p), rep)
    # Gaussian smoothing - sigma is 0.3 dex, Piette & Madhusudhan (2020)
    sig = 0.3 / np.abs(np.log10(p)[0] - np.log10(p)[1])
    return gaussian_filter1d(T, sigma=sig, mode='nearest')


def PT_generator(p, free_params, PTfunc, PTargs=None):
  '''
  Wrapper to generate an inverted or non-inverted temperature and pressure
  profile.

  Parameters:
  -----------
  p: 1D float ndarray
     Atmospheric pressure profile (in bar).
  free_params: 1D array of floats
     Set of 6 PT parameters (a1, a2, p1, p2, p3, T3) as
     in Madhusudhan & Seager (2009).
  PTfunc: pointer to function
     Determines the method of evaluating the PT profile's temperature
  PTargs: list
     If not None, passed to `PTfunc`

  Returns
  -------
  T_smooth: 1D floats array
     Temperature array (in K).

  Modification History:
  ---------------------
  2013-11-23  Jasmina   Written by.
  2014-04-05  Jasmina   Set T3 as free parameter instead of T0.
  2014-04-11  Patricio  Modified to work with the BART code.
  2014-12-27  Patricio  Added profile from Line et al. (2013).
  2019-02-13  mhimes    Refactored to take a pointer to a function.
  '''
  # Pass extra args if given
  if PTargs is None:
    Temp = PTfunc(p, *free_params)
  else:
    Temp = PTfunc(p, *(list(free_params) + PTargs))

  # madhu profiles have 2 returns, the params and the temps
  if type(Temp) == tuple:
    Temp = Temp[-1] # only take the temps
    
  return Temp


# plots PT profiles
def plot_PT(p, PT, T_smooth, MadhuPT):
     '''
     This function plots two figures:
     1.
     Madhu's PT profile for inversion or non-inversion case. It uses returned arrays
     from the PT generator, to plot the PT curves part by part.
     2. 
     Smoothed PT profile without kinks on layer transitions.

     Parameters
     ----------
     p: 1D array of floats 
         Pressure in the atmosphere.
     PT: tuple of 1D arrays of floats
         Tuple containing temperatures and pressures for every layer in
         the atmosphere.
     T_smooth: 1D array of floats
         Temperatures smoothed with Gaussian.
     MadhuPT: string
         Defines inversion or non-inversion case.
      
     Returns
     -------
     None

     Revisions
     ---------
     2013-11-20 Jasmina   Written by.
     2014-04-05 Jasmina   Excluded free_params as arguments. Added T0 to be read by 
                          PT generator.
     2014-07-24 Jasmina   Integrated general plotting function.
     2014-09-24 Jasmina   Updated documentation.
     ''' 
    
     if MadhuPT == 'MadhuPT_Inv':
          # Takes temperatures from PT generator
          T, T0, T1, T2, T3 = PT[8], PT[9], PT[10], PT[11], PT[12]

          # Sets plots in the middle 
          minT= min(T0, T2)*0.75
          maxT= max(T1, T3)*1.25

          # Plots raw PT profile
          plt.figure(1)
          plt.clf()
          plt.semilogy(PT[0], PT[1], '.', color = 'r'     )
          plt.semilogy(PT[2], PT[3], '.', color = 'b'     )
          plt.semilogy(PT[4], PT[5], '.', color = 'orange')
          plt.semilogy(PT[6], PT[7], '.', color = 'g'     )
          plt.title('Thermal Inversion Raw', fontsize=14)
          plt.xlabel('T [K]'               , fontsize=14)
          plt.ylabel('logP [bar]'          , fontsize=14)
          plt.xlim(minT  , maxT)
          plt.ylim(max(p), min(p))
          #plt.savefig('ThermInverRaw.png', format='png')
          #plt.savefig('ThermInverRaw.ps' , format='ps' )

          # Plots Gaussian smoothing
          plt.figure(2)
          plt.clf()
          plt.semilogy(T_smooth, p, '-', color = 'b', linewidth=1)
          plt.title('Thermal Inversion Smoothed with Gaussian', fontsize=14)
          plt.xlabel('T [K]'     , fontsize=14)
          plt.ylabel('logP [bar]', fontsize=14)
          plt.xlim(minT  , maxT)
          plt.ylim(max(p), min(p))

          # Overplots with the raw PT
          plt.semilogy(T, p, color = 'r')

          return

     elif MadhuPT == 'MadhuPT_NoInv':  
          # Takes temperatures from PT generator
          T, T0, T1, T3 = PT[6], PT[7], PT[8], PT[9]

          # Sets plots in the middle 
          minT= T0*0.75
          maxT= max(T1, T3)*1.25

          # Plots raw PT profile
          plt.figure(3)
          plt.clf()
          plt.semilogy(PT[0], PT[1], '.', color = 'r'     )
          plt.semilogy(PT[2], PT[3], '.', color = 'b'     )
          plt.semilogy(PT[4], PT[5], '.', color = 'orange')
          plt.title('No Thermal Inversion Raw', fontsize=14)
          plt.xlabel('T [K]'                  , fontsize=14)
          plt.ylabel('logP [bar]'             , fontsize=14)
          plt.xlim(minT  , maxT)
          plt.ylim(max(p), min(p))
          #plt.savefig('NoThermInverRaw.png', format='png')
          #plt.savefig('NoThermInverRaw.ps' , format='ps' )

          # Plots Gaussian smoothing
          plt.figure(4)
          plt.clf()
          plt.semilogy(T_smooth, p, '-', color = 'b', linewidth=1)
          plt.title('No Thermal Inversion Smoothed with Gaussian', fontsize=14)
          plt.xlabel('T [K]'     , fontsize=14)
          plt.ylabel('logP [bar]', fontsize=14)
          plt.xlim(minT  , maxT)
          plt.ylim(max(p), min(p))

          # Overplots with the raw PT
          plt.semilogy(T, p, color = 'r')

          return


