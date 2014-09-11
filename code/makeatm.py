#! /usr/bin/env python

# ******************************* END LICENSE *******************************
# Thermal Equilibrium Abundances (TEA), a code to calculate gaseous molecular
# abundances for hot-Jupiter atmospheres under thermochemical equilibrium
# conditions.
#
# This project was completed with the support of the NASA Earth and Space 
# Science Fellowship Program, grant NNX12AL83H, held by Jasmina Blecic, 
# PI Joseph Harrington. Lead scientist and coder Jasmina Blecic, 
# assistant coder for the first pre-release Oliver M. Bowman.  
# 
# Copyright (C) 2014 University of Central Florida.  All rights reserved.
# 
# This is a test version only, and may not be redistributed to any third
# party.  Please refer such requests to us.  This program is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.
# 
# We welcome your feedback, but do not guarantee support.  Please send
# feedback or inquiries to both:
# 
# Jasmina Blecic <jasmina@physics.ucf.edu>
# Joseph Harrington <jh@physics.ucf.edu>
# 
# or alternatively,
# 
# Jasmina Blecic and Joseph Harrington
# UCF PSB 441
# 4000 Central Florida Blvd
# Orlando, FL 32816-2385
# USA
# 
# Thank you for testing TEA!
# ******************************* END LICENSE *******************************

import os
import numpy as np
import scipy.constants as sc

import PT as pt
import reader as rd

"""
  This code produces a pre-atm file in the format that TEA can read. 

  Functions:
  ----------
  get_g:
     Extract the surface gravity from the given TEP file.
  radpress:
      Calculate the radii for each layer given planetary surface gravity and
      the pressure, temperature, and mean-molecular-mass arrays.
  make_preatm:
      Write a pre-atm file.

  Notes:
  ------
  Possible user errors in making pre-atm file that conflict with TEA:
   - H, and He elements as in_elem must be included
   - H_g, He_ref and H2_ref species in out_spec must be included
   - all in_elem must be included in the list of out_species with their states
   - use elements names as they appear in the periodic table
   - use species names as readJANAF.py produces them. See sorted version of the
     conversion-record.txt for the correct names of the species. 
   - If the code stalls at the first iteration of the first temperature, check 
     if all elements that appear in the species list are included with their 
     correct names.

  Developers:
  -----------
  Jasmina Blecic     jasmina@physics.ucf.edu
  Oliver Bowman      email@oliver.ucf               FINDME
  Patricio Cubillos  pcubillos@fulbrightmail.org

  Modification History:
  ---------------------
  2014-06-01  Oliver, Jasmina  First releaste as part of the TEA package.
  2014-07-14  Jasmina          Added get_g, radpress, and read_press_file
                               functions.
  2014-08-15  Patricio         Removed read_press_file, use function from PT.py
"""

# reads the tep file and calculates surface gravity
def get_g(tepfile):
    '''
    Calculates planetary surface gravity. Calls tep reader and 
    gets data needed for calculation (g = G*M/r^2). Returns
    surface gravity and surface radius.

    Parameters
    ----------
    tepfile: tep file, ASCII file

    Returns
    -------
    g: Float
       The planet surface gravity in m/s^2.
    Rp: Float
       The planet radius in km.

    Revisions
    ---------
    2014-06-11  Jasmina   Initial version
    2014-08-15  Patricio  Updated docstring.  Get NASA Jupiter values.
    '''

    # Jupiter mass and radius:
    # Source: http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
    Mjup = 1.8983e27      # kg
    Rjup = 71492.0 * 1e3  # m

    # Open tepfile to read and get data:
    tep = rd.File(tepfile)

    # Get planet mass in Mjup:
    planet_mass = np.float(tep.getvalue('Mp')[0])

    # Get planet radius in units of Rjup:
    planet_rad = np.float(tep.getvalue('Rp')[0])

    # Mass and radius of the planet in MKS units:
    Mp = planet_mass * Mjup  # kg
    Rp = planet_rad  * Rjup  # m

    # Calculate the planet surface gravity:
    # g = G*Mp/ Rp^2 in m/s^2
    g = sc.G * Mp / (Rp**2)

    # convert Rp back to km
    Rp = Rp / 1000.

    return g, Rp


# calculate radii for each pressure in the atmosphere
def radpress(tepfile, temp, mu, pres):
    '''
    Given a pressure in bar, temperature in K, mean molecular
    mass in g/mol, the function calculates radii in km for each 
    layer in the atmosphere. It declares constants, calls get_g()
    function, allocates array for radii, and calculates radii for
    each pressure. The input pressure and temperature arrays must
    in descending order, because reference radius Rp is placed on
    the bottom of the atmosphere where pressure is largest. 
    All upper radii are calculated based on the reference level. 
    The final radii array is given in km, and converted back to 
    ascending order for pre-atm file.

    Parameters
    ----------
    g: float
       Surface gravity in m/s^2.
    R0: float
       Surface radius in km.
    T: array of floats
       Array containing temperatures in K for each layer in the atmosphere.
    mu: array of floats
       Array containing mean molecular mass in g/mol for each layer 
       from in the atmosphere.
    P: array of floats
       Array containing pressures in bar for each layer in the atmosphere.

    Returns
    -------
    rad: array of floats
        Array containing radii in km for each pressure layer in the atmosphere.

    Revisions
    ---------
    2014-07-04 0.1  Jasmina Blecic, jasmina@physics.ucf.edu   Original version

    '''

    # Define physical constants
    Boltzmann = 1.38065e-23 # J/K == (N*m) / K = kg m/s^2 * m / K
    Avogadro  = 6.022e23    # 1/mol

    # Calculate surface gravity
    g, Rp = get_g(tepfile)

    # Number of layers in the atmosphere
    n = len(pres)

    # Allocate array of radii
    rad = np.zeros(n)

    # Reverse the order of pressure and temperature array, so it starts from 
    #         the bottom of the atmosphere where first radii is defined as Rp
    pres = pres[::-1]
    temp = temp[::-1]

    # Andrews:"Introduction to Atmospheric Physics" page 26
    # rad2 - rad1 = (R * T1)/g * ln(p1/p2)
    # R = R*/mu, R* - universal gas constant, mu - mean molecular mass
    # R* = Avogadro * Boltzmann
    # rad2 = rad1 + (Avogadro * Bolztmann * T1) / (g * mu1) * ln(p1/p2)
    # to convert from g to kg, mu[i] need to be multiplied with 1000
    for i in np.arange(n):
        if i == 0:
            rad[0] = Rp
        else:
            rad[i] = rad[i-1] + (Avogadro * Boltzmann * temp[i-1]  * np.log(pres[i-1]/pres[i])) / (g * mu[i-1])

    # Reverse the order of calculated radii to write them in the right order
    #         in pre-atm file
    rad = rad[::-1]

    return rad


def make_preatm(tepfile, press_file, abun_file, in_elem, out_spec,
                pre_atm, Temp):
  '''
  This code produces a pre-atm file in the format that TEA can read it. It  
  defines the directory with user inputs, then it reads the pressure file and
  elemental dex abundance data file (default: abundances.txt). It trims the 
  abundance data to the elements of interest, and takes the column with the
  weights information to calculate mean molecular weight, based on an assumption
  that 85% of the atmosphere is filled with H2 and 15% with He. Then, it converts
  dex abundances of all elements of interest to number density and divides them 
  by the sum of all number densities in the mixture to get fractional abundances.
  It calls the get_Teff() function to calculate effective temperature of the 
  planer, produces initial PT free parameters for DEMC, plots the figures for
  user to check the output, call radpress() to calculate radii, and writes
  the data (radii, pressure, temperature, elemental abundances) into a pre-atm
  file.

  Parameters
  ----------
  tepfile: String
  press_file: String
  abun_file: String
  in_elem: String
  out_spec: String
  pre_atm: String
  Temp: 1D float ndarray

  Returns
  -------
  None

  Revisions
  ---------
  2014-07-12  Jasmina   Initial version
  2014-08-15  Patricio  Removed PT re-calculation, use input T array.
  '''

  # Read pressure data:
  pres = pt.read_press_file(press_file)

  # Number of layers in the atmosphere
  n_layers = len(pres)

  # Read abundance data and convert to array:
  f = open(abun_file, 'r')
  abundata = []
  for line in f.readlines():
      if line.startswith('#'):
          continue
      else:
          l = [value for value in line.split()]
          abundata.append(l)
  abundata = np.asarray(abundata)
  f.close()

  # Trim abundata to elements we need
  in_elem_split = in_elem.split(" ")
  nelem  = np.size(in_elem_split)
  data_slice = np.zeros(abundata.shape[0], dtype=bool)
  for i in np.arange(nelem):
    data_slice += (abundata[:,1] == in_elem_split[i])

  # List of elements of interest and their corresponding data
  abun_trim = abundata[data_slice]

  # Take data and create list:
  out_elem = abun_trim[:,1].tolist()

  # Convert strings to floats:
  out_dex  = map(float, abun_trim[:,2])

  # Convert dex exponents to elemental counts:
  out_num  = 10**np.array(out_dex)

  # Get fractions of elemental counts to the 
  #     total sum of elemental counts in the mixture
  out_abn  = (out_num / np.sum(out_num)).tolist()

  # Calculate mean molecular weight guess from H2 and He only:
  H_weight  = np.float(abundata[1][4])
  He_weight = np.float(abundata[2][4])
  mu = 2.0*0.85*H_weight + 0.15*He_weight
  mu_array = np.linspace(mu, mu, n_layers) 

  # Call radpress function to calculate radii:
  rad = radpress(tepfile, Temp, mu_array, pres)
  
  # Write pre-atm file:
  f = open(pre_atm, 'w')

  # Pre-atm header with basic instructions:
  header = ("# This is a TEA pre-atmosphere input file.\n"
     "# TEA accepts a file in this format to produce species abundances as\n"
     "# a function of pressure and temperature.\n"
     "# Output species must be added in the line immediately following the \n"
     "# FINDSPEC marker and must be named to match JANAF converted names.\n"
     "# Units: radius (km), pressure (bar), temperature (K), abundance\n"
     "# (number density)\n\n")
  f.write(header)

  f.write('#FINDSPEC\n' + out_spec + '\n\n')

  f.write('#FINDTEA\n')
  # Column's header:
  f.write("#Radius        Pressure      Temp     " +
          "".join(["{:<18s}".format(elem) for elem in out_elem]) + "\n")

  # For each layer:
  for i in np.arange(n_layers):
    # Radius, Pressure, and Temperature:
    f.write("{:12.3f} {:14.6e} {:8.2f}  ".format(rad[i], pres[i], Temp[i]))
    # Elemental abundance list:
    f.write("  ".join(["{:16.10e}".format(abun) for abun in out_abn]) + "\n")
  f.close()


def get_press(atmfile):
  """
  Reads an atmospheric file made in Transit (Patricio Rojo) format.
  The function takes the column with pressures and converts it to floats.

  Parameters:
  -----------
  atmfile: String
     Name of TEA (pre)atmospheric ASCII file.

  Returns:
  --------
  pressure: 1D double ndarrya
     Array of pressures listed in the atmospheric file.

  Notes:
  ------
  Per-layer atmospheric data starts two lines after the header line:
    '#FINDTEA\n'
  The pressure is given in the second column.

  Modification History:
  ---------------------
  2013-11-17  Jasmina   Initial version.
  2014-08-19  Patricio  Adapted to TEA output.
  """
  # Open the atmospheric file and read:
  f = open(atmfile, 'r')
  lines = np.asarray(f.readlines())
  f.close()

  # Find the line where the layers info begins:
  start = np.where(lines == "#FINDTEA\n")[0][0] + 2
  # Allocates array of pressures:
  ndata = len(lines) - start
  pressure = np.zeros(ndata, np.double)

  # Extract pressure data:
  for i in np.arange(ndata):
    pressure[i] = lines[start + i].strip().split()[1]

  return pressure


def readatm(atmfile):
  """
  Reads an atmospheric file made in Transit (Patricio Rojo) format.
  The function takes the column with pressures and converts it to floats.

  Parameters:
  -----------
  atmfile: String
     Name of TEA (pre)atmospheric ASCII file.

  Returns:
  --------
  pressure: 1D double ndarrya
     Array of pressures listed in the atmospheric file.

  Notes:
  ------
  Per-layer atmospheric data starts two lines after the header line:
    '#FINDTEA\n'
  The pressure is given in the second column.

  Modification History:
  ---------------------
  2014-09-05  Patricio  Modified from getpressure.
  """
  # Open the atmospheric file and read:
  f = open(atmfile, 'r')
  lines = np.asarray(f.readlines())
  f.close()

  # Get the molecules:
  imol = np.where(lines == "#FINDSPEC\n")[0][0] + 1
  molecules = lines[imol].split()

  # Find the line where the layers info begins:
  start = np.where(lines == "#FINDTEA\n")[0][0] + 2
  # Number of layers:
  ndata = len(lines) - start
  # Number of abundances (elements per line except Rad, Press, and T).
  nabun = len(lines[start].split()) - 3
  abundances = np.zeros((ndata, nabun), np.double)
  pressure = np.zeros(ndata, np.double)

  # Extract pressure data:
  for i in np.arange(ndata):
    line = lines[start+i].strip().split()
    pressure[i]   = line[1]
    abundances[i] = line[3:]

  return molecules, pressure, abundances


def reformat(atmfile):
  """
  Re-format the TEA output file to work for transit.

  Parameters:
  -----------
  atmfile: String
     Name of TEA atmospheric ASCII file.

  Modification History:
  ---------------------
  2014-08-25  patricio  Initial implementation
  """
  # Open the atmospheric file and read:
  f = open(atmfile, 'r')
  lines = np.asarray(f.readlines())
  f.close()

  # Find the Molecule names and remove their suffix:
  imol = np.where(lines == "#FINDSPEC\n")[0][0] + 1
  molecules = lines[imol].split()
  for m in np.arange(len(molecules)):
    molecules[m] = molecules[m].partition('_')[0]
  lines[imol] = " ".join(molecules) + "\n"

  # Repeat for the column headers:
  start = np.where(lines == "#FINDTEA\n")[0][0] + 1
  molecules = lines[start].split()
  for m in np.arange(len(molecules)):
    molecules[m] = "{:10s}".format(molecules[m].partition('_')[0])
  lines[start] = " ".join(molecules) + "\n"

  lines = list(lines)
  # Reverse order to have layers from bottom to top:
  datalines = lines[start+1:]
  datalines.reverse()
  lines = lines[:start+1] + datalines

  # Add values' units:
  # Abundance by number density:
  lines.insert(imol-2, "q number\n")
  # Pressure in bars:
  lines.insert(imol-2, "up 1e6\n")
  # Radius in kilometers:
  lines.insert(imol-2, "\n#Values units:\nur 1e5\n")

  # Save file with the updated lines:
  f = open(atmfile + "_m", 'w')
  f.writelines(lines)
  f.close()
