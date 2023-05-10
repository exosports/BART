# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

import os
import six
import shutil
import re
import numpy as np
import scipy.constants as sc
from scipy.interpolate import interp1d

import PT        as pt
import reader    as rd
import constants as c

"""
    This code produces a pre-atm file in the format that TEA can read it.
    It, then, reads the final TEA output, adds radius array, and reformats
    the file for Transit to use it.

    Functions
    ---------
    get_g:
          Extracts the surface gravity from the given TEP file.
    readatm:
          Reads atmospheric file made by TEA.
    radpress:
          Calculates the radii for each layer given planetary surface gravity
          and the pressure, temperature, and mean-molecular-mass arrays.
    readAbun:
          Trims the elemental data of interest from 'abundances.txt' file.
    stoich:
          Calculates stoichiometric values of output species.
    mean_molar_mass:
          Calculates mean molecular mass of all output species.
    makeRadius:
          Rewrites TEA atmospheric file with new column for radii.
    make_preatm:
          Writes a pre-atmospheric file.
    reformat:
          Reformats final atmospheric file for Transit code.

    Notes
    -----
    - The user should have unique T-P values for each line in the
      PT-profile file before running makeatm.py, otherwise
      each repeated layer will be overwritten and the final atmosphere
      file will end up with fewer lines.
    - input_elem must have names as they appear in the
      conversion_record_sorted.txt file for the correct names of
      the species.
    - Should the code stall on the first iteration of the first
      temperature, check if all elements that appear in the species
      list are included with their correct names.
    - Elemental hydrogen (H) and helium (He) must be included in
      in_elem for hot-Jupiter atmosphere calculations. Similarly,
      the H_g, He_ref and H2_ref species must also appear in out_spec
      for these calculations.

    Developers
    ----------
    Jasmina Blecic     jasmina@physics.ucf.edu
    Oliver Bowman      obowman@knights.ucf.edu
    Patricio Cubillos  pcubillos@fulbrightmail.org

    Revisions
    ---------
    2014-07-14 Jasmina   Added get_g, radpress, and read_press_file
                         functions, and rewrote and renames
                         make_preatm() for BART project.
    2014-08-15 Patricio  Removed read_press_file, use function from PT.py.
                         Modified read_press_file to readatm().
                         Added reformat()
    2014-09-15 Jasmina   Added makeRadius(), readAbun(),
                         mean_molar_mass(), stoich(). Modified readatm(),
                         make_preatm(), reformat().
    2014-09-27 Jasmina   Modified radpress() to calculate radius from the
                         referenced surface pressure of 1bar.
                         Cleaned documentation.
    2014-10-02  Jasmina  Added option to read atmfile with or without radius
                         array in it.
    2014-11-06  Jasmina  Minor adjustments for proper TEA execution.
    2015-05-03  jasmina  Corrected atm header.
"""

def read_eabun(solabun):
  """
  Extract the Solar elemental-abundances information from file.

  Parameters:
  -----------
  efile: String
    Name of the elemental abundances file.

  Returns:
  --------
  index: 1D integer ndarray
     Ordinal index.
  symbol: 1D string ndarray
     Elemental chemical symbol.
  dex: 1D float ndarray
     Logarithmic number-abundance, scaled to log(H) = 12.
  name: 1D string ndarray
     Element names.
  mass: 1D float ndarray
     Elemental mass in amu.

  Modification History:
  --------------------
  2014-07-12  Jasmina   Written by.
  2014-08-15  Patricio  Rewrote data handling. Updated data strings.
  2014-09-24  Jasmina   Updated documentation.
  2015-03-06  patricio  Reworked code from makeAbun function.
  """

  # Read the elemental-abundances file:
  f = open(solabun, 'r')
  lines = f.readlines()
  f.close()

  # Count the number of elements:
  nelements = len(lines)
  for line in lines:
    if line.startswith("#"):
      nelements -= 1

  # Allocate arrays to put information:
  index  = np.zeros(nelements, int)
  symbol = np.zeros(nelements, '|S2' if six.PY2 else 'U2')
  dex    = np.zeros(nelements, np.double)
  name   = np.zeros(nelements, '|S20' if six.PY2 else 'U20')
  mass   = np.zeros(nelements, np.double)

  # Store data into the arrays:
  i = 0
  for line in lines:
    if not line.startswith("#"):
      index[i], symbol[i], dex[i], name[i], mass[i] = line.strip().split()
      i += 1
  return index, symbol, dex, name, mass


# reads the tep file and calculates surface gravity
def get_g(tepfile):
  '''
  Calculates planetary surface gravity. Calls tep reader and
  gets data needed for calculation (g = G*M/r^2). Returns
  surface gravity and surface radius.

  Parameters
  ----------
  tepfile: String
     Name of the tep ASCII file.

  Returns
  -------
  g: Float
     The planet surface gravity in m/s^2.
  Rp: Float
     The planet radius in km.

  Revisions
  ---------
  2014-06-11  Jasmina   Written by.
  2014-08-15  Patricio  Updated docstring.  Got NASA Jupiter values.
  '''
  # Open tepfile to read and get data:
  tep = rd.File(tepfile)

  # Get planet mass in kg:
  Mplanet = np.float(tep.getvalue('Mp')[0]) * c.Mjup
  # Get planet radius in m:
  Rplanet = np.float(tep.getvalue('Rp')[0]) * c.Rjup

  # Calculate the planet surface gravity in m/s^2:
  g = sc.G * Mplanet / (Rplanet**2)
  # Return Rp in km:
  Rp = Rplanet / 1000.0

  return g, Rp


def radpress(pressure, temperature, mu, p0, R0, g0):
  """
  Calculate the radii for the atmospheric layers applying the
  hidrostatic-equilibrium equation with the constraint that:
  radius(p0) = R0.

  Parameters:
  -----------
  pressure: 1D float ndarray
     Atmospheric layers' pressure in bars.
  temperature: 1D float ndarray
     Atmospheric layers' temperature in K.
  mu: 1D float ndarray
     Atmospheric layers' mean molecular mass.
  p0: Float
     Reference pressure level, i.e. R(p0) = R0, in bars.
  R0: Float
     Reference radius level in km.
  g0: Float
     Atmospheric gravity in m s-2.
  """
  # Number of layers in the atmosphere:
  n = len(pressure)

  # Allocate array of radii:
  rad = np.zeros(n)
  g   = np.zeros(n)

  # Interpolate temp and mu in lin-log space (1bar)
  interPT = interp1d(np.log10(pressure), temperature)
  intermu = interp1d(np.log10(pressure), mu)

  # Get temp and mu at surface pressure (1bar)
  try:
    temp0 = interPT(np.log10(p0))
    mu0   = intermu(np.log10(p0))
  except IOError:
    print("Referenced surface pressure of {:.3e} bar is not in the "
          "range of pressures: [{}, {}] bar.".
           format(p0, np.amin(pressure), np.amax(pressure)))

  # Return back to desending order for radius calculation
  press = pressure   [::-1]
  temp  = temperature[::-1]
  mu    = mu         [::-1]

  # Find the index of the closest pressure point to p0:
  idx = np.argmin(np.abs(press - p0))

  # If p0 is not in press:
  if press[idx] != p0:
      # If the point is above p0:
      if press[idx] > p0:
          rad[idx] = R0 + 0.5 * (temp[idx] / mu[idx] + temp0 / mu0)    \
                     * (sc.Avogadro * sc.k * np.log(p0/press[idx]) / g0)
      # If the point is below p0:
      else:
          rad[idx] = R0 - 0.5 * (temp[idx] / mu[idx] + temp0 / mu0)    \
                     * (sc.Avogadro * sc.k * np.log(press[idx]/p0) / g0)
      g[idx] = g0 * R0**2 / rad[idx]**2

  else:
      rad[idx] = R0
      g[idx]   = g0

  # Calculate radius below p0:
  for i in reversed(np.arange(idx)):
      rad[i] = rad[i+1] - 0.5 * (temp[i] / mu[i] + temp[i+1] / mu[i+1]) * \
               (sc.Avogadro * sc.k * np.log(press[i]/press[i+1]) / g[i+1])
      g[i] = g[i+1] * rad[i+1]**2 / rad[i]**2
  # Calculate radius above p0:
  for i in np.arange(idx+1, n):
      rad[i] = rad[i-1] + 0.5 * (temp[i] / mu[i] + temp[i-1] / mu[i-1]) * \
               (sc.Avogadro * sc.k * np.log(press[i-1]/press[i]) / g[i-1])
      g[i] = g[i-1] * rad[i-1]**2 / rad[i]**2

  # Reverse the order of calculated radii to write them in the right order
  # in pre-atm file:
  rad = rad[::-1]

  return rad


def makeAbun(solar_abun, abun_file, solar_times=1, COswap=False):
    """
    This function makes the abundaces file to be used by BART.
    The function uses Asplund et al (2009) elemental abundances file
    http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A, (abudances.txt),
    and multiplies the abundances by the number user wants, or swaps the
    C/O ratio.

    Parameters
    ----------
    solar_abun: String
       Input Solar abundances filename.
    abun_file: String
       Output filename to store the modified elemental abundances.

    Optional parameters
    -------------------
    solar_times: Integer
       Multiplication factor for metallic elemental abundances (everything
       except H and He).
    COswap: Boolean
       If True, swap the abundances of C and O.

    Returns
    -------
    None

    Developers
    ----------
    Jasmina Blecic     UCF  jasmina@physics.ucf.edu
    Patricio Cubillos  UCF  pcubillos@fulbrightmail.org

    Revisions
    ---------
    2014-07-12  Jasmina   Written by.
    2014-08-15  Patricio  Rewrote data handling. Updated data strings.
    2014-09-24  Jasmina   Updated documentation.
    2015-03-06  Patricio  Updated code to read the solar abundances.
    """

    # Read the solar abundances file:
    index, symbol, dex, name, mass = read_eabun(solar_abun)
    # Count the number of elements:
    nelements = len(symbol)

    # Scale the metals aundances:
    imetals = np.where((symbol != "H") & (symbol != "He"))
    dex[imetals] += np.log10(solar_times)

    # Swap C and O abundances if requested:
    if COswap:
      Cdex = dex[np.where(symbol == "C")]
      dex[np.where(symbol == "C")] = dex[np.where(symbol == "O")]
      dex[np.where(symbol == "O")] = Cdex

    # Save data to file
    f = open(abun_file, "w")
    # Write header
    f.write("# Elemental abundances:\n"
            "# Columns: ordinal, symbol, dex abundances, name, molar mass.\n")
    # Write data
    for i in np.arange(nelements):
      f.write("{:3d}  {:2s}  {:5.2f}  {:10s}  {:12.8f}\n".format(
              index[i], symbol[i], dex[i], name[i], mass[i]))
    f.close()


# calculates species stoichiometric values
def stoich(species):
    '''
    Species counting function. Counts the number of each element in a chemical
    species. Takes in a string of a chemical species (i.e., "H2O") and returns
    an array containing every element with corresponding counts found in that
    species.

    Parameters
    ----------
    species : string
             Chemical species name. MUST include elemental species listed in
             the order they appear in the 'abundances.txt'.

    Returns
    -------
    stoich_info : 2D array
               Array containing two columns of equal length: the first
               column is a list of atomic symbols of the elements present
               in a species and the second column counts of each of the
               elements found in the species.

    Notes
    ------
    "Weight" in the code is the count of each element occurrence in the
    species, and the sum of all weights for that element is the stoichiometric
    coefficient (i.e., ClSSCl that appears in JANAF tables has weight 1 for
    first occurrence of Cl, weight 1 for first occurrence of S, and the final
    stoichiometric values of Cl is 2, and for S is 2).

    Capitals imply the beginning of an atomic symbol (i.e., He, Be,  Mg, etc)
    and the digit indicates the count of the element (weight) preceding it in
    the species (i.e, H2 has 2 H's and Li4 has 4 Li's).

    Revisions
    ---------
    2014-06-01  Oliver/Jasmina Written by.
    2014-09-16  Jasmina        Modified for BART project to return only data
                               of interest.

    2018-10-22  Michael        Modified to work for ions
    '''

    # Allocate string length and array of booleans to indicate if characters
    #          are capitals, digits, or +/-
    chars   = len(species)
    iscaps  = np.empty(chars, dtype=np.bool)
    isdigit = np.empty(chars, dtype=np.bool)
    isplus  = np.empty(chars, dtype=np.bool)
    isminus = np.empty(chars, dtype=np.bool)

    # Check each character in string to fill in boolean arrays for capitals,
    #       digits, and +/-
    for i in np.arange(len(species)):
        iscaps[i]  = (re.findall('[A-Z]', species[i]) != [])
        isdigit[i] = species[i].isdigit()
        isplus[i]  = species[i] == '+'
        isminus[i] = species[i] == '-'

    # Indicator for ending each count and blank stoich_info array
    endele = True
    stoich_info = []

    # Loop over all characters in species string
    for i in np.arange(len(species)):
        # Start tracking new element
        if endele == True:
            ele = ''
            weight = 0
            endele = False

        # Check if character is a letter, if so add to element name
        if (isdigit[i] == False and isplus[i] == False and isminus[i] == False):
            ele += species[i]

        # Check if character is a digit, if so, make this the element's weight
        if isdigit[i] == True:
            weight = np.int(species[i])

        # Check if charcater is a + or -, if so, make it an electron
        if isplus[i] == True:
            ele    += 'e'
            weight  = -1
            endele  = True

        if isminus[i] == True:
            ele    += 'e'
            weight  = 1
            endele  = True

        # Check if element name ends (next capital is reached)
        #       and if no weight (count) is found, set it to 1
        if (isdigit[i] == False and                \
           (iscaps[i+1:i+2] == True or i == chars-1)):
            weight = 1

        # If next element is found or if end of species name is reached
        #    (end of string), stop tracking
        if (iscaps [i+1:i+2] == True or isplus[i+1:i+2] == True or 
            isminus[i+1:i+2] == True or i == chars-1):
            endele = True

        # End of element has been reached, so output weights of element
        #     into elemental array
        if endele == True:
            stoich_info.append([ele,weight])

    # Return full array of elements and stoichiometry
    return stoich_info


# calculates mean molecular mass
def mean_molar_mass(abun_file, atmfile=None, spec=None, pressure=None,
                    temp=None, abundances=None):
    """
    This function calculates mean molar mass at each layer in the atmosphere.
    For input elements it trims the data from the abundances file, and
    extracts elemental molar mass. Then, it reads the final TEA
    output atmospheric file to get all the data, trims the names of the
    output species and makes a stoichiometric array to store the values of
    all output species. It calls the stoich() function to get each species
    stoichiometric values, multiplies elemental weights with each element
    number in a species and sum them for all output species at each layer
    in the atmosphere. It stores the values in the mu array for every layer.

    Parameters
    ----------
    abun_file: string
       Name of the file carrying abundance information
    atmfile: String
       Name of TEA atmospheric ASCII file.

    Returns
    -------
    mu: 1D array of floats
       Array containing mean molecular mass for each layer in the atmosphere.

    Revisions
    ---------
    2014-09-29  Jasmina   Written by.
    2015-03-05  Patricio  Simplified a few calculations.
    """

    # Read the elemental abundances file:
    index, element, dex, name, weights = read_eabun(abun_file)

    # Read the atmospheric file:
    if atmfile is not None:
      spec, pressure, temp, abundances = readatm(atmfile)

    # Number of molecules:
    nspec   = len(spec)
    # Number of layers in the atmosphere:
    nLayers = len(abundances)

    spec_weight = np.zeros(nspec)
    # Get the mass of each species:
    for i in np.arange(nspec):
      # Remove the JANAF extension from species name and get the
      #  stoichiometric data:
      spec_stoich = stoich(spec[i].partition('_')[0])
      # Add the mass from each element in this species:
      for j in np.arange(len(spec_stoich)):
        # Find the element:
        elem_idx = np.where(element == spec_stoich[j][0])
        # Add the weighted mass:
        spec_weight[i] += weights[elem_idx][0] * float(spec_stoich[j][1])

    # Allocate array (for each layer) of mean molar mass:
    mu = np.zeros(nLayers)

    # Sum of all species weight in the layer:
    for i in np.arange(nLayers):
        mu[i] = sum(spec_weight * abundances[i])

    return mu


def makeRadius(out_spec, atmfile, abun_file, tepfile, p0):
  """
  Add radius array into the final TEA output atmospheric file.
  It opens a new file to write, adds headers, reads the final TEA output
  atmospheric file to take the species, pressure, temperature, and abundances,
  calls the mean_molar_mass() function to calculate mu, calls radpress() to
  calculate radius, and then write all the data into a new file. Radius array
  is added as the first column in the file, the rest of the TEA format is
  preserved.

  Parameters
  ----------
  out_spec: String
      String containing all output molecular species.
  atmfile: String
      Name of TEA atmospheric ASCII file.
  abun_file: String
      Name of the abundances file.
      (default: 'abundances.txt', Asplund et al 2009)
  tepfile: String
      Name of the tepfile.
  p0: Float
      Reference pressure level (corresponding to Rplanet from the tepfile).

  Revisions
  ---------
  2014-09-20 Jasmina   Written by.
  2015-05-03 Jasmina   Corrected atm header.
  """

  # Read the final atmospheric file
  molecules, pressure, temperature, abundances = readatm(atmfile)

  # Calculate the mean molecular mass of each layer:
  mu = mean_molar_mass(abun_file, atmfile)

  # Open atmfile to overwrite:
  fout = open(atmfile, 'w')

  # Write a header file
  header = (
  "# This is a final TEA output file with calculated abundances (mixing fractions) for all listed species.\n"
  "# Units: pressure (bar), temperature (K), abundance (unitless).")
  fout.write(header + '\n\n')

  # Retrieve planet name and surface radius
  g, Rp = get_g(tepfile)

  # Write species names
  fout.write('#SPECIES\n' + out_spec + '\n\n')

  # Write TEA data
  fout.write('#TEADATA\n')

  # Calculate the radius of each layer:
  rad = radpress(pressure, temperature, mu, p0, Rp, g)

  # Number of layers in the atmosphere
  nLayers = len(abundances)

  # Number of molecules
  nspec = len(molecules)

  # Make a list of labels
  label = ['#Radius'.ljust(11)] + ['Pressure'.ljust(11)] + ['Temp'.ljust(8)]
  for i in np.arange(nspec):
       label = label + [molecules[i].ljust(14)]
  label = ''.join(label)

  # Write new label
  fout.write(label + '\n')

  # Write atm file for each run
  for i in np.arange(nLayers):
      # Radius, pressure, and temp for the current line
      radi = str('%10.3f'%rad[i])
      presi = str('%10.4e'%pressure[i])
      tempi = str('%7.2f'%temperature[i])

      # Insert radii array
      fout.write(radi.ljust(10) + ' ')

      # Insert results from the current line (T-P) to atm file
      fout.write(presi.ljust(10) + ' ')
      fout.write(tempi.ljust(7) + ' ')

      # Write current abundances
      for j in np.arange(nspec):
          fout.write('%1.4e'%abundances[i][j] + ' ')
      fout.write('\n')

  # Close atm file
  fout.close()


def make_preatm(tepfile, press_file, abun_file, in_elem, out_spec,
                pre_atm, Temp):
  """
  This code produces a pre-atm file in the format that TEA can read it.
  It reads the pressure file and elemental dex abundance data, trims to
  the selected elements.  It converts dex abundances to number density
  and divide them by hydrogen number density to get fractional abundances.
  It writes the pressure, temperature, and elemental-abundances data
  into a pre-atm file.

  Parameters
  ----------
  tepfile: String
     Name of the tepfile.
  press_file: String
     Name of the pressure file.
  abun_file: String
     Name of the abundances file.
  in_elem: String
     String containing input elemental species.
  out_spec: String
     String containing output molecular species.
  pre_atm: String
     Pre-atmospheric filename.
  Temp: 1D float array
     Array containing temperatures for each layer in the atmosphere (in K).

  Revisions
  ---------
  2014-07-12  Jasmina   Written by.
  2014-08-15  Patricio  Removed PT re-calculation, use input T array.
  2014-09-15  Jasmina   Added call to readAbun() function.
  2014-11-06  Jasmina   Adjusted to work properly with TEA.
  2015-03-05  Patricio  Adapted to use read_eabun function, reworked some bits.
  """

  # Read pressure data
  pres = pt.read_press_file(press_file)

  # Number of layers in the atmosphere
  n_layers = len(pres)

  # Get the elemental abundace data:
  index, symbol, dex, name, mass = read_eabun(abun_file)

  # Take only the elements we need:
  in_arg = np.in1d(symbol, in_elem.split())
  in_symbol = symbol[in_arg]
  in_dex    = dex   [in_arg]

  # Hydrogen number density:
  H_num = 10**12

  # Get fractional number density relative to hydrogen:
  out_abn  = (10.0**in_dex / H_num).tolist()

  # Write pre-atm file
  f = open(pre_atm, 'w')

  # Pre-atm header with basic instructions
  header = ("# This is a TEA pre-atmosphere input file.\n"
   "# TEA accepts a file in this format to produce species abundances as\n"
   "# a function of pressure and temperature.\n"
   "# Output species must be added in the line immediately following the \n"
   "# SPECIES marker and must be named to match JANAF converted names.\n"
   "# Units: pressure (bar), temperature (K), abundance (unitless).\n\n")

  # Write header
  f.write(header)

  # Write species names
  f.write('#SPECIES\n' + out_spec + '\n\n')

  # Write header of TEA data
  f.write('#TEADATA\n')
  # Column's header:
  f.write("#Pressure   Temp          " +
        "".join(["{:<18s}".format(elem) for elem in in_symbol]) + "\n")

  # Write data for each layer:
  for i in np.arange(n_layers):
    # Pressure, and Temperature:
    f.write("{:10.4e} {:8.2f}  ".format(pres[i], Temp[i]))
    # Elemental abundance list:
    f.write("  ".join(["{:16.10e}".format(abun) for abun in out_abn]) + "\n")
  f.close()


def uniform(atmfile, press_file, abun_file, tepfile, species, abundances, temp, p0):
  """
  Generate an atmospheric file with uniform abundances.

  Parameters:
  -----------
  atmfile: String
     Name of output atmospheric file.
  press_file: String
     Input pressure-array filename.
  abun_file: String
     Input elemental-abundances filename.
  tepfile: String
     Transiting extrasolar planet filename.
  species: String
     String with list of atmospheric species (blank space separated).
  abundances: String
     String with list of species mole-mixing-ratio (blank space separated).
  temp: 1D float ndarray
     Temperature profile.
  p0: Float
     Reference pressure level (corresponding to Rp from the tepfile).

  Modification History:
  ---------------------
  2015-03-04  patricio  Initial implementation.
  2018-10-22  mhimes    Altered header writing to have wider fields to avoid 
                        bugs with longer ion names
  """

  # Read pressure array:
  press = pt.read_press_file(press_file)
  # Put species names into an array:
  spec = np.asarray(species.split())
  # Put abundance values into an array:
  abun = np.asarray(abundances, np.double)

  # Write to file:
  f = open(atmfile, 'w')

  # Write species names:
  f.write('#SPECIES\n' + species + '\n\n')

  # Write header of TEA data:
  f.write('#TEADATA\n')
  # Column's header:
  f.write("#Pressure   Temp     " +
        "".join(["{:<18s}".format(mol) for mol in spec]) + "\n")

  # For each layer write TEA data
  for i in np.arange(len(press)):
      # Pressure, and Temperature
      f.write("{:10.4e} {:8.2f}  ".format(press[i], temp[i]))
      # Elemental abundance list
      f.write("  ".join(["{:16.10e}".format(ab) for ab in abun]) + "\n")
  f.close()

  # Calculate the radius of each layer and put it into the atmfile:
  makeRadius(species, atmfile, abun_file, tepfile, p0)
  # Reformat to use in transit:
  reformat(atmfile)


# reads final TEA atmospheric file
def readatm(atmfile):
    """
    Reads TEA atmospheric file.

    Parameters
    ----------
    atmfile: String
               Name of TEA atmospheric ASCII file.

    Returns
    -------
    molecules: 1D string array
               Array of output molecular species.
    pressure: 1D float array
               Array of pressures listed in the atmospheric file (in bar).
    temp: 1D float array
               Array of pressures listed in the atmospheric file (in K)
    abundances: 2D float array
               Array containing abundances (mixing fractions) for
               each species at each atmospheric layer (unitelss).

    Notes
    -----
    Atmospheric data starts two lines after the header line: #TEADATA

    Revisions
    ---------
    2013-11-17  Jasmina   Written by.
    2014-09-05  Patricio  Modified from getpressure.
    2014-09-12  Jasmina   Added new markers, excluded radii, returned
                          temperature, added documentation.
    2014-10-02  Jasmina   Added option to read atmfile with or without
                          radius array in it.
    """

    # Open the atmospheric file and read
    f = open(atmfile, 'r')
    lines = np.asarray(f.readlines())
    f.close()

    # Get the molecules
    imol = np.where(lines == "#SPECIES\n")[0][0] + 1
    molecules = lines[imol].split()

    # Find the line where the layers info begins
    start = np.where(lines == "#TEADATA\n")[0][0] + 2

    # Number of columns
    ncol = len(lines[start].split())

    # Number of layers
    ndata = len(lines) - start

    # Read atmfile without radius array in it or with it
    pressure = np.zeros(ndata, np.double)
    temp     = np.zeros(ndata, np.double)
    if ncol == len(molecules) + 2:
        # Number of abundances (elements per line except Press and T)
        nabun = len(lines[start].split()) - 2
        abundances = np.zeros((ndata, nabun), np.double)

        # Extract pressure, temperature, and abundance data
        for i in np.arange(ndata):
            line = lines[start+i].strip().split()
            pressure[i]   = line[0]
            temp[i]       = line[1]
            abundances[i] = line[2:]
    else:
        # Number of abundances (elements per line except Radius, Press and T)
        nabun = len(lines[start].split()) - 3
        abundances = np.zeros((ndata, nabun), np.double)

        # Extract pressure, temperature, and abundance data
        for i in np.arange(ndata):
            line = lines[start+i].strip().split()
            pressure[i]   = line[1]
            temp[i]       = line[2]
            abundances[i] = line[3:]

    return molecules, pressure, temp, abundances


# reformats final TEA atmospheric file for Transit
def reformat(atmfile):
    """
    Re-format the TEA output file to work for transit.

    Parameters
    ----------
    atmfile: String
       Name of TEA atmospheric ASCII file.

    Revisions
    ---------
    2014-08-25  Patricio  Written by.
    2014-09-20  Jasmina   Fixed issue with the line break.
    2018-10-22  Michael   Added support for ions
    """

    # Open the atmospheric file to read and write
    f = open(atmfile, 'r')
    lines = np.asarray(f.readlines())
    f.close()

    # Find the Molecule names and remove their suffix
    imol = np.where(lines == "#SPECIES\n")[0][0] + 1
    molecules = lines[imol].split()
    for m in np.arange(len(molecules)):
        molecules[m] = molecules[m].replace('_ion_p', '+')
        molecules[m] = molecules[m].replace('_ion_n', '-')
        molecules[m] = molecules[m].partition('_')[0]
    lines[imol] = " ".join(molecules) + "\n"

    # Repeat for the column headers
    start = np.where(lines == "#TEADATA\n")[0][0] + 1
    molecules = lines[start].split()
    for m in np.arange(len(molecules) - 1):
        molecules[m] = "{:10s}".format(molecules[m].replace('_ion_p', '+').replace('_ion_n', '-').partition('_')[0])
    molecules[len(molecules) - 1] = molecules[len(molecules) - 1].replace('_ion_p', '+').replace('_ion_n', '-').partition('_')[0]
    lines[start] = " ".join(molecules) + "\n"

    lines = list(lines)
    # Reverse order to have layers from bottom to top
    datalines = lines[start+1:]
    datalines.reverse()
    lines = lines[:start+1] + datalines

    # Add values' units
    # Abundance by number density
    lines.insert(imol-2, "q number\n")
    # Pressure in bars:
    lines.insert(imol-2, "up 1e6\n")
    # Radius in kilometers:
    lines.insert(imol-2, "\n#Values units:\nur 1e5\n")

    # Save file with the updated lines
    f = open(atmfile, 'w')
    f.writelines(lines)
    f.close()
