#! /usr/bin/env python

# ****************************** START LICENSE ******************************
# Thermal Equilibrium Abundances (TEA), a code to calculate gaseous molecular
# abundances under thermochemical equilibrium conditions.
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

import numpy as np
import scipy.constants as sc
from scipy.interpolate import interp1d
import os
import shutil

import PT as pt
import reader as rd
import re

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
    Possible user errors in making pre-atm file that conflict with TEA:
     - The user should have unique T-P values for each line in the
     PT-profile file before running \texttt{makeatm.py}, otherwise
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
    2014-06-01 Oliver/Jasmina Made initial makeatm() function for TEA.
    2014-07-14 Jasmina        Added get_g, radpress, and read_press_file
                              functions, and rewrote and renames 
                              make_preatm() for BART project.
    2014-08-15 Patricio       Removed read_press_file, use function from PT.py.
                              Modified read_press_file to readatm(). 
                              Added reformat()
    2014-09-15 Jasmina        Added makeRadius(), readAbun(), 
                              mean_molar_mass(), stoich(). Modified readatm(),
                              make_preatm(), reformat().
    2014-09-27 Jasmina        Modified radpress() to calculate radius from the
                              referenced surface pressure of 1bar. 
                              Cleaned documentation.
    2014-10-02  Jasmina       Added option to read atmfile with or without radius 
                              array in it.
    2014-11-06  Jasmina       Minor adjustments for proper TEA execution.
"""

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


# calculates radii for each pressure in the atmosphere
def radpress(tepfile, temp, mu, pres):
    '''
    Given a pressure in bar, temperature in K, mean molecular
    mass in g/mol, the function calculates radii in km for each 
    layer in the atmosphere. It declares constants, calls get_g()
    function, allocates array for radii, and calculates radii for
    each pressure. The input pressure and temperature arrays MUST be
    in ascending order before interpolation is done, and then in
    in descending order because reference radius Rp is placed on
    the bottom of the atmosphere where pressure is largest. 
    All radii are calculated based on the reference level (p = 1bar). 
    The final radii array is calculated based on the referenced surface
    radius, above and below the surface radius. Above values are positive
    below values are negative relative to the value of the surface radius.
    Radii array is given in km, and converted back to ascending order for
    pre-atm file. 

    Parameters
    ----------
    tepfile: String
          Name of the tep ASCII file.
    temp: 1D array of floats
          Array containing temperatures for each layer in the atmosphere (in K).
    mu: 1D array of floats
          Array containing mean molecular mass for each layer from in the atmosphere
          (in g/mol).
    pres: 1D array of floats
          Array containing pressures for each layer in the atmosphere (in bar).

    Returns
    -------
    rad: 1D array of floats
        Array containing radii for each pressure layer in the atmosphere (in km).

    Revisions
    ---------
    2014-07-04 Jasmina   Written by.
    2014-09-29 Jasmina   Added radii calculation based on the referenced pressure 
                         level of 1bar and surface radius defined in the tepfile.
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

    # Referenced surface pressure
    p0 = 1.0  # bar

    # Interpolate temp and mu in lin-log space (1bar)
    interPT = interp1d(np.log10(pres), temp)
    intermu = interp1d(np.log10(pres),   mu)

    # Get temp and mu at surface pressure (1bar)
    try:
        temp_1bar = interPT(np.log10(p0))
        mu_1bar   = intermu(np.log10(p0))
    except IOError:
        print("\nReferenced surface pressure of 1bar is not in the range of pressures: " + \
               "[" + str(min(pres), max(pres)) + "]\n")

    # Return back to desending order for radius calculation
    pres = pres[::-1]
    temp = temp[::-1]
    mu   =   mu[::-1]

    # Find the index of the closest pressure point to surface (1bar)
    idx = (np.abs(pres - p0)).argmin()

    # Andrews:"Introduction to Atmospheric Physics" page 26
    # rad2 - rad1 = (R * T1)/g * ln(p1/p2)
    # R = R*/mu, R* - universal gas constant, mu - mean molecular mass
    # R* = Avogadro * Boltzmann
    # rad2 = rad1 + (Avogadro * Boltzmann * T1) / (g * mu1) * ln(p1/p2)
    # to convert from g to kg, mu[i] need to be multiplied with 1000

    # ====== If pres = 1 bar is not in the list of pressure points ====== #
    if float(pres[idx]) != p0:
        # If the point is above the surface (1bar)
        if pres[idx] > p0:
            rad[idx] = Rp + 1./2. * (temp[idx] / mu[idx] + temp_1bar / mu_1bar) * \
                       (Avogadro * Boltzmann * np.log(p0/pres[idx]) / g)
        # If the point is below the surface (1 bar)
        else:
            rad[idx] = Rp - 1./2. * (temp[idx] / mu[idx] + temp_1bar / mu_1bar) * \
                       (Avogadro * Boltzmann * np.log(pres[idx]/p0) / g)
    # ====== If pres = 1 bar is in list of pressure points ====== #
    else:
        rad[idx] = Rp
 
    # Calculate radius below the surface (1 bar)
    for i in reversed(np.arange(idx)):
        rad[i] = rad[i+1] - 1./2. * (temp[i] / mu[i] + temp[i+1] / mu[i+1]) * \
                 (Avogadro * Boltzmann * np.log(pres[i]/pres[i+1]) / g)
    # Calculate radius above the surface (1bar)  
    for i in np.arange(idx+1, n):
        rad[i] = rad[i-1] + 1./2. * (temp[i] / mu[i] + temp[i-1] / mu[i-1]) * \
                 (Avogadro * Boltzmann * np.log(pres[i-1]/pres[i]) / g)  

    # Reverse the order of calculated radii to write them in the right order
    #         in pre-atm file
    rad = rad[::-1]

    # Subtract radius array from the referenced surface radius
    rad = rad #- Rp

    return rad


# reads abundances for elements of interest
def readAbun(in_elem, abun_file):
    """
    Read the basic abundances file,'abundances.txt', that carries
    Asplund et al (2009) solar photosphere elemental abundances
    information.

    Parameters:
    in_elem: string
             String with all input elemental species.
    abun_file: string
             Name of the abundances file (default: abundances.txt)

    Returns:
    abun_trim: list of strings
             List containing information from abundances.txt file only
             for elements of interest

    Notes: 
             in_elem MUST be in the order thay appear in 'abundances.txt'

    Revisions
    ---------
    2014-09-10  Jasmina   Written by.
    """

    # Read abundance data and convert to array (skip rows with comments):
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

    return abun_trim


# calculates species stoichiometric values
def stoich(specie):
    '''
    Species counting function. Counts the number of each element in a chemical
    species. Takes in a string of a chemical species (i.e., "H2O") and returns
    an array containing every element with corresponding counts found in that
    species. 

    Parameters
    ----------
    specie : string
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
    '''

    # Allocate string length and array of booleans to indicate if characters
    #          are capitals or digits
    chars   = len(specie)
    iscaps  = np.empty(chars, dtype=np.bool)
    isdigit = np.empty(chars, dtype=np.bool)
    
    # Check each character in string to fill in boolean arrays for capitals
    #       or digits; 
    for i in np.arange(len(specie)):
        iscaps[i] = (re.findall('[A-Z]', specie[i]) != [])
        isdigit[i] = specie[i].isdigit()
    
    # Indicator for ending each count and blank stoich_info array
    endele = True
    stoich_info = [[]]
    
    # Loop over all characters in species string
    for i in np.arange(len(specie)):  
        # Start tracking new element
        if endele == True:
            ele = ''
            weight = 0
            endele = False
        
        # Check if character is a letter, if so add to element name
        if (isdigit[i] == False): 
            ele += specie[i]
        
        # Check if character is a digit, if so, make this the element's weight
        if isdigit[i] == True:
            weight = np.int(specie[i])
        
        # Check if element name ends (next capital is reached) 
        #       and if no weight (count) is found, set it to 1
        if (isdigit[i] == False and                \
           (iscaps[i+1:i+2] == True or i == chars-1)):
            weight = 1
        
        # If next element is found or if end of species name is reached 
        #    (end of string), stop tracking
        if (iscaps[i+1:i+2] == True or i == chars-1): 
            endele = True
        
        # End of element has been reached, so output weights of element 
        #     into elemental array
        if endele == True:

            # Create array containing only the elements used in this run 
            if stoich_info == [[]]:
                stoich_info = np.append(stoich_info, [[ele, weight]], axis=1)
            else:
                stoich_info = np.append(stoich_info, [[ele, weight]], axis=0)
    
    # Return full array of elements and stoichiometry
    return stoich_info


# calculates mean molecular mass
def mean_molar_mass(in_elem, abun_file, atmfile):
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
    in_elem: string
             String containing input elements. Input elements MUST be in the
             order they appear in the 'abundances.txt' file.
    abun_file: string
             Name of the file carrying abundance information
             (default: 'abundances.txt', Asplund et al. 2009)
    atmfile: String
             Name of TEA atmospheric ASCII file.

    Returns
    -------
    mu: 1D array of floats
        Array containing mean molecular mass for each layer in the atmosphere.

    Revisions
    ---------
    2014-09-29 Jasmina   Written by.
    """

    # Take the abundance data of interest
    abun_trim = readAbun(in_elem, abun_file)

    # Take only weights from abundance data and convert them to floats
    weights = map(float, abun_trim[:,4])

    # Read the final atmospheric file
    out_spec, pressure, temp, abundances = readatm(atmfile)

    # Number of molecules
    nspec = len(out_spec)

    # Take the JANAF extension from out_spec names
    for m in np.arange(len(out_spec)):
        out_spec[m] = out_spec[m].partition('_')[0]

    # Take elements from abun_trim to get the correct elements order
    elements = abun_trim[:,1]

    # Initiate a list of stoichiometric values
    stoich_val = []
    
    # Loop over all species to make 2D array of stoichiometry
    for i in np.arange(nspec):
        # Call stoich function and get each species full info 
        #      each element name and its stoichiometric value
        spec_info = stoich(out_spec[i])
        
        # Initiate list of each species stoichiometry
        spec_stoich = np.zeros(len(elements), int).tolist()

        # Loop over all species information
        for j in np.arange(len(spec_info)):

            # Take each element info
            elem_info = spec_info[j]
 
            # Take the name of the element
            atom = elem_info[0]

            # Loop over all elements
            for k in np.arange(len(elements)):

               # If name of the atom equal to the name of the input element
               if  atom == elements[k]:
                   
                   # Fill out species stoichiometry
                   spec_stoich[k] = int(elem_info[1])

        # To the final stoichiometric values append each species stoichiometry
        stoich_val.append(spec_stoich)

        # Multiple elemental molar mass with stoichiometric value
        elem_weight = np.asarray(weights) * np.asarray(stoich_val)  

    # Allocate space for the sum of elemental weights
    spec_weight = np.zeros(nspec)

    # Sum all elemental weights to get the molar mass of the species
    for i in np.arange(nspec):
            spec_weight[i] = sum(elem_weight[i])

    # Number of layers in the atmosphere
    nLayers = len(abundances)

    # Allocate space for mu, final mean molar mass for each T-P
    mu = np.zeros(nLayers)

    # Sum of all species weight fraction for one T-P
    for i in np.arange(nLayers):
        mu[i] = sum(spec_weight * abundances[i])

    return mu


# adds radius array to the final TEA atmospheric file
def makeRadius(in_elem, out_spec, atmfile, abun_file, tepfile):
    """
    This function adds radius array into the final TEA output atmospheric file.
    It opens a new file to write, adds headers, reads the final TEA output
    atmospheric file to take the species, pressure, temperature, and abundances, 
    calls the mean_molar_mass() function to calculate mu, calls radpress() to
    calculate radius, and then write all the data into a new file. Radius array
    is added as the first column in the file, the rest of the TEA format is
    preserved. 

    Parameters
    ----------
    in_elem: string
             String containing input elements.
    out_spec: String
             String containing all output molecular species.
    atmfile: String
             Name of TEA atmospheric ASCII file.
    abun_file: String
             Name of the abundances file.
             (default: 'abundances.txt', Asplund et al 2009)
    tepfile: String
             Name of the tepfile.

    Return
    ------
    None

    Revisions
    ---------
    2014-09-20 Jasmina   Written by.
    2014-09-26 Jasmina   Added correct atm header.
    """

    # Make a copy of the original TEA atmfile
    #shutil.copy2(atmfile, atmfile[:-4] + '_genuine.tea')

    # Open atmfile to read and write
    fout = open(atmfile, 'r+')

    # Write a header file
    header = ("# This is a final TEA output file with calculated abundances (mixing fractions) for all listed species. \n"
              "# Units: radius (km), pressure (bar), temperature (K), abundance (unitless).")

    # Write header
    fout.write(header + '\n\n')
    
    # Retrieve planet name and surface radius
    head, tepname = os.path.split(tepfile)
    g, Rp = get_g(tepfile)

    # Write planet name and surface radius
    fout.write('#' + tepname[:-4] +'\n\n')
    #fout.write('Rp ' + str(Rp) + '\n\n')

    # Write species names
    fout.write('#SPECIES\n' + out_spec + '\n\n')

    # Write TEA data
    fout.write('#TEADATA\n')

    # Read the final atmospheric file
    molecules, pressure, temperature, abundances = readatm(atmfile)

    # Call mean_molar_mass() function to calculate mu for each T-P
    mu = mean_molar_mass(in_elem, abun_file, atmfile)

    # Call radpress() to calculate radii
    rad = radpress(tepfile, temperature, mu, pressure)

    # Number of layers in the atmosphere
    nLayers = len(abundances)

    # Number of molecules
    nspec = len(molecules)

    # Make a list of labels
    label = ['#Radius'.ljust(11)] + ['Pressure'.ljust(11)] + ['Temp'.ljust(8)]
    #label = ['#Radius-Rp'.ljust(11)] + ['Pressure'.ljust(11)] + ['Temp'.ljust(8)]
    for i in np.arange(nspec):
         label = label + [molecules[i].ljust(11)]
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


# Make pre-atmospheric file for TEA
def make_preatm(tepfile, press_file, abun_file, in_elem, out_spec,
                pre_atm, Temp):
    '''
    This code produces a pre-atm file in the format that TEA can read it. 
    It reads the pressure file and trimmed elemental dex abundance data 
    by calling readAbun() function. It converts dex abundances of all 
    elements of interest to number density and divides them by the sum
    of all number densities in the mixture to get fractional abundances.
    It calls radpress() to calculate radii, and writes the data 
    (pressure, temperature, elemental abundances) into a pre-atm file.

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

    Returns
    -------
    None

    Revisions
    ---------
    2014-07-12  Jasmina   Written by.
    2014-08-15  Patricio  Removed PT re-calculation, use input T array.
    2014-09-15  Jasmina   Added call to readAbun() function.
    2014-11-06  Jasmina   Adjusted to work properly with TEA.
    '''

    # Read pressure data
    pres = pt.read_press_file(press_file)

    # Number of layers in the atmosphere
    n_layers = len(pres)

    # Take the abundace data of interest
    abun_trim = readAbun(in_elem, abun_file)

    # Take data and create list
    out_elem = abun_trim[:,1].tolist()

    # Convert strings to floats
    out_dex  = map(float, abun_trim[:,2])

    # Convert dex exponents to number density
    out_num  = 10**np.array(out_dex)

    # Hydrogen number density
    H_num = 10**12

    # Get fractions of element number density to hydrogen number density
    out_abn  = (out_num / H_num).tolist()
 
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
          "".join(["{:<18s}".format(elem) for elem in out_elem]) + "\n")

    # For each layer write TEA data
    for i in np.arange(n_layers):
        # Pressure, and Temperature
        f.write("{:10.4e} {:8.2f}  ".format(pres[i], Temp[i]))
        # Elemental abundance list
        f.write("  ".join(["{:16.10e}".format(abun) for abun in out_abn]) + "\n")
    f.close()


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
    """

    # Open the atmospheric file to read and write
    f = open(atmfile, 'r')
    lines = np.asarray(f.readlines())
    f.close() 

    # Find the Molecule names and remove their suffix
    imol = np.where(lines == "#SPECIES\n")[0][0] + 1
    molecules = lines[imol].split()
    for m in np.arange(len(molecules)):
        molecules[m] = molecules[m].partition('_')[0]
    lines[imol] = " ".join(molecules) + "\n"

    # Repeat for the column headers
    start = np.where(lines == "#TEADATA\n")[0][0] + 1
    molecules = lines[start].split()
    for m in np.arange(len(molecules) - 1):
        molecules[m] = "{:10s}".format(molecules[m].partition('_')[0]) 
    molecules[len(molecules) - 1] = molecules[len(molecules) - 1].partition('_')[0]
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



