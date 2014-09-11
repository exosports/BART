import numpy as np
import os

def makeAbun(abun_basic, filename, solar_times=1, COswap=False): 
    '''
    Make the abundaces file to be used by BART.
    The function uses Asplund et al (2009) elemental abundances file
    (abudances.txt) and multiplies the abundances by the number user wants,
    or swaps the C/O ratio. It is called by BART_config.py. 
  
    Parameters:
    -----------
    abun_basic: String
       Name of the basic abundances file (default:abundances.txt)
    filename: String
       Name of the file.

    Optional parameters:
    --------------------
    solar_times: Integer
       Multiplication factor for metallic elemental abundances (everything
       except H and He).
    COswap: Boolean
       If True, swap the abundances of C and O.

    Returns:
    --------
    None

    Developers:
    -----------
    Jasmina Blecic     UCF  jasmina@physics.ucf.edu
    Patricio Cubillos  UCF  pcubillos@fulbrightmail.org

    Revisions:
    ----------
    2014-07-12  Jasmina   Initial version
    2014-08-15  Patricio  Updated docstring.  Reworked data handling.
    '''

    # Read the basic solar abundances file - abundances.txt
    f = open(abun_basic, 'r')
    lines = f.readlines()
    f.close()

    # Count the number of elements:
    nelements = len(lines)
    for line in lines:
      if line.startswith("#"):
        nelements -= 1
    
    # Allocate arrays to put info:
    index  = np.zeros(nelements, int)
    symbol = np.zeros(nelements, '|S2')
    dex    = np.zeros(nelements, np.double)
    name   = np.zeros(nelements, '|S20')
    mass   = np.zeros(nelements, np.double)

    # Store elements info into arrays:
    i = 0
    for line in lines:
      if not line.startswith("#"):
        index[i], symbol[i], dex[i], name[i], mass[i] = line.strip().split()
        i += 1

    # Scale metals aundance:
    imetals = np.where((symbol != "H") & (symbol != "He"))
    dex[imetals] += np.log10(solar_times)

    # Swap C and O abundances:
    if COswap:
      Cdex = dex[np.where(symbol == "C")]
      dex[np.where(symbol == "C")] = dex[np.where(symbol == "O")]
      dex[np.where(symbol == "O")] = Cdex

    # Save data to file:
    f = open(filename, "w")
    # Write header:
    f.write("# Elemental abundances:\n"
            "# Columns: ordinal, symbol, dex abundances, name, molar mass.\n")
    # Write data:
    for i in np.arange(nelements):
      f.write("{:3d}  {:2s}  {:5.2f}  {:10s}  {:12.8f}\n".format(
              index[i], symbol[i], dex[i], name[i], mass[i]))
    f.close()
