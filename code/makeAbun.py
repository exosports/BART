# ****************************** START LICENSE *******************************
# Bayesian Atmospheric Radiative Transfer (BART), a code to infer
# properties of planetary atmospheres based on observed spectroscopic
# information.
# 
# This project was completed with the support of the NASA Planetary
# Atmospheres Program, grant NNX12AI69G, held by Principal Investigator
# Joseph Harrington. Principal developers included graduate students
# Patricio E. Cubillos and Jasmina Blecic, programmer Madison Stemm, and
# undergraduates M. Oliver Bowman and Andrew S. D. Foster.  The included
# 'transit' radiative transfer code is based on an earlier program of
# the same name written by Patricio Rojo (Univ. de Chile, Santiago) when
# he was a graduate student at Cornell University under Joseph
# Harrington.  Statistical advice came from Thomas J. Loredo and Nate
# B. Lust.
# 
# Copyright (C) 2014 University of Central Florida.  All rights reserved.
# 
# This is a test version only, and may not be redistributed to any third
# party.  Please refer such requests to us.  This program is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.
# 
# Our intent is to release this software under an open-source,
# reproducible-research license, once the code is mature and the first
# research paper describing the code has been accepted for publication
# in a peer-reviewed journal.  We are committed to development in the
# open, and have posted this code on github.com so that others can test
# it and give us feedback.  However, until its first publication and
# first stable release, we do not permit others to redistribute the code
# in either original or modified form, nor to publish work based in
# whole or in part on the output of this code.  By downloading, running,
# or modifying this code, you agree to these conditions.  We do
# encourage sharing any modifications with us and discussing them
# openly.
# 
# We welcome your feedback, but do not guarantee support.  Please send
# feedback or inquiries to:
# 
# Joseph Harrington <jh@physics.ucf.edu>
# Patricio Cubillos <pcubillos@fulbrightmail.org>
# Jasmina Blecic <jasmina@physics.ucf.edu>
# 
# or alternatively,
# 
# Joseph Harrington, Patricio Cubillos, and Jasmina Blecic
# UCF PSB 441
# 4111 Libra Drive
# Orlando, FL 32816-2385
# USA
# 
# Thank you for testing BART!
# ******************************* END LICENSE *******************************

import numpy as np
import os

def makeAbun(abun_basic, filename, solar_times=1, COswap=False): 
    '''
    This function makes the abundaces file to be used by BART.
    The function uses Asplund et al (2009) elemental abundances file
    http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A, (abudances.txt), 
    and multiplies the abundances by the number user wants, or swaps the
    C/O ratio. 
  
    Parameters
    ----------
    abun_basic: String
       Name of the basic abundances file (default:abundances.txt)
    filename: String
       Name of the ASCII file.

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
    '''

    # Read the basic solar elemental abundances file - abundances.txt
    f = open(abun_basic, 'r')
    lines = f.readlines()
    f.close()

    # Count the number of elements
    nelements = len(lines)
    for line in lines:
      if line.startswith("#"):
        nelements -= 1
    
    # Allocate arrays to put info
    index  = np.zeros(nelements, int)
    symbol = np.zeros(nelements, '|S2')
    dex    = np.zeros(nelements, np.double)
    name   = np.zeros(nelements, '|S20')
    mass   = np.zeros(nelements, np.double)

    # Store elements info into arrays
    i = 0
    for line in lines:
      if not line.startswith("#"):
        index[i], symbol[i], dex[i], name[i], mass[i] = line.strip().split()
        i += 1

    # Scale metals aundance
    imetals = np.where((symbol != "H") & (symbol != "He"))
    dex[imetals] += np.log10(solar_times)

    # Swap C and O abundances
    if COswap:
      Cdex = dex[np.where(symbol == "C")]
      dex[np.where(symbol == "C")] = dex[np.where(symbol == "O")]
      dex[np.where(symbol == "O")] = Cdex

    # Save data to file
    f = open(filename, "w")
    # Write header
    f.write("# Elemental abundances:\n"
            "# Columns: ordinal, symbol, dex abundances, name, molar mass.\n")
    # Write data
    for i in np.arange(nelements):
      f.write("{:3d}  {:2s}  {:5.2f}  {:10s}  {:12.8f}\n".format(
              index[i], symbol[i], dex[i], name[i], mass[i]))
    f.close()
