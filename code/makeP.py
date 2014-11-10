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

def makeP(n_layers, p_top, p_bottom, filename, log=True): 
    '''
    Function to make a logarithmic or linearly equally-spaced pressure
    array and store it in a file.
  
    Parameters:
    -----------
    n_layers: Integer
       Number of layers in the atmosphere.
    p_top: Float
       Pressure at the top of the atmosphere (in bars).
    p_bottom: Float
       Pressure at the bottom of the atmopshere (in bars).
    filename: String
       Name of the ASCII file to be produced.
    log: Boolean
       If True, logarithmic sampling of the pressure range is chosen,
       else linear.

    Returns:
    --------
    None

    Developers:
    -----------
    Jasmina Blecic      jasmina@physics.ucf.edu
    Patricio Cubillos   pcubillos@fulbrightmail.org

    Revisions
    ---------
    2014-06-22  Jasmina   Written by.
    2014-08-15  Patricio  Complemented doc-string and updated string
                          formatting.
    2014-09-24  Jasmina  Updated documentation.
    '''

    # Make logarithmic or linearly equispaced array
    if log:
      pres = np.logspace(np.log10(p_top), np.log10(p_bottom), n_layers)   
    else:
      pres = np.linspace(p_top, p_bottom, n_layers)   

    # Header line
    header = "Layer  P (bar)\n"

    # Open file to write
    f = open(str(filename), 'w')
    f.write(header)

    # Write layer index number and pressure
    for i in np.arange(n_layers):
      f.write("{:5d}  {:.4e}\n".format(i+1, pres[i]))
    f.close()

