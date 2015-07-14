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
# Copyright (C) 2015 University of Central Florida.  All rights reserved.
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


# List of constants used in the BART project:

# Source: http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
Mjup = 1.8983e+27     # Jupiter mass in kg:
Rjup = 7.1492e+07     # Jupiter radius in m:

# Source: http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
Rsun = 6.96e+08       # Solar radius in m:

# Source: transit module, transit/include/constrants_tr.h
H  = 6.6260755e-27    # Planck's constant (erg * s)
LS = 2.99792458e10    # Light Speed (cm / s)
KB = 1.380658e-16     # Boltzmann constant (erg / K) 
