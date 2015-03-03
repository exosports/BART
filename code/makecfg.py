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

import os, sys
import argparse, ConfigParser
import numpy as np

filedir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(filedir + "/../modules/MCcubed/src/")
import mcutils as mu


def makecfg(cfg, directory, atmfile):
  """
  Make a transit configuration file.

  Parameters:
  -----------
  cfg: String
     BART configuration file
  directory:  String
     Directory where to store the transit configuration file.
  atmfile: String
     Name of the atmospheric filename.

  Modification History:
  ---------------------
  2014-09-09  patricio  Initial implementation.
  """

  # Known (single-value) transit arguments:
  known_args = ["radlow",  "radhigh",  "raddelt",  "radfct",
                "wllow",   "wlhigh",   "wlfct",
                "wnlow",   "wnhigh",   "wndelt",  "wnfct", "wnosamp",
                "tlow",    "thigh",    "tempdelt",
                "allowq",  "nwidth",
                "opacityfile",
                #"molfile",
                "toomuch", "tauiso", "outtau", "taulevel", "modlevel",
                "starrad", "transparent",
                "solution", "raygrid",
                "gsurf", "refpress", "refradius",
                "cloudrad", "cloudfct", "cloudext",
                "verb",
                "outtoomuch", "outsample", "output",
               ]
  # Known multiple-value transit arguments:
  known_tuple_args = ["linedb", "cia", "orbpars", "orbparsfct"]

  # Get all BART arguments in a dict:
  config = ConfigParser.SafeConfigParser()
  config.read([cfg])
  defaults = dict(config.items("MCMC"))

  # Key names of the arguments in the BART configuration file:
  cfg_args = defaults.keys()

  # transit configuration filename:
  # FINDME: Add file-not-found exception
  transitcfg = defaults["config"]
  tcfg = open(directory + transitcfg, "w")

  # Print the atmospheric name:
  tcfg.write("atm {:s}\n".format(atmfile))
  # Default the molfile path:
  tcfg.write("molfile {:s}\n".format(
      os.path.normpath(filedir + "/../modules/transit/inputs/molecules.dat")))

  # Print the known arguments:
  for key in known_tuple_args:
    if key in cfg_args:
      values = defaults[key].strip().split()
      for val in values:
        tcfg.write("{:s} {:s}\n".format(key, val))

  for key in known_args:
    if key in cfg_args:
      tcfg.write("{:s} {:s}\n".format(key, defaults[key]))

  tcfg.close()


def makeTEA(cfg, TEAdir):
  """
  Make a TEA configuration file.

  Parameters:
  -----------
  cfg: String
     BART configuration file
  TEAdir: String
     Default TEA directory.

  Modification History:
  ---------------------
  2014-12-08  patricio  Initial implementation.
  """
  # Open New ConfigParser:
  config = ConfigParser.SafeConfigParser()
  config.add_section('TEA')

  # Open BART ConfigParser:
  Bconfig = ConfigParser.SafeConfigParser()
  Bconfig.read([cfg])

  # List of known TEA arguments:
  keys = ["maxiter",      "save_headers", "save_outputs", "doprint", "times",
          "location_TEA", "abun_file",    "location_out"]
  # TEA default values:
  defaults = ["100", "False", "False", "False", "False",
              "",    "None",  "None"]

  # Set TEA default arguments:
  for i in np.arange(len(keys)):
    if Bconfig.has_option("MCMC", keys[i]):
      config.set("TEA", keys[i], Bconfig.get("MCMC", keys[i]))
    else:
      config.set("TEA", keys[i], defaults[i])

  # Aliases:
  if Bconfig.has_option("MCMC", "abun_basic"):
    config.set("TEA", "abun_file", Bconfig.get("MCMC", "abun_basic"))
  if Bconfig.has_option("MCMC", "loc_dir"):
    config.set("TEA", "location_out", Bconfig.get("MCMC", "loc_dir"))
  # Default TEA dir:
  if config.get("TEA", "location_TEA") == "":
    config.set("TEA", "location_TEA", TEAdir)

  # For completion:
  config.add_section('PRE-ATM')
  config.set("PRE-ATM", "PT_file",        "None")
  config.set("PRE-ATM", "pre_atm_name",   "None")
  config.set("PRE-ATM", "input_elem",     "None")
  config.set("PRE-ATM", "output_species", "None")

  # Write TEA configuration file:
  with open("TEA.cfg", 'w') as configfile:
    config.write(configfile)
