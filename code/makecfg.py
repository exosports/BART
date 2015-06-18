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

import os, sys
import argparse, ConfigParser
import numpy as np
import scipy.constants as sc

import reader as rd
import constants as c

filedir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(filedir + "/../modules/MCcubed/src/")
import mcutils as mu


def makeTransit(cfile, tepfile):
  """
  Make the transit configuration file.

  Parameters:
  -----------
  cfile: String
     BART configuration file.
  tepfile: String
     A TEP file.
  """

  # Known transit arguments:
  known_args = ["radlow",  "radhigh",  "raddelt",  "radfct",
                "wllow",   "wlhigh",   "wlfct",
                "wnlow",   "wnhigh",   "wndelt",  "wnfct", "wnosamp",
                "tlow",    "thigh",    "tempdelt",
                "allowq",  "nwidth",
                "opacityfile",
                "molfile",
                "toomuch", "tauiso", "outtau", "taulevel", "modlevel",
                "starrad", "transparent",
                "solution", "raygrid",
                "gsurf", "refpress", "refradius",
                "cloudrad", "cloudfct", "cloudext",
                "verb",
                "savefiles",
                "outtoomuch", "outsample", "outmod", "outflux", "outintens",
                "linedb", "cia", "orbpars", "orbparsfct"]

  # Name of the configuration-file section:
  section = "MCMC"
  # Read BART configuration file:
  Bconfig = ConfigParser.SafeConfigParser()
  Bconfig.read([cfile])

  # Keyword names of the arguments in the BART configuration file:
  args = Bconfig.options(section)

  # transit configuration filename:
  tcfile = open(Bconfig.get(section, "tconfig"), "w")

  # Keyword for the atmospheric file is different in transit:
  tcfile.write("atm {:s}\n".format(Bconfig.get(section, "atmfile")))

  # Default molfile path:
  if "molfile" not in args:
    tcfile.write("molfile {:s}\n".format(
      os.path.realpath(filedir + "/../modules/transit/inputs/molecules.dat")))

  # Calculate gsurf and refradius from the tepfile:
  tep = rd.File(tepfile)
  mplanet = float(tep.getvalue('Mp')[0]) * c.Mjup
  rplanet = float(tep.getvalue('Rp')[0]) * c.Rjup

  # Planetary radius reference level in km:
  Bconfig.set(section, "refradius", "{:.2f}".format(rplanet * 1e-3))
  # Planetary surface gravity in cm/s2:
  Bconfig.set(section, "gsurf", "{:.1f}".format(100*sc.G*mplanet/rplanet**2))
  # Add these keywords:
  args = np.union1d(args, ["refradius", "gsurf"])

  #CIA flag
  ciawritten = 0

  # Print the known arguments to file:
  for key in np.intersect1d(args, known_args):
    values = Bconfig.get(section, key).split("\n")
    for val in values:
      #if the first cia entry, flip flag, print normally but without final \n
      if key == "cia" and ciawritten == 0:
        tcfile.write("{:s} {:s}".format(key, val))
        ciawritten = 1
      #if not first cia entry, print comma then new file
      elif key == "cia" and ciawritten == 1:
        tcfile.write(",{:s}".format(val))
      #if first not cia entry after cia entries, print \n to end cia entries, then print normally
      elif key != "cia" and ciawritten == 1:
        tcfile.write("\n{:s} {:s}\n".format(key, val))
        ciawritten = 0
      #print normally
      else:
        tcfile.write("{:s} {:s}\n".format(key, val))

  tcfile.close()


def makeMCMC(cfile, MCMC_cfile):
  """
  Reformat configuration file to remove relative paths.  This output 
  configuration file is used by the BART's MCMC program.

  Parameters:
  -----------
  cfile: String
     BART configuration file.
  MCMC_cfile: String
     Reformated configuration file.
  """

  # Name of the configuration-file section:
  section = "MCMC"
  # Open input BART configuration file:
  Bconfig = ConfigParser.SafeConfigParser()
  Bconfig.optionxform = str
  Bconfig.read([cfile])
  # Keyword names of the arguments in the BART configuration file:
  args = Bconfig.options(section)

  # Known arguments that may have a path:
  input_args = ["tep_name", "kurucz", "molfile", "filter", "linedb",
                "cia", "loc_dir"]
  output_args = ["tconfig", "atmfile", "opacityfile", "press_file",
               "abun_basic", "abun_file", "preatm_file", "output", "savemodel"]

  # Inputs should always exist:
  for arg in np.intersect1d(args, input_args):
    # Split multiple values (if more than one), get full path, join back:
    values = Bconfig.get(section, arg).split("\n")
    for v in np.arange(len(values)):
      values[v] = os.path.realpath(values[v])
    Bconfig.set(section, arg, "\n".join(values))

  # Outputs are stored/copied into loc_dir:
  for arg in np.intersect1d(args, output_args):
    Bconfig.set(section, arg,
                Bconfig.get(section, "loc_dir") + "/" +
                os.path.basename(Bconfig.get(section, arg)))

  # Func is a special case:
  funcvalues = Bconfig.get(section, "func").split()
  funcvalues[2] = os.path.realpath(funcvalues[2])
  Bconfig.set(section, "func", " ".join(funcvalues))

  # Params is another special case:
  params = Bconfig.get(section, "params")
  # It may or not be a file path:
  if os.path.isfile(params):
    Bconfig.set(section, "params", os.path.realpath(params))

  # Write the configuration file for use by MC3:
  with open(MCMC_cfile, 'w') as configfile:
    Bconfig.write(configfile)


def makeTEA(cfile, TEAdir):
  """
  Make a TEA configuration file.

  Parameters:
  -----------
  cfile: String
     BART configuration file
  TEAdir: String
     Default TEA directory.
  """
  # Open New ConfigParser:
  config = ConfigParser.SafeConfigParser()
  config.add_section('TEA')

  # Open BART ConfigParser:
  Bconfig = ConfigParser.SafeConfigParser()
  Bconfig.read([cfile])

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
