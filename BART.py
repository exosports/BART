#! /usr/bin/env python

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

import sys, os, re, shutil, time, subprocess
import argparse, ConfigParser
import numpy as np

# Directory of BART.py file:
BARTdir = os.path.dirname(os.path.realpath(__file__))
TEAdir     = BARTdir + "/modules/TEA/"
MC3dir     = BARTdir + "/modules/MCcubed/src"
Transitdir = BARTdir + "/modules/transit/"

# Add path to submodules:
sys.path.append(BARTdir + "/code")
sys.path.append(MC3dir)

# Import submodules:
import mcutils   as mu
import makeP     as mp
import makeAbun  as ma
import InitialPT as ipt
#import PT        as pt
import makeatm   as mat
import makecfg   as mc

def main():
  """
  One function to run them all.

  Developer Team:
  ---------------
  Patricio Cubillos   pcubillos@fulbrightmail.org
  Jasmina Blecic      jasmina@physics.ucf.edu
  Joseph Harrington   jh@physics.ucf.edu
  Madison Stemm       astromaddie@gmail.com  (FINDME)

  Modification History:
  ---------------------
  2014-07-25  Jasmina   Initial version.
  2014-08-15  Patricio  put code into main() function.
  2014-08-18  Patricio  Merged with MC3 module.  Added flag to sort
                        read/execute steps.
  2014-09-20  Jasmina   Made call to makeRadius() function. Added progress
                        statements.
  2014-10-12  Jasmina   Updated to new TEA structure.
  2014-12-13  patricio  Added Opacity calculation step (through Transit), 
                        added flags to break after TEA or Opacity calculation.
  """

  mu.msg(1,
     "\n======= Bayesian Atmospheric Radiative Transfer (BART) ==============="
     "\nA code to infer planetary atmospheric properties based on observed  "
     "\nspectroscopic information."
   "\n\nCopyright (C) 2014 University of Central Florida. All rights reserved."
   "\n\nDevelopers contact:  Patricio Cubillos  pcubillos@fulbrightmail.org"
     "\n                     Jasmina Blecic     jasmina@physics.ucf.edu"
     "\n                     Joseph Harrington  jh@physics.ucf.edu"
     "\n======================================================================")

  mu.msg(1, "\nInitialization:")

  # Parse the config file from the command line:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--config_file",
                       help="Configuration file", metavar="FILE")
  cparser.add_argument("--justTEA",               action='store_true',
                       help="Run only TEA.")
  cparser.add_argument("--justOpacity",           action='store_true',
                       help="Run only Transit to generate the Opacity table.")
  cparser.add_argument("--resume",                action='store_true',
                       help="Resume a previous run.")
  # Remaining_argv contains all other command-line-arguments:
  args, remaining_argv = cparser.parse_known_args()

  # Take configuration file from command-line:
  cfile = args.config_file

  if cfile is not None and not os.path.isfile(cfile):
    mu.error("Configuration file: '{:s}' not found.".format(cfile))
  else: # Search for a .cfg file in current directory:
    files = os.listdir('.')
    try:
      #cfile = re.search("[\n]?(.+\.cfg)[\n]?", "\n".join(files)).groups(1)[0]
      cfile = 'BART.cfg'
    except:
      mu.error("No configuration file specified.")
  mu.msg(1, "The configuration file is: '{:s}'.".format(cfile), 2)

  # FINDME: Put this into an argparse object
  # Open config file:
  config = ConfigParser.RawConfigParser({})
  config.read(cfile)

  # Read BART.cfg parameters:
  cfgsec = "MCMC"
  # Directories and files:
  loc_dir     = config.get(cfgsec, 'loc_dir')
  tep_name    = config.get(cfgsec, 'tep_name')

  # STEP 1: Pressure array:
  n_layers    = config.getint(cfgsec,     'n_layers')
  p_top       = config.getfloat(cfgsec,   'p_top')
  p_bottom    = config.getfloat(cfgsec,   'p_bottom')
  log         = config.getboolean(cfgsec, 'log')
  press_file  = config.get(cfgsec,        'press_file')

  # STEP 2 Elemental abundances:
  abun_file   = config.get(cfgsec, 'abun_file')
  abun_basic  = config.get(cfgsec, 'abun_basic')
  solar_times = config.getfloat(cfgsec, 'solar_times')
  COswap      = config.getboolean(cfgsec, 'COswap')

  # STEP 3: Temperature profile:
  PTtype = config.get(cfgsec, 'PTtype')
  PTinit = np.asarray(config.get(cfgsec, 'PTinit').split(), np.double)

  # STEP 4: Elemental-abundances profile:
  in_elem     = config.get(cfgsec, 'in_elem')
  out_spec    = config.get(cfgsec, 'out_spec')
  preatm_file = config.get(cfgsec, 'preatm_file')

  # STEP 5: Atmospheric (species) file:
  #output_dir = config.get(cfgsec, 'output_dir')
  atmfile = config.get(cfgsec, 'atmfile')

  # Flag to break after TEA:
  TEAbreak = args.justTEA
  # Flag to break after the opacity calculation:
  Transitbreak = args.justOpacity

  # Transit variables:
  linedb  = config.get(cfgsec, 'linedb')
  cia     = config.get(cfgsec, 'cia')
  tconfig = config.get(cfgsec, 'config')
  opacity = config.get(cfgsec, 'opacityfile')

  # Make output directory:
  # Make a subdirectory with the date and time
  dirfmt = loc_dir + "%4d-%02d-%02d_%02d:%02d:%02d"
  date_dir = dirfmt % time.localtime()[0:6]
  # FINDME: Temporary hack:
  date_dir = os.path.normpath(loc_dir) + "/"
  if not os.path.isabs(date_dir):
    date_dir = os.getcwd() + "/" + date_dir
  mu.msg(1, "Output folder: '{:s}'".format(date_dir), 2)
  try:
    os.mkdir(date_dir)
  except OSError, e:
    if e.errno == 17: # FINDME: Allow overwritting while we debug
      pass
    else:
      mu.error("Cannot create folder '{:s}'. {:s}.".format(date_dir,
                                                     os.strerror(e.errno)))
  # Copy files to date dir:
  # BART configuration file:
  shutil.copy2(cfile, date_dir)
  # TEP file:
  if not os.path.isfile(tep_name):
    mu.error("Tepfile ('{:s}') Not found.".format(tep_name))
  else:
    shutil.copy2(tep_name, date_dir + os.path.basename(tep_name))


  # Check if files already exist:
  runMCMC = 0  # Flag that indicate which steps to run
  # Atmospheric file:
  if os.path.isfile(atmfile):
    shutil.copy2(atmfile, date_dir + os.path.basename(atmfile))
    mu.msg(1, "Atmospheric file copied from: '{:s}'.".format(atmfile),2)
    runMCMC |= 8
  # Pre-atmospheric file:
  if os.path.isfile(preatm_file):
    shutil.copy2(preatm_file, date_dir + os.path.basename(preatm_file))
    mu.msg(1, "Pre-atmospheric file copied from: '{:s}'.".format(preatm_file),2)
    runMCMC |= 4
  # Elemental-abundances file:
  if os.path.isfile(abun_file):
    shutil.copy2(abun_file, date_dir + os.path.basename(abun_file))
    mu.msg(1, "Elemental abundances file copied from: '{:s}'.".format(
                                                                 abun_file), 2)
    runMCMC |= 2
  # Pressure file:
  if os.path.isfile(press_file):
    shutil.copy2(press_file, date_dir + os.path.basename(press_file))
    mu.msg(1, "Pressure file copied from: '{:s}'.".format(press_file), 2)
    runMCMC |= 1


  # Generate files as needed:
  if runMCMC < 1:  # Pressure file
    press_file = date_dir + press_file
    mp.makeP(n_layers, p_top, p_bottom, press_file, log)
    mu.msg(1, "Created new pressure file.", 2)

  if runMCMC < 2:  # Elemental-abundances file
    abun_file = date_dir + abun_file
    ma.makeAbun(abun_basic, abun_file, solar_times, COswap)
    mu.msg(1, "Created new elemental abundances file.", 2)

  if runMCMC < 4:  # Pre-atmospheric file
    # Calculate the temperature profile:
    temp = ipt.initialPT2(PTinit, press_file, PTtype, tep_name)
    # Choose a pressure-temperature profile
    mu.msg(1, "\nChoose temperature and pressure profile:", 2)
    raw_input("  Press enter to continue, or quit and choose other initial "
              "PT parameters.")
    preatm_file = date_dir + preatm_file
    mat.make_preatm(tep_name, press_file, abun_file, in_elem, out_spec,
                  preatm_file, temp)
    mu.msg(1, "Created new pre-atmospheric file.", 2)

  if runMCMC < 8:  # Atmospheric file
    # Generate the TEA configuration file:
    mc.makeTEA(cfile, TEAdir)
    # Call TEA to calculate the atmospheric file:
    TEAcall = TEAdir + "tea/runatm.py"
    TEAout  = os.path.splitext(atmfile)[0]  # Remove extension
    # Execute TEA:
    mu.msg(1, "\nExecute TEA:")
    proc = subprocess.Popen([TEAcall, preatm_file, 'TEA'])
    proc.communicate()

    shutil.copy2(date_dir+"TEA/results/TEA.tea", date_dir+atmfile) 
    atmfile = date_dir + atmfile
    # Add radius array:
    mat.makeRadius(in_elem, out_spec, atmfile, abun_file, tep_name)
    mu.msg(1, "Added radius column to TEA atmospheric file.", 2)
    # Re-format file for use with transit:
    mat.reformat(atmfile)
    mu.msg(1, "Atmospheric file reformatted for Transit.", 2)

  if TEAbreak:
    mu.msg(1, "~~ BART End (after TEA) ~~")
    return

  # Make transit configuration file:
  mc.makecfg(cfile, date_dir, atmfile)

  # Generate the opacity file if it doesn't exist:
  if (not os.path.isfile(opacity) and 
      not os.path.isfile(date_dir + os.path.normpath(opacity))):
    mu.msg(1, "Transit call to generate the Opacity grid table.")
    Tcall = Transitdir + "/transit/transit"
    subprocess.call(["{:s} -c {:s} --justOpacity".format(Tcall, tconfig)],
                    shell=True, cwd=date_dir)

  if Transitbreak:
    mu.msg(1, "~~ BART End (after Transit) ~~")
    return

  # Run the MCMC:
  mu.msg(1, "\nStart MCMC:")
  MC3call = MC3dir + "/mccubed.py"
  subprocess.call(["mpiexec {:s} -c {:s}".format(MC3call, cfile)], shell=True,
                  cwd=date_dir)
  mu.msg(1, "~~ BART End ~~")


if __name__ == "__main__":
  main()
