#! /usr/bin/env python

import sys, os, re, shutil, time, subprocess
import argparse, ConfigParser
import numpy as np

# Directory of BART.py file:
BARTdir = os.path.dirname(os.path.realpath(__file__))
TEAdir  = BARTdir + "/modules/TEA/"
MC3dir  = BARTdir + "/modules/MCcubed/src"

# Add path to submodules:
sys.path.append(BARTdir + "/code")
sys.path.append(TEAdir)
sys.path.append(MC3dir)

# Import submodules:
import mcutils   as mu
import makeP     as mp
import makeAbun  as ma
import InitialPT as ipt
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
  """
  mu.msg(1, "\nThis is BART!")
  mu.msg(1, "\nInitialization:")

  # Parse the config file from the command line:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--config_file",
                       help="Configuration file", metavar="FILE")
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

  # STEP 1 pressure:
  n_layers    = config.getint(cfgsec,     'n_layers')
  p_top       = config.getfloat(cfgsec,   'p_top')
  p_bottom    = config.getfloat(cfgsec,   'p_bottom')
  log         = config.getboolean(cfgsec, 'log')
  press_file  = config.get(cfgsec,        'press_file')

  # STEP 2 abundances file:
  abun_file   = config.get(cfgsec, 'abun_file')
  abun_basic  = config.get(cfgsec, 'abun_basic')
  solar_times = config.getfloat(cfgsec, 'solar_times')
  COswap      = config.getboolean(cfgsec, 'COswap')

  # STEP 3 initial PT profile:
  p3     = config.getfloat(cfgsec, 'p3')
  p1     = config.getfloat(cfgsec, 'p1')
  a2     = config.getfloat(cfgsec, 'a2')
  a1     = config.getfloat(cfgsec, 'a1')
  T3_fac = config.getfloat(cfgsec, 'T3_fac')

  # STEP 4 makeatm:
  in_elem     = config.get(cfgsec, 'in_elem')
  out_spec    = config.get(cfgsec, 'out_spec')
  preatm_file = config.get(cfgsec, 'preatm_file')

  # STEP 5 run TEA:
  #output_dir = config.get(cfgsec, 'output_dir')
  atmfile = config.get(cfgsec, 'atmfile')

  # Transit variables:
  linedb = config.get(cfgsec, 'linedb')
  cia    = config.get(cfgsec, 'cia')

  # Make output directory:
  # Make a subdirectory with the date and time
  dirfmt = loc_dir + "%4d-%02d-%02d_%02d:%02d:%02d"
  date_dir = dirfmt % time.localtime()[0:6]
  # FINDME: Temporary hack:
  date_dir = os.path.normpath(loc_dir) + "/"
  if not os.path.isabs(date_dir):
    date_dir = os.getcwd() + "/" + date_dir + "/"
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
  shutil.copy2(cfile, date_dir)     # Configuration file
  if not os.path.isfile(tep_name):  # TEP file
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
    temp = ipt.initialPT(date_dir, tep_name, press_file, a1, a2, p1, p3, T3_fac)
    # Choose a pressure-temperature profile
    mu.msg(1, "\nChoose temperature and pressure profile:", 2)
    raw_input("  Press enter to continue, or quit and choose other initial "
              "PT parameters.")
    preatm_file = date_dir + preatm_file
    mat.make_preatm(tep_name, press_file, abun_file, in_elem, out_spec,
                  preatm_file, temp)
    mu.msg(1, "Created new pre-atmospheric file.", 2)

  if runMCMC < 8:  # Atmospheric file
    # Call TEA to calculate the atmospheric file:
    TEAcall = TEAdir + "tea/runatm.py"
    TEAout  = os.path.splitext(atmfile)[0]  # Remove extension
    #preatm_file = '/home/jasmina/BART-test/modules/TEA/doc/examples/multiTP/inputs/multiTP_Example.atm'    # temporal test case
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

  # Make transit configuration file:
  mc.makecfg(cfile, date_dir, atmfile)

  # Run the MCMC:
  mu.msg(1, "\nStart MCMC:")
  MC3call = MC3dir + "/mccubed.py"
  subprocess.call(["mpiexec {:s} -c {:s}".format(MC3call, cfile)], shell=True,
                  cwd=date_dir)
  mu.msg(1, "~~  The End  ~~", 2)


if __name__ == "__main__":
  main()
