#! /usr/bin/env python

# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

import sys, os, re, shutil, time, subprocess
import argparse, ConfigParser
import numpy as np

# Directory of BART.py file:
BARTdir = os.path.dirname(os.path.realpath(__file__))
TEAdir     = BARTdir + "/modules/TEA/"
MC3dir     = BARTdir + "/modules/MCcubed/"
Transitdir = BARTdir + "/modules/transit/"

# Add path to submodules and import:
sys.path.append(BARTdir + "/code")
import makeP     as mp
import InitialPT as ipt
import makeatm   as mat
import makecfg   as mc
import bestFit   as bf
import cf        as cf

sys.path.append(MC3dir)
import MCcubed.utils as mu

def main():
  """
  One function to run them all.
  """

  mu.msg(1,
     "\n======= Bayesian Atmospheric Radiative Transfer (BART) ==============="
     "\nA code to infer planetary atmospheric properties based on observed  "
     "\nspectroscopic information."
   "\n\nCopyright (C) 2015-2016 University of Central Florida."
     "\nAll rights reserved."
   "\n\nContact:  Patricio Cubillos  patricio.cubillos[at]oeaw.ac.at"
     "\n          Jasmina Blecic     jasmina[at]physics.ucf.edu"
     "\n          Joseph Harrington  jh[at]physics.ucf.edu"
     "\n======================================================================")

  mu.msg(1, "\nInitialization:")

  # Parse the config file from the command line:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--config_file",
                       help="Configuration file", metavar="FILE")

  # Parser for the MCMC arguments:
  parser = argparse.ArgumentParser(parents=[cparser])
  parser.add_argument("--justTEA",               action='store_true',
                       help="Run only TEA.")
  parser.add_argument("--justOpacity",           action='store_true',
                       help="Run only Transit to generate the Opacity table.")
  parser.add_argument("--resume",                action='store_true',
                       help="Resume a previous run.")
  # Directories and files options:
  group = parser.add_argument_group("Directories and files")
  group.add_argument("--loc_dir", dest="loc_dir",
           help="Output directory to store results [default: %(default)s]",
           type=str, action="store", default="outdir")
  group.add_argument("--tep_name", dest="tep_name",
           help="Transiting exoplanet file name.",
           type=str, action="store", default=None)
  group.add_argument("--logfile", dest="logfile",
           help="MCMC log file [default: %(default)s]",
           type=str, action="store", default="MCMC.log")
  # Pressure layers options:
  group = parser.add_argument_group("Layers pressure sampling")
  group.add_argument("--n_layers", dest="n_layers",
           help="Number of atmospheric layers [default: %(default)s]",
           type=int, action="store", default=100)
  group.add_argument("--p_top", dest="p_top",
           help="Pressure at the top of the atmosphere (bars) "
                "[default: %(default)s]",
           type=np.double, action="store", default=1.0e-5)
  group.add_argument("--p_bottom", dest="p_bottom",
           help="Pressure at the botom of the atmosphere (bars) "
                "[default: %(default)s]",
           type=np.double, action="store", default=100.0)
  group.add_argument("--log", dest="log",
           help="Use log (True) or linear (False) scale sampling "
                "[default: %(default)s]",
           type=eval, action="store", default=True)
  group.add_argument("--press_file", dest="press_file",
           help="Input/Output file with pressure array.",
           type=str, action="store", default=None)

  # Elemental abundance options:
  group = parser.add_argument_group("Elemental abundances")
  group.add_argument("--abun_basic", dest="abun_basic",
           help="Input elemental abundances file [default: %(default)s]",
           type=str, action="store", default="../BART/inputs/abundances_Asplund2009.txt")
  group.add_argument("--abun_file", dest="abun_file",
           help="Input/Output modified elemental abundances file",
           type=str, action="store", default=None)
  group.add_argument("--solar_times", dest="solar_times",
           help="Multiplication factor for metal-element abundances",
           type=int, action="store", default=1.0)
  group.add_argument("--COswap", dest="COswap",
           help="Swap C and O abundances if True [default: %(default)s]",
           type=eval, action="store", default=False)

  # Temperature profile options:
  group = parser.add_argument_group("Temperature profile")
  group.add_argument("--PTtype", dest="PTtype",
           help="Temperature profile model [default: %(default)s]",
           type=str, action="store", default="line", choices=("line","madhu"))
  group.add_argument("--PTinit", dest="PTinit",
           help="Temperature profile model parameters",
           type=mu.parray, action="store", default=None)

  # Atmospheric model options:
  group = parser.add_argument_group("Atmospheric model")
  group.add_argument("--in_elem", dest="in_elem",
           help="Input elements to consider in TEA [default: %(default)s]",
           type=str, action="store", default='H He C N O')
  group.add_argument("--out_spec", dest="out_spec",
           help="Output species to include in the atmospheric model "
                "[default: %(default)s]",
           type=str, action="store",
           default='H_g He_ref C_g N_g O_g H2_ref CO_g CO2_g CH4_g H2O_g')
  group.add_argument("--preatm_file", dest="preatm_file",
           help="Pre-atmospheric file with elemental abundances per layer "
                "[default: %(default)s]",
           type=str, action="store", default="elem.atm")
  group.add_argument("--atmfile", dest="atmfile",
           help="Atmospheric model file [default: %(default)s]",
           type=str, action="store", default="")
  group.add_argument("--uniform", dest="uniform",
           help="If not None, set uniform abundances with the specified "
                "values for each species in out_spec [default: %(default)s]",
           type=mu.parray, action="store", default=None)
  group.add_argument("--refpress", dest="refpress",
           help="Reference pressure level (bar) corresponding to the pressure"
                "at the planet radius [default: %(default)s]",
           type=float, action="store", default=0.1)

  # MCMC options:
  group = parser.add_argument_group("MCMC")
  group.add_argument("--params",  dest="params",
           help="Model-fitting parameters [default: %(default)s]",
           type=mu.parray, action="store", default=None)
  group.add_argument("--molfit",  dest="molfit", 
           help="Molecules fit [default: %(default)s]",
           type=mu.parray, action="store", default=None)
  group.add_argument("--Tmin",    dest="Tmin", 
           help="Lower Temperature boundary [default: %(default)s]",
           type=float,     action="store", default=400.0)
  group.add_argument("--Tmax",    dest="Tmax",
           help="Higher Temperature boundary [default: %(default)s]",
           type=float,     action="store", default=3000.0)
  group.add_argument("--quiet",   dest="quiet",
           help="Set verbosity level to minimum",
           action="store_true")
  group.add_argument("--stepsize", dest="stepsize",
           help="Parameters stepsize",
           type=mu.parray, action="store", default=None)
  group.add_argument("--burnin", dest="burnin",
           help="Number of burn-in iterations per chain",
           type=mu.parray, action="store", default=None)
  group.add_argument("--data", dest="data",
           help="Transit or eclipse depths",
           type=mu.parray, action="store", default=None)
  group.add_argument("--uncert", dest="uncert",
           help="Uncertanties on transit or eclipse depths",
           type=mu.parray, action="store", default=None)

  # Input converter options:
  group = parser.add_argument_group("Input Converter Options")
  group.add_argument("--tint", dest="tint",
           help="Internal temperature of the planet [default: %(default)s].",
           type=float, action="store", default=100.0)

  # Output-Converter Options:
  group = parser.add_argument_group("Output Converter Options")
  group.add_argument("--filter",                 action="store",
           help="Waveband filter name [default: %(default)s]",
           dest="filter",   type=mu.parray, default=None)
  group.add_argument("--kurucz_file",           action="store",
           help="Stellar Kurucz file [default: %(default)s]",
           dest="kurucz",   type=str,       default=None)
  group.add_argument("--solution",                    action="store",
           help="Solution geometry [default: %(default)s]",
           dest="solution", type=str,       default="None",
           choices=('transit', 'eclipse'))

  # Transit options:
  group = parser.add_argument_group("Transit variables")
  group.add_argument("--tconfig", dest="tconfig",
           help="Transit configuration file [default: %(default)s]",
           type=str, action="store", default="transit.cfg")
  group.add_argument("--opacityfile", dest="opacityfile",
           help="Opacity table file [default: %(default)s]",
           type=str, action="store", default=None)
  group.add_argument("--outflux", dest="outflux",
           help="Output with flux values [default: %(default)s]",
           type=str, action="store", default=None)
  group.add_argument("--outmod", dest="outmod",
           help="Output with modulation values [default: %(default)s]",
           type=str, action="store", default=None)
  group.add_argument("--shareOpacity", dest="shareOpacity",
           help="If True, use shared memory for the Transit opacity file "
                "[default: %(default)s]",
           type=eval, action="store", default=True)


  # Remaining_argv contains all other command-line-arguments:
  cargs, remaining_argv = cparser.parse_known_args()
  # Get only the arguments defined above:
  known, unknown = parser.parse_known_args(remaining_argv)


  # Get configuration file from command-line:
  cfile = cargs.config_file
  # Default:
  if cfile is None:
    cfile = "./BART.cfg"
  # Always require a configuration file:
  if not os.path.isfile(cfile):
    mu.error("Configuration file: '{:s}' not found.".format(cfile))

  # Read values from configuration file:
  config = ConfigParser.SafeConfigParser()
  config.optionxform = str  # This one enable Uppercase in arguments
  config.read([cfile])
  defaults = dict(config.items("MCMC"))
  mu.msg(1, "The configuration file is: '{:s}'.".format(cfile), indent=2)

  # Set the defaults from the configuration file:
  parser.set_defaults(**defaults)
  # Set values from command line:
  args, unknown = parser.parse_known_args(remaining_argv)

  # Unpack the variables from args:
  variables = dir(args)
  for var in dir(known):
    if not var.startswith("_"):
      exec("{:s} = args.{:s}".format(var, var))

  # Make output directory:
  # Make a subdirectory with the date and time
  dirfmt = loc_dir + "%4d-%02d-%02d_%02d:%02d:%02d"
  date_dir = dirfmt % time.localtime()[0:6]
  # FINDME: Temporary hack (temporary?):
  date_dir = os.path.normpath(loc_dir) + "/"
  if not os.path.isabs(date_dir):
    date_dir = os.getcwd() + "/" + date_dir
  mu.msg(1, "Output folder: '{:s}'".format(date_dir), indent=2)
  try:
    os.mkdir(date_dir)
  except OSError, e:
    if e.errno == 17: # Allow overwritting while we debug
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
    atmfile = os.path.realpath(atmfile)
    shutil.copy2(atmfile, date_dir + os.path.basename(atmfile))
    mu.msg(1, "Atmospheric file copied from: '{:s}'.".format(atmfile),indent=2)
    runMCMC |= 8
  # Pre-atmospheric file:
  if os.path.isfile(preatm_file):
    preatm_file = os.path.realpath(preatm_file)
    shutil.copy2(preatm_file, date_dir + os.path.basename(preatm_file))
    mu.msg(1, "Pre-atmospheric file copied from: '{:s}'.".format(preatm_file),
           indent=2)
    runMCMC |= 4
  # Elemental-abundances file:
  if abun_file is not None and os.path.isfile(abun_file):
    shutil.copy2(abun_file, date_dir + os.path.basename(abun_file))
    mu.msg(1, "Elemental abundances file copied from: '{:s}'.".
              format(abun_file), indent=2)
    runMCMC |= 2
  # Pressure file:
  if press_file is not None and os.path.isfile(press_file):
    shutil.copy2(press_file, date_dir + os.path.basename(press_file))
    mu.msg(1, "Pressure file copied from: '{:s}'.".format(press_file), indent=2)
    runMCMC |= 1


  # Generate files as needed:
  if runMCMC < 1:  # Pressure file
    press_file = date_dir + press_file
    mp.makeP(n_layers, p_top, p_bottom, press_file, log)
    mu.msg(1, "Created new pressure file.", indent=2)

  # Make uniform-abundance profiles if requested:
  if uniform is not None and runMCMC < 8:
    # Calculate the temperature profile:
    temp = ipt.initialPT2(date_dir, PTinit, press_file, PTtype, tep_name)
    # Generate the uniform-abundance profiles file:
    mat.uniform(date_dir + atmfile, press_file, abun_basic, tep_name,
               out_spec, uniform, temp, refpress)
    atmfile = date_dir + atmfile
    # Update the runMCMC flag to skip upcoming steps:
    runMCMC |= 8

  if runMCMC < 2:  # Elemental-abundances file
    abun_file = date_dir + abun_file
    mu.msg(1, "CO swap: {}".format(COswap), indent=2)
    mat.makeAbun(abun_basic, abun_file, solar_times, COswap)
    mu.msg(1, "Created new elemental abundances file.", indent=2)

  if runMCMC < 4:  # Pre-atmospheric file
    # Calculate the temperature profile:
    temp = ipt.initialPT2(date_dir, PTinit, press_file, PTtype, tep_name)
    # Choose a pressure-temperature profile
    mu.msg(1, "\nChoose temperature and pressure profile:", indent=2)
    raw_input("  open Initial PT profile figure and\n" 
              "  press enter to continue or quit and choose other initial "
              "PT parameters.")
    preatm_file = date_dir + preatm_file
    mat.make_preatm(tep_name, press_file, abun_file, in_elem, out_spec,
                  preatm_file, temp)
    mu.msg(1, "Created new pre-atmospheric file.", indent=2)

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
    mat.makeRadius(out_spec, atmfile, abun_file, tep_name, refpress)
    mu.msg(1, "Added radius column to TEA atmospheric file.", indent=2)
    # Re-format file for use with transit:
    mat.reformat(atmfile)
    mu.msg(1, "Atmospheric file reformatted for Transit.", indent=2)

  if justTEA:
    mu.msg(1, "~~ BART End (after TEA) ~~")
    return

  # Make the MC3 configuration file:
  MCMC_cfile = os.path.realpath(loc_dir) + "/MCMC_" + os.path.basename(cfile)
  mc.makeMCMC(cfile, MCMC_cfile, logfile)
  # Make transit configuration file:
  mc.makeTransit(MCMC_cfile, tep_name, shareOpacity)

  # Generate the opacity file if it doesn't exist:
  if not os.path.isfile(opacityfile):
    mu.msg(1, "Transit call to generate the Opacity grid table.")
    Tcall = Transitdir + "/transit/transit"
    subprocess.call(["{:s} -c {:s} --justOpacity".format(Tcall, tconfig)],
                    shell=True, cwd=date_dir)
  else:
    mu.msg(1, "\nTransit copies the existing opacity file from:\n '{:s}'.".
                 format(opacityfile), indent=2)
    shutil.copy2(opacityfile, date_dir + os.path.basename(opacityfile))

  if justOpacity:
    mu.msg(1, "~~ BART End (after Transit opacity calculation) ~~")
    return


  # Run the MCMC:
  MC3call = MC3dir + "/MCcubed/mccubed.py"
  subprocess.call(["mpiexec {:s} -c {:s}".format(MC3call, MCMC_cfile)],
                  shell=True, cwd=date_dir)

  # Run best-fit Transit call
  mu.msg(1, "\nTransit call with the best-fitting values.")
 
  # MCcubed output file
  MCfile = date_dir + logfile
  
  # Call bestFit submodule and make new bestFit_tconfig.cfg
  bf.callTransit(atmfile, tep_name, MCfile, stepsize, molfit, solution,
                 refpress, tconfig, date_dir, params, burnin, abun_basic)

  # Best-fit tconfig
  bestFit_tconfig = date_dir + 'bestFit_tconfig.cfg'

  # Call Transit with the best-fit tconfig
  Tcall = Transitdir + "/transit/transit"
  subprocess.call(["{:s} -c {:s}".format(Tcall, bestFit_tconfig)],
                   shell=True, cwd=date_dir)

  # Plot best-fit eclipse or modulation spectrum, depending on solution:
  if solution == 'eclipse':
    # Plot best-fit eclipse spectrum
    bf.plot_bestFit_Spectrum(filter, kurucz, tep_name, solution, outflux,
                             data, uncert, date_dir)
  elif solution == 'transit':
    # Plot best-fit transit spectrum
    bf.plot_bestFit_Spectrum(filter, kurucz, tep_name, solution, outmod,
                             data, uncert, date_dir)


  # Plot abundance profiles
  bf.plotabun(date_dir, 'bestFit.atm', molfit)
  
  # Run Transit with unlimited 'toomuch' argument for contribution
  # function calculation:
  mu.msg(1, "\nTransit call for contribution functions calculation.")

  # Make cf_tconfig
  cf.cf_tconfig(date_dir)

  # Contribution functions tconfig
  cf_tconfig = date_dir + 'cf_tconfig.cfg'

  # Call Transit with the cf_tconfig
  Tcall = Transitdir + "/transit/transit"
  subprocess.call(["{:s} -c {:s}".format(Tcall, cf_tconfig)],
                    shell=True, cwd=date_dir)

  # Calculate contribution functions and plot them
  mu.msg(1, "Calculating contribution functions ...", indent=2)
  bestFit_atmfile = date_dir + 'bestFit.atm'
  cf.cf(date_dir, bestFit_atmfile, filter)

  mu.msg(1, "~~ BART End ~~")


if __name__ == "__main__":
  main()
