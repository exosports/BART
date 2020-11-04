#! /usr/bin/env python

# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

import sys, os, re, shutil, time, subprocess
import argparse
from six.moves import configparser
import six
if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser
import numpy as np

# Directory of BART.py file:
BARTdir = os.path.dirname(os.path.realpath(__file__))
TEAdir     = os.path.join(BARTdir, "modules", "TEA", "")
MC3dir     = os.path.join(BARTdir, "modules", "MCcubed", "")
Transitdir = os.path.join(BARTdir, "modules", "transit", "")

# Add path to submodules and import:
sys.path.append(os.path.join(BARTdir, "code"))
import makeP       as mp
import InitialPT   as ipt
import        PT   as  pt
import makeatm     as mat
import makecfg     as mc
import bestFit     as bf
import cf          as cf
import mc3plots    as mcp

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
  parser.add_argument("--justTEA", dest="justTEA", action='store_true',
                      help="Run only TEA.", default=False)
  parser.add_argument("--justOpacity", dest="justOpacity", action='store_true',
                      help="Run only Transit to generate the Opacity table.", 
                      default=False)
  parser.add_argument("--justPlots", dest="justPlots", action='store_true',
                      help="Remakes plots of BART output.", default=False)
  parser.add_argument("--resume", dest="resume", action='store_true',
                      help="Resume a previous run.", default=False)
  # Directories and files options:
  group = parser.add_argument_group("Directories and files")
  group.add_argument("--fext",    dest="fext", 
           help="File extension for plots [default: %(default)s]", 
           type=str, action="store", default=".png")
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
           type=str, action="store", default="line", 
           choices=("line","madhu_noinv","madhu_inv","iso"))
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
                " at the planet radius [default: %(default)s]",
           type=float, action="store", default=0.1)
  group.add_argument("--cloudtop",           action="store",
           help="Cloud deck top pressure [default: %(default)s]",
                     dest="cloudtop",   type=float, default=None)
  group.add_argument("--scattering",           action="store",
           help="Rayleigh scattering [default: %(default)s]",
                     dest="scattering", type=float, default=None)

  # MCMC options:
  group = parser.add_argument_group("MCMC")
  group.add_argument("--params",  dest="params",
           help="Model-fitting parameters [default: %(default)s]",
           type=mu.parray, action="store", default=None)
  group.add_argument("--parnames", dest="parnames", 
           help="Labels for model-fitting parameters [default: %(default)s]", 
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
  group.add_argument("--nchains",  dest="nchains", 
           help="Number of parallel chains for MCMC", 
           type=int,       action="store", default=10)
  group.add_argument("--walk",  dest="walk", 
           help="MCMC algorithm", 
           type=str, action="store", default="snooker", 
           choices=('snooker', 'mrw', 'demc', 'unif'))
  group.add_argument("--stepsize", dest="stepsize",
           help="Parameters stepsize",
           type=mu.parray, action="store", default=None)
  group.add_argument("--burnin", dest="burnin",
           help="Number of burn-in iterations per chain",
           type=int,       action="store", default=None)
  group.add_argument("--thinning", dest="thinning",
           help="Thinning factor of the chains (use every thinning-th "
                 "iteration) used in the GR test and plots", 
           type=int,       action="store", default=1)
  group.add_argument("--data", dest="data",
           help="Transit or eclipse depths",
           type=mu.parray, action="store", default=None)
  group.add_argument("--uncert", dest="uncert",
           help="Uncertanties on transit or eclipse depths",
           type=mu.parray, action="store", default=None)
  group.add_argument("--savemodel", dest="savemodel", 
           help="Filename to save out models.",
           type=str, action="store", default=None)
  group.add_argument("--modelper", dest="modelper", 
           help="Determines how to split MC3's `savemodel`. " + \
                "0 makes no split, >0 sets the # of iterations per split. " + \
                "If nchains=10 and modelper=5, it will save every 50 " + \
                "models to a new .NPY file.",
           type=int,       action="store", default=0)
  group.add_argument("--plots", dest="plots", 
                     help="Determines whether to produce plots.", 
                     type=bool, action="store", default=True)

  # Input converter options:
  group = parser.add_argument_group("Input Converter Options")
  group.add_argument("--tint", dest="tint",
           help="Internal temperature of the planet [default: %(default)s].",
           type=float, action="store", default=100.0)
  group.add_argument("--tint_type", dest="tint_type",
           help="Method to evaluate `tint`. Options: const or thorngren. " + \
                "[default: %(default)s].",
           type=str,   action="store", default='const', 
           choices=("const","thorngren"))

  # Output-Converter Options:
  group = parser.add_argument_group("Output Converter Options")
  group.add_argument("--filters",                action="store",
           help="Waveband filter names [default: %(default)s]",
           dest="filters",   type=mu.parray, default=None)
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
  group.add_argument("--outspec", dest="outspec",
           help="Output spectrum filename [default: %(default)s]",
           type=str, action="store", default="outspec.dat")
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
    cfile = os.path.join(os.getcwd(), "BART.cfg")
  # Always require a configuration file:
  if not os.path.isfile(cfile):
    mu.error("Configuration file: '{:s}' not found.".format(cfile))

  # Read values from configuration file:
  config = ConfigParser()
  config.optionxform = str  # This one enable Uppercase in arguments
  config.read([cfile])
  defaults = dict(config.items("MCMC"))
  mu.msg(1, "The configuration file is: '{:s}'.".format(cfile), indent=2)

  # Set the defaults from the configuration file:
  parser.set_defaults(**defaults)
  # Set values from command line:
  args, unknown = parser.parse_known_args(remaining_argv)

  # Unpack configuration-file/command-line arguments:
  justTEA     = args.justTEA
  justOpacity = args.justOpacity
  justPlots   = args.justPlots
  resume      = args.resume

  loc_dir  = args.loc_dir
  fext     = args.fext
  tep_name = args.tep_name
  logfile  = args.logfile

  n_layers   = args.n_layers
  p_top      = args.p_top
  p_bottom   = args.p_bottom
  log        = args.log
  press_file = args.press_file

  abun_basic  = args.abun_basic
  abun_file   = args.abun_file
  solar_times = args.solar_times
  COswap      = args.COswap
  cloud    = args.cloudtop
  rayleigh = args.scattering

  PTtype = args.PTtype
  PTinit = args.PTinit

  in_elem     = args.in_elem
  out_spec    = args.out_spec
  preatm_file = args.preatm_file
  atmfile     = args.atmfile
  uniform     = args.uniform
  refpress    = args.refpress

  params    = args.params
  parnames  = args.parnames
  molfit    = args.molfit
  Tmin      = args.Tmin
  Tmax      = args.Tmax
  quiet     = args.quiet
  nchains   = args.nchains
  walk      = args.walk
  stepsize  = args.stepsize
  burnin    = args.burnin
  thinning  = args.thinning
  data      = args.data
  uncert    = args.uncert
  savemodel = args.savemodel
  modelper  = args.modelper
  plots     = args.plots

  tint      = args.tint
  tint_type = args.tint_type

  filters  = args.filters
  kurucz   = args.kurucz
  solution = args.solution

  tconfig      = args.tconfig
  opacityfile  = args.opacityfile
  outspec      = args.outspec
  shareOpacity = args.shareOpacity

  # Unpack the variables from args:
  '''
  argd = {}
  for key, val in vars(args).items():
    if type(val) == str and val in ['True', 'False', 'None']:
      if val == 'True':
        argd.update({key:True})
      elif val == 'False':
        argd.update({key:False})
      elif val == 'None':
        argd.update({key:None})
    else:
      argd.update({key:val})
  vars(sys.modules[__name__]).update(argd)
  '''

  # Dictionary of functions to calculate temperature for PTtype
  PTfunc = {'iso'         : pt.PT_iso,
            'line'        : pt.PT_line, 
            'madhu_noinv' : pt.PT_NoInversion,
            'madhu_inv'   : pt.PT_Inversion}

  # Check that the user gave a valid PTtype:
  if PTtype not in PTfunc.keys():
    print("The specified 'PTtype' is not valid. Options are 'line', " + \
          "'madhu_noinv', 'madhu_inv', or 'iso'. Please try again.")
    sys.exit()

  # Check that out_spec and uniform are valid specifications
  if uniform is not None and len(uniform) != len(out_spec.split(' ')):
    print('The inputs for out_spec and uniform are not compatible.')
    diffuniout = len(uniform) - len(out_spec.split(' '))
    if diffuniout > 0:
      if diffuniout == 1:
        print('uniform has ' + str(diffuniout) + 'extra entry.')
      else:
        print('uniform has ' + str(diffuniout) + 'extra entries.')
    else:
      if diffuniout == -1:
        print('out_spec has ' + str(-1*diffuniout) + 'extra entry.')
      else:
        print('out_spec has ' + str(-1*diffuniout) + 'extra entries.')
    print('Please correct this and run again.')
    sys.exit()

  # Make output directory:
  # Make a subdirectory with the date and time
  dirfmt = loc_dir + "%4d-%02d-%02d_%02d:%02d:%02d"
  date_dir = dirfmt % time.localtime()[0:6]
  # FINDME: Temporary hack (temporary?):
  date_dir = os.path.join(os.path.normpath(loc_dir), "")
  if not os.path.isabs(date_dir):
    date_dir = os.path.join(os.getcwd(), date_dir)
  mu.msg(1, "Output folder: '{:s}'".format(date_dir), indent=2)
  try:
    os.mkdir(date_dir)
  except OSError as e:
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
  if justPlots:
    mu.msg(1, "\nRe-making output plots.", indent=0)
    runMCMC |= 16
  # Atmospheric file:
  if os.path.isfile(atmfile):
    fatmfile = os.path.realpath(atmfile)
    shutil.copy2(fatmfile, date_dir + os.path.basename(fatmfile))
    mu.msg(1, "Atmospheric file copied from: '{:s}'.".format(fatmfile),indent=2)
    runMCMC |= 8
  atmfile = date_dir + os.path.basename(atmfile)
  # Pre-atmospheric file:
  if os.path.isfile(preatm_file):
    fpreatm_file = os.path.realpath(preatm_file)
    shutil.copy2(fpreatm_file, date_dir + os.path.basename(fpreatm_file))
    mu.msg(1, "Pre-atmospheric file copied from: '{:s}'.".format(fpreatm_file),
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

  press_file = date_dir + os.path.basename(press_file)
  # Generate files as needed:
  if runMCMC < 1:  # Pressure file
    mp.makeP(n_layers, p_top, p_bottom, press_file, log)
    mu.msg(1, "Created new pressure file.", indent=2)

  # Make uniform-abundance profiles if requested:
  if uniform is not None and runMCMC < 8:
    # Calculate the temperature profile:
    temp = ipt.initialPT2(date_dir, PTinit,         press_file, 
                          PTtype,   PTfunc[PTtype], tep_name)
    # Generate the uniform-abundance profiles file:
    mat.uniform(atmfile,  press_file, abun_basic, tep_name,
                out_spec, uniform,    temp,       refpress)
    # Update the runMCMC flag to skip upcoming steps:
    runMCMC |= 8

  if runMCMC < 2:  # Elemental-abundances file
    mu.msg(1, "CO swap: {}".format(COswap), indent=2)
    mat.makeAbun(abun_basic, date_dir+abun_file, solar_times, COswap)
    mu.msg(1, "Created new elemental abundances file.", indent=2)

  abun_file = date_dir + abun_file

  if runMCMC < 4:  # Pre-atmospheric file
    # Calculate the temperature profile:
    temp = ipt.initialPT2(date_dir, PTinit,         press_file, 
                          PTtype,   PTfunc[PTtype], tep_name, 
                          tint_type=tint_type)
    # Choose a pressure-temperature profile
    mu.msg(1, "\nChoose temperature and pressure profile:", indent=2)
    raw_input("  open Initial PT profile figure and\n" 
              "  press enter to continue or quit and choose other initial "
              "PT parameters.")
    mat.make_preatm(tep_name, press_file, abun_file, 
                    in_elem, out_spec, preatm_file, temp)
    mu.msg(1, "Created new pre-atmospheric file.", indent=2)

  if runMCMC < 8:  # Atmospheric file
    # Generate the TEA configuration file:
    mc.makeTEA(cfile, TEAdir)
    # Call TEA to calculate the atmospheric file:
    TEAcall = os.path.join(TEAdir, "tea", "runatm.py")
    TEAout  = os.path.splitext(atmfile)[0]  # Remove extension
    # Execute TEA:
    mu.msg(1, "\nExecute TEA:")
    proc = subprocess.Popen([TEAcall, preatm_file, 'TEA'])
    proc.communicate()

    TEAres = os.path.join("TEA", "results", "TEA.tea")
    shutil.copy2(os.path.join(date_dir, TEAres), atmfile) 
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
  if runMCMC < 16: # MCMC
    MCMC_cfile = os.path.join(os.path.realpath(loc_dir), 
                              "MCMC_" + os.path.basename(cfile))
    mc.makeMCMC(cfile, MCMC_cfile, logfile)
    # Make transit configuration file:
    mc.makeTransit(MCMC_cfile, tep_name, shareOpacity)

    # Generate the opacity file if it doesn't exist:
    if not os.path.isfile(opacityfile):
      mu.msg(1, "Transit call to generate the Opacity grid table.")
      Tcall = os.path.join(Transitdir, "transit", "transit")
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
  if runMCMC < 16:
    MC3call = os.path.join(MC3dir, "MCcubed", "mccubed.py")
    subprocess.call(["mpiexec {:s} -c {:s}".format(MC3call, MCMC_cfile)],
                    shell=True, cwd=date_dir)

  if walk=='unif' and modelper > 0:
    # Clean up the output directory
    model_dir = os.path.join(date_dir, savemodel.replace('.npy', ''), '')
    # Make directory
    try:
      os.mkdir(model_dir)
    except OSError as e:
      if e.errno == 17: # Already exists
        pass
      else:
        mu.error("Cannot create folder '{:s}'. {:s}.".format(model_dir,
                                                       os.strerror(e.errno)))
    # Move model files to subdirectory
    subprocess.call(['mv {:s} {:s}'.format('*'.join(savemodel.split('.')), 
                                           model_dir)                     ], 
                    shell=True, cwd=date_dir)

  if plots and walk != 'unif':
    # Re-plot MCMC results in prettier format
    mcp.mc3plots('output.npy', burnin,   thinning, nchains, uniform, molfit, 
                 out_spec,     parnames, stepsize, date_dir, 
                 ["output_trace", "output_pairwise", 
                  "output_posterior"], fext)

    # Run best-fit Transit call
    mu.msg(1, "\nTransit call with the best-fitting values.")

    # MCcubed output file
    MCfile = date_dir + logfile
  
    # Call bestFit submodule: make new bestFit_tconfig.cfg, run best-fit Transit
    bf.callTransit(atmfile, tep_name, MCfile,  stepsize, molfit, 
                   cloud, rayleigh,
                   solution, refpress, tconfig, date_dir, burnin, 
                   abun_basic, PTtype, PTfunc[PTtype], 
                   tint, tint_type, filters, fext=fext)
    # Plot best-fit eclipse or modulation spectrum, depending on solution:
    bf.plot_bestFit_Spectrum(filters, kurucz, tep_name, solution, outspec,
                             data, uncert, date_dir, fext)

    bestFit_atmfile = 'bestFit.atm'

    # Plot abundance profiles
    bf.plotabun(date_dir, bestFit_atmfile, molfit, fext)
  
    mu.msg(1, "\nTransit call for contribution functions/transmittance.")
    # Run Transit with unlimited 'toomuch' argument:
    cf.cf_tconfig(date_dir)
    # Call Transit with the cf_tconfig
    cf_tconfig = date_dir + 'cf_tconfig.cfg'
    Tcall = os.path.join(Transitdir, "transit", "transit")
    subprocess.call(["{:s} -c {:s}".format(Tcall, cf_tconfig)],
                    shell=True, cwd=date_dir)

    # Calculate and plot contribution functions:
    if solution == "eclipse":
      # Compute contribution fucntions if this is a eclipse run:
      mu.msg(1, "Calculating contribution functions.", indent=2)
      ctfraw, ctf = cf.cf(date_dir, bestFit_atmfile, filters, fext)
    else:
      # Compute transmittance if this is a transmission run:
      mu.msg(1, "Calculating transmittance.", indent=2)
      ctf = cf.transmittance(date_dir, bestFit_atmfile, filters, fext)

    # Make a plot of MCMC profiles with contribution functions/transmittance
    bf.callTransit(atmfile, tep_name, MCfile, stepsize, molfit, 
                   cloud, rayleigh,
                   solution, refpress, tconfig, date_dir, burnin, 
                   abun_basic, PTtype, PTfunc[PTtype], 
                   tint, tint_type, filters, ctf, fext=fext)

  mu.msg(1, "~~ BART End ~~")


if __name__ == "__main__":
  main()
