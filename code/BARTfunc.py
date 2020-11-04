#! /usr/bin/env python

# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

import sys, os
import argparse
from six.moves import configparser
import six
if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser
import numpy as np
import scipy.constants as sc
from mpi4py import MPI

import makeatm   as mat
import PT        as pt
import wine      as w
import reader    as rd
import constants as c

BARTdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
sys.path.append(os.path.join(BARTdir, "modules", "MCcubed", ""))
import MCcubed.utils as mu

Tdir = os.path.join(BARTdir, "modules", "transit", "transit", "python")
sys.path.append(Tdir)
import transit_module as trm


def main(comm):
  """
  This is a hacked version of MC3's func.py.
  This function directly call's the modeling function for the BART project.
  """
  # Parse arguments:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--config_file", 
                       help="Configuration file", metavar="FILE")
  # Remaining_argv contains all other command-line-arguments:
  args, remaining_argv = cparser.parse_known_args()

  # Get parameters from configuration file:
  cfile = args.config_file
  if cfile:
    config = ConfigParser()
    config.optionxform = str
    config.read([cfile])
    defaults = dict(config.items("MCMC"))
  else:
    defaults = {}
  parser = argparse.ArgumentParser(parents=[cparser])
  parser.add_argument("--func",      dest="func",      type=mu.parray, 
                                     action="store",  default=None)
  parser.add_argument("--indparams", dest="indparams", type=mu.parray, 
                                     action="store",   default=[])
  parser.add_argument("--params",    dest="params",    type=mu.parray,
                                     action="store",   default=None,
                      help="Model-fitting parameters [default: %(default)s]")
  parser.add_argument("--molfit",    dest="molfit",    type=mu.parray,
                                     action="store",   default=None,
                      help="Molecules fit [default: %(default)s]")
  parser.add_argument("--Tmin",      dest="Tmin",      type=float,
                      action="store",  default=400.0,
                      help="Lower Temperature boundary [default: %(default)s]")
  parser.add_argument("--Tmax",      dest="Tmax",      type=float,
                      action="store",  default=3000.0,
                      help="Higher Temperature boundary [default: %(default)s]")
  parser.add_argument("--quiet",             action="store_true",
                      help="Set verbosity level to minimum",
                      dest="quiet")
  # Input-Converter Options:
  group = parser.add_argument_group("Input Converter Options")
  group.add_argument("--atmospheric_file",  action="store",
                     help="Atmospheric file [default: %(default)s]",
                     dest="atmfile", type=str,    default=None)
  group.add_argument("--PTtype",            action="store",
                     help="PT profile type.",
                     dest="PTtype",  type=str,    default="none")
  group.add_argument("--tint",              action="store",
                     help="Internal temperature of the planet [default: "
                     "%(default)s].",
                     dest="tint",    type=float,  default=100.0)
  group.add_argument("--tint_type", dest="tint_type",
                     help="Method to evaluate `tint`, const or thorngren " +\
                          "[default: %(default)s].",
                     type=str, action="store", default='const')
  # transit Options:
  group = parser.add_argument_group("transit Options")
  group.add_argument("--config",  action="store",
                     help="transit configuration file [default: %(default)s]",
                     dest="config", type=str,    default=None)
  # Output-Converter Options:
  group = parser.add_argument_group("Output Converter Options")
  group.add_argument("--filters",                action="store",
                     help="Waveband filter name [default: %(default)s]",
                     dest="filters",  type=mu.parray, default=None)
  group.add_argument("--tep_name",          action="store",
                     help="A TEP file [default: %(default)s]",
                     dest="tep_name", type=str,    default=None)
  group.add_argument("--kurucz_file",           action="store",
                     help="Stellar Kurucz file [default: %(default)s]",
                     dest="kurucz",   type=str,       default=None)
  group.add_argument("--cloudtop",           action="store",
                     help="Cloud deck top pressure [default: %(default)s]",
                     dest="cloudtop",   type=float, default=None)
  group.add_argument("--scattering",           action="store",
                     help="Rayleigh scattering [default: %(default)s]",
                     dest="scattering", type=float, default=None)
  group.add_argument("--solution",                    action="store",
                     help="Solution geometry [default: %(default)s]",
                     dest="solution", type=str,       default="None",
                     choices=('transit', 'eclipse'))

  parser.set_defaults(**defaults)
  args2, unknown = parser.parse_known_args(remaining_argv)

  # Quiet all threads except rank 0:
  rank = comm.Get_rank()
  verb = rank == 0

  # Get (Broadcast) the number of parameters and iterations from MPI:
  array1 = np.zeros(2, np.int)
  mu.comm_bcast(comm, array1)
  npars, niter = array1

  # :::::::  Initialize the Input converter ::::::::::::::::::::::::::
  atmfile   = args2.atmfile
  molfit    = args2.molfit
  PTtype    = args2.PTtype
  params    = args2.params
  tepfile   = args2.tep_name
  tint      = args2.tint
  tint_type = args2.tint_type
  Tmin      = args2.Tmin
  Tmax      = args2.Tmax
  cloudtop  = args2.cloudtop
  scattering = args2.scattering
  solution  = args2.solution  # Solution type

  # Dictionary of functions to calculate temperature for PTtype
  PTfunc = {'iso'         : pt.PT_iso,
            'line'        : pt.PT_line, 
            'madhu_noinv' : pt.PT_NoInversion,
            'madhu_inv'   : pt.PT_Inversion}

  # Extract necessary values from the TEP file:
  tep = rd.File(tepfile)
  # Stellar temperature in K:
  tstar = float(tep.getvalue('Ts')[0])
  # Stellar radius (in meters):
  rstar = float(tep.getvalue('Rs')[0]) * c.Rsun
  # Semi-major axis (in meters):
  sma   = float(tep.getvalue( 'a')[0]) * sc.au
  # Planetary radius (in meters):
  rplanet = float(tep.getvalue('Rp')[0]) * c.Rjup
  # Planetary mass (in kg):
  mplanet = float(tep.getvalue('Mp')[0]) * c.Mjup

  # Number of fitting parameters:
  nfree   = len(params)                 # Total number of free parameters
  nmolfit = len(molfit)                 # Number of molecular free parameters
  ncloud  = int(cloudtop is not None)
  nray    = int(scattering is not None)
  nradfit = int(solution == 'transit')  # 1 for transit, 0 for eclipse
  nPT     = nfree - nmolfit - ncloud - nray - nradfit   # Number of PT pars

  # Read atmospheric file to get data arrays:
  species, pressure, temp, abundances = mat.readatm(atmfile)
  # Reverse pressure order (for PT to work):
  pressure = pressure[::-1]
  nlayers  = len(pressure)   # Number of atmospheric layers
  nspecies = len(species)    # Number of species in the atmosphere
  mu.msg(verb, "There are {:d} layers and {:d} species.".format(nlayers,
                                                                nspecies))
  # Find index for Hydrogen and Helium:
  species = np.asarray(species)
  iH2     = np.where(species=="H2")[0]
  iHe     = np.where(species=="He")[0]
  # Get H2/He abundance ratio:
  ratio = (abundances[:,iH2] / abundances[:,iHe]).squeeze()
  # Find indices for the metals:
  imetals = np.where((species != "He") & (species != "H2") & \
                     (species != "H-") & (species != 'e-'))[0]
  # Index of molecular abundances being modified:
  imol = np.zeros(nmolfit, dtype='i')
  for i in np.arange(nmolfit):
    imol[i] = np.where(np.asarray(species) == molfit[i])[0]
  
  # Pressure-Temperature profile:
  if PTtype == "line":
    # Planetary surface gravity (in cm s-2):
    gplanet = 100.0 * sc.G * mplanet / rplanet**2
    # Additional PT arguments:
    PTargs  = [rstar, tstar, tint, sma, gplanet, tint_type]
  else:
    PTargs  = None

  # Allocate arrays for receiving and sending data to master:
  freepars = np.zeros(nfree,                 dtype='d')
  profiles = np.zeros((nspecies+1, nlayers), dtype='d')
  # This are sub-sections of profiles, containing just the temperature and
  # the abundance profiles, respectively:
  tprofile  = profiles[0, :]
  aprofiles = profiles[1:,:]

  # Store abundance profiles:
  for i in np.arange(nspecies):
    aprofiles[i] = abundances[:, i]

  # :::::::  Spawn transit code  :::::::::::::::::::::::::::::::::::::
  # # transit configuration file:
  transitcfile = args2.tconfig
 
  # Initialize the transit python module:
  transit_args = ["transit", "-c", transitcfile]
  trm.transit_init(len(transit_args), transit_args)

  # Get wavenumber array from transit:
  nwave  = trm.get_no_samples()
  specwn = trm.get_waveno_arr(nwave)

  # :::::::  Output Converter  :::::::::::::::::::::::::::::::::::::::
  ffile    = args2.filters    # Filter files
  kurucz   = args2.kurucz     # Kurucz file

  # Log10(stellar gravity)
  gstar = float(tep.getvalue('loggstar')[0])
  # Planet-to-star radius ratio:
  rprs  = rplanet / rstar

  nfilters = len(ffile)  # Number of filters:

  # FINDME: Separate filter/stellar interpolation?
  # Get stellar model:
  starfl, starwn, tmodel, gmodel = w.readkurucz(kurucz, tstar, gstar)
  # Read and resample the filters:
  nifilter  = [] # Normalized interpolated filter
  istarfl   = [] # interpolated stellar flux
  wnindices = [] # wavenumber indices used in interpolation
  for i in np.arange(nfilters):
    # Read filter:
    filtwaven, filttransm = w.readfilter(ffile[i])
    # Check that filter boundaries lie within the spectrum wn range:
    if filtwaven[0] < specwn[0] or filtwaven[-1] > specwn[-1]:
      mu.exit(message="Wavenumber array ({:.2f} - {:.2f} cm-1) does not "
              "cover the filter[{:d}] wavenumber range ({:.2f} - {:.2f} "
              "cm-1).".format(specwn[0], specwn[-1], i, filtwaven[0],
                                                        filtwaven[-1]))

    # Resample filter and stellar spectrum:
    nifilt, strfl, wnind = w.resample(specwn, filtwaven, filttransm,
                                              starwn,    starfl)
    nifilter.append(nifilt)
    istarfl.append(strfl)
    wnindices.append(wnind)

  # Allocate arrays for receiving and sending data to master:
  spectrum = np.zeros(nwave,    dtype='d')
  bandflux = np.zeros(nfilters, dtype='d')

  # Allocate array to receive parameters from MPI:
  params = np.zeros(npars, np.double)

  # ::::::  Main MCMC Loop  ::::::::::::::::::::::::::::::::::::::::::
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  while niter >= 0:
    niter -= 1
    # Receive parameters from MCMC:
    mu.comm_scatter(comm, params)

    # Check for the MCMC-end flag:
    if params[0] == np.inf:
      break

    # Input converter calculate the profiles:
    try:
      tprofile[:] = pt.PT_generator(pressure,       params[0:nPT], 
                                    PTfunc[PTtype], PTargs       )[::-1]
    except ValueError:
      mu.msg(verb, 'Input parameters give non-physical profile.')
      # FINDME: what to do here?

    # If the temperature goes out of bounds:
    if np.any(tprofile < Tmin) or np.any(tprofile > Tmax):
      mu.comm_gather(comm, -np.ones(nfilters), MPI.DOUBLE)
      continue
    # Scale abundance profiles:
    for i in np.arange(nmolfit):
      m = imol[i]
      # Use variable as the log10:
      aprofiles[m] = abundances[:, m] * 10.0**params[nPT+nradfit+ncloud+nray+i]

    # Update H2, He abundances so sum(abundances) = 1.0 in each layer:
    q = 1.0 - np.sum(aprofiles[imetals], axis=0)
    if np.any(q < 0.0):
      mu.comm_gather(comm, -np.ones(nfilters), MPI.DOUBLE)
      #mu.comm_gather(comm, -np.ones(nwave), MPI.DOUBLE)
      continue
    aprofiles[iH2] = ratio * q / (1.0 + ratio)
    aprofiles[iHe] =         q / (1.0 + ratio)

    # Set the 'surface' level:
    if solution == "transit":
      trm.set_radius(params[nPT])

    if ncloud == 1:
      trm.set_cloudtop(params[nPT+nradfit])

    if nray == 1:
      trm.set_scattering(params[nPT+nradfit+ncloud])

    # Let transit calculate the model spectrum:
    spectrum = trm.run_transit(profiles.flatten(), nwave)

    # Calculate the band-integrated intensity per filter:
    for i in np.arange(nfilters):
      if   solution == "eclipse":
        fluxrat = (spectrum[wnindices[i]]/istarfl[i]) * rprs*rprs
        bandflux[i] = w.bandintegrate(fluxrat, specwn,
                                      nifilter[i], wnindices[i])
      elif solution == "transit":
        bandflux[i] = w.bandintegrate(spectrum[wnindices[i]], specwn,
                                      nifilter[i], wnindices[i])

    # Send resutls back to MCMC:
    mu.comm_gather(comm, bandflux, MPI.DOUBLE)

  # ::::::  End main Loop  :::::::::::::::::::::::::::::::::::::::::::
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # Close communications and disconnect:
  mu.comm_disconnect(comm)
  trm.free_memory()


if __name__ == "__main__":
  # Open communications with the master:
  comm = MPI.Comm.Get_parent()
  main(comm)
