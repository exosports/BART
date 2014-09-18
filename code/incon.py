#! /usr/bin/env python
import sys, os
import argparse, ConfigParser
import numpy as np
from mpi4py import MPI

# Input converted dir:
ICdir = os.path.dirname(os.path.realpath(__file__))

import makeatm as mat
import PT as pt
sys.path.append(ICdir + "/../modules/MCcubed/src/")
import mcutils as mu


def main(comm):
  """
  Modification History:
  ---------------------
  2014-03-24  Madison   First implementation completed.
                        Madison Stemm, UCF.  user@email.com
  2014-04-13  patricio  Adapted for BART project.
                        Patricio Cubillos, UCF.  pcubillos@fulbrightmail.org
  2014-04-24  patricio  
  """
  # Initialize parser to process a configuration file:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--config_file",
                       help="Specify config file", metavar="FILE")
  # Remaining_argv contains all other command-line-arguments:
  args, remaining_argv = cparser.parse_known_args()

  # Get parameters from configuration file (if exists):
  cfile = args.config_file # The configuration file
  if cfile:
    config = ConfigParser.SafeConfigParser()
    config.read([cfile])
    defaults = dict(config.items("MCMC"))
  else:
    defaults = {}
  # Now, parser for the Input-Converter arguments:
  parser = argparse.ArgumentParser(parents=[cparser])

  # IC Options:
  group = parser.add_argument_group("Input Converter Options")
  group.add_argument("-a", "--atmospheric_file",  action="store",
                     help="Atmospheric file [default: %(default)s]",
                     dest="atmfile", type=str,    default=None)
  group.add_argument("-t", "--tep_file",          action="store",
                     help="A TEP file [default: %(default)s]",
                     dest="tepfile", type=str,    default=None)
  group.add_argument("-q", "--quiet",             action="store_true",
                     help="Set verbosity level to minimum",
                     dest="quiet")
  # MCMC Options:
  group = parser.add_argument_group("MCMC Options")
  group.add_argument("-f", "--params",         action="store",
                     help="Model fitting parameter [default: %(default)s]",
                     dest="params", type=mu.parray,   default=None)
  group.add_argument("-m", "--molfit",         action="store",
                     help="Molecule fit [default: %(default)s]",
                     dest="molfit", type=mu.parray,   default=None)

  # Set the defaults from the configuration file:
  parser.set_defaults(**defaults)
  # Set values from command line:
  args2 = parser.parse_args(remaining_argv)
  rank = comm.Get_rank()

  # Now start working:
  atmfile = args2.atmfile
  tepfile = args2.tepfile
  molfit  = args2.molfit
  params = args2.params

  # Number of PT parameters:
  nfree   = len(params)
  nmolfit = len(molfit)
  nPT     = len(params) - len(molfit)

  # Verbosity level:
  if args2.quiet:
    verb = 0
  else:
    verb = 1

  mu.msg(verb, "ICON FLAG 50")
  # Read atmospheric file to get data arrays:
  species, pressure, abundances = mat.readatm(atmfile)
  nlayers  = len(pressure)
  nspecies = len(species)
  mu.msg(verb, "There are {:d} layers and {:d} species.".format(nlayers,
                                                                nspecies))
  # Send nlayers + nspecies to master:
  mu.comm_gather(comm, np.array([nlayers, nspecies], dtype='i'), MPI.INT)
  mu.msg(verb, "ICON FLAG 55")

  # Get the number of PT params and iterations:
  array1 = np.zeros(1, dtype='i')
  mu.comm_bcast(comm, array1)
  niter = array1
  mu.msg(verb, "ICON FLAG 60: nPT=%d,  niter=%d"%(nPT, niter))

  # Determine inversion or non-inversion (by looking at parameter P2):
  if params[3] == -1.0:
    MadhuPT = False
  else:
    MadhuPT = True

  # Allocate arrays for receiving and sending data to master:
  freepars = np.zeros(nfree,                dtype='d')
  profiles = np.zeros(nlayers*(nspecies+1), dtype='d')

  # Index of molecular abundances being modified:
  imol = np.zeros(nmolfit, dtype='i')
  for i in np.arange(nmolfit):
    imol[i] = np.where(np.asarray(species) == molfit[i])[0]
  mu.msg(verb, "ICON FLAG 62: i-molecs {}".format(imol))

  # Store abundance profiles:
  for i in np.arange(nspecies):
    profiles[(i+1)*nlayers:(i+2)*nlayers] = abundances[:, i]

  # ::::::  MAIN MCMC LOOP  ::::::::::::::::::::::::::::::::::::::::::
  mu.msg(verb, "ICON FLAG 70: Enter main loop")
  while niter >= 0:
    # Receive parameters from master:
    mu.comm_scatter(comm, freepars)
    mu.msg(verb, "ICON FLAG 72: incon pars: {}".format(freepars))

    # Get temperature profile and check for non-physical outputs:
    try:
      profiles[0:nlayers] = pt.PT_generator(pressure, freepars[0:nPT], MadhuPT)
    except ValueError:
      mu.msg(verb, 'Input parameters give non-physical profile.')
      # FINDME: what to do here?

    # Scale abundance profiles:
    for i in np.arange(nmolfit):
      m = imol[i]
      profiles[(m+1)*nlayers:(m+2)*nlayers] = abundances[:,m]*freepars[nPT+i]
    
    mu.msg(verb, "ICON FLAG 76: profiles size: {:d}".format(len(profiles)))
    # Gather (send) the temperature and abundance profile:
    mu.comm_gather(comm, profiles, MPI.DOUBLE)
    niter -= 1

  # Close communications and disconnect:
  if rank == 0:
    mu.msg(verb, "ICON FLAG 99: Input Con is out")
  mu.exit(comm)

if __name__ == "__main__":
  comm = MPI.Comm.Get_parent()
  main(comm)
