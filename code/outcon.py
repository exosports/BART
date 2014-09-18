#! /usr/bin/env python
import sys, os
import argparse, ConfigParser
import numpy as np
from mpi4py import MPI

# Output converter dir:
OCdir = os.path.dirname(os.path.realpath(__file__))

import wine   as w
import reader as rd
sys.path.append(OCdir + "/../modules/MCcubed/src/")
import mcutils as mu


def main(comm):
  """
  Calculate the waveband-integrated values for a model spectrum.

  Developers:
  -----------
  Madison Stemm      UCF  FINDME
  Patricio Cubillos  UCF  pcubillos@fulbrightmail.org

  Modification History:
  ---------------------
  2014-03-24  Madison   First implementation completed.
  2014-04-13  patricio  Adapted for BART project.
  """
  # Initialise parser to process a configuration file:
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

  # Now, parser for the Output-Converter arguments:
  parser = argparse.ArgumentParser(parents=[cparser])
  
  # Output-Converter Options:
  group = parser.add_argument_group("Output Converter Options")
  group.add_argument(     "--filter",                 action="store",
                     help="Waveband filter name [default: %(default)s]",
                     dest="filter",   type=mu.parray, default=None)
  group.add_argument("--tep_name",         action="store",
                     help="Transiting Exo-Planet file [default: %(default)s]",
                     dest="tep_name", type=str,       default=None)
  group.add_argument("-k", "--kurucz_file",           action="store",
                     help="Stellar Kurucz file [default: %(default)s]",
                     dest="kurucz",   type=str,       default=None)
  group.add_argument("-q", "--quiet",                 action="store_true",
                     help="Set verbosity level to minimum",
                     dest="quiet")

  # Set the defaults from the configuration file:
  parser.set_defaults(**defaults)
  # Set values from command line:
  args2 = parser.parse_args(remaining_argv)

  rank  = comm.Get_rank()
  size  = comm.Get_size()

  # Retrieve options and arguments:
  ffile   = args2.filter    # Filter files
  kurucz  = args2.kurucz    # Kurucz file
  tepfile = args2.tep_name  # TEP file

  # Verbosity level:
  if args2.quiet:
    verb = 0
  else:
    verb = 1

  # Some constants:
  # http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
  # http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
  Rjup =  71492.0 * 1e3 # m
  Rsun = 696000.0 * 1e3 # m

  # Extract necessary values from the TEP file:
  tep = rd.File(tepfile)
  # Stellar temperature in K.
  tstar = float(tep.getvalue('Ts')[0])
  # Log10(stellar gravity)
  gstar = float(tep.getvalue('loggstar')[0])
  # Planet-to-star radius ratio:
  rprs  = float(tep.getvalue('Rp')[0])*Rjup/(float(tep.getvalue('Rs')[0])*Rsun)
  mu.msg(verb, "OCON FLAG 10: {}, {}, {}".format(tstar, gstar, rprs))

  # Send the number of filters:
  nfilters = len(ffile)
  mu.comm_gather(comm, np.array([nfilters], dtype='i'), MPI.INT)
  mu.msg(verb,"OCON FLAG 55: Start")

  # Get the number of iterations and wavenumbers:
  array1 = np.zeros(2, dtype='i')
  mu.comm_bcast(comm, array1)
  nspec, niter = array1
  mu.msg(verb, "OCON FLAG 65: nspec={:d}, niter={:d}".format(nspec, niter))

  # Get the wavenumber array:
  specwn = np.zeros(nspec, dtype="d")
  mu.comm_scatter(comm, specwn)
  mu.msg(verb, "OCON FLAG 65.5: wn=[{:.2f}, {:.2f}, ..., {:.2f}]".format(
              specwn[0], specwn[1], specwn[-1]))

  # Get stellar model:
  starfl, starwn, tmodel, gmodel = w.readkurucz(kurucz, tstar, gstar)
  # Read and resample the filters:
  nifilter  = [] # Normalized interpolated filter
  istarfl   = [] # interpolated stellar flux
  wnindices = [] # wavenumber indices used in interpolation
  mu.msg(verb, "OCON FLAG 66: Prepare!")
  for i in np.arange(nfilters):
    # Read filter:
    filtwaven, filttransm = w.readfilter(ffile[i])
    # Check that filter boundaries lie within the spectrum wn range:
    if filtwaven[0] < specwn[0] or filtwaven[-1] > specwn[-1]:
      mu.exit(message="Wavenumber array ({:.2f} - {:.2f} cm-1) does not cover "
            "the filter[{:d}] wavenumber range ({:.2f} - {:.2f} cm-1).".format(
            specwn[0], specwn[-1], i, filtwaven[0], filtwaven[-1]))

    # Resample filter and stellar spectrum:
    nifilt, strfl, wnind = w.resample(specwn, filtwaven, filttransm,
                                              starwn,    starfl)
    mu.msg(verb, "OCON FLAG 67: mean star flux: %.3e"%np.mean(strfl))
    nifilter.append(nifilt)
    istarfl.append(strfl)
    wnindices.append(wnind)

  # Allocate arrays for receiving and sending data to master:
  array3   = np.zeros(nspec,    dtype='d')
  array4   = np.zeros(nfilters, dtype='d')

  # Worker loop, communicating with Master
  mu.msg(verb, "OCON FLAG 70: enter loop")
  while niter >= 0:
    # Receive planet spectrum:
    mu.comm_scatter(comm, array3)
    mu.msg(verb, "OCON FLAG 72: receive spectum")

    # Calculate the band-integrated intensity per filter:
    for i in np.arange(nfilters):
      fluxrat = (array3[wnindices[i]]/istarfl[i])*rprs*rprs
      array4[i] = w.bandintegrate(fluxrat, specwn, nifilter[i], wnindices[i])
    #print("gather send: " + str(array4))

    # Gather (send) the band-integrated fluxes:
    mu.comm_gather(comm, array4, MPI.DOUBLE)
    niter -= 1

  # Close communications and disconnect:
  if rank == 0:
    mu.msg(verb, "OCON FLAG 99: OutputCon is out.")
  mu.exit(comm)


if __name__ == "__main__":
  # Open communications with the master:
  comm = MPI.Comm.Get_parent()
  main(comm)
