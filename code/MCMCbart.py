import os, sys
import numpy as np
from mpi4py import MPI

BARTdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(BARTdir + "/../modules/MCcubed/src/")
import mcutils as mu

def init(incomm, transitcomm, outcomm, nparams, niter, verb=0):
  """
  BART initialization routine (everything that needs to be done before
  entering the MCMC loop).

  Parameters:
  -----------
  incomm: MPI communicator
  transitcomm: MPI communicator
  outcomm: MPI communicator
  nparams: Integer
  niter: Integer

  Modification History:
  ---------------------
  2014-06-25  patricio  Initial implementation.  pcubillos@fulbrightmail.org
  2014-07-24  patricio  Added documentation.
  """
  # Number of sub-functions that comprise the BART code:
  nfuncs = 3

  # Number of PT params:
  nPT = 6

  # Get the output-array sizes:
  arrsize = np.zeros(2, dtype="i")
  # Get the number of layers from incomm:
  mu.comm_gather(incomm,      arrsize)
  nlayers  = arrsize[0]
  nspecies = arrsize[1]

  # Get the number of spectral samples from transit:
  arrsize = np.zeros(1, dtype="i")
  mu.comm_gather(transitcomm, arrsize)
  nwave   = arrsize[0]
  # Get the number of data points from outcomm:
  mu.comm_gather(outcomm,     arrsize)
  nbands  = arrsize[0]

  # Transits receives the temperature and abundance profiles:
  nprofile = nlayers * (nspecies + 1)

  mu.msg(verb, "BART FLAG 50:")
  # Send the input-array sizes to incomm, transit, and outcomm:
  mu.comm_bcast(incomm,      np.array([          niter], dtype="i"), MPI.INT)
  mu.comm_bcast(transitcomm, np.array([nprofile, niter], dtype="i"), MPI.INT)
  mu.comm_bcast(outcomm,     np.array([nwave,    niter], dtype="i"), MPI.INT)

  mu.msg(verb, "BART FLAG 55:")
  # Get wavenumber array from transit and sent to outcomm:
  specwn = np.zeros(nwave, dtype="d")
  mu.comm_gather(transitcomm, specwn)
  mu.msg(verb, "BART FLAG 57:")
  mu.comm_scatter(outcomm,    specwn, MPI.DOUBLE)

  # Allocate arrays received form each communicator:
  profiles = np.empty(nprofile, np.double)
  spectrum = np.empty(nwave,    np.double)
  bandflux = np.empty(nbands,   np.double)

  # Return:
  return [profiles, spectrum, bandflux]
  

def bart(params, incomm, transitcomm, outcomm, profiles, spectrum, bandflux,
         verb=1):
  """
  Function to loop between input - transit - and output

  Parameters:
  -----------
  params:  1D float ndarray
     Array of fitting parameters.
  incomm:  MPI communicator
  transitcomm:  MPI communicator
  outcomm:  MPI communicator
  profiles:  1D float ndarray
     Array of data returned by incom.
  spectrum:  1D float ndarray
     Array of data returned by transit.
  bandflux:  1D float ndarray
     Array of data returned by outcom.

  Return:
  -------
  bandflux: 1D ndarray

  Modification History:
  ---------------------
  2014-05-28  patricio  Initial implementation
  2014-09-05  patricio  Get returned arrays as arguments instead of allocating
                        each time.
  """
  mu.msg(verb, "BART FLAG 73: Start iteration")
  # Scatter (send) free parameter to IC:
  mu.msg(verb, "BART FLAG 74: free pars: {}".format(params))
  mu.comm_scatter(incomm, params, MPI.DOUBLE)
  # Gather (receive) the profiles:
  mu.msg(verb, "BART FLAG 76: Tprofile size: %d"%len(profiles))
  mu.comm_gather(incomm, profiles)

  mu.msg(verb, "BART FLAG 77: profiles={}".format(profiles[0:10]))
  mu.msg(verb, "BART FLAG 78: intransit size: {:d}.".format(len(profiles)))
  # Scatter (send) Tprofile and abundance factors to transit:
  mu.comm_scatter(transitcomm, profiles, MPI.DOUBLE)
  mu.msg(verb, "BART FLAG 79: sent data to transit!")
  # Gather (receive) spectrum from transit:
  mu.comm_gather(transitcomm, spectrum)

  mu.msg(verb, "BART FLAG 80: Send spectrum to OC.")
  # Scatter (send) spectrum to output converter:
  mu.comm_scatter(outcomm, spectrum, MPI.DOUBLE)
  # Gather (receive) band-integrated flux:
  mu.comm_gather(outcomm, bandflux)

  mu.msg(verb, "BART FLAG 85: OC out: %.6e, %.6e"%(bandflux[0], bandflux[1]))
  return bandflux


def main():
  """

  Returns:
  --------
  submodules: List
     List with the names of the submodules to run.
  bart: callable
     Function that initializes the submodules.

  Modification History:
  ---------------------
  2014-08-19  patricio  Added documentation.
  """
  return ["incon.py",
          BARTdir + "/../modules/transit/transit/MPItransit",
          "outcon.py"], bart

if __name__ == "__main__":
  # Open communications with the master:
  main()
