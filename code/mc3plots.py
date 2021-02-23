# Copyright (c) 2015-2018 Patricio Cubillos and contributors.
# MC3 is open-source software under the MIT license (see LICENSE).

import sys, os
import six
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                '..', 'modules', 'MCcubed', 'MCcubed', 'plots', ''))
from mcplots import trace, pairwise, histogram

__all__ = ["mc3plots"]


def mc3plots(output,   burnin,   thinning, nchains, uniform, molfit, 
             out_spec, parnames, stepsize, date_dir, fnames, fext='.png'):
  """
  Reformats the MC3 output file so that the log(abundance) factor is with 
  respect to molar fraction, rather than the initial values (as MC3 does). 
  Calls trace(), pairwise(), and histogram() using these values.

  Parameters
  ----------
  output  : string. Path to MC3 output.npy file.
  burnin  : int. Number of burn-in iterations.
  thinning: int. Thinning factor of the chains (use every thinning-th 
                 iteration) used for plotting.
  uniform : array-like. If not None, set uniform abundances with the 
                        specified values for each species.
  nchains : int. Number of parallel chains in MCMC.
  molfit  : list, strings. Molecules to be fit by the MCMC.
  out_spec: list, strings. Molecules included in atmospheric file.
  parnames: list, strings. Parameter names.
  stepsize: array, floats.  Initial stepsizes of MCMC parameters.
  date_dir: string. Path to directory where plots are to be saved.
  fnames  : list, strings. File names for the trace, pairwise, and histogram 
                           plots, in that order.
  fext    : string.        File extension for the plots to be saved.
                           Options: .png, .pdf
                           Default: .png
  """
  # Load and stack results, excluding burn-in
  allparams = np.load(date_dir + output)
  allstack  = allparams[0, :, burnin:]
  for c in np.arange(1, allparams.shape[0]):
    allstack = np.hstack((allstack, allparams[c, :, burnin:]))

  # Subtract initial abundances if uniform, so that plots are log(abundance)
  if uniform is not None:
    molind = []
    for imol in range(len(molfit)):
      for j  in range(len(out_spec.split())):
        if molfit[imol]+'_' in out_spec.split()[j] and \
           stepsize[-len(molfit):][imol] > 0:
          molind.append(j)
    allstack[-len(molfit):, :] += \
                               np.log10(uniform[molind]).reshape(len(molind),1)

  # Slice only params that are varied (remove static params)
  ipar = stepsize != 0
  if parnames is None:
    nparams  = stepsize.size
    namelen  = int(2+np.log10(np.amax([nparams-1,1])))
    parnames = np.zeros(nparams, "|S%d"%namelen if six.PY2 else "<U%d"%namelen)
    for i in np.arange(nparams):
      parnames[i] = "P" + str(i).zfill(namelen-1)
  parnames = [parnames[i] for i in range(len(parnames)) if ipar[i]]

  # Trace plot:
  trace(    allstack, parname=parnames, thinning=thinning,
            savefile=date_dir + fnames[0] + fext,
            sep=np.size(allstack[0])/nchains)
  # Pairwise posteriors:
  pairwise( allstack, parname=parnames, thinning=thinning,
            savefile=date_dir + fnames[1] + fext)
  # Histograms:
  histogram(allstack, parname=parnames, thinning=thinning,
            savefile=date_dir + fnames[2] + fext)


