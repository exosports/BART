# This is a modified version of mcplots.py from MC3 
# for improved plotting in BART.
# Copyright (c) 2015-2018 Patricio Cubillos and contributors.
# MC3 is open-source software under the MIT license (see LICENSE).

"""
    This code makes trace, pairwise, and histogram plots of the MCMC results.

    Functions
    ---------
    mcplots: 
        Makes trace, pairwise, and histogram plots as MC3 does, but 
        with features unique to BART (plots with respect to the molar 
        fraction instead of initial abundances (if uniform), and 
        reformats parameter names in plots)
    trace: 
        Plot parameter trace MCMC sampling
    pairwise: 
        Plot parameter pairwise posterior distributions
    histogram: 
        Plot parameter marginal posterior distributions
    RMS: 
        Plot the RMS vs binsize
    modelfit: 
        Plot the model and (binned) data arrays, and their residuals
    reformatparname:
        Reformats special-case parameter names for plotting purposes

    Revisions
    ---------
    2018-12-21  Michael  Adapted from MC3 code. Added mcplots() and 
                         reformatparname() functions. Added/updated 
                         documentation.
"""

import sys, os
import numpy as np
import matplotlib as mpl
#mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/../modules/MCcubed/MCcubed/lib")
import binarray as ba

__all__ = ["mcplots", "trace", "pairwise", "histogram", "RMS", "modelfit"]

def mcplots(output,   burnin,   thinning, nchains, uniform, molfit, 
            out_spec, parnames, stepsize, date_dir, fnames):
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
  """
  # Load and stack results, excluding burn-in
  allparams = np.load(date_dir + output)
  allstack  = allparams[0, :, burnin:]
  for c in np.arange(1, allparams.shape[0]):
    allstack = np.hstack((allstack, allparams[c, :, burnin:]))

  # Subtract initial abundances if uniform, so that plots are log(abundance)
  if uniform is not None and np.all(stepsize > 0):
    molind = []
    for imol in range(len(molfit)):
      for j  in range(len(out_spec.split())):
        if molfit[imol]+'_' in out_spec.split()[j] and \
           stepsize[-len(molfit):][imol] > 0:
          molind.append(j)
    allstack[-len(molfit):, :] += \
                               np.log10(uniform[molind]).reshape(len(molind),1)

  # Slice only params that are varied (remove static params)
  ipar     = stepsize != 0
  # Note 'parnames' is a list, so cannot index using an array/list
  parnames = [parnames[i] for i in range(len(parnames)) if ipar[i]]

  # Trace plot:
  trace(    allstack, parname=parnames, thinning=thinning,
            savefile=date_dir + fnames[0],
            sep=np.size(allstack[0])/nchains)
  # Pairwise posteriors:
  pairwise( allstack, parname=parnames, thinning=thinning,
            savefile=date_dir + fnames[1])
  # Histograms:
  histogram(allstack, parname=parnames, thinning=thinning,
            savefile=date_dir + fnames[2])


def trace(allparams, title=None, parname=None, thinning=1,
          fignum=-10, savefile=None, fmt=".", sep=None):
  """
  Plot parameter trace MCMC sampling

  Parameters
  ----------
  allparams: 2D ndarray
     An MCMC sampling array with dimension (number of parameters,
     sampling length).
  title: String
     Plot title.
  parname: Iterable (strings)
     List of label names for parameters.  If None use ['P0', 'P1', ...].
  thinning: Integer
     Thinning factor for plotting (plot every thinning-th value).
  fignum: Integer
     The figure number.
  savefile: Boolean
     If not None, name of file to save the plot.
  fmt: String
     The format string for the line and marker.
  sep: Integer
     Number of samples per chain. If not None, draw a vertical line
     to mark the separation between the chains.

  Uncredited developers
  ---------------------
  Kevin Stevenson (UCF)
  """
  # Get number of parameters and length of chain:
  npars, niter = np.shape(allparams)
  fs = 10

  # Set default parameter names:
  if parname is None:
    namelen = int(2+np.log10(np.amax([npars-1,1])))
    parname = np.zeros(npars, "|S%d"%namelen)
    for i in np.arange(npars):
      parname[i] = "P" + str(i).zfill(namelen-1)

  # Reformat special-case parameter names
  reformatpar = reformatparname(parname)

  # Get location for chains separations:
  xmax = len(allparams[0,0::thinning])
  if sep is not None:
    xsep = np.arange(sep/thinning, xmax, sep/thinning)

  # Make the trace plot:
  plt.figure(fignum, figsize=(8,8))
  plt.clf()
  if title is not None:
    plt.suptitle(title, size=16)

  plt.subplots_adjust(left=0.15, right=0.95, bottom=0.10, top=0.90,
                      hspace=0.3)

  for i in np.arange(npars):
    a = plt.subplot(npars, 1, i+1)
    plt.plot(allparams[i, 0::thinning], fmt)
    yran = a.get_ylim()
    if sep is not None:
      plt.vlines(xsep, yran[0], yran[1], "0.3")
    plt.xlim(0, xmax)
    plt.ylim(yran)
    plt.ylabel(reformatpar[i], size=fs+4, multialignment='center')
    plt.yticks(size=fs)
    # Make sure ticks are read-able
    if   len(a.get_yticks()[::2]) > 5:
      a.set_yticks(a.get_yticks()[::3])
    elif len(a.get_yticks()[::2]) > 2:
      a.set_yticks(a.get_yticks()[::2])
    if i == npars - 1:
      plt.xticks(size=fs)
      if thinning > 1:
        plt.xlabel('MCMC (thinned) iteration', size=fs+4)
      else:
        plt.xlabel('MCMC iteration', size=fs+4)
    else:
      plt.xticks(visible=False)
    # Align labels
    a.yaxis.set_label_coords(-0.1, 0.5)

  if savefile is not None:
    plt.savefig(savefile, bbox_inches='tight')


def pairwise(allparams, title=None, parname=None, thinning=1,
             fignum=-11, savefile=None, style="hist"):
  """
  Plot parameter pairwise posterior distributions

  Parameters
  ----------
  allparams: 2D ndarray
     An MCMC sampling array with dimension (number of parameters,
     sampling length).
  title: String
     Plot title.
  parname: Iterable (strings)
     List of label names for parameters.  If None use ['P0', 'P1', ...].
  thinning: Integer
     Thinning factor for plotting (plot every thinning-th value).
  fignum: Integer
     The figure number.
  savefile: Boolean
     If not None, name of file to save the plot.
  style: String
     Choose between 'hist' to plot as histogram, or 'points' to plot
     the individual points.

  Uncredited developers
  ---------------------
  Kevin Stevenson  (UCF)
  Ryan Hardy  (UCF)
  """
  # Get number of parameters and length of chain:
  npars, niter = np.shape(allparams)

  # Don't plot if there are no pairs:
  if npars == 1:
    return

  # Set default parameter names:
  if parname is None:
    namelen = int(2+np.log10(np.amax([npars-1,1])))
    parname = np.zeros(npars, "|S%d"%namelen)
    for i in np.arange(npars):
      parname[i] = "P" + str(i).zfill(namelen-1)
  fs = 10

  # Reformat parameter names for special cases when plotting
  reformatpar = reformatparname(parname)

  # Set palette color:
  palette = plt.matplotlib.colors.LinearSegmentedColormap('YlOrRd2',
                                               plt.cm.datad['YlOrRd'], 256)
  palette.set_under(alpha=0.0)
  palette.set_bad(alpha=0.0)

  fig = plt.figure(fignum, figsize=(8,8))
  plt.clf()
  if title is not None:
    plt.suptitle(title, size=16)

  h = 1 # Subplot index
  plt.subplots_adjust(left  =0.15, right =0.85, bottom=0.15, top=0.85,
                      hspace=0.3,  wspace=0.3)

  for   j in np.arange(npars): # Rows
    for i in np.arange(npars):  # Columns
      if j > i or j == i:
        a = plt.subplot(npars, npars, h)
        # Y labels:
        if i == 0 and j != 0:
          plt.yticks(size=fs)
          plt.ylabel(reformatpar[j], size=fs+4, multialignment='center')
        elif i == 0 and j == 0:
          plt.yticks(visible=False)
          plt.ylabel(reformatpar[j], size=fs+4, multialignment='center')
        else:
          a = plt.yticks(visible=False)
        # X labels:
        if j == npars-1:
          plt.xticks(size=fs, rotation=90)
          plt.xlabel(reformatpar[i], size=fs+4)
        else:
          a = plt.xticks(visible=False)
        # The plot:
        if style=="hist":
          if j > i:
            hist2d, xedges, yedges = np.histogram2d(allparams[i, 0::thinning],
                                                    allparams[j, 0::thinning], 
                                                    20, normed=False)
            vmin = 0.0
            hist2d[np.where(hist2d == 0)] = np.nan
            a = plt.imshow(hist2d.T, extent=(xedges[0], xedges[-1], yedges[0],
                           yedges[-1]), cmap=palette, vmin=vmin, aspect='auto',
                           origin='lower', interpolation='bilinear')
          else:
            a = plt.hist(allparams[i,0::thinning], 20, normed=False)
          a = plt.gca()

        elif style=="points":
          if j > i:
            a = plt.plot(allparams[i], allparams[j], ",")
          else:
            a = plt.hist(allparams[i,0::thinning], 20, normed=False)
          a = plt.gca()
        # Make sure ticks are read-able
        if   len(a.get_xticks()[::2]) > 4:
          a.set_xticks(a.get_xticks()[::3])
        elif len(a.get_xticks()[::2]) > 2:
          a.set_xticks(a.get_xticks()[::2])
        if   len(a.get_yticks()[::2]) > 4:
          a.set_yticks(a.get_yticks()[::3])
        elif len(a.get_yticks()[::2]) > 2:
          a.set_yticks(a.get_yticks()[::2])
        # Align labels
        if j == npars-1 and i == npars-1:
          axs = fig.get_axes()
          for ax in axs:
            ss = ax.get_subplotspec()
            nrows, ncols, start, stop = ss.get_geometry()
            if start//nrows == nrows-1:
              ax.xaxis.set_label_coords(0.5, -0.155 * npars)
            if start%ncols == 0:
              ax.yaxis.set_label_coords(-0.155 * npars, 0.5)
    
      h += 1
  # The colorbar:
  if style == "hist":
    if npars > 2:
      a = plt.subplot(2, 6, 5, frameon=False)
      a.yaxis.set_visible(False)
      a.xaxis.set_visible(False)
    bounds = np.linspace(0, 1.0, 64)
    norm = mpl.colors.BoundaryNorm(bounds, palette.N)
    ax2 = fig.add_axes([0.85, 0.535, 0.025, 0.36])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=palette, norm=norm,
          spacing='proportional', boundaries=bounds, format='%.1f')
    cb.set_label("Normalized Point Density", fontsize=fs)
    cb.set_ticks(np.linspace(0, 1, 5))
    plt.draw()

  # Save file:
  if savefile is not None:
    plt.savefig(savefile, bbox_inches='tight')


def histogram(allparams, title=None, parname=None, thinning=1,
              fignum=-12, savefile=None):
  """
  Plot parameter marginal posterior distributions

  Parameters
  ----------
  allparams: 2D ndarray
     An MCMC sampling array with dimension (number of parameters,
     sampling length).
  title: String
     Plot title.
  parname: Iterable (strings)
     List of label names for parameters.  If None use ['P0', 'P1', ...].
  thinning: Integer
     Thinning factor for plotting (plot every thinning-th value).
  fignum: Integer
     The figure number.
  savefile: Boolean
     If not None, name of file to save the plot.

  Uncredited developers
  ---------------------
  Kevin Stevenson  (UCF)
  """
  # Get number of parameters and length of chain:
  npars, niter = np.shape(allparams)
  fs = 14  # Fontsize

  # Set default parameter names:
  if parname is None:
    namelen = int(2+np.log10(np.amax([npars-1,1])))
    parname = np.zeros(npars, "|S%d"%namelen)
    for i in np.arange(npars):
      parname[i] = "P" + str(i).zfill(namelen-1)

  # Reformat special cases of parameter names
  reformatpar = reformatparname(parname)

  # Set number of rows:
  if npars < 10:
    nrows = (npars - 1)/3 + 1
  else:
    nrows = (npars - 1)/4 + 1
  # Set number of columns:
  if   npars > 9:
    ncolumns = 4
  elif npars > 4:
    ncolumns = 3
  else:
    ncolumns = (npars+2)/3 + (npars+2)%3  # (Trust me!)

  histheight = np.amin((2 + 2*(nrows), 8))
  if nrows == 1:
    bottom = 0.25
  else:
    bottom = 0.15

  plt.figure(fignum, figsize=(8, histheight))
  plt.clf()
  plt.subplots_adjust(left=0.1, right=0.95, bottom=bottom, top=0.9,
                      hspace=0.1, wspace=0.1)

  if title is not None:
    a = plt.suptitle(title, size=16)

  maxylim = 0  # Max Y limit
  for i in np.arange(npars):
    ax = plt.subplot(nrows, ncolumns, i+1)
    a  = plt.xticks(size=fs, rotation=90)
    if i%ncolumns == 0:
      a = plt.yticks(size=fs)
    else:
      a = plt.yticks(visible=False)
    plt.xlabel(reformatpar[i], size=fs)
    a = plt.hist(allparams[i,0::thinning], 20, normed=False)
    maxylim = np.amax((maxylim, ax.get_ylim()[1]))

  # Set uniform height:
  for i in np.arange(npars):
    ax = plt.subplot(nrows, ncolumns, i+1)
    ax.set_ylim(0, maxylim)

  if savefile is not None:
    plt.tight_layout()
    plt.savefig(savefile)


def RMS(binsz, rms, stderr, rmserr, cadence=None, binstep=1,
        timepoints=[], ratio=False, fignum=-20,
        yran=None, xran=None, savefile=None):
  """
  Plot the RMS vs binsize

  Parameters
  ----------
  binsz: 1D ndarray
     Array of bin sizes.
  rms: 1D ndarray
     RMS of dataset at given binsz.
  stderr: 1D ndarray
     Gaussian-noise rms Extrapolation
  rmserr: 1D ndarray
     RMS uncertainty
  cadence: Float
     Time between datapoints in seconds.
  binstep: Integer
     Plot every-binstep point.
  timepoints: List
     Plot a vertical line at each time-points.
  ratio: Boolean
     If True, plot rms/stderr, else, plot both curves.
  fignum: Integer
     Figure number
  yran: 2-elements tuple
     Minimum and Maximum y-axis ranges.
  xran: 2-elements tuple
     Minimum and Maximum x-axis ranges.
  savefile: String
     If not None, name of file to save the plot.

  Uncredited developers
  ---------------------
  Kevin Stevenson  (UCF)
  """

  if np.size(rms) <= 1:
    return

  # Set cadence:
  if cadence is None:
    cadence = 1.0
    xlabel = "Bin size"
  else:
    xlabel = "Bin size  (sec)"

  # Set plotting limits:
  if yran is None:
    #yran = np.amin(rms), np.amax(rms)
    yran = [np.amin(rms-rmserr), np.amax(rms+rmserr)]
    yran[0] = np.amin([yran[0],stderr[-1]])
    if ratio:
      yran = [0, np.amax(rms/stderr) + 1.0]
  if xran is None:
    xran = [cadence, np.amax(binsz*cadence)]

  fs = 14 # Font size
  if ratio:
    ylabel = r"$\beta =$ RMS / std. error"
  else:
    ylabel = "RMS"

  plt.figure(fignum, (8,6))
  plt.clf()
  ax = plt.subplot(111)

  if ratio: # Plot the residuals-to-Gaussian RMS ratio:
    a = plt.errorbar(binsz[::binstep]*cadence, (rms/stderr)[::binstep],
                     (rmserr/stderr)[::binstep], fmt='k-', ecolor='0.5',
                     capsize=0, label="__nolabel__")
    a = plt.semilogx(xran, [1,1], "r-", lw=2)
  else:     # Plot residuals and Gaussian RMS individually:
    # Residuals RMS:
    a = plt.errorbar(binsz[::binstep]*cadence, rms[::binstep],
                     rmserr[::binstep], fmt='k-', ecolor='0.5',
                     capsize=0, label="RMS")
    # Gaussian noise projection:
    a = plt.loglog(binsz*cadence, stderr, color='red', ls='-',
                   lw=2, label="Gaussian std.")
    a = plt.legend()
  for time in timepoints:
    a = plt.vlines(time, yran[0], yran[1], 'b', 'dashed', lw=2)

  a = plt.yticks(size=fs)
  a = plt.xticks(size=fs)
  a = plt.ylim(yran)
  a = plt.xlim(xran)
  a = plt.ylabel(ylabel, fontsize=fs)
  a = plt.xlabel(xlabel, fontsize=fs)

  if savefile is not None:
    plt.savefig(savefile)


def modelfit(data, uncert, indparams, model, nbins=75, title=None,
             fignum=-22, savefile=None):
  """
  Plot the model and (binned) data arrays, and their residuals.

  Parameters
  ----------
  data: 1D float ndarray
     The data array.
  uncert: 1D float ndarray
     Uncertainties of the data-array values.
  indparams: 1D float ndarray
     X-axis values of the data-array values.
  model: 1D ndarray
     The model of data (evaluated at indparams values).
  nbins: Integer
     Output number of data binned values.
  title: String
     Plot title.
  fignum: Integer
     The figure number.
  savefile: Boolean
     If not None, name of file to save the plot.
  """

  # Bin down array:
  binsize = (np.size(data)-1)/nbins + 1
  bindata, binuncert, binindp = ba.binarray(data, uncert, indparams, binsize)
  binmodel = ba.weightedbin(model, binsize)
  fs = 14 # Font-size

  p = plt.figure(fignum, figsize=(8,6))
  p = plt.clf()

  # Residuals:
  a = plt.axes([0.15, 0.1, 0.8, 0.2])
  p = plt.errorbar(binindp, bindata-binmodel, binuncert, fmt='ko', ms=4)
  p = plt.plot([indparams[0], indparams[-1]], [0,0],'k:',lw=1.5)
  p = plt.xticks(size=fs)
  p = plt.yticks(size=fs)
  p = plt.xlabel("x", size=fs)
  p = plt.ylabel('Residuals', size=fs)

  # Data and Model:
  a = plt.axes([0.15, 0.35, 0.8, 0.55])
  if title is not None:
    p = plt.title(title, size=fs)
  p = plt.errorbar(binindp, bindata, binuncert, fmt='ko', ms=4,
                   label='Binned Data')
  p = plt.plot(indparams, model, "b", lw=2, label='Best Fit')
  p = plt.setp(a.get_xticklabels(), visible = False)
  p = plt.yticks(size=13)
  p = plt.ylabel('y', size=fs)
  p = plt.legend(loc='best')

  if savefile is not None:
      p = plt.savefig(savefile)


def reformatparname(parname):
  """
  This function reformats special-case parameter names for plotting purposes.
  Cases: kappa, gamma, alpha, beta, and molecules considered by Transit.

  Inputs
  ------
  parname: List of parameter names.

  Outputs
  -------
  reformatpar: List of reformatted names.
  """
  # Reformat parameter names for special cases when plotting
  reformatpar = []
  for i in range(len(parname)):
    # kappa, gamma1, gamma2, alpha, beta params for Line PT profile
    if   'kappa' in parname[i]:
      reformatpar.append(parname[i].replace('kappa', '$\kappa$'))
    elif 'g1' in parname[i]:
      reformatpar.append(parname[i].replace('g1', '$\gamma_1$'))
    elif 'g2' in parname[i]:
      reformatpar.append(parname[i].replace('g2', '$\gamma_2$'))
    elif parname[i] == 'alpha':
      reformatpar.append('$\\alpha$')
    elif parname[i] == 'beta':
      reformatpar.append('$\\beta$')
    # Radius of the planet, for transit geometry cases
    elif parname[i] == 'Rp' or parname[i] == 'R_p':
      reformatpar.append('R$_\mathrm{p}$')
    # Molecule names -- this needs to be updated whenever new species data 
    # is added
    elif np.any(parname[i] == np.array(['H2O',   'CO2',  'CH4',  'C2H2', 
                                        'C2H4',  'C2H6', 'H2',   'NH3',   
                                        'H2S',   'NO2',  'N2O',  'H2O2', 
                                        'CH3OH', 'H2CO', 'HNO3', 'N2', 'O2'])):
      reformatpar.append(parname[i].replace('2', '$_{2}$')\
                                   .replace('3', '$_{3}$')\
                                   .replace('4', '$_{4}$')\
                                   .replace('6', '$_{6}$'))
    else:
      reformatpar.append(parname[i])

  return reformatpar
