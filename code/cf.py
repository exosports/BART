# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

"""
This code produced contribution functions.

Functions
---------
cf_tconfig:
      Write cf_tconfig file for cf Transit run.
readTauDat:
      Read the tau.dat file that carries cf tau values.
Planck:
      Calculate Planch function.
cf_eq:
  Apply the contribution function equation as in Knutson et al 2008
      equation (2).
filter_cf:
      Band-averaged (filters) contribution functions.
cf:
  Call above functions to calculate cf and plot them
"""

import numpy as np
from scipy.interpolate import interp1d
import sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import makeatm as mat
import wine as w
import constants as c


def cf_tconfig(date_dir):
  '''
  Write cf_tconfig file for cf Transit run
  '''
  # Open atmfile to read
  f = open(date_dir + 'bestFit_tconfig.cfg', 'r')
  lines = np.asarray(f.readlines())
  f.close()

  # Find where toomuch argument is written and put 500.0
  data = [[] for x in range(len(lines))]
  for i in np.arange(len(lines)):
    data[i] = lines[i].split()
    if data[i][0] == 'toomuch':
      lines[i] = 'toomuch 1e100\n'
    elif data[i][0] == 'verb':
      lines[i] = 'verb 0\n'
    elif data[i][0] == 'outspec':
      lines[i] = 'outspec ./cf-flux.dat\n'
    elif data[i][0] == 'outintens':
      lines[i] = 'outintens ./cf-intens.dat\n'
    elif data[i][0] == 'outtoomuch':
      lines[i] = 'outtoomuch ./cf-toom.dat\n'

  # Write lines into the bestFit config file
  f = open(date_dir + 'cf_tconfig.cfg', 'w')
  f.writelines(lines)
  f.writelines('savefiles yes')
  f.close()


def readTauDat(file, nlayers):
  """
  Read the tau.dat file that carries cf tau values.
  """
  # Open and read the filter file:
  data = open(file, "r")
  lines = data.readlines()
  data.close()

  # Remove header comments and empty lines:
  while lines[0].startswith("#") or not lines[0].strip():
      comment = lines.pop(0)

  # Take only lines with tau  and wns data,
  #      skipps the lines with other information
  tau_lines = lines[1:-1:3]
  wn_lines  = lines[0:-1:3]
  tau = np.zeros((len(tau_lines), nlayers))
  wns = np.zeros(len(wn_lines))
  for i in np.arange(len(tau_lines)):
      tau[i] = tau_lines[i].split()
      wns[i] = np.float(wn_lines[i].split()[1])

  # Transpose the order of tau array
  tau = tau.T

  return tau, wns


def Planck(T, wns):
  """
  Calculate Planch function.
  """
  # Number of layers
  nlayers = len(T)

  # Calculate black body function for every layer and every wavenumber
  BB = np.zeros((nlayers, len(wns)))
  for i in np.arange(nlayers):
    for j in np.arange(len(wns)):
      BB[i, j] = (2.0 * c.H * wns[j]**3.0 * c.LS **2.0) / \
                  (np.exp((c.H * wns[j] * c.LS) / (c.KB * T[i])) - 1.0)

  return BB


def cf_eq(BB, p, tau, nlayers, wns):
  """
  Apply the contribution function equation as in Knutson et al. (2008)
  equation (2).
  """
  # Calculate change in exp(-tau) and log(p), and cf
  de_tau = np.zeros((nlayers, len(wns)))
  d_logp = np.zeros(nlayers)
  cf = np.zeros((nlayers, len(wns)))
  for i in np.arange(nlayers):
    # Outermost layer
    if i==0:
      de_tau[i, :] = 0.0
      d_logp[i]    = 0.0
      cf[i]        = 0.0
    else:
      de_tau[i, :] = np.exp(-tau[i-1, :]) - np.exp(-tau[i, :])
      d_logp[i] = np.log(p[i]*1e6) - np.log(p[i-1]*1e6)
      cf[i] = BB[i, :] * de_tau[i, :]/d_logp[i]

  return cf


def filter_cf(filters, nlayers, wns, cf, normalize=False):
  """
  Band-averaged (filters) contribution functions.
  """
  # Number of filters
  nfilters = len(filters)
  # Allocate arrays for filter cf and normalize cf
  filt_cf = np.zeros((nfilters, nlayers))
  filt_cf_norm = np.zeros((len(filters), nlayers))

  for i in np.arange(nfilters):
    # Read filter
    wn, response = w.readfilter(filters[i])

    # Find where filters starts and ends
    wn_filt = wns[(wns>min(wn)) * (wns<max(wn))]
    start_filt, stop_filt = np.where(wns==min(wn_filt))[0][0], \
                            np.where(wns==max(wn_filt))[0][0]

    # Interpolate filter response functions
    inter_filt = interp1d(wn, response)
    resp_filt = inter_filt(wn_filt)

    # Extract filter contribution functions
    cf_filt = cf[:, start_filt:stop_filt+1]

    # Extract filter response functions
    cf_filt_resp = np.zeros((nlayers, len(wn_filt)))
    for j in np.arange(len(wn_filt)):
      cf_filt_resp[:,j] = cf_filt[:,j] * resp_filt[j]

    # Integrate cf across bandpass (filter)
    for k in np.arange(nlayers):
      filt_cf[i,k] = np.trapz(cf_filt_resp[k,:]) / np.trapz(resp_filt)

    if normalize:
      # Normalize to 1
      filt_cf_norm[i] = (filt_cf[i] - min(filt_cf[i])) / \
                        (max(filt_cf[i]) - min(filt_cf[i]))

  if normalize:
    return filt_cf, filt_cf_norm

  return filt_cf


def transmittance(date_dir, atmfile, filters, fext='.png', plot=True):
  """
  """
  # Read atmfile
  molecules, p, T, abundances = mat.readatm(date_dir + atmfile)
  nlayers = len(p)
  p = p[::-1]  # top to bottom of the atmosphere
  # Read tau.dat
  foo      = date_dir + 'tau.dat'
  tau, wns = readTauDat(foo, nlayers)
  # Transmittance:
  transmit = np.exp(-tau)
  # Band intgrate it:
  filt_tr  = filter_cf(filters, nlayers, wns, transmit)
  nfilters = len(filters)

  meanwl = np.zeros(nfilters)
  for i in np.arange(nfilters):
    filtwaven, filttransm = w.readfilter(filters[i])
    meanwn    = np.sum(filtwaven*filttransm)/np.sum(filttransm)
    meanwl[i] = 1e4/meanwn
  maxmeanwl   = np.max(meanwl)
  minmeanwl   = np.min(meanwl)
  colors      = (meanwl-minmeanwl) / (maxmeanwl-minmeanwl)

  if plot:
    print("  Plotting contribution functions.\n")
    # Not normalized cf
    plt.figure(4, figsize=(8,6.5))
    plt.clf()
    gs       = gridspec.GridSpec(1, 2, width_ratios=[20, 1], wspace=0.05)
    ax0      = plt.subplot(gs[0])
    ax1      = plt.subplot(gs[1])
    for i in np.arange(nfilters):
      (head, tail) = os.path.split(filters[i])
      lbl         = tail[:-4]
      ax0.semilogy(filt_tr[i], p, '-', color=plt.cm.rainbow(colors[i]), 
                   linewidth = 1, label=lbl)
    norm = matplotlib.colors.Normalize(vmin=minmeanwl, vmax=maxmeanwl)
    cbar = matplotlib.colorbar.ColorbarBase(ax1, cmap=plt.cm.rainbow,
                                            norm=norm,
                                            orientation='vertical')
    cbar.set_label("Mean Wavelength (\u00b5m)", fontsize=14)
    ax0.set_ylim(max(p), min(p))
    ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax0.set_xlabel('Transmittance', fontsize=14)
    ax0.set_ylabel('Pressure (bar)' , fontsize=14)
    plt.savefig(date_dir + 'Transmittance' + fext, bbox_inches='tight')
    plt.close()

  return filt_tr[:,::-1]


def cf(date_dir, atmfile, filters, fext='.png', plot=True, fs=15):
  """
  Call above functions to calculate cf and plot them
  """
  # Read atmfile
  molecules, p, T, abundances = mat.readatm(date_dir + atmfile)
  nlayers = len(p)

  # Read tau.dat
  foo      = date_dir + 'tau.dat'
  tau, wns = readTauDat(foo, nlayers)

  # Calculate BB, reverse the order of p and T
  p  = p[::-1]
  T  = T[::-1]
  BB = Planck(T, wns)

  # Call cf_eq() to calculate cf
  cf = cf_eq(BB, p, tau, nlayers, wns)

  # Call filter_cf() to calculate cf
  filt_cf, filt_cf_norm = filter_cf(filters, nlayers, wns, cf, normalize=True)
  nfilters = len(filters)

  meanwl = np.zeros(nfilters)
  for i in np.arange(nfilters):
    filtwaven, filttransm = w.readfilter(filters[i])
    meanwn    = np.sum(filtwaven*filttransm)/np.sum(filttransm)
    meanwl[i] = 1e4/meanwn
  maxmeanwl   = np.max(meanwl)
  minmeanwl   = np.min(meanwl)
  colors      = (meanwl-minmeanwl) / (maxmeanwl-minmeanwl)

  if plot:
    print("  Plotting contribution functions.\n")
    # Not normalized cf
    plt.figure(4)
    plt.clf()
    gs       = gridspec.GridSpec(1, 2, width_ratios=[20, 1], wspace=0.05)
    ax0      = plt.subplot(gs[0])
    ax1      = plt.subplot(gs[1])
    for i in np.arange(len(filt_cf)):
      (head, tail) = os.path.split(filters[i])
      lbl          = tail[:-4]
      ax0.semilogy(filt_cf[i], p, '-', color=plt.cm.rainbow(colors[i]), 
                   linewidth = 1, label=lbl)
    norm = matplotlib.colors.Normalize(vmin=minmeanwl, vmax=maxmeanwl)
    cbar = matplotlib.colorbar.ColorbarBase(ax1, cmap=plt.cm.rainbow,
                                            norm=norm,
                                            orientation='vertical')
    cbar.set_label("Mean Wavelength (\u00b5m)", fontsize=fs)
    cbar.ax.tick_params(labelsize=fs-4)
    ax0.set_ylim(max(p), min(p))
    ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax0.set_xlabel('Contribution Functions', fontsize=fs)
    ax0.set_ylabel('Pressure (bar)' , fontsize=fs)
    plt.xticks(size=fs-4)
    plt.yticks(size=fs-4)
    plt.savefig(date_dir + 'ContrFuncs' + fext, bbox_inches='tight')
    plt.close()

    # Normalized cf
    plt.figure(5, figsize=(8,6.5))
    plt.clf()
    gs       = gridspec.GridSpec(1, 2, width_ratios=[20, 1], wspace=0.05)
    ax0      = plt.subplot(gs[0])
    ax1      = plt.subplot(gs[1])
    for i in np.arange(len(filt_cf_norm)):
      (head, tail) = os.path.split(filters[i])
      lbl          = tail[:-4]
      ax0.semilogy(filt_cf_norm[i], p, '--', color=plt.cm.rainbow(colors[i]), 
                   linewidth = 1, label=lbl)
    norm = matplotlib.colors.Normalize(vmin=minmeanwl, vmax=maxmeanwl)
    cbar = matplotlib.colorbar.ColorbarBase(ax1, cmap=plt.cm.rainbow,
                                            norm=norm,
                                            orientation='vertical')
    cbar.set_label("Mean Wavelength (\u00b5m)", fontsize=fs)
    cbar.ax.tick_params(labelsize=fs-4)
    ax0.set_ylim(max(p), min(p))
    ax0.set_xlim(0, 1.0)
    ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax0.set_xlabel('Normalized Contribution Functions', fontsize=fs)
    ax0.set_ylabel('Pressure (bar)' , fontsize=fs)
    plt.xticks(size=fs-4)
    plt.yticks(size=fs-4)
    plt.savefig(date_dir + 'NormContrFuncs' + fext, bbox_inches='tight')
    plt.close()

  return filt_cf[:,::-1], filt_cf_norm[:,::-1]


