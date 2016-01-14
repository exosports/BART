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

    Revisions
    ---------
    2015-07-14  Jasmina  Original implementation

"""

import numpy as np
from scipy.interpolate import interp1d
import os
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
		elif data[i][0] == 'outflux':
			lines[i] = 'outflux ./cf-flux.dat\n'
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
			BB[i, j] = (2. * c.H * wns[j]**3 * c.LS **(2)) / \
						(np.exp((c.H * wns[j] * c.LS) / (c.KB * T[i])) - 1)

	return BB


def cf_eq(BB, p, tau, nlayers, wns):
	"""
	Apply the contribution function equation as in Knutson et al 2008
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


def filter_cf(filters, nlayers, wns, cf):
	"""
	Band-averaged (filters) contribution functions.
	""" 
	# Allocate arrays for filter cf and normalize cf
	filt_cf      = np.zeros((len(filters), nlayers))
	filt_cf_norm = np.zeros((len(filters), nlayers))

	# Number of filters
	nfilters = len(filters)
	for i in np.arange(nfilters):

		# Read filter
		wn, response = w.readfilter(filters[i])

        # Find where filters starts and ends
		wn_filt = np.asarray(filter(lambda x: x>min(wn) and x<max(wn), wns))
		start_filt, stop_filt = np.where(wns==min(wn_filt))[0][0], np.where(wns==max(wn_filt))[0][0]

		# Interpolate filter response functions
		inter_filt = interp1d(wn, response)
		resp_filt = inter_filt(wn_filt)

		# Extract filter contribution functions
		cf_filt = cf[:, start_filt:stop_filt+1]

        # Extract filter response functions
		cf_filt_resp = np.zeros((nlayers, len(wn_filt)))
		for j in np.arange(len(wn_filt)):
			cf_filt_resp[:, j] = cf_filt[:, j]*resp_filt[j]

        # Integrate cf across bandpass (filter)
		for k in np.arange(nlayers):
			filt_cf[i, k] = np.trapz(cf_filt_resp[k, :])

		# Normalize to 1
		filt_cf_norm[i] = (filt_cf[i] - min(filt_cf[i])) / (max(filt_cf[i]) - min(filt_cf[i]))

	return filt_cf, filt_cf_norm 
	

def cf(date_dir, atmfile, filters):
	"""
	Call above functions to calculate cf and plot them
	"""
	# Read atmfile
	molecules, p, T, abundances = mat.readatm(atmfile)
	nlayers = len(p)

	# Read tau.dat
	file = date_dir + 'tau.dat'
	tau, wns = readTauDat(file, nlayers)

	# Calculate BB, reverse the order of p and T
	p = p[::-1]
	T = T[::-1]
	BB = Planck(T, wns)

	# Call cf_eq() to calculate cf
	cf = cf_eq(BB, p, tau, nlayers, wns)

	# Call filter_cf() to calculate cf
	filt_cf, filt_cf_norm = filter_cf(filters, nlayers, wns, cf)

	############ Plotting ################

	print("  Plotting contribution functions.\n")
	# Not normalized cf
	plt.figure(4)
	plt.clf()
	gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1]) 
	ax0 = plt.subplot(gs[0])
	for i in np.arange(len(filt_cf)):
		(head, tail) = os.path.split(filters[i])
		lbl = tail[:-4]
		ax0.semilogy(filt_cf[i], p, '-', linewidth = 1, label=lbl)
	ax0.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), prop={'size':8})
	ax0.set_ylim(max(p), min(p))
	ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	ax0.set_xlabel('Contribution functions', fontsize=14)
	ax0.set_ylabel('Pressure [bar]' , fontsize=14)
	plt.savefig(date_dir + 'ContrFuncs.png')

	# Normalized cf
	plt.figure(5)
	plt.clf()
	gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1]) 
	ax0 = plt.subplot(gs[0])
	for i in np.arange(len(filt_cf_norm)):
		(head, tail) = os.path.split(filters[i])
		lbl = tail[:-4]
		ax0.semilogy(filt_cf_norm[i], p, '--', linewidth = 1, label=lbl)

	ax0.legend(loc='center left', bbox_to_anchor=(1,0.5), prop={'size':8})
	ax0.set_ylim(max(p), min(p))
	ax0.set_xlim(0, 1.0)
	ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	ax0.set_xlabel('Normalized contribution functions', fontsize=14)
	ax0.set_ylabel('Pressure [bar]' , fontsize=14)
	plt.savefig(date_dir + 'NormContrFuncs.png')


