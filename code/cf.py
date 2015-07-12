# ****************************** START LICENSE *******************************
# Bayesian Atmospheric Radiative Transfer (BART), a code to infer
# properties of planetary atmospheres based on observed spectroscopic
# information.
#
# This project was completed with the support of the NASA Planetary
# Atmospheres Program, grant NNX12AI69G, held by Principal Investigator
# Joseph Harrington. Principal developers included graduate students
# Patricio E. Cubillos and Jasmina Blecic, programmer Madison Stemm, and
# undergraduates M. Oliver Bowman and Andrew S. D. Foster.  The included
# 'transit' radiative transfer code is based on an earlier program of
# the same name written by Patricio Rojo (Univ. de Chile, Santiago) when
# he was a graduate student at Cornell University under Joseph
# Harrington.  Statistical advice came from Thomas J. Loredo and Nate
# B. Lust.
#
# Copyright (C) 2015 University of Central Florida.  All rights reserved.
#
# This is a test version only, and may not be redistributed to any third
# party.  Please refer such requests to us.  This program is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.
#
# Our intent is to release this software under an open-source,
# reproducible-research license, once the code is mature and the first
# research paper describing the code has been accepted for publication
# in a peer-reviewed journal.  We are committed to development in the
# open, and have posted this code on github.com so that others can test
# it and give us feedback.  However, until its first publication and
# first stable release, we do not permit others to redistribute the code
# in either original or modified form, nor to publish work based in
# whole or in part on the output of this code.  By downloading, running,
# or modifying this code, you agree to these conditions.  We do
# encourage sharing any modifications with us and discussing them
# openly.
#
# We welcome your feedback, but do not guarantee support.  Please send
# feedback or inquiries to:
#
# Joseph Harrington <jh@physics.ucf.edu>
# Patricio Cubillos <pcubillos@fulbrightmail.org>
# Jasmina Blecic <jasmina@physics.ucf.edu>
#
# or alternatively,
#
# Joseph Harrington, Patricio Cubillos, and Jasmina Blecic
# UCF PSB 441
# 4111 Libra Drive
# Orlando, FL 32816-2385
# USA
#
# Thank you for testing BART!
# ******************************* END LICENSE *******************************

"""
    This code produced contribution functions.
    
    Functions
    ---------
    read_MCMC_out:
          Read the MCMC output log file. Extract the best fitting parameters.
    get_params:
	      Get correct number of all parameters from stepsize
    get_starData:
          Extract stellar temperature, radius, and mass from TEP file
    write_atmfile:
          Write best-fit atm file with scaled H2 and He to abundances sum of 1
    bestFit_tconfig:
          Write best-fit config file for best-fit Transit run
    callTransit:
          Call Transit to produce best-fit outputs. Plot MCMC posterior PT plot.
    plot_bestFit_Spectrum:
          Plot BART best-model spectrum

    Revisions
    ---------
    2015-07-12  Jasmina  Initial version

"""

import numpy as np
from scipy.interpolate import interp1d
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import makeatm as mat
import wine as w
import constants as c



def readTauDat(file, nlayers):
    """
    Read the tau.dat file that carries best-fit tau values.
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


def Plunck(T, wns):
	"""
	Calculate Plunch function.
	"""

	H  = 6.6260755e-27       # Planck's constant (erg * s)
	LS = 2.99792458e10       # Light Speed (cm / s)
	KB = 1.380658e-16        # Boltzmann constant (erg / K) 

    # Number of layers
	nlayers = len(T)

	# Calculate black body function for every layer and every wavenumber
	BB = np.zeros((nlayers, len(wns)))
	for i in np.arange(nlayers):
		for j in np.arange(len(wns)):
			BB[i, j] = (2. * H * wns[j]**3 * LS **(2)) / \
						(np.exp((H * wns[j] * LS) / (KB * T[i])) - 1)

	return BB



def cf_eq(BB, p, tau, nlayers, wns):
	"""
	Applies the contribution function equation as in Knutson et al 2008
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

	# Compensates for Transit's tau=0 after reaching toomuch
	huck = np.linspace(11, 110, 100)
	for i in np.arange(nlayers):
		for j in np.arange(len(wns)):
			if i != 0 and tau[i, j] == 0.0:
				tau[i, j] = huck[i]

	# Calculate BB, reverse the order of p and T
	p = p[::-1]
	T = T[::-1]
	BB = Plunck(T, wns)

	# Call cf_eq() to calculate cf
	cf = cf_eq(BB, p, tau, nlayers, wns)

	# Call filter_cf() to calculate cf
	filt_cf, filt_cf_norm = filter_cf(filters, nlayers, wns, cf)

	############ Plotting ################

	print("  Plotting contribution functions.")
	# Not normalized cf
	plt.figure(1)
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
	plt.figure(2)
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


