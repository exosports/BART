# ::: Frequently-Used Scripts :::
# -------------------------------

# Index:
# ( 0) Make a log-scaled pressure profile.
# ( 1) Make a Temperature profile.
# ( 2) Make an atmospheric file with uniform abundances.
# ( 3) Read and plot a transit spectrum.
# ( 4) Calculate the minimum and maximum line-profile widths.
# ( 5) Plot the posterior TP profiles.

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Preamble:
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

# Assuming that the current working directory is /home/.../BART/scripts/
sys.path.append("../code/")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ( 0) Make a log-scaled pressure profile:
import makeP as mp

# Output pressure file:
pfile = "layers.press"
# Pressure variables:
nlayers = 100
ptop    = 1e-5 # bar
pbottom = 1e2  # bar

# Write the pressure to file:
mp.makeP(nlayers, ptop, pbottom, pfile, log=True)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ( 1) Make a Temperature profile using model from Line et al. (2013):
import PT as PT

# [Run script ( 0) to make a pressure profile file]

# Read the pressure file to an array:
press = PT.read_press_file(pfile)

# System parameters:
Rplanet = 1.0*6.995e8 # m
Tstar   = 5700.0      # K
Tint    =  100.0      # K
gplanet = 2200.0      # cm s-2
smaxis  = 0.050*sc.au # m

# Fitting parameters: [log10(kappa), log10(g1), log10(g2), alpha, beta]
params = np.asarray(  [-2.0,         -0.55,     -0.8,      0.5,   1.0])

# Calculate the temperature profile:
temp = PT.PT_line(press, params, Rplanet, Tstar, Tint, smaxis, gplanet)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ( 2) Make an atmospheric file with uniform abundances:
import makeatm as ma

# Output filename:
atmfile  = "uniform.atm"
# Elemental abundances file:
elemabun = "../inputs/abundances_Asplund2009.txt"
# Transiting extrasolar planet filename:
tep      = "../inputs/tep/HD209458b.tep"

# Atmospheric species:
species    = "He    H2    CO    CO2   CH4   H2O   NH3    C2H2   C2H4"
# Abundances (mole mixing ratio):
abundances = "0.15  0.85  1e-4  1e-4  1e-4  1e-4  1e-10  1e-10  1e-10"

# [Run script ( 2) to make a temperature profile]

# Make the atmospheric file:
ma.uniform(atmfile, pfile, elemabun, tep, species, abundances, temp)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ( 3) Read and plot a transit spectrum:
import readtransit as rt

wl, spectrum = rt.readspectrum("eclipse_out.dat.-Flux", 0)

plt.figure(1)
plt.clf()
plt.semilogx(wl, spectrum, "b", label="Planet spectrum")
plt.xlim(wl[-1], wl[0])
plt.legend(loc="best")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ( 4) Calculate the minimum and maximum line-profile widths:

# TBD

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ( 5) Plot the posterior TP profiles:

# TBD
