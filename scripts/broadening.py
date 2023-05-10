import sys, os
import numpy as np
import scipy.constants as sc
import ConfigParser

scriptsdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(scriptsdir + "/../code")
import makeatm as ma

def get_widths(cfile):
  """
  Calculate the max and min Lorentz and Doppler broadening HWHM for
  the given configuration file.

  Parameters:
  -----------
  cfile: String
     A BART configuration file.
  """

  # Read config:
  config = ConfigParser.SafeConfigParser()
  config.optionxform = str  # This one enable Uppercase in arguments
  config.read([cfile])
  defaults = dict(config.items("MCMC"))

  # Get min-max wavenumber:
  try:
    wnmin = 1.0/(float(defaults["wlhigh"])*float(defaults["wlfct"]))
  except:
    wnmin = float(defaults["wnlow"]) * float(defaults["wnfct"])
  try:
    wnmax = 1.0/(float(defaults["wllow"])*float(defaults["wlfct"]))
  except:
    wnmax = float(defaults["wnhigh"]) * float(defaults["wnfct"])

  # Read atmospheric file:
  molecs, pressure, temps, abun = ma.readatm(defaults["atmfile"])

  # Get min-max temperatures:
  try:
    tmin = float(defaults["tlow"])
  except:
    tmin = np.amin(temps)
  try:
    tmax = float(defaults["thigh"])
  except:
    tmax = np.amax(temps)

  # Get min-max pressures: 
  pmin = np.amin(pressure) * 1e6
  pmax = np.amax(pressure) * 1e6

  # Get masses:
  molfile = scriptsdir + "/../modules/transit/inputs/molecules.dat"
  ID, mol, mass, diam = readmol(molfile)

  # Keep only molecules from the atmospheric file:
  mask = np.zeros(len(mol), int)
  for i in np.arange(len(mol)):
    mask[i] = mol[i] in molecs

  mol  = mol [np.where(mask)]
  mass = mass[np.where(mask)] * sc.u * 1e3  # grams
  diam = diam[np.where(mask)] * sc.angstrom * 100  # cm

  # Sort molecules according to the order in the atmospheric file:
  isort = np.zeros(len(mol), int)
  for i in np.arange(len(mol)):
    isort[i] = np.where(np.asarray(mol) == molecs[i])[0][0]
  mol  = mol [isort]
  mass = mass[isort]
  diam = diam[isort]

  # Calculate minimum and maximum Doppler widths:
  dmin = Doppler(wnmin, tmin, np.amin(mass))
  dmax = Doppler(wnmax, tmax, np.amax(mass))

  iH2 = np.where(mol == 'H2')[0]
  iHe = np.where(mol == 'He')[0]

  # Calculate minimum and maximum Lorentz widths:
  lmin = Lorentz(pmin, tmax, mass, iH2, iHe, abun[-1], diam, True)
  lmax = Lorentz(pmax, tmin, mass, iH2, iHe, abun[ 0], diam, False)

  print("Doppler minimum and maximum HWHM (cm-1): {:.3e}, {:.3e}\n"
        "Lorentz minimum and maximum HWHM (cm-1): {:.3e}, {:.3e}".
        format(dmin, dmax, lmin, lmax))


def Lorentz(pressure, temperature, mass, iH2, iHe, abundance,
            diameter, min=True):
  """
  Calculate the Lorentz HWHM (in cm-1) for the given parameters:

  Parameters:
  -----------
  pressure: Float
     Pressure in Barye.
  temperature: Float
     Temperature in K.
  mass: 1D float array
     Array of species masses in grams.
  iH2: Integer
     Index in mass correspoding to H2.
  iHe: Integer
     Index in mass correspoding to helium.
  abundance: 1D float array
     Mole mixing ratio of the species.
  diameter: 1D float array
     Collisional diameter of the species in cm.
  min: Bool
     Flag to calculate the minimum (True) or maximum (False) from the species.
  """
  if min:
    func = np.amin
  else:
    func = np.amax

  # Take the minimum or maximum from the species:
  flim = func(abundance[iH2] * ((diameter+diameter[iH2])*0.5)**2.0 *
                np.sqrt((1/mass + 1/mass[iH2]))              +
              abundance[iHe] * ((diameter+diameter[iHe])*0.5)**2.0 *
                np.sqrt((1/mass + 1/mass[iHe]))              )

  # Multiply by the other factors and return:
  return np.sqrt(2)/(sc.c*100) / np.sqrt(temperature*np.pi*(sc.k*1e7)) * pressure * flim


def Doppler(wavenumber, temp, mass):
  """
  Calculate the Doppler HWHM (in cm-1) for the given paramteres.

  Parameters:
  -----------
  wavenumber: Float
    The wavenumber in cm-1.
  temp: Float
    Temperature in Kelvin degrees.
  mass: Float
    Mass of the species in gr.
  """
  return wavenumber/(sc.c*100) * np.sqrt(2*np.log(2)*(sc.k*1e7) * temp/mass)


def readmol(molfile):
  """
  Read and extract species info from molfile.

  Parameters:
  -----------
  molfile: String
     Path to the molecular info file.

  Returns:
  --------
  ID: 1D integer array
     Universal ID for each species.
  mol: 1D string array
     Names of the species.
  mass: 1D float array
     Mass of the species in amu.
  diam: 1D float array
     Diameter of the species in angstroms.
  """
  # Open read the file:
  f = open(molfile, "r")
  lines = f.readlines()
  f.close()

  # Find the line where the species are listed.
  for start in np.arange(len(lines)):
    if lines[start].startswith("# ID"):
      break
  start += 2

  # Extract the species info:
  ID, mol, mass, diam = [], [], [], []
  while lines[start].strip() != "":
    line = lines[start].split()
    ID  .append(line[0])
    mol .append(line[1])
    mass.append(line[2])
    diam.append(line[3])
    start += 1

  return (np.asarray(ID, int), np.asarray(mol),
          np.asarray(mass, np.double), np.asarray(diam, np.double))


