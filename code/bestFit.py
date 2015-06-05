import numpy as np
import reader as rd
import scipy.constants as sc
import scipy.special   as sp
import matplotlib.pyplot as plt

import makeatm   as mat
import PT        as pt
import constants as c

def read_MCMC_out(MCfile):
    """
    Read the MCMC output log file. Extract the best fitting parameters.
    """
    # Open file to read
    f = open(MCfile, 'r')
    lines = np.asarray(f.readlines())
    f.close() 

    # Find where the data starts and ends:
    for ini in np.arange(len(lines)):
        if lines[ini].startswith(' Best-fit params'):
            break
    ini += 1
    end = ini
    for end in np.arange(ini, len(lines)):
        if lines[end].strip() == "":
            break

    # Read data:
    data = np.zeros((end-ini, 4))
    for i in np.arange(ini, end):
        data[i-ini] = lines[i].split()
    bestP = data[:, 0]
    uncer = data[:, 1]
    SN    = data[:, 2] 
    mean  = data[:, 3] 

    return bestP, uncer, SN, mean


def get_params(bestP, stepsize, params):
    """
    Get correct number of all parameters from stepsize
    """
    j = 0
    allParams = np.zeros(len(stepsize))
    for i in np.arange(len(stepsize)):
        if stepsize[i] != 0.0:
            allParams[i] = bestP[j]
            j +=1
        else:
            allParams[i] = params[i]
            
    return allParams


def get_starData(tepfile):
  """
  Extract the Stellar temperature, radius, and mass from a TEP file.
  """
  # Open tepfile to read and get data:
  tep = rd.File(tepfile)

  # Get star mass in Mjup:
  Tstar = np.float(tep.getvalue('Ts')[0])
  # Get star radius in MKS units:
  Rstar = np.float(tep.getvalue('Rs')[0]) * c.Rsun
  # Get semi major axis in meters:
  sma = np.float(tep.getvalue('a')[0]) * sc.au

  return Rstar, Tstar, sma


def write_atmfile(atmfile, molfit, T_line, allParams, date_dir):
    """
    Write best-fit atm file with scaled H2 and He to abundances sum of 1.
    """
    # Open atmfile to read
    f = open(atmfile, 'r')
    lines = np.asarray(f.readlines())
    f.close()

    # Get the molecules
    imol = np.where(lines == "#SPECIES\n")[0][0] + 1
    molecules = lines[imol].split()

    # Find the line where the layers info begins
    start = np.where(lines == "#TEADATA\n")[0][0] + 2 
    headers = lines[start-1].split()
    datalines = lines[start:]

    # Number of columns
    ncol = len(lines[start].split())

    # Number of layers
    ndata = len(datalines)  

    # Allocate space for rad and pressure
    rad = np.zeros(ndata, np.double)
    pressure = np.zeros(ndata, np.double) 

    # Number of abundances (elements per line except Radius, Press and T) 
    nabun = len(lines[start].split()) - 3  
    abundances = np.zeros((ndata, nabun), np.double)

    # Extract radius, pressure, and abundance data
    data = np.zeros((ndata, len(headers)))
    for i in np.arange(ndata):
            data[i] = datalines[i].split()
    rad      = data[:, 0]
    pressure = data[:, 1] 

    # fill out the abundances array
    abundances = np.zeros((len(molecules), ndata))
    for i in np.arange(len(molecules)):
        for j in np.arange(ndata):
            abundances[i] = data[:, i+3] 

    # recognize which columns to take from the atmospheric file
    headers = lines[start-1].split()
    columns = np.zeros(len(molfit))
    for i in np.arange(len(molfit)):
        for j in np.arange(len(headers)):
            if molfit[i] == headers[j]:
                columns[i] = j

    # number of molecules to fit
    nfit = len(molfit)
    abun_fact = allParams[5:]

    # multiply the abundances of molfit molecules
    for i in np.arange(len(columns)):
       abundances[columns[i]-3] = abundances[columns[i]-3] * 10**abun_fact[i]

    # ===== Scale H2 and He if sum abundances > 1 ===== #
    # Find index for Hydrogen and Helium
    molecules = np.asarray(molecules)
    iH2     = np.where(molecules=="H2")[0][0]
    iHe     = np.where(molecules=="He")[0][0]

    # Get H2/He abundance ratio:
    ratio = (abundances[iH2,:] / abundances[iHe,:])

    # Find which level has sum(abundances)>1 and get the difference
    q = np.sum(abundances, axis=0) - 1

    # Correct H2, and He abundances respecting their ration
    for i in np.arange(ndata):
        if q[i]>0:
            abundances[iH2, i] -= ratio[i] * q[i] / (1.0 + ratio[i])
            abundances[iHe, i] -=            q[i] / (1.0 + ratio[i])

    # open best fit atmospheric file
    fout = open(date_dir + 'bestFit.atm', 'w')
    fout.writelines(lines[:start])

    # Write atm file for each run
    for i in np.arange(ndata): 
        # Radius, pressure, and temp for the current line
        radi = str('%10.3f'%rad[i])
        presi = str('%10.4e'%pressure[i])
        tempi = str('%7.2f'%T_line[i])

        # Insert radii array
        fout.write(radi.ljust(10) + ' ')

        # Insert results from the current line (T-P) to atm file
        fout.write(presi.ljust(10) + ' ')
        fout.write(tempi.ljust(7) + ' ')

        # Write current abundances
        for j in np.arange(len(headers) - 3):
            fout.write('%1.4e'%abundances[j][i] + ' ')
        fout.write('\n')

    # Close atm file
    fout.close()

# write best-fit config file for new Transit run
def bestFit_tconfig(tconfig, date_dir):

    # Open atmfile to read
    f = open(date_dir + tconfig, 'r')
    lines = np.asarray(f.readlines())
    f.close()

    atm_line = 'atm ' + date_dir + 'bestFit.atm' + '\n' 
    lines[0] = atm_line
    f = open(date_dir + 'bestFit_tconfig.cfg', 'w')
    f.writelines(lines)
    f.writelines('savefiles yes')
    f.close()


# call above functions to prepare Transit for best-fit execution, plot MCMC PT profiles
def callTransit(atmfile, tepfile, MCfile, stepsize, molfit, tconfig, date_dir, params, burnin):
    '''
    Call Transit to produce best-fit outputs.

    '''
    # read atmfile
    molecules, pressure, temp, abundances = mat.readatm(atmfile)

    # get surface gravity
    grav, Rp = mat.get_g(tepfile)

    # convert gravity to cm/s^2
    grav = grav*1e2

    # get star data
    R_star, T_star, sma = get_starData(tepfile)

    # get best parameters
    bestP, uncer, SN, mean = read_MCMC_out(MCfile)

    # get all params
    allParams = get_params(bestP, stepsize, params)

    # get PTparams and abundances factors
    nparams = len(allParams)
    nmol = len(molfit)
    nPTparams = nparams - nmol
    PTparams  = allParams[:nPTparams]

    # HARDCODED !
    T_int = 100 # K

    # call PT line profile to calculate temperature
    best_T = pt.PT_line(pressure, PTparams, R_star, T_star, T_int, sma, grav)

    # Plot best PT profile
    plt.figure(1)
    plt.semilogy(best_T, pressure, '-', color = 'r')
    plt.xlim(0.9*min(best_T), 1.1*max(best_T))
    plt.ylim(max(pressure), min(pressure))
    plt.title('Best PT', fontsize=14)
    plt.xlabel('T [K]'     , fontsize=14)
    plt.ylabel('logP [bar]', fontsize=14)

    # Save plot to current directory
    plt.savefig(date_dir + 'Best_PT.png') 

    # write best-fit atmospheric file
    write_atmfile(atmfile, molfit, best_T, allParams, date_dir)

    # bestFit atm file
    bestFit_atm = date_dir + 'bestFit.atm'

    # write new bestFit Transit config
    bestFit_tconfig(tconfig, date_dir)

    # plot MCMC PT profiles ===================================================

    # get MCMC data:
    MCMCdata = date_dir + "/output.npy"
    data = np.load(MCMCdata)
    nchains, npars, niter = np.shape(data)

    # stuck chains:
    data_stack = data[0,:,burnin:]
    for c in np.arange(1, nchains):
        data_stack = np.hstack((data_stack, data[c, :, burnin:]))

    # create array of PT profiles
    PTprofiles = np.zeros((np.shape(data_stack)[1], len(pressure)))

    # current PT parameters for each chain, iteration
    curr_PTparams = PTparams

    # fill-in PT profiles array
    print("  Plotting MCMC PT profile figure.")
    for k in np.arange(0, np.shape(data_stack)[1]):
        j = 0
        for i in np.arange(len(PTparams)):
            if stepsize[i] != 0.0:
                curr_PTparams[i] = data_stack[j,k]
                j +=1
            else:
                pass
        PTprofiles[k] = pt.PT_line(pressure, curr_PTparams, R_star, T_star, T_int, sma, grav)
        if k %10 == 0:
            print(k)

    # get percentiles (for 1,2-sigma boundaries):
    low1 = np.percentile(PTprofiles, 16.0, axis=0)
    hi1  = np.percentile(PTprofiles, 84.0, axis=0)
    low2 = np.percentile(PTprofiles,  2.5, axis=0)
    hi2  = np.percentile(PTprofiles, 97.5, axis=0)
    median = np.median(PTprofiles, axis=0)

    # plot figure
    plt.figure(2)
    ax=plt.subplot(111)
    ax.fill_betweenx(pressure, low2, hi2, facecolor="#62B1FF", edgecolor="0.5")
    ax.fill_betweenx(pressure, low1, hi1, facecolor="#1873CC", edgecolor="#1873CC")
    plt.semilogy(median, pressure, "-", lw=2, label='Median',   color="limegreen")
    plt.semilogy(best_T, pressure, "-", lw=2, label="Best fit", color="r")
    plt.ylim(pressure[0], pressure[-1])
    plt.legend(loc="best")
    plt.xlabel("Temperature  (K)", size=15)
    plt.ylabel("Pressure  (bar)",  size=15)

    # save figure
    savefile = date_dir + "MCMC_PTprofiles.ps" 
    plt.savefig(savefile, transparent=True)
