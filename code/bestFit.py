import numpy as np
import reader as rd
import scipy.constants as sc
import scipy.special   as sp

import makeatm as mat
import PT as pt

# read MCMC out
def read_MCMC_out(MCfile):

    # open file to read
    f = open(MCfile, 'r')
    lines = np.asarray(f.readlines())
    f.close() 

    # find where the data starts and ends
    datalines = lines[1:]
    ndata = len(datalines)
    data = np.zeros((ndata, 4))
    for i in np.arange(ndata):
            data[i] = datalines[i].split()
    bestP = data[:, 0]
    uncer = data[:, 1]
    SN    = data[:, 2] 
    mean  = data[:, 3] 

    return bestP, uncer, SN, mean

# get correct number of all parameters from stepsize
def get_params(bestP, stepsize):

    j = 0
    params = np.zeros(len(stepsize))
    for i in np.arange(len(stepsize)):
        if stepsize[i] != 0.0:
            params[i] = bestP[j]
            j +=1

    return params


# reads the tep file and calculates surface gravity
def get_starData(tepfile):

    # Sun mass, radius and AU:
    # Source: http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
    #         http://neo.jpl.nasa.gov/glossary/au.html
    Rsun = 696.0 * 1e3          # m
    AU   = 149597870.7 * 1e3    # m

    # Open tepfile to read and get data:
    tep = rd.File(tepfile)

    # Get star mass in Mjup:
    T_star = np.float(tep.getvalue('Ts')[0])

    # Get star radius in units of Rjup:
    star_rad = np.float(tep.getvalue('Rs')[0])

    # Get semi major axis in units of AU:
    a = np.float(tep.getvalue('a')[0])

    # Mass, radius, semi major axis in MKS units:
    R_star = star_rad  * Rsun  # m
    sma    = a * AU            # m

    return R_star, T_star, sma

# write best-fit atm file with scaled H2 and He to abundances sum of 1
def write_atmfile(atmfile, molfit, T_line, params, date_dir):

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
    abun_fact = params[5:]

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

    #print T_line
    #print type(T_line)
    # TEMPORARY HACK!
    T_line = np.array([ 1592.61,  1592.49,  1592.22,  1591.62,  1590.43,  1588.31,
        1584.9 ,  1579.83,  1572.82,  1563.69,  1552.41,  1539.06,
        1523.84,  1507.05,  1489.04,  1470.21,  1450.94,  1431.61,
        1412.56,  1394.09,  1376.43,  1359.8 ,  1344.32,  1330.07,
        1317.1 ,  1305.41,  1294.96,  1285.7 ,  1277.56,  1270.44,
        1264.25,  1258.9 ,  1254.31,  1250.37,  1247.01,  1244.15,
        1241.73,  1239.68,  1237.95,  1236.49,  1235.27,  1234.24,
        1233.38,  1232.66,  1232.06,  1231.56,  1231.14,  1230.79,
        1230.51,  1230.26,  1230.06,  1229.9 ,  1229.76,  1229.65,
        1229.55,  1229.47,  1229.41,  1229.36,  1229.31,  1229.28,
        1229.25,  1229.22,  1229.2 ,  1229.19,  1229.17,  1229.16,
        1229.15,  1229.15,  1229.14,  1229.14,  1229.13,  1229.13,
        1229.13,  1229.12,  1229.12,  1229.12,  1229.12,  1229.12,
        1229.12,  1229.12,  1229.12,  1229.12,  1229.12,  1229.12,
        1229.12,  1229.12,  1229.12,  1229.12,  1229.12,  1229.12,
        1229.12,  1229.12,  1229.12,  1229.12,  1229.12,  1229.12,
        1229.12,  1229.12,  1229.12,  1229.12])


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
    #atm_line = 'atm /home/jasmina/BART_WASP43b/BART-new/run/WASP-43b-short/WASP43b_noHCN_atm.tea\n'
    lines[0] = atm_line
    f = open(date_dir + 'bestFit_tconfig.cfg', 'w')
    f.writelines(lines)
    f.writelines('savefiles yes')
    f.close()


# call above functions to prepare Transit for best-fit execution
def callTransit(atmfile, tepfile, MCfile, stepsize, molfit, tconfig, date_dir):

    # read atmfile
    molecules, pressure, temp, abundances = mat.readatm(atmfile)

    # get surface gravity
    grav, Rp = mat.get_g(tepfile)

    # get star data
    R_star, T_star, sma = get_starData(tepfile)

    # get best parameters
    bestP, uncer, SN, mean = read_MCMC_out(MCfile)

    # get all params
    params = get_params(bestP, stepsize)

    # get PTparams and abundances factors
    nparams = len(params)
    nmol = len(molfit)
    nPTparams = nparams - nmol
    PTparams  = params[:nPTparams]
    AbunFact  = params[nPTparams:]

    # HARDCODED !!!!!!!!!!!!
    T_int = 100 # K

    # call PT line profile to calculate temperature
    T_line = pt.PT_line(pressure, PTparams, R_star, T_star, T_int, sma, grav)

    # write best-fit atmospheric file
    write_atmfile(atmfile, molfit, T_line, params, date_dir)

    # bestFit atm file
    bestFit_atm = date_dir + 'bestFit.atm'

    # write new bestFit Transit config
    bestFit_tconfig(tconfig, date_dir)


