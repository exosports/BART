# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# BART is under an open-source, reproducible-research license (see LICENSE).

import numpy as np

class File:
  """
  Loads a data parameters file and let you querry the parameters and values.
  Useful to read tep files and run files.

  SYNTAX:
      filedata = File(file)

  PURPOSE:
      To read, parse, and return the contents of a tep file
      in a 2D list form.

  RESTRICTIONS:
      This class automatically interprets the type of the values. If they
      can be cast into a numeric value retuns a numeric value, otherwise
      returns a string.

      The parameter 'tepfile' must be of base class FILE with READ privledges!
      This can be obtained by using the built in function 'open'.

      Input MUST be a text file in tep format (.tep).  A similar algorithm
      can be used to read any text file, only requiring a few modest changes
      to this function.


  Examples:
  ---------
  >>> import reader as rd
  >>> reload(rd)
  >>>
  >>> t = rd.File('wasp18b.ev')
  >>> t.listparams()
  >>> t.getvalue('planet')
  >>> t.getvalue('npos')
  >>> t.getvalue('photap')
  >>> t.getvalue('sigma')
  >>> t.getvalue('aorname')
  >>> t.getstr('aorname')
  >>>
  >>> t.getstr('nnod')
  >>> t.getvalue('nnod')

  >>> tep = rd.File('wasp18b.tep')
  >>> tep.getvalue('Rs')
  >>> tep.getvalue('RA')

  Developers:
  -----------
  Christopher Campo  UCF  ccampo@gmail.com
  Patricio Cubillos  UCF  pcubillos@fulbrightmail.org

  Modification History:
  ---------------------
  2009-01-02  Chris      Initial Version
  2010-03-08  Patricio   Modified from CCampo version
  """


  def __init__(self, file):
    # List for the parameter and values
    self.params = np.array([])       
    self.values = []

    # Read the file
    file = open(file, 'r')
    lines = file.readlines()
    file.close()

    for line in lines:
      try:
        line = line[0:line.index('#')].strip()
      except:
        line = line.strip()

      if len(line) > 0:
        self.params = np.append(self.params, line.split()[0])
        self.values.append( np.array(line.split()[1:]) )


  def evaluate(self, value):
    '''
    determines if the value is a numeric expression,
    if not, return a string.
    if it is, returns the numeric value.
    '''
    try:
      v = eval(value)
      return v
    except:
      return value


  def checkpar(self, par):
    ''' 
    check if the input parameter exists. 
    If it does, returns the reference to the values, 
    If not, returns NaN.
    '''
    try:
      id = np.where(self.params == par)[0][0]
      value = self.values[id]
      return value[0] if value.size == 1 else value
    except:
      return np.nan


  def getvalue(self, par):
    '''
    Get the values of a parameter, if it has more than one value, 
    returns an array, if not, returns the value.
    '''
    val = self.checkpar(par)
    if val is np.nan:
      return np.nan

    if val.size > 1:
      value = []   # list that contains the values 
      for i in np.arange(val.size):
        value.append( self.evaluate(val[i]) )
      return np.array(value)

    else:
      value = self.getvalue(val)
      return self.evaluate(val) if np.isnan(value)  else value


  def getstr(self, par):
    '''
    Get the values of a parameter as strings.
    '''
    val = self.checkpar(par)
    return val
