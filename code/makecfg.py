import os, sys
import argparse, ConfigParser
import numpy as np

filedir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(filedir + "/../modules/MCcubed/src/")
import mcutils as mu


def makecfg(cfg, directory, atmfile):
  """
  Make a transit configuration file.

  Parameters:
  -----------
  cfg: String
     BART configuration file
  directory:  String
     Directory where to store the transit configuration file.
  atmfile: String
     Name of the atmospheric filename.

  Modification History:
  ---------------------
  2014-09-09  patricio  Initial implementation.
  """

  # Known (single-value) transit arguments:
  known_args = ["numberabund", "allowq",
                "rad-low", "rad-high", "rad-delt",  "rad-fct",
                "wl-low",  "wl-high",  "wl-delt",   "wl-fct", "wl-osamp",
                "wn-low",  "wn-high",  "wn-delt",   "wn-fct", "wn-osamp",
                "tlow",    "thigh",    "tempdelt",
                "finebin", "nwidth", "blowex", "minelow",
                "opacityfile",
                "toomuch", "tauiso", "outtau", "taulevel", "modlevel",
                "starrad", "transparent",
                "solution-type", "raygrid",
                "cloudrad", "cloudfct", "cloudext",
                "verb",
                "outtoomuch", "outsample", "output",
               ]
  # Known multiple-value transit arguments:
  known_tuple_args = ["linedb", "cia", "orbpars", "orbparsfct"]

  # Get all BART arguments in a dict:
  config = ConfigParser.SafeConfigParser()
  config.read([cfg])
  defaults = dict(config.items("MCMC"))

  # Key names of the arguments in the BART configuration file:
  cfg_args = defaults.keys()

  # transit configuration filename:
  # FINDME: Add file-not-found exception
  transitcfg = defaults["config"]
  tcfg = open(directory + transitcfg, "w")

  # Print the atmospheric name:
  tcfg.write("atm {:s}\n".format(atmfile))

  # Print the known arguments:
  for key in known_tuple_args:
    if key in cfg_args:
      values = defaults[key].strip().split()
      for val in values:
        tcfg.write("{:s} {:s}\n".format(key, val))

  for key in known_args:
    if key in cfg_args:
      tcfg.write("{:s} {:s}\n".format(key, defaults[key]))

  tcfg.close()
