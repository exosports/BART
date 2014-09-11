import sys, traceback, textwrap

# Warning separator:
sep = 70*":"

def msg(verblevel, message, indent=0):
  """
  Conditional message printing to screen.

  Modification History:
  ---------------------
  2014-06-15  patricio  Added Documentation.
  2014-08-18  patricio  Copied to BART project.
  """
  sentences = message.splitlines()
  indspace = " "*indent
  if verblevel > 0:
    for s in sentences:
      msg = textwrap.fill(s, replace_whitespace=True,
                          initial_indent=indspace, subsequent_indent=indspace)
      print(msg)


def warning(message):
  """
  Print message surrounded by colon bands.

  Modification History:
  ---------------------
  2014-06-15  patricio  Initial implementation.
  2014-08-18  patricio  Copied to BART project.
  """
  print("{:s}\n  Warning:".format(sep))
  msg(1, message, 4)
  print("{:s}".format(sep))


def error(message):
  """
  Pretty print error message.

  Modification History:
  ---------------------
  2014-06-15  patricio  Initial implementation.
  2014-08-18  patricio  Copied to BART project.
  """
  # Trace back the file, function, and line where the error source:
  t = traceback.extract_stack()
  # Extract fields:
  efile = t[-2][0]
  efile = efile[efile.rfind('/')+1:]
  efunc = t[-2][2]
  eline = t[-2][1]
  # Indent and wrap message to 70 characters:
  msg = textwrap.fill(message, initial_indent   ="    ",
                               subsequent_indent="    ")
  # Print it out:
  print("{:s}\n  Error in module: '{:s}', function: '{:s}', line: {:d}\n"
        "{:s}\n{:s}".format(sep, efile, efunc, eline, msg, sep))
  sys.exit(0)
