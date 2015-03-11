#! /usr/bin/python


if __name__ == '__main__':

  import os
  import sys
  from optparse import OptionParser
  parser = OptionParser()
  parser.usage = ""
  parser.add_option("-f","--filter", dest="filter",
                    default="",
                    help="skim only selected dataset")

  (options,args) = parser.parse_args()

  dir = args[0]

  signalDirs = os.listdir(dir)

  for i in signalDirs:

    if i.startswith("datacard_"): continue

    signal = i.split("datacards_")[1]

    logFileName = dir + "/" + i + "/limit.log"
    if not os.path.isfile(logFileName): continue

    logFile = open(logFileName, 'r')

    for line in logFile:

      if not line.startswith( "Expected 50.0%" ): continue

      limit = line.split( "Expected 50.0%: r < " )[1]
      spaces = ""
      if signal.startswith("T2"): spaces = "    "
      print signal + spaces + " \t" + limit

      
