#! /usr/bin/python


if __name__ == '__main__':

  import os
  import sys
  from optparse import OptionParser
  parser = OptionParser()
  parser.usage = ""
  parser.add_option("-s","--significance", dest="significance",
                    default="",
                    help="compute significance instead of limit")
  parser.add_option("-t","--toys", dest="toys",
                    default="",
                    help="use toys instead of asymptotic")


  (options,args) = parser.parse_args()

  dir = args[0]

  signalDirs = os.listdir(dir)

  for i in signalDirs:

    if( i=="datacard_templates" ): continue
    print "Signal: " + i
    os.system("combineCards.py " + dir + "/" + i + "/datacard_HT*txt > " + dir + "/" + i + "/datacard_combine.txt")
    print "Computing limit"
    os.system("combine -M Asymptotic " + dir + "/" + i + "/datacard_combine.txt | tee " + dir + "/" + i + "/limit.log")

  os.system("rm roostats*")

