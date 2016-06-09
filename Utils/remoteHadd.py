#! /usr/bin/python



import os
import sys
import getopt


def main(argv):
    
        
  if len(argv) < 1:
      print 'Need at least one arguments...\nUsage: python remoteHadd.py <inputPath> <outputPath (default = inputPath/merged)>'
      sys.exit()

  inputdir = argv[0]
  outputdir = inputdir + "/merged"
  if len(argv)>1: 
    outputdir = argv[1]


  inputPath  = "/pnfs/psi.ch/cms/trivcat/store/user/" + inputdir
  outputPath = "/pnfs/psi.ch/cms/trivcat/store/user/" + outputdir
  #os.system("gfal-mkdir srm://t3se01.psi.ch" + outputPath) # old way
  os.system("xrdfs t3dcachedb.psi.ch mkdir " + outputPath)

  subdirs = os.listdir( inputPath )

  for s in subdirs:

    if s=="merged" : continue

    subdirPath = inputPath + "/" + s 

    print "hadding: " + subdirPath

    outfile = s + ".root"
    command = "hadd -f " + outfile


    files = os.listdir( subdirPath )

    for f in files:

      filePath = subdirPath + "/" + f
      command += " dcap://t3se01.psi.ch:22125" + filePath

    os.system( command )

    #os.system( "gfal-copy file:///`pwd`/" + outfile + " srm://t3se01.psi.ch" + outputPath ) # old way
    os.system( "xrdcp -d 1 `pwd`/" + outfile + " root://t3dcachedb.psi.ch:1094" + outputPath )
  
    os.system("rm " + outfile )
      
 

if __name__ == "__main__":
    main(sys.argv[1:])
