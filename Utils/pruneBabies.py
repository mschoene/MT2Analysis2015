#! /usr/bin/python

from ROOT import TFile, TTree
import commands




def pruneBaby( fname, dir, prunedir, pruneBranches ):

  fullname = dir+"/"+fname
  if "pnfs" in fullname : 
    fullname = "dcap://t3se01.psi.ch:22125//"+fullname
  file = TFile.Open(fullname)
  tree = file.Get("mt2")

  print "-> Pruning " + fullname

  
  for i in pruneBranches:
    tree.SetBranchStatus( i+"_*", 0 )
    tree.SetBranchStatus( "n"+i, 0 )

  newfilename = fname.split(".root")[0] + "_prune.root"
  newfile = TFile(prunedir+"/"+newfilename, "recreate")
  print "newfile: " + newfile.GetName()
  newtree = tree.CloneTree()
  newtree.Write("", 5)
  newfile.Close()



if __name__ == '__main__':

   import os
   import sys

   from optparse import OptionParser

   parser = OptionParser()
   parser.usage = ""
   parser.add_option("-f","--filter", dest="filter",
                      default="",
                      help="prune only selected dataset")
   parser.add_option("-x","--useXRD", dest="useXRD",
                      default="false",
                      help="useXRD (default=false -> gfal)")
   parser.add_option("-p","--gfalProtocol", dest="gfalProtocol",
                      default="gsiftp",
                      help="gfal protocol (default and recommended: gsiftp; supported alternative: srm)")

   (options,args) = parser.parse_args()
   if len(args)==0:
     print "ERROR! Script needs at least one argument to run:"
     print "   python pruneBabies.py [input_directory] [output_dir] [configFile=\"pruneBranches.txt\"]"
     exit()

   dir = args[0]
   outdir = args[1]
   pruneBranches = args[2]
   if len(args)>3:
     print "WARNING! Passed more than three arguments to script! Will ignore the additional ones."


   if options.filter !="" :
     print "-> Pruning only files containing: " + str(options.filter)

   if options.useXRD == "true":
     print "-> chosen to use xrootd"
   else:
     print "-> chosen to use gfal via "+gfalProtocol+" protocol"

   pruneBranches = pruneBranches.replace(",", " ")
   pruneBranches = pruneBranches.split()

   print "-> Will prune these branches: " + str(pruneBranches)


   prunedir = outdir
   os.system("mkdir -p " + prunedir)

   if "pnfs/psi.ch" in dir : 
     if options.useXRD == "true":
       status,files = commands.getstatusoutput("xrdfs t3dcachedb.psi.ch ls "+dir)
     else:
       status,files = commands.getstatusoutput("env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-ls "+gfalProtocol+"://t3se01.psi.ch"+dir)

     files=files.splitlines()
   else :
     files = os.listdir(dir)



   for f in files:
     f = f.split("/")[-1] # xrootd takes the full path, truncate
     if ".root" in f:
       if options.filter in f:
         pruneBaby(f, dir, prunedir, pruneBranches)

   #print "Find pruned babies in " + prunedir

