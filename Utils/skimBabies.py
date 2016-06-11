#! /usr/bin/python

from ROOT import TFile, TTree
import commands




def skimBaby( fname, dir, skimdir, cuts ):

  fullname = dir+"/"+fname
  if "pnfs" in fullname : 
    fullname = "dcap://t3se01.psi.ch:22125//"+fullname
  file = TFile.Open(fullname)
  tree = file.Get("mt2")

  print "-> Skimming " + fullname
  newfilename = fname.split(".root")[0] + "_skim.root"
  newfile = TFile(skimdir+"/"+newfilename, "recreate")
  newtree = tree.CopyTree(cuts)
  newfile.cd()
  newtree.Write()
  newfile.Close()



if __name__ == '__main__':

   import os
   import sys

   from optparse import OptionParser

   parser = OptionParser()
   parser.usage = ""
   parser.add_option("-f","--filter", dest="filter",
                      default="",
                      help="skim only selected dataset")
   parser.add_option("-x","--useXRD", dest="useXRD",
                      default="false",
                      help="useXRD (default=false -> gfal)")
   parser.add_option("-p","--gfalProtocol", dest="gfalProtocol",
                      default="gsiftp",
                      help="gfal protocol (default and recommended: gsiftp; supported alternative: srm)")

   (options,args) = parser.parse_args()
   if len(args)==0:
     print "ERROR! Script needs at least one argument to run:"
     print "   python skimBabies.py [input_directory] [output_directory=\"mt2_200_ht_450\"] [selection=\"mt2>200. && ht>450.\"]"
     exit()

   dir = args[0]
   subdir = "mt2_200_ht_450"
   if len(args)>1 : outdir = args[1]
   cuts = "mt2>200. && ht>450."
   if len(args)>2 : cuts = args[2]

   print "-> Will apply following skim: " + cuts

   if options.filter !="" :
     print "-> Skimming only files containing: " + str(options.filter)
   if options.useXRD == "true":
     print "-> chosen to use xrootd"
   else:
     print "-> chosen to use gfal via "+gfalProtocol+" protocol"


   skimdir = outdir
   os.system("mkdir -p " + skimdir)


   
   if "pnfs/psi.ch" in dir : 
     if options.useXRD == "true":
       status,files = commands.getstatusoutput("xrdfs t3dcachedb.psi.ch ls "+dir)
     else:
       status,files = commands.getstatusoutput("gfal-ls "+options.gfalProtocol+"://t3se01.psi.ch"+dir)
     files=files.splitlines()
   else :
     files = os.listdir(dir)

   for f in files:
     f = f.split("/")[-1] # xrootd takes the full path, truncate
     if ".root" in f:
       if options.filter in f:
         skimBaby(f, dir, skimdir, cuts)

   #print "Find skimmed babies in " + skimdir

