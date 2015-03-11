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


   skimdir = outdir
   os.system("mkdir -p " + skimdir)



   if "pnfs/psi.ch" in dir : 
     status,files = commands.getstatusoutput("gfal-ls "+"srm://t3se01.psi.ch"+dir)
     files=files.splitlines()
   else :
     files = os.listdir(dir)

   for f in files:
     if ".root" in f:
       if options.filter in f:
         skimBaby(f, dir, skimdir, cuts)

   #print "Find skimmed babies in " + skimdir

