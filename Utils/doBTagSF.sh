#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/PostProcessed/150415/skimAndPrune"
# --------------------------



# initialization
outputFolder="/scratch/`whoami`/BTagSF/"
objName="jet"
#scaleFile="scaleFactors.root"

# here I compile the root macro only once
echo "gROOT->LoadMacro(\"btagSF.C+\"); gSystem->Exit(0);" |root.exe -b -l ;


while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    name=`echo $line |awk '{print $2}'`
    name=$name"_post_skim_prune"
    outputFile=${outputFolder}/${name}".root";

    
    mkdir -p $outputFolder
    
    echo "btagSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objName\")"
    echo "gROOT->LoadMacro(\"btagSF.C\"); btagSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objName\"); gSystem->Exit(0);" |root.exe -b -l ;
    
done < postProcessing.cfg

rm -f btagSF_C.d btagSF_C.so;
