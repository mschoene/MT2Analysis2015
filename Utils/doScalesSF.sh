#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/PostProcessed/160415/skimAndPrune"
# --------------------------



# initialization
outputFolder="/scratch/`whoami`/scalesSF/"
#scaleFile="scaleFactors.root"

# here I compile the root macro only once
echo "gROOT->LoadMacro(\"scalesSF.C+\"); gSystem->Exit(0);" |root.exe -b -l ;


while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    name=`echo $line |awk '{print $2}'`
    name=$name"_post_skim_prune"
    outputFile=${outputFolder}/${name}".root";

    
    mkdir -p $outputFolder
    
    echo "scalesSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\")"
    echo "gROOT->LoadMacro(\"scalesSF.C\"); scalesSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\"); gSystem->Exit(0);" |root.exe -b -l ;
    
done < postProcessing.cfg

rm -f scalesSF_C.d scalesSF_C.so;
