#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/PostProcessed/150415/skimAndPrune"
# --------------------------



# initialization
outputFolder="/scratch/`whoami`/Scaled_1504/"
objName="lep"
scaleFile="scaleFactors.root"

# here I compile the root macro only once
echo "gROOT->LoadMacro(\"scaleFactors.C+\"); gSystem->Exit(0);" |root.exe -b -l ;


while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    name=`echo $line |awk '{print $2}'`
    name=$name"_post_skim_prune"
    outputFile=${outputFolder}/${name}".root";

    
    mkdir -p $outputFolder
    
    echo "scaleFactors(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objName\",\"$scaleFile\")"
    echo "gROOT->LoadMacro(\"scaleFactors.C\"); scaleFactors(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objName\",\"$scaleFile\"); gSystem->Exit(0);" |root.exe -b -l ;
    
done < postProcessing.cfg

rm -f scaleFactors_C.d scaleFactors_C.so;
