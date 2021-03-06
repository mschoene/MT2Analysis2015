#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/PostProcessed/150415/skimAndPrune"
# --------------------------



# initialization
outputFolder="/scratch/`whoami`/allSF/"
scaleFile="lepSF.root"

# here I compile the root macro only once
echo "gROOT->LoadMacro(\"allSF.C+\"); gSystem->Exit(0);" |root.exe -b -l ;


while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    name=`echo $line |awk '{print $2}'`
    name=$name"_post_skim"
    outputFile=${outputFolder}/${name}".root";

    mkdir -p $outputFolder
    
    echo "allSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$scaleFile\")"
    echo "gROOT->LoadMacro(\"allSF.C\"); allSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$scaleFile\"); gSystem->Exit(0);" |root.exe -b -l ;
    
done < postProcessing.cfg

rm -f allSF_C.d allSF_C.so;
