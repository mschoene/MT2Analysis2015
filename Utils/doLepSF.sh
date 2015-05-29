#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/PostProcessed/150415/skimAndPrune"
# --------------------------



# initialization
outputFolder="/scratch/`whoami`/lepSF/"
scaleFile="lepSF.root"

# here I compile the root macro only once
echo "gROOT->LoadMacro(\"lepSF.C+\"); gSystem->Exit(0);" |root.exe -b -l ;


while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    name=`echo $line |awk '{print $2}'`
    name=$name"_post_skim"
    outputFile=${outputFolder}/${name}".root";

    mkdir -p $outputFolder
    
    echo "lepSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$scaleFile\")"
    echo "gROOT->LoadMacro(\"lepSF.C\"); lepSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$scaleFile\"); gSystem->Exit(0);" |root.exe -b -l ;
    
done < postProcessing.cfg

rm -f lepSF_C.d lepSF_C.so;
