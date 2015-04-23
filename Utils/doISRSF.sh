#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/PostProcessed/160415/skimAndPrune"
# --------------------------



# initialization
outputFolder="/scratch/`whoami`/isrSF_2204/"
objName="GenPart"

# here I compile the root macro only once
echo "gROOT->LoadMacro(\"isrSF.C+\"); gSystem->Exit(0);" |root.exe -b -l ;


while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    name=`echo $line |awk '{print $2}'`
    name=$name"_post_skim"
    outputFile=${outputFolder}/${name}".root";

    
    mkdir -p $outputFolder
    
    echo "isrSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objName\")"
    echo "gROOT->LoadMacro(\"isrSF.C\"); isrSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objName\"); gSystem->Exit(0);" |root.exe -b -l ;
    
done < postProcessing.cfg

rm -f isrSF_C.d isrSF_C.so;
