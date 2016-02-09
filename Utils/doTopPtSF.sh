#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
#inputFolder="/pnfs/psi.ch/cms/trivcat/store/user//mmasciov/MT2production/74X/Spring15_forCommissioningPAS/testBatch/TTJets_LO_50n/"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/74X/Spring15/PostProcessed/18Jan2016_jecV6/skimAndPrune/"

# --------------------------



# initialization
outputFolder="/scratch/`whoami`/sfTests/"
objectName="GenPart"

# here I compile the root macro only once
echo "gROOT->LoadMacro(\"topPtSF.C+\"); gSystem->Exit(0);" |root.exe -b -l ;


while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    name=`echo $line |awk '{print $2}'`
    name=$name"_post_skim"
    outputFile=${outputFolder}/${name}".root";

    mkdir -p $outputFolder
    
    echo "topPtSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objectName\")"
    echo "gROOT->LoadMacro(\"topPtSF.C\"); topPtSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objectName\"); gSystem->Exit(0);" |root.exe -b -l ;
    
done < postProcessing.cfg

rm -f topPtSF_C.d topPtSF_C.so;
