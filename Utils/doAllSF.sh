#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/74X/Spring15/PostProcessed/18Jan2016_jecV6/skimAndPrune/"

# --------------------------

# initialization
outputFolder="/scratch/`whoami`/allSF/"
objectName=""

args=("$@")
if [ ${#args[@]}>0 ]; then
    inputFolder=$1
fi

# here I compile the root macro only once
echo "gROOT->LoadMacro(\"allSF.C+\"); gSystem->Exit(0);" |root.exe -b -l ;

while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    name=`echo $line |awk '{print $2}'`
    name=$name"_post"
    outputFile=${outputFolder}/${name}".root";

    mkdir -p $outputFolder
    
    echo "allSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objectName\")"
    echo "gROOT->LoadMacro(\"allSF.C\"); allSF(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",\"$objectName\"); gSystem->Exit(0);" |root.exe -b -l ;
    
done < postProcessing.cfg

rm -f allSF_C.d allSF_C.so allSF_C_ACLiC_dict_rdict.pcm;
