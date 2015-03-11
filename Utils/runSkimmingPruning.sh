#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "=== ERROR ==="
    echo "You need to provide one single argument to the script, ie the name of the input configuration file"
    exit;
fi

if [ ! -f $1 ]; then
    echo "File "$1 "not found! Quitting..."
    exit;
fi


rand1=$RANDOM
rand2=$RANDOM
rand3=$RANDOM
rand4=$RANDOM

echo "Reading config file"$1"..." 
configfile=$1
configfile_secured=$rand1.cfg

# check if the file contains something we don't want
if egrep -q -v '^#|^[^ ]*=[^;]*' "$configfile"; then
  echo "Config file is unclean, cleaning it..." 
  # filter the original to a new file
  egrep '^#|^[^ ]*=[^;&]*'  "$configfile" > "$configfile_secured"
  configfile="$configfile_secured"
fi

source $configfile; rm -f $rand1.cfg;
echo "Following configuration will be used:"  
echo "====================================="  
echo "doSkimming"=$doSkimming 
echo "doPruning"=$doPruning 
echo "inputDir"=$inputDir 
echo "inputFilter"=$inputFilter 
echo "outputDir"=$outputDir
echo "skimmingSelection"=$skimmingSelection
echo "branchesToPrune"=$branchesToPrune
echo "====================================="  
echo "" 


if [ -d $outputDir ]; then
    echo "====================================="  
    echo "Output directory "$outputDir" already exists. "
    echo "Be careful that you are not overwriting existing files by mistake";
    echo "====================================="  
    echo "" 
fi


# --- here I run the skimming python code
outputSkimming=/scratch/`whoami`/dirOutSkimming_$rand2
if [ "$doSkimming" = true ]; then
    echo "Running skimming... "
    mkdir $outputSkimming
    python skimBabies.py $inputDir $outputSkimming "$skimmingSelection"  --filter="$inputFilter"
fi



# --- here I run the pruning python code
outputPruning=/scratch/`whoami`/dirOutPruning_$rand3
if [ "$doPruning" = true ]; then
    if [ "$doSkimming" = true ]; then
      inputDir=$outputSkimming
    fi
    python pruneBabies.py $inputDir $outputPruning "$branchesToPrune"  --filter="$inputFilter"
fi



# --- creating destinatinon folder, copying files, cleaning of tmp folders in scratch
if [[ "$outputDir" == *"/pnfs/psi.ch/"* ]]; then
    gfal-mkdir -p srm://t3se01.psi.ch/$outputDir
else
    mkdir -p $outputDir
fi


echo "cleaning/moving temp folders...";
if [[ "$doSkimming" = true && ! "$doPruning" = true ]]; then
    if [[ "$outputDir" == *"/pnfs/psi.ch/"* ]]; then
	for x in $outputSkimming/*; do 
	    gfal-copy file://$x "srm://t3se01.psi.ch"$outputDir/
	done;
    else
	for x in $outputSkimming/*; do 
	    cp $x $outputDir/
	done;
    fi
    rm $outputSkimming/*; rmdir $outputSkimming;
elif [[ ! "$doSkimming" = true &&  "$doPruning" = true ]]; then
    if [[ "$outputDir" == *"/pnfs/psi.ch/"* ]]; then
	for x in $outputPruning/*; do 
	    gfal-copy file://$x "srm://t3se01.psi.ch"$outputDir/
	done;
    else
	for x in $outputPruning/*; do 
	    cp $x $outputDir/
	done;
    fi
    rm $outputPruning/*; rmdir $outputPruning;
elif [[ "$doSkimming" = true &&  "$doPruning" = true ]]; then
    if [[ "$outputDir" == *"/pnfs/psi.ch/"* ]]; then
	for x in $outputSkimming/*; do 
	    gfal-copy file://$x "srm://t3se01.psi.ch"$outputDir/
	done;
	for x in $outputPruning/*; do 
	    gfal-copy file://$x "srm://t3se01.psi.ch"$outputDir/
	done;
    else
	for x in $outputSkimming/*; do 
	    cp $x $outputDir/
	done;
	for x in $outputPruning/*; do 
	    cp $x $outputDir/
	done;
    fi
    rm $outputPruning/*; rmdir $outputPruning;
    rm $outputSkimming/*; rmdir $outputSkimming;
fi

echo ""
echo "Find your files in: " $outputDir
