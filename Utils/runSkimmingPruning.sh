#!/bin/bash

if [ $CMSSW_BASE ]; then
    myCMSSW=$CMSSW_BASE  
else
    myCMSSW=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw-patch/CMSSW_7_4_12_patch4
fi

gfalProtocol="gsiftp" # if useXRD disabled, use gfal via the given protocol
#gfalProtocol="srm" # alternative to gsiftp (gsiftp supposed to be more stable)

# workaround to use gfal with 80X release
shopt -s expand_aliases
alias semkdir="env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-mkdir -p"
alias secp="env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-copy"


source $VO_CMS_SW_DIR/cmsset_default.sh
#source /swshare/glite/external/etc/profile.d/grid-env.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export LD_LIBRARY_PATH=/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH
echo "Loading your CMSSW release or CMSSW_7_4_12/"
echo "from $myCMSSW"
cd $myCMSSW
eval `scramv1 runtime -sh`
cd -

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
    python skimBabies.py $inputDir $outputSkimming "$skimmingSelection"  --filter="$inputFilter" --useXRD="$useXRD" --gfalProtocol="$gfalProtocol"
fi



# --- here I run the pruning python code
outputPruning=/scratch/`whoami`/dirOutPruning_$rand3
if [ "$doPruning" = true ]; then
    if [ "$doSkimming" = true ]; then
      inputDir=$outputSkimming
    fi
    python pruneBabies.py $inputDir $outputPruning "$branchesToPrune"  --filter="$inputFilter" --useXRD="$useXRD" --gfalProtocol="$gfalProtocol"
fi



# --- creating destinatinon folder, copying files, cleaning of tmp folders in scratch
if [[ "$outputDir" == *"/pnfs/psi.ch/"* ]]; then
    semkdir ${gfalProtocol}://t3se01.psi.ch/$outputDir
else
    mkdir -p $outputDir
fi


echo "cleaning/moving temp folders...";
if [[ "$doSkimming" = true && ! "$doPruning" = true ]]; then
    if [[ "$outputDir" == *"/pnfs/psi.ch/"* ]]; then
	for x in $outputSkimming/*; do 
	    secp file://$x ${gfalProtocol}://t3se01.psi.ch$outputDir/
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
	    secp file://$x ${gfalProtocol}://t3se01.psi.ch$outputDir/
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
	    secp file://$x ${gfalProtocol}://t3se01.psi.ch$outputDir/
	done;
	for x in $outputPruning/*; do 
	    secp file://$x ${gfalProtocol}://t3se01.psi.ch$outputDir/
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
