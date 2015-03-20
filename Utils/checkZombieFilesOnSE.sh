#!/bin/bash


#example:
#./checkZombie.sh /pnfs/psi.ch/cms/trivcat/store/user/mangano/babies/chunks/productionForFrancescoTry5/


if [[ "$#" -eq 0 ]]; then
    echo "Relunch the script as: "
    echo "'checkZombieFilesOnSE.sh inputPnfsFolder'"
    echo "'checkZombieFilesOnSE inputPnfsFolder datasetName'"
    exit;
fi;

dir="$1"
sample=$2

if [[ "$#" -lt 2 ]]; then
    nSamplesToProcess=`ls -1 $dir | wc -l`
    if [ "$nSamplesToProcess" -gt 1 ]; then
	echo "There are " $nSamplesToProcess " to process in the input folder. It may take some time";
	read -r -p "Do you want to process only a subset of them? [y/N] " response
	if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
	then
	    echo "Ok, then rerun the script as 'checkZombieFilesOnSE.sh inputFolder datasetName'"
	    echo "where 'dataset' is one of the following entries: "
	    ls -1 $dir;
	    exit;
	fi
    fi;
fi;



echo "gROOT->LoadMacro(\"checkZombie.C+\"); gSystem->Exit(0);" 
echo "gROOT->LoadMacro(\"checkZombie.C+\"); gSystem->Exit(0);" |root.exe -b -l ;


for x in `ls -1 $dir `; do 
    if [[ (-n "$sample") && ("$x" != "$sample") ]]; then continue;fi; 
    echo "checking sample: " $x;    
    for y in `ls $dir/$x`; do
	echo "check " $dir/$x/$y
	echo "gROOT->LoadMacro(\"checkZombie.C+\"); checkZombie(\"$dir/$x/$y\"); gSystem->Exit(0);"| root.exe -b -l; 
    done;    
done;

rm checkZombie_C.d checkZombie_C.so;