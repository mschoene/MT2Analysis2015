#!bin/bash

### HOW TO USE:
### sh checkBash.sh 'dataset_name'
### for instance: sh checkBash.sh 'TTJets', after replacing workDir and remoteDir inside the script

### local output dir of heppy
workDir=/afs/cern.ch/user/m/mmasciov/CMSSW_7_2_3_PhysicsTools/src/CMGTools/TTHAnalysis/cfg/06Feb2015/

### remote output dir of heppy
remoteDir=/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/PHYS14_production/06Feb2015/

### dataset name
sample=$1

for i in $(ls $workDir); 
do 
    s="ZJetsToNuNu_HT600toInf"
    if [[ $i == $sample* ]]
    then 
	n=${i#*k} 
	Y=0 
	for m in $(ls $remoteDir$s/mt2_*) 
	do 
	    m=${m#*mt2_} 
	    m=${m%.root} 
	    if [[ $m == $n ]] 
	    then 
		Y=1 
	    fi 
	done 
	if [[ $Y != 1 ]] 
	then 
	    echo $n 
	fi
    fi
done