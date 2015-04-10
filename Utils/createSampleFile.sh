#!/bin/bash

#echo $#
if [ $# != '2' ]
then
    echo "Usage: ./createSamplesFile.sh inputDir outputFile"
    echo "e.g.: ./createSamplesFile.sh /pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/PostProcessed/080415/skimAndPrune/ PHYS14_v4"
    echo "This will create the file ../samples/samples_PHYS14_v4_skimprune.dat"
    exit
fi

inputDir=$1
outputFile=`pwd`/../samples/samples_$2_skimprune.dat

for i in $(ls $inputDir/*prune.root)
do
    echo $i >> $outputFile
done

sed -i 's/\/pnfs\//srm:\/\/t3se01.psi.ch:8443\/srm\/managerv2?SFN=\/pnfs\//g' $outputFile