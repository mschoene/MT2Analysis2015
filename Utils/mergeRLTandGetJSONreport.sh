#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: merge.sh path"
    echo "e.g.: merge.sh /pnfs/psi.ch/cms/trivcat/store/user/casal/babies/MT2_CMGTools-from-CMSSW_7_4_7/prod747data_Run2015B_golden_v2/"
    exit 1
fi

path=$1


for i in `ls -1 $path`; do
    mkdir $i
    hadd $i/RLTInfo.root `ls $path/$i/*/0000/RLTInfo_* | awk 'BEGIN{ORS=" "}{print "dcap://t3se01.psi.ch:22125/"$1}'`
    heppy_report.py $i -a "" -o lumiSummary_$i.txt # this works if you have 'cmsenv' in a CMSSW release with CMG
done
