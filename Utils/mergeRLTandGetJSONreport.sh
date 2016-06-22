#!/bin/bash

goldenJson=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt

echo "calculating missing lumis using json: " $goldenJson

if [ $# -ne 1 ]; then
    echo "usage: merge.sh path"
    echo "e.g.: merge.sh /pnfs/psi.ch/cms/trivcat/store/user/casal/babies/MT2_CMGTools-from-CMSSW_7_4_7/prod747data_Run2015B_golden_v2/"
    exit 1
fi

path=$1

prod_label=`basename $path`
mkdir $prod_label

for i in `ls -1 $path`; do
    mkdir $prod_label/$i
    hadd $prod_label/$i/RLTInfo.root `ls $path/$i/*/0000/RLTInfo_* | awk 'BEGIN{ORS=" "}{print "dcap://t3se01.psi.ch:22125/"$1}'`
    heppy_report.py $prod_label/$i -a "" -o lumiSummary_$i.txt # this works if you have 'cmsenv' in a CMSSW release with CMG
    compareJSON.py --sub $goldenJson $prod_label/$i/lumiSummary_$i.txt > $prod_label/$i/missingRuns_$i.txt
done


## to filter out some run ranges
#for x in *.txt ; do  filterJSON.py $x --min=274422 --output=$x.new ; done;

## to calculate the lumi for all files
#for x in *.new; do echo $x; brilcalc lumi -b "STABLE BEAMS" -i $x -u /pb ;done;
