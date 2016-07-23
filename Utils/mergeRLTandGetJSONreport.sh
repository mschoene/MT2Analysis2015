#!/bin/bash

# --- Start of configuration ----------------------------------

#path=/pnfs/psi.ch/cms/trivcat/store/user/mangano/crab/MT2_8_0_11/prodJuly15_runC_all_v1/
path=/pnfs/psi.ch/cms/trivcat/store/user/mangano/crab/MT2_8_0_11/prodJuly15_runC_all_forQCD_v1/
goldenJson=gold_runC.txt

#path=/pnfs/psi.ch/cms/trivcat/store/user/mangano/crab/MT2_8_0_11/prodJuly15_runD_276311-276384_v1
#path=/pnfs/psi.ch/cms/trivcat/store/user/mangano/crab/MT2_8_0_11/prodJuly15_runD_276311-276384_forQCD_v1
#goldenJson=gold_runD.txt


# --- End of configuration ----------------------------------



echo "calculating missing lumis using json: " $goldenJson

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
