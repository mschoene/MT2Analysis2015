#!/bin/bash

# This script checks that both RLT and MT2 files corresponding to the same crab job
# have been succesfully copied to the SE. 
# This script does NOT check that these files correspond to ALL the crab jobs originally created/submitted.

# EDIT THIS:
input=/pnfs/psi.ch/cms/trivcat/store/user/mangano/crab/MT2_8_0_11/prodJuly15_runD_276311-276384_forQCD_v1


# ----------
for x in $input/*; 
do 
    prova=$(echo $x |awk -F/ '{print $NF}'|awk -F_ '{print $1}'); 
    ls $x/*/*/mt2*.root| awk -F/ '{print $NF}'| sed 's#mt2_##'|sed 's#.root##'|sort -n > $prova.mt2.txt; 
    ls $x/*/*/RLT*.root| awk -F/ '{print $NF}'| sed 's#RLTInfo_##'|sed 's#.root##'|sort -n > $prova.rlt.txt; 
    echo "=== $prova:"; 
    echo "#mt2, #rlt: " $(cat $prova.mt2.txt|wc -l) " , " $(cat $prova.rlt.txt|wc -l); 
    echo "diff: ";
    diff $prova.mt2.txt $prova.rlt.txt;
    rm $prova.mt2.txt $prova.rlt.txt;
    echo "";
done;