#!/bin/bash                                                                                                                                                                          
echo $#;
if [ $# != 5 ]; then
    echo "USAGE: ${0} CFG M1 M2 M11 M22";
    exit;
fi

CFG=$1
M1=$2
M2=$3
M11=$4
M22=$5

source $VO_CMS_SW_DIR/cmsset_default.sh
source /swshare/glite/external/etc/profile.d/grid-env.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export LD_LIBRARY_PATH=/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH

echo "Loading CMSSW_7_4_7"
cd /shome/mmasciov/CMSSW_7_4_7_MT2/src/
echo $PWD
eval `scramv1 runtime -sh`
eval `cmsenv`

DIRECTORY=/scratch/mmasciov/analysisCode_forMerge_replica_T2bb/
JOBDIR=/scratch/mmasciov/datacardCreation_$JOB_ID

###rm -rf /scratch/mmasciov/*

###if [ -d "$DIRECTORY" ]; then
###
###    echo "Directory already exists. Not copying again"
###   
###else
###    
###    mkdir -p /scratch/mmasciov/analysisCode_forMerge_replica_T2bb/analysis/
###    cp -r /shome/mmasciov/analysisCode_forMerge_replica/analysis/signalScansFromDominick /scratch/mmasciov/analysisCode_forMerge_replica_T2bb/analysis/
###    cp -r /shome/mmasciov/analysisCode_forMerge_replica/analysis/EventYields_data_Run2015D_25nsGolden_miniAODv2_fullStat_final_2p155_T2bb /scratch/mmasciov/analysisCode_forMerge_replica_T2bb/analysis/
###    cp -r /shome/mmasciov/analysisCode_forMerge_replica/analysis/cfgs /scratch/mmasciov/analysisCode_forMerge_replica_T2bb/analysis/
###    cp -r /shome/mmasciov/analysisCode_forMerge_replica/samples/ /scratch/mmasciov/analysisCode_forMerge_replica_T2bb/
###    
###    cp /shome/mmasciov/analysisCode_forMerge_replica/analysis/createDatacards_T2bb /scratch/mmasciov/analysisCode_forMerge_replica_T2bb/analysis/
###    
###fi

echo "Copying all needed stuff..."

mkdir -p $JOBDIR/analysis/
cp -r /shome/mmasciov/analysisCode_forMerge_replica/analysis/signalScansFromDominick $JOBDIR/analysis/
cp -r /shome/mmasciov/analysisCode_forMerge_replica/analysis/EventYields_data_Run2015_25nsGolden_2p3ifb_T2bb $JOBDIR/analysis/
cp -r /shome/mmasciov/analysisCode_forMerge_replica/analysis/cfgs $JOBDIR/analysis/
cp -r /shome/mmasciov/analysisCode_forMerge_replica/samples/ $JOBDIR

cp /shome/mmasciov/analysisCode_forMerge_replica/analysis/createDatacards_T2bb $JOBDIR/analysis/

cd $JOBDIR/analysis/
echo $PWD

echo "Starting to create datacards..."
./createDatacards_T2bb $1 $2 $3 $4 $5

cd /scratch/mmasciov/
rm -rf $JOBDIR