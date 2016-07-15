#!/bin/bash                                                                                                                                                                          
echo $#;
if [ $# -lt 6 ]; then
    echo "USAGE: ${0} MODEL CFG M1 M2 M11 M22 [LABEL]";
    exit;
fi

MODEL=$1
CFG=$2
M1=$3
M2=$4
M11=$5
M22=$6
if [ $# -ge 7 ]; then
    LABEL=$7;
else
    LABEL="";
fi

source $VO_CMS_SW_DIR/cmsset_default.sh
#source /mnt/t3nfs01/data01/swshare/glite/external/etc/profile.d/grid-env.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export LD_LIBRARY_PATH=/mnt/t3nfs01/data01/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH

echo "Loading CMSSW_7_4_7"
cd /mnt/t3nfs01/data01/shome/casal/CMSSW_7_4_12_patch4_MT2/src/
echo $PWD
eval `scramv1 runtime -sh`

JOBDIR=/scratch/casal/datacardCreation_$JOB_ID

echo "Copying all needed stuff..."

mkdir -p $JOBDIR/analysis/
cp -r /mnt/t3nfs01/data01/shome/casal/CMSSW_7_4_12_patch4_MT2/src/analysisJul/analysis/signalScansFromDominick $JOBDIR/analysis/
cp -r /mnt/t3nfs01/data01/shome/casal/CMSSW_7_4_12_patch4_MT2/src/analysisJul/analysis/EventYields_data_Run2016_7p7ifb $JOBDIR/analysis/
cp -r /mnt/t3nfs01/data01/shome/casal/CMSSW_7_4_12_patch4_MT2/src/analysisJul/analysis/cfgs $JOBDIR/analysis/
cp -r /mnt/t3nfs01/data01/shome/casal/CMSSW_7_4_12_patch4_MT2/src/analysisJul/samples/ $JOBDIR

cp /mnt/t3nfs01/data01/shome/casal/CMSSW_7_4_12_patch4_MT2/src/analysisJul/analysis/createDatacards_general $JOBDIR/analysis/

cd $JOBDIR/analysis/
echo $PWD

echo "Starting to create datacards..."
./createDatacards_general $1 $2 $3 $4 $5 $6 $7

cd /scratch/casal/
rm -rf $JOBDIR