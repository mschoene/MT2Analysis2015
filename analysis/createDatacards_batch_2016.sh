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

echo "Loading 80X"
cd /mnt/t3nfs01/data01/shome/mschoene/80X/src/
echo $PWD
eval `scramv1 runtime -sh`

JOBDIR=/scratch/mschoene/datacardCreation_$JOB_ID

echo "Copying all needed stuff..."

mkdir -p $JOBDIR/analysis/
cp -r /mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/signalScansFromDominick $JOBDIR/analysis/
cp -r /mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_dataETH_SnTMC_35p9ifb/ $JOBDIR/analysis/
#cp -r /mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_2016_SnT_36p8_FixedWJets/ $JOBDIR/analysis/
cp -r /mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/cfgs $JOBDIR/analysis/
cp -r /mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/samples/ $JOBDIR

cp /mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/createDatacards_general_zllZinvEst  $JOBDIR/analysis/
#cp /mnt/t3nfs01/data01/shome/mschoene/80X/src/myMT2Analysis/analysis/createDatacards_general $JOBDIR/analysis/

cd $JOBDIR/analysis/
echo $PWD

echo "Starting to create datacards..."
./createDatacards_general_zllZinvEst $1 $2 $3 $4 $5 $6 $7
#./createDatacards_general $1 $2 $3 $4 $5 $6 $7

mkdir /mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/testCopy/

ls $JOBDIR/analysis/EventYields_dataETH_SnTMC_35p9ifb/
ls $JOBDIR/analysis/
ls $JOBDIR/

#if no datacards produced, don't copy anything (to be tested)
if [[ -s $JOBDIR/analysis/EventYields_dataETH_SnTMC_35p9ifb/datacards_$2/datacard* ]]; then
    #if ls $JOBDIR/analysis/EventYields_dataETH_SnTMC_35p9ifb/datacards_$2/*.txt 1>/dev/null 2>&1; then
    #if $JOBDIR/analysis/EventYields_dataETH_SnTMC_35p9ifb/datacards_$2/; then
    
    cd $JOBDIR/analysis/EventYields_dataETH_SnTMC_35p9ifb/datacards_$2/
    tar -czvf $JOBDIR/tared_${M1}_${M11}.tar.gz datacard*${M1}_${M11}*
#fi
    
    #copyOnSE( Form("xrdcp -v %s root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/`whoami`/%s/datacards_%.0f_%.0f/datacard_%s_%s_%.0f_%.0f.txt", newDatacard.c_str(), pathSE.c_str(), mParent, mLSP, binName.c_str(), sigName.c_str(), mParent, mLSP) );
    xrdcp -v  $JOBDIR/tared_${M1}_${M11}.tar.gz  root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/`whoami`/EventYields_${1}/datacards_${2}_${7}/

fi;


cd /scratch/mschoene/
rm -rf $JOBDIR
