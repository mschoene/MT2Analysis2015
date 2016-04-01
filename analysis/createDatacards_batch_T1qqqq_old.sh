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
cd /shome/mmasciov/analysisCode_forMerge_replica/analysis
echo $PWD

./createDatacards_T1qqqq $1 $2 $3 $4 $5

