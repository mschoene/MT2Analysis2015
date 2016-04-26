#!/bin/bash

echo $#;
if [ $# -le 4 ]; then
    echo "USAGE: ${0} cfg mc sampleID ijob Njobs ";
    exit;
fi

WDIR=$SGE_O_WORKDIR

CFG=$1
MC=$2
SID=$3
JOB=$4
NJOB=$5
TOSE=$6

source $VO_CMS_SW_DIR/cmsset_default.sh
cd $WDIR
eval `scramv1 runtime -sh`
#cmsenv

./doRebalancing $CFG $MC $SID $JOB $NJOB $TOSE