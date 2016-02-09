#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/babies/MT2_CMGTools-from-CMSSW_7_4_12/fullMC_Spring15_miniAODv2_03Nov2015/"
#"/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/babies/MT2_CMGTools-from-CMSSW_7_4_12/fullMC_miniAODv2_06Nov2015/"
#"/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/babies/MT2_CMGTools-from-CMSSW_7_4_12/fullMC_miniAODv2_09Nov2015/"


#"/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/babies/MT2_CMGTools-from-CMSSW_7_4_12/ZG_20Dec2015/"
#"/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/babies/MT2_CMGTools-from-CMSSW_7_4_12/fullData_miniAODv2_15Nov2015_jecV6/"
productionName="09Feb2016_fullMC_mAODv2_WithSF"
fileExt="_post.root"
isCrab=1
inputPU="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/74X/Spring15/PostProcessed/23Oct2015_data_noSkim/JetHT_Run2015D_post.root"
PUvar="nVert"
#GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt"
#GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-257599_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
#GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
#GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
#GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260426_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
SilverJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt"
applyJSON=0
doSilver=0
doFilterTxt=0
doAllSF=1
# --------------------------


# initialization
jobsLogsFolder="./${productionName}"
outputFolder="/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/MT2production/74X/Spring15/PostProcessed/"$productionName"/"
workingFolder="/scratch/`whoami`/"$productionName


if [[ "$#" -eq 0 ]]; then
    echo "Relunch the script with one of the following options: "
    echo "./doTreeProduction.sh post      # post-processing"
    echo "./doTreeProduction.sh postCheck # check post-processing"
    echo "./doTreeProduction.sh mergeData # merge data and remove duplicates (not implemented yet)"
    echo "./doTreeProduction.sh addAllSF  # add all scale factor weights"
    echo "./doTreeProduction.sh addISR    # add isr weights variables(not fully implemented yet for the standalone version)"
    echo "./doTreeProduction.sh addBtag   # add b-tagg weights variables"
    echo "./doTreeProduction.sh addLepSF  # add lepton weights variables"
    echo "./doTreeProduction.sh clean     # clean (not implemented yet)"
    exit;
fi;


 

#######################################
if [[ "$1" = "post" ]]; then

# --- check the existence of outputFolder on SE ---
gfal-ls srm://t3se01.psi.ch$outputFolder &> /tmp/checkOutputDir
if [ -n "`cat /tmp/checkOutputDir|grep 'No such file or directory'`"  ]; then
    :
else
    echo "WARNING: output directory " $outputFolder "already exists."
    read -r -p "Are you sure you want to write there ? [y/N] " response
    if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
    then
	echo "Ok. Will proceed then ... "
    else
	echo "Exiting ..."
	exit;
    fi
fi
# --- 

if [ -d "$jobsLogsFolder" ]; then 
    echo "ERROR: the logFolder" $jobsLogsFolder " already exists."
    echo "Delete it and start from a clean area, or redirect the logs in a different place."
    echo "Exiting ..."
    exit
else
    mkdir  $jobsLogsFolder
fi

python $PWD/convertGoodRunsList_JSON.py $GoldenJSON >& goodruns_golden.txt
python $PWD/convertGoodRunsList_JSON.py $SilverJSON >& goodruns_silver.txt
gfal-mkdir -p srm://t3se01.psi.ch/$outputFolder 
gfal-copy file://$GoldenJSON srm://t3se01.psi.ch/$outputFolder/ 
gfal-copy file://$SilverJSON srm://t3se01.psi.ch/$outputFolder/ 

echo "Location of log files is: " $jobsLogsFolder
echo "Location of final files on SE is: " $outputFolder
echo "Working folder on working-node is: " $workingFolder


if [ $CMSSW_BASE ]; then
    myCMSSW=$CMSSW_BASE  
else
    myCMSSW=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw-patch/CMSSW_7_4_12_patch4
fi

#get the pileUp histogram
if [ ! -f MyDataPileupHistogram.root ]; then
    pileupCalc.py -i $GoldenJSON --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 80000 --maxPileupBin 50 --numPileupBins 50  MyDataPileupHistogram.root
fi

### here I compile the root macro only once
### Uncomment for ROOT v5
#echo "gROOT->LoadMacro(\"goodrun.cc+\"); gSystem->Exit(0);" |root.exe -b -l ;
###echo "gROOT->LoadMacro(\"postProcessing.C+\"); gSystem->Exit(0);" |root.exe -b -l ;

while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    id=`echo $line |awk '{print $1}'`
    name=`echo $line |awk '{print $2}'`
    xsec=`echo $line |awk '{print $3}'`
    filter=`echo $line |awk '{print $4}'`
    kfactor=`echo $line |awk '{print $5}'`

    doPruning="true"
    if [ $id -lt 10 ]; then
	doPruning="false"
    fi;

    if [ $id -lt 10 ]; then
	doAllSF="false"
    fi;


    if [[ $id -lt 10 && $doFilterTxt == 1 ]]; then
	# this is the merged (and duplicate removed) file including Nov15 txts and Oct 15 txts for all datasets
	#filterTxt=/shome/casal/eventlist_Nov14/filter_cscNov15_ecalscnNov15_cscOct15_sortu.txt
	# latest one from december
	filterTxt=/shome/casal/eventlist_Dec4/allfilters_Dec01.txt
	outputFilteredFile=${workingFolder}/${name}_filtered$fileExt;
    fi;
    
    outputFile=${workingFolder}/${name}$fileExt;

    if [ ${isCrab} = 1 ]; then
	crabExt=$(ls $inputFolder/$name/)
    else
	crabExt=""
    fi;

    cat <<EOF > batchScript_${name}.sh
#!/bin/bash


#### The following configurations you should not need to change
# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N postProcessing_${name}_`whoami`

### Specify the queue on which to run
#$ -q short.q

# Change to the current working directory from which the job got
# submitted . This will also result in the job report stdout/stderr being
# written to this directory, if you do not override it (below).
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
#$ -o $jobsLogsFolder/${name}.out
#$ -e $jobsLogsFolder/${name}.err

#source /swshare/psit3/etc/profile.d/cms_ui_env.sh
#source $VO_CMS_SW_DIR/cmsset_default.sh
#source /swshare/ROOT/thisroot.sh 
#eval \`scramv1 runtime -sh\`



source $VO_CMS_SW_DIR/cmsset_default.sh
source /swshare/glite/external/etc/profile.d/grid-env.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export LD_LIBRARY_PATH=/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH
echo "Loading your CMSSW release or CMSSW_7_4_12/"
echo "from $myCMSSW"
cd $myCMSSW
eval `scramv1 runtime -sh`
cd -


mkdir -p $workingFolder
gfal-mkdir -p srm://t3se01.psi.ch/$outputFolder

echo "postProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\", \"$inputPU\", \"$PUvar\", $applyJSON);"

### Uncomment for ROOT v5
#echo "gSystem->Load(\"goodrun_cc\"); gROOT->LoadMacro(\"postProcessing.C\"); postProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\",\"$inputPU\",\"$PUvar\",$applyJSON); gSystem->Exit(0);" |root.exe -b -l ;
### Comment for ROOT v5

echo "gROOT->LoadMacro(\"postProcessing.C+\"); postProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\",\"$inputPU\",\"$PUvar\",$applyJSON,$doAllSF,$doSilver); gSystem->Exit(0);" |root.exe -b -l ;

if [[ $id -lt 10 && $doFilterTxt == 1 ]]; then 
   echo "filterFromTxt(\"$filterTxt\",\"$outputFile\",1,\"$outputFilteredFile\",\"mt2\",\"\")"
   echo "gROOT->LoadMacro(\"filterFromTxt.C\"); filterFromTxt(\"$filterTxt\",\"$outputFile\",1,\"$outputFilteredFile\",\"mt2\",\"\"); gSystem->Exit(0);" |root.exe -b -l ;

   mv $outputFilteredFile $outputFile;
fi;


#mv $outputFile $outputFolder
gfal-copy file://$outputFile srm://t3se01.psi.ch/$outputFolder
rm $outputFile

skimmingPruningCfg="${workingFolder}/skimmingPruning_${name}.cfg"
    cat skimmingPruning.cfg |grep -v \# | sed  "s#INPUTDIR#${outputFolder}#" |sed "s#INPUTFILTER#${name}#" \
	| sed "s#OUTPUTDIR#${outputFolder}/skimAndPrune#" | sed "s#DOPRUNING#${doPruning}#" > \$skimmingPruningCfg

./runSkimmingPruning.sh \$skimmingPruningCfg
rm \$skimmingPruningCfg

#echo "is anything left in working folder? workingFolder: " 
#ls $workingFolder


EOF

qsub -q long.q batchScript_${name}.sh;
rm batchScript_${name}.sh;

done < postProcessing.cfg

###rm -f postProcessing_C.d postProcessing_C.so;

fi;

if [[ "$1" = "postCheck" ]]; then
    sizeLogsErr=`ls -l ${jobsLogsFolder}/*.err |awk '{sum+=$5} END {print sum}'`
    if [[ $sizeLogsErr -eq 0 ]]; then
	echo "there were no errors. Zipping all logs and copying them to the SE"	 
	cd $jobsLogsFolder
	tar -czvf logs.tgz  *
	gfal-copy file://`pwd`/logs.tgz srm://t3se01.psi.ch/$outputFolder
	cd ..
	rm $jobsLogsFolder/*
	rmdir $jobsLogsFolder
	rm postProcessing_C.d;
	rm postProcessing_C.so;
    else
	echo "ERROR: something went wrong. Check your logs in " $jobsLogsFolder
    fi

fi

if [[ "$1" = "mergeData" ]]; then
    echo "need to implement the data merging + duplicate removal"
fi

if [[ "$1" = "addLepSF" ]]; then
#the outputfolder of postprocessing is the input file for the adding of the scale factors
    ./doLepSF.sh $outputFolder 
fi

if [[ "$1" = "addAllSF" ]]; then
#the outputfolder of postprocessing is the input file for the adding of the scale factors
    ./doAllSF.sh $outputFolder 
fi

if [[ "$1" = "addBtag" ]]; then
    ./doBTagSF.sh
fi

if [[ "$1" = "addISR" ]]; then
    ./doISRSF.sh
fi

if [[ "$1" = "clean" ]]; then
    echo "INFO: option 'clean' is not implemented yet."

fi