#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
#inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/casal/babies/MT2_CMGTools-from-CMSSW_7_4_12/data_Run2015D_round2/"
#inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/babies/MT2_CMGTools-from-CMSSW_7_4_7/prod747mc_Spring15/"
#inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/babies/MT2_CMGTools-from-CMSSW_7_4_12/MC_forZGratio_05Oct2015/"
#inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/74X/Spring15/22Sep2015_singleTop/"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/babies/MT2_CMGTools-from-CMSSW_7_4_12/fullData_miniAODv2_04Nov2015/"
productionName="05Nov2015_dataRunD_goldenJSON_1p25ifb_dataSkim"
#productionName="05Oct2015_forMonojetBins"
fileExt="_post.root"
isCrab=1
inputPU="/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/74X/firstData2015/PostProcessed/04Sep2015_GoldenJSON_v2/JetHT_Run2015C_post.root"
PUvar="nVert"
#GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt"
#GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-257599_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
#GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
#GoldenJSON="/afs/cern.ch/user/g/gzevi/public/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON_Photon.txt"
applyJSON=1
doFilterTxt=1
# --------------------------


# initialization
jobsLogsFolder="./${productionName}"
outputFolder="/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/MT2production/74X/Spring15/PostProcessed/"$productionName"/"
workingFolder="/scratch/`whoami`/"$productionName


if [[ "$#" -eq 0 ]]; then
    echo "Relunch the script with one of the following options: "
    echo "./doTreeProduction.sh post"
    echo "./doTreeProduction.sh postCheck"
    echo "./doTreeProduction.sh clean"
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

python $PWD/convertGoodRunsList_JSON.py $GoldenJSON >& goodruns.txt
gfal-mkdir -p srm://t3se01.psi.ch/$outputFolder 
gfal-copy file://$GoldenJSON srm://t3se01.psi.ch/$outputFolder/ 

echo "Location of log files is: " $jobsLogsFolder
echo "Location of final files on SE is: " $outputFolder
echo "Working folder on working-node is: " $workingFolder

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

    if [[ $id -lt 10 && $doFilterTxt == 1 ]]; then	
	filterTxt=/shome/casal/eventlist_csc2015/eventlist_`echo $name | cut -d _ -f 1`_csc2015.txt
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

source /swshare/psit3/etc/profile.d/cms_ui_env.sh
source $VO_CMS_SW_DIR/cmsset_default.sh
source /swshare/ROOT/thisroot.sh 
eval \`scramv1 runtime -sh\`

mkdir -p $workingFolder
gfal-mkdir -p srm://t3se01.psi.ch/$outputFolder

echo "postProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\", \"$inputPU\", \"$PUvar\", $applyJSON);"
### Uncomment for ROOT v5
#echo "gSystem->Load(\"goodrun_cc\"); gROOT->LoadMacro(\"postProcessing.C\"); postProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\",\"$inputPU\",\"$PUvar\",$applyJSON); gSystem->Exit(0);" |root.exe -b -l ;
### Comment for ROOT v5
echo "gROOT->LoadMacro(\"postProcessing.C\"); postProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\",\"$inputPU\",\"$PUvar\",$applyJSON); gSystem->Exit(0);" |root.exe -b -l ;

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



if [[ "$1" = "clean" ]]; then
    echo "INFO: option 'clean' is not implemented yet."

fi