#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="mt2"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mangano/crab/MT2_8_0_5/data05June_all_v1"
listOfSamplesFile="postProcessing2016-Data.cfg"
#listOfSamplesFile="postProcessing2016-MC.cfg"

productionName="$(basename $inputFolder)" 


#outputFolder="/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/ZGproduction/80X/PostProcessed/"$productionName"_pu/"
outputFolder="/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/MT2production/80X/PostProcessed/"$productionName"/"
#outputFolder="/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/babies/80X/MT2/PostProcessed/"$productionName"/"

# in current implementation one also needs to change this in runSkimmingPruning.sh
useXRD="false"
gfalProtocol="gsiftp" # if useXRD disabled, use gfal via the given protocol
#gfalProtocol="srm" # alternative to gsiftp (gsiftp supposed to be more stable)

fileExt="_post.root"
isCrab=1
inputPU="MyDataPileupHistogram.root"
PUvar="nTrueInt"
GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt"
SilverJSON=$GoldenJSON
applyJSON=1     #0 for MC
doSilver=0      #0 for MC
doFilterTxt=0   #0 for MC
doAllSF=0       #1 for MC
doPreProc=0     #0 (only 1 for TTJets or if you want to split MC samples, then run ./doTreeProduction pre first)

# ------------------------------------------------------------------------------
# NOTA BENE: the following files should also be (may have to) edited according 
# to the type of production (for example MT2 vs ZGamma, data vs MC):
#
# postProcessing2016-MC.cfg/postProcessing2016-Data.cfg
# skimmingPruning.cfg, skimmingPruningQCD.cfg, skimmingPruningMonoJet.cfg
#
# TO DO: add switches here in doTreeProduction.sh instead of force editing by hand
# ------------------------------------------------------------------------------








# initialization
jobsLogsFolder="./${productionName}"
workingFolder="/scratch/`whoami`/"$productionName



if [[ "$#" -eq 0 ]]; then
    echo "Relaunch the script with one of the following options: "
    echo "./doTreeProduction.sh pre       # pre-processing"
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


##### beginning pre processing
if [[ "$1" = "pre" ]]; then
    if [ -d "$jobsLogsFolder" ]; then 
	echo "ERROR: the logFolder" $jobsLogsFolder " already exists."
	echo "Delete it and start from a clean area, or redirect the logs in a different place."
	echo "Exiting ..."
	exit
    else
	mkdir  $jobsLogsFolder
    fi


    while read line; 
    do 
	case "$line" in \#*) continue ;; esac; #skip commented lines
	case "$line" in *"_ext"*) continue ;; esac; #skip extensions
	if [ -z "$line" ]; then continue; fi;  #skip empty lines
	id=`echo $line |awk '{print $1}'`
	name=`echo $line |awk '{print $2}'`

	fileList=inputChunkList.txt
	
	if [ -e inputChunkList.txt ]; then
	    echo "deleting the old file list"
	    rm $fileList
	fi;

	crabExt=""
	if [ ${isCrab} = 1 ]; then
	    crabExt=$(ls $inputFolder/$name/)
	fi;

	if [ ${isCrab} = 1 ]; then
	    for ((i=0; i<10; i++)); do
		if [ -d $inputFolder/${name}/$crabExt/000${i}/ ]; then
		    for f in $inputFolder/$name/$crabExt/000${i}/mt2*.root; do
			echo $f>>$fileList
		    done;
		fi;
	    done;
	else
	    for f in $inputFolder/$name/mt2*.root; do
		echo $f>>$fileList
	    done;
	fi;
	

	numFiles=$(wc -l inputChunkList.txt | awk '{print $1}')
	echo "number of files = " $numFiles
	
	if [ -d $inputFolder/${name}_ext/ ]; then
	    if [ ${isCrab} = 1 ]; then
    		crabExt=$(ls $inputFolder/${name}_ext/)
		for ((i=0; i<10; i++)); do
		    if [ -d $inputFolder/${name}_ext/$crabExt/000${i}/ ]; then
			for f in $inputFolder/${name}_ext/${crabExt}/000${i}/mt2*.root; do
			    echo $f>>$fileList
			done;
		    fi;
		done;
	    else
		for f in $inputFolder/${name}_ext_${i}/mt2*.root; do
		    echo $f>>$fileList
		done;
	    fi;
	fi;
	    

	numFiles=$(wc -l inputChunkList.txt | awk '{print $1}')
	echo "number of files = " $numFiles

	awk '$0="dcap://t3se01.psi.ch:22125/"$0' $fileList > inputChunkList_dcap.txt 
	mv inputChunkList_dcap.txt inputChunkList.txt
	cp  inputChunkList.txt  inputChunkList_${name}.txt
	fileListTemp=inputChunkList_${name}.txt
    #   fileList=inputChunkList_dcap.txt
	
	scriptName=batchScript_${name}.sh

	outputFile=${name}_pre.cfg

	cat <<EOF > $scriptName
#!/bin/bash

#### The following configurations you should not need to change
# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N pP_${name}_`whoami`

### Specify the queue on which to run
#$ -q short.q

# Change to the current working directory from which the job go
# submitted . This will also result in the job report stdout/stderr being
# written to this directory, if you do not override it (below).
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
####v$ -o $jobsLogsFolder/${name}.out
####c$ -e $jobsLogsFolder/${name}.err
#$ -o $jobsLogsFolder/${name}.out
#$ -e $jobsLogsFolder/${name}.err


source $VO_CMS_SW_DIR/cmsset_default.sh
source /mnt/t3nfs01/data01/swshare/glite/external/etc/profile.d/grid-env.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export LD_LIBRARY_PATH=/mnt/t3nfs01/data01/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH
echo "Loading your CMSSW release or CMSSW_7_4_12/"
echo "from $myCMSSW"
cd $myCMSSW
eval `scramv1 runtime -sh`
cd -
echo "preProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$id,$fileList);"
echo "gROOT->LoadMacro(\"preProcessing.C+\"); preProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$id,\"$fileListTemp\"); gSystem->Exit(0);" |root.exe -b -l ;


EOF

	qsub -q long.q $scriptName;
	rm $scriptName;
	
    done < $listOfSamplesFile

fi;
##### done pre processing






#######################################
if [[ "$1" = "post" ]]; then

if [[ $useXRD == "true" ]]; then
    xrdfs t3dcachedb.psi.ch ls $outputFolder &> ./checkOutputDir
else
    gfal-ls ${gfalProtocol}://t3se01.psi.ch$outputFolder &> ./checkOutputDir
fi
# --- check the existence of outputFolder on SE ---
if [ -n "`cat ./checkOutputDir|grep 'No such file or directory'`"  ]; then
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


if [[ $useXRD == "true" ]]; then
    #xrdfs t3dcachedb.psi.ch mkdir -p $outputFolder # recursive mkdir does not work via xrootd
    gfal-mkdir -p ${gfalProtocol}://t3se01.psi.ch/$outputFolder # use gfal instead
else
    gfal-mkdir -p ${gfalProtocol}://t3se01.psi.ch/$outputFolder
fi

python $PWD/convertGoodRunsList_JSON.py $GoldenJSON >& goodruns_golden.txt
if [[ $useXRD == "true" ]]; then
    xrdcp -d 1 $GoldenJSON root://t3dcachedb.psi.ch:1094/$outputFolder/ 
else
    gfal-copy file://$GoldenJSON ${gfalProtocol}://t3se01.psi.ch/$outputFolder/
fi


if [ $doSilver -eq 1 ]; then
    if [[ $useXRD == "true" ]]; then
	xrdcp -d 1 $SilverJSON root://t3dcachedb.psi.ch:1094/$outputFolder/ 
    else
	gfal-copy file://$SilverJSON ${gfalProtocol}://t3se01.psi.ch/$outputFolder/
    fi
    python $PWD/convertGoodRunsList_JSON.py $SilverJSON >& goodruns_silver.txt
fi

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
    pileupCalc.py -i $GoldenJSON --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50  MyDataPileupHistogram.root
fi


echo "gROOT->LoadMacro(\"goodrunClass.cc+\"); gSystem->Exit(0);" |root.exe -b -l ;
echo "gSystem->Load(\"goodrunClass_cc.so\"); gROOT->LoadMacro(\"postProcessing.C+\"); gSystem->Exit(0);" |root.exe -b -l ;

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
	# filterTxt=/shome/casal/eventlist_Nov14/filter_cscNov15_ecalscnNov15_cscOct15_sortu.txt
	# latest one from december
	filterTxt=/shome/mmasciov/CMSSW_7_4_7_MT2PostProcessing/src/analysisCode/Utils/allFilters_19Jan.txt
	outputFilteredFile=${workingFolder}/${name}_filtered$fileExt;
    fi;
    
    crabExt=""
    if [ ${isCrab} = 1 ]; then
	crabExt=$(ls $inputFolder/$name/)
    fi;

    #default is to not do the splitting into 10 files
    #also this should NEVER be done for MC as some scale factors need normalization
    #unless you did the preprocessing

    fileList=inputChunkList.txt

    if [ -e  inputChunkList.txt ]; then
	echo "deleting the old file list"
	rm $fileList
    fi;

    if [ ${isCrab} = 1 ]; then
	for ((i=0; i<10; i++)); do
	    if [ -d $inputFolder/${name}/$crabExt/000${i}/ ]; then
		for f in $inputFolder/$name/$crabExt/000${i}/mt2*.root; do
		    echo $f>>$fileList
		done;
	    fi;
	done;
    else
	for f in $inputFolder/$name/mt2*.root; do
	    echo $f>>$fileList
	done;
    fi;
    

    numFiles=$(wc -l inputChunkList.txt | awk '{print $1}')
    echo "number of files = " $numFiles

    if [ -d $inputFolder/${name}_ext/ ]; then
	if [ ${isCrab} = 1 ]; then
    	    crabExt=$(ls $inputFolder/${name}_ext/)
	    for ((i=0; i<10; i++)); do
		if [ -d $inputFolder/${name}_ext/$crabExt/000${i}/ ]; then
		    for f in $inputFolder/${name}_ext/${crabExt}/000${i}/mt2*.root; do
			echo $f>>$fileList
		    done;
		fi;
	    done;
	else
	    for f in $inputFolder/${name}_ext_${i}/mt2*.root; do
		echo $f>>$fileList
	    done;
	fi;
    fi;

    numFiles=$(wc -l inputChunkList.txt | awk '{print $1}')
    echo "number of files = " $numFiles


    maxNfiles=1000
    counter=-1
    preProcFile=""

    if [[ (( $doPreProc -eq 1 && (( $numFiles -gt $maxNfiles )) ))  || (( $id -lt 10 )) ]]; then
	echo "File will be split into multiple files for speed and memory limit purposes"
	counter=0;
	maxNfiles=100;
	if [[ $id > 10 ]]; then
	    preProcFile=${name}_pre.cfg;
	fi;
    fi;
    

    while (( (( (( $numFiles + 1 )) > (($counter * $maxNfiles)) )) || $(($counter < 0 )) )); do 

        counterFile=chunkPart_${name}_${counter}.txt

        ##the input file list for the current range
	if [[ $counter -gt -1 ]]; then
	    echo "running on chunks $((counter*maxNfiles+1)) to $(((counter+1)*maxNfiles))"
 	    sed -n $((counter*maxNfiles+1)),$(((counter+1)*maxNfiles))p  inputChunkList.txt > chunkPart_${name}_${counter}.txt
	    awk '$0="dcap://t3se01.psi.ch:22125/"$0' chunkPart_${name}_$counter.txt > chunkPart_${name}_${counter}_dcap.txt
	    mv chunkPart_${name}_${counter}_dcap.txt chunkPart_${name}_${counter}.txt
	fi;

	echo "Submitting the job for part" $counter; 

	#if [[ ${doPreProc} == 0 ]]; then
	if [[ $counter == -1 ]]; then
    	    awk '$0="dcap://t3se01.psi.ch:22125/"$0' inputChunkList.txt > temp_${name}_${counter}_dcap.txt
	    mv temp_${name}_${counter}_dcap.txt chunkPart_${name}_${counter}.txt
	fi;

	counterName=${name}_${counter};

	if [[ $counter == -1 ]]; then
	    outputFile=${workingFolder}/${name}$fileExt;
	    counterName=${name}
	    scriptName=batchScript_${name}.sh
	else
	    scriptName=batchScript_${counter}_${name}.sh;
	    outputFile=${workingFolder}/${counterName}$fileExt;
	    if [[ $doFilterTxt == 1 ]]; then
		outputFilteredFile=${workingFolder}/${counterName}_filtered$fileExt;
	    fi;
	fi;


	cat <<EOF > $scriptName
#!/bin/bash


#### The following configurations you should not need to change
# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N pP_${name}_`whoami`

### Specify the queue on which to run
#$ -q short.q

# Change to the current working directory from which the job go
# submitted . This will also result in the job report stdout/stderr being
# written to this directory, if you do not override it (below).
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
####v$ -o $jobsLogsFolder/${name}.out
####c$ -e $jobsLogsFolder/${name}.err
#$ -o $jobsLogsFolder/${counterName}.out
#$ -e $jobsLogsFolder/${counterName}.err

#source /mnt/t3nfs01/data01/swshare/psit3/etc/profile.d/cms_ui_env.sh
#source $VO_CMS_SW_DIR/cmsset_default.sh
#source //mnt/t3nfs01/data01/swshare/ROOT/thisroot.sh 
#eval \`scramv1 runtime -sh\`



source $VO_CMS_SW_DIR/cmsset_default.sh
#source /mnt/t3nfs01/data01/swshare/glite/external/etc/profile.d/grid-env.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export LD_LIBRARY_PATH=/mnt/t3nfs01/data01/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH
echo "Loading your CMSSW release or CMSSW_7_4_12/"
echo "from $myCMSSW"
cd $myCMSSW
eval `scramv1 runtime -sh`
cd -


mkdir -p $workingFolder
if [[ $useXRD == "true" ]]; then
    #xrdfs t3dcachedb.psi.ch mkdir -p $outputFolder # recursive mkdir does not work via xrootd
    gfal-mkdir -p ${gfalProtocol}://t3se01.psi.ch/$outputFolder # let's use gfal then
else
    gfal-mkdir -p ${gfalProtocol}://t3se01.psi.ch/$outputFolder
fi


echo "postProcessing(\"$name\",\"$counterFile\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\",\"$inputPU\",\"$PUvar\",$applyJSON,$doAllSF,$doSilver,\"$preProcFile\"); gSystem->Exit(0);"

echo "gSystem->Load(\"goodrunClass_cc.so\"); gSystem->Load(\"libCondFormatsBTauObjects.so\"); gROOT->LoadMacro(\"postProcessing.C\"); postProcessing(\"$name\",\"$counterFile\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\",\"$inputPU\",\"$PUvar\",$applyJSON,$doAllSF,$doSilver,\"$preProcFile\"); gSystem->Exit(0);" |root.exe -b -l ;



if [[ $id -lt 10 && $doFilterTxt == 1 ]]; then 
   echo "filterFromTxt(\"$filterTxt\",\"$outputFile\",1,\"$outputFilteredFile\",\"mt2\",\"\")"
   echo "gROOT->LoadMacro(\"filterFromTxt.C\"); filterFromTxt(\"$filterTxt\",\"$outputFile\",1,\"$outputFilteredFile\",\"mt2\",\"\"); gSystem->Exit(0);" |root.exe -b -l ;
   mv $outputFilteredFile $outputFile;
fi;

#######mv $outputFile $outputFolder
if [[ $useXRD == "true" ]]; then
    xrdcp -d 1 $outputFile root://t3dcachedb.psi.ch:1094/$outputFolder
else
    gfal-copy file://$outputFile ${gfalProtocol}://t3se01.psi.ch/$outputFolder
fi
rm $outputFile

#Normal skim
skimmingPruningCfg="${workingFolder}/skimmingPruning_${counterName}.cfg"
     cat skimmingPruning.cfg |grep -v \# | sed  "s#INPUTDIR#${outputFolder}#" |sed "s#INPUTFILTER#${counterName}#" \
 	| sed "s#OUTPUTDIR#${outputFolder}/skimAndPrune#" | sed "s#DOPRUNING#${doPruning}#" > \$skimmingPruningCfg
./runSkimmingPruning.sh \$skimmingPruningCfg
rm \$skimmingPruningCfg

#qcd skim
skimmingPruningCfgQCD="${workingFolder}/skimmingPruningQCD_${counterName}.cfg"
    cat skimmingPruningQCD.cfg |grep -v \# | sed  "s#INPUTDIR#${outputFolder}#" |sed "s#INPUTFILTER#${counterName}#" \
	| sed "s#OUTPUTDIR#${outputFolder}/QCDskimAndPrune#" | sed "s#DOPRUNING#${doPruning}#" > \$skimmingPruningCfgQCD
./runSkimmingPruning.sh \$skimmingPruningCfgQCD
rm \$skimmingPruningCfgQCD

#qcd skim for monojet
skimmingPruningCfgMonoJet="${workingFolder}/skimmingPruningMonoJet_${counterName}.cfg"
    cat skimmingPruningMonoJet.cfg |grep -v \# | sed  "s#INPUTDIR#${outputFolder}#" |sed "s#INPUTFILTER#${counterName}#" \
	| sed "s#OUTPUTDIR#${outputFolder}/QCDMonoJetSkimAndPrune#" | sed "s#DOPRUNING#${doPruning}#" > \$skimmingPruningCfgMonoJet
./runSkimmingPruning.sh \$skimmingPruningCfgMonoJet
rm \$skimmingPruningCfgMonoJet


#echo "is anything left in working folder? workingFolder: " 
#ls $workingFolder


EOF

        #if you have a big file and no time to change the code to be smoother: qsub  -q long.q -l h_vmem=5g batchScript_${name}.sh;

	qsub -q long.q $scriptName;
	rm $scriptName;

	if (($counter < 0)); then
	    break;
	else
	counter=$[$counter+1]; #	    counter=$[$counter+1]; 
	fi;
	
	#rm $counterFile

    done;

done < $listOfSamplesFile




###rm -f postProcessing_C.d postProcessing_C.so;

fi;

if [[ "$1" = "postCheck" ]]; then
    sizeLogsErr=`ls -l ${jobsLogsFolder}/*.err |awk '{sum+=$5} END {print sum}'`
    if [[ $sizeLogsErr -eq 0 ]]; then
	echo "there were no errors. Zipping all logs and copying them to the SE"	 
	cd $jobsLogsFolder
	tar -czvf logs.tgz  *
	if [[ $useXRD == "true" ]]; then
	    xrdcp -d 1 `pwd`/logs.tgz root://t3dcachedb.psi.ch:1094/$outputFolder
	else
	    gfal-copy file://`pwd`/logs.tgz ${gfalProtocol}://t3se01.psi.ch/$outputFolder
	fi
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
    rm -f inputChunkList.txt;
    rm -f postProcessing_C*;
    rm -f chunkPart_*.txt;
    rm -f inputChunkList.txt;
    rm -f goodruns_golden.txt;
    rm -f goodrun_cc*;
    rm -f goodrunClass_cc*;
    rm -f BTagCalibrationStandalone_cc*;
    rm -f checkOutputDir
fi


