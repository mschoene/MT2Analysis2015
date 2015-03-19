#!/bin/bash

# --- configuration (consider to move this into a separate file) ---
treeName="tree"
#inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/casal/babies/PHYS14_Production_QCDpt_noSietaieta/"
inputFolder="/pnfs/psi.ch/cms/trivcat/store/user/mangano/babies/testMarch7attempt8/"
productionName="prodMarch9"
fileExt="_post.root"
# --------------------------



# initialization
jobsLogsFolder="./$productionName"
outputFolder="/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/babies/postprocessed/"$productionName/
workingFolder="/scratch/`whoami`/"$productionName


if [[ "$#" -eq 0 ]]; then
    echo "Relunch the script with one of the following options: "
    echo "./doTreeProduction.sh post"
    echo "./doTreeProduction.sh postCheck"
    echo "./doTreeProduction.sh clean"
fi;

if [[ "$1" = "post" ]]; then


if [ -d "$jobsLogsFolder" ]; then 
    echo "ERROR: the logFolder" $jobsLogsFolder " already exists."
    echo "Delete it and start from a clean area, or redirect the logs in a different place."
    echo "Exiting ..."
    exit
else
    mkdir  $jobsLogsFolder
fi

echo "Location of log files is: " $jobsLogsFolder
echo "Location of final files on SE is: " $outputFolder
echo "Working folder on working-node is: " $workingFolder

# here I compile the root macro only once
echo "gROOT->LoadMacro(\"postProcessing.C+\"); gSystem->Exit(0);" |root.exe -b -l ;


while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    id=`echo $line |awk '{print $1}'`
    name=`echo $line |awk '{print $2}'`
    xsec=`echo $line |awk '{print $3}'`
    filter=`echo $line |awk '{print $4}'`
    kfactor=`echo $line |awk '{print $5}'`
    
    outputFile=${workingFolder}/${name}$fileExt;


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

echo "postProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id);"
echo "gROOT->LoadMacro(\"postProcessing.C\"); postProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id); gSystem->Exit(0);" |root.exe -b -l ;

#mv $outputFile $outputFolder
gfal-copy file://$outputFile srm://t3se01.psi.ch/$outputFolder
rm $outputFile

skimmingPruningCfg="${workingFolder}/skimmingPruning_${name}.cfg"
    cat skimmingPruning.cfg |grep -v \# | sed  "s#INPUTDIR#${outputFolder}#" |sed "s#INPUTFILTER#${name}#" \
	| sed "s#OUTPUTDIR#${outputFolder}/skimAndPrune#" > \$skimmingPruningCfg

./runSkimmingPruning.sh \$skimmingPruningCfg
rm \$skimmingPruningCfg

#echo "is anything left in working folder? workingFolder: " 
#ls $workingFolder


EOF

qsub batchScript_${name}.sh;
rm batchScript_${name}.sh;

done < postProcessing.cfg
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