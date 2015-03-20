#!/bin/bash



if [[ "$#" -lt 2 ]]; then
    echo "Relunch the script with one of the following options: "
    echo "./checkLxbatchLogs.sh inputFolder checkLogs"
    echo "./checkLxbatchLogs.sh inputFolder checkRemoteFiles"
    echo "./checkLxbatchLogs.sh inputFolder checkTiming"
    echo "./checkLxbatchLogs.sh inputFolder checkRate"
    echo "./checkLxbatchLogs.sh ... "
    echo "./checkLxbatchLogs.sh inputFolder debug"
    exit;
fi;


inputFolder=$1

# --- This parses the log file and printout a message if there were issues in copying the file to the SE
# --- This script should be extended to verify if there were error messages in general
if [ $2 = "checkLogs" ]; then 
    #for x in `ls -1 $inputFolder| grep _Chunk`; do
    for x in `ls -1 $inputFolder| grep _Chunk|awk -F_Chunk '{print $1}'|uniq`; do
	echo "processing logs for sample " $x
	for z in `ls -1 $inputFolder|grep $x|grep _Chunk`; do
	    cd $inputFolder/$z
	    if [ -e std_output.txt ]; then
		if  [ -z "`tail -n 50 std_output.txt |grep -v echo|grep "remote copy succeeded"`" ]; 
		then 
		    if [ -e mt2.root ]; then
			echo "ERROR: problem in copying output to SE for job: " $x ", but file was saved locally" 
			echo "copy file by hand and everything should be fine"
			echo " "
		    else
			echo "ERROR: problem in copying output to SE for job: " $x 
			echo "File was not saved locally and job probably needs to be resubmitted"
			echo ""
		    fi		
		fi;
	    else
		echo "ERROR: no log output for job " $x;
	    fi
	    cd $OLDPWD
	done #done loop on z
    done #done loop on x
fi


# --- This is to compare directly the list of files in the SE with the number of jobs that were submitted
if [ $2 = "checkRemoteFiles" ]; then  
    if [[ "$#" -lt 3 ]]; then
	echo "You need to specify also the remote folder of the SE where the chunk files are stored"
	echo "./checkLxbatchLogs.sh /afs/../ProductionLogs checkRemoteFiles /pnfs/../productionFolder"
	exit;
    fi;

    remoteFolderOnSE=$3

    listJobs="$PWD/listJobs_$RANDOM.txt"
    listFiles="$PWD/listFiles_$RANDOM.txt"

    #echo "listJobs :" $listJobs
    #echo "listFiles :" $listFiles

    ls -1 $inputFolder/|grep _Chunk|sort -n > $listJobs
    
    cd $remoteFolderOnSE
    for x in `ls -1`; do for y in `ls -1 $x`; do echo $x/$y| sed 's#/mt2_#_Chunk#'|sed 's#.root##'; \
	done;  done |sort -n > $listFiles;
    cd $OLDPWD

    echo "diff listJobs listFilesOnSE:"
    diff $listJobs  $listFiles | grep "<";

    rm $listJobs $listFiles;
 
    
fi

# --- This produces distributions of the time/job in hours
if [ $2 = "checkTiming" ]; then  
    cat <<EOF > plotTiming.C
#include "TTree.h"
#include "TH1F.h"
#include "TPad.h"

using namespace std;

void plotTiming(string folder,string name){
 TTree tree;
 tree.ReadFile((folder+"/"+name+"_timing.txt").c_str(),"n");
 TH1F histo("histo","histo",40,0,8);
 tree.Draw("n/60./60.>>histo");
 histo.SetXTitle("time/job [hours]");
 histo.SetYTitle("#jobs ");
 histo.SetTitle(name.c_str());
 histo.SetLineColor(2);
 histo.SetLineWidth(2);
 gPad->Print((folder+"/"+name+".png").c_str());
}
EOF
    echo "gROOT->LoadMacro(\"plotTiming.C+\"); gSystem->Exit(0);" 
    echo "gROOT->LoadMacro(\"plotTiming.C+\"); gSystem->Exit(0);" |root.exe -b -l ;

    outputFolder=$inputFolder
    for sample in `ls -1 $inputFolder | grep -v txt| grep -v png| awk -F_Chunk '{print $1}'|uniq`; do
	echo $sample;
	for y in `ls -d $inputFolder/${sample}_Chunk*/std_output.txt`; do 
	    cat $y|grep "This process used approx"|tail -n 1|awk -F\( '{print $2 }'|awk '{print $1}';
	done > $outputFolder/${sample}_timing.txt; 		
	echo "gROOT->LoadMacro(\"plotTiming.C+\"); plotTiming(\"$outputFolder\",\"$sample\"); gSystem->Exit(0);" |root.exe -b -l ;
    done

rm plotTiming.C
rm plotTiming_C.d
rm plotTiming_C.so

fi


# --- This produces tables about the average rate to process events
if [ $2 = "checkRate" ]; then  
    outputFolder=$inputFolder
    for sample in `ls -1 $inputFolder | grep -v txt| grep -v png| awk -F_Chunk '{print $1}'|uniq`; do
    #for sample in TTJets; do
	#echo $sample;
	for y in `ls -d $inputFolder/${sample}_Chunk*/std_output.txt`; do 
	    cat $y|grep "event"|grep "ev/s" | awk -F\( '{print $2}'|awk '{print $1}';
	done > $outputFolder/${sample}_rates.txt;
 	cat $outputFolder/${sample}_rates.txt| awk '{sum+=$1}{n+=1}END{print sum/n}' > $outputFolder/${sample}_meanRate.txt;
	rm $outputFolder/${sample}_rates.txt;
	echo "Average rate for "$sample" is: " `cat $outputFolder/${sample}_meanRate.txt`;
    done


fi



# Debugging
if [ $2 = "debug" ]; then 
    for x in `ls -1 $inputFolder`; do
	cd $inputFolder/$x
	if [ -e std_output.err ]; then
	    echo "----- printing log for " $x " -----------"
	    tail -n 20 std_output.err
	    echo ""
	fi
	cd $OLDPWD
    done
fi


