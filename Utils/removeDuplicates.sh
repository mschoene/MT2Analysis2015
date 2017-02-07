for i in /scratch/mmasciov/SnT_Moriond2017_Feb01_data/mergedData*.root;
do
    echo $i
    n=${i%.root*}
    n=${n#*data/}
    n="/scratch/mmasciov/SnT_Moriond2017_Feb01_dataNoDuplicates/"$n".root"
    echo $n
    echo "gROOT->LoadMacro(\"removeDuplicates.C\"); removeDuplicates(\"$i\",1,\"$n\"); gSystem->Exit(0);" | root.exe -b -l;
#    echo "gROOT->LoadMacro(\"removeDuplicates.C\"); removeDuplicates(\"$i\",1,\"$n\"); gSystem->Exit(0);"
done
