###inputdir="/scratch/mmasciov/SnT_Moriond2017_Feb01_forQCD/"
inputdir="/scratch/mmasciov/SnT_Moriond2017_Feb01_wjetsFix_forQCD/"
outputdir="/scratch/mmasciov/SnT_Moriond2017_Feb01_updatedLepSF_forQCD/"
for i in $(ls /scratch/mmasciov/SnT_Moriond2017_Feb01_wjetsFix_forQCD/)
do
    name=${i%.root*}
    echo $name
    newname=$outputdir$name".root"
    echo "gROOT->LoadMacro(\"leptonSF.C\"); setLeptonSF(\"$name\",\"$inputdir\",\"$newname\"); gSystem->Exit(0);" | root.exe -b -l;
done
