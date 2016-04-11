# usage: ./runFullAnalysis [configFileName]
make all
doPlots=0


#first you will want to run and create the babies baby tree:
./regionEventYields $1 $

./gammaControlRegion $1 $
./makeZinvGammaTemplates $1 data RC $
# ./makeZinvGammaTemplates $1 MC MC

./zllControlRegion $1

./llepControlRegion $1 $


./qcdControlRegion $1 #data/MC monojet
./qcdControlRegion $1 all monojet #check all

echo "Wait till all the baby-trees have been produced, this will take a couple of minutes"
wait 
echo "Finished all the background processes"

#the neccessary steps before creating the DCs
./binIsoAlongAxes $1
./fitPurityGamma $1 data RC axes

#QCD fits and estimate
./computeQCDFromDeltaPhi $1 MC
./computeQCDFromDeltaPhi $1 data
./drawQCDForMonoJet $1


if [ doPlots == 1 ]; then

./drawPurityFits $1 data RC axes

./drawGammaControlRegionDataMC $1
./drawGammaControlRegionDataMC $1 shape

./drawQCDFromDeltaPhi $1

fi


./computeZinvFromGamma $1

#./computeZllGammaRatio $1


./computeZllGammaRatio $1

./createDatacards $1
