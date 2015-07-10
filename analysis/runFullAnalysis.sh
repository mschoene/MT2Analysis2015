# usage: ./runFullAnalysis [configFileName]
make all
./regionEventYields $1
./gammaControlRegion $1
./makeZinvGammaTemplates $1
./fitPurityGamma $1
# ./makeZinvGammaTemplates $1 MC
# ./fitPurityGamma $1 MC
# ./drawPurityFits $1 
