#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "THStack.h"
#include "TGraphErrors.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2DrawTools.h"




int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./drawGammaControlRegionDataMC [configFileName] [lumi/shape]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  MT2DrawTools::setStyle();
  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  std::string mcFile = cfg.getEventYieldDir() + "/gammaControlRegion/mc.root";
  std::string dataFile = cfg.getEventYieldDir() + "/gammaControlRegion/data.root";

  MT2Analysis<MT2EstimateTree>* mc_   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "gammaCRtree_loose");
  // MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(dataFile, "gammaCRtree_loose");

  Double_t bins_HT[6] = {200,450, 575, 1000, 1500, 2000};

  MT2Analysis<MT2EstimateZinvGamma>* mc_ht = MT2EstimateZinvGamma::makeInclusiveEstimateFromInclusiveTree( "purity_ht", mc_ ,"", "ht", 5, bins_HT ); 



  std::string outFile = cfg.getEventYieldDir() + "/gammaControlRegion/purity_HT.root";
  mc_ht->writeToFile(outFile, "recreate");

  return 0;

}

