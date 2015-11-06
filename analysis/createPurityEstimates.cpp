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
    std::cout << "USAGE: ./createPurityEstimates [configFileName]" << std::endl;
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

  Double_t bins_ht[] = {200,450, 575, 1000, 1500};
  int size_ht = sizeof(bins_ht)/sizeof(double)-1;
  
  Double_t bins_njets[] = {2,4,7,12};
  int size_njets = sizeof(bins_njets)/sizeof(double)-1;
  
  Double_t bins_nbjets[] = {0,1,2,3,6};
  int size_nbjets = sizeof(bins_nbjets)/sizeof(double)-1;

  
  MT2Analysis<MT2EstimateZinvGamma>* mc_ht = MT2EstimateZinvGamma::makeInclusiveEstimateFromInclusiveTree( "purity_ht", mc_ ,"", "ht", size_ht, bins_ht ); 

  MT2Analysis<MT2EstimateZinvGamma>* mc_njets = MT2EstimateZinvGamma::makeInclusiveEstimateFromInclusiveTree( "purity_njets", mc_ ,"", "nJets", size_njets, bins_njets ); 

  MT2Analysis<MT2EstimateZinvGamma>* mc_nbjets = MT2EstimateZinvGamma::makeInclusiveEstimateFromInclusiveTree( "purity_nbjets", mc_ ,"", "nBJets", size_nbjets, bins_nbjets );
 

  std::string outFile = cfg.getEventYieldDir() + "/gammaControlRegion/purity_HT.root";
  mc_ht->writeToFile(outFile, "recreate");

  std::string outFile_njets = cfg.getEventYieldDir() + "/gammaControlRegion/purity_nJets.root";
  mc_njets->writeToFile(outFile_njets, "recreate");

  std::string outFile_nbjets = cfg.getEventYieldDir() + "/gammaControlRegion/purity_nBJets.root";
  mc_nbjets->writeToFile(outFile_nbjets, "recreate");

  return 0;

}

