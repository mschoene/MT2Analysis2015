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
#include "../interface/MT2EstimateSyst.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2DrawTools.h"




int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./binIsoAlongAxes [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  MT2DrawTools::setStyle();
  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  std::string mcFile = cfg.getEventYieldDir() + "/gammaControlRegion/mc.root";
  std::string dataFile = cfg.getEventYieldDir() + "/gammaControlRegion/data.root";

  MT2Analysis<MT2EstimateTree>* mc_   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "gammaCRtree_loose");
  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(dataFile, "gammaCRtree_loose");

  Double_t bins_ht[] = {200,450, 575, 1000, 1500,2000};
  int size_ht = sizeof(bins_ht)/sizeof(double)-1;
  
  Double_t bins_njets[] = {2,4,7,12};
  int size_njets = sizeof(bins_njets)/sizeof(double)-1;
  
  Double_t bins_nbjets[] = {0,1,2,3,6};
  int size_nbjets = sizeof(bins_nbjets)/sizeof(double)-1;

  
  MT2Analysis<MT2EstimateZinvGamma>* mc_ht = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_ht", data ,"", "ht" , size_ht, bins_ht ); 

  MT2Analysis<MT2EstimateZinvGamma>* mc_njets = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_njets", data ,"", "nJets", size_njets, bins_njets ); 

  MT2Analysis<MT2EstimateZinvGamma>* mc_nbjets = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_nbjets", data ,"", "nBJets", size_nbjets, bins_nbjets );
 

  std::string outFile = cfg.getEventYieldDir() + "/gammaControlRegion/iso_ht.root";
  mc_ht->writeToFile(outFile, "recreate");

  std::string outFile_njets = cfg.getEventYieldDir() + "/gammaControlRegion/iso_nJets.root";
  mc_njets->writeToFile(outFile_njets, "recreate");

  std::string outFile_nbjets = cfg.getEventYieldDir() + "/gammaControlRegion/iso_nBJets.root";
  mc_nbjets->writeToFile(outFile_nbjets, "recreate");


  //need to also create the purity from MC binned in ht, njets & nbjets
  //just so that it is binned correctly, could also just use mc_ otherwise
  MT2Analysis<MT2EstimateTree>* all_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_ht, bins_ht, "weight" , "ht" ) ;
  MT2Analysis<MT2EstimateTree>* prompt_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_ht, bins_ht, "weight*(prompt==2 || prompt==1)" , "ht" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_ht = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)prompt_ht, (MT2Analysis<MT2Estimate>*)all_ht);
  std::string outFile_purity_ht = cfg.getEventYieldDir() + "/gammaControlRegion/purityMC_ht.root";

  MT2Analysis<MT2EstimateTree>* all_pass_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_ht, bins_ht, "(iso<2.5)","ht");
  MT2Analysis<MT2EstimateTree>* prompt_tight_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_ht, bins_ht, "((prompt==2 || prompt==1) && iso<2.5)" , "ht" ) ;
 MT2Analysis<MT2EstimateSyst>* purityTight_ht = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)prompt_tight_ht, (MT2Analysis<MT2Estimate>*)all_pass_ht);

  purityLoose_ht->writeToFile( outFile_purity_ht );
  purityTight_ht->addToFile( outFile_purity_ht );
 

  MT2Analysis<MT2EstimateTree>* all_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree("all", cfg.regionsSet(), mc_ , size_njets, bins_njets, "" , "nJets");
  MT2Analysis<MT2EstimateTree>* prompt_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_njets, bins_njets, "prompt==2 || prompt==1" , "nJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_njets = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)prompt_njets, (MT2Analysis<MT2Estimate>*)all_njets);
  std::string outFile_purity_njets = cfg.getEventYieldDir() + "/gammaControlRegion/purityMC_njets.root";

  MT2Analysis<MT2EstimateTree>* all_pass_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_njets, bins_njets, "(iso<2.5)","nJets");
  MT2Analysis<MT2EstimateTree>* prompt_tight_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_njets, bins_njets, "((prompt==2 || prompt==1) && iso<2.5)" , "nJets" ) ;
 MT2Analysis<MT2EstimateSyst>* purityTight_njets = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)prompt_tight_njets, (MT2Analysis<MT2Estimate>*)all_pass_njets);

  purityLoose_njets->writeToFile( outFile_purity_njets );
  purityTight_njets->addToFile( outFile_purity_njets );



  MT2Analysis<MT2EstimateTree>* all_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree("all", cfg.regionsSet(), mc_ , size_nbjets, bins_nbjets, "" , "nBJets");
  MT2Analysis<MT2EstimateTree>* prompt_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_nbjets, bins_nbjets, "prompt==2 || prompt==1" , "nBJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_nbjets = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)prompt_nbjets, (MT2Analysis<MT2Estimate>*)all_nbjets);
  std::string outFile_purity_nbjets = cfg.getEventYieldDir() + "/gammaControlRegion/purityMC_nbjets.root";

  MT2Analysis<MT2EstimateTree>* all_pass_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_nbjets, bins_nbjets, "(iso<2.5)","nBJets");
  MT2Analysis<MT2EstimateTree>* prompt_tight_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all", cfg.regionsSet(), mc_ , size_nbjets, bins_nbjets, "((prompt==2 || prompt==1) && iso<2.5)" , "nBJets" ) ;
 MT2Analysis<MT2EstimateSyst>* purityTight_nbjets = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", cfg.regionsSet(), (MT2Analysis<MT2Estimate>*)prompt_tight_nbjets, (MT2Analysis<MT2Estimate>*)all_pass_nbjets);

  purityLoose_nbjets->writeToFile( outFile_purity_nbjets );
  purityTight_nbjets->addToFile( outFile_purity_nbjets );


  return 0;

}

