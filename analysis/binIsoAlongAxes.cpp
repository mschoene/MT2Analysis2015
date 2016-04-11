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

  std::string mcFile = cfg.getGammaCRdir() + "/mc.root";
  std::string dataFile = cfg.getGammaCRdir() + "/data.root";

  MT2Analysis<MT2EstimateTree>* mc   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "gammaCRtree_loose");
  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(dataFile, "gammaCRtree_loose");

  std::string outputDir = cfg.getGammaCRdir() + "/isolationAndPurity/";
  system(Form("mkdir -p %s", outputDir.c_str()));

  Double_t bins_central[] = {200,3000};
  int size_central = sizeof(bins_central)/sizeof(double)-1;

  MT2Analysis<MT2EstimateZinvGamma>* mc_central = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_central", data, cfg.crRegionsSet(), "nJets>1" , "ht" , size_central, bins_central ); 
  std::string outFile_central = outputDir + "/iso_central.root";
  mc_central->writeToFile(outFile_central, "recreate");



  Double_t bins_ht[] = {200,450, 575, 1000, 1500,3000};
  int size_ht = sizeof(bins_ht)/sizeof(double)-1;
  Double_t bins_njets[] = {2,4,7,12};
  int size_njets = sizeof(bins_njets)/sizeof(double)-1;
  Double_t bins_nbjets[] = { 0, 1, 2, 3, 6};
  int size_nbjets = sizeof(bins_nbjets)/sizeof(double)-1;
  Double_t bins_mono_nbjets[] = { 0, 1, 2};
  int size_mono_nbjets = sizeof(bins_mono_nbjets)/sizeof(double)-1;

  Double_t bins_mono_ht[] = {200,250,350,450, 575, 700, 1000, 1500};
  int size_mono_ht = sizeof(bins_mono_ht)/sizeof(double)-1;

  Double_t bins_mt2[] = { 200, 1500};
  int size_mt2 = sizeof(bins_mt2)/sizeof(double)-1;
  
  Double_t bins_incl_njets[] = {1,2,4,7,12};
  int size_incl_njets = sizeof(bins_incl_njets)/sizeof(double)-1;



  MT2Analysis<MT2EstimateZinvGamma>* mc_mt2 = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_mt2", data, "zurichPlus" ,"nJets>0", "mt2" , size_mt2, bins_mt2 ); 
  std::string outFile_mt2 = outputDir + "/iso_mt2.root";
  mc_mt2->writeToFile(outFile_mt2, "recreate");


  MT2Analysis<MT2EstimateZinvGamma>* mc_mono_ht = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_mono_ht", data, cfg.crRegionsSet() , "nJets>1" , "ht" , size_mono_ht, bins_mono_ht ); 
  std::string outFile_mono_ht = outputDir + "/iso_mono_ht.root";
  mc_mono_ht->writeToFile(outFile_mono_ht, "recreate");
  

  // MT2Analysis<MT2EstimateZinvGamma>* mc_ht = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_ht", data, cfg.crRegionsSet() , "nJets>1" , "ht" , size_ht, bins_ht ); 
  // std::string outFile = outputDir + "/iso_ht.root";
  // mc_ht->writeToFile(outFile, "recreate");

  // MT2Analysis<MT2EstimateZinvGamma>* mc_njets = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_njets", data, cfg.crRegionsSet() ,"nJets>1", "nJets", size_njets, bins_njets ); 
  // std::string outFile_njets = outputDir + "/iso_nJets.root";
  // mc_njets->writeToFile(outFile_njets, "recreate");

  // MT2Analysis<MT2EstimateZinvGamma>* mc_nbjets = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_nbjets", data, cfg.crRegionsSet() ,"nJets>1", "nBJets", size_nbjets, bins_nbjets );
  // std::string outFile_nbjets = outputDir + "/iso_nBJets.root";
  // mc_nbjets->writeToFile(outFile_nbjets, "recreate");
 
  // MT2Analysis<MT2EstimateZinvGamma>* mc_mono_nbjets = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_mono_nbjets", data, cfg.crRegionsSet() ,"nJets==1", "nBJets", size_mono_nbjets, bins_mono_nbjets );
  // std::string outFile_mono_nbjets = outputDir + "/iso_mono_nBJets.root";
  // mc_mono_nbjets->writeToFile(outFile_mono_nbjets, "recreate");
  
   

  MT2Analysis<MT2EstimateZinvGamma>* incl_ht = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_incl_ht", data, cfg.crRegionsSet() , "" , "ht" , size_ht, bins_ht ); 
  std::string outFile_incl_ht = outputDir + "/iso_incl_ht.root";
  incl_ht->writeToFile(outFile_incl_ht, "recreate");

  MT2Analysis<MT2EstimateZinvGamma>* incl_njets = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_incl_njets", data, cfg.crRegionsSet() ,"", "nJets", size_incl_njets, bins_incl_njets ); 
  std::string outFile_incl_njets = outputDir + "/iso_incl_nJets.root";
  incl_njets->writeToFile(outFile_incl_njets, "recreate");

  MT2Analysis<MT2EstimateZinvGamma>* incl_nbjets = MT2EstimateZinvGamma::makeInclusiveAnalysisFromInclusiveTree( "iso_incl_nbjets", data, cfg.crRegionsSet() ,"" , "nBJets", size_nbjets, bins_nbjets );
  std::string outFile_incl_nbjets = outputDir + "/iso_incl_nBJets.root";
  incl_nbjets->writeToFile(outFile_incl_nbjets, "recreate");
 



  std::string multijet = "(nJets>1)";
  std::string monojet = "(nJets==1)";
  std::string iso = "(iso<2.5)";
  std::string prompt = "(prompt==1 || prompt==2)";
  
  
  //need to also create the purity from MC binned in ht, njets & nbjets
////  //just so that it is binned correctly, could also just use mc otherwise
////  MT2Analysis<MT2EstimateTree>* all_mt2 =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_mt2", "zurichPlus", mc , "", size_mt2, bins_mt2, "mt2" ) ;
////  MT2Analysis<MT2EstimateTree>* prompt_mt2 =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_mt2", "zurichPlus", mc , prompt, size_mt2, bins_mt2, "mt2" ) ;
////
////  MT2Analysis<MT2EstimateSyst>* purityLoose_mt2 = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", "zurichPlus", (MT2Analysis<MT2Estimate>*)prompt_mt2, (MT2Analysis<MT2Estimate>*)all_mt2);
//// 
////  MT2Analysis<MT2EstimateTree>* all_pass_mt2 =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_pass_mt2", "zurichPlus", mc , iso, size_mt2, bins_mt2,"mt2");
////  MT2Analysis<MT2EstimateTree>* prompt_tight_mt2 =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_tight_mt2", "zurichPlus", mc , iso+"&&"+prompt , size_mt2, bins_mt2, "mt2" ) ;
////  MT2Analysis<MT2EstimateSyst>* purityTight_mt2 = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", "zurichPlus", (MT2Analysis<MT2Estimate>*)prompt_tight_mt2, (MT2Analysis<MT2Estimate>*)all_pass_mt2);
////
////  std::string outFile_purity_mt2 = cfg.getEventYieldDir() + "/gammaControlRegion/purityMC_mt2.root";
////  purityLoose_mt2->writeToFile( outFile_purity_mt2,"recreate" );
////  purityTight_mt2->addToFile( outFile_purity_mt2 );
  
  
  //NEW HT BINNING FOR THE BINNED MONO JET
  MT2Analysis<MT2EstimateTree>* all_mono_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_mono_ht", cfg.crRegionsSet(), mc , "(id>151.)", size_mono_ht, bins_mono_ht, "ht" ) ;
  
  MT2Analysis<MT2EstimateTree>* prompt_mono_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_mono_ht", cfg.crRegionsSet(), mc ,  prompt+"&& (id>151.)", size_mono_ht, bins_mono_ht, "ht" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_mono_ht = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose",   (MT2Analysis<MT2Estimate>*)prompt_mono_ht, (MT2Analysis<MT2Estimate>*)all_mono_ht);
 
  MT2Analysis<MT2EstimateTree>* all_pass_mono_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_pass_mono_ht", cfg.crRegionsSet(), mc , iso+"&& (id>151.)", size_mono_ht, bins_mono_ht,"ht");
  MT2Analysis<MT2EstimateTree>* prompt_tight_mono_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_tight_mono_ht", cfg.crRegionsSet(), mc , iso+"&&"+prompt+"&& (id>151.)" , size_mono_ht, bins_mono_ht, "ht" ) ;

  //  std::string outFile_yield_mono_ht = cfg.getEventYieldDir() + "/gammaControlRegion/yieldMC_mono_ht.root";
  //  all_pass_mono_ht->writeToFile( outFile_yield_mono_ht,"RECREATE" );
  //  prompt_tight_mono_ht->addToFile( outFile_yield_mono_ht );
 

  MT2Analysis<MT2EstimateSyst>* purityTight_mono_ht = MT2EstimateSyst::makeEfficiencyAnalysis( "purity",   (MT2Analysis<MT2Estimate>*)prompt_tight_mono_ht, (MT2Analysis<MT2Estimate>*)all_pass_mono_ht);

  std::string outFile_purity_mono_ht = outputDir + "/purityMC_mono_ht.root";
  purityLoose_mono_ht->writeToFile( outFile_purity_mono_ht,"recreate" );
  purityTight_mono_ht->addToFile( outFile_purity_mono_ht );
  

  /*
  //MULTIJET HT
  MT2Analysis<MT2EstimateTree>* all_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_ht", cfg.crRegionsSet(), mc , multijet+"&& (id>151.)", size_ht, bins_ht, "ht" ) ;

  MT2Analysis<MT2EstimateTree>* prompt_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_ht", cfg.crRegionsSet(), mc , multijet +"&&"+ prompt+"&& (id>151.)", size_ht, bins_ht, "ht" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_ht = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", cfg.crRegionsSet(), (MT2Analysis<MT2Estimate>*)prompt_ht, (MT2Analysis<MT2Estimate>*)all_ht);
 
  MT2Analysis<MT2EstimateTree>* all_pass_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_pass_ht", cfg.crRegionsSet(), mc , multijet +"&&"+iso+"&& (id>151.)", size_ht, bins_ht,"ht");
  MT2Analysis<MT2EstimateTree>* prompt_tight_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_tight_ht", cfg.crRegionsSet(), mc , multijet+"&&"+iso+"&&"+prompt+"&& (id>151.)" , size_ht, bins_ht, "ht" ) ;

  std::string outFile_yield_ht = cfg.getEventYieldDir() + "/gammaControlRegion/yieldMC_ht.root";
  all_pass_ht->writeToFile( outFile_yield_ht,"RECREATE" );
  prompt_tight_ht->addToFile( outFile_yield_ht );
 

  MT2Analysis<MT2EstimateSyst>* purityTight_ht = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", cfg.crRegionsSet(), (MT2Analysis<MT2Estimate>*)prompt_tight_ht, (MT2Analysis<MT2Estimate>*)all_pass_ht);

  std::string outFile_purity_ht = outputDir + "/purityMC_ht.root";
  purityLoose_ht->writeToFile( outFile_purity_ht,"recreate" );
  purityTight_ht->addToFile( outFile_purity_ht );
  


  //NJETS
  MT2Analysis<MT2EstimateTree>* all_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree("all_njets", cfg.crRegionsSet(), mc , multijet + " && (id>151.)" , size_njets, bins_njets, "nJets");
  MT2Analysis<MT2EstimateTree>* prompt_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_njets", cfg.crRegionsSet(), mc , multijet +"&&"+ prompt+"&&(id>151.)", size_njets, bins_njets, "nJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_njets = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", cfg.crRegionsSet(), (MT2Analysis<MT2Estimate>*)prompt_njets, (MT2Analysis<MT2Estimate>*)all_njets);
 
  MT2Analysis<MT2EstimateTree>* all_pass_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_pass_njets", cfg.crRegionsSet(), mc ,  multijet + " && " + iso + " && (id>151.)", size_njets, bins_njets,"nJets");
  MT2Analysis<MT2EstimateTree>* prompt_tight_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_tight_njets", cfg.crRegionsSet(), mc ,  multijet + " && " + iso + " && "+ prompt + " && (id>151.)", size_njets, bins_njets, "nJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityTight_njets = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", cfg.crRegionsSet(), (MT2Analysis<MT2Estimate>*)prompt_tight_njets, (MT2Analysis<MT2Estimate>*)all_pass_njets);


  std::string outFile_purity_njets = outputDir + "/purityMC_njets.root";
  purityLoose_njets->writeToFile( outFile_purity_njets ,"recreate");
  purityTight_njets->addToFile( outFile_purity_njets );



  MT2Analysis<MT2EstimateTree>* all_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree("all_nbjets", cfg.crRegionsSet(), mc,  multijet+"&& (id>151.)" , size_nbjets, bins_nbjets, "nBJets");
  MT2Analysis<MT2EstimateTree>* prompt_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_nbjets", cfg.crRegionsSet(), mc , multijet+"&&"+prompt+"&& (id>151.)" , size_nbjets, bins_nbjets, "nBJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_nbjets = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", cfg.crRegionsSet(), (MT2Analysis<MT2Estimate>*)prompt_nbjets, (MT2Analysis<MT2Estimate>*)all_nbjets);
 
  MT2Analysis<MT2EstimateTree>* all_pass_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_pass_nbjets", cfg.crRegionsSet(), mc ,  multijet +"&&"+ iso+"&& (id>151.)" , size_nbjets, bins_nbjets, "nBJets");
  MT2Analysis<MT2EstimateTree>* prompt_tight_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_tight_nbjets", cfg.crRegionsSet(), mc ,  multijet+"&&"+iso+"&&"+ prompt+"&& (id>151.)" , size_nbjets, bins_nbjets, "nBJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityTight_nbjets = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", cfg.crRegionsSet(), (MT2Analysis<MT2Estimate>*)prompt_tight_nbjets, (MT2Analysis<MT2Estimate>*)all_pass_nbjets);


  std::string outFile_purity_nbjets = outputDir + "/purityMC_nbjets.root";
  purityLoose_nbjets->writeToFile( outFile_purity_nbjets ,"recreate");
  purityTight_nbjets->addToFile( outFile_purity_nbjets );
  */








  //////INCLUSIVE REGION//////////////

  MT2Analysis<MT2EstimateTree>* all_incl_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_incl_ht", cfg.crRegionsSet(), mc , "(id>151.)", size_ht, bins_ht, "ht" ) ;

  MT2Analysis<MT2EstimateTree>* prompt_incl_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_incl_ht", cfg.crRegionsSet(), mc ,  prompt+"&& (id>151.)", size_ht, bins_ht, "ht" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_incl_ht = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose",   (MT2Analysis<MT2Estimate>*)prompt_incl_ht, (MT2Analysis<MT2Estimate>*)all_incl_ht);
 
  MT2Analysis<MT2EstimateTree>* all_pass_incl_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_pass_incl_ht", cfg.crRegionsSet(), mc , iso+"&& (id>151.)", size_ht, bins_ht,"ht");
  MT2Analysis<MT2EstimateTree>* prompt_tight_incl_ht =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_tight_incl_ht", cfg.crRegionsSet(), mc , iso+"&&"+prompt+"&& (id>151.)" , size_ht, bins_ht, "ht" ) ;
  /*
  std::string outFile_yield_incl_ht = cfg.getEventYieldDir() + "/gammaControlRegion/yieldMC_incl_ht.root";
  all_pass_incl_ht->writeToFile( outFile_yield_incl_ht,"RECREATE" );
  prompt_tight_incl_ht->addToFile( outFile_yield_incl_ht );
  */
  MT2Analysis<MT2EstimateSyst>* purityTight_incl_ht = MT2EstimateSyst::makeEfficiencyAnalysis( "purity",   (MT2Analysis<MT2Estimate>*)prompt_tight_incl_ht, (MT2Analysis<MT2Estimate>*)all_pass_incl_ht);

  std::string outFile_purity_incl_ht = outputDir + "/purityMC_incl_ht.root";
  purityLoose_incl_ht->writeToFile( outFile_purity_incl_ht,"recreate" );
  purityTight_incl_ht->addToFile( outFile_purity_incl_ht );
  


  //NJETS
  MT2Analysis<MT2EstimateTree>* all_incl_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree("all_incl_njets", cfg.crRegionsSet(), mc , "(id>151.)" , size_incl_njets, bins_incl_njets, "nJets");
  MT2Analysis<MT2EstimateTree>* prompt_incl_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_incl_njets", cfg.crRegionsSet(), mc , prompt+"&&(id>151.)", size_incl_njets, bins_incl_njets, "nJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_incl_njets = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose",   (MT2Analysis<MT2Estimate>*)prompt_incl_njets, (MT2Analysis<MT2Estimate>*)all_incl_njets);
 
  MT2Analysis<MT2EstimateTree>* all_pass_incl_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_pass_incl_njets", cfg.crRegionsSet(), mc , iso + " && (id>151.)", size_incl_njets, bins_incl_njets,"nJets");
  MT2Analysis<MT2EstimateTree>* prompt_tight_incl_njets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_tight_incl_njets", cfg.crRegionsSet(), mc , iso + " && "+ prompt + " && (id>151.)", size_incl_njets, bins_incl_njets, "nJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityTight_incl_njets = MT2EstimateSyst::makeEfficiencyAnalysis( "purity",   (MT2Analysis<MT2Estimate>*)prompt_tight_incl_njets, (MT2Analysis<MT2Estimate>*)all_pass_incl_njets);


  std::string outFile_purity_incl_njets = outputDir + "/purityMC_incl_njets.root";
  purityLoose_incl_njets->writeToFile( outFile_purity_incl_njets ,"recreate");
  purityTight_incl_njets->addToFile( outFile_purity_incl_njets );



  MT2Analysis<MT2EstimateTree>* all_incl_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree("all_incl_nbjets", cfg.crRegionsSet(), mc,  "(id>151.)" , size_nbjets, bins_nbjets, "nBJets");
  MT2Analysis<MT2EstimateTree>* prompt_incl_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_incl_nbjets", cfg.crRegionsSet(), mc , prompt+"&& (id>151.)" , size_nbjets, bins_nbjets, "nBJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityLoose_incl_nbjets = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose",   (MT2Analysis<MT2Estimate>*)prompt_incl_nbjets, (MT2Analysis<MT2Estimate>*)all_incl_nbjets);
 
  MT2Analysis<MT2EstimateTree>* all_pass_incl_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "all_pass_incl_nbjets", cfg.crRegionsSet(), mc , iso+"&& (id>151.)" , size_nbjets, bins_nbjets, "nBJets");
  MT2Analysis<MT2EstimateTree>* prompt_tight_incl_nbjets =  MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "prompt_tight_incl_nbjets", cfg.crRegionsSet(), mc , iso+"&&"+ prompt+"&& (id>151.)" , size_nbjets, bins_nbjets, "nBJets" ) ;
  MT2Analysis<MT2EstimateSyst>* purityTight_incl_nbjets = MT2EstimateSyst::makeEfficiencyAnalysis( "purity",   (MT2Analysis<MT2Estimate>*)prompt_tight_incl_nbjets, (MT2Analysis<MT2Estimate>*)all_pass_incl_nbjets);


  std::string outFile_purity_incl_nbjets = outputDir + "/purityMC_incl_nbjets.root";
  purityLoose_incl_nbjets->writeToFile( outFile_purity_incl_nbjets ,"recreate");
  purityTight_incl_nbjets->addToFile( outFile_purity_incl_nbjets );



  return 0;

}

