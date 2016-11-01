#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>



#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2DrawTools.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2EstimateSyst.h"

#define mt2_cxx
#include "../interface/mt2.h"


#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLorentzVector.h"



bool use_extrapolation = true;
//bool use_extrapolation = false;
bool do_dummyMC = true;

int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeLostLepton                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc!=2 ) {
    std::cout << "USAGE: ./computeLostLepton [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  TH1::AddDirectory(kTRUE);


  std::string llepControlRegionDir = cfg.getEventYieldDir() + "/llepControlRegion"; 


  // Data Control Region
  MT2Analysis<MT2Estimate>* llepCR;
  if( !do_dummyMC )
    llepCR  = MT2Analysis<MT2Estimate>::readFromFile(llepControlRegionDir + "/data.root", "llepCR");
  else{
    llepCR  = MT2Analysis<MT2Estimate>::readFromFile(llepControlRegionDir + "/mc.root", "llepCR");
    (*llepCR) *= cfg.lumi();
  }

  if( llepCR==0 ) {
    std::cout << "-> Please run llepControlRegion first. I need to get the llepCR yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(193);
  }

  // MC Control Region
  MT2Analysis<MT2Estimate>* MCcr_ = MT2Analysis<MT2Estimate>::readFromFile(llepControlRegionDir + "/mc.root", "llepCR");
  if( MCcr_==0 ) {
    std::cout << "-> Please run llepControlRegion first. I need to get the MC llepCR yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(193);
  }

  MT2Analysis<MT2Estimate>* MCcr = MT2Analysis<MT2Estimate>::readFromFile(llepControlRegionDir + "/mc.root", "llepCR");
  MCcr->setName("llepCR_mc");
  
  MT2Analysis<MT2Estimate>* MCcr_integral;
  if( use_extrapolation ){
    MCcr_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "MCcr_integral", cfg.regionsSet(), MCcr_ );
    MCcr_integral->setName("llepCR_mc_integral");
  }
  
  // MC Signal Region
  MT2Analysis<MT2Estimate>* Top   = MT2Analysis<MT2Estimate>::readFromFile(cfg.getEventYieldDir() + "/analyses.root", "Top");
  MT2Analysis<MT2Estimate>* WJets = MT2Analysis<MT2Estimate>::readFromFile(cfg.getEventYieldDir() + "/analyses.root", "WJets");
  if( Top==0 || WJets==0 ) {
    std::cout << "-> Please run regionEventYields on MC first. I need to get the Top and W+jets MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }

  MT2Analysis<MT2Estimate>* MCsr = new MT2Analysis<MT2Estimate>( *(Top) );
  (*MCsr) += (*(WJets));
  MCsr->setName("Top + W+jets");

  // MC ratio SR/CR
  MT2Analysis<MT2Estimate>* RatioMC = new MT2Analysis<MT2Estimate>( "llepRatioMC", cfg.regionsSet() );
  if( !use_extrapolation )
    (*RatioMC) = ( (*MCsr) / (*MCcr) );
  else
    (*RatioMC) = ( (*MCsr) / (*MCcr_integral) );

  MT2Analysis<MT2EstimateSyst>* llep_est_ = MT2EstimateSyst::makeAnalysisFromEstimate( "llep_est_", cfg.regionsSet(), llepCR );

  MT2Analysis<MT2EstimateSyst>* llep_est;
  if( !use_extrapolation )
    llep_est = MT2EstimateSyst::makeAnalysisFromEstimate( "llep_est", cfg.regionsSet(), llepCR );

  MT2Analysis<MT2EstimateSyst>* llep_est_integral;
  if( use_extrapolation )
    llep_est_integral = MT2EstimateSyst::makeIntegralAnalysisFromEstimate( "llep_est_integral", cfg.regionsSet(), llep_est_ );


  // Lost lepton estimate
  MT2Analysis<MT2Estimate>* llepEstimate = new MT2Analysis<MT2Estimate>( "llepEstimate", cfg.regionsSet() );
  if( !use_extrapolation )
    (*llepEstimate) = (* (MT2Analysis<MT2Estimate>*)llep_est) * (*RatioMC);
  else
    (*llepEstimate) = (* (MT2Analysis<MT2Estimate>*)llep_est_integral) * (*RatioMC);


  std::string outFile = cfg.getEventYieldDir() + "/llepEstimate";
  outFile += ".root";

  llepEstimate->writeToFile( outFile, "recreate" );
  RatioMC->addToFile( outFile );
  MCsr->addToFile( outFile );
  if( !use_extrapolation ){
    llep_est->addToFile( outFile );
    MCcr->addToFile( outFile );
  }
  else{
    llep_est_integral->addToFile( outFile );
    MCcr_integral->addToFile( outFile );
  }

  return 0;

}
