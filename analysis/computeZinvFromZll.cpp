#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>


#include "../interface/MT2Analysis.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2EstimateSyst.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Sample.h"
#include "../interface/MT2Config.h"

#define mt2_cxx
#include "../interface/mt2.h"


#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"



MT2Analysis<MT2Estimate>* combineDataAndMC( MT2Analysis<MT2Estimate>* data, MT2Analysis<MT2Estimate>* mc );
MT2Analysis<MT2Estimate>* getInclusiveRatioMC( const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* Zinv, MT2Analysis<MT2EstimateTree>* gammaCRtree );



int main(int argc, char* argv[]) {

  
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeZinvFromGamma             |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc!=2 ) {
    std::cout << "USAGE: ./computeZinvFromGamma [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string zllControlRegionDir = cfg.getEventYieldDir() + "/zllControlRegion";

  MT2Analysis<MT2Estimate>* zllCR = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/data.root", "zllCR");
  MT2Analysis<MT2Estimate>* zllCR_mc = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc.root", "zllCR");                                                                     

  MT2Analysis<MT2Estimate>* zllCR_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "zllCR_integral", cfg.regionsSet(), zllCR );
  MT2Analysis<MT2Estimate>* zllCR_mc_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "zllCR_mc_integral", cfg.regionsSet(), zllCR_mc );

  if( zllCR==0 || zllCR_mc==0 ) {
    std::cout << "-> Please run zllControlRegion first. I need to get the zllCR yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(193);
  }

  MT2Analysis<MT2EstimateTree>* Zinv = MT2Analysis<MT2EstimateTree>::readFromFile(cfg.getEventYieldDir() + "/analyses.root", "ZJets");
  if( Zinv==0 ) {
    std::cout << "-> Please run regionEventYields on MC first. I need to get the Z->vv MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }

  MT2Analysis<MT2Estimate>* ZZllRatioMC = new MT2Analysis<MT2Estimate>( "ZZllRatioMC", cfg.regionsSet() );
  //(*ZZllRatioMC) = ( (* (MT2Analysis<MT2Estimate>*)Zinv) / (*zllCR_mc) );
  (*ZZllRatioMC) = ( (* (MT2Analysis<MT2Estimate>*)Zinv) / (*zllCR_mc_integral) );

  MT2Analysis<MT2Estimate>* ZZllRatio = new MT2Analysis<MT2Estimate>( "ZZllRatio", cfg.regionsSet() );
  (*ZZllRatio) = (*ZZllRatioMC);

  MT2Analysis<MT2Estimate>* ZinvEstimateFromZll = new MT2Analysis<MT2Estimate>( "ZinvEstimateFromZll", cfg.regionsSet() );
  //(*ZinvEstimateFromZll) = (*zllCR) * (*ZZllRatio);
  (*ZinvEstimateFromZll) = (*zllCR_integral)*(*ZZllRatio);


  //MT2Analysis<MT2EstimateSyst>* ZinvEstimateFromGamma = MT2EstimateSyst::makeAnalysisFromEstimate( "ZinvEstimateFromGamma", regionsSet, gammaCR_times_ZgammaRatio );

  MT2Analysis<MT2Estimate>* ZinvEstimate = combineDataAndMC( ZinvEstimateFromZll, (MT2Analysis<MT2Estimate>*)Zinv );

  std::string outFile = cfg.getEventYieldDir() + "/zinvFromZll";
  outFile += ".root";

  ZinvEstimate->writeToFile( outFile, "recreate" );
  ZZllRatio->addToFile( outFile );
  Zinv->setName("Zinv");
  Zinv->addToFile( outFile );
  zllCR->addToFile( outFile );

  return 0;

}



MT2Analysis<MT2Estimate>* combineDataAndMC( MT2Analysis<MT2Estimate>* data, MT2Analysis<MT2Estimate>* mc ) {

  std::string dataname = data->getName();
  std::string mcname = mc->getName();

  // temporarily set all names to the output name so that returned MT2Analysis has consistent naming in all regions:
  std::string estimateName = "ZinvEstimate";
  data->setName( estimateName );
  mc->setName( estimateName );

  std::set<MT2Region> regions = data->getRegions();

  std::set<MT2Estimate*> newData;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* dataEst = data->get(*iR);
    MT2Estimate* mcEst = mc->get(*iR);

    MT2Estimate* thisNewEstimate;
    if( iR->nBJetsMin()>1 ) {
      thisNewEstimate =  new MT2Estimate(*mcEst);
    } else {
      thisNewEstimate =  new MT2Estimate(*dataEst);
    }
    newData.insert( thisNewEstimate );

  }

  MT2Analysis<MT2Estimate>* analysis = new MT2Analysis<MT2Estimate>( estimateName, newData );

  // set names back to original:
  data->setName( dataname );
  mc->setName( mcname );


  return analysis;

}

















































