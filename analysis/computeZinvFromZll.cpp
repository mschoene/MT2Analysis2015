#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>


#include "../interface/MT2Analysis.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2EstimateTree.h"
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





//float lumi = 0.1;
float lumi = 20.;

MT2Analysis<MT2Estimate>* combineDataAndMC( MT2Analysis<MT2Estimate>* data, MT2Analysis<MT2Estimate>* mc );


MT2Analysis<MT2Estimate>* getInclusiveRatioMC( const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* Zinv, MT2Analysis<MT2EstimateTree>* gammaCRtree );



int main(int argc, char* argv[]) {


   std::string regionsSet = "zurich";
  //std::string regionsSet = "13TeV_inclusive";
 if( argc>2 ) {
    regionsSet = std::string(argv[2]);
  }

  std::cout << "-> Using regions: " << regionsSet << std::endl;

  
  if( argc<2 ) {
    std::cout << "USAGE: ./regionEventYields [configFileName] regionSet" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  

  std::string configFileName(argv[1]);


  MT2Config cfg("cfgs/" + configFileName + ".txt");


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
  //  std::string samplesFileName = "../samples/samples_samples_PHYS14_skimprune.dat";

  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;



  std::string outputdir = "Zll_CR_" + cfg.mcSamples() + "_" + regionsSet;
  //  std::string outputdir = "Zll_" + configFileName;
    double intpart;
    double fracpart = modf(lumi, &intpart);
    std::string suffix;
    if( fracpart>0. )
      suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
    else
      suffix = std::string( Form("_%.0ffb", intpart ) );
    outputdir += suffix;
  
  system(Form("mkdir -p %s", outputdir.c_str()));
 

  MT2Analysis<MT2Estimate>* Zll_data = MT2Analysis<MT2Estimate>::readFromFile(Form("%s/data.root",outputdir.c_str()), "Zll");

  MT2Analysis<MT2Estimate>* Zll_mc = MT2Analysis<MT2Estimate>::readFromFile(Form("%s/mc.root",outputdir.c_str()), "Zll");



  MT2Analysis<MT2EstimateTree>* Zinv = MT2Analysis<MT2EstimateTree>::readFromFile(Form("EventYields_mc_PHYS14_v5_dummy_%.0ffb/analyses.root", lumi), "ZJets");
  if( Zinv==0 ) {
    std::cout << "-> Please run regionEventYields on MC first. I need to get the Z->vv MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }




 
  MT2Analysis<MT2Estimate>* ZllRatioMC = new MT2Analysis<MT2Estimate>( "ZllRatioMC", regionsSet );
  (*ZllRatioMC) = ( (* (MT2Analysis<MT2Estimate>*)Zinv) / (*Zll_mc) );
  


  MT2Analysis<MT2Estimate>* ZllRatio = new MT2Analysis<MT2Estimate>( "ZllRatio", regionsSet );
  (*ZllRatio) = (*ZllRatioMC);


  //"Data" estimate from the rounded data times the mc ratio of the SR/CR
  MT2Analysis<MT2Estimate>* ZinvEstimateFromZll = new MT2Analysis<MT2Estimate>( "ZinvEstimateFromZll", regionsSet );
  (*ZinvEstimateFromZll) = (*Zll_data) * (*ZllRatio);


  MT2Analysis<MT2Estimate>* ZinvEstimate = combineDataAndMC( ZinvEstimateFromZll, (MT2Analysis<MT2Estimate>*)Zinv );



  std::string outFile = outputdir + "/MT2ZinvEstimate.root";


  ZinvEstimate->writeToFile( outFile );
  ZllRatio->addToFile( outFile );
  Zinv->setName("Zinv");

  Zinv->addToFile( outFile );

  Zll_data->addToFile( outFile );

 





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

















































