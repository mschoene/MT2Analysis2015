#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>



#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2EstimateSyst.h"

#define mt2_cxx
#include "../interface/mt2.h"


#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLorentzVector.h"




float lumi = 4.; // fb-1



int type = 1;

// 0: use bare MC for ratio (pure GJet)
// 1: use GJet+QCD and multiply by fitted purity




MT2Analysis<MT2EstimateSyst>* combineDataAndMC( MT2Analysis<MT2EstimateSyst>* data, MT2Analysis<MT2Estimate>* mc );


int main( int argc, char* argv[] ) {


  std::string samples = "PHYS14_v4_skimprune";
  if( argc>1 ) {
    std::string samples_tmp(argv[1]); 
    samples = samples_tmp;
  }


  std::string regionsSet = "zurich";

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string outputdir( Form("ZinvEstimateFromGamma_%s_%s_%.0ffb_type%d", samples.c_str(), regionsSet.c_str(), lumi, type) );
  system(Form("mkdir -p %s", outputdir.c_str()));


  std::string gammaControlRegionDir(Form("GammaControlRegion_%s_%s_%.0ffb", samples.c_str(), regionsSet.c_str(), lumi));

  MT2Analysis<MT2Estimate>* gammaCR = MT2Analysis<MT2Estimate>::readFromFile(gammaControlRegionDir + "/data.root", "gammaCR");
  MT2Analysis<MT2Estimate>* gamma_prompt = MT2Analysis<MT2Estimate>::readFromFile(gammaControlRegionDir + "/mc.root", "prompt");

  if( gammaCR==0 || gamma_prompt==0 ) {
    std::cout << "-> Please run gammaControlRegion first. I need to get the gammaCR yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(193);
  }

  

  MT2Analysis<MT2EstimateTree>* Zinv = MT2Analysis<MT2EstimateTree>::readFromFile(Form("EventYields_mc_PHYS14_v4_dummy_%.0ffb/analyses.root", lumi), "ZJets");
  if( Zinv==0 ) {
    std::cout << "-> Please run regionEventYields on MC first. I need to get the Z->vv MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }


  MT2Analysis<MT2Estimate>* ZgammaRatio = new MT2Analysis<MT2Estimate>( "ZgammaRatio", regionsSet );
  (*ZgammaRatio) = (* (MT2Analysis<MT2Estimate>*)Zinv) / (*gamma_prompt);


  MT2Analysis<MT2Estimate>* gammaCR_times_ZgammaRatio = new MT2Analysis<MT2Estimate>( "gammaCR_times_ZgammaRatio", regionsSet );
  if( type==0 )
    (*gammaCR_times_ZgammaRatio) = (*gamma_prompt) * (*ZgammaRatio);
  else
    (*gammaCR_times_ZgammaRatio) = (*gammaCR) * (*ZgammaRatio);


  //MT2Analysis<MT2EstimateSyst>* purity = 0;

  MT2Analysis<MT2EstimateSyst>* ZinvEstimateFromGamma = MT2EstimateSyst::makeAnalysisFromEstimate( "ZinvEstimateFromGamma", regionsSet, gammaCR_times_ZgammaRatio );

  if( type!=0 ) {

    MT2Analysis<MT2EstimateSyst>* purity = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsDataRC/purityFit.root", "purity" );

    if( purity==0 ) {
      std::cout << "-> Please run fitPurityGamma first. I need to get the purity from there." << std::endl;
      std::cout << "-> Thank you for your cooperation." << std::endl;
      exit(195);
    }
    ZinvEstimateFromGamma = new MT2Analysis<MT2EstimateSyst>( "ZinvEstimateFromGamma", regionsSet );
    (*ZinvEstimateFromGamma) = (*purity);
    (*ZinvEstimateFromGamma) *= 0.92; // flat f
    (*ZinvEstimateFromGamma) *= (*gammaCR_times_ZgammaRatio);

  }


  MT2Analysis<MT2EstimateSyst>* ZinvEstimate = combineDataAndMC( ZinvEstimateFromGamma, (MT2Analysis<MT2Estimate>*)Zinv );

  std::string outFile = outputdir + "/MT2ZinvEstimate.root";

  ZinvEstimate->writeToFile( outFile );
  ZgammaRatio->addToFile( outFile );
  //if( purity!=0 )
  //  purity->addToFile( outFile );
  Zinv->setName("Zinv");
  Zinv->addToFile( outFile );

  return 0;

}





MT2Analysis<MT2EstimateSyst>* combineDataAndMC( MT2Analysis<MT2EstimateSyst>* data, MT2Analysis<MT2Estimate>* mc ) {

  std::string dataname = data->getName();
  std::string mcname = mc->getName();

  // temporarily set all names to the output name so that returned MT2Analysis has consistent naming in all regions:
  std::string estimateName = "ZinvEstimate";
  data->setName( estimateName );
  mc->setName( estimateName );

  std::set<MT2Region> regions = data->getRegions();

  std::set<MT2EstimateSyst*> newData;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateSyst* dataEst = data->get(*iR);
    MT2Estimate* mcEst = mc->get(*iR);

    MT2EstimateSyst* thisNewEstimate;
    if( iR->nBJetsMin()>1 ) {
      thisNewEstimate =  new MT2EstimateSyst(*mcEst);
      for( int ibin=1; ibin<thisNewEstimate->yield->GetNbinsX()+1; ++ibin ) {
        thisNewEstimate->yield_systUp->SetBinContent( ibin, 2.*thisNewEstimate->yield->GetBinContent(ibin) );
        thisNewEstimate->yield_systDown->SetBinContent( ibin, 0. );
      }
    } else {
      thisNewEstimate =  new MT2EstimateSyst(*dataEst);
    }
    newData.insert( thisNewEstimate );

  }

  MT2Analysis<MT2EstimateSyst>* analysis = new MT2Analysis<MT2EstimateSyst>( estimateName, newData );

  // set names back to original:
  data->setName( dataname );
  mc->setName( mcname );


  return analysis;

}

