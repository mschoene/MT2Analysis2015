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


MT2Analysis<MT2Estimate>* getInclusiveRatioMC( const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* Zinv, MT2Analysis<MT2EstimateTree>* gammaCRtree );
MT2Analysis<MT2EstimateSyst>* combineDataAndMC( MT2Analysis<MT2EstimateSyst>* data, MT2Analysis<MT2Estimate>* mc );


int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeZinvFromZll               |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc!=2 ) {
    std::cout << "USAGE: ./computeZinvFromZll [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this




  std::string zllControlRegionDir = cfg.getEventYieldDir() + "/zllControlRegion"; //(Form("ZllControlRegion_%s_%s_%.0ffb", samples.c_str(), regionsSet.c_str(), lumi));

  MT2Analysis<MT2Estimate>* zllData = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/data_forZinvEst.root", "data");
  MT2Analysis<MT2Estimate>* zllMC_ = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_forZinvEst.root", "zllCR");

 

  
  MT2Analysis<MT2Estimate>* zllMC;
  if ( !use_extrapolation )  {
    zllMC = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_forZinvEst.root", "");
    (*zllMC) = (*zllMC) * cfg.lumi(); 
  }else{
    (*zllMC_) = (*zllMC_) * cfg.lumi();
  }

  MT2Analysis<MT2Estimate>* zllMC_integral;
  if( use_extrapolation )
    zllMC_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "zllMC_integral", cfg.regionsSet(), zllMC_ );

  if( zllData==0 || zllMC_==0  ) {
    std::cout << "-> Please run gammaControlRegion first. I need to get the zllData yields from there." << std::endl;
    std::cout << "->Thank you for your cooperation." << std::endl;
    exit(193);
  }

  

  MT2Analysis<MT2EstimateTree>* Zinv = MT2Analysis<MT2EstimateTree>::readFromFile(cfg.getEventYieldDir() + "/analyses.root", "ZJets");
  if( Zinv==0 ) {
    std::cout << "-> Please run regionEventYields on MC first. I need to get the Z->vv MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }
  (* (MT2Analysis<MT2Estimate>*) Zinv) = (* (MT2Analysis<MT2Estimate>*)Zinv) * cfg.lumi();


  //MT2Analysis<MT2EstimateTree>* gammaCRtree = MT2Analysis<MT2EstimateTree>::readFromFile(gammaControlRegionDir + "/data.root", "gammaCRtree");
  //MT2Analysis<MT2Estimate>* ZgammaRatioMC = getInclusiveRatioMC( regionsSet, Zinv, gammaCRtree );

  MT2Analysis<MT2Estimate>* ZinvZllRatioMC = new MT2Analysis<MT2Estimate>( "ZinvZllRatioMC", cfg.regionsSet() );
  if( !use_extrapolation )
    (*ZinvZllRatioMC) = ( (* (MT2Analysis<MT2Estimate>*)Zinv) / (*zllMC) );
  else
    (*ZinvZllRatioMC) = ( (* (MT2Analysis<MT2Estimate>*)Zinv) / (*zllMC_integral) );
  //  (*ZinvZllRatioMC) = ( (* (MT2Analysis<MT2Estimate>*)Zinv) / (*zllMC) );

  //MT2Analysis<MT2Estimate>* ZinvZllRatio = MT2EstimateSyst::makeAnalysisFromEstimate( "ZinvZllRatio", regionsSet, ZinvZllRatioMC );
  MT2Analysis<MT2Estimate>* ZinvZllRatio = new MT2Analysis<MT2Estimate>( "ZinvZllRatio", cfg.regionsSet() );
  (*ZinvZllRatio) = (*ZinvZllRatioMC);



  MT2Analysis<MT2EstimateSyst>* zll_est_ = MT2EstimateSyst::makeAnalysisFromEstimate( "zll_est_", cfg.regionsSet(), zllData );

  MT2Analysis<MT2EstimateSyst>* zll_est;
  if( !use_extrapolation ) {
    zll_est = MT2EstimateSyst::makeAnalysisFromEstimate( "zll_est", cfg.regionsSet(), zllData );
  }



  MT2Analysis<MT2Estimate>* TopMC_ = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_Top_forZinvEst.root", "zllCR");
  (*TopMC_) = (*TopMC_) * cfg.lumi();
  //  MT2Analysis<MT2Estimate>* TopMC_integral;
  //  TopMC_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "TopMC_integral", cfg.regionsSet(), TopMC_ );



  MT2Analysis<MT2Estimate>* diLep = new MT2Analysis<MT2Estimate>( "diLep", cfg.regionsSet() ); 
  (*diLep) = (*TopMC_) + (*zllMC_);
  
  MT2Analysis<MT2EstimateSyst>* purity = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", (MT2Analysis<MT2Estimate>*)zllMC_, (MT2Analysis<MT2Estimate>*)diLep);


  (*zll_est_) = (*zll_est_ ) *  (*purity);




  MT2Analysis<MT2EstimateSyst>* zll_est_integral;
  if( use_extrapolation ){
    zll_est_integral = MT2EstimateSyst::makeIntegralAnalysisFromEstimate( "zll_est_integral", cfg.regionsSet(), zll_est_ );
  }



 


  MT2Analysis<MT2EstimateSyst>* ZinvEstimateFromZll = new MT2Analysis<MT2EstimateSyst>( "ZinvEstimateFromZll", cfg.regionsSet() );
  if( !use_extrapolation )
    (*ZinvEstimateFromZll) = (*zll_est) * (*ZinvZllRatio);
  else
    (*ZinvEstimateFromZll) = (*zll_est_integral) * (*ZinvZllRatio);

  MT2Analysis<MT2EstimateSyst>* ZinvEstimate = combineDataAndMC( ZinvEstimateFromZll, (MT2Analysis<MT2Estimate>*)Zinv );

  std::string outFile = cfg.getEventYieldDir() + "/zinvFromZll";
  outFile += ".root";

  ZinvEstimate->writeToFile( outFile, "recreate" );
  ZinvZllRatio->addToFile( outFile );
  Zinv->setName("Zinv");
  Zinv->addToFile( outFile );
  purity->addToFile( outFile );

  if( !use_extrapolation ){
    zllMC->addToFile( outFile );
    zll_est->addToFile( outFile );
  }
  else{
    zllMC_integral->addToFile( outFile );
    zll_est_integral->addToFile( outFile );
  }

  return 0;

}






MT2Analysis<MT2Estimate>* getInclusiveRatioMC( const MT2Config& cfg, MT2Analysis<MT2EstimateTree>* Zinv, MT2Analysis<MT2EstimateTree>* gammaCRtree ) {


  std::cout << "THIS FUNCTION IS BUGGED!!! DID YOU FIX THE BUG???" << std::endl;

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this


  MT2Analysis<MT2Estimate>* inclusiveZinv  = new MT2Analysis<MT2Estimate>( "inclusiveZinv" , cfg.regionsSet() );
  MT2Analysis<MT2Estimate>* inclusiveGamma = new MT2Analysis<MT2Estimate>( "inclusiveGamma", cfg.regionsSet() );

  std::set<MT2Region> regions = Zinv->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()>1 ) continue;

    TTree* thisTree_Zinv  = Zinv->get(*iR)->tree;
    TTree* thisTree_gamma = gammaCRtree->get(*iR)->tree;


    // loop on all regions and fill all histograms:
    for( std::set<MT2Region>::iterator iR2 = regions.begin(); iR2!=regions.end(); ++iR2 ) {

      TH1D* thisZinv  = inclusiveZinv ->get(*iR2)->yield;
      TH1D* thisGamma = inclusiveGamma->get(*iR2)->yield;

      TH1D* tmp_gamma = new TH1D(*thisGamma); // so that it gets the same binning
      tmp_gamma->SetName("tmp_gamma");
      TH1D* tmp_Zinv  = new TH1D(*thisZinv); // so that it gets the same binning
      tmp_Zinv->SetName("tmp_Zinv");

      tmp_gamma->Fill( 500., 10000.);

      thisTree_Zinv ->Project( "tmp_Zinv" , "mt2", "weight" );
      thisTree_gamma->Project( "tmp_gamma", "mt2", "weight*(prompt>1.5)" );

      thisZinv ->Add( tmp_Zinv  );
      thisGamma->Add( tmp_gamma );

      delete tmp_gamma;
      delete tmp_Zinv;

    } // for regions 2

  } // for regions 


  MT2Analysis<MT2Estimate>* inclusiveRatio = new MT2Analysis<MT2Estimate>( "ZgammaRatioMC", cfg.regionsSet() );
  (*inclusiveRatio) = (*inclusiveZinv) / (*inclusiveGamma);


  return inclusiveRatio;

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
    if( (!(use_extrapolation) && iR->nBJetsMin()>1) || iR->nBJetsMin()>3 /*|| iR->nJetsMin()==1*/ ) {
    //    if( (!(use_extrapolation) && iR->nBJetsMin()>1) || iR->nBJetsMin()>2 /*|| iR->nJetsMin()==1*/ ) {
      
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

