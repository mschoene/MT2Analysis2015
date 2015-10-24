#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TVirtualFitter.h"
#include "TLegend.h"


#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateQCD.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"




bool closureTest = false;

float scaleMC = 0.; // simulate stats for the given lumi [0 means MC stats]

void projectFromInclusive( MT2Analysis<MT2Estimate>* analysis, MT2Analysis<MT2EstimateTree>* ana_inclusive, const std::string& selection );
void fillFromTreeAndRatio( MT2Estimate* estimate, MT2Estimate* nCR, MT2Estimate* r_effective, TTree* tree, TF1* f1_ratio, TH1D* h_band );
void get_rHat( MT2Analysis<MT2Estimate>* rHat  , MT2Analysis<MT2EstimateTree>* analysis, MT2Analysis<MT2EstimateTree>* ana_rest=NULL );
void get_fJets( MT2Analysis<MT2Estimate>* fJets, MT2Analysis<MT2EstimateTree>* analysis, MT2Analysis<MT2EstimateTree>* ana_rest=NULL );
void drawSingleFit( const MT2Config& cfg, bool useMC, const std::string& outdir, MT2EstimateQCD* qcd, MT2EstimateQCD* all, TF1* thisFitQCD, TH1D* h_band, float xMin_fit, float xMax_fit );
void computePurity( TH1D* purity, TH1D* nonQCD, TH1D* all );
void multiplyHisto( TH1D* histo, TH1D* other );
void scaleHisto( TH1D* histo, float val, float err );
void drawClosure( const std::string& outputdir, MT2Analysis<MT2Estimate>* estimate, MT2Analysis<MT2Estimate>* mcTruth );



int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "----------------------------------------------------" << std::endl;
  std::cout << "|                                                  |" << std::endl;
  std::cout << "|                                                  |" << std::endl;
  std::cout << "|          Running computeQCDFromDeltaPhi          |" << std::endl;
  std::cout << "|                                                  |" << std::endl;
  std::cout << "|                                                  |" << std::endl;
  std::cout << "----------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./fitDeltaPhiQCD [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  MT2DrawTools::setStyle();


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool useMC = true;

  if( argc>2 ) {

    std::string mc_or_data = std::string(argv[2]); 
    if( mc_or_data=="mc" ) mc_or_data="MC";
    if( mc_or_data=="MC" ) useMC = true;
    else useMC=false;

  } 

  TH1D::AddDirectory(kTRUE);



  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion/";
  std::string outputdir = qcdCRdir;
  if( closureTest ) {
    outputdir = outputdir + "/test/";
  }
  std::string fitsDir = outputdir;
  if( useMC ) fitsDir = fitsDir + "/fitsMC";
  else        fitsDir = fitsDir + "/fitsData";
  system( Form("mkdir -p %s", fitsDir.c_str() ));



  MT2Analysis<MT2EstimateTree>* qcdTree_mc   = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/mc.root",   "qcdCRtree" );
  MT2Analysis<MT2EstimateTree>* qcdTree_data = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/data.root", "qcdCRtree" );
  

  std::string regionsSet = cfg.regionsSet();
  std::string regionsSet_fJets = "zurich_onlyHT";
  std::string regionsSet_rHat  = "zurich_onlyJets_noB";

  std::cout << "-> Making MT2EstimateTrees from inclusive tree...";
  MT2Analysis<MT2EstimateTree>* qcd_4rHat ;
  MT2Analysis<MT2EstimateTree>* qcd_4fJets;
  MT2Analysis<MT2EstimateTree>* qcd_4rHat_rest ;
  MT2Analysis<MT2EstimateTree>* qcd_4fJets_rest;
  if( useMC ) {
    qcd_4rHat  = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mc_4rHat"   , regionsSet_rHat  , qcdTree_mc, "id>=153 && id<200 && mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3" ); // invert deltaPhi
    qcd_4fJets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mc_4fJets"  , regionsSet_fJets , qcdTree_mc, "id>=153 && id<200 && mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3" ); // invert deltaPhi 
  } else {
    qcd_4rHat  = MT2EstimateTree::makeAnalysisFromInclusiveTree( "data_4rHat" , regionsSet_rHat  , qcdTree_data, "id==1 && ht>1000. &&  mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3" ); // invert deltaPhi
    qcd_4fJets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "data_4fJets", regionsSet_fJets , qcdTree_data, "id==1 &&              mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3" ); // invert deltaPhi 
    qcd_4rHat_rest  = MT2EstimateTree::makeAnalysisFromInclusiveTree( "rest_4rHat" , regionsSet_rHat  , qcdTree_mc, "id>=300 && ht>1000. &&  mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3" ); // invert deltaPhi
    qcd_4fJets_rest = MT2EstimateTree::makeAnalysisFromInclusiveTree( "rest_4fJets", regionsSet_fJets , qcdTree_mc, "id>=300 &&              mt2>100. && mt2<200. && nJets>1 && deltaPhiMin<0.3" ); // invert deltaPhi 
  }

  //MT2Analysis<MT2EstimateTree>* mc_4rHat    = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mc_4rHat"   , regionsSet_rHat  , qcdTree_mc, "id>=153 && id<200 && mt2>100. && mt2<200. && deltaPhiMin<0.3" ); // invert deltaPhi
  //MT2Analysis<MT2EstimateTree>* mc_4fJets   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mc_4fJets"  , regionsSet_fJets , qcdTree_mc, "id>=153 && id<200 && mt2>100. && mt2<200. && deltaPhiMin<0.3" ); // invert deltaPhi 

  std::cout << " Done." << std::endl;


  std::cout << "-> Creating the MT2Estimates...";
  MT2Analysis<MT2Estimate>* estimate     = new MT2Analysis<MT2Estimate>("qcdEstimate", regionsSet);
  MT2Analysis<MT2Estimate>* nCR          = new MT2Analysis<MT2Estimate>("nCR"        , regionsSet);
  MT2Analysis<MT2Estimate>* r_effective  = new MT2Analysis<MT2Estimate>("r_effective", regionsSet);

  MT2Analysis<MT2Estimate>* est_mcRest   = new MT2Analysis<MT2Estimate>("est_mcRest"  , regionsSet);
  MT2Analysis<MT2Estimate>* nCR_mcRest   = new MT2Analysis<MT2Estimate>("nCR_mcRest"  , regionsSet);
  MT2Analysis<MT2Estimate>* r_eff_mcRest = new MT2Analysis<MT2Estimate>("r_eff_mcRest", regionsSet);

  MT2Analysis<MT2Estimate>* qcdPurity    = new MT2Analysis<MT2Estimate>("qcdPurity"  , regionsSet);

  MT2Analysis<MT2Estimate>* r_hat        = new MT2Analysis<MT2Estimate>("r_hat"      , regionsSet_rHat);
  MT2Analysis<MT2Estimate>* f_jets       = new MT2Analysis<MT2Estimate>("f_jets"     , regionsSet_fJets);
  std::cout << " Done." << std::endl;


  if( closureTest ) {

    MT2Estimate::rebinYields( estimate    , 1, 150., 200. );
    MT2Estimate::rebinYields( nCR         , 1, 150., 200. );
    MT2Estimate::rebinYields( r_effective , 1, 150., 200. );

    MT2Estimate::rebinYields( est_mcRest  , 1, 150., 200. );
    MT2Estimate::rebinYields( nCR_mcRest  , 1, 150., 200. );
    MT2Estimate::rebinYields( r_eff_mcRest, 1, 150., 200. );

    MT2Estimate::rebinYields( qcdPurity   , 1, 150., 200. );

  }

  std::cout << "-> Getting fJets...";
  if ( useMC ) 
    get_fJets( f_jets, qcd_4fJets );
  else
    get_fJets( f_jets, qcd_4fJets, qcd_4fJets_rest );    
  std::cout << " Done." << std::endl;
  std::cout << "-> Getting rHat...";
  if ( useMC ) 
    get_rHat ( r_hat , qcd_4rHat );
  else
    get_rHat ( r_hat , qcd_4rHat, qcd_4rHat_rest );
  std::cout << " Done." << std::endl;



  std::cout << "-> Making MT2EstimateQCD from inclusive tree...";
  MT2Analysis<MT2EstimateQCD>* est_all;
  if( useMC ) {
    est_all  = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "est", cfg.qcdRegionsSet(), qcdTree_mc  , "" );
  } else {
    est_all  = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "est", cfg.qcdRegionsSet(), qcdTree_data, "", "id==1" ); //use HT-only triggers for dphi-ratio
  }
  MT2Analysis<MT2EstimateQCD>* mc_rest = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "mc_rest", cfg.qcdRegionsSet(), qcdTree_mc  , "id>=300" );
  //MT2Analysis<MT2EstimateQCD>* data    = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "data"   , cfg.qcdRegionsSet(), qcdTree_data, "(ht>1000. && id==1) || (ht>450 && ht<1000. && id==2)" );

  MT2Analysis<MT2EstimateQCD>* est_minus_nonQCD = new MT2Analysis<MT2EstimateQCD>("est_minus_nonQCD"  , cfg.qcdRegionsSet() );
  //MT2Analysis<MT2EstimateQCD>* data_minus_nonQCD  = new MT2Analysis<MT2EstimateQCD>("data_minus_nonQCD", cfg.qcdRegionsSet() );
  std::cout << " Done." << std::endl;





  std::set<MT2Region> regions = estimate->getRegions();


  std::string lastRegionName = "";

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* this_estimate    = estimate    ->get( *iR );
    MT2Estimate* this_nCR         = nCR         ->get( *iR );
    MT2Estimate* this_r_effective = r_effective ->get( *iR );

    MT2Estimate* this_est_mcRest   = est_mcRest  ->get( *iR );
    MT2Estimate* this_nCR_mcRest   = nCR_mcRest  ->get( *iR );
    MT2Estimate* this_r_eff_mcRest = r_eff_mcRest->get( *iR );

    MT2Estimate* this_qcdPurity   = qcdPurity   ->get( *iR );

    MT2Region* regionToMatch;
    if( iR->nBJetsMin()==3 && iR->nJetsMin()==2 ) 
      regionToMatch = new MT2Region( iR->htMin(), iR->htMax(), 4, 6, iR->nBJetsMin(), iR->nBJetsMax() );
    else
      regionToMatch = new MT2Region( *iR );

    MT2Estimate* this_r_hat  = r_hat ->getWithMatch( *regionToMatch );
    MT2Estimate* this_f_jets = f_jets->getWithMatch( *regionToMatch );


    MT2EstimateQCD* matchedEstimate      = est_all          ->getWithMatch( *iR );
    //MT2EstimateQCD* matchedEstimate_qcd  = (useMC) ? mc_qcd ->getWithMatch( *iR ) : 0;
    MT2EstimateQCD* matchedEstimate_rest = mc_rest          ->getWithMatch( *iR );
    MT2EstimateQCD* matchedEstimate_mnQ  = est_minus_nonQCD ->getWithMatch( *iR );


    float lumiScale = 1.;
    if( !useMC || scaleMC!=0. ) {
      lumiScale = useMC ? scaleMC : cfg.lumi();
      if     ( iR->htMin() < 300. ) lumiScale /= 7000.; // prescale
      else if( iR->htMin() < 500. ) lumiScale /=  180.; // prescale
      else if( iR->htMin() < 600. ) lumiScale /=   60.;
    }


    matchedEstimate_mnQ->lDphi = (TH1D*) matchedEstimate->lDphi->Clone(matchedEstimate_mnQ->lDphi->GetName());
    if (useMC && scaleMC!=0.)  matchedEstimate_mnQ->lDphi->Scale(lumiScale);
    matchedEstimate_mnQ->lDphi->Add( matchedEstimate_rest->lDphi, -lumiScale );

    matchedEstimate_mnQ->hDphi = (TH1D*) matchedEstimate->hDphi->Clone(matchedEstimate_mnQ->hDphi->GetName());
    if (useMC && scaleMC!=0.)  matchedEstimate_mnQ->hDphi->Scale(lumiScale);
    matchedEstimate_mnQ->hDphi->Add( matchedEstimate_rest->hDphi, -lumiScale );

    if (useMC && scaleMC!=0.)  matchedEstimate_mnQ->sqrtErrors();

    float xMin_fit = (iR->htMin()>=1000.) ? 70. : 60.;
    float xMax_fit = 100.;
    TF1* f1_ratio = matchedEstimate_mnQ->getFit( "pow", xMin_fit, xMax_fit );
    TH1D* h_band = new TH1D("band", "", 500, matchedEstimate->lDphi->GetXaxis()->GetXmin(), matchedEstimate->lDphi->GetXaxis()->GetXmax());
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(h_band, 0.68);

    //std::cout << iR->htRegion()->getNiceName() << std::endl
    //          << "par0 : " << f1_ratio->GetParameter(0) << " +/- " << f1_ratio->GetParError(0) << std::endl
    //	        << "par1 : " << f1_ratio->GetParameter(1) << " +/- " << f1_ratio->GetParError(1) << std::endl;

    if( matchedEstimate->region->getName()!=lastRegionName ) // draw only one per HT region
      drawSingleFit( cfg, useMC, fitsDir, matchedEstimate_mnQ, matchedEstimate, f1_ratio, h_band, xMin_fit, xMax_fit );
    lastRegionName = matchedEstimate->region->getName();


    fillFromTreeAndRatio( this_estimate  , this_nCR       , this_r_effective , matchedEstimate->tree     , f1_ratio, h_band );
    fillFromTreeAndRatio( this_est_mcRest, this_nCR_mcRest, this_r_eff_mcRest, matchedEstimate_rest->tree, f1_ratio, h_band );
    //fillFromTreeAndRatio( this_estimate, this_nCR, this_r_effective, matchedEstimate_qcd->tree, f1_ratio, h_band );

    
    computePurity( this_qcdPurity->yield, this_est_mcRest->yield, this_estimate->yield );
    multiplyHisto( this_estimate->yield, this_qcdPurity->yield );
    multiplyHisto( this_r_effective->yield, this_qcdPurity->yield );

    int bin_bJets = this_r_hat ->yield->FindBin(iR->nBJetsMin());
    float thisRhatValue = this_r_hat->yield->GetBinContent( bin_bJets );
    float thisRhatError = this_r_hat->yield->GetBinError  ( bin_bJets );

    int bin_jets = this_f_jets->yield->FindBin(iR->nJetsMin() );
    float thisFjetsValue = this_f_jets->yield->GetBinContent( bin_jets );
    float thisFjetsError = this_f_jets->yield->GetBinError  ( bin_jets );

    if( iR->nBJetsMin()==3 && iR->nJetsMin()==2 ) {
      int bin_jets_2  = this_f_jets->yield->FindBin(4);
      thisFjetsValue += this_f_jets->yield->GetBinContent( bin_jets_2 );
      thisFjetsError *= thisFjetsError;
      thisFjetsError += this_f_jets->yield->GetBinError( bin_jets_2 )*this_f_jets->yield->GetBinError( bin_jets_2 );
      thisFjetsError  = sqrt(thisFjetsError);
    }

    scaleHisto( this_estimate->yield, thisRhatValue , thisRhatError  );
    scaleHisto( this_estimate->yield, thisFjetsValue, thisFjetsError );

    delete h_band;

  }  // for regions
      
    
  std::string outfileName = outputdir;
  if( useMC ) outfileName = outfileName + "/qcdEstimateMC.root";
  else        outfileName = outfileName + "/qcdEstimateData.root";

  estimate   ->writeToFile( outfileName, "recreate" );
  nCR        ->writeToFile( outfileName );
  r_effective->writeToFile( outfileName );
  r_hat      ->writeToFile( outfileName );
  f_jets     ->writeToFile( outfileName );
  qcdPurity  ->writeToFile( outfileName );




  std::string outputFile_fits = fitsDir;
  if( useMC ) outputFile_fits = outputFile_fits + "/qcdFitsMC.root";
  else        outputFile_fits = outputFile_fits + "/qcdFitsData.root";

  est_all          ->writeToFile( outputFile_fits, "recreate" );
  if( useMC ) {
    MT2Analysis<MT2EstimateQCD>* mc_qcd  = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "mc_qcd" , cfg.qcdRegionsSet(), qcdTree_mc  , "id>=100 && id<200" );
    mc_qcd         ->writeToFile( outputFile_fits );
  }
  mc_rest          ->writeToFile( outputFile_fits );
  est_minus_nonQCD ->writeToFile( outputFile_fits );



  return 0;

}





void projectFromInclusive( MT2Analysis<MT2Estimate>* analysis, MT2Analysis<MT2EstimateTree>* ana_inclusive, const std::string& selection ) {

  MT2EstimateTree* treeInclusive = ana_inclusive->get( MT2Region("HT450toInf_j2toInf_b0toInf") );
  if( treeInclusive==0 ) {
    std::cout << "ERROR!! You need to pass an inclusive MT2EstimateTree Analysis to use this function!" << std::endl;
    exit(19191);
  }

  std::set<MT2Region> regions = analysis->getRegions();


  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* thisEstimate= analysis->get( *iR );

    std::string fullSelection = iR->getRegionCuts();
    if( selection!="" ) fullSelection = fullSelection + " && " + selection;
    treeInclusive->tree->Project( thisEstimate->yield->GetName(), "mt2", Form("weight*(%s)", fullSelection.c_str()) );

  } // for regions


}



void get_rHat( MT2Analysis<MT2Estimate>* rHat, MT2Analysis<MT2EstimateTree>* analysis, MT2Analysis<MT2EstimateTree>* ana_rest ) {


  std::vector<float> uncert;
  uncert.push_back( 0.08 ); // from bruno
  uncert.push_back( 0.20 ); 
  uncert.push_back( 0.35 ); 
  uncert.push_back( 0.70 ); 

  std::set<MT2Region> regions = rHat->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* thisEst = rHat->get(*iR);

    std::string name(thisEst->yield->GetName());
    int nBins = 4;
    Double_t bins[nBins+1];
    bins[0] = 0.;
    bins[1] = 1.;
    bins[2] = 2.;
    bins[3] = 3.;
    bins[4] = 6.;

    delete thisEst->yield;
    thisEst->yield = new TH1D( name.c_str(), "", nBins, bins );
    thisEst->yield->Sumw2();

    MT2EstimateTree* thisTree = analysis->get(*iR);
    thisTree->tree->Project( name.c_str(), "nBJets", "weight" );
    MT2DrawTools::addOverflowSingleHisto(thisEst->yield);

    if ( ana_rest!=NULL ) {
      TH1D *toSubtract = new TH1D( Form("%s_rest",name.c_str()), "", nBins, bins );
      MT2EstimateTree* anotherTree = ana_rest->get(*iR);
      anotherTree->tree->Project( Form("%s_rest",name.c_str()), "nBJets", "weight" );
      MT2DrawTools::addOverflowSingleHisto(toSubtract);
      
      thisEst->yield->Add(toSubtract,-1);
    }

    for( int iBin=1; iBin<nBins+1; ++iBin ) {

      // add error in quadrature to estimate
      float val_est    = thisEst->yield->GetBinContent(iBin);
      float error_est  = thisEst->yield->GetBinError(iBin);
      float errorRel_est = error_est/val_est;

      float errorRel_tot = sqrt( errorRel_est*errorRel_est + uncert[iBin-1]*uncert[iBin-1] );
      thisEst->yield->SetBinError( iBin, errorRel_tot*val_est );

    }

    thisEst->yield->Scale( 1./thisEst->yield->Integral(1, nBins+1) );

  } // for regions

    
}


void get_fJets( MT2Analysis<MT2Estimate>* fJets, MT2Analysis<MT2EstimateTree>* analysis, MT2Analysis<MT2EstimateTree>* ana_rest ) {

  std::vector<float> uncert;
  uncert.push_back( 0.25 ); // from bruno
  uncert.push_back( 0.07 ); 
  uncert.push_back( 0.20 ); 

  std::set<MT2Region> regions = fJets->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* thisEst = fJets->get(*iR);

    std::string name(thisEst->yield->GetName());
    int nBins = 3;
    Double_t bins[nBins+1];
    bins[0] = 2.;
    bins[1] = 4.;
    bins[2] = 7.;
    bins[3] = 11.;

    delete thisEst->yield;
    thisEst->yield = new TH1D( name.c_str(), "", nBins, bins );
    thisEst->yield->Sumw2();

    MT2EstimateTree* thisTree = analysis->get(*iR);
    thisTree->tree->Project( name.c_str(), "nJets", "weight" );

    MT2DrawTools::addOverflowSingleHisto(thisEst->yield);

    if ( ana_rest!=NULL ) {
      TH1D *toSubtract = new TH1D( Form("%s_rest",name.c_str()), "", nBins, bins );
      MT2EstimateTree* anotherTree = ana_rest->get(*iR);
      anotherTree->tree->Project( Form("%s_rest",name.c_str()), "nJets", "weight" );
      MT2DrawTools::addOverflowSingleHisto(toSubtract);
      
      float ps = 1.;
      if      (iR->htMin() < 500.) ps = 180.;
      else if (iR->htMin() < 600.) ps =  60.;
      
      thisEst->yield->Add(toSubtract,-1/ps);
    }

    for( int iBin=1; iBin<nBins+1; ++iBin ) {

      // add error in quadrature to estimate
      float val_est    = thisEst->yield->GetBinContent(iBin);
      float error_est  = thisEst->yield->GetBinError(iBin);
      float errorRel_est = error_est/val_est;

      float errorRel_tot = sqrt( errorRel_est*errorRel_est + uncert[iBin-1]*uncert[iBin-1] );
      thisEst->yield->SetBinError( iBin, errorRel_tot*val_est );

    }
    

    thisEst->yield->Scale( 1./thisEst->yield->Integral(1, nBins+1) );

  } // for regions

}




void fillFromTreeAndRatio( MT2Estimate* estimate, MT2Estimate* nCR, MT2Estimate* r_effective, TTree* tree, TF1* f1_ratio, TH1D* h_band ) {


  int nBins;
  double* bins;
  estimate->region->getBins(nBins, bins);

  TProfile* hp_r = new TProfile( "r", "", nBins, bins );
  TProfile* hp_rErr = new TProfile( "rErr", "", nBins, bins );

  float weight;
  tree->SetBranchAddress( "weight", &weight );
  float mt2;
  tree->SetBranchAddress( "mt2", &mt2 );

  int nentries = tree->GetEntries();

  // this tree has only lowDeltaPhi events
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    float r = f1_ratio->Eval( mt2 );

    nCR     ->yield->Fill( mt2, weight   );
    estimate->yield->Fill( mt2, weight*r );

    hp_r->Fill( mt2, r, weight );
    hp_rErr->Fill( mt2, h_band->GetBinError(h_band->FindBin(mt2)), weight );

  } // for entries



  for( int iBin=1; iBin<nBins+1; ++iBin ) {

    float r = hp_r->GetBinContent(iBin);
    r_effective->yield->SetBinContent( iBin, hp_r   ->GetBinContent(iBin) );

    float error_fit  = hp_rErr->GetBinContent(iBin);

    // add fit error in quadrature to R
    float error_mean = hp_r   ->GetBinError  (iBin);
    float error_r    = sqrt(error_fit*error_fit + error_mean*error_mean);
    r_effective->yield->SetBinError(iBin, error_r);

    // add R error in quadrature to estimate
    float errorRel_r = error_r/r;
    float val_est    = estimate->yield->GetBinContent(iBin);
    float error_est  = estimate->yield->GetBinError(iBin);
    float errorRel_est = error_est/val_est;

    float errorRel_tot = sqrt( errorRel_est*errorRel_est + errorRel_r*errorRel_r );
    estimate->yield->SetBinError( iBin, errorRel_tot*val_est );

  }

  delete hp_r;
  delete hp_rErr;

}







void drawSingleFit( const MT2Config& cfg, bool useMC, const std::string& outdir, MT2EstimateQCD* qcd, MT2EstimateQCD* all, TF1* thisFitQCD, TH1D* h_band, float xMin_fit, float xMax_fit ) {


  TH1D* thisRatioAll = all->getRatio();
  TH1D* thisRatioQCD = qcd->getRatio();

  TCanvas* c1 = new TCanvas( "c2", "", 600, 600 );
  c1->cd();
  c1->SetLogx();
  c1->SetLogy();

  
  float xMin = thisRatioAll->GetXaxis()->GetXmin();
  float xMax = useMC ? thisRatioAll->GetXaxis()->GetXmax() : 200;  // we are blind to MT2>200 in data

  float yMax    = thisRatioAll->GetMaximum()*5.;
  float yMinAll = thisRatioQCD->GetMinimum()/2.;
  float yMin    = thisRatioQCD->GetMinimum()/2.;
  if( yMin < 0.03 ) yMin = 0.03;
  if( yMin > yMinAll && yMinAll>0.001 ) yMin = yMinAll;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax );
  h2_axes->SetXTitle( "M_{T2} [GeV]" );
  h2_axes->SetYTitle( "Ratio");
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->GetXaxis()->SetMoreLogLabels();
  h2_axes->GetYaxis()->SetNoExponent();
  h2_axes->Draw();

  std::vector<std::string> regionNiceNames = qcd->region->getNiceNames();

  TPaveText* regionName = new TPaveText( 0.4, 0.81, 0.9, 0.9, "brNDC" );
  regionName->SetTextAlign( 11 );
  regionName->SetTextSize( 0.035 );
  regionName->SetFillColor( 0 );
  regionName->AddText( regionNiceNames[0].c_str() );
  regionName->AddText( regionNiceNames[1].c_str() );
  //regionName->Draw("same");

  TLegend* legend = new TLegend( 0.4, 0.7, 0.8, 0.92, regionNiceNames[0].c_str() );
  //TLegend* legend = new TLegend( 0.4, 0.65, 0.8, 0.82 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  if( useMC ) {
    legend->AddEntry( thisRatioAll, "MC (all)", "P" );
    legend->AddEntry( thisRatioQCD, "MC (QCD Only)", "P" );
  } else {
    legend->AddEntry( thisRatioAll, "Data (all)", "P" );
    legend->AddEntry( thisRatioQCD, "Data (non-QCD subtracted)", "P" );
  }
  legend->AddEntry( thisFitQCD, "Fit", "L" );
  legend->Draw("same");

  TLine* lineLeft = new TLine( xMin_fit, yMin, xMin_fit, yMax );
  lineLeft->SetLineStyle(2);
  lineLeft->Draw("same");

  TLine* lineRight = new TLine( 100., yMin, 100., yMax );
  lineRight->SetLineStyle(2);
  lineRight->Draw("same");


  h_band->SetMarkerSize(0);
  h_band->SetFillColor(18); 
  h_band->SetFillStyle(3001);
  h_band->Draw("C E3 same");

  thisFitQCD->SetLineColor(46); 
  thisFitQCD->SetLineWidth(2); 
  thisFitQCD->Draw("L same");

  thisRatioAll->SetMarkerStyle(20);
  thisRatioAll->SetMarkerSize(1.3);
  thisRatioAll->SetLineColor(kBlack);
  thisRatioAll->Draw("P same");

  thisRatioQCD->SetMarkerStyle(24);
  thisRatioQCD->SetMarkerSize(1.3);
  thisRatioQCD->SetLineColor(kBlack);
  thisRatioQCD->Draw("P same");

  gPad->RedrawAxis();

  TPaveText* labelTop;
  if( useMC && scaleMC==0.) labelTop = MT2DrawTools::getLabelTopSimulation();
  else if( useMC )          labelTop = MT2DrawTools::getLabelTopSimulation(scaleMC);
  else        labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
  labelTop->Draw("same");

  c1->SaveAs( Form("%s/ratio_%s%s.pdf", outdir.c_str(), qcd->region->getName().c_str(), (useMC && scaleMC!=0) ? Form("_%.1ffb",scaleMC) : "" ));
  c1->SaveAs( Form("%s/ratio_%s%s.eps", outdir.c_str(), qcd->region->getName().c_str(), (useMC && scaleMC!=0) ? Form("_%.1ffb",scaleMC) : "" ));
  

  delete c1;
  delete h2_axes;

}





void computePurity( TH1D* purity, TH1D* nonQCD, TH1D* all ) {


  for( int iBin=1; iBin<nonQCD->GetXaxis()->GetNbins()+1; ++iBin ) {

    float nonQCD_val = nonQCD->GetBinContent(iBin);
    float nonQCD_err = nonQCD->GetBinError  (iBin);
    if( nonQCD_err>nonQCD_val ) nonQCD_err = nonQCD_val;

    float allCR_val  = all->GetBinContent(iBin);
    float allCR_err  = all->GetBinError  (iBin);
    if( allCR_err>allCR_val ) allCR_err = allCR_val;
    //float allCR_err  = sqrt( allCR_val );

    float qcdPurity = (allCR_val>0.) ? (allCR_val-nonQCD_val)/allCR_val : 0.;
    if( qcdPurity<0. ) qcdPurity = 0.;
    
    float qcdPurityErr = (qcdPurity>0.) ? sqrt( nonQCD_err*nonQCD_err/(allCR_val*allCR_val) + allCR_err*allCR_err*nonQCD_val*nonQCD_val/(allCR_val*allCR_val*allCR_val*allCR_val) ) : 0.;

    purity->SetBinContent( iBin, qcdPurity    );
    purity->SetBinError  ( iBin, qcdPurityErr );

  }

}



void multiplyHisto( TH1D* histo, TH1D* other ) {

  for( int iBin=1; iBin<histo->GetXaxis()->GetNbins()+1; ++iBin ) {

    float val1     = histo->GetBinContent(iBin);
    float val1_err = histo->GetBinError  (iBin);
    float val1_errRel = val1_err/val1;

    float val2     = other->GetBinContent(iBin);
    float val2_err = other->GetBinError  (iBin);
    float val2_errRel = val2_err/val2;

    float errRel = sqrt( val1_errRel*val1_errRel + val2_errRel*val2_errRel );

    histo->SetBinContent( iBin, val1*val2 );
    histo->SetBinError  ( iBin, errRel*val1*val2 );

  }

}



void scaleHisto( TH1D* histo, float val, float err ) {

  float relErr = (val>0.) ? err/val : 0.;

  for( int iBin=1; iBin<histo->GetXaxis()->GetNbins()+1; ++iBin ) {

    float thisVal = histo->GetBinContent(iBin);
    float newVal = thisVal*val;
    histo->SetBinContent( iBin, newVal );

    float thisError = histo->GetBinError(iBin);
    float thisRelErr = (thisVal>0.) ? thisError/thisVal : 0.;

    float newRelErr = sqrt( thisRelErr*thisRelErr + relErr*relErr );
    float newErr = newRelErr*newVal;

    histo->SetBinError( iBin, newErr );

  }


}


