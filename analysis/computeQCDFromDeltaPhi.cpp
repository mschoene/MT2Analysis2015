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







void projectFromInclusive( MT2Analysis<MT2Estimate>* analysis, MT2Analysis<MT2EstimateTree>* ana_inclusive, const std::string& selection );
void fillFromTreeAndRatio( MT2Estimate* estimate, MT2Estimate* nCR, MT2Estimate* r_effective, TTree* tree, TF1* f1_ratio, TH1D* h_band );
void get_rHat( MT2Analysis<MT2Estimate>* rHat, MT2Analysis<MT2EstimateTree>* analysis );
void get_fJets( MT2Analysis<MT2Estimate>* fJets, MT2Analysis<MT2EstimateTree>* analysis );
void subtractNonQCD( MT2EstimateQCD* dmnQ, MT2EstimateQCD* data, MT2EstimateQCD* rest );
void subtractNonQCDSingleHisto( TH1D* h1_dmnQ, TH1D* h1_data, TH1D* h1_rest, float prescale );
void drawSingleFit( const std::string& outdir, MT2EstimateQCD* qcd, MT2EstimateQCD* all, TF1* thisFitQCD, TH1D* h_band, float xMin_fit, float xMax_fit );
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
  std::string fitsDir = qcdCRdir + "/fits";
  system( Form("mkdir -p %s", fitsDir.c_str() ));



  MT2Analysis<MT2EstimateTree>* qcdTree_mc   = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/mc.root",   "qcdCRtree" );
  MT2Analysis<MT2EstimateTree>* qcdTree_data = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/data.root", "qcdCRtree" );
  

  std::string regionsSet = cfg.regionsSet();
  std::string regionsSet_fJets = "zurich_onlyHT";
  std::string regionsSet_rHat  = "zurich_onlyJets_noB";

  std::cout << "-> Making MT2EstimateTrees from inclusive tree...";
  MT2Analysis<MT2EstimateTree>* mcTruth     = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mcTruth"    , regionsSet       , qcdTree_mc, "id>=153 && id<200 && mt2>200. && deltaPhiMin>0.3" ); // signal region for mcTruth

  MT2Analysis<MT2EstimateTree>* qcd_4rHat ;
  MT2Analysis<MT2EstimateTree>* qcd_4fJets;
  if( useMC ) {
    qcd_4rHat  = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mc_4rHat"   , regionsSet_rHat  , qcdTree_mc, "id>=153 && id<200 && mt2>100. && mt2<200. && deltaPhiMin<0.3" ); // invert deltaPhi
    qcd_4fJets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mc_4fJets"  , regionsSet_fJets , qcdTree_mc, "id>=153 && id<200 && mt2>100. && mt2<200. && deltaPhiMin<0.3" ); // invert deltaPhi 
  } else {
    qcd_4rHat  = MT2EstimateTree::makeAnalysisFromInclusiveTree( "data_4rHat" , regionsSet_rHat  , qcdTree_data, "id==1 && ht>1000. &&  mt2>100. && mt2<200. && deltaPhiMin<0.3" ); // invert deltaPhi
    qcd_4fJets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "data_4fJets", regionsSet_fJets , qcdTree_data, "id==1 &&              mt2>100. && mt2<200. && deltaPhiMin<0.3" ); // invert deltaPhi 
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


  std::cout << "-> Getting fJets...";
  get_fJets( f_jets, qcd_4fJets );
  std::cout << " Done." << std::endl;
  std::cout << "-> Getting rHat...";
  get_rHat ( r_hat , qcd_4rHat );
  std::cout << " Done." << std::endl;



  std::cout << "-> Making MT2EstimateQCD from inclusive tree...";
  MT2Analysis<MT2EstimateQCD>* mc_all  = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "mc"     , cfg.qcdRegionsSet(), qcdTree_mc  , "" );
  MT2Analysis<MT2EstimateQCD>* mc_rest = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "mc_rest", cfg.qcdRegionsSet(), qcdTree_mc  , "id>=300" );
  MT2Analysis<MT2EstimateQCD>* mc_qcd  = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "mc_qcd" , cfg.qcdRegionsSet(), qcdTree_mc  , "id>=100 && id<200" );
  MT2Analysis<MT2EstimateQCD>* data    = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "data"   , cfg.qcdRegionsSet(), qcdTree_data, "(ht>1000. && id==1) || (ht>450 && ht<1000. && id==2)" );

  MT2Analysis<MT2EstimateQCD>* mc_minus_nonQCD    = new MT2Analysis<MT2EstimateQCD>("mc_minus_nonQCD"  , cfg.qcdRegionsSet() );
  MT2Analysis<MT2EstimateQCD>* data_minus_nonQCD  = new MT2Analysis<MT2EstimateQCD>("data_minus_nonQCD", cfg.qcdRegionsSet() );
  std::cout << " Done." << std::endl;





  std::set<MT2Region> regions = estimate->getRegions();


  std::string lastRegionName = "";

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->htMin()<400. ) continue;
    if( iR->nJetsMax()==1 ) continue;

std::cout << "region: " << iR->getName() << std::endl;
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


    MT2EstimateQCD* matchedEstimate_mc   = mc_all           ->getWithMatch( *iR );
    MT2EstimateQCD* matchedEstimate_qcd  = mc_qcd           ->getWithMatch( *iR );
    MT2EstimateQCD* matchedEstimate_rest = mc_rest          ->getWithMatch( *iR );
    MT2EstimateQCD* matchedEstimate_mmnQ = mc_minus_nonQCD  ->getWithMatch( *iR );
    MT2EstimateQCD* matchedEstimate_data = data             ->getWithMatch( *iR );
    MT2EstimateQCD* matchedEstimate_dmnQ = data_minus_nonQCD->getWithMatch( *iR );

    subtractNonQCD( matchedEstimate_mmnQ, matchedEstimate_mc  , matchedEstimate_rest );
    subtractNonQCD( matchedEstimate_dmnQ, matchedEstimate_data, matchedEstimate_rest );

    float xMin_fit = (iR->htMin()>=1000.) ? 70. : 60.;
    float xMax_fit = 100.;
    TF1* f1_ratio = matchedEstimate_qcd->getFit( "pow", xMin_fit, xMax_fit );
    TH1D* h_band = new TH1D("band", "", 500, matchedEstimate_qcd->yield->GetXaxis()->GetXmin(), matchedEstimate_qcd->yield->GetXaxis()->GetXmax());
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(h_band, 0.68);

    if( matchedEstimate_qcd->region->getName()!=lastRegionName ) // draw only one per HT region
      drawSingleFit( fitsDir, matchedEstimate_qcd, matchedEstimate_mc, f1_ratio, h_band, xMin_fit, xMax_fit );
    lastRegionName = matchedEstimate_qcd->region->getName();

std::cout << "k8" << std::endl;

    fillFromTreeAndRatio( this_estimate  , this_nCR       , this_r_effective , matchedEstimate_mc ->tree, f1_ratio, h_band );
    fillFromTreeAndRatio( this_est_mcRest, this_nCR_mcRest, this_r_eff_mcRest, matchedEstimate_rest->tree, f1_ratio, h_band );
    //fillFromTreeAndRatio( this_estimate, this_nCR, this_r_effective, matchedEstimate_qcd->tree, f1_ratio, h_band );

std::cout << "k9" << std::endl;
    
    computePurity( this_qcdPurity->yield, this_est_mcRest->yield, this_estimate->yield );
    multiplyHisto( this_estimate->yield, this_qcdPurity->yield );
    multiplyHisto( this_r_effective->yield, this_qcdPurity->yield );

std::cout << "k10" << std::endl;
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

std::cout << "k11" << std::endl;
    scaleHisto( this_estimate->yield, thisRhatValue , thisRhatError  );
    scaleHisto( this_estimate->yield, thisFjetsValue, thisFjetsError );

std::cout << "k12" << std::endl;
    delete h_band;

  }  // for regions
      
    
  estimate   ->writeToFile( outputdir + "/qcdEstimate.root", "recreate" );
  nCR        ->writeToFile( outputdir + "/qcdEstimate.root" );
  mcTruth    ->writeToFile( outputdir + "/qcdEstimate.root" );
  r_effective->writeToFile( outputdir + "/qcdEstimate.root" );
  r_hat      ->writeToFile( outputdir + "/qcdEstimate.root" );
  f_jets     ->writeToFile( outputdir + "/qcdEstimate.root" );
  qcdPurity  ->writeToFile( outputdir + "/qcdEstimate.root" );


  mcTruth->setColor(kQCD);
  estimate->setColor(kBlack);
  drawClosure( qcdCRdir, estimate, (MT2Analysis<MT2Estimate>*)mcTruth );


  mc_all           ->writeToFile( outputdir + "/mcFits.root", "recreate" );
  mc_qcd           ->writeToFile( outputdir + "/mcFits.root" );
  mc_rest          ->writeToFile( outputdir + "/mcFits.root" );
  mc_minus_nonQCD  ->writeToFile( outputdir + "/mcFits.root" );
  data             ->writeToFile( outputdir + "/mcFits.root" );
  data_minus_nonQCD->writeToFile( outputdir + "/mcFits.root" );



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



void get_rHat( MT2Analysis<MT2Estimate>* rHat, MT2Analysis<MT2EstimateTree>* analysis ) {

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
    bins[4] = 10.;

    delete thisEst->yield;
    thisEst->yield = new TH1D( name.c_str(), "", nBins, bins );
    thisEst->yield->Sumw2();

    MT2EstimateTree* thisTree = analysis->get(*iR);
    thisTree->tree->Project( name.c_str(), "nBJets", "weight" );

    thisEst->yield->Scale( 1./thisEst->yield->Integral(1, nBins+1) );

  } // for regions

    
}


void get_fJets( MT2Analysis<MT2Estimate>* fJets, MT2Analysis<MT2EstimateTree>* analysis ) {

  std::set<MT2Region> regions = fJets->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* thisEst = fJets->get(*iR);

    std::string name(thisEst->yield->GetName());
    int nBins = 3;
    Double_t bins[nBins+1];
    bins[0] = 2.;
    bins[1] = 4.;
    bins[2] = 7.;
    bins[3] = 25.;

    delete thisEst->yield;
    thisEst->yield = new TH1D( name.c_str(), "", nBins, bins );
    thisEst->yield->Sumw2();

    MT2EstimateTree* thisTree = analysis->get(*iR);
    thisTree->tree->Project( name.c_str(), "nJets", "weight" );

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
    float error_est  = estimate->yield->GetBinError(iBin);
    float errorRel_est = error_est/estimate->yield->GetBinContent(iBin);

    float errorRel_tot = sqrt( errorRel_est*errorRel_est + errorRel_r*errorRel_r );
    estimate->yield->SetBinError( iBin, errorRel_tot );

  }

  delete hp_r;
  delete hp_rErr;

}




void subtractNonQCD( MT2EstimateQCD* dmnQ, MT2EstimateQCD* data, MT2EstimateQCD* rest ) {

  MT2Region* region = data->region;

  float prescale = 1.;
  if( region->htMin() < 500. ) prescale = 180.;
  else if( region->htMin() < 600. ) prescale = 60.;

  subtractNonQCDSingleHisto( dmnQ->lDphi, data->lDphi, rest->lDphi, prescale );
  subtractNonQCDSingleHisto( dmnQ->hDphi, data->hDphi, rest->hDphi, prescale );

}


void subtractNonQCDSingleHisto( TH1D* h1_dmnQ, TH1D* h1_data, TH1D* h1_rest, float prescale ) {

  std::string oldName(h1_dmnQ->GetName());

  h1_dmnQ = new TH1D( *h1_data );
  h1_dmnQ->Add( h1_rest, -prescale );

  h1_dmnQ->SetName( oldName.c_str() );

}



void drawSingleFit( const std::string& outdir, MT2EstimateQCD* qcd, MT2EstimateQCD* all, TF1* thisFitQCD, TH1D* h_band, float xMin_fit, float xMax_fit ) {


std::cout << "mmm1" << std::endl;
  TH1D* thisRatioAll = all->getRatio();
  TH1D* thisRatioQCD = qcd->getRatio();

  TCanvas* c1 = new TCanvas( "c2", "", 600, 600 );
  c1->cd();
  c1->SetLogx();
  c1->SetLogy();

  
  float xMin = thisRatioAll->GetXaxis()->GetXmin();
  float xMax = thisRatioAll->GetXaxis()->GetXmax();

  float yMax    = thisRatioAll->GetMaximum()*5.;
  float yMinAll = thisRatioQCD->GetMinimum()/2.;
  float yMin    = thisRatioQCD->GetMinimum()/2.;
  if( yMin < 0.03 ) yMin = 0.03;
  if( yMin > yMinAll && yMinAll>0.001 ) yMin = yMinAll;

std::cout << "mmm2" << std::endl;
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
  regionName->Draw("same");

std::cout << "mmm3" << std::endl;
  TLegend* legend = new TLegend( 0.4, 0.65, 0.8, 0.82 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( thisRatioAll, "MC (all)", "P" );
  legend->AddEntry( thisRatioQCD, "MC (QCD Only)", "P" );
  legend->AddEntry( thisFitQCD, "Fit", "L" );
  legend->Draw("same");

  TLine* lineLeft = new TLine( xMin_fit, yMin, xMin_fit, yMax );
  lineLeft->SetLineStyle(2);
  lineLeft->Draw("same");

  TLine* lineRight = new TLine( 100., yMin, 100., yMax );
  lineRight->SetLineStyle(2);
  lineRight->Draw("same");

std::cout << "mmm4" << std::endl;

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

std::cout << "mmm5" << std::endl;
  gPad->RedrawAxis();

  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
  labelTop->Draw("same");

  c1->SaveAs( Form("%s/ratio_%s.pdf", outdir.c_str(), qcd->region->getName().c_str()) );
  c1->SaveAs( Form("%s/ratio_%s.eps", outdir.c_str(), qcd->region->getName().c_str()) );
  

std::cout << "mmm6" << std::endl;
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


void drawClosure( const std::string& outputdir, MT2Analysis<MT2Estimate>* estimate, MT2Analysis<MT2Estimate>* mcTruth ) {


  system(Form("mkdir -p %s", outputdir.c_str()));

  
  std::set<MT2Region> MT2Regions = estimate->getRegions();
  
  TH1D* h_estimate_tot = new TH1D("h_estimate_tot", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  h_estimate_tot->Sumw2();
  h_estimate_tot->GetYaxis()->SetTitle("Events");
  h_estimate_tot->SetMarkerStyle(20);
  h_estimate_tot->SetMarkerSize(1.6);
  h_estimate_tot->SetLineColor( estimate->getColor() );
  h_estimate_tot->SetMarkerColor( estimate->getColor() );
  
  TH1D* h_mcTruth_tot = new TH1D("h_mcTruth_tot", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  h_mcTruth_tot->Sumw2();
  h_mcTruth_tot->GetYaxis()->SetTitle("Events");
  h_mcTruth_tot->SetFillColor(0);
  h_mcTruth_tot->SetLineColor( mcTruth->getColor() );
  h_mcTruth_tot->SetMarkerColor( mcTruth->getColor() );
  h_mcTruth_tot->SetMarkerStyle(20);
  h_mcTruth_tot->SetMarkerSize(1.6);
  
  TH1D* hPull = new TH1D("hPull", "", 20, -5, 5);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Data Driven - MC)/#sigma");
  hPull->GetYaxis()->SetTitle("Events");
  
  
  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::string fullPath = outputdir;

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      TH1D* h_estimate = estimate->get(*iMT2)->yield;
      h_estimate->SetMarkerStyle(20);
      h_estimate->SetMarkerSize(1.6);
      h_estimate->SetLineColor( estimate->getColor() );
      h_estimate->SetMarkerColor( estimate->getColor() );


      int nBins = h_estimate->GetXaxis()->GetNbins();
      double err_estimate;
      double int_estimate = h_estimate->IntegralAndError(1, nBins+1, err_estimate);
 
      h_estimate_tot->SetBinContent(iRegion, int_estimate);
      h_estimate_tot->SetBinError(iRegion, err_estimate);
      

      h_estimate_tot->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
      
      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.15);
      pad1->Draw();
      pad1->cd();

      TH1D* h_mcTruth = mcTruth->get(*iMT2)->yield;

      h_mcTruth->SetLineColor( mcTruth->getColor() );
      h_mcTruth->SetMarkerColor( mcTruth->getColor() );
      h_mcTruth->SetMarkerStyle(20);
      h_mcTruth->SetMarkerSize(1.6);
      
   
      double err_int;
      double int_mcTruth = h_mcTruth->IntegralAndError(1, nBins+1, err_int);
      h_mcTruth_tot->SetBinContent(iRegion, h_mcTruth->Integral());
      h_mcTruth_tot->SetBinError(iRegion, err_int);
      h_mcTruth_tot->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );

      if(int_mcTruth>0)
        hPull->Fill((int_estimate-int_mcTruth)/sqrt(err_estimate*err_estimate+err_int*err_int));


      float xMin = h_estimate->GetXaxis()->GetXmin();
      float xMax = h_estimate->GetXaxis()->GetXmax();
      float yMax_1 = h_estimate->GetMaximum()*1.5;
      float yMax_2 = 1.2*(h_estimate->GetMaximum() + h_estimate->GetBinError(h_estimate->GetMaximumBin()));
      float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
      float yMax_3 = h_mcTruth->GetMaximum()*1.5;
      float yMax_4 = 1.2*(h_mcTruth->GetMaximum() + h_mcTruth->GetBinError(h_mcTruth->GetMaximumBin()));
      float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
      float yMax = (yMax1>yMax2) ? yMax1 : yMax2;

      TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
      h2_axes->SetXTitle("M_{T2} [GeV]");
      h2_axes->SetYTitle("Events");

      h2_axes->Draw();
 
      //std::vector<std::string> niceNames = iMT2->getNiceNames();
      for( unsigned i=0; i<niceNames.size(); ++i ) {

        float yMax = 0.9-(float)i*0.05;
        float yMin = yMax - 0.05;
        TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
        regionText->SetTextSize(0.035);
        regionText->SetTextFont(42);
        regionText->SetFillColor(0);
        regionText->SetTextAlign(11);
        regionText->AddText( niceNames[i].c_str() );
        regionText->Draw("same");

      }


      TLegend* legend = new TLegend( 0.6, 0.9-2.*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
      legend->AddEntry( h_estimate, "data-driven", "P" );
      legend->AddEntry( h_mcTruth, "MC QCD", "P" );

      legend->Draw("same");

      h_estimate->Draw("P same");
      h_mcTruth->Draw("P same");
      //      bgStack.Draw("histoE, same");

      TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
      labelTop->Draw("same");

      gPad->RedrawAxis();

      c1->cd();
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.1);
      pad2->Draw();
      pad2->cd();

      std::string thisName = Form("%s_ratio", h_estimate->GetName());
      TH1D* h_ratio = (TH1D*) h_estimate->Clone(thisName.c_str());
      h_ratio->Divide(h_mcTruth);
      h_ratio->SetStats(0);	    
      h_ratio->SetMarkerStyle(20);
      h_ratio->SetLineColor(1);
      //      h_ratio->SetMarkerSize(0.02);
      h_ratio->GetXaxis()->SetLabelSize(0.00);
      h_ratio->GetXaxis()->SetTickLength(0.09);
      h_ratio->GetYaxis()->SetNdivisions(5,5,0);
      h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
      h_ratio->GetYaxis()->SetTitleSize(0.17);
      h_ratio->GetYaxis()->SetTitleOffset(0.4);
      h_ratio->GetYaxis()->SetLabelSize(0.17);
      h_ratio->GetYaxis()->SetTitle("Ratio");
      
      

      h_ratio->SetLineWidth(2);
      //h_ratio->Draw("PE");
      
      TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, 0.0, 2.0 );
      
      TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
      lineCentral->SetLineColor(1);
      

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      h_ratio->Draw("pe,same");

      gPad->RedrawAxis();
      
      c1->cd();

      c1->SaveAs( Form("%s/closure_%s.eps", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/closure_%s.pdf", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;
      delete h2_axes;
      delete h2_axes_ratio;
      delete h_ratio;
      delete h_estimate;
      delete h_mcTruth;
      
      ++iRegion;

  } // for MT2 regions

  
  TCanvas* c2 = new TCanvas("c2", "", 1200, 600);
  c2->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  pad1->SetLogy();
  
  float yMax_1 = h_estimate_tot->GetMaximum();
  float yMax_2 = h_estimate_tot->GetMaximum() + h_estimate_tot->GetBinError(h_mcTruth_tot->GetMaximumBin());
  float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
  float yMax_3 = h_mcTruth_tot->GetMaximum();
  float yMax_4 = h_mcTruth_tot->GetMaximum() + h_mcTruth_tot->GetBinError(h_mcTruth_tot->GetMaximumBin());
  float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
  float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
  yMax*=20.;
  
  float yMin = 1e-1;
  for( int iBin=1; iBin<h_estimate_tot->GetXaxis()->GetNbins()+1; ++iBin ) {
    if( h_estimate_tot    ->GetBinContent(iBin)>0. && h_estimate_tot    ->GetBinContent(iBin)<yMin ) yMin = h_estimate_tot    ->GetBinContent(iBin);
    if( h_mcTruth_tot->GetBinContent(iBin)>0. && h_mcTruth_tot->GetBinContent(iBin)<yMin ) yMin = h_mcTruth_tot->GetBinContent(iBin);
  }
  yMin /= 3.;
  
  h_mcTruth_tot->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  h_mcTruth_tot->GetYaxis()->SetRangeUser(yMin, yMax);
  h_mcTruth_tot->GetXaxis()->LabelsOption("v");
  h_mcTruth_tot->Draw("PE");


  h_estimate_tot->SetMarkerStyle(20);
  h_estimate_tot->SetMarkerSize(1.6);
  h_estimate_tot->SetLineColor( estimate->getColor() );
  h_estimate_tot->SetMarkerColor( estimate->getColor() );

  h_estimate_tot->Draw("pe,same");

  TLegend* legend = new TLegend( 0.18, 0.7, 0.32, 0.82 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( h_estimate_tot, "data-driven", "PL" );
  legend->AddEntry( h_mcTruth_tot, "QCD MC", "PL" );

  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
  labelTop->Draw("same");
  
  TLine* lHT[3];
  for( int iHT=1; iHT < 4; iHT++ ){
    lHT[iHT-1] = new TLine(11*iHT, -3., 11*iHT, yMax );
    lHT[iHT-1]->SetLineColor(kBlack);
    lHT[iHT-1]->SetLineStyle(3);
    lHT[iHT-1]->SetLineWidth(2);

    lHT[iHT-1]->Draw("same");
  }

  int nHTRegions = 4;
  std::vector< std::string > htRegions;
  htRegions.push_back("low H_{T}");
  htRegions.push_back("medium H_{T}");
  htRegions.push_back("high H_{T}");
  htRegions.push_back("extreme H_{T}");
  
  TPaveText* htBox[nHTRegions];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){
    
    htBox[iHT] = new TPaveText(0.16+0.2*iHT, 0.9-0.06, 0.34+0.2*iHT, 0.9, "brNDC");
    htBox[iHT]->AddText( htRegions[iHT].c_str() );
    
    htBox[iHT]->SetBorderSize(0);
    htBox[iHT]->SetFillColor(kWhite);
    htBox[iHT]->SetTextSize(0.038);
    htBox[iHT]->SetTextAlign(21); // align centered
    htBox[iHT]->SetTextFont(62);
    htBox[iHT]->Draw("same");

  }

  gPad->RedrawAxis();
  
  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();

  std::string thisName = Form("%s_ratio", h_estimate_tot->GetName());
  TH1D* h_Ratio = (TH1D*) h_estimate_tot->Clone(thisName.c_str());
  for( int iBin=1; iBin<h_Ratio->GetXaxis()->GetNbins()+1; ++iBin ) {
    float mc = h_estimate_tot->GetBinContent(iBin);
    float mc_err = h_estimate_tot->GetBinError(iBin);
    float est = h_mcTruth_tot->GetBinContent(iBin);
    float est_err = h_mcTruth_tot->GetBinError(iBin);
    float denom = sqrt( mc_err*mc_err + est_err*est_err );
    if( denom!=0. && mc>0. ) {
      h_Ratio->SetBinContent( iBin, (mc-est)/denom );
      h_Ratio->SetBinError( iBin, 1. );
    }
  }
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->SetLineWidth(2);
  

  TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( 0, MT2Regions.size(), -3., 3.);
  h2_axes_ratio->SetYTitle("Pull");

  TLine* LineCentral = new TLine(0, 0, MT2Regions.size(), 0);
  LineCentral->SetLineColor(1);


  h2_axes_ratio->Draw("");
  LineCentral->Draw("same");
  h_Ratio->Draw("pe,same");

  for( int iHT=1; iHT < 4; iHT++ ){
    lHT[iHT-1]->Draw("same");
  }


  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/closure_allRegions_pull.pdf", outputdir.c_str()) );
  c2->SaveAs( Form("%s/closure_allRegions_pull.eps", outputdir.c_str()) );

  pad2->cd();
  pad2->Clear();

  delete h2_axes_ratio;
  h2_axes_ratio = MT2DrawTools::getRatioAxes( 0, MT2Regions.size(), 0., 2.);
  h2_axes_ratio->SetYTitle("Data / MC");
  h2_axes_ratio->Draw("");

  TLine* lineOne = new TLine(0, 1., MT2Regions.size(), 1.);
  lineOne->SetLineColor(1);
  lineOne->Draw("same");

  delete h_Ratio;
  h_Ratio = (TH1D*) h_estimate_tot->Clone(thisName.c_str());
  h_Ratio->Divide( h_mcTruth_tot );
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->SetLineWidth(2);

  h_Ratio->Draw("pe,same");

  for( int iHT=1; iHT < 4; iHT++ ){
    lHT[iHT-1]->Draw("same");
  }


  c2->SaveAs( Form("%s/closure_allRegions_ratio.pdf", outputdir.c_str()) );
  c2->SaveAs( Form("%s/closure_allRegions_ratio.eps", outputdir.c_str()) );


  gStyle->SetOptStat(1110);
  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  hPull->SetStats(1110);
  TF1* f1_gaus = new TF1("f1_pull", "gaus", -2., 2.);
  f1_gaus->SetLineColor(kRed);
  hPull->Fit( f1_gaus, "QRL" );
  TPaveText* fitPars = new TPaveText( 0.2, 0.7, 0.5, 0.9, "brNDC" );
  fitPars->SetTextSize(0.03);
  fitPars->SetTextAlign(11);
  fitPars->SetFillColor(0);
  fitPars->AddText("Gaussian Fit:");
  fitPars->AddText(Form("Mean : %.2f +/- %.2f", f1_gaus->GetParameter(1), f1_gaus->GetParError(1) ));
  fitPars->AddText(Form("Sigma: %.2f +/- %.2f", f1_gaus->GetParameter(2), f1_gaus->GetParError(2) ));
  fitPars->Draw("same");
  hPull->Draw("hist same");
  f1_gaus->Draw("l same");
  c3->SaveAs( Form("%s/closure_pull.pdf", outputdir.c_str()) );
  c3->SaveAs( Form("%s/closure_pull.eps", outputdir.c_str()) );

  //TFile* file_4snt = TFile::Open("qcd_histo.root", "recreate");
  //file_4snt->cd();
  //h_mcTruth_tot->Write();
  //file_4snt->Close();

  delete c2;
  delete c3;
  delete h2_axes_ratio;
  
}


