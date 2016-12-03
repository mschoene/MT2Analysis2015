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
#include "../interface/MT2DrawTools.h"


bool do_dummyMC = false;


void getQCDMonojet( TCanvas* c1, const std::string& regionName, int &nCR, float &qcdFraction, float &qcdFractionError );



int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./drawQCDForMonojet [configFileName] [lumi/shape]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  MT2DrawTools::setStyle();

  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  bool shapeNorm = false;
  if( argc>2 ) {
    std::string normType(argv[2]);
    if( normType=="lumi" ) shapeNorm=false;
    else if( normType=="shape" ) shapeNorm=true;
    else {
      std::cout << "-> Only 'lumi' and 'shape' are supported normTypes." << std::endl;
      exit(17);
    }
  }


  if( shapeNorm )
    std::cout << "-> Using shape normalization." << std::endl;
  else
    std::cout << "-> Using lumi normalization." << std::endl;




  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion";

  MT2Analysis<MT2EstimateTree>* mcTree = MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/mc_forMonojet.root"  , "qcdCRtree");
  //MT2Analysis<MT2EstimateTree>* mcTree2= MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/mc_zinv.root"  , "qcdCRtree");

  std::cout << "-> Making analyses from inclusive tree..." << std::endl;
  MT2Analysis<MT2EstimateTree>* qcd   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "QCD"  , "13TeV_2016_inclusive", mcTree, "id>=100 && id<200" ); 
  std::cout << "    QCD done." << std::endl;
  MT2Analysis<MT2EstimateTree>* wjets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "WJets", "13TeV_2016_inclusive", mcTree, "id>=500 && id<600" ); 
  std::cout << "    WJets done." << std::endl;
  //MT2Analysis<MT2EstimateTree>* top   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "Top"  , "13TeV_2016_inclusive", mcTree, "id>=300 && id<500" ); 
  //std::cout << "    Top done." << std::endl;
  MT2Analysis<MT2EstimateTree>* zjets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "ZJets", "13TeV_2016_inclusive", mcTree, "id>=600 && id<700" ); 
  std::cout << "    ZJets done." << std::endl;



  MT2Analysis<MT2EstimateTree>* data;
  if( !do_dummyMC )
    data   = MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/data_forMonojet.root", "qcdCRtree");
  else{
    data   = MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/mc_forMonojet.root", "qcdCRtree");
    (*data) *= cfg.lumi(); //This is not exactly right
    //(*data) = cfg.lumi() * ((*qcd)+(*wjets)+(*zjets));
  }


  wjets->setFullName("W+Jets");
  wjets->setColor(kWJets);
  //wjets->setWeight(1.17);

  //top->setFullName("Top");
  //top->setColor(kTop);

  qcd->setFullName("QCD");
  qcd->setColor(kQCD);

  zjets->setFullName("Z+Jets");
  zjets->setColor(kZJets);

  std::vector< MT2Analysis<MT2EstimateTree>* > mc;
  mc.push_back(qcd);
  mc.push_back(wjets);
  mc.push_back(zjets);
  //mc.push_back(top);

  std::string plotsDir = qcdCRdir + "/plotsMonojet";
  if( shapeNorm ) plotsDir += "_shape";


  MT2DrawTools dt(plotsDir, cfg.lumi() );
  if( shapeNorm ) {
    dt.set_shapeNorm( shapeNorm );
    dt.set_lumiErr(0.);
  }


  dt.set_data( data );
  dt.set_mc( &mc );

  dt.set_addOverflow( false );
  dt.set_displaySF( false );
  //dt.set_mcSF( 1.3 );



  std::vector<TCanvas*> canvases;


  //std::string selection = "nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200.";
  std::string selection = "(id<100 || id>=151) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>250. && met>250.";
  //std::string selection = "(id<100 || id>152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200.";
  canvases = dt.drawRegionYields_fromTree( "jet2_pt" , "jet2_pt" , selection, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 250 GeV", "N(j) = 2" );

  float mcSF = MT2DrawTools::getDataMCSF( canvases[0] );
  //dt.set_mcSF( mcSF );



  MT2Analysis<MT2Estimate>* qcdMonojet = new MT2Analysis<MT2Estimate>( "monojet_qcdEstimate", cfg.regionsSet() );
  MT2Analysis<MT2Estimate>* nCRMonojet = new MT2Analysis<MT2Estimate>( "monojet_nCR"        , cfg.regionsSet() );
  MT2Analysis<MT2Estimate>* rMonojet   = new MT2Analysis<MT2Estimate>( "monojet_r"          , cfg.regionsSet() );



  std::set<MT2Region> regions = qcdMonojet->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nJetsMin()!=1 || iR->nJetsMax()!=1 ) continue; // only monojet

    int nBJets = iR->nBJetsMin();

    float ptMin = iR->htMin();
    float ptMax = iR->htMax();

    std::cout<< ptMin << "   " << ptMax << std::endl;

    MT2Estimate* thisEst = qcdMonojet->get( *iR );
    MT2Estimate* thisNCR = nCRMonojet->get( *iR );
    MT2Estimate* thisR   = rMonojet  ->get( *iR );


    TString fullSelection = fullSelection.Format("(%s && %s && jet1_pt>%f && jet1_pt<%f)", selection.c_str(), iR->sigRegion()->getBJetCuts().c_str(), ptMin, ptMax );
    TString fullSelection_end = fullSelection_end.Format("(%s && %s && jet1_pt>%f)", selection.c_str(), iR->sigRegion()->getBJetCuts().c_str(), ptMin ) ;

    std::cout << fullSelection << std::endl;

    std::string bJetsLabel = (nBJets==0) ? "b = 0" : "b #geq 1";


    
    if(ptMax<0)
      canvases = dt.drawRegionYields_fromTree( Form("jet2_pt_%s", iR->getName().c_str()), "jet2_pt", (std::string)fullSelection_end, 20, 0., 300., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", bJetsLabel );
    else 
      canvases = dt.drawRegionYields_fromTree( Form("jet2_pt_%s", iR->getName().c_str()), "jet2_pt", (std::string)fullSelection, 20, 0., 300., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", bJetsLabel );

    int nCR;
    float r, r_err;
    getQCDMonojet( canvases[0], iR->getName(), nCR, r, r_err );

    // should be one bin in "mt2" anyways:
    thisNCR->yield->SetBinContent( 1, nCR   );
    thisR  ->yield->SetBinContent( 1, r     );
    thisR  ->yield->SetBinError  ( 1, r_err );
    thisEst->yield->SetBinContent( 1, r*nCR );

    //for( int iBin=1; iBin<thisEst->yield->GetXaxis()->GetNbins(); ++iBin ) {

    //  float ptMin = thisEst->yield->GetXaxis()->GetBinLowEdge(iBin);
    //  float ptMax = thisEst->yield->GetXaxis()->GetBinLowEdge(iBin+1);
    //  std::string fullSelection(Form("%s && %s && jet1_pt>%f && jet1_pt<%f", selection.c_str(), iR->sigRegion()->getBJetCuts().c_str(), ptMin, ptMax ) );

    //  std::string bJetsLabel = (nBJets==0) ? "b = 0" : "b #geq 1";
    //  canvases = dt.drawRegionYields_fromTree( Form("jet2_pt_bin%d_b%d", iBin, nBJets) , "jet2_pt" , fullSelection, 20, 0., 300., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", bJetsLabel );

    //  int nCR;
    //  float r, r_err;
    //  getQCDMonojet( canvases[0], iR->getName(), nCR, r, r_err );

    //  thisNCR->yield->SetBinContent( iBin, nCR   );
    //  thisR  ->yield->SetBinContent( iBin, r     );
    //  thisR  ->yield->SetBinError  ( iBin, r_err );
    //  thisEst->yield->SetBinContent( iBin, r*nCR );

    //}

  }
    

  std::string outfileName( Form( "%s/qcdEstimateMonojet.root", qcdCRdir.c_str() ) );
  qcdMonojet->writeToFile( outfileName, "recreate" ); 
  nCRMonojet->writeToFile( outfileName );
  rMonojet  ->writeToFile( outfileName );
  

  return 0;

}






void getQCDMonojet( TCanvas* c1, const std::string& regionName, int &nCR, float &qcdFraction, float &qcdFractionError ) {

  TList* list = MT2DrawTools::getCorrectList( c1 );
  TGraphAsymmErrors* gr_data  = (TGraphAsymmErrors*)list->FindObject("Graph");

  nCR = MT2DrawTools::graphIntegral( gr_data, 30., 60. );


  // decided to take flat conservative QCD purity values
  TString regionName_tstr(regionName);
  if( regionName_tstr.Contains("b0") ) {
    qcdFraction = 0.3;
  } else {
    qcdFraction = 0.4;
  }
  qcdFractionError = qcdFraction*0.5;

//THStack* bgStack = (THStack*)list->FindObject( "bgStack" );

//TList* mcHists = bgStack->GetHists();

//TH1D* h1_qcd   = (TH1D*)mcHists->FindObject( "h1_QCD_HT200toInf_j1toInf_b0toInf" );
//TH1D* h1_wjets = (TH1D*)mcHists->FindObject( "h1_WJets_HT200toInf_j1toInf_b0toInf" );
//TH1D* h1_zjets = (TH1D*)mcHists->FindObject( "h1_ZJets_HT200toInf_j1toInf_b0toInf" );


//int bin30 = h1_qcd->FindBin(30);
//int bin60 = h1_qcd->FindBin(60);

//double qcd3060_err;
//double qcd3060 = h1_qcd->IntegralAndError( bin30, bin60, qcd3060_err );
//double wjets3060 = h1_wjets->Integral( bin30, bin60 );
//double zjets3060 = h1_zjets->Integral( bin30, bin60 );

//double mcTot = qcd3060 + wjets3060 + zjets3060;

//qcdFraction = qcd3060 / mcTot;
//qcdFractionError = qcd3060_err / mcTot;

}
