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



void getQCDMonojet( TCanvas* c1, const std::string& regionName, int &nCR, float &qcdFraction );



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

  MT2Analysis<MT2EstimateTree>* data   = MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/data.root", "qcdCRtree");

  MT2Analysis<MT2EstimateTree>* mcTree = MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/mc.root"  , "qcdCRtree");
  MT2Analysis<MT2EstimateTree>* mcTree2= MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/mc_zinv.root"  , "qcdCRtree");

  std::cout << "-> Making analyses from inclusive tree..." << std::endl;
  MT2Analysis<MT2EstimateTree>* qcd   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "QCD"  , "13TeV_inclusive", mcTree, "id>=100 && id<200" ); 
  std::cout << "    QCD done." << std::endl;
  MT2Analysis<MT2EstimateTree>* wjets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "WJets", "13TeV_inclusive", mcTree, "id>=500 && id<600" ); 
  std::cout << "    WJets done." << std::endl;
  //MT2Analysis<MT2EstimateTree>* top   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "Top"  , "13TeV_inclusive", mcTree, "id>=300 && id<500" ); 
  //std::cout << "    Top done." << std::endl;
  MT2Analysis<MT2EstimateTree>* zjets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "ZJets", "13TeV_inclusive", mcTree2, "id>=600 && id<700" ); 
  std::cout << "    ZJets done." << std::endl;

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



  std::vector<TCanvas*> canvases;


  std::string selection = "((id<100 && nJets==2) || (id>=100 && nJets<=2)) && deltaPhiMin<0.3 && jet1_pt>200. && met>200. ";
  canvases = dt.drawRegionYields_fromTree( "jet2_pt" , "jet2_pt" , selection, 20, 30., 330., "Subleading Jet p_{T}", "GeV" );

  float mcSF = MT2DrawTools::getDataMCSF( canvases[0] );
  dt.set_mcSF( mcSF );



  MT2Analysis<MT2Estimate>* qcdMonojet = new MT2Analysis<MT2Estimate>( "monojet_qcdEstimate", cfg.regionsSet() );
  MT2Analysis<MT2Estimate>* nCRMonojet = new MT2Analysis<MT2Estimate>( "monojet_nCR"        , cfg.regionsSet() );
  MT2Analysis<MT2Estimate>* rMonojet   = new MT2Analysis<MT2Estimate>( "monojet_r"          , cfg.regionsSet() );



  std::set<MT2Region> regions = qcdMonojet->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nJetsMin()!=1 || iR->nJetsMax()!=1 ) continue; // only monojet

    int nBJets = iR->nBJetsMin();

    MT2Estimate* thisEst = qcdMonojet->get( *iR );
    MT2Estimate* thisNCR = nCRMonojet->get( *iR );
    MT2Estimate* thisR   = rMonojet  ->get( *iR );

    for( int iBin=1; iBin<thisEst->yield->GetXaxis()->GetNbins(); ++iBin ) {

      float ptMin = thisEst->yield->GetXaxis()->GetBinLowEdge(iBin);
      float ptMax = thisEst->yield->GetXaxis()->GetBinLowEdge(iBin+1);
      std::string fullSelection(Form("%s && deltaPhiMin<0.3 && nJets==2 && jet1_pt>%f && jet1_pt<%f && met>200.", iR->sigRegion()->getBJetCuts().c_str(), ptMin, ptMax ) );

      std::string bJetsLabel = (nBJets==0) ? "b = 0" : "b #geq 1";
      canvases = dt.drawRegionYields_fromTree( Form("jet2_pt_bin%d_b%d", iBin, nBJets) , "jet2_pt" , fullSelection, 20, 0., 300., "Subleading Jet p_{T}", "GeV", "H_{T} > 200 GeV", bJetsLabel );

      int nCR;
      float r;
      getQCDMonojet( canvases[0], iR->getName(), nCR, r );

      thisNCR->yield->SetBinContent( iBin, nCR   );
      thisR  ->yield->SetBinContent( iBin, r     );
      thisEst->yield->SetBinContent( iBin, r*nCR );

    }

  }
    

  std::string outfileName( Form( "%s/qcdEstimateMonojet.root", qcdCRdir.c_str() ) );
  qcdMonojet->writeToFile( outfileName, "recreate" ); 
  nCRMonojet->writeToFile( outfileName );
  rMonojet  ->writeToFile( outfileName );
  

  return 0;

}






void getQCDMonojet( TCanvas* c1, const std::string& regionName, int &nCR, float &qcdFraction ) {

  TList* list = MT2DrawTools::getCorrectList( c1 );
  TGraphAsymmErrors* gr_data  = (TGraphAsymmErrors*)list->FindObject("Graph");

  nCR = MT2DrawTools::graphIntegral( gr_data, 30., 60. );

  THStack* bgStack = (THStack*)list->FindObject( "bgStack" );

  TList* mcHists = bgStack->GetHists();

  TH1D* h1_qcd   = (TH1D*)mcHists->FindObject( "h1_QCD_HT200toInf_j1toInf_b0toInf" );
  TH1D* h1_wjets = (TH1D*)mcHists->FindObject( "h1_WJets_HT200toInf_j1toInf_b0toInf" );
  TH1D* h1_zjets = (TH1D*)mcHists->FindObject( "h1_ZJets_HT200toInf_j1toInf_b0toInf" );


  int bin30 = h1_qcd->FindBin(30);
  int bin60 = h1_qcd->FindBin(60);

  float qcd3060   = h1_qcd  ->Integral( bin30, bin60 );
  float wjets3060 = h1_wjets->Integral( bin30, bin60 );
  float zjets3060 = h1_zjets->Integral( bin30, bin60 );
  float mcTot = qcd3060 + wjets3060 + zjets3060;

  qcdFraction = qcd3060 / mcTot;

}
