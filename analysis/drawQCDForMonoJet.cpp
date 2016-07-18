#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"



void getQCDMonojet( TCanvas* c1, const std::string& regionName, int &nCR, float &qcdFraction, float &qcdFractionError );

void fitQCDMonojet( TCanvas* c1, const std::string& regionName, int &nCR, float &qcdFraction, float &qcdFractionError, MT2Config cfg );



void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* data_of, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields_of, const std::string& saveName, const std::string& varName, const std::string& selection, const std::string& selection_of, int nBins, double* bins, std::string name = "", std::string name_of = "", std::string axisName="", const std::string& units="", const std::string& kinCuts="", const std::string& topoCuts="", bool drawData=0 );


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

  MT2Analysis<MT2EstimateTree>* data   = MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/data_forMonojet.root", "qcdCRtree");

  MT2Analysis<MT2EstimateTree>* mcTree = MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/mc_forMonojet.root"  , "qcdCRtree");
  //MT2Analysis<MT2EstimateTree>* mcTree2= MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir+"/mc_zinv.root"  , "qcdCRtree");

  std::cout << "-> Making analyses from inclusive tree..." << std::endl;
  MT2Analysis<MT2EstimateTree>* qcd   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "QCD"  , "13TeV_inclusive", mcTree, "id>=100 && id<200" ); 
  std::cout << "    QCD done." << std::endl;
  MT2Analysis<MT2EstimateTree>* wjets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "WJets", "13TeV_inclusive", mcTree, "id>=500 && id<600" ); 
  std::cout << "    WJets done." << std::endl;
  //MT2Analysis<MT2EstimateTree>* top   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "Top"  , "13TeV_inclusive", mcTree, "id>=300 && id<500" ); 
  //std::cout << "    Top done." << std::endl;
  MT2Analysis<MT2EstimateTree>* zjets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "ZJets", "13TeV_inclusive", mcTree, "id>=600 && id<700" ); 
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

  dt.set_addOverflow( false );
  // dt.set_displaySF( false );
  //dt.set_mcSF( 1.3 );



  std::string qcdCRdir_old = "/mnt/t3nfs01/data01/shome/mschoene/CMSSW_7_4_12_patch4/src/analysisCode/analysis/EventYields_data_Run2015_qcdMonojetTest/qcdControlRegion";

 MT2Analysis<MT2EstimateTree>* data_old   = MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir_old+"/data_forMonojet.root", "qcdCRtree");

  MT2Analysis<MT2EstimateTree>* mcTree_old = MT2Analysis<MT2EstimateTree>::readFromFile(qcdCRdir_old+"/mc_forMonojet.root"  , "qcdCRtree");

  std::cout << "-> Making analyses from inclusive tree..." << std::endl;
  MT2Analysis<MT2EstimateTree>* qcd_old   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "QCD"  , "13TeV_inclusive", mcTree_old, "id>=100 && id<200" ); 
  std::cout << "    QCD done." << std::endl;
  MT2Analysis<MT2EstimateTree>* wjets_old = MT2EstimateTree::makeAnalysisFromInclusiveTree( "WJets", "13TeV_inclusive", mcTree_old, "id>=500 && id<600" ); 
  std::cout << "    WJets done." << std::endl;
  MT2Analysis<MT2EstimateTree>* zjets_old = MT2EstimateTree::makeAnalysisFromInclusiveTree( "ZJets", "13TeV_inclusive", mcTree_old, "id>=600 && id<700" ); 
  std::cout << "    ZJets done." << std::endl;
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_old; 
  
  bgYields_old.push_back( zjets_old );
  bgYields_old.push_back( wjets_old );
  bgYields_old.push_back( qcd_old );


  std::vector< MT2Analysis<MT2EstimateTree>* > mc_rev;
  mc_rev.push_back(qcd); 
  mc_rev.push_back(zjets);
  mc_rev.push_back(wjets);


 


  double bins_jet2_pt[] = {30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330};
  int size_jet2_pt = sizeof(bins_jet2_pt)/sizeof(double)-1;
  double bins_fraction[] = {0,0.0333333,0.0666667,0.1,0.133333,0.166667,0.2,0.233333,0.266667,0.3,0.333333,0.366667,0.4,0.433333,0.466667,0.5,0.533333,0.566667,0.6,0.633333,0.666667,0.7,0.733333,0.766667,0.8,0.833333,0.866667,0.9,0.933333,0.966667,1.};
  int size_fraction = sizeof(bins_fraction)/sizeof(double)-1;


  double bins_jet2_phi[] = {-3.2,-2.98667,-2.77333,-2.56,-2.34667,-2.13333,-1.92,-1.70667,-1.49333,-1.28,-1.06667,-0.853333,-0.64,-0.426666,-0.213333,2.63055e-07,0.213334,0.426667,0.64,0.853334,1.06667,1.28,1.49333,1.70667,1.92,2.13333,2.34667,2.56,2.77333,2.98667,3.2};
  int size_jet2_phi = sizeof(bins_jet2_phi)/sizeof(double)-1;


  double bins_jet2_eta[] = {-2.5,-2.33333,-2.16667,-2,-1.83333,-1.66667,-1.5,-1.33333,-1.16667,-1,-0.833333,-0.666667,-0.5,-0.333333,-0.166667,-6.45717e-08,0.166667,0.333333,0.5,0.666667,0.833333,1,1.16667,1.33333,1.5,1.66667,1.83333,2,2.16667,2.33333,2.5 };
  int size_jet2_eta = sizeof(bins_jet2_eta)/sizeof(double)-1;


  double bins_energy[] = {0,3.3333,6.6666,9.9999,13.3332,16.6665,19.9998,23.3331,26.6664,29.9997,33.333,36.6663,39.9996,43.3329,46.6662,49.9995,53.3328,56.6661,59.9994,63.3327,66.666,69.9993,73.3326,76.6659,79.9992,83.3325,86.6658,89.9991,93.3324,96.6657, 100.};
  int size_energy = sizeof(bins_energy)/sizeof(double)-1;

  std::string selection_2016 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) && jet2_pt>=30.";
  std::string selection_2015 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200.&& jet2_pt>=30.";

  /*
    std::string selection_2016_lowPU = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) && jet2_pt>=30.&& nVert < 12";
    std::string selection_2015_lowPU = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. && jet2_pt>=30. && nVert < 12";
    drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPU_jet2_pt" , "jet2_pt" , selection_2016_lowPU, selection_2015_lowPU, size_jet2_pt, bins_jet2_pt,  "", "","Subleading Jet p_{T}", "GeV" , "p_{T}(jet1) > 200 GeV, nVert<12", "j=2, #geq0b");

    std::string selection_2016_highPU = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter)&& jet2_pt>=30. && nVert >18";
    std::string selection_2015_highPU = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200.&& && jet2_pt>=30.nVert >18";
    drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_highPU_jet2_pt" , "jet2_pt" , selection_2016_highPU, selection_2015_highPU, size_jet2_pt, bins_jet2_pt,  "", "","Subleading Jet p_{T}", "GeV" , "p_{T}(jet1) > 200 GeV, nVert>18", "j=2, #geq0b");
  */

  std::string selection_2016_lowPt2 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter)&& jet2_pt>=30. && jet2_pt<60";
  std::string selection_2015_lowPt2 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. && jet2_pt>=30 && jet2_pt<60";
  

  //  std::string selection_phChHEFfilter_2016 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) &&( abs(met_phi-jet2_phi)<0.3 || abs(met_phi-jet2_phi)>5.9) && (jet2_phEF+jet2_neHEF)<0.8 && (jet2_eEF+jet2_muEF)<0.7 && (jet2_chHEF + jet2_neHEF)>0.5 && jet2_chHEF>0.1";
  //  std::string selection_phChHEFfilter_2015 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( abs(met_phi-jet2_phi)<0.3 || abs(met_phi-jet2_phi)>5.9) && (jet2_phEF+jet2_neHEF)<0.8 && (jet2_eEF+jet2_muEF)<0.7 && (jet2_chHEF + jet2_neHEF)>0.5 && jet2_chHEF>0.1";

  std::string selection_phChHEFfilter_2016 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) &&( abs(met_phi-jet2_phi)<0.3 || abs(met_phi-jet2_phi)>5.9)";

  std::string selection_phChHEFfilter_2015 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( abs(met_phi-jet2_phi)<0.3 || abs(met_phi-jet2_phi)>5.9) ";


drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_EFMetFilter_jet2_pt" , "jet2_pt" , selection_phChHEFfilter_2016, selection_phChHEFfilter_2015, size_jet2_pt, bins_jet2_pt,  "", "","Subleading Jet p_{T}", "GeV" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");



 double bins_jet2_deltaPhi[] = {-6.3,-6.1,-5.9,-5.7,-5.5,-5.3,-5.1,-4.9,-4.7,-4.5,-4.3,-4.1,-3.9,-3.7,-3.5,-3.3,-3.1,-2.9,-2.7,-2.5,-2.3,-2.1,-1.9,-1.7,-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1, 6.3};
  int size_jet2_deltaPhi = sizeof(bins_jet2_deltaPhi)/sizeof(double)-1;

  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_Metjet2_phi" , "met_phi-jet2_phi" , selection_2016, selection_2015, size_jet2_deltaPhi, bins_jet2_deltaPhi,  "", "","#Delta#phi(ME_{T} - jet2)", "" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");

  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_Metjet2_phi" , "met_phi-jet2_phi" , selection_2016_lowPt2, selection_2015_lowPt2, size_jet2_deltaPhi, bins_jet2_deltaPhi,  "", "","#Delta#phi(ME_{T} - jet2)", "" , "p_{T}(jet1)>200 GeV, p_{T}(jet2)<60 GeV", "j=2, #geq0b");


   drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_chPlusNeHEF" , "jet2_neHEF+jet2_chHEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Neutral+Charged Hadron Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
 drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_phPlusEEF" , "jet2_phEF+jet2_eEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Photon + ElectronEnergy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");

  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_nePlusPhHEF" , "jet2_neHEF+jet2_phEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Neutral Hadron+Photon Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");

  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_ePlusMuHEF" , "jet2_eEF+jet2_muEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Electron+Muon Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
 
  // drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_eEF" , "jet2_eEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Electron Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  // drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_muEF" , "jet2_muEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Muon Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");



  /*
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_jet2_pt" , "jet2_pt" , selection_2016, selection_2015, size_jet2_pt, bins_jet2_pt,  "", "","Subleading Jet p_{T}", "GeV" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");
  

   drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_jet2_phi" , "jet2_phi" , selection_2016, selection_2015, size_jet2_phi, bins_jet2_phi,  "", "","#phi (jet2)", "GeV" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_jet2_eta" , "jet2_eta" , selection_2016, selection_2015, size_jet2_eta, bins_jet2_eta,  "", "","#eta (jet2)", "GeV" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");
 
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_jet2_neHEF" , "jet2_neHEF" , selection_2016, selection_2015, size_fraction, bins_fraction,  "", "","Neutral Hadron Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_jet2_chHEF" , "jet2_chHEF" , selection_2016, selection_2015, size_fraction, bins_fraction,  "", "","Charged Hadron Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_jet2_phEF" , "jet2_phEF" , selection_2016, selection_2015, size_fraction, bins_fraction,  "", "","Photon Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_jet2_eEF" , "jet2_eEF" , selection_2016, selection_2015, size_fraction, bins_fraction,  "", "","Electron Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_jet2_muEF" , "jet2_muEF" , selection_2016, selection_2015, size_fraction, bins_fraction,  "", "","Muon Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");
  

 
 drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_phi" , "jet2_phi" , selection_2016_lowPt2, selection_2015_lowPt2, size_jet2_phi, bins_jet2_phi,  "", "","#phi (jet2)", "GeV" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_eta" , "jet2_eta" , selection_2016_lowPt2, selection_2015_lowPt2, size_jet2_eta, bins_jet2_eta,  "", "","#eta (jet2)", "GeV" , "p_{T}(jet1) > 200 GeV", "j=2, #geq0b");
   drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_neHEF" , "jet2_neHEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Neutral Hadron Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_chHEF" , "jet2_chHEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Charged Hadron Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_phEF" , "jet2_phEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Photon Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_eEF" , "jet2_eEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Electron Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_muEF" , "jet2_muEF" , selection_2016_lowPt2, selection_2015_lowPt2, size_fraction, bins_fraction,  "", "","Muon Energy Fraction (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  


  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_neHEnergy" , "jet2_neHEF*jet2_energy" , selection_2016_lowPt2, selection_2015_lowPt2, size_energy, bins_energy,  "", "","Neutral Hadron Energy (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_chHEnergy" , "jet2_chHEF*jet2_energy" , selection_2016_lowPt2, selection_2015_lowPt2, size_energy, bins_energy,  "", "","Charged Hadron Energy (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_phEnergy" , "jet2_phEF*jet2_energy" , selection_2016_lowPt2, selection_2015_lowPt2, size_energy, bins_energy,  "", "","Photon Energy (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_eEnergy" , "jet2_eEF*jet2_energy" , selection_2016_lowPt2, selection_2015_lowPt2, size_energy, bins_energy,  "", "","Electron Energy (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  drawYields( cfg, data, data_old,  mc_rev, bgYields_old, "comparison_lowPt2_jet2_muEnergy" , "jet2_muEF*jet2_energy" , selection_2016_lowPt2, selection_2015_lowPt2, size_energy, bins_energy,  "", "","Muon Energy (jet2)", "" , "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "j=2, #geq0b");
  

  */

  std::vector<TCanvas*> canvases;


  //std::string selection = "nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200.";


  std::string selection = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200.";




  //std::string selection = "(id<100 || id>152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200.";
   std::string selection_filter = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) && jet2_pt>=30.";

   std::string selection_filter_runNum = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) && jet2_pt>=30. && run<=275125.";


   std::string selection_filter_phi = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) && jet2_pt>=30. &&( abs(met_phi-jet2_phi)<1.0 || abs(met_phi-jet2_phi)>5.0) ";

  std::string selection_filter_eta2 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. && jet2_eta>-2.0 && jet2_eta<2.0 &&( id>10 || badFilter)";

  std::string selection_id = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( jet2_id || (id>=600 && id <= 700) )";

  std::string selection_filterAndId = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( jet2_id || (id>=600 && id <= 700) ) &&  (id>10 ||(id>=600 && id <= 700) || badFilter)";

  std::string selection_filterAndIdLowPt = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. && jet2_pt<60. &&( jet2_id || (id>=600 && id <= 700) ) && (id>10 ||(id>=600 && id <= 700) || badFilter)";

 std::string selection_filterAndLowPt = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. && jet2_pt<60. && (id>10 ||(id>=600 && id <= 700) || badFilter)";


 
  //dt.drawRegionYields_fromTree( "metOcaloMet" , "met/caloMet" , selection, 40, 0.0, 5.0, "ME_{T}/caloME_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

  std::string selection_lowpt2 = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. && jet2_pt<60";
 // dt.drawRegionYields_fromTree( "metOcaloMet_lowJet2pt" , "met/caloMet" , selection_lowpt2, 40, 0.5, 10.5, "ME_{T}/caloME_{T}", "GeV", "p_{T}(jet1) > 200 GeV, p_{T}(jet2) < 60 GeV", "N(j) = 2" );


  dt.drawRegionYields_fromTree( "filter_jet1_chHEF" , "jet1_chHEF" , selection_filter, 30, 0., 1.   , "charged Hadron Energy Fraction (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet1_neHEF" , "jet1_neHEF" , selection_filter, 30, 0., 1.   , "neutral Hadron Energy Fraction (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet1_phEF" , "jet1_phEF" , selection_filter, 30, 0., 1.   , "Photon Energy Fraction (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet1_muEF" , "jet1_muEF" , selection_filter, 30, 0., 1.   , "Muon Energy Fraction (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet1_eEF" , "jet1_eEF" , selection_filter, 30, 0., 1.   , "Electron Energy Fraction (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

  dt.drawRegionYields_fromTree( "filter_jet2_chHEF" , "jet2_chHEF" , selection_filter, 30, 0., 1.   , "charged Hadron Energy Fraction (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet2_neHEF" , "jet2_neHEF" , selection_filter, 30, 0., 1.   , "neutral Hadron Energy Fraction (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet2_phEF" , "jet2_phEF" , selection_filter, 30, 0., 1.   , "Photon Energy Fraction (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet2_muEF" , "jet2_muEF" , selection_filter, 30, 0., 1.   , "Muon Energy Fraction (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "filter_jet2_eEF" , "jet2_eEF" , selection_filter, 30, 0., 1.   , "Electron Energy Fraction (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );


  dt.drawRegionYields_fromTree( "filter_jet1_chHEnergy" , "jet1_chHEF*jet1_energy" , selection_filter, 30, 0., 600.   , "charged Hadron Energy (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet1_neHEnergy" , "jet1_neHEF*jet1_energy" , selection_filter, 30, 0., 600.   , "neutral Hadron Energy (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet1_phEnergy" , "jet1_phEF*jet1_energy" , selection_filter, 30, 0., 600.   , "Photon Energy (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet1_muEnergy" , "jet1_muEF*jet1_energy" , selection_filter, 30, 0., 600.   , "Muon Energy (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet1_eEnergy" , "jet1_eEF*jet1_energy" , selection_filter, 30, 0., 600.   , "Electron Energy (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

  dt.drawRegionYields_fromTree( "filter_jet2_chHEnergy" , "jet2_chHEF*jet2_energy" , selection_filter, 30, 0., 100.   , "charged Hadron Energy (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet2_neHEnergy" , "jet2_neHEF*jet2_energy" , selection_filter, 30, 0., 100.   , "neutral Hadron Energy (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet2_phEnergy" , "jet2_phEF*jet2_energy" , selection_filter, 30, 0., 100.   , "Photon Energy (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filter_jet2_muEnergy" , "jet2_muEF*jet2_energy" , selection_filter, 30, 0., 100.   , "Muon Energy (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
dt.drawRegionYields_fromTree( "filter_jet2_eEnergy" , "jet2_eEF*jet2_energy" , selection_filter, 30, 0., 100.   , "Electron Energy (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
/*
 std::string selection_filter_lowChHE = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) && (jet2_chHEF*jet2_energy < 20)";

 dt.drawRegionYields_fromTree( "lowChHE_filter_jet2_neHEF" , "jet2_neHEF" , selection_filter_lowChHE, 30, 0., 1.   , "neutral Hadron Energy Fraction (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "lowChHE_filter_jet2_phEF" , "jet2_phEF" , selection_filter_lowChHE, 30, 0., 1.   , "Photon Energy Fraction (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "lowChHE_filter_jet2_muEF" , "jet2_muEF" , selection_filter_lowChHE, 30, 0., 1.   , "Muon Energy Fraction (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "lowChHE_filter_jet2_eEF" , "jet2_eEF" , selection_filter_lowChHE, 30, 0., 1.   , "Electron Energy Fraction (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );


 dt.drawRegionYields_fromTree( "lowChHE_filter_jet2_neHEnergy" , "jet2_neHEF*jet2_energy" , selection_filter_lowChHE, 30, 0., 100.   , "neutral Hadron Energy (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "lowChHE_filter_jet2_phEnergy" , "jet2_phEF*jet2_energy" , selection_filter_lowChHE, 30, 0., 100.   , "Photon Energy (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "lowChHE_filter_jet2_muEnergy" , "jet2_muEF*jet2_energy" , selection_filter_lowChHE, 30, 0., 100.   , "Muon Energy (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "lowChHE_filter_jet2_eEnergy" , "jet2_eEF*jet2_energy" , selection_filter_lowChHE, 30, 0., 100.   , "Electron Energy (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );


*/







 dt.drawRegionYields_fromTree( "jet1_eta" , "jet1_eta" , selection_filter, 30, -2.5  , 2.5   , "#eta (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "jet2_eta" , "jet2_eta" , selection_filter, 30, -2.5  , 2.5   , "#eta (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );


 dt.drawRegionYields_fromTree( "jet1_phi" , "jet1_phi" , selection_filter, 32, -3.2  , 3.2   , "#phi (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "jet2_phi" , "jet2_phi" , selection_filter, 32, -3.2 , 3.2   , "#phi (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );


 dt.drawRegionYields_fromTree( "jet1_pt" , "jet1_pt" , selection_filter, 40, 200., 800., "Leading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

 //dt.drawRegionYields_fromTree( "lowChHE_filter_metOcaloMet" , "met/caloMet" , selection_filter_lowChHE, 30, 0., 2.   , "ME_{T} / calo-ME{T}" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

 /*
   dt.drawRegionYields_fromTree( "jet1_eta" , "jet1_eta" , selection, 30, -2.5  , 2.5   , "#eta (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
   dt.drawRegionYields_fromTree( "jet2_eta" , "jet2_eta" , selection, 30, -2.5  , 2.5   , "#eta (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );


 dt.drawRegionYields_fromTree( "jet1_phi" , "jet1_phi" , selection, 32, -3.2  , 3.2   , "#phi (jet1)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "jet2_phi" , "jet2_phi" , selection, 32, -3.2 , 3.2   , "#phi (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

  dt.drawRegionYields_fromTree( "filterAndId_jet2_phi" , "jet2_phi" , selection_filterAndId, 32, -3.2 , 3.2   , "#phi (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
dt.drawRegionYields_fromTree( "filterAndId_lowPt_jet2_phi" , "jet2_phi" , selection_filterAndIdLowPt, 32, -3.2 , 3.2   , "#phi (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

 dt.drawRegionYields_fromTree( "filterAndId_lowPt_jet2_eta" , "jet2_eta" , selection_filterAndIdLowPt, 30, -2.5  , 2.5   , "#eta (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

dt.drawRegionYields_fromTree( "filterAndlowPt_jet2_phi" , "jet2_phi" , selection_filterAndLowPt, 32, -3.2 , 3.2   , "#phi (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

 dt.drawRegionYields_fromTree( "filterAndlowPt_jet2_eta" , "jet2_eta" , selection_filterAndLowPt, 30, -2.5  , 2.5   , "#eta (jet2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
*/

//  dt.drawRegionYields_fromTree( "delta_phi_jets" , "jet1_phi-jet2_phi" , selection, 60, -1. , 1.0   , "#Delta#phi (jet1-2)" , ""    , "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
/*
 

  dt.drawRegionYields_fromTree( "badFilter_jet2_pt" , "jet2_pt" , selection_filter, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

  dt.drawRegionYields_fromTree( "filter_eta2_jet2_pt" , "jet2_pt" , selection_filter_eta2, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

 dt.drawRegionYields_fromTree( "id_jet2_pt" , "jet2_pt" , selection_id, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
  dt.drawRegionYields_fromTree( "filterAndId_jet2_pt" , "jet2_pt" , selection_filterAndId, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );
*/
 //  return jet_id[j]>=3 && jet_chHEF[j]>0.05 && jet_neHEF[j]<0.8 && jet_phEF[j]<0.7;


 //meh, doesn't work that well... 
 std::string selection_deltaPhiMetJet2= "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) &&( abs(met_phi-jet2_phi)<0.3 || abs(met_phi-jet2_phi)>5.5) ";
 std::string selection_phChHEFfilter = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) &&( abs(met_phi-jet2_phi)<0.3 || abs(met_phi-jet2_phi)>5.9) && (jet2_phEF+jet2_neHEF)<0.8 && (jet2_eEF+jet2_muEF)<0.7 && (jet2_chHEF + jet2_neHEF)>0.5 && jet2_chHEF>0.1";
 //std::string selection_phChHEFfilter = "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) && ( (jet2_phEF + jet2_eEF + jet2_muEF) < 0.8) && ((jet2_chHEF + jet2_neHEF) >0.2 && jet2_neHEF>0.01 && jet2_muEF<0.5 )";
 dt.drawRegionYields_fromTree( "phChHEF_jet2_pt" , "jet2_pt" , selection_phChHEFfilter, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV, chHEF>0.2", "N(j) = 2" );
 dt.drawRegionYields_fromTree( "deltaPhiMetJet2_jet2_pt" , "jet2_pt" , selection_deltaPhiMetJet2, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV, chHEF>0.2", "N(j) = 2" );

 



 std::string selection_deltaPhi_mildFilter= "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) &&( abs(met_phi-jet2_phi)<0.3 || abs(met_phi-jet2_phi)>5.5) && jet2_muEF<0.9 && jet2_phEF<0.9 && jet2_eEF<0.9 ";

 dt.drawRegionYields_fromTree( "mildFilter_jet2_pt" , "jet2_pt" , selection_deltaPhi_mildFilter, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );



 std::string selection_deltaPhi_milderFilter= "(id<100 || id>=152) && nJets==2 && deltaPhiMin<0.3 && jet1_pt>200. && met>200. &&( id>10 || badFilter) &&( abs(met_phi-jet2_phi)<0.3 || abs(met_phi-jet2_phi)>5.5) && jet2_muEF<0.95 && jet2_phEF<0.95 && jet2_eEF<0.95 ";

 dt.drawRegionYields_fromTree( "milderFilter_jet2_pt" , "jet2_pt" , selection_deltaPhi_milderFilter, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );




 dt.drawRegionYields_fromTree( "runNum_jet2_pt" , "jet2_pt" , selection_filter_runNum , 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );



 canvases = dt.drawRegionYields_fromTree( "jet2_pt" , "jet2_pt" , selection_filter, 20, 30., 330., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", "N(j) = 2" );

  float mcSF = MT2DrawTools::getDataMCSF( canvases[0] );
  //dt.set_mcSF( mcSF );

  int nControlRegion;
  float boh, boh_err;
  fitQCDMonojet( canvases[0], "boh" , nControlRegion, boh, boh_err, cfg );



  MT2Analysis<MT2Estimate>* qcdMonojet = new MT2Analysis<MT2Estimate>( "monojet_qcdEstimate", cfg.regionsSet() );
  MT2Analysis<MT2Estimate>* nCRMonojet = new MT2Analysis<MT2Estimate>( "monojet_nCR"        , cfg.regionsSet() );
  MT2Analysis<MT2Estimate>* rMonojet   = new MT2Analysis<MT2Estimate>( "monojet_r"          , cfg.regionsSet() );



  std::set<MT2Region> regions = qcdMonojet->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nJetsMin()!=1 || iR->nJetsMax()!=1 ) continue; // only monojet

    int nBJets = iR->nBJetsMin();

    float ptMin = iR->htMin();
    float ptMax = iR->htMax();

    MT2Estimate* thisEst = qcdMonojet->get( *iR );
    MT2Estimate* thisNCR = nCRMonojet->get( *iR );
    MT2Estimate* thisR   = rMonojet  ->get( *iR );

    std::string fullSelection(Form("%s && %s && jet1_pt>%f && jet1_pt<%f", selection.c_str(), iR->sigRegion()->getBJetCuts().c_str(), ptMin, ptMax ) );

    std::string bJetsLabel = (nBJets==0) ? "b = 0" : "b #geq 1";
    canvases = dt.drawRegionYields_fromTree( Form("jet2_pt_%s", iR->getName().c_str()), "jet2_pt", fullSelection, 20, 0., 300., "Subleading Jet p_{T}", "GeV", "p_{T}(jet1) > 200 GeV", bJetsLabel );

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
    

  std::string outfileName( Form( "%s/qcdEstimateMonojet.root", cfg.getEventYieldDir().c_str() ) );
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











void fitQCDMonojet( TCanvas* c1, const std::string& regionName, int &nCR, float &qcdFraction, float &qcdFractionError, MT2Config cfg ) {

  TList* list = MT2DrawTools::getCorrectList( c1 );
  TGraphAsymmErrors* gr_data  = (TGraphAsymmErrors*)list->FindObject("Graph");

  nCR = MT2DrawTools::graphIntegral( gr_data, 30., 60. );

  THStack* bgStack = (THStack*)list->FindObject( "bgStack" );

  TList* mcHists = bgStack->GetHists();

  TH1D* h1_qcd   = (TH1D*)mcHists->FindObject( "h1_QCD_HT200toInf_j1toInf_b0toInf" );
  TH1D* h1_wjets = (TH1D*)mcHists->FindObject( "h1_WJets_HT200toInf_j1toInf_b0toInf" );
  TH1D* h1_zjets = (TH1D*)mcHists->FindObject( "h1_ZJets_HT200toInf_j1toInf_b0toInf" );


  int nPoints = gr_data->GetN();


  TGraphAsymmErrors* gr_data_sub = new TGraphAsymmErrors(nPoints);

  for(int i=0; i < nPoints; i++){

    Double_t x, y;
    gr_data->GetPoint(i, x, y);

    int iBin = h1_wjets->FindBin(x);
    float wJets = h1_wjets->GetBinContent(iBin);
    float zJets = h1_zjets->GetBinContent(iBin);
    float ewk = wJets + zJets;

    gr_data_sub->SetPoint( i, x, y - ewk );
  }

  gr_data->SetMarkerStyle(24);
  gr_data_sub->SetMarkerStyle(20);

  h1_qcd->SetFillColor(kQCD);
  h1_qcd->SetLineWidth(2);


  TCanvas* c2 = new TCanvas("c2", "", 600, 600);


  float yMaxScale = 1.1;
  // float yMax = gr_data->GetHistogramm()->GetMaximum()*yMaxScale;
  // if( yMax <= 0 ) std::cout << "You fucked up" << std::endl;

  double * y = gr_data->GetY();
  int locmax= TMath::LocMax(nPoints,y);
  double yMax = y[locmax];
  yMax*=1.25;

  TLegend* legend = new TLegend( 0.59, 0.92-(3)*0.06, 0.89, 0.92 );
  legend->SetTextSize(0.04);
  // legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( h1_qcd ,"QCD", "F" );
  legend->AddEntry( gr_data ,"Data", "PL" );
  legend->AddEntry( gr_data_sub ,"Data - EWK MC", "PL" );



  TH2D* h2_axes = new TH2D("axes", "", 20, 30, 330, 10, 0., yMax );
  h2_axes->SetYTitle("Entries/ 30 GeV");
  h2_axes->SetXTitle("Subleading Jet p_{T} [GeV]");
  h2_axes->Draw("");

  h1_qcd->Draw("histo same");
  gr_data->Draw("P same");
  gr_data_sub->Draw("P same");

  legend->Draw("same");

  gPad->RedrawAxis();

  c2->SaveAs(Form( "%s/test.eps", cfg.getEventYieldDir().c_str() ) );


  /*
    int bin30 = h1_qcd->FindBin(30);
    int bin60 = h1_qcd->FindBin(60);

    double qcd3060_err;
    double qcd3060 = h1_qcd->IntegralAndError( bin30, bin60, qcd3060_err );
    double wjets3060 = h1_wjets->Integral( bin30, bin60 );
    double zjets3060 = h1_zjets->Integral( bin30, bin60 );

    double mcTot = qcd3060 + wjets3060 + zjets3060;

    qcdFraction = qcd3060 / mcTot;
    qcdFractionError = qcd3060_err / mcTot;

  */


}

















void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* data_of , std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields,  std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields_of, const std::string& saveName, const std::string& varName, const std::string& selection,  const std::string& selection_of, int nBins, double* bins, std::string name, std::string name_of, std::string axisName, const std::string& units , const std::string& kinCuts, const std::string& topoCuts, bool drawData ) {

  float xMin = bins[0];
  float xMax = bins[nBins];

  //  float binWidth = (xMax-xMin)/nBins;
  //  if( axisName=="" ) axisName = varName;

  std::vector<int> colors;
  if( bgYields.size()==1 ) { // estimates
    colors.push_back(855); 
    // colors.push_back(430); 
    // colors.push_back(418); 
  } else { // mc
    //colors.push_back(430); // other=zll
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    colors.push_back(419); // z+jets
    // colors.push_back(855); // top
    // colors.push_back(); // other
  }


  bool shapeNorm = true;

  std::string fullPathPlots = cfg.getEventYieldDir() + "/plots2016vs2015";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_data = data->get(thisRegion)->tree;
    TH1D* h1_data = new TH1D("h1_data", "", nBins, bins );
    tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );
    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data, false);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.2);

    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<bgYields.size(); ++i ) { 
      TTree* tree_mc = (bgYields[i]->get(thisRegion)->tree);
      std::string thisName = "h1_" + bgYields[i]->getName();
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, bins );
      h1_mc->Sumw2();
      tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", cfg.lumi(), selection.c_str()) ); 
      // MT2DrawTools::addOverflowSingleHisto(h1_mc);
      histos_mc.push_back(h1_mc);
    }

    TH1D* mc_sum;
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      if( i==0 ) {
        mc_sum = new TH1D( *histos_mc[i] );
        mc_sum->SetName("mc_sum");
      } else {
        mc_sum->Add( histos_mc[i] );
      }
    }
    /*
    std::vector< TH1D* > histos_mc_of;
    for( unsigned i=0; i<bgYields_of.size(); ++i ) { 
      TTree* tree_mc = (bgYields_of[i]->get(thisRegion)->tree);
      std::string thisName = "h1_of_" + bgYields_of[i]->getName();
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, bins );
      h1_mc->Sumw2();
      tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s)", cfg.lumi(), selection_of.c_str()) ); 
      // MT2DrawTools::addOverflowSingleHisto(h1_mc);
      histos_mc_of.push_back(h1_mc);
    }

    TH1D* mc_sum_of;
    for( unsigned i=0; i<histos_mc_of.size(); ++i ) { 
      if( i==0 ) {
        mc_sum_of = new TH1D( *histos_mc_of[i] );
        mc_sum_of->SetName("mc_sum");
      } else {
        mc_sum_of->Add( histos_mc_of[i] );
      }
    }
    */


    std::cout << "Integrals: " << h1_data->Integral(1, nBins) << "\t" << mc_sum->Integral(1, nBins) << std::endl;
    float scaleFactor = h1_data->Integral(1, nBins)/mc_sum->Integral(1, nBins);   
    // if( shapeNorm ) 
      std::cout << "SF: " << scaleFactor << std::endl;

    TObject* oldStack = gROOT->FindObject("bgStack");
    if( oldStack ) delete oldStack;
    TH1D* histo_mc;
    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = bgYields.size() - i - 1;
      histos_mc[index]->SetFillColor( colors[index] );
      histos_mc[index]->SetLineColor( kBlack );
      //histos_mc[index]->SetLineWidth( 0 );
      //if( shapeNorm ) {
      histos_mc[index]->Scale( scaleFactor );
      //}
      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack.Add(histos_mc[index]);
    }

    /*
    TObject* oldStack_of = gROOT->FindObject("bgStack_of");
    if( oldStack_of ) delete oldStack_of;
    TH1D* histo_mc_of;
    THStack bgStack_of("bgStack_of", "");
    for( unsigned i=0; i<histos_mc_of.size(); ++i ) { 
      int index = bgYields_of.size() - i - 1;
      // histos_mc_of[index]->SetFillColor( colors[index] );
      // histos_mc_of[index]->SetFillStyle( 3444 );
      histos_mc_of[index]->SetLineColor( kBlack );
      histos_mc_of[index]->SetLineWidth( 2 );
      histos_mc_of[index]->SetLineStyle( i+1 );
      if( shapeNorm ) {
        histos_mc_of[index]->Scale( scaleFactor );
      }
      if(i==0) histo_mc_of = (TH1D*) histos_mc_of[index]->Clone("histo_mc_of");
      else histo_mc_of->Add(histos_mc_of[index]);
      bgStack_of.Add(histos_mc_of[index]);
    }
    */

    TTree* tree_data_of = data_of->get(thisRegion)->tree;
    TH1D* h1_data_of = new TH1D("h1_data_of", "", nBins, bins );
    tree_data_of->Project( "h1_data_of", varName.c_str(), selection_of.c_str() );
    // if( shapeNorm )   
    h1_data_of->Scale( 5.4/2.3 );
    //    if( shapeNorm )   h1_data_of->Scale( scaleFactor);
    TH1D* gr_data_of = new TH1D("gr_data_of", "", nBins, bins);
    gr_data_of->Add(h1_data_of);
    //TGraphAsymmErrors* gr_data_of = MT2DrawTools::getPoissonGraph(h1_data_of);
    gr_data_of->SetMarkerStyle(24);
    gr_data_of->SetMarkerSize(1);


    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph( h1_data_of, h1_data);
    // TH1D* g_ratio = new TH1D("g_ratio", "", nBins, bins);
    // g_ratio->Add(h1_data_of);
    //   g_ratio->Divide( h1_data );
    g_ratio->SetMarkerStyle( 20 ); g_ratio->SetMarkerSize(1.2);
    /*
    //    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(histo_mc, histo_mc_of);
    TH1D* h_ratio = new TH1D("h_ratio","",nBins, bins);
    h_ratio->Add(h1_data_of);
    h_ratio->Divide(h1_data);
    h_ratio->SetMarkerSize(1.4);
    h_ratio->SetMarkerStyle(20);
    */
    //    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
    // TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data_of, histo_mc);
    
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    float lumiErr=0;
    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);
    
    // TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TF1* fSF = 0;// MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TGraphErrors* SFFitBand = (fSF) ? MT2DrawTools::getSFFitBand(fSF, xMin, xMax) : 0;

    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );


    TCanvas* c1 = new TCanvas(Form("c1_%s", iMT2->getName().c_str()), "", 600, 600);
    c1->cd();
    
    TCanvas* c1_log = new TCanvas(Form("c1_log_%s", iMT2->getName().c_str()), "", 600, 600);

    float yMaxScale = 1.1;
    float yMax1 = histo_mc->GetMaximum()*yMaxScale ;
    //    float yMax1 = h1_data->GetMaximum()*yMaxScale ;
    float yMax2 = yMaxScale*(histo_mc->GetMaximum() + sqrt(histo_mc->GetMaximum()));
    //  float yMax2 = yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
    float yMax3 = yMaxScale*(bgStack.GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    if( histo_mc->GetNbinsX()<2 ) yMax *=3.;
    yMax*=1.25;


    std::string xAxisTitle;
    if( units!="" ) 
      xAxisTitle = (std::string)(Form("%s [%s]", axisName.c_str(), units.c_str()) );
    else
      xAxisTitle = (std::string)(Form("%s", axisName.c_str()) );


    std::string yAxisTitle;
    yAxisTitle = (std::string)("Events");
    /*
    if(binWidth>0.99){
      if( units!="" ) 
	yAxisTitle = (std::string)(Form("Events / (%.0f %s)", binWidth, units.c_str()));
      else
	yAxisTitle = (std::string)(Form("Events / (%.0f)", binWidth));
    }
    else{
      if( units!="" ) 
	yAxisTitle = (std::string)(Form("Events / (%.2f %s)", binWidth, units.c_str()));
      else
	yAxisTitle = (std::string)(Form("Events / (%.2f)", binWidth));
	}*/

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();  

    TPad* pad1 = 0;
    pad1 = MT2DrawTools::getCanvasMainPad();
    pad1->Draw();
    pad1->cd();
   

    h2_axes->Draw();

   
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.1, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();
    
    TPad* pad1_log = 0;
    pad1_log = MT2DrawTools::getCanvasMainPad( true );
    pad1_log->Draw();
    pad1_log->cd();

  
    h2_axes_log->Draw();
  
  
    std::vector<std::string> niceNames = thisRegion.getNiceNames();

    for( unsigned i=0; i<niceNames.size(); ++i ) {
      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.04);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);

      if( i==0 ) {

        if(kinCuts!="") {
          regionText->AddText( kinCuts.c_str() );
        } else {
          regionText->AddText( niceNames[i].c_str() );
        }

      } else if( i==1 ) {
    
        if(topoCuts!="") {
          regionText->AddText( topoCuts.c_str() );
        } else {
          regionText->AddText( niceNames[i].c_str() );
        }

      }
      pad1->cd();
      regionText->Draw("same");
      
      pad1_log->cd();
      regionText->Draw("same");
    }


    /*
    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      if( scaleFactor!=1. ) {
	normText->AddText( Form("#splitline{2015 Data scaled}{by %.2f}", scaleFactor) );
      }
      pad1->cd();
      normText->Draw("same");
      pad1_log->cd();
      normText->Draw("same");
    }
    */

    TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+2)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( gr_data, "Data 2016", "EP" );
    //if(drawData == true) 
    legend->AddEntry( gr_data_of, "Data 2015", "EP" );
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
     legend->AddEntry( histos_mc[i], Form("%s %s", bgYields[i]->getFullName().c_str(), name.c_str()  ) , "F" );
     //  legend->AddEntry( histos_mc_of[i], Form("%s %s", bgYields_of[i]->getFullName().c_str(), name_of.c_str() ) , "L" );
    }
    //legend->AddEntry( mcBand, "MC Uncert.", "F" );


    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    
    TPaveText* fitText = (fSF) ? MT2DrawTools::getFitText( fSF ) : 0;


    float yMinR=0.0;
    float yMaxR=2.0;

    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
 
    h2_axes_ratio->SetYTitle( "2015/2016" );
    // if( drawData == false ) h2_axes_ratio->SetYTitle( "#frac{ee+#mu#mu}{e#mu+#mue}" );

    c1->cd();
    pad1->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    //bgStack.Draw("L same");
    //   bgStack_of.Draw("L same");
    gr_data->Draw("Ep same");
    //if(drawData==true) 
    gr_data_of->Draw("EP same");
    if( drawData == false ) {
      MT2DrawTools::addLabels( (TCanvas*)pad1, cfg.lumi(), "CMS Simulation" );
    } else {
      MT2DrawTools::addLabels( (TCanvas*)pad1, cfg.lumi(), "CMS Preliminary" );
    }


    //  if( !shapeNorm )
    //   fitText->Draw("same");
    // ratioText->Draw("same");
  
    gPad->RedrawAxis();

    c1_log->cd();
    pad1_log->cd();
    legend->Draw("same");
    bgStack.Draw("L same");
    bgStack.Draw("histo same");
    //   bgStack_of.Draw("L same");
    gr_data->Draw("p same");
    //if(drawData==true)
    gr_data_of->Draw("PE same");
    //  labelTop->Draw("same");
    //if( !shapeNorm )
    // fitText->Draw("same");
    //  ratioText->Draw("same");

    //if( drawData == false ) {
    //  MT2DrawTools::addLabels( (TCanvas*)pad1_log, cfg.lumi(), "CMS Simulation" );
    //} else {
    MT2DrawTools::addLabels( (TCanvas*)pad1_log, cfg.lumi(), "CMS Preliminary" );
      //}
    gPad->RedrawAxis();

    /*
      TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
      line->SetLineColor(1);
      TLine* lineSF = MT2DrawTools::getSFLine(integral_data, integral_mc, xMin, xMax);
      TGraphErrors* SFband = MT2DrawTools::getSFBand(integral_data, error_data, integral_mc, error_mc, xMin, xMax);
    */

    c1->cd();
    TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
    pad2->Draw();
    pad2->cd();

    h2_axes_ratio->Draw("");
 
    /*  line->Draw("same");
	SFband->Draw("3,same");
	lineSF->Draw("same");
    */
    lineCentral->Draw("same");
    if( !shapeNorm ){


      systBand->Draw("3,same");
      lineCentral->Draw("same");

      //         SFFitBand->Draw("3,same");
      //   fSF->Draw("same");

    }


    // h_ratio->Draw("PE,same"); 
    //if(drawData==true)
      g_ratio->Draw("PE,same");    
    gPad->RedrawAxis();



    c1_log->cd();
    TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
    pad2_log->Draw();
    pad2_log->cd();

    h2_axes_ratio->Draw(""); 
    
    lineCentral->Draw("same");
    if( !shapeNorm ){

      systBand->Draw("3,same");
      lineCentral->Draw("same");

      //  SFFitBand->Draw("3,same");
      //  fSF->Draw("same");
    }
    /*
    line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same"); */

    // h_ratio->Draw("PE,same");
    // if(drawData==true)
    g_ratio->Draw("PE,same");
   
    gPad->RedrawAxis();


    c1->SaveAs( Form("%s/%s.eps", fullPathPlots.c_str(), saveName.c_str() ) );
    c1->SaveAs( Form("%s/%s.png", fullPathPlots.c_str(), saveName.c_str() ) );
    c1->SaveAs( Form("%s/%s.pdf", fullPathPlots.c_str(), saveName.c_str() ) );

    c1_log->SaveAs( Form("%s/%s_log.eps", fullPathPlots.c_str(), saveName.c_str() ) );
    c1_log->SaveAs( Form("%s/%s_log.png", fullPathPlots.c_str(), saveName.c_str() ) );
    c1_log->SaveAs( Form("%s/%s_log.pdf", fullPathPlots.c_str(), saveName.c_str() ) );

    delete c1;
    delete h2_axes;

    delete c1_log;
    delete h2_axes_log;

    delete h2_axes_ratio;

    delete h1_data; delete h1_data_of;// delete h_ratio;
    delete g_ratio; delete gr_data_of;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

    //  for( unsigned i=0; i<histos_mc_of.size(); ++i )
    //     delete histos_mc_of[i];

  }// for MT2 regions

}













