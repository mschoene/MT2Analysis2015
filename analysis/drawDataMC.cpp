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



std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units="" );



int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./drawDataMC [configFileName] [lumi/shape]" << std::endl;
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


  std::string mcFile = cfg.getEventYieldDir() + "/analyses.root";
  std::string dataFile = cfg.getEventYieldDir() + "/analyses.root";

  MT2Analysis<MT2EstimateTree>* zjets = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "ZJets");
  zjets->setFullName("Z+Jets");
  zjets->setColor(kZJets);

  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "WJets");
  wjets->setFullName("W+Jets");
  wjets->setColor(kWJets);

  MT2Analysis<MT2EstimateTree>* top   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "Top");
  top->setFullName("Top");
  top->setColor(kTop);

  MT2Analysis<MT2EstimateTree>* qcd   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "QCD");
  qcd->setFullName("QCD");
  qcd->setColor(kQCD);

  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(dataFile, "data");
  data->setFullName("Data");

  std::vector< MT2Analysis<MT2EstimateTree>* > mc;
  mc.push_back(qcd);
  mc.push_back(wjets);
  mc.push_back(zjets);
  mc.push_back(top);

  std::string plotsDir = cfg.getEventYieldDir() + "/plotsDataMC";
  if( shapeNorm ) plotsDir += "_shape";


  MT2DrawTools dt(plotsDir, cfg.lumi() );
  dt.set_shapeNorm( shapeNorm );


  // +++++++++++++++++++++++++
  // +++      high-HT      +++
  // +++++++++++++++++++++++++

  dt.set_lumi( cfg.lumi_JetHT() );

  float htMin=1000, htMax=-1;
  std::string cutsLabel = getCutLabel(htMin, htMax, "H_T", "GeV");

  std::string selection = "weight*(ht>1000. && nJets>1 && met>30. && mt2>10. && deltaPhiMin>0.3 && diffMetMht<0.5*met)/puWeight";
  dt.drawRegionYields_fromTree( data, mc, "nVert_noPU" , "nVert" , selection, 50, 0.5, 50.5, "Number of Vertices", "", cutsLabel );


  //selection = "weight*(ht>1000. && nJets>1 && met>30. && mt2>10. && deltaPhiMin>0.3 && diffMetMht<0.5*met && nJetHF30==0)";
  selection = "weight*(ht>=1000. && nJets>1 && met>30. && mt2>10. && deltaPhiMin>0.3 && diffMetMht<0.5*met)";
  
  //drawYields( cfg, data, mc, "nVert" , "nVert" , selection, 50, 0.5, 50.5, "Number of Vertices", "", htMin, htMax );
  dt.drawRegionYields_fromTree( data, mc, "mt2"   , "mt2"   , selection, 60, 10., 310., "M_{T2}", "GeV", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "met"   , "met"   , selection, 80, 30., 430., "Missing E_{T}", "GeV", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "mht"   , "mht"   , selection, 80, 30., 430., "Missing H_{T}", "GeV", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "ht"    , "ht"    , selection, 50, 1000., 3500., "H_{T}", "GeV", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "nJets" , "nJets" , selection, 12, 1.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "nJetHF" , "nJetHF" , selection, 7, -0.5, 6.5, "N(jets, p_{T} > 30 GeV & |#eta|>3.0)", "", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel );

  dt.drawRegionYields_fromTree( data, mc, "mt2_tail"   , "mt2"   , selection, 100, 0., 5000., "M_{T2}", "GeV", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "met_tail"   , "met"   , selection, 100,  0., 5000., "Missing E_{T}", "GeV", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "ht_tail"    , "ht"    , selection, 120, 1000., 13000., "H_{T}", "GeV", cutsLabel );

  selection = "weight*(ht>=1000. && nJets>1 && met>30. && mt2>10. && diffMetMht<0.5*met)";

  dt.drawRegionYields_fromTree( data, mc, "deltaPhiMin", "deltaPhiMin", selection, 32, 0.0, 3.2, "min #Delta#phi(jets, ME_{T})", "", cutsLabel );

  selection = "weight*(ht>=1000. && nJets>1 && met>30. && mt2>10. && deltaPhiMin>0.3)";

  dt.drawRegionYields_fromTree( data, mc, "diffMetMht_overMet", "diffMetMht/met", selection, 30, 0., 3.0, "|ME_{T}-MH_{T}|/ME_{T}", "", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "diffMetMht", "diffMetMht", selection, 100, 0., 500., "|ME_{T}-MH_{T}|", "GeV", cutsLabel );





  // +++++++++++++++++++++++++
  // +++      low-HT       +++
  // +++++++++++++++++++++++++

  dt.set_lumi( cfg.lumi_HTMHT() );

  htMin=450; 
  htMax=1000;

  cutsLabel = std::string(Form("%.0f < H_{T} < %.0f GeV", htMin, htMax) );

  //selection = "weight*(ht>450. && ht<1000. && met>200. && nJets>1 && mt2>10. && deltaPhiMin>0.3 && diffMetMht<0.5*met && nJetHF30==0)";
  selection = "weight*(ht>450. && ht<1000. && met>200. && nJets>1 && mt2>10. && deltaPhiMin>0.3 && diffMetMht<0.5*met)";

  

  dt.drawRegionYields_fromTree( data, mc, "lowHT_nVert" , "nVert" , selection, 25, 0.5, 50.5, "Number of Vertices", "", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "lowHT_mt2"   , "mt2"   , selection, 18, 10., 910., "M_{T2}", "GeV", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "lowHT_met"   , "met"   , selection, 20, 200., 1200., "Missing E_{T}", "GeV", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "lowHT_ht"    , "ht"    , selection, 11, 450., 1000., "H_{T}", "GeV", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "lowHT_nJets" , "nJets" , selection, 12, 1.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "lowHT_nJetHF", "nJetHF", selection, 7, -0.5, 6.5, "N(jets, p_{T} > 30 GeV & |#eta|>3.0)", "", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "lowHT_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel );

  selection = "weight*(ht>450. && ht<1000. && met>200. && nJets>1 && mt2>10. &&  diffMetMht<0.5*met)";

  dt.drawRegionYields_fromTree( data, mc, "lowHT_deltaPhiMin", "deltaPhiMin", selection, 32, 0.0, 3.2, "min #Delta#phi(jets, ME_{T})", "", cutsLabel );

  selection = "weight*(ht>450. && ht<1000. && met>200. && nJets>1 && mt2>10. && deltaPhiMin>0.3)";

  dt.drawRegionYields_fromTree( data, mc, "lowHT_diffMetMht_overMet", "diffMetMht/met", selection, 30, 0.0, 3.0, "|ME_{T}-MH_{T}|/ME_{T}", "", cutsLabel );
  dt.drawRegionYields_fromTree( data, mc, "lowHT_diffMetMht", "diffMetMht", selection, 24, 0., 1200., "|ME_{T}-MH_{T}|", "GeV", cutsLabel );
  
  return 0;

}



std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units ) {

  std::string cutLabel;
  if( theMax>theMin ) cutLabel = std::string(Form("%.0f < %s < %.0f %s", theMin, name.c_str(), theMax, units.c_str()) );
  else                cutLabel = std::string(Form("%s > %.0f %s", name.c_str(), theMin, units.c_str()) );

  return cutLabel;

}
