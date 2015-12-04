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
    std::cout << "USAGE: ./drawGammaControlRegionDataMC [configFileName] [lumi/shape]" << std::endl;
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


  std::string mcFile = cfg.getEventYieldDir() + "/gammaControlRegion/mc.root";
  std::string dataFile = cfg.getEventYieldDir() + "/gammaControlRegion/data.root";

  MT2Analysis<MT2EstimateTree>* mc_   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "gammaCRtree_loose");
  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(dataFile, "gammaCRtree_loose");

  MT2Analysis<MT2EstimateTree>* mc_prompt = MT2EstimateTree::makeAnalysisFromInclusiveTree( "Prompt" , cfg.crRegionsSet(), mc_, "prompt==2" ); 
  MT2Analysis<MT2EstimateTree>* mc_nip    = MT2EstimateTree::makeAnalysisFromInclusiveTree( "Fragm." , cfg.crRegionsSet(), mc_, "prompt==1" ); 
  MT2Analysis<MT2EstimateTree>* mc_fake   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "Fakes"  , cfg.crRegionsSet(), mc_, "prompt==0" ); 
  mc_prompt->setColor(18);
  mc_nip   ->setColor(38);
  mc_fake  ->setColor(46);

  std::vector< MT2Analysis<MT2EstimateTree>* > mc;
  mc.push_back(mc_prompt);
  mc.push_back(mc_nip);
  mc.push_back(mc_fake);

  std::string plotsDir = cfg.getEventYieldDir() + "/gammaControlRegion/plotsDataMC";
  if( shapeNorm ) plotsDir += "_shape";


  MT2DrawTools dt(plotsDir, cfg.lumi() );
  dt.set_shapeNorm( shapeNorm );

  dt.set_data( data );
  dt.set_mc( &mc );

  // +++++++++++++++++++++++++
  // +++     multijet      +++
  // +++++++++++++++++++++++++

  dt.set_lumi( cfg.lumi_SinglePhoton() );

  float htMin=200, htMax=-1;
  std::string cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  std::string selection = "ptGamma>180. && ht>450. && nJets>1 && iso<2.5 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200. && met>200";

  dt.drawRegionYields_fromTree( "nVert"      , "nVert"              , selection, 25, 0.5   , 50.5  , "Number of Vertices"               , ""    , cutsLabel, "#geq 2 j" );
  dt.drawRegionYields_fromTree( "mt2"        , "mt2"                , selection, 40, 0.    , 1000. , "M_{T2} (Photon Removed)"          , "GeV" , cutsLabel, "#geq 2 j" );
  dt.drawRegionYields_fromTree( "met"        , "met"                , selection, 36, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, "#geq 2 j" );
  dt.drawRegionYields_fromTree( "ht"         , "ht"                 , selection, 64, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, "#geq 2 j" );
  dt.drawRegionYields_fromTree( "nJets"      , "nJets"              , selection, 10, 1.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, "#geq 2 j" );
  dt.drawRegionYields_fromTree( "nBJets"     , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, "#geq 2 j" );
  dt.drawRegionYields_fromTree( "ptGamma"    , "ptGamma"            , selection, 36, 180.  , 1080  , "Photon p_{T}"                     , "GeV" , cutsLabel, "#geq 2 j" );
  dt.drawRegionYields_fromTree( "etaGamma"   , "etaGamma"           , selection, 50, -2.5  , 2.5   , "Photon #eta"                      , ""    , cutsLabel, "#geq 2 j" );
  dt.drawRegionYields_fromTree( "sietaieta"  , "sietaieta"          , selection, 12, 0.0075, 0.0111, "Photon #sigma_{i#eta i#eta}"      , ""    , cutsLabel, "#geq 2 j" );
  dt.drawRegionYields_fromTree( "isoRC"      , "gamma_chHadIsoRC[0]", selection, 40, 0.    , 10.   , "Random Cone Isolation"            , "GeV" , cutsLabel, "#geq 2 j" );

  selection = "ptGamma>180. && ht>450. && nJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.";
  dt.drawRegionYields_fromTree("loose_iso"         , "iso"                , selection, 40, 0.    , 10.   , "Photon Charged Isolation"         , "GeV" , cutsLabel, "#ge 2 j" );



  // +++++++++++++++++++++++++
  // +++      Inclusive      +++
  // +++++++++++++++++++++++++


  selection = "ptGamma>180. && ht>200 && met>200. && nJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.";
  dt.set_addOverflow(false); // one sec
  dt.drawRegionYields_fromTree( "inclusive_sietaietaEB", "sietaieta"          , selection, 30, 0.0075, 0.015 , "Photon #sigma_{i#eta i#eta}"      , ""    , "Barrel", " "  );
  dt.drawRegionYields_fromTree( "inclusive_sietaietaEE", "sietaieta"          , selection, 30, 0.02  , 0.035 , "Photon #sigma_{i#eta i#eta}"      , ""    , "Endcap", " "  );
  dt.set_addOverflow(true);


  htMin=200, htMax=-1;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");
  
  selection = "ptGamma>180. && ht>200 && met>200. && nJets>1 && iso<2.5 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.";
  dt.drawRegionYields_fromTree( "inclusive_nVert"      , "nVert"              , selection, 25, 0.5   , 50.5  , "Number of Vertices"               , ""    , cutsLabel );
  dt.drawRegionYields_fromTree( "inclusive_mt2"        , "mt2"                , selection, 40, 0.    , 1000. , "M_{T2} (Photon Removed)"          , "GeV" , cutsLabel );
  dt.drawRegionYields_fromTree( "inclusive_met"        , "met"                , selection, 36, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel );
  dt.drawRegionYields_fromTree( "inclusive_ht"         , "ht"                 , selection, 72, 200.  , 2000. , "H_{T}"                            , "GeV" , cutsLabel );
  dt.drawRegionYields_fromTree( "inclusive_nJets"      , "nJets"              , selection, 10, 1.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel );
  dt.drawRegionYields_fromTree( "inclusive_nBJets"     , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel );
  dt.drawRegionYields_fromTree( "inclusive_ptGamma"    , "ptGamma"            , selection, 36, 180.  , 1080  , "Photon p_{T}"                     , "GeV" , cutsLabel );
  dt.drawRegionYields_fromTree( "inclusive_etaGamma"   , "etaGamma"           , selection, 50, -2.5  , 2.5   , "Photon #eta"                      , ""    , cutsLabel );
  dt.drawRegionYields_fromTree( "inclusive_sietaieta"  , "sietaieta"          , selection, 24, 0.0075, 0.0111, "Photon #sigma_{i#eta i#eta}"      , ""    , cutsLabel );
  dt.drawRegionYields_fromTree( "inclusive_isoRC"      , "gamma_chHadIsoRC[0]", selection, 40, 0.    , 10.   , "Random Cone Isolation"            , "GeV" , cutsLabel );

  selection = "ptGamma>180. && ht>200. && met>200. && nJets>=1 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.";
  dt.drawRegionYields_fromTree("inclusive_loose_iso"         , "iso"                , selection, 40, 0.    , 10.   , "Photon Charged Isolation"         , "GeV" , cutsLabel );


  // +++++++++++++++++++++++++
  // +++      Monojet      +++
  // +++++++++++++++++++++++++

  htMin=200, htMax=-1;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  selection = "ptGamma>180. && ht>200 && met>200. && nJets==1 && iso<2.5 && deltaPhiMin>0.3 && diffMetMht<0.5*met";

  dt.drawRegionYields_fromTree( "monojet_nVert"     , "nVert"              , selection, 25, 0.5   , 50.5  , "Number of Vertices"               , ""    , cutsLabel );
  //dt.drawRegionYields_fromTree( "monojet_mt2"       , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2} (Photon Removed)"          , "GeV" , cutsLabel );
  dt.drawRegionYields_fromTree( "monojet_met"       , "met"                , selection, 36, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel );
  dt.drawRegionYields_fromTree( "monojet_ht"        , "ht"                 , selection, 72, 200.  , 2000. , "H_{T}"                            , "GeV" , cutsLabel );
  dt.drawRegionYields_fromTree( "monojet_nJets"     , "nJets"              , selection, 12, -0.5  , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel );
  dt.drawRegionYields_fromTree( "monojet_nBJets"    , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel );
  dt.drawRegionYields_fromTree( "monojet_ptGamma"   , "ptGamma"            , selection, 36, 180.  , 1080  , "Photon p_{T}"                     , "GeV" , cutsLabel );
  dt.drawRegionYields_fromTree( "monojet_etaGamma"  , "etaGamma"           , selection, 50, -2.5  , 2.5   , "Photon #eta"                      , ""    , cutsLabel );
  dt.drawRegionYields_fromTree( "monojet_sietaieta" , "sietaieta"          , selection, 25, 0.0075, 0.0111, "Photon #sigma_{i#eta i#eta}"      , ""    , cutsLabel );
  dt.drawRegionYields_fromTree( "monojet_isoRC"     , "gamma_chHadIsoRC[0]", selection, 40, 0.    , 10.   , "Random Cone Charged Isolation"    , "GeV" , cutsLabel );

  selection = "ptGamma>180. && ht>200 && met>200. && nJets==1 && deltaPhiMin>0.3 && diffMetMht<0.5*met";
  dt.drawRegionYields_fromTree("monojet_loose_iso"  , "iso"                , selection, 40, 0.    , 10.   , "Photon Charged Isolation"         , "GeV" , cutsLabel );

  return 0;

}

std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units ) {

  std::string cutLabel;
  if( theMax>theMin ) cutLabel = std::string(Form("%.0f < %s < %.0f %s", theMin, name.c_str(), theMax, units.c_str()) );
  else                cutLabel = std::string(Form("%s > %.0f %s", name.c_str(), theMin, units.c_str()) );

  return cutLabel;

}
