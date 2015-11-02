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
std::string getJetCutLabel( int jMin, int jMax, int bMin, int bMax);



int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./drawLostLeptonControlRegionDataMC [configFileName] [lumi/shape]" << std::endl;
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


  std::string mcFile = cfg.getEventYieldDir() + "/llepControlRegion/mc.root";
  std::string dataFile = cfg.getEventYieldDir() + "/llepControlRegion/data.root";

  MT2Analysis<MT2EstimateTree>* mc_   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "llepCR");
  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(dataFile, "llepCR");

  MT2Analysis<MT2EstimateTree>* mc_top      = MT2EstimateTree::makeAnalysisFromInclusiveTree( "Top"   , cfg.regionsSet(), mc_, "id>=300 && id<500" ); 
  MT2Analysis<MT2EstimateTree>* mc_wjets    = MT2EstimateTree::makeAnalysisFromInclusiveTree( "W+jets", cfg.regionsSet(), mc_, "id>=500 && id<600" ); 
  mc_top   ->setColor(855);
  mc_wjets ->setColor(417);

  std::vector< MT2Analysis<MT2EstimateTree>* > mc;
  mc.push_back(mc_top);
  mc.push_back(mc_wjets);

  std::string plotsDir = cfg.getEventYieldDir() + "/llepControlRegion/plotsDataMC";
  if( shapeNorm ) plotsDir += "_shape";


  MT2DrawTools dt(plotsDir, cfg.lumi() );
  dt.set_shapeNorm( shapeNorm );

  dt.set_data( data );
  dt.set_mc( &mc );

  // +++++++++++++++++++++++++
  // +++      Inclusive    +++
  // +++++++++++++++++++++++++

  dt.set_lumi( cfg.lumi() );

  float htMin=200, htMax=-1;
  std::string cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");
  
  int jetMin=2, jetMax=-1;
  int bMin=0, bMax=-1;
  std::string jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);
  
  std::string selection = "weight*(ht>200. && met>200  && nJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.)";

  dt.drawRegionYields_fromTree( "nVert"              , "nVert"                , selection, 100, -0.5    , 99.5  , "Number of Vertices", "" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "nJets"            , "nJets"              , selection, 10, 1.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );



  // +++++++++++++++++++++++++
  // +++      b-veto       +++
  // +++++++++++++++++++++++++

  htMin=200, htMax=-1;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  jetMin=2, jetMax=-1;
  bMin=0, bMax=0;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>200 && met>200. && nJets>1 && nBJets==0 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.)";
  
  dt.drawRegionYields_fromTree( "0b_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "0b_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "0b_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "0b_nJets"            , "nJets"              , selection, 10, 1.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "0b_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );



  // +++++++++++++++++++++++++
  // +++     b-enriched    +++
  // +++++++++++++++++++++++++

  htMin=200, htMax=-1;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");
  
  jetMin=2, jetMax=-1;
  bMin=2, bMax=-1;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>200 && met>200. && nJets>1 && nBJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met &&  mt2>200.)";

  dt.drawRegionYields_fromTree( "2b_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "2b_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "2b_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "2b_nJets"            , "nJets"              , selection, 12, -0.5  , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "2b_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );


  // +++++++++++++++++++++++++
  // +++       Monojet     +++
  // +++++++++++++++++++++++++

  dt.set_lumi( cfg.lumi() );

  htMin=200, htMax=-1;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  jetMin=1, jetMax=1;
  bMin=0, bMax=-1;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>200. && met>200  && nJets==1 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.)";

  //  dt.drawRegionYields_fromTree( "monojet_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_nJets"            , "nJets"              , selection, 12, -0.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );



  // +++++++++++++++++++++++++
  // +++      b-veto       +++
  // +++++++++++++++++++++++++

  htMin=200, htMax=-1;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  jetMin=1, jetMax=1;
  bMin=0, bMax=0;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>200 && met>200. && nJets==1 && nBJets==0 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.)";
  
  //  dt.drawRegionYields_fromTree( "monojet_0b_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_0b_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_0b_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_0b_nJets"            , "nJets"              , selection, 12, -0.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_0b_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );



  // +++++++++++++++++++++++++
  // +++     b-enriched    +++
  // +++++++++++++++++++++++++

  htMin=200, htMax=-1;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  jetMin=1, jetMax=1;
  bMin=2, bMax=-1;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>200 && met>200. && nJets==1 && nBJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met &&  mt2>200.)";

  //dt.drawRegionYields_fromTree( "monojet_2b_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_2b_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_2b_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_2b_nJets"            , "nJets"              , selection, 12, -0.5  , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "monojet_2b_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );


  // +++++++++++++++++++++++++
  // +++     Very Low HT   +++
  // +++++++++++++++++++++++++

  dt.set_lumi( cfg.lumi() );

  htMin=200, htMax=450;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  jetMin=2, jetMax=-1;
  bMin=0, bMax=-1;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>200. && ht<450. && met>200  && nJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.)";

  dt.drawRegionYields_fromTree( "veryLowHT_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_nJets"            , "nJets"              , selection, 10, 1.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );



  // +++++++++++++++++++++++++
  // +++      b-veto       +++
  // +++++++++++++++++++++++++

  htMin=200, htMax=450;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  jetMin=2, jetMax=-1;
  bMin=0, bMax=0;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>200. && ht<450. && met>200  && nJets>1 && nBJets==0 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.)";
  
  dt.drawRegionYields_fromTree( "veryLowHT_0b_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_0b_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_0b_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_0b_nJets"            , "nJets"              , selection, 10, 1.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_0b_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );



  // +++++++++++++++++++++++++
  // +++     b-enriched    +++
  // +++++++++++++++++++++++++

  htMin=200, htMax=450;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  jetMin=2, jetMax=-1;
  bMin=2, bMax=-1;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>200. && ht<450. && met>200  && nJets>1 && nBJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met &&  mt2>200.)";

  dt.drawRegionYields_fromTree( "veryLowHT_2b_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_2b_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_2b_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_2b_nJets"            , "nJets"              , selection, 12, -0.5  , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "veryLowHT_2b_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );


  // +++++++++++++++++++++++++
  // +++     High HT   +++
  // +++++++++++++++++++++++++

  dt.set_lumi( cfg.lumi() );

  htMin=450, htMax=-1;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  jetMin=2, jetMax=-1;
  bMin=0, bMax=-1;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>450. && met>200  && nJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.)";

  dt.drawRegionYields_fromTree( "LowHT_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_nJets"            , "nJets"              , selection, 10, 1.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );



  // +++++++++++++++++++++++++
  // +++      b-veto       +++
  // +++++++++++++++++++++++++


  jetMin=2, jetMax=-1;
  bMin=0, bMax=0;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>450. && met>200  && nJets>1 && nBJets==0 && deltaPhiMin>0.3 && diffMetMht<0.5*met && mt2>200.)";
  
  dt.drawRegionYields_fromTree( "LowHT_0b_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_0b_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_0b_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_0b_nJets"            , "nJets"              , selection, 10, 1.5   , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_0b_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );



  // +++++++++++++++++++++++++
  // +++     b-enriched    +++
  // +++++++++++++++++++++++++

  jetMin=2, jetMax=-1;
  bMin=2, bMax=-1;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>450. && met>200  && nJets>1 && nBJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met &&  mt2>200.)";

  dt.drawRegionYields_fromTree( "LowHT_2b_mt2"              , "mt2"                , selection, 12, 0.    , 600.  , "M_{T2}"          , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_2b_met"              , "met"                , selection, 18, 0.    , 900.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_2b_ht"               , "ht"                 , selection, 37, 200.  , 2050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_2b_nJets"            , "nJets"              , selection, 12, -0.5  , 11.5  , "Number of Jets (p_{T} > 30 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "LowHT_2b_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 20 GeV)", ""    , cutsLabel, jetCutsLabel );


  // +++++++++++++++++++++++++
  // +++     EAs common    +++
  // +++++++++++++++++++++++++

  dt.set_lumi( cfg.lumi() );

  htMin=500, htMax=-1;
  cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  jetMin=2, jetMax=-1;
  bMin=0, bMax=-1;
  jetCutsLabel = getJetCutLabel(jetMin, jetMax, bMin, bMax);

  selection = "weight*(ht>500. && met>250.  && nJets40>1)";

  dt.drawRegionYields_fromTree( "EAs_lepPt"              , "lep_pt"                , selection, 24, 0.    , 600.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "EAs_met"              , "met"                , selection, 27, 200.    , 1550.  , "Missing E_{T}"                    , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "EAs_ht"               , "ht"                 , selection, 52, 450.  , 3050. , "H_{T}"                            , "GeV" , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "EAs_nJets"            , "nJets40"              , selection, 10, 1.5   , 11.5  , "Number of Jets (p_{T} > 40 GeV)"  , ""    , cutsLabel, jetCutsLabel );
  dt.drawRegionYields_fromTree( "EAs_nBJets"           , "nBJets"             , selection, 6 , -0.5  , 5.5   , "Number of b-Jets (p_{T} > 40 GeV)", ""    , cutsLabel, jetCutsLabel );


  return 0;

}

std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units ) {

  std::string cutLabel;
  if( theMax>theMin ) cutLabel = std::string(Form("%.0f < %s < %.0f %s", theMin, name.c_str(), theMax, units.c_str()) );
  else                cutLabel = std::string(Form("%s > %.0f %s", name.c_str(), theMin, units.c_str()) );

  return cutLabel;

}

std::string getJetCutLabel( int jetMin, int jetMax, int bMin, int bMax ) {

  std::string jetCutLabel;
  if( jetMax>jetMin ) jetCutLabel = std::string(Form("%d-%dj, ", jetMin, jetMax));
  else if( jetMax==jetMin ) jetCutLabel = std::string(Form("%dj, ", jetMin));
  else                jetCutLabel = std::string(Form("#geq%dj, ", jetMin));
  
  if( bMax>bMin )     jetCutLabel = jetCutLabel + (std::string)(Form("%d-%db", bMin, bMax));
  else if (bMax<bMin) jetCutLabel = jetCutLabel + (std::string)(Form("#geq%db", bMin));
  else                jetCutLabel = jetCutLabel + (std::string)(Form("%db", bMin));
  
  return jetCutLabel;

}
