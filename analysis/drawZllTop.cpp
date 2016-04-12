#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGaxis.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"

#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"


#include <iostream>
#include "string.h"


#define mt2_cxx
#include "interface/mt2.h"


double lumiErr = 0.046;
bool shapeNorm = false;


//void drawMll( const std::string& outputdir,  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, MT2Analysis<MT2EstimateTree>* data, bool of, float lumi ); 

//void drawStacks(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  MT2Analysis<MT2EstimateTree>* data,const MT2Region thisRegion, std::string cut, float lumi);

MT2Analysis<MT2Estimate>* getEstimate( MT2Config cfg, const std::string& saveName, MT2Analysis<MT2EstimateTree>* top, const std::string& varName, const std::string& selection,  int nBins, double *bins );


void drawYieldsFromHisto( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* DY, MT2Analysis<MT2Estimate>* top,  const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, double* bins, std::string name = "", std::string name_of = "", std::string axisName="", const std::string& units="", const std::string& kinCuts="", const std::string& topoCuts="", bool drawData=0 );


void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* data_of, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields_of, const std::string& saveName, const std::string& varName, const std::string& selection, const std::string& selection_of, int nBins, double* bins, std::string name = "", std::string name_of = "", std::string axisName="", const std::string& units="", const std::string& kinCuts="", const std::string& topoCuts="", bool drawData=0 );

//void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* data_of, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields_of, const std::string& saveName, const std::string& varName, const std::string& selection, const std::string& selection_of, int nBins, float xMin, float xMax, std::string name = "", std::string name_of = "", std::string axisName="", const std::string& units="", const std::string& kinCuts="", const std::string& topoCuts="", bool drawData=0 );

void drawYieldsTopSplit( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* data_of , std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields,  std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields_of, const std::string& saveName, const std::string& varName, const std::string& selection,  const std::string& selection_of, int nBins, float xMin, float xMax,  std::string name, std::string name_of, std::string axisName, const std::string& units , const std::string& kinCuts, const std::string& topoCuts, bool drawData=0 ); 

std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units );

MT2Analysis<MT2Estimate>* getTopBG( MT2Config cfg,  MT2Analysis<MT2EstimateTree>* data_of , MT2Analysis<MT2EstimateTree>* tree_top, MT2Analysis<MT2EstimateTree>*  tree_top_of, const std::string& saveName, const std::string& varName, const std::string& selection,  const std::string& selection_of, int nBins, double *bins, std::string name, std::string name_of, std::string axisName, const std::string& units , const std::string& kinCuts, const std::string& topoCuts ); 

 


int main(int argc, char* argv[]){

  std::string regionsSet = "zurich";


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

  if( argc<2 ) {
    std::cout << "USAGE: ./drawZllControlRegion [configFileName] (lumi/shape)" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  std::string configFileName(argv[1]);
  MT2Config cfg( configFileName);
  regionsSet = cfg.regionsSet();

  std::string outputdir = cfg.getEventYieldDir() + "/zllTop";
  //  std::string outputdir_of = cfg.getEventYieldDir() + "/zllTop_of";

  std::cout << "-> Using regions: " << regionsSet << std::endl;
  MT2DrawTools::setStyle();
 
  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  double intpart;
  double fracpart = modf(cfg.lumi(), &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("_%.0ffb", intpart ) );
  
  system(Form("mkdir -p %s", outputdir.c_str()));


  std::string ZllDir = cfg.getEventYieldDir() + "/zllControlRegion";
  std::string ZllDir_of = cfg.getEventYieldDir() + "/zllControlRegion";

  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "DYJets");
  if( Zll==0 ) {
    std::cout << "-> Please run zllControlRegion first. I need to get the yields from there." << std::endl;   std::cout << "-> Thank you for your cooperation." << std::endl;   exit(197);
  } 

  Zll->setColor(kZinv);
  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str()  ), "QCD");
  qcd->setColor(kQCD);
  MT2Analysis<MT2EstimateTree>* top = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "Top");
  top->setColor(kTop);
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "WJets");
  wjets->setColor(kWJets);

  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ZllDir.c_str()) , "data");
 
  data->setFullName("Data ee/#mu#mu");
  Zll->setFullName("Z+jets ee/#mu#mu");
  wjets->setFullName("W+jets ee/#mu#mu");
  top->setFullName("Top ee/#mu#mu");

  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
  bgYields.push_back( Zll );
  //  bgYields.push_back( qcd );
  //  bgYields.push_back( wjets );
  bgYields.push_back( top );


  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_top; 
  bgYields_top.push_back( top );


  //OPPOSITE FLAVOR TREES
  MT2Analysis<MT2EstimateTree>* Zll_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "DYJets");
  Zll_of->setColor(kZinv);
  MT2Analysis<MT2EstimateTree>* qcd_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str()  ), "QCD");
  qcd_of->setColor(kQCD);
  MT2Analysis<MT2EstimateTree>* top_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "Top");
  top_of->setColor(kTop);
  MT2Analysis<MT2EstimateTree>* wjets_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "WJets");
  wjets_of->setColor(kWJets);
  MT2Analysis<MT2EstimateTree>* data_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data_of.root", ZllDir_of.c_str() ) , "data_of");

  Zll_of->setFullName("Z+jets e#mu");
  wjets_of->setFullName("W+jets e#mu");
  data_of->setFullName("Data e#mu");
  top_of->setFullName("Top e#mu");

  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_of; 
  bgYields_of.push_back( Zll_of );
  // bgYields_of.push_back( qcd_of );
  // bgYields_of.push_back( wjets_of );
  bgYields_of.push_back( top_of );

  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_top_of; 
  bgYields_top_of.push_back( top_of );


  std::string plotsDir = cfg.getEventYieldDir() + "/zllControlRegion/plotsDataMC";
  if( shapeNorm ) plotsDir += "_shape";

  std::string plotsDir_of = cfg.getEventYieldDir() + "/zllControlRegion/plotsDataMC_of";
  if( shapeNorm ) plotsDir_of += "_shape";

  
  MT2DrawTools dt(plotsDir, cfg.lumi() );
  dt.set_shapeNorm( shapeNorm );
  dt.set_data( data );
  dt.set_mc( &bgYields );
  dt.set_lumi( cfg.lumi() );

  MT2DrawTools dt_of(plotsDir_of, cfg.lumi() );
  dt_of.set_shapeNorm( shapeNorm );
  dt_of.set_data( data_of );
  dt_of.set_mc( &bgYields_of );
  dt_of.set_lumi( cfg.lumi() );
  

 
  float htMin=200, htMax=-1;
  std::string cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");

  std::string selection = "( Z_mass>0 &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )*weight_lep0*weight_lep1 ";

  dt_of.drawRegionYields_fromTree( "mll_topOnly"   , "Z_mass"   , selection, 25, 50., 300., "M_{ll}", "GeV", cutsLabel, "#geq1j, #geq0b" );

  //selection = "( abs(Z_mass-91.19)<10 && ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";

  // dt.drawRegionYields_fromTree( "ht"   , "ht"   , selection, 25, 50., 300., "M_{ll}", "GeV", cutsLabel, "#geq1j, #geq0b" );
  //"ht" , selection_peak, selection_of,  28, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b"


  selection = "weight*( Z_mass>0 &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) *weight_lep0*weight_lep1";


  // std::string    selection = "weight*(abs(Z_mass-91.19)<20 && ht>200)";
  // std::string selection_mass_el = "weight*(Z_mass>50 && Z_pt>180 && Z_lepId==11)";
  // std::string selection_mass_mu = "weight*(Z_mass>50 && Z_pt>180 && Z_lepId==13)";


  //Mll with data
  drawYieldsTopSplit( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_split" , "Z_mass" , selection, selection, 20, 50, 150, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b", 1);

  drawYieldsTopSplit( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_split_simulation" , "Z_mass" , selection, selection, 90, 0, 180, "", "","M_{ll}", "GeV", cutsLabel, "#geq1j, #geq0b", 0);
  

  std::string selection_peak = "weight*( ID!=412  && abs(Z_mass-91.19)<10 &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )*weight_lep0*weight_lep1 ";
  std::string selection_of = "weight*( ID!=412 && abs(Z_mass-91.19)<10  &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )*weight_lep0*weight_lep1 ";
  // std::string selection_of = "weight*( ID!=412 && Z_mass > 81.19  &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )*weight_lep0*weight_lep1 ";
  

  double bins_mass[] =  {50,60,70,80,90,100,110,120,130,140,150};
  int size_mass = sizeof(bins_mass)/sizeof(double)-1;

  double bins_ht[] =  {200,450,575,1000,1500,3000};
  int size_ht = sizeof(bins_ht)/sizeof(double)-1;
  double bins_nJets[] = {1,2,4,7,12}; 
  int size_nJets = sizeof(bins_nJets)/sizeof(double)-1;
  double bins_nBJets[] = {0,1,2,3,6}; 
  int size_nBJets = sizeof(bins_nBJets)/sizeof(double)-1;

  //double bins_mt2[] ={200,250,300,350,400,500,800, 1500 };
  double bins_mt2[] ={200,300,400,500, 600, 1500 };
  int size_mt2 = sizeof(bins_mt2)/sizeof(double)-1;
  double bins_mono_ht[] = {200, 250, 350,450, 575, 700, 1000, 1500}; 
  int size_mono_ht = sizeof(bins_mono_ht)/sizeof(double)-1;
  double bins_central[] = {200, 3000}; 
  int size_central = sizeof(bins_central)/sizeof(double)-1;


  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "Z_mass_comp_noTTZll" , "Z_mass" , selection_peak, selection_of, size_mass, bins_mass,  "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "mt2_comp_noTTZll" , "mt2" , selection_peak, selection_of,  size_mt2, bins_mt2, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "ht_comp_noTTZll" , "ht" , selection_peak, selection_of,  size_ht, bins_ht, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "nJets_comp_noTTZll" , "nJets" , selection_peak, selection_of,  size_nJets, bins_nJets ,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "nBJets_comp_noTTZll" , "nBJets" , selection_peak, selection_of,  size_nBJets, bins_nBJets, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
 

  //Old version with equidistant bins
  /*  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "Z_mass_comp_noTTZll" , "Z_mass" , selection_peak, selection_of,  24, 50, 400 , "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "mt2_comp_noTTZll" , "mt2" , selection_peak, selection_of,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "ht_comp_noTTZll" , "ht" , selection_peak, selection_of,  28, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "nJets_comp_noTTZll" , "nJets" , selection_peak, selection_of,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "nBJets_comp_noTTZll" , "nBJets" , selection_peak, selection_of,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
 
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "lep_pt0_comp_noTTZll", "lep_pt0" , selection_peak, selection_of,  20, 0, 500, "", "","Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_top, bgYields_top_of, "lep_pt1_comp_noTTZll", "lep_pt1" , selection_peak, selection_of,  20, 0, 500, "", "","Sub-Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  */





  selection_peak = "weight*( Z_mass >81.19 && Z_mass < 101.19 &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )*weight_lep0*weight_lep1 ";
  selection_of = "weight*( Z_mass > 81.19 &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )*weight_lep0*weight_lep1 ";
  
  MT2Analysis<MT2Estimate>* top_central = getTopBG( cfg, data_of, top, top_of, "central", "ht", selection_peak, selection_of, size_central, bins_central, "", "", "H_{T}", "GeV" , cutsLabel , "#geq1j, #geq0b" ); 



  MT2Analysis<MT2Estimate>* top_mt2 = getTopBG( cfg, data_of, top, top_of, "mt2", "mt2", selection_peak, selection_of, size_mt2, bins_mt2, "", "", "H_{T}", "GeV" , cutsLabel , "#geq1j, #geq0b" ); 
  drawYieldsFromHisto( cfg, data, Zll, top_mt2 , "mt2_comp" , "mt2", selection_peak,  size_mt2, bins_mt2, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");


  MT2Analysis<MT2Estimate>* top_ht = getTopBG( cfg, data_of, top, top_of, "ht", "ht", selection_peak, selection_of, size_ht, bins_ht, "", "", "H_{T}", "GeV" , cutsLabel , "#geq1j, #geq0b" ); 
  drawYieldsFromHisto( cfg, data, Zll, top_ht , "ht_comp" , "ht", selection_peak,  size_ht, bins_ht, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");

  MT2Analysis<MT2Estimate>* top_mono_ht = getTopBG( cfg, data_of, top, top_of, "ht_mono", "ht", selection_peak, selection_of, size_mono_ht, bins_mono_ht, "", "", "H_{T}", "GeV" , cutsLabel , "#geq1j, #geq0b" ); 
  drawYieldsFromHisto( cfg, data, Zll, top_mono_ht , "ht_mono_comp" , "ht", selection_peak,  size_mono_ht, bins_mono_ht, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");

  MT2Analysis<MT2Estimate>* top_nJets = getTopBG( cfg, data_of, top, top_of, "nJets", "nJets", selection_peak, selection_of, size_nJets, bins_nJets, "", "", "Jet Multiplicity", "GeV" , cutsLabel , "#geq1j, #geq0b" ); 
  drawYieldsFromHisto( cfg, data, Zll, top_nJets , "nJets_comp" , "nJets", selection_peak,  size_nJets, bins_nJets, "", "","Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");

 MT2Analysis<MT2Estimate>* top_nBJets = getTopBG( cfg, data_of, top, top_of, "nBJets", "nBJets", selection_peak, selection_of, size_nBJets, bins_nBJets, "", "", "Jet Multiplicity", "GeV" , cutsLabel , "#geq1j, #geq0b" ); 
  drawYieldsFromHisto( cfg, data, Zll, top_nBJets , "nBJets_comp" , "nBJets", selection_peak,  size_nBJets, bins_nBJets, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");



 
  //With background purely from Simulation:

  MT2Analysis<MT2Estimate>* topMC_mt2 = getEstimate( cfg, "topMC_mt2", top, "mt2", selection_peak, size_mt2, bins_mt2 );
  drawYieldsFromHisto( cfg, data, Zll, topMC_mt2 , "mt2_comp_sim" , "mt2" , selection_peak,  size_mt2, bins_mt2, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");

  MT2Analysis<MT2Estimate>* topMC_ht = getEstimate( cfg, "topMC_ht", top, "ht", selection_peak, size_ht, bins_ht );
  drawYieldsFromHisto( cfg, data, Zll, topMC_ht , "ht_comp_sim" , "ht" , selection_peak,  size_ht, bins_ht, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");

  MT2Analysis<MT2Estimate>* topMC_mono_ht = getEstimate( cfg, "topMC_mono_ht", top, "ht", selection_peak, size_mono_ht, bins_mono_ht );
  drawYieldsFromHisto( cfg, data, Zll, topMC_mono_ht , "ht_mono_comp_sim" , "ht" , selection_peak,  size_mono_ht, bins_mono_ht, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");

 

  MT2Analysis<MT2Estimate>* topMC_nJets = getEstimate( cfg, "topMC_nJets", top, "nJets", selection_peak, size_nJets, bins_nJets );
  drawYieldsFromHisto( cfg, data, Zll, topMC_nJets , "nJets_comp_sim" , "nJets" , selection_peak,  size_nJets, bins_nJets, "", "","Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");

  MT2Analysis<MT2Estimate>* topMC_nBJets = getEstimate( cfg, "topMC_nBJets", top, "nBJets", selection_peak, size_nBJets, bins_nBJets );
  drawYieldsFromHisto( cfg, data, Zll, topMC_nBJets , "nBJets_comp_sim" , "nBJets" , selection_peak,  size_nBJets, bins_nBJets, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");

 


  /*
 MT2Analysis<MT2Estimate>* top_nJets = getTopBG( cfg, data_of, top, top_of, "nJets", "nJets", selection_peak, selection_of, size_nJets, bins_nJets, "", "", "Jet Multiplicity", "" , cutsLabel , "#geq1j, #geq0b" ); 
 MT2Analysis<MT2Estimate>* top_nBJets = getTopBG( cfg, data_of, top, top_of, "nBJets", "nBJets", selection_peak, selection_of, size_nBJets, bins_nBJets, "", "", "b Jet Multiplicity", "" , cutsLabel , "#geq1j, #geq0b" ); 


  MT2Analysis<MT2Estimate>* top_mono_ht = getTopBG( cfg, data_of, top, top_of, "ht", "ht", selection_peak, selection_of, size_mono_ht, bins_mono_ht, "", "", "H_{T}", "GeV" , cutsLabel , "#geq1j, #geq0b" ); 



 
  //////*/

  //MT2Analysis<MT2Estimate>* getTopBG( MT2Config cfg,  MT2Analysis<MT2EstimateTree>* data_of , MT2Analysis<MT2EstimateTree>* tree_top, MT2Analysis<MT2EstimateTree>*  tree_top_of, const std::string& saveName, const std::string& varName, const std::string& selection,  const std::string& selection_of, int nBins, double *bins, std::string name, std::string name_of, std::string axisName, const std::string& units , const std::string& kinCuts, const std::string& topoCuts ){

  //drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_comp_noTTZll" , "Z_mass" , selection_peak, selection_of,  24, 50, 400 , "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");



  /*
  selection = "weight*( id==412 &&  Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_TTZll" , "Z_mass" , selection, selection,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_TTZll" , "nJets" , selection, selection,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
   drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_TTZll" , "mt2" , selection, selection,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_TTZll" , "ht" , selection, selection,  24, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_TTZll" , "nBJets" , selection, selection,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  


  selection = "weight*( id==413 &&  Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_TTZqq" , "Z_mass" , selection, selection,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_TTZqq" , "nJets" , selection, selection,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
   drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_TTZqq" , "mt2" , selection, selection,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_TTZqq" , "ht" , selection, selection,  24, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_TTZqq" , "nBJets" , selection, selection,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  


  selection = "weight*( id>=410 && id<=416 && id!=412 && id!=413&&   Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_TTWH" , "Z_mass" , selection, selection,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_TTWH" , "nJets" , selection, selection,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
 
   drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_TTWH" , "mt2" , selection, selection,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_TTWH" , "ht" , selection, selection,  24, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_TTWH" , "nBJets" , selection, selection,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  

  // selection = "weight*( id>=410 && id<=416 &&  Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  // drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_TTV" , "Z_mass" , selection, selection,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  // drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_TTV" , "nJets" , selection, selection,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
 

  selection = "weight*( id==304 &&  Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_TTtoLL" , "Z_mass" , selection, selection,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_TTtoLL" , "nJets" , selection, selection,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");


   drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_TTtoLL" , "mt2" , selection, selection,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_TTtoLL" , "ht" , selection, selection,  24, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_TTtoLL" , "nBJets" , selection, selection,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
   

 selection = "weight*( id>=302 && id<=303 &&  Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_TsingleL" , "Z_mass" , selection, selection,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_TsingleL" , "nJets" , selection, selection,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  
   drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_TsingleL" , "mt2" , selection, selection,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_TsingleL" , "ht" , selection, selection,  24, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_TsingleL" , "nBJets" , selection, selection,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
   

 selection = "weight*( id>=400 && id<=403  &&  Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_Tsingle" , "Z_mass" , selection, selection,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
    drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_Tsingle" , "nJets" , selection, selection,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
 
  
   drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_Tsingle" , "mt2" , selection, selection,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_Tsingle" , "ht" , selection, selection,  24, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_Tsingle" , "nBJets" , selection, selection,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  */




  /*  selection = "weight*( Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  

  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass" , "Z_mass" , selection, selection,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2" , "mt2" , selection, selection,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht" , "ht" , selection, selection,  24, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets" , "nJets" , selection, selection,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets" , "nBJets" , selection, selection,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt0", "lep_pt0" , selection, selection,  24, 0, 600, "", "","Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt1", "lep_pt1" , selection, selection,  24, 0, 600, "", "","Sub-Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  */
  
  /*

  std::string selection_peak = "weight*( abs(Z_mass-91.19)<10 &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  
  std::string selection_of = "weight*( Z_mass>100 && Z_mass < 140 &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  

  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_comp" , "Z_mass" , selection_peak, selection_of,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_comp" , "mt2" , selection_peak, selection_of,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_comp" , "ht" , selection_peak, selection_of,  24, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_comp" , "nJets" , selection_peak, selection_of,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_comp" , "nBJets" , selection_peak, selection_of,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt0_comp", "lep_pt0" , selection_peak, selection_of,  24, 0, 600, "", "","Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt1_comp", "lep_pt1" , selection_peak, selection_of,  24, 0, 600, "", "","Sub-Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");


  selection_of = "weight*( Z_mass>0 && Z_mass < 80 &&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";
  

  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_comp_left" , "Z_mass" , selection_peak, selection_of,  24, 0, 600, "", "","M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_comp_left" , "mt2" , selection_peak, selection_of,  24, 200, 600, "", "","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_comp_left" , "ht" , selection_peak, selection_of,  24, 200, 1600, "", "","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_comp_left" , "nJets" , selection_peak, selection_of,  11, 1, 12,"", "", "Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_comp_left" , "nBJets" , selection_peak, selection_of,  6, 0, 6, "", "","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt0_comp_left", "lep_pt0" , selection_peak, selection_of,  24, 0, 600, "", "","Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt1_comp_left", "lep_pt1" , selection_peak, selection_of,  24, 0, 600, "", "","Sub-Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  */





  /*

  selection = "weight*( Z_lepId==11 && Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) ) ";

  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_el" , "Z_mass" , selection, selection,  24, 0, 600,  "ee", "e#mu", "M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_el" , "mt2" , selection, selection,  24, 200, 600, "ee", "e#mu", "M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_el" , "ht" , selection, selection,  24, 200, 1600, "ee", "e#mu", "H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_el" , "nJets" , selection, selection,  11, 1, 12,  "ee", "e#mu","Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_el" , "nBJets" , selection, selection,  6, 0, 6,  "ee", "e#mu", "b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt0_el", "lep_pt0" , selection, selection,  24, 0, 600,  "ee", "e#mu", "Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt1_el", "lep_pt1" , selection, selection,  24, 0, 600,  "ee", "e#mu", "Sub-Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");


  selection = "weight*( Z_lepId==13 && Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )";

  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "Z_mass_mu" , "Z_mass" , selection, selection,  24, 0, 600, "#mu#mu", "#mue", "M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "mt2_mu" , "mt2" , selection, selection,  24, 200, 600, "#mu#mu", "#mue","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "ht_mu" , "ht" , selection, selection,  24, 200, 1600, "#mu#mu", "#mue","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nJets_mu" , "nJets" , selection, selection,  11, 1, 12, "#mu#mu", "#mue","Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "nBJets_mu" , "nBJets" , selection, selection,  6, 0, 6, "#mu#mu", "#mue","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt0_mu", "lep_pt0" , selection, selection,  24, 0, 600, "#mu#mu", "#mue","Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields_of, "lep_pt1_mu", "lep_pt1" , selection, selection,  24, 0, 600, "#mu#mu", "#mue", "Sub-Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");




  
  //std::string selection_el = "weight*( Z_lepId==11 && Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met  )";
  //std::string selection_mu = "weight*( Z_lepId==13 && Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met  )";
  std::string selection_el = "weight*( Z_lepId==11 && Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met&& (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )";
  std::string selection_mu = "weight*( Z_lepId==13 && Z_mass>0&&  ht>200. && nJets>=1 && met>200. && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met&& (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )";

  drawYields( cfg, data, data_of,  bgYields, bgYields, "Z_mass_sf" , "Z_mass" , selection_el, selection_mu,  24, 0, 600, "ee", "#mu#mu", "M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields, "mt2_sf" , "mt2" , selection_el, selection_mu,  24, 200, 600, "ee", "#mu#mu","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields, "ht_sf" , "ht" , selection_el, selection_mu,  24, 200, 1600, "ee", "#mu#mu","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields, "nJets_sf" , "nJets" , selection_el, selection_mu,  11, 1, 12, "ee", "#mu#mu","Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields, "nBJets_sf" , "nBJets" , selection_el, selection_mu,  6, 0, 6, "ee", "#mu#mu","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields, bgYields, "lep_pt0_sf", "lep_pt0" , selection_el, selection_mu,  24, 0, 600, "ee", "#mu#mu","Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields, bgYields, "lep_pt1_sf", "lep_pt1" , selection_el, selection_mu,  24, 0, 600, "ee", "#mu#mu", "Sub-Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");



  //OF compare emu and mue
  drawYields( cfg, data, data_of,  bgYields_of, bgYields_of, "Z_mass_of" , "Z_mass" , selection_el, selection_mu, 24, 0, 600, "e#mu", "#mue", "M_{ll}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields_of, bgYields_of, "mt2_of" , "mt2" , selection_el, selection_mu,  24, 200, 600, "e#mu", "#mue","M_{T2}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_of, bgYields_of, "ht_of" , "ht" , selection_el, selection_mu,  24, 200, 1600, "e#mu", "#mue","H_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_of, bgYields_of, "nJets_of" , "nJets" , selection_el, selection_mu,  11, 1, 12, "e#mu", "#mue","Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_of, bgYields_of, "nBJets_of" , "nBJets" , selection_el, selection_mu,  6, 0, 6, "e#mu", "#mue","b Jet Multiplicity", "" , cutsLabel, "#geq1j, #geq0b");
  
  drawYields( cfg, data, data_of,  bgYields_of, bgYields_of, "lep_pt0_of", "lep_pt0" , selection_el, selection_mu,  24, 0, 600, "e#mu", "#mue","Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");
  drawYields( cfg, data, data_of,  bgYields_of, bgYields_of, "lep_pt1_of", "lep_pt1" , selection_el, selection_mu,  24, 0, 600, "e#mu", "#mue", "Sub-Leading Lepton p_{T}", "GeV" , cutsLabel, "#geq1j, #geq0b");

  
  */

  return 0;

}





std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units) {

  std::string cutLabel;
  if( theMax>theMin ) cutLabel = std::string(Form("%.0f < %s < %.0f %s", theMin, name.c_str(), theMax, units.c_str()) );
  else                cutLabel = std::string(Form("%s > %.0f %s", name.c_str(), theMin, units.c_str()) );

  return cutLabel;

}


MT2Analysis<MT2Estimate>* getEstimate( MT2Config cfg, const std::string& saveName, MT2Analysis<MT2EstimateTree>* top, const std::string& varName, const std::string& selection,  int nBins, double *bins ){

  TH1::AddDirectory(kTRUE);
  //Get the first region (NB inclusive!)
  std::set<MT2Region> MT2Regions = top->getRegions();
  std::set<MT2Region>::iterator iMT2 = MT2Regions.begin();
  MT2Region thisRegion( (*iMT2) );


  MT2Analysis<MT2Estimate>* estimate = new MT2Analysis<MT2Estimate>( saveName.c_str(), cfg.regionsSet() ); 
  MT2Estimate::rebinYields( estimate, nBins, bins);

  TTree* treeTop = top->get(thisRegion)->tree;
  TH1D* h1_top = new TH1D("h1_top", "", nBins, bins );
  treeTop->Project( "h1_top", varName.c_str(), Form("%f*(%s)", cfg.lumi(), selection.c_str())  );

  for(int i=1; i<=nBins; i++){
    double val = h1_top->GetBinContent(i);
    double err = val * 0.5;
    estimate->get(thisRegion)->yield->SetBinContent(i, val); 
    estimate->get(thisRegion)->yield->SetBinError(i, err);
    std::cout << "Error in bin " << i  << " with value = " << val << " is = " << err << std::endl;
  }

  return estimate;
}


MT2Analysis<MT2Estimate>* getTopBG( MT2Config cfg,  MT2Analysis<MT2EstimateTree>* data_of , MT2Analysis<MT2EstimateTree>* tree_top, MT2Analysis<MT2EstimateTree>*  tree_top_of, const std::string& saveName, const std::string& varName, const std::string& selection,  const std::string& selection_of, int nBins, double *bins, std::string name, std::string name_of, std::string axisName, const std::string& units , const std::string& kinCuts, const std::string& topoCuts ){

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  //The correctly binned estimate of the top background
  //Estimated from opposite flavor events, scaled by the ratio of SF/OF (due to different ranges) and added up to it the TTbarZ events
  MT2Analysis<MT2Estimate>* estimate = new MT2Analysis<MT2Estimate>( saveName.c_str(), cfg.regionsSet() ); 
  MT2Estimate::rebinYields( estimate, nBins, bins);

  //Get the first region (NB inclusive!)
  std::set<MT2Region> MT2Regions = data_of->getRegions();
  std::set<MT2Region>::iterator iMT2 = MT2Regions.begin();
  MT2Region thisRegion( (*iMT2) );

  //SF histogram
  TTree* treeTop = tree_top->get(thisRegion)->tree;
  TH1D* h1_top = new TH1D("h1_top", "", nBins, bins );
  treeTop->Project( "h1_top", varName.c_str(), Form("%f*(%s)*(ID!=412)", cfg.lumi(), selection.c_str())  );

  //OF histogram
  TTree* treeTop_of = tree_top_of->get(thisRegion)->tree;
  TH1D* h1_top_of = new TH1D("h1_top_of", "", nBins, bins );
  treeTop_of->Project( "h1_top_of", varName.c_str(), Form("%f*(%s)*(ID!=412)", cfg.lumi(), selection_of.c_str()) );

  //Ratio of SF/OF for scaling of OF data to SF
  float scaleFactor = h1_top->Integral(1, nBins+1)/h1_top_of->Integral(1, nBins+1);   
  std::cout << "Scale factor = SF/OF = " << h1_top->Integral(1, nBins+1) << " / " << h1_top_of->Integral(1, nBins+1) << " = " << scaleFactor << std::endl;

  //Get data histogram & then scale to area at peak
  TTree* tree_data_of = data_of->get(thisRegion)->tree;
  TH1D* h1_data_of = new TH1D("h1_data_of", "", nBins, bins ); h1_data_of->Sumw2();
  tree_data_of->Project( "h1_data_of", varName.c_str(), selection_of.c_str() );

  for(int i=1; i<=nBins; i++){
    double data = h1_data_of->GetBinContent(i);
    std::cout << "Data Events in bin " << i << " is = " << data << std::endl;
  }

  h1_data_of->Scale( scaleFactor );

  //Get ttZ contribution and add it to the top background
  TH1D* h1_topZ = new TH1D("h1_topZ", "", nBins, bins );
  treeTop->Project( "h1_topZ", varName.c_str(),  Form("%f*(%s)*(ID==412) ", cfg.lumi(), selection.c_str()) );

  h1_top->Sumw2();
  h1_top->Add( h1_topZ );
  h1_top->SetLineColor(kBlue);

  for(int i=1; i<=nBins; i++){
    double data = h1_data_of->GetBinContent(i);
    h1_data_of->SetBinError( i, sqrt( data + 0.1*0.1*data*data) );

    double uncertZ = h1_topZ->GetBinContent(i);
    h1_topZ->SetBinError( i, uncertZ*0.5 );

    double uncert = h1_top->GetBinContent(i);
    h1_top->SetBinError( i, uncert*0.5 );
  }

  h1_data_of->Add( h1_topZ );


  TCanvas* c1 = new TCanvas();
  h1_top->Draw("same");
  h1_data_of->Draw("same");
  // h1_top_of->Draw("same");
  c1->SaveAs("meh.png");
  
  // for(int i=0; i<=nBins; i++){
  //   double val = h1_data_of->GetBinContent(i);
  //   estimate->get(thisRegion)->yield->SetBinContent(i, val);
  // }


  for(int i=1; i<=nBins; i++){
    double val = h1_data_of->GetBinContent(i);
    double err = h1_data_of->GetBinError(i);
    estimate->get(thisRegion)->yield->SetBinContent(i, val); 
    estimate->get(thisRegion)->yield->SetBinError(i, err);
    std::cout << "Error in bin " << i  << " with value = " << val << " is = " << err << std::endl;
  }

  std::string outputdir = cfg.getEventYieldDir() + "/zllTop";
  estimate->writeToFile( outputdir + "/top.root");

  delete h1_top; delete h1_top_of; delete h1_data_of;
  delete h1_topZ;

  return estimate;

  
}




void drawYieldsTopSplit( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* data_of , std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields,  std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields_of, const std::string& saveName, const std::string& varName, const std::string& selection,  const std::string& selection_of, int nBins, float xMin, float xMax,  std::string name, std::string name_of, std::string axisName, const std::string& units, const std::string& kinCuts, const std::string& topoCuts , bool drawData ) {



  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;

  std::vector<int> colors;
  if( bgYields.size()==1 ) { // estimates
    //colors.push_back(855);  
    //   colors.push_back(865); // top 
    colors.push_back(860); // top
    colors.push_back(855); // top
    // colors.push_back(430); 
    //  colors.push_back(418); 
  } else { // mc
    colors.push_back(430); // other=zll
    // colors.push_back(401); // qcd
    // colors.push_back(417); // w+jets
    // colors.push_back(419); // z+jets
    // colors.push_back(865); // top
    colors.push_back(860); // top
    colors.push_back(855); // top
   
    // colors.push_back(855); // top
    // colors.push_back(); // other
  }

  std::string fullPathPlots = cfg.getEventYieldDir() + "/plotsFlavorComparisson";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_data = data->get(thisRegion)->tree;
    TH1D* h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
    tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );
    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.2);

    TTree* tree_data_of = data_of->get(thisRegion)->tree;
    TH1D* h1_data_of = new TH1D("h1_data_of", "", nBins, xMin, xMax );
    tree_data_of->Project( "h1_data_of", varName.c_str(), selection.c_str() );
    TGraphAsymmErrors* gr_data_of = MT2DrawTools::getPoissonGraph(h1_data_of);
    gr_data_of->SetMarkerStyle(20);
    gr_data_of->SetMarkerSize(1.2);


    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<bgYields.size(); ++i ) { 
      if(i==1){
	for(int k=0; k<2; k++){	
	  std::string extension = "ttbar";
	  if( k==0) extension = "ttZ";
	  TTree* tree_mc = (bgYields[i]->get(thisRegion)->tree);
	  std::string thisName = "h1_" + bgYields[i]->getName() + "_" + extension.c_str() ;
	  TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
	  h1_mc->Sumw2();
	  if( k == 1 ) {
	    tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s * (ID>=300&&ID<=500 && ID!=412) )", cfg.lumi(), selection.c_str()) ); 
	  }else {
	    tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s * (ID==412))", cfg.lumi(), selection.c_str()) ); 
	  }
	  // MT2DrawTools::addOverflowSingleHisto(h1_mc);
	  histos_mc.push_back(h1_mc);
	}
      }else {
	TTree* tree_mc = (bgYields[i]->get(thisRegion)->tree);
	std::string thisName = "h1_" + bgYields[i]->getName();
	TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
	h1_mc->Sumw2();
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s)", cfg.lumi(), selection.c_str()) ); 
	// MT2DrawTools::addOverflowSingleHisto(h1_mc);
	histos_mc.push_back(h1_mc);
      }
    } //end histo_mc

    
    TH1D* mc_sum;
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      if( i==0 ) {
        mc_sum = new TH1D( *histos_mc[i] );
        mc_sum->SetName("mc_sum");
      } else {
        mc_sum->Add( histos_mc[i] );
      }
    }


    std::vector< TH1D* > histos_mc_of;
    for( unsigned i=0; i<bgYields_of.size(); ++i ) { 
      TTree* tree_mc = ( bgYields_of[i]->get(thisRegion)->tree );
      std::string thisName = "h1_of_" + bgYields_of[i]->getName();
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_mc->Sumw2();
      tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s)", cfg.lumi(), selection_of.c_str()) ); 
      //	MT2DrawTools::addOverflowSingleHisto(h1_mc);
      histos_mc_of.push_back(h1_mc);
    }//end histos_of

    TH1D* mc_sum_of;
    for( unsigned i=0; i<histos_mc_of.size(); ++i ) { 
      if( i==0 ) {
        mc_sum_of = new TH1D( *histos_mc_of[i] );
        mc_sum_of->SetName("mc_sum");
      } else {
        mc_sum_of->Add( histos_mc_of[i] );
      }

    }

   
    std::cout << "Integrals: " << mc_sum->Integral(0, nBins+1) << "\t" << mc_sum_of->Integral(0, nBins+1) << std::endl;
    float scaleFactor = mc_sum->Integral(0, nBins+1)/mc_sum_of->Integral(0, nBins+1);   
    // if( shapeNorm ) 
    std::cout << "SF: " << scaleFactor << std::endl;

  
    TObject* oldStack = gROOT->FindObject("bgStack");
    if( oldStack ) delete oldStack;
    TH1D* histo_mc;
    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = histos_mc.size() - i - 1;
      histos_mc[index]->SetFillColor( colors[index] );
      // histos_mc[index]->SetLineColor( kBlack );
      histos_mc[index]->SetLineWidth( 0 );
      // if( shapeNorm ) {
      //   histos_mc[index]->Scale( scaleFactor );
      // }
      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack.Add(histos_mc[index]);
    }

    TObject* oldStack_of = gROOT->FindObject("bgStack_of");
    if( oldStack_of ) delete oldStack_of;
    TH1D* histo_mc_of;
    THStack bgStack_of("bgStack_of", "");
    for( unsigned i=0; i<histos_mc_of.size(); ++i ) { 
      //int index = histos_mc_of.size() - i - 1;
      int index = i;
      // histos_mc_of[index]->SetFillColor( colors[index] );
      // histos_mc_of[index]->SetFillStyle( 3444 );
      histos_mc_of[index]->SetLineColor( kBlack );
      histos_mc_of[index]->SetLineWidth( 2 );
      //histos_mc_of[index]->SetLineStyle( i+1 );
      if( shapeNorm ) {
        histos_mc_of[index]->Scale( scaleFactor );
      }
      if(i==0) histo_mc_of = (TH1D*) histos_mc_of[index]->Clone("histo_mc_of");
      else histo_mc_of->Add(histos_mc_of[index]);
      bgStack_of.Add(histos_mc_of[index]);
    }

     histos_mc_of[0]->SetLineColor( kGray );

    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(histo_mc, histo_mc_of);
    TH1D* h_ratio = new TH1D("h_ratio","",nBins, xMin, xMax);
    h_ratio->Add(histo_mc);
    h_ratio->Divide(histo_mc_of);
    h_ratio->SetMarkerSize(1.4);
    h_ratio->SetMarkerStyle(20);

    // TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
    // TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data_of, histo_mc);
    
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);
    
    TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TGraphErrors* SFFitBand = (fSF) ? MT2DrawTools::getSFFitBand(fSF, xMin, xMax) : 0;

    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );


    TCanvas* c1 = new TCanvas(Form("c1_%s", iMT2->getName().c_str()), "", 600, 600);
    c1->cd();
    
    TCanvas* c1_log = new TCanvas(Form("c1_log_%s", iMT2->getName().c_str()), "", 600, 600);

    float yMaxScale = 1.1;
    float yMax1 = histo_mc_of->GetMaximum()*yMaxScale ;
    //  float yMax1 = h1_data->GetMaximum()*yMaxScale ;
    float yMax2 = yMaxScale*(histo_mc_of->GetMaximum() + sqrt(histo_mc_of->GetMaximum()));
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
    }

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();  

    // TPad* pad1 = 0;
    // pad1 = MT2DrawTools::getCanvasMainPad();
    // pad1->Draw();
    // pad1->cd();

    h2_axes->Draw();

   
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.1, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();
    c1_log->SetLogy();   
    //  TPad* pad1_log = 0;
    //  pad1_log = MT2DrawTools::getCanvasMainPad( true );
    //  pad1_log->Draw();
    //  pad1_log->cd();

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
      //  pad1->cd();
      c1->cd();
      regionText->Draw("same");
      c1_log->cd();
      // pad1_log->cd();
      regionText->Draw("same");
    }


    
    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      if( scaleFactor!=1. ) {
	normText->AddText( Form("#splitline{MC scaled}{by %.2f}", scaleFactor) );
      }
      //  pad1->cd();
      c1->cd();
      normText->Draw("same");
      // pad1_log->cd();
      c1_log->cd();
      normText->Draw("same");
    }

    TLegend* legend = new TLegend( 0.65, 0.9-(bgYields.size()+3)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    if( drawData==true )   legend->AddEntry( gr_data, "Data ee/#mu#mu", "P" );

    for( unsigned i=0; i<=histos_mc.size(); ++i ) {
      if( i == 0){
	TH1D* h1_bg1 = bgYields[i]->get(thisRegion)->yield;
	h1_bg1->SetFillColor( colors[i] );
	legend->AddEntry( h1_bg1, bgYields[i]->getFullName().c_str(), "F" );
      } if( i == 1) {
	TH1D* h1_bg1 = new TH1D("h1_bg1","",10,0,2);
	h1_bg1->SetFillColor( colors[2] );
	legend->AddEntry( h1_bg1, "Top ee/#mu#mu", "F" );
      } else if( i == 2){
	TH1D* h1_bg1 = new TH1D("h1_bg1","",10,0,2);
	h1_bg1->SetFillColor( colors[1]);
	legend->AddEntry( h1_bg1, "tt+Z ee/#mu#mu", "F" );
      } //else if( i ==3){
      /// TH1D* h1_bg1 = new TH1D("h1_bg1","",10,0,2);
      // h1_bg1->SetFillColor( colors[i]);
      // legend->AddEntry( h1_bg1, "other Top SF", "F" );
      // }
      //TH1D* h1_bg1 = bgYields[index]->get(thisRegion)->yield;
      // legend->AddEntry( h1_bg1, bgYields[index]->getFullName().c_str(), "F" );
    }

    if(drawData==true) {
      legend->AddEntry( gr_data_of, "Data e#mu", "P" );
      legend->AddEntry( histos_mc_of[1] , "MC e#mu" , "L" );
    } else {
    legend->AddEntry( histos_mc_of[1], "Top e#mu" , "L" );
    legend->AddEntry( histos_mc_of[0], "Z+jets e#mu" , "L" );
    }
   
    // for( unsigned i=0; i<bgYields.size(); ++i ) {
    // std::cout <<  bgYields_of[i]->getFullName() << std::endl;
    //      legend->AddEntry( histos_mc[i], Form("%s %s", bgYields[i]->getFullName().c_str(), name.c_str()  ) , "F" );
    //  legend->AddEntry( histos_mc_of[i], Form("%s %s", bgYields_of[i]->getFullName().c_str(), name_of.c_str() ) , "L" );
   // }
    // legend->AddEntry( mcBand, "MC Uncert.", "F" );


    TPaveText* fitText = (fSF) ? MT2DrawTools::getFitText( fSF ) : 0;


    float yMinR=0.0;
    float yMaxR=2.0;

    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
  
    gr_data->SetMarkerSize(1);
    gr_data_of->SetMarkerSize(1);
    gr_data_of->SetMarkerStyle(24);


    c1->cd();
    // pad1->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    bgStack.Draw("L same");
 
    // bgStack_of.GetStack()->Last()->Draw("L same");
    if( drawData == false ) {
      histos_mc_of[0]->Draw("L same");
      histos_mc_of[1]->Draw("L same");
      MT2DrawTools::addLabels( gPad, cfg.lumi(), "CMS Simulation" );
    } else {
      bgStack_of.GetStack()->Last()->Draw("L same");
      gr_data->Draw("p same");
      gr_data_of->Draw("p same");
      MT2DrawTools::addLabels( gPad, cfg.lumi(), "CMS Preliminary" );
    }
    // if( !shapeNorm )
    //     fitText->Draw("same");
    // ratioText->Draw("same");
  
    gPad->RedrawAxis();

    c1_log->cd();
    // pad1_log->cd();
    legend->Draw("same");
    bgStack.Draw("L same");
    bgStack.Draw("histo same");
    if( drawData == false ) {
      histos_mc_of[0]->Draw("L same");
      histos_mc_of[1]->Draw("L same");
      MT2DrawTools::addLabels( gPad, cfg.lumi(), "CMS Simulation" );
    } else {
      bgStack_of.GetStack()->Last()->Draw("L same");
      gr_data->Draw("p same");
      gr_data_of->Draw("p same");
      MT2DrawTools::addLabels( gPad, cfg.lumi(), "CMS Preliminary" );
    }

   //  if( !shapeNorm )
    // fitText->Draw("same");
    //  ratioText->Draw("same");

    gPad->RedrawAxis();

   /*
    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    TLine* lineSF = MT2DrawTools::getSFLine(integral_data, integral_mc, xMin, xMax);
    TGraphErrors* SFband = MT2DrawTools::getSFBand(integral_data, error_data, integral_mc, error_mc, xMin, xMax);
    */

    /*
    c1->cd();
    TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
    pad2->Draw();
    pad2->cd();

    h2_axes_ratio->Draw("");
 
    /*  line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same");
    */
    /*
    lineCentral->Draw("same");
    if( !shapeNorm ){

      systBand->Draw("3,same");
      lineCentral->Draw("same");

      SFFitBand->Draw("3,same");
      fSF->Draw("same");
    }

    h_ratio->Draw("PE,same"); 
    //    g_ratio->Draw("PE,same");    
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

      SFFitBand->Draw("3,same");
      fSF->Draw("same");
    }
    */

    /*
    line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same"); */

    //  h_ratio->Draw("PE,same");
    //    g_ratio->Draw("PE,same");
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

    delete h1_data; delete h1_data_of; delete h_ratio;
    delete g_ratio; delete gr_data_of;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

    for( unsigned i=0; i<histos_mc_of.size(); ++i )
      delete histos_mc_of[i];

  }// for MT2 regions

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
    colors.push_back(430); // other=zll
    // colors.push_back(401); // qcd
    // colors.push_back(417); // w+jets
    // colors.push_back(419); // z+jets
    colors.push_back(855); // top
    // colors.push_back(); // other
  }


  std::string fullPathPlots = cfg.getEventYieldDir() + "/plotsFlavorComparisson";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_data = data->get(thisRegion)->tree;
    TH1D* h1_data = new TH1D("h1_data", "", nBins, bins );
    tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );
    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.2);

    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<bgYields.size(); ++i ) { 
      TTree* tree_mc = (bgYields[i]->get(thisRegion)->tree);
      std::string thisName = "h1_" + bgYields[i]->getName();
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, bins );
      h1_mc->Sumw2();
      tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s)", cfg.lumi(), selection.c_str()) ); 
      MT2DrawTools::addOverflowSingleHisto(h1_mc);
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

    std::vector< TH1D* > histos_mc_of;
    for( unsigned i=0; i<bgYields_of.size(); ++i ) { 
      TTree* tree_mc = (bgYields_of[i]->get(thisRegion)->tree);
      std::string thisName = "h1_of_" + bgYields_of[i]->getName();
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, bins );
      h1_mc->Sumw2();
      tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s)", cfg.lumi(), selection_of.c_str()) ); 
      MT2DrawTools::addOverflowSingleHisto(h1_mc);
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



    std::cout << "Integrals: " << mc_sum->Integral(0, nBins+1) << "\t" << mc_sum_of->Integral(0, nBins+1) << std::endl;
    float scaleFactor = mc_sum->Integral(0, nBins+1)/mc_sum_of->Integral(0, nBins+1);   
    // if( shapeNorm ) 
      std::cout << "SF: " << scaleFactor << std::endl;

    TObject* oldStack = gROOT->FindObject("bgStack");
    if( oldStack ) delete oldStack;
    TH1D* histo_mc;
    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = bgYields.size() - i - 1;
      histos_mc[index]->SetFillColor( colors[index] );
      // histos_mc[index]->SetLineColor( kBlack );
      histos_mc[index]->SetLineWidth( 0 );
      // if( shapeNorm ) {
      //   histos_mc[index]->Scale( scaleFactor );
      // }
      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack.Add(histos_mc[index]);
    }


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

    TTree* tree_data_of = data_of->get(thisRegion)->tree;
    TH1D* h1_data_of = new TH1D("h1_data_of", "", nBins, bins );
    tree_data_of->Project( "h1_data_of", varName.c_str(), selection_of.c_str() );
    if( shapeNorm )   h1_data_of->Scale( scaleFactor);
    TH1D* gr_data_of = new TH1D("gr_data_of", "", nBins, bins);
    gr_data_of->Add(h1_data_of);
    //TGraphAsymmErrors* gr_data_of = MT2DrawTools::getPoissonGraph(h1_data_of);
    gr_data_of->SetMarkerStyle(24);
    gr_data_of->SetMarkerSize(1);


    //TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph( h1_data_of, histo_mc_of);
    TH1D* g_ratio = new TH1D("g_ratio", "", nBins, bins);
    g_ratio->Add(h1_data_of);
    g_ratio->Divide( histo_mc );
    g_ratio->SetMarkerStyle( 24 ); g_ratio->SetMarkerSize(1.4);

    //    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(histo_mc, histo_mc_of);
    TH1D* h_ratio = new TH1D("h_ratio","",nBins, bins);
    h_ratio->Add(histo_mc);
    h_ratio->Divide(histo_mc_of);
    h_ratio->SetMarkerSize(1.4);
    h_ratio->SetMarkerStyle(20);

    //    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
    // TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data_of, histo_mc);
    
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);
    
    // TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TF1* fSF = 0;// MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TGraphErrors* SFFitBand = (fSF) ? MT2DrawTools::getSFFitBand(fSF, xMin, xMax) : 0;

    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );


    TCanvas* c1 = new TCanvas(Form("c1_%s", iMT2->getName().c_str()), "", 600, 600);
    c1->cd();
    
    TCanvas* c1_log = new TCanvas(Form("c1_log_%s", iMT2->getName().c_str()), "", 600, 600);

    float yMaxScale = 1.1;
    float yMax1 = histo_mc_of->GetMaximum()*yMaxScale ;
    //    float yMax1 = h1_data->GetMaximum()*yMaxScale ;
    float yMax2 = yMaxScale*(histo_mc_of->GetMaximum() + sqrt(histo_mc_of->GetMaximum()));
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


    
    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      if( scaleFactor!=1. ) {
	normText->AddText( Form("#splitline{e#mu scaled}{by %.2f}", scaleFactor) );
      }
      pad1->cd();
      normText->Draw("same");
      pad1_log->cd();
      normText->Draw("same");
    }

    TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+2)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    // legend->AddEntry( gr_data, "Data", "P" );
    if(drawData == true) legend->AddEntry( gr_data_of, "Data e#mu", "EP" );
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
     legend->AddEntry( histos_mc[i], Form("%s %s", bgYields[i]->getFullName().c_str(), name.c_str()  ) , "F" );
      legend->AddEntry( histos_mc_of[i], Form("%s %s", bgYields_of[i]->getFullName().c_str(), name_of.c_str() ) , "L" );
    }
    //legend->AddEntry( mcBand, "MC Uncert.", "F" );


    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    
    TPaveText* fitText = (fSF) ? MT2DrawTools::getFitText( fSF ) : 0;


    float yMinR=0.0;
    float yMaxR=2.0;

    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
 
    if( drawData == false ) h2_axes_ratio->SetYTitle( "#frac{ee+#mu#mu}{e#mu+#mue}" );

    c1->cd();
    pad1->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    bgStack.Draw("L same");
    bgStack_of.Draw("L same");
    // gr_data->Draw("p same");
    if(drawData==true) gr_data_of->Draw("EP same");
    if( drawData == false ) {
      MT2DrawTools::addLabels( pad1, cfg.lumi(), "CMS Simulation" );
    } else {
      MT2DrawTools::addLabels( pad1, cfg.lumi(), "CMS Preliminary" );
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
    bgStack_of.Draw("L same");
    //  gr_data->Draw("p same");
    if(drawData==true) gr_data_of->Draw("PE same");
    //  labelTop->Draw("same");
    //if( !shapeNorm )
      // fitText->Draw("same");
    //  ratioText->Draw("same");

    if( drawData == false ) {
      MT2DrawTools::addLabels( pad1_log, cfg.lumi(), "CMS Simulation" );
    } else {
      MT2DrawTools::addLabels( pad1_log, cfg.lumi(), "CMS Preliminary" );
    }
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


    h_ratio->Draw("PE,same"); 
    if(drawData==true) g_ratio->Draw("PE,same");    
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

    h_ratio->Draw("PE,same");
    if(drawData==true) g_ratio->Draw("PE,same");
   
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

    delete h1_data; delete h1_data_of; delete h_ratio;
    delete g_ratio; delete gr_data_of;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

    for( unsigned i=0; i<histos_mc_of.size(); ++i )
      delete histos_mc_of[i];

  }// for MT2 regions

}













void drawYieldsFromHisto( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* DY , MT2Analysis<MT2Estimate>* top , const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, double* bins,  std::string name, std::string name_of, std::string axisName, const std::string& units , const std::string& kinCuts, const std::string& topoCuts, bool drawData ) {

  float xMax = bins[nBins];
  float xMin = bins[0];

  //  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;

  std::vector<int> colors;
  colors.push_back(430); // other=zll
  // colors.push_back(401); // qcd
  // colors.push_back(417); // w+jets
  // colors.push_back(419); // z+jets
  colors.push_back(855); // top
  // colors.push_back(); // other
 

  std::string fullPathPlots = cfg.getEventYieldDir() + "/plotsComparisson";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_data = data->get(thisRegion)->tree;
    TH1D* h1_data = new TH1D("h1_data", "", nBins, bins );
    tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );
    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.2);

    std::vector< TH1D* > histos_mc;
  
    TTree* tree_mc = (DY->get(thisRegion)->tree);
    std::string thisName = "h1_" + DY->getName();
    TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, bins );
    h1_mc->Sumw2();
    tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s)", cfg.lumi(), selection.c_str()) ); 
    MT2DrawTools::addOverflowSingleHisto(h1_mc);
    histos_mc.push_back(h1_mc);

    histos_mc.push_back( top->get(thisRegion)->yield );
 

    TH1D* mc_sum;
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      if( i==0 ) {
        mc_sum = new TH1D( *histos_mc[i] );
        mc_sum->SetName("mc_sum");
      } else {
        mc_sum->Add( histos_mc[i] );
      }
    }


    std::cout << "Integrals: " << h1_data->Integral(0, nBins+1) << "\t" << mc_sum->Integral(0, nBins+1) << std::endl;
    float scaleFactor = h1_data->Integral(0, nBins+1)/mc_sum->Integral(0, nBins+1);   
    // if( shapeNorm ) 
      std::cout << "SF: " << scaleFactor << std::endl;

    TObject* oldStack = gROOT->FindObject("bgStack");
    if( oldStack ) delete oldStack;
    TH1D* histo_mc;
    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = histos_mc.size() - i - 1;
      histos_mc[index]->SetFillColor( colors[index] );
      // histos_mc[index]->SetLineColor( kBlack );
      histos_mc[index]->SetLineWidth( 0 );
      // if( shapeNorm ) {
      //   histos_mc[index]->Scale( scaleFactor );
      // }
      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack.Add(histos_mc[index]);
    }




    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
    // TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data_of, histo_mc);
    
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);
    
    // TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TF1* fSF = 0;// MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TGraphErrors* SFFitBand = (fSF) ? MT2DrawTools::getSFFitBand(fSF, xMin, xMax) : 0;

    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );


    TCanvas* c1 = new TCanvas(Form("c1_%s", iMT2->getName().c_str()), "", 600, 600);
    c1->cd();
    
    TCanvas* c1_log = new TCanvas(Form("c1_log_%s", iMT2->getName().c_str()), "", 600, 600);

    float yMaxScale = 1.1;
    float yMax1 = h1_data->GetMaximum()*yMaxScale ;
    float yMax2 = yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
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


    std::string yAxisTitle = (std::string)(Form("Events"));
 

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


    
    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      if( scaleFactor!=1. ) {
	normText->AddText( Form("#splitline{e#mu scaled}{by %.2f}", scaleFactor) );
      }
      pad1->cd();
      normText->Draw("same");
      pad1_log->cd();
      normText->Draw("same");
    }

    TLegend* legend = new TLegend( 0.7, 0.9-(histos_mc.size()+2)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( gr_data, "Data", "P" ); 
    legend->AddEntry( histos_mc[0], "Z+jets", "F" );
    legend->AddEntry( histos_mc[1], "Top", "F" );


    //legend->AddEntry( mcBand, "MC Uncert.", "F" );


    TPaveText* fitText = (fSF) ? MT2DrawTools::getFitText( fSF ) : 0;


    float yMinR=0.0;
    float yMaxR=2.0;

    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
 
 
    c1->cd();
    pad1->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    MT2DrawTools::addLabels( pad1, cfg.lumi(), "CMS Preliminary" );

    //  if( !shapeNorm )
    //   fitText->Draw("same");
    // ratioText->Draw("same");
  
    gPad->RedrawAxis();

    c1_log->cd();
    pad1_log->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    MT2DrawTools::addLabels( pad1_log, cfg.lumi(), "CMS Preliminary" );
    //if( !shapeNorm )
    // fitText->Draw("same");
    //  ratioText->Draw("same");
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

    delete h1_data;
    delete g_ratio; 
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

  }// for MT2 regions

}






























