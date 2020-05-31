// ROOT includes
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TAxis.h"
#include "TString.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TKey.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TChain.h"
#include "TList.h"
#include "TObject.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"

// RooFit includes
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"
#include "RooCmdArg.h"

// RooStats includes
#include "RooWorkspace.h"

// standard includes
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>


#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"

#define mt2_cxx
#include "interface/mt2.h"

using namespace std;
using namespace RooFit;

bool doMT2binning = true;

bool doSingleLbin = true;
bool doLLbin = true;
bool doBBbin = true;



int main(int argc, char* argv[]){

  std::string regionsSet = "zurich";

  std::string configFileName(argv[1]);
  MT2Config cfg( configFileName);
  regionsSet = cfg.regionsSet();

  std::string ggDir = cfg.getEventYieldDir() + "/diPhotonControlRegion/";

  std::string outputdir = cfg.getEventYieldDir() + "/WS";
  system(Form("mkdir -p %s", outputdir.c_str()));

  bool doWH = false;
  bool doHZ = false;
  bool doHH = false;
  bool doT2bH = false;

  if( argc>2 ) {
   std::string inputName(argv[2]);
   if( inputName == "T2bH" )
      doT2bH = true;
    else if( inputName == "HH" )
      doHH = true;
    else if( inputName == "HZ" )
      doHZ = true;
    else if( inputName == "WH" )
      doWH = true;
    else 
      std::cout << "Unknown new physics you are looking for is not here" << std::endl;
  }

  bool  doISRsyst = false;
  if( doT2bH )
    doISRsyst = true;


  vector<std::string> signalList_T2bH = { "SMS_T2bH_mSbottom250_mLSP100", "SMS_T2bH_mSbottom250_mLSP1", "SMS_T2bH_mSbottom250_mLSP50", "SMS_T2bH_mSbottom300_mLSP100", "SMS_T2bH_mSbottom300_mLSP150", "SMS_T2bH_mSbottom300_mLSP1", "SMS_T2bH_mSbottom300_mLSP50", "SMS_T2bH_mSbottom350_mLSP100", "SMS_T2bH_mSbottom350_mLSP150", "SMS_T2bH_mSbottom350_mLSP1", "SMS_T2bH_mSbottom350_mLSP200", "SMS_T2bH_mSbottom350_mLSP50", "SMS_T2bH_mSbottom400_mLSP100", "SMS_T2bH_mSbottom400_mLSP150", "SMS_T2bH_mSbottom400_mLSP1", "SMS_T2bH_mSbottom400_mLSP200", "SMS_T2bH_mSbottom400_mLSP250", "SMS_T2bH_mSbottom400_mLSP50", "SMS_T2bH_mSbottom450_mLSP1", "SMS_T2bH_mSbottom450_mLSP100", "SMS_T2bH_mSbottom450_mLSP150", "SMS_T2bH_mSbottom450_mLSP200", "SMS_T2bH_mSbottom450_mLSP250", "SMS_T2bH_mSbottom450_mLSP300", "SMS_T2bH_mSbottom450_mLSP50", "SMS_T2bH_mSbottom500_mLSP100", "SMS_T2bH_mSbottom500_mLSP150", "SMS_T2bH_mSbottom500_mLSP1", "SMS_T2bH_mSbottom500_mLSP200", "SMS_T2bH_mSbottom500_mLSP250", "SMS_T2bH_mSbottom500_mLSP300", "SMS_T2bH_mSbottom500_mLSP50", "SMS_T2bH_mSbottom600_mLSP1", "SMS_T2bH_mSbottom600_mLSP100", "SMS_T2bH_mSbottom600_mLSP200", "SMS_T2bH_mSbottom600_mLSP300" };

  vector<std::string> signalList_WH = { "SMS_TChiWH_HToGG_127_1","SMS_TChiWH_HToGG_150_1","SMS_TChiWH_HToGG_150_24","SMS_TChiWH_HToGG_175_1","SMS_TChiWH_HToGG_175_25","SMS_TChiWH_HToGG_175_49","SMS_TChiWH_HToGG_200_1","SMS_TChiWH_HToGG_200_25","SMS_TChiWH_HToGG_200_50","SMS_TChiWH_HToGG_200_74" };

  vector<std::string> signalList_HZ = { "SMS_TChiHZ_HToGG_m127", "SMS_TChiHZ_HToGG_m150", "SMS_TChiHZ_HToGG_m175", "SMS_TChiHZ_HToGG_m200" };
  
//  vector<std::string> signalList_HZ = { "SMS_TChiHZ_HToGG_m127", "SMS_TChiHZ_HToGG_m150", "SMS_TChiHZ_HToGG_m175", "SMS_TChiHZ_HToGG_m200", "SMS_TChiHZ_HToGG_m225", "SMS_TChiHZ_HToGG_m250", "SMS_TChiHZ_HToGG_m275", "SMS_TChiHZ_HToGG_m300", "SMS_TChiHZ_HToGG_m325", "SMS_TChiHZ_HToGG_m350", "SMS_TChiHZ_HToGG_m375", "SMS_TChiHZ_HToGG_m400", "SMS_TChiHZ_HToGG_m425", "SMS_TChiHZ_HToGG_m450", "SMS_TChiHZ_HToGG_m475", "SMS_TChiHZ_HToGG_m500", "SMS_TChiHZ_HToGG_m525", "SMS_TChiHZ_HToGG_m550", "SMS_TChiHZ_HToGG_m575", "SMS_TChiHZ_HToGG_m600", "SMS_TChiHZ_HToGG_m625", "SMS_TChiHZ_HToGG_m650", "SMS_TChiHZ_HToGG_m675", "SMS_TChiHZ_HToGG_m700", "SMS_TChiHZ_HToGG_m725", "SMS_TChiHZ_HToGG_m750", "SMS_TChiHZ_HToGG_m775", "SMS_TChiHZ_HToGG_m800", "SMS_TChiHZ_HToGG_m825", "SMS_TChiHZ_HToGG_m850", "SMS_TChiHZ_HToGG_m875", "SMS_TChiHZ_HToGG_m900", "SMS_TChiHZ_HToGG_m925", "SMS_TChiHZ_HToGG_m950", "SMS_TChiHZ_HToGG_m975", "SMS_TChiHZ_HToGG_m1000" };

  vector<std::string> signalList_HH_HH0p25 = { "SMS_TChiHH_HToGG_m127_HH0p25", "SMS_TChiHH_HToGG_m150_HH0p25", "SMS_TChiHH_HToGG_m175_HH0p25", "SMS_TChiHH_HToGG_m200_HH0p25"};

  //  vector<std::string> signalList_HH_HH0p25 = { "SMS_TChiHH_HToGG_m127_HH0p25", "SMS_TChiHH_HToGG_m150_HH0p25", "SMS_TChiHH_HToGG_m175_HH0p25", "SMS_TChiHH_HToGG_m200_HH0p25", "SMS_TChiHH_HToGG_m225_HH0p25", "SMS_TChiHH_HToGG_m250_HH0p25", "SMS_TChiHH_HToGG_m275_HH0p25", "SMS_TChiHH_HToGG_m300_HH0p25", "SMS_TChiHH_HToGG_m325_HH0p25", "SMS_TChiHH_HToGG_m350_HH0p25", "SMS_TChiHH_HToGG_m375_HH0p25", "SMS_TChiHH_HToGG_m400_HH0p25", "SMS_TChiHH_HToGG_m425_HH0p25", "SMS_TChiHH_HToGG_m450_HH0p25", "SMS_TChiHH_HToGG_m475_HH0p25", "SMS_TChiHH_HToGG_m500_HH0p25", "SMS_TChiHH_HToGG_m525_HH0p25", "SMS_TChiHH_HToGG_m550_HH0p25", "SMS_TChiHH_HToGG_m575_HH0p25", "SMS_TChiHH_HToGG_m600_HH0p25", "SMS_TChiHH_HToGG_m625_HH0p25", "SMS_TChiHH_HToGG_m650_HH0p25", "SMS_TChiHH_HToGG_m675_HH0p25", "SMS_TChiHH_HToGG_m700_HH0p25", "SMS_TChiHH_HToGG_m725_HH0p25", "SMS_TChiHH_HToGG_m750_HH0p25", "SMS_TChiHH_HToGG_m775_HH0p25", "SMS_TChiHH_HToGG_m800_HH0p25", "SMS_TChiHH_HToGG_m825_HH0p25", "SMS_TChiHH_HToGG_m850_HH0p25", "SMS_TChiHH_HToGG_m875_HH0p25", "SMS_TChiHH_HToGG_m900_HH0p25", "SMS_TChiHH_HToGG_m925_HH0p25", "SMS_TChiHH_HToGG_m950_HH0p25", "SMS_TChiHH_HToGG_m975_HH0p25", "SMS_TChiHH_HToGG_m1000_HH0p25" };


  vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127", "SMS_TChiHH_HToGG_m150", "SMS_TChiHH_HToGG_m175", "SMS_TChiHH_HToGG_m200"};

  //vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127", "SMS_TChiHH_HToGG_m150", "SMS_TChiHH_HToGG_m175", "SMS_TChiHH_HToGG_m200", "SMS_TChiHH_HToGG_m225", "SMS_TChiHH_HToGG_m250", "SMS_TChiHH_HToGG_m275", "SMS_TChiHH_HToGG_m300", "SMS_TChiHH_HToGG_m325", "SMS_TChiHH_HToGG_m350", "SMS_TChiHH_HToGG_m375", "SMS_TChiHH_HToGG_m400", "SMS_TChiHH_HToGG_m425", "SMS_TChiHH_HToGG_m450", "SMS_TChiHH_HToGG_m475", "SMS_TChiHH_HToGG_m500", "SMS_TChiHH_HToGG_m525", "SMS_TChiHH_HToGG_m550", "SMS_TChiHH_HToGG_m575", "SMS_TChiHH_HToGG_m600", "SMS_TChiHH_HToGG_m625", "SMS_TChiHH_HToGG_m650", "SMS_TChiHH_HToGG_m675", "SMS_TChiHH_HToGG_m700", "SMS_TChiHH_HToGG_m725", "SMS_TChiHH_HToGG_m750", "SMS_TChiHH_HToGG_m775", "SMS_TChiHH_HToGG_m800", "SMS_TChiHH_HToGG_m825", "SMS_TChiHH_HToGG_m850", "SMS_TChiHH_HToGG_m875", "SMS_TChiHH_HToGG_m900", "SMS_TChiHH_HToGG_m925", "SMS_TChiHH_HToGG_m950", "SMS_TChiHH_HToGG_m975", "SMS_TChiHH_HToGG_m1000" };





  //BACKGROUNDS
  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_qcd.root", ggDir.c_str()  ), "qcd");
  MT2Analysis<MT2EstimateTree>* diPhoton = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diPhoton.root", ggDir.c_str() ), "diPhoton");
  MT2Analysis<MT2EstimateTree>* gjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_gjets.root", ggDir.c_str() ), "gjets");
  MT2Analysis<MT2EstimateTree>* higgs = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir.c_str() ), "higgs");
  
  //DATA
  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ggDir.c_str() ), "diPhoton_data");

  //SIGNAL
  string sigName;
  vector<std::string> signalList;

  if( doT2bH ){        signalList = signalList_T2bH;     sigName = "SMS_T2bH_mSbottom";
  }else if( doHZ ){    signalList = signalList_HZ;       sigName = "SMS_TChiHZ_HToGG";
  }else if( doHH ){    signalList = signalList_HH;       sigName = "SMS_TChiHH_HToGG";     
  }else if( doWH ){    signalList = signalList_WH;       sigName = "SMS_TChiWH_HToGG";     }


  vector< MT2Analysis<MT2EstimateTree>* > signals;
  for( unsigned int i=0; i<signalList.size(); i++){
    MT2Analysis<MT2EstimateTree>* sigTemp = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), signalList[i]  );
    signals.push_back( sigTemp );
  }

  vector< MT2Analysis<MT2EstimateTree>* > signals_HH_HH0p25;
  if( doHZ )
    for( unsigned int i=0; i<signalList_HH_HH0p25.size(); i++){
      MT2Analysis<MT2EstimateTree>* sigTemp_HH_HH0p25 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), signalList_HH_HH0p25[i]  );
      signals_HH_HH0p25.push_back( sigTemp_HH_HH0p25 );
  }


  std::string outFileName( Form("%s/ws_%s", outputdir.c_str(), sigName.c_str() ));
  if( doLLbin )    outFileName += "_ll";
  if( doMT2binning )    outFileName += "_mt2";
  if( doSingleLbin )    outFileName += "_onelPT";
  if( doBBbin )    outFileName += "_bb";
  
  TFile* dataFile = new TFile(  Form("%s.root",outFileName.c_str()) , "RECREATE" );
  RooWorkspace ws_sig ("ws_sig");
  RooWorkspace ws_bg  ("ws_bg");
  RooWorkspace ws_data("ws_data");

  std::string dataAmountName( Form("%s/dataAmount",  outputdir.c_str() ));
  if( doLLbin )    dataAmountName += "_ll";
  if( doMT2binning )    dataAmountName += "_mt2";
  if( doSingleLbin )    dataAmountName += "_onelPT";
  if( doBBbin )    dataAmountName += "_bb";
  std::ofstream dataAmount( Form("%s.txt", dataAmountName.c_str() ) );

  std::string isrName( Form("%s/isr",  outputdir.c_str() ));
  if( doLLbin )    isrName += "_ll";
  if( doMT2binning )    isrName += "_mt2";
  if( doSingleLbin )    isrName += "_onelPT";
  if( doBBbin )    isrName += "_bb";
  std::ofstream isr( Form("%s.txt", isrName.c_str()) );

  std::string sigAmountName( Form("%s/sigAmount_%s",  outputdir.c_str(), sigName.c_str() ));
  if( doLLbin )    sigAmountName += "_ll";
  if( doMT2binning )    sigAmountName += "_mt2";
  if( doSingleLbin )    sigAmountName += "_onelPT";
  if( doBBbin )    sigAmountName += "_bb";
  std::ofstream sigAmount( Form("%s.txt", sigAmountName.c_str()) );
        
  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this


  if( doMT2binning )
    std::cout << "doing more than 1 MT2 bins, whooop" << std::endl;

  std::string lowMT2sel  = " && (mt2 < 30 ) ";
  std::string highMT2sel = " && (mt2 >= 30 ) ";
  std::string lowLowPtSel  = " && ( (h_pt/h_mass) <   0.4 ) ";
  std::string highLowPtSel = " && ( (h_pt/h_mass) >=  0.4 ) ";
  std::string lowHighPtSel  = " && ( (h_pt/h_mass) <   1. ) ";
  std::string highHighPtSel = " && ( (h_pt/h_mass) >=  1. ) ";

  // std::string lowMT2sel  = " && (gg_mt2 < 60 ) ";
  // std::string highMT2sel = " && (gg_mt2 >= 60 ) ";

  std::string lowMT2name = "loMT2";
  std::string highMT2name = "hiMT2";
 
  std::string lowPtname = "loPt";
  std::string highPtname = "hiPt";



  RooRealVar h_mass_var("h_mass","Invariant mass", 100,180); 
  RooRealVar weight("weight", "weight", -99,99 );
  //dummy variables to be pushed through
  RooRealVar IntLumi("IntLumi", "IntLumi", 1. );
  RooRealVar dZ("dZ", "dZ", 0. );
   
  RooArgSet vars;
  vars.add(h_mass_var);
  vars.add(weight);
  vars.add(IntLumi);
  vars.add(dZ);


  RooFormulaVar hgg_mass_120( "hgg_mass", "(h_mass - 5.)",  vars );
  RooFormulaVar hgg_mass_125( "hgg_mass", "(h_mass )",      vars );
  RooFormulaVar hgg_mass_130( "hgg_mass", "(h_mass + 5.)",  vars );


  //DILEP REGION
  TTree* tree_data_ll; 
  TTree* tree_higgs_ll; 
  TTree* tree_qcd_ll; 
  TTree* tree_diPhoton_ll;
  TTree* tree_gjets_ll;
  TTree* tree_sig_ll;

  std::string llRegionName = "diLepZ";

  RooDataSet higgs_ds_130_ll ( Form("higgs_130_13TeV_%s", llRegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_130_ll    ( Form("bg_130_13TeV_%s", llRegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_130_ll  ( Form("data_130_13TeV_%s", llRegionName.c_str() ), "Data",   vars, WeightVar(weight)); ;
  
  RooDataSet higgs_ds_125_ll ( Form("higgs_125_13TeV_%s", llRegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_125_ll    ( Form("bg_125_13TeV_%s", llRegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_125_ll  ( Form("data_125_13TeV_%s", llRegionName.c_str() ), "Data",   vars, WeightVar(weight));
  
  RooDataSet higgs_ds_120_ll ( Form("higgs_120_13TeV_%s", llRegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_120_ll    ( Form("bg_120_13TeV_%s", llRegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_120_ll  ( Form("data_120_13TeV_%s", llRegionName.c_str() ), "Data",   vars, WeightVar(weight));

  vector< RooDataSet > sig_ds_130_ll;
  vector< RooDataSet > sig_ds_125_ll;
  vector< RooDataSet > sig_ds_120_ll;
  for( unsigned int i=0; i<signalList.size(); i++){
    RooDataSet sig_ds_130_ll_iter( Form("%s_130_13TeV_%s", signalList[i].c_str(), llRegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_130_ll.push_back( sig_ds_130_ll_iter );
    RooDataSet sig_ds_125_ll_iter( Form("%s_125_13TeV_%s", signalList[i].c_str(), llRegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_125_ll.push_back( sig_ds_125_ll_iter );
    RooDataSet sig_ds_120_ll_iter( Form("%s_120_13TeV_%s", signalList[i].c_str(), llRegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_120_ll.push_back( sig_ds_120_ll_iter );
  }

  TH2D* h2D_isr_ll     = new TH2D( Form("h2D_isr_%s", llRegionName.c_str() ),     "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
  TH2D* h2D_isr_err_ll = new TH2D( Form("h2D_isr_err_%s", llRegionName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );

  vector<TH1D*> h_isr_ll_vec;
  vector<TH1D*> h_isr_UP_ll_vec;
  vector<TH1D*> h_isr_DN_ll_vec;
  for( unsigned int i=0; i<signalList.size(); i++){
    TH1D* h_isr_ll    = new TH1D( Form("h_isr_ll_%s_%s", signalList[i].c_str(), llRegionName.c_str()),    "", 100, 0, 2);
    TH1D* h_isr_UP_ll = new TH1D( Form("h_isr_UP_ll_%s_%s", signalList[i].c_str(), llRegionName.c_str()), "", 100, 0, 2);
    TH1D* h_isr_DN_ll = new TH1D( Form("h_isr_DN_ll_%s_%s", signalList[i].c_str(), llRegionName.c_str()), "", 100, 0, 2);

    h_isr_ll_vec.push_back( h_isr_ll );
    h_isr_UP_ll_vec.push_back( h_isr_UP_ll );
    h_isr_DN_ll_vec.push_back( h_isr_DN_ll );
  }

  //SINGLE ELECTRON REGION
  TTree* tree_data_el;   TTree* tree_higgs_el;   TTree* tree_qcd_el;   TTree* tree_diPhoton_el;  TTree* tree_gjets_el;  TTree* tree_sig_el;
  std::string elRegionName = "is1El";
  //PT0 bin
  std::string el_pT0RegionName = "is1El_pT0";

  RooDataSet higgs_ds_130_el_pT0 ( Form("higgs_130_13TeV_%s", el_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_130_el_pT0    ( Form("bg_130_13TeV_%s", el_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_130_el_pT0  ( Form("data_130_13TeV_%s", el_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight)); ;
  
  RooDataSet higgs_ds_125_el_pT0 ( Form("higgs_125_13TeV_%s", el_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_125_el_pT0    ( Form("bg_125_13TeV_%s", el_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_125_el_pT0  ( Form("data_125_13TeV_%s", el_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight));
  
  RooDataSet higgs_ds_120_el_pT0 ( Form("higgs_120_13TeV_%s", el_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_120_el_pT0    ( Form("bg_120_13TeV_%s", el_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_120_el_pT0  ( Form("data_120_13TeV_%s", el_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight));

  vector< RooDataSet > sig_ds_130_el_pT0;
  vector< RooDataSet > sig_ds_125_el_pT0;
  vector< RooDataSet > sig_ds_120_el_pT0;
  for( unsigned int i=0; i<signalList.size(); i++){
    RooDataSet sig_ds_130_el_pT0_iter( Form("%s_130_13TeV_%s", signalList[i].c_str(), el_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_130_el_pT0.push_back( sig_ds_130_el_pT0_iter );
    RooDataSet sig_ds_125_el_pT0_iter( Form("%s_125_13TeV_%s", signalList[i].c_str(), el_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_125_el_pT0.push_back( sig_ds_125_el_pT0_iter );
    RooDataSet sig_ds_120_el_pT0_iter( Form("%s_120_13TeV_%s", signalList[i].c_str(), el_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_120_el_pT0.push_back( sig_ds_120_el_pT0_iter );
  }

  TH2D* h2D_isr_el_pT0     = new TH2D( Form("h2D_isr_%s", el_pT0RegionName.c_str() ),     "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
  TH2D* h2D_isr_err_el_pT0 = new TH2D( Form("h2D_isr_err_%s", el_pT0RegionName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );

  vector<TH1D*> h_isr_el_pT0_vec;
  vector<TH1D*> h_isr_UP_el_pT0_vec;
  vector<TH1D*> h_isr_DN_el_pT0_vec;
  for( unsigned int i=0; i<signalList.size(); i++){
    TH1D* h_isr_el_pT0    = new TH1D( Form("h_isr_el_pT0_%s_%s", signalList[i].c_str(), el_pT0RegionName.c_str()),    "", 100, 0, 2);
    TH1D* h_isr_UP_el_pT0 = new TH1D( Form("h_isr_UP_el_pT0_%s_%s", signalList[i].c_str(), el_pT0RegionName.c_str()), "", 100, 0, 2);
    TH1D* h_isr_DN_el_pT0 = new TH1D( Form("h_isr_DN_el_pT0_%s_%s", signalList[i].c_str(), el_pT0RegionName.c_str()), "", 100, 0, 2);

    h_isr_el_pT0_vec.push_back( h_isr_el_pT0 );
    h_isr_UP_el_pT0_vec.push_back( h_isr_UP_el_pT0 );
    h_isr_DN_el_pT0_vec.push_back( h_isr_DN_el_pT0 );
  }

  //PT1 bin
  std::string el_pT1RegionName = "is1El_pT1";

  RooDataSet higgs_ds_130_el_pT1 ( Form("higgs_130_13TeV_%s", el_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_130_el_pT1    ( Form("bg_130_13TeV_%s", el_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_130_el_pT1  ( Form("data_130_13TeV_%s", el_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight)); ;
  
  RooDataSet higgs_ds_125_el_pT1 ( Form("higgs_125_13TeV_%s", el_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_125_el_pT1    ( Form("bg_125_13TeV_%s", el_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_125_el_pT1  ( Form("data_125_13TeV_%s", el_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight));
  
  RooDataSet higgs_ds_120_el_pT1 ( Form("higgs_120_13TeV_%s", el_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_120_el_pT1    ( Form("bg_120_13TeV_%s", el_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_120_el_pT1  ( Form("data_120_13TeV_%s", el_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight));

  vector< RooDataSet > sig_ds_130_el_pT1;
  vector< RooDataSet > sig_ds_125_el_pT1;
  vector< RooDataSet > sig_ds_120_el_pT1;
  for( unsigned int i=0; i<signalList.size(); i++){
    RooDataSet sig_ds_130_el_pT1_iter( Form("%s_130_13TeV_%s", signalList[i].c_str(), el_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_130_el_pT1.push_back( sig_ds_130_el_pT1_iter );
    RooDataSet sig_ds_125_el_pT1_iter( Form("%s_125_13TeV_%s", signalList[i].c_str(), el_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_125_el_pT1.push_back( sig_ds_125_el_pT1_iter );
    RooDataSet sig_ds_120_el_pT1_iter( Form("%s_120_13TeV_%s", signalList[i].c_str(), el_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_120_el_pT1.push_back( sig_ds_120_el_pT1_iter );
  }

  TH2D* h2D_isr_el_pT1     = new TH2D( Form("h2D_isr_%s", el_pT1RegionName.c_str() ),     "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
  TH2D* h2D_isr_err_el_pT1 = new TH2D( Form("h2D_isr_err_%s", el_pT1RegionName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );

  vector<TH1D*> h_isr_el_pT1_vec;
  vector<TH1D*> h_isr_UP_el_pT1_vec;
  vector<TH1D*> h_isr_DN_el_pT1_vec;
  for( unsigned int i=0; i<signalList.size(); i++){
    TH1D* h_isr_el_pT1    = new TH1D( Form("h_isr_el_pT1_%s_%s", signalList[i].c_str(), el_pT1RegionName.c_str()),    "", 100, 0, 2);
    TH1D* h_isr_UP_el_pT1 = new TH1D( Form("h_isr_UP_el_pT1_%s_%s", signalList[i].c_str(), el_pT1RegionName.c_str()), "", 100, 0, 2);
    TH1D* h_isr_DN_el_pT1 = new TH1D( Form("h_isr_DN_el_pT1_%s_%s", signalList[i].c_str(), el_pT1RegionName.c_str()), "", 100, 0, 2);

    h_isr_el_pT1_vec.push_back( h_isr_el_pT1 );
    h_isr_UP_el_pT1_vec.push_back( h_isr_UP_el_pT1 );
    h_isr_DN_el_pT1_vec.push_back( h_isr_DN_el_pT1 );
  }
 
  //SINGLE MUON REGION
  TTree* tree_data_mu;  TTree* tree_higgs_mu;   TTree* tree_qcd_mu;   TTree* tree_diPhoton_mu;  TTree* tree_gjets_mu;  TTree* tree_sig_mu;
  std::string muRegionName = "is1Mu";
  //PT0
  std::string mu_pT0RegionName = "is1Mu_pT0";

  RooDataSet higgs_ds_130_mu_pT0 ( Form("higgs_130_13TeV_%s", mu_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_130_mu_pT0    ( Form("bg_130_13TeV_%s", mu_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_130_mu_pT0  ( Form("data_130_13TeV_%s", mu_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight)); ;
  
  RooDataSet higgs_ds_125_mu_pT0 ( Form("higgs_125_13TeV_%s", mu_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_125_mu_pT0    ( Form("bg_125_13TeV_%s", mu_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_125_mu_pT0  ( Form("data_125_13TeV_%s", mu_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight));
  
  RooDataSet higgs_ds_120_mu_pT0 ( Form("higgs_120_13TeV_%s", mu_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_120_mu_pT0    ( Form("bg_120_13TeV_%s", mu_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_120_mu_pT0  ( Form("data_120_13TeV_%s", mu_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight));

  vector< RooDataSet > sig_ds_130_mu_pT0;
  vector< RooDataSet > sig_ds_125_mu_pT0;
  vector< RooDataSet > sig_ds_120_mu_pT0;
  for( unsigned int i=0; i<signalList.size(); i++){
    RooDataSet sig_ds_130_mu_pT0_iter( Form("%s_130_13TeV_%s", signalList[i].c_str(), mu_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_130_mu_pT0.push_back( sig_ds_130_mu_pT0_iter );
    RooDataSet sig_ds_125_mu_pT0_iter( Form("%s_125_13TeV_%s", signalList[i].c_str(), mu_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_125_mu_pT0.push_back( sig_ds_125_mu_pT0_iter );
    RooDataSet sig_ds_120_mu_pT0_iter( Form("%s_120_13TeV_%s", signalList[i].c_str(), mu_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_120_mu_pT0.push_back( sig_ds_120_mu_pT0_iter );
  }

  TH2D* h2D_isr_mu_pT0     = new TH2D( Form("h2D_isr_%s", mu_pT0RegionName.c_str() ),     "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
  TH2D* h2D_isr_err_mu_pT0 = new TH2D( Form("h2D_isr_err_%s", mu_pT0RegionName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );

  vector<TH1D*> h_isr_mu_pT0_vec;
  vector<TH1D*> h_isr_UP_mu_pT0_vec;
  vector<TH1D*> h_isr_DN_mu_pT0_vec;
  for( unsigned int i=0; i<signalList.size(); i++){
    TH1D* h_isr_mu_pT0    = new TH1D( Form("h_isr_mu_pT0_%s_%s", signalList[i].c_str(), mu_pT0RegionName.c_str()),    "", 100, 0, 2);
    TH1D* h_isr_UP_mu_pT0 = new TH1D( Form("h_isr_UP_mu_pT0_%s_%s", signalList[i].c_str(), mu_pT0RegionName.c_str()), "", 100, 0, 2);
    TH1D* h_isr_DN_mu_pT0 = new TH1D( Form("h_isr_DN_mu_pT0_%s_%s", signalList[i].c_str(), mu_pT0RegionName.c_str()), "", 100, 0, 2);

    h_isr_mu_pT0_vec.push_back( h_isr_mu_pT0 );
    h_isr_UP_mu_pT0_vec.push_back( h_isr_UP_mu_pT0 );
    h_isr_DN_mu_pT0_vec.push_back( h_isr_DN_mu_pT0 );
  }

  //PT0
  std::string mu_pT1RegionName = "is1Mu_pT1";

  RooDataSet higgs_ds_130_mu_pT1 ( Form("higgs_130_13TeV_%s", mu_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_130_mu_pT1    ( Form("bg_130_13TeV_%s", mu_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_130_mu_pT1  ( Form("data_130_13TeV_%s", mu_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight)); ;
  
  RooDataSet higgs_ds_125_mu_pT1 ( Form("higgs_125_13TeV_%s", mu_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_125_mu_pT1    ( Form("bg_125_13TeV_%s", mu_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_125_mu_pT1  ( Form("data_125_13TeV_%s", mu_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight));
  
  RooDataSet higgs_ds_120_mu_pT1 ( Form("higgs_120_13TeV_%s", mu_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_120_mu_pT1    ( Form("bg_120_13TeV_%s", mu_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_120_mu_pT1  ( Form("data_120_13TeV_%s", mu_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight));

  vector< RooDataSet > sig_ds_130_mu_pT1;
  vector< RooDataSet > sig_ds_125_mu_pT1;
  vector< RooDataSet > sig_ds_120_mu_pT1;
  for( unsigned int i=0; i<signalList.size(); i++){
    RooDataSet sig_ds_130_mu_pT1_iter( Form("%s_130_13TeV_%s", signalList[i].c_str(), mu_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_130_mu_pT1.push_back( sig_ds_130_mu_pT1_iter );
    RooDataSet sig_ds_125_mu_pT1_iter( Form("%s_125_13TeV_%s", signalList[i].c_str(), mu_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_125_mu_pT1.push_back( sig_ds_125_mu_pT1_iter );
    RooDataSet sig_ds_120_mu_pT1_iter( Form("%s_120_13TeV_%s", signalList[i].c_str(), mu_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_120_mu_pT1.push_back( sig_ds_120_mu_pT1_iter );
  }

  TH2D* h2D_isr_mu_pT1     = new TH2D( Form("h2D_isr_%s", mu_pT1RegionName.c_str() ),     "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
  TH2D* h2D_isr_err_mu_pT1 = new TH2D( Form("h2D_isr_err_%s", mu_pT1RegionName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );

  vector<TH1D*> h_isr_mu_pT1_vec;
  vector<TH1D*> h_isr_UP_mu_pT1_vec;
  vector<TH1D*> h_isr_DN_mu_pT1_vec;
  for( unsigned int i=0; i<signalList.size(); i++){
    TH1D* h_isr_mu_pT1    = new TH1D( Form("h_isr_mu_pT1_%s_%s", signalList[i].c_str(), mu_pT1RegionName.c_str()),    "", 100, 0, 2);
    TH1D* h_isr_UP_mu_pT1 = new TH1D( Form("h_isr_UP_mu_pT1_%s_%s", signalList[i].c_str(), mu_pT1RegionName.c_str()), "", 100, 0, 2);
    TH1D* h_isr_DN_mu_pT1 = new TH1D( Form("h_isr_DN_mu_pT1_%s_%s", signalList[i].c_str(), mu_pT1RegionName.c_str()), "", 100, 0, 2);

    h_isr_mu_pT1_vec.push_back( h_isr_mu_pT1 );
    h_isr_UP_mu_pT1_vec.push_back( h_isr_UP_mu_pT1 );
    h_isr_DN_mu_pT1_vec.push_back( h_isr_DN_mu_pT1 );
  }
 


  //DIB region

  //Z BIN
  TTree* tree_data_bbZ; TTree* tree_higgs_bbZ; TTree* tree_qcd_bbZ; TTree* tree_diPhoton_bbZ; TTree* tree_gjets_bbZ; TTree* tree_sig_bbZ;
  std::string bbZRegionName = "diBBZ";
  //LOW PT BIN
  std::string bbz_pT0RegionName = "diBBZ_pT0";

  RooDataSet higgs_ds_130_bbz_pT0 ( Form("higgs_130_13TeV_%s", bbz_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_130_bbz_pT0    ( Form("bg_130_13TeV_%s", bbz_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_130_bbz_pT0  ( Form("data_130_13TeV_%s", bbz_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight)); ;
  
  RooDataSet higgs_ds_125_bbz_pT0 ( Form("higgs_125_13TeV_%s", bbz_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_125_bbz_pT0    ( Form("bg_125_13TeV_%s", bbz_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_125_bbz_pT0  ( Form("data_125_13TeV_%s", bbz_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight));
  
  RooDataSet higgs_ds_120_bbz_pT0 ( Form("higgs_120_13TeV_%s", bbz_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_120_bbz_pT0    ( Form("bg_120_13TeV_%s", bbz_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_120_bbz_pT0  ( Form("data_120_13TeV_%s", bbz_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight));

  vector< RooDataSet > sig_ds_130_bbz_pT0;
  vector< RooDataSet > sig_ds_125_bbz_pT0;
  vector< RooDataSet > sig_ds_120_bbz_pT0;
  for( unsigned int i=0; i<signalList.size(); i++){
    RooDataSet sig_ds_130_bbz_pT0_iter( Form("%s_130_13TeV_%s", signalList[i].c_str(), bbz_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_130_bbz_pT0.push_back( sig_ds_130_bbz_pT0_iter );
    RooDataSet sig_ds_125_bbz_pT0_iter( Form("%s_125_13TeV_%s", signalList[i].c_str(), bbz_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_125_bbz_pT0.push_back( sig_ds_125_bbz_pT0_iter );
    RooDataSet sig_ds_120_bbz_pT0_iter( Form("%s_120_13TeV_%s", signalList[i].c_str(), bbz_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_120_bbz_pT0.push_back( sig_ds_120_bbz_pT0_iter );
  }

  TH2D* h2D_isr_bbz_pT0     = new TH2D( Form("h2D_isr_%s", bbz_pT0RegionName.c_str() ),     "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
  TH2D* h2D_isr_err_bbz_pT0 = new TH2D( Form("h2D_isr_err_%s", bbz_pT0RegionName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );

  vector<TH1D*> h_isr_bbz_pT0_vec;
  vector<TH1D*> h_isr_UP_bbz_pT0_vec;
  vector<TH1D*> h_isr_DN_bbz_pT0_vec;
  for( unsigned int i=0; i<signalList.size(); i++){
    TH1D* h_isr_bbz_pT0    = new TH1D( Form("h_isr_bbz_pT0_%s_%s", signalList[i].c_str(), bbz_pT0RegionName.c_str()),    "", 100, 0, 2);
    TH1D* h_isr_UP_bbz_pT0 = new TH1D( Form("h_isr_UP_bbz_pT0_%s_%s", signalList[i].c_str(), bbz_pT0RegionName.c_str()), "", 100, 0, 2);
    TH1D* h_isr_DN_bbz_pT0 = new TH1D( Form("h_isr_DN_bbz_pT0_%s_%s", signalList[i].c_str(), bbz_pT0RegionName.c_str()), "", 100, 0, 2);
    h_isr_bbz_pT0_vec.push_back( h_isr_bbz_pT0 );
    h_isr_UP_bbz_pT0_vec.push_back( h_isr_UP_bbz_pT0 );
    h_isr_DN_bbz_pT0_vec.push_back( h_isr_DN_bbz_pT0 );
  }

  //HIGH PT BIN////////////////////////////////////////////////////////
  std::string bbz_pT1RegionName = "diBBZ_pT1";

  RooDataSet higgs_ds_130_bbz_pT1 ( Form("higgs_130_13TeV_%s", bbz_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_130_bbz_pT1    ( Form("bg_130_13TeV_%s", bbz_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_130_bbz_pT1  ( Form("data_130_13TeV_%s", bbz_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight)); ;
  
  RooDataSet higgs_ds_125_bbz_pT1 ( Form("higgs_125_13TeV_%s", bbz_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_125_bbz_pT1    ( Form("bg_125_13TeV_%s", bbz_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_125_bbz_pT1  ( Form("data_125_13TeV_%s", bbz_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight));
  
  RooDataSet higgs_ds_120_bbz_pT1 ( Form("higgs_120_13TeV_%s", bbz_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_120_bbz_pT1    ( Form("bg_120_13TeV_%s", bbz_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_120_bbz_pT1  ( Form("data_120_13TeV_%s", bbz_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight));

  vector< RooDataSet > sig_ds_130_bbz_pT1;
  vector< RooDataSet > sig_ds_125_bbz_pT1;
  vector< RooDataSet > sig_ds_120_bbz_pT1;
  for( unsigned int i=0; i<signalList.size(); i++){
    RooDataSet sig_ds_130_bbz_pT1_iter( Form("%s_130_13TeV_%s", signalList[i].c_str(), bbz_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_130_bbz_pT1.push_back( sig_ds_130_bbz_pT1_iter );
    RooDataSet sig_ds_125_bbz_pT1_iter( Form("%s_125_13TeV_%s", signalList[i].c_str(), bbz_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_125_bbz_pT1.push_back( sig_ds_125_bbz_pT1_iter );
    RooDataSet sig_ds_120_bbz_pT1_iter( Form("%s_120_13TeV_%s", signalList[i].c_str(), bbz_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_120_bbz_pT1.push_back( sig_ds_120_bbz_pT1_iter );
  }

  TH2D* h2D_isr_bbz_pT1     = new TH2D( Form("h2D_isr_%s", bbz_pT1RegionName.c_str() ),     "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
  TH2D* h2D_isr_err_bbz_pT1 = new TH2D( Form("h2D_isr_err_%s", bbz_pT1RegionName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );

  vector<TH1D*> h_isr_bbz_pT1_vec;
  vector<TH1D*> h_isr_UP_bbz_pT1_vec;
  vector<TH1D*> h_isr_DN_bbz_pT1_vec;
  for( unsigned int i=0; i<signalList.size(); i++){
    TH1D* h_isr_bbz_pT1    = new TH1D( Form("h_isr_bbz_pT1_%s_%s", signalList[i].c_str(), bbz_pT1RegionName.c_str()),    "", 100, 0, 2);
    TH1D* h_isr_UP_bbz_pT1 = new TH1D( Form("h_isr_UP_bbz_pT1_%s_%s", signalList[i].c_str(), bbz_pT1RegionName.c_str()), "", 100, 0, 2);
    TH1D* h_isr_DN_bbz_pT1 = new TH1D( Form("h_isr_DN_bbz_pT1_%s_%s", signalList[i].c_str(), bbz_pT1RegionName.c_str()), "", 100, 0, 2);
    h_isr_bbz_pT1_vec.push_back( h_isr_bbz_pT1 );
    h_isr_UP_bbz_pT1_vec.push_back( h_isr_UP_bbz_pT1 );
    h_isr_DN_bbz_pT1_vec.push_back( h_isr_DN_bbz_pT1 );
  }



  //H BIN
  TTree* tree_data_bbH; TTree* tree_higgs_bbH; TTree* tree_qcd_bbH; TTree* tree_diPhoton_bbH; TTree* tree_gjets_bbH; TTree* tree_sig_bbH;
  std::string bbHRegionName = "diBBH";
  //LOW PT BIN/////////////////////////////////////////////////////////////////////////
  std::string bbh_pT0RegionName = "diBBH_pT0";

  RooDataSet higgs_ds_130_bbh_pT0 ( Form("higgs_130_13TeV_%s", bbh_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_130_bbh_pT0    ( Form("bg_130_13TeV_%s", bbh_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_130_bbh_pT0  ( Form("data_130_13TeV_%s", bbh_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight)); ;
  
  RooDataSet higgs_ds_125_bbh_pT0 ( Form("higgs_125_13TeV_%s", bbh_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_125_bbh_pT0    ( Form("bg_125_13TeV_%s", bbh_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_125_bbh_pT0  ( Form("data_125_13TeV_%s", bbh_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight));
  
  RooDataSet higgs_ds_120_bbh_pT0 ( Form("higgs_120_13TeV_%s", bbh_pT0RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_120_bbh_pT0    ( Form("bg_120_13TeV_%s", bbh_pT0RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_120_bbh_pT0  ( Form("data_120_13TeV_%s", bbh_pT0RegionName.c_str() ), "Data",   vars, WeightVar(weight));

  vector< RooDataSet > sig_ds_130_bbh_pT0;
  vector< RooDataSet > sig_ds_125_bbh_pT0;
  vector< RooDataSet > sig_ds_120_bbh_pT0;
  for( unsigned int i=0; i<signalList.size(); i++){
    RooDataSet sig_ds_130_bbh_pT0_iter( Form("%s_130_13TeV_%s", signalList[i].c_str(), bbh_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_130_bbh_pT0.push_back( sig_ds_130_bbh_pT0_iter );
    RooDataSet sig_ds_125_bbh_pT0_iter( Form("%s_125_13TeV_%s", signalList[i].c_str(), bbh_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_125_bbh_pT0.push_back( sig_ds_125_bbh_pT0_iter );
    RooDataSet sig_ds_120_bbh_pT0_iter( Form("%s_120_13TeV_%s", signalList[i].c_str(), bbh_pT0RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_120_bbh_pT0.push_back( sig_ds_120_bbh_pT0_iter );
  }


  TH2D* h2D_isr_bbh_pT0     = new TH2D( Form("h2D_isr_%s", bbh_pT0RegionName.c_str() ),     "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
  TH2D* h2D_isr_err_bbh_pT0 = new TH2D( Form("h2D_isr_err_%s", bbh_pT0RegionName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );

  vector<TH1D*> h_isr_bbh_pT0_vec;
  vector<TH1D*> h_isr_UP_bbh_pT0_vec;
  vector<TH1D*> h_isr_DN_bbh_pT0_vec;
  for( unsigned int i=0; i<signalList.size(); i++){
    TH1D* h_isr_bbh_pT0    = new TH1D( Form("h_isr_bbh_pT0_%s_%s", signalList[i].c_str(), bbh_pT0RegionName.c_str()),    "", 100, 0, 2);
    TH1D* h_isr_UP_bbh_pT0 = new TH1D( Form("h_isr_UP_bbh_pT0_%s_%s", signalList[i].c_str(), bbh_pT0RegionName.c_str()), "", 100, 0, 2);
    TH1D* h_isr_DN_bbh_pT0 = new TH1D( Form("h_isr_DN_bbh_pT0_%s_%s", signalList[i].c_str(), bbh_pT0RegionName.c_str()), "", 100, 0, 2);
    h_isr_bbh_pT0_vec.push_back( h_isr_bbh_pT0 );
    h_isr_UP_bbh_pT0_vec.push_back( h_isr_UP_bbh_pT0 );
    h_isr_DN_bbh_pT0_vec.push_back( h_isr_DN_bbh_pT0 );
  }

  //LOW PT BIN/////////////////////////////////////////////////////////////////////////
  std::string bbh_pT1RegionName = "diBBH_pT1";

  RooDataSet higgs_ds_130_bbh_pT1 ( Form("higgs_130_13TeV_%s", bbh_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_130_bbh_pT1    ( Form("bg_130_13TeV_%s", bbh_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_130_bbh_pT1  ( Form("data_130_13TeV_%s", bbh_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight)); ;
  
  RooDataSet higgs_ds_125_bbh_pT1 ( Form("higgs_125_13TeV_%s", bbh_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_125_bbh_pT1    ( Form("bg_125_13TeV_%s", bbh_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_125_bbh_pT1  ( Form("data_125_13TeV_%s", bbh_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight));
  
  RooDataSet higgs_ds_120_bbh_pT1 ( Form("higgs_120_13TeV_%s", bbh_pT1RegionName.c_str() ), "Higgs", vars, WeightVar(weight));
  RooDataSet bg_ds_120_bbh_pT1    ( Form("bg_120_13TeV_%s", bbh_pT1RegionName.c_str() ), "Background",  vars, WeightVar(weight));
  RooDataSet data_ds_120_bbh_pT1  ( Form("data_120_13TeV_%s", bbh_pT1RegionName.c_str() ), "Data",   vars, WeightVar(weight));

  vector< RooDataSet > sig_ds_130_bbh_pT1;
  vector< RooDataSet > sig_ds_125_bbh_pT1;
  vector< RooDataSet > sig_ds_120_bbh_pT1;
  for( unsigned int i=0; i<signalList.size(); i++){
    RooDataSet sig_ds_130_bbh_pT1_iter( Form("%s_130_13TeV_%s", signalList[i].c_str(), bbh_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_130_bbh_pT1.push_back( sig_ds_130_bbh_pT1_iter );
    RooDataSet sig_ds_125_bbh_pT1_iter( Form("%s_125_13TeV_%s", signalList[i].c_str(), bbh_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_125_bbh_pT1.push_back( sig_ds_125_bbh_pT1_iter );
    RooDataSet sig_ds_120_bbh_pT1_iter( Form("%s_120_13TeV_%s", signalList[i].c_str(), bbh_pT1RegionName.c_str() ), "Signal", vars, WeightVar(weight) );
    sig_ds_120_bbh_pT1.push_back( sig_ds_120_bbh_pT1_iter );
  }


  TH2D* h2D_isr_bbh_pT1     = new TH2D( Form("h2D_isr_%s", bbh_pT1RegionName.c_str() ),     "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
  TH2D* h2D_isr_err_bbh_pT1 = new TH2D( Form("h2D_isr_err_%s", bbh_pT1RegionName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );

  vector<TH1D*> h_isr_bbh_pT1_vec;
  vector<TH1D*> h_isr_UP_bbh_pT1_vec;
  vector<TH1D*> h_isr_DN_bbh_pT1_vec;
  for( unsigned int i=0; i<signalList.size(); i++){
    TH1D* h_isr_bbh_pT1    = new TH1D( Form("h_isr_bbh_pT1_%s_%s", signalList[i].c_str(), bbh_pT1RegionName.c_str()),    "", 100, 0, 2);
    TH1D* h_isr_UP_bbh_pT1 = new TH1D( Form("h_isr_UP_bbh_pT1_%s_%s", signalList[i].c_str(), bbh_pT1RegionName.c_str()), "", 100, 0, 2);
    TH1D* h_isr_DN_bbh_pT1 = new TH1D( Form("h_isr_DN_bbh_pT1_%s_%s", signalList[i].c_str(), bbh_pT1RegionName.c_str()), "", 100, 0, 2);
    h_isr_bbh_pT1_vec.push_back( h_isr_bbh_pT1 );
    h_isr_UP_bbh_pT1_vec.push_back( h_isr_UP_bbh_pT1 );
    h_isr_DN_bbh_pT1_vec.push_back( h_isr_DN_bbh_pT1 );
  }










  ///////////////////////////////////////////////////////////////////////////////////
  ///// BEGIN LOOP OVER REGIONS /////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  std::set<MT2Region> MT2Regions = higgs->getRegions();

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    // std::cout << *iMT2 <<std::endl;

    MT2Region thisRegion( (*iMT2) );

    int mt2BinsIt = 0;
    if( !doMT2binning )
      mt2BinsIt = 1;

    for( ; mt2BinsIt<2; mt2BinsIt++){

      //INITIAL SELECTION 
      string treeSel = "((( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )" ;

      string treeSelWideEta = "((( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=2.4  && fabs(etaGamma1)<=2.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )" ;
      
      string treeSel_ll  = treeSelWideEta + " && (isDiLepZ)";
      if(doLLbin){
	treeSel += " && !(isDiLepZ) ";
      }

      string treeSel_el  = treeSel + " && (is1El)";
      string treeSel_mu  = treeSel + " && (is1Mu)";
      if(doSingleLbin){
	treeSel += " && !(is1El) && !(is1Mu) ";
      }

      string treeSel_bbZ = treeSel + " && isDiBZ";
      string treeSel_bbH = treeSel + " && isDiBH";
      if(doBBbin){
	treeSel += " && (!isDiBH && !isDiBZ) ";
      }

      string regionSaveName =  thisRegion.getName();

      if( doMT2binning && !(regionSaveName=="HT0toInf_j4toInf_b2toInf_pT1")   && !(regionSaveName=="HT0toInf_j4toInf_b2toInf_pT0") && !(regionSaveName=="HT0toInf_j1to3_b2toInf_pT1")   && !(regionSaveName=="HT0toInf_j0_b0toInf_pT0")  && !(regionSaveName=="HT0toInf_j0_b0toInf_pT1") ){
	if( mt2BinsIt==0){
	  treeSel += lowMT2sel;
	  regionSaveName += "_"+lowMT2name;
	}else{
	  treeSel += highMT2sel;
	  regionSaveName += "_"+highMT2name;	
	}

      }else if(doMT2binning && ((regionSaveName=="HT0toInf_j0_b0toInf_pT0")  || (regionSaveName=="HT0toInf_j0_b0toInf_pT1")) ){

	if( (regionSaveName=="HT0toInf_j0_b0toInf_pT0")){
	  if( mt2BinsIt==0){
	    treeSel += lowLowPtSel;
	    regionSaveName += "_"+lowPtname;
	  }else{
	    treeSel += highLowPtSel;
	    regionSaveName += "_"+highPtname;	
	  }

	}else{
	  if( mt2BinsIt==0){
	    treeSel += lowHighPtSel;
	    regionSaveName += "_"+lowPtname;
	  }else{
	    treeSel += highHighPtSel;
	    regionSaveName += "_"+highPtname;	
	  }
	}

      }


      std::string isHighPt = regionSaveName.substr( regionSaveName.find("_pT")+1);
      std::cout << isHighPt << std::endl;
      std::cout << (isHighPt=="pT1") << std::endl;

      if( isHighPt == "pT1_loMT2" || isHighPt == "pT1_loPt" || isHighPt == "pT1_hiMT2" || isHighPt == "pT1_hiPt" )
	isHighPt = "pT1";

      if( isHighPt == "pT0_loMT2" || isHighPt == "pT0_loPt" || isHighPt == "pT0_hiMT2" || isHighPt == "pT0_hiPt" )
	isHighPt = "pT0";


      TH2D* h2D_isr    = new TH2D( Form("h2D_isr_%s", regionSaveName.c_str() ),    "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );
      TH2D* h2D_isr_err = new TH2D(Form("h2D_isr_err_%s", regionSaveName.c_str() ), "", 15, 250-12.5, 612.5, 13 , -12.5, 312.5 );


      TTree* tree_data_ini = data->get(thisRegion)->tree;
      TTree* tree_higgs_ini = higgs->get(thisRegion)->tree;
      TTree* tree_qcd_ini = qcd->get(thisRegion)->tree;
      TTree* tree_diPhoton_ini = diPhoton->get(thisRegion)->tree;
      TTree* tree_gjets_ini = gjets->get(thisRegion)->tree;

      TTree* tree_data = tree_data_ini->CopyTree(treeSel.c_str() );
      TTree* tree_higgs = tree_higgs_ini->CopyTree(treeSel.c_str() );
      TTree* tree_qcd = tree_qcd_ini->CopyTree(treeSel.c_str() );
      TTree* tree_diPhoton = tree_diPhoton_ini->CopyTree(treeSel.c_str() );
      TTree* tree_gjets = tree_gjets_ini->CopyTree(treeSel.c_str() );


      if( doLLbin && mt2BinsIt==1 ){
	tree_data_ll = tree_data_ini->CopyTree(treeSel_ll.c_str() );
	tree_higgs_ll = tree_higgs_ini->CopyTree(treeSel_ll.c_str() );
	tree_qcd_ll = tree_qcd_ini->CopyTree(treeSel_ll.c_str() );
	tree_diPhoton_ll = tree_diPhoton_ini->CopyTree(treeSel_ll.c_str() );
	tree_gjets_ll = tree_gjets_ini->CopyTree(treeSel_ll.c_str() );

	RooDataSet higgs_ds_templl(Form("higgs_temp_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_ll) );
	RooDataSet bg_ds_templl(Form("bg_temp_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton_ll) );
	RooDataSet gjets_ds_templl("gjets_temp_ds", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets_ll) );
	RooDataSet qcd_ds_templl("qcd_temp_ds", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd_ll) );
	RooDataSet data_ds_templl( Form("data_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data_ll) );
	bg_ds_templl.append( qcd_ds_templl );   //append bg data into one
	bg_ds_templl.append( gjets_ds_templl ); //append bg data into one

	RooDataSet higgs_ds_130_templl( higgs_ds_templl, Form("higgs_temp_130_13TeV_%s", llRegionName.c_str() ) );
	RooDataSet bg_ds_130_templl   ( bg_ds_templl,    Form("bg_temp_130_13TeV_%s",    llRegionName.c_str() ) );
	RooDataSet data_ds_130_templl ( data_ds_templl,  Form("data_temp_130_13TeV_%s",  llRegionName.c_str() ) );

	RooDataSet higgs_ds_125_templl( higgs_ds_templl, Form("higgs_temp_125_13TeV_%s", llRegionName.c_str() ) );
	RooDataSet bg_ds_125_templl   ( bg_ds_templl,    Form("bg_temp_125_13TeV_%s",    llRegionName.c_str() ) );
	RooDataSet data_ds_125_templl ( data_ds_templl,  Form("data_temp_125_13TeV_%s",  llRegionName.c_str() ) );

	RooDataSet higgs_ds_120_templl( higgs_ds_templl, Form("higgs_temp_120_13TeV_%s", llRegionName.c_str() ) );
	RooDataSet bg_ds_120_templl   ( bg_ds_templl,    Form("bg_temp_120_13TeV_%s",    llRegionName.c_str() ) );
	RooDataSet data_ds_120_templl ( data_ds_templl,  Form("data_temp_120_13TeV_%s",  llRegionName.c_str() ) );

	higgs_ds_130_ll.append( higgs_ds_130_templl);
	data_ds_130_ll.append( data_ds_130_templl);
	higgs_ds_125_ll.append( higgs_ds_125_templl);
	data_ds_125_ll.append( data_ds_125_templl);
	higgs_ds_120_ll.append( higgs_ds_120_templl);
	data_ds_120_ll.append( data_ds_120_templl);
	bg_ds_130_ll.append(bg_ds_130_templl);
	bg_ds_125_ll.append(bg_ds_125_templl);
	bg_ds_120_ll.append(bg_ds_120_templl);

	std::cout << "Appended a region for the di-Lep region" << std::endl;
	std::cout << regionSaveName << "    " << data_ds_125_templl.sumEntries() << std::endl;
      }
      
      if( doSingleLbin && mt2BinsIt==1 ){
	//SINGLE ELECTRON REGION
	tree_data_el = tree_data_ini->CopyTree(treeSel_el.c_str() );
	tree_higgs_el = tree_higgs_ini->CopyTree(treeSel_el.c_str() );
	tree_qcd_el = tree_qcd_ini->CopyTree(treeSel_el.c_str() );
	tree_diPhoton_el = tree_diPhoton_ini->CopyTree(treeSel_el.c_str() );
	tree_gjets_el = tree_gjets_ini->CopyTree(treeSel_el.c_str() );

	RooDataSet higgs_ds_tempel(Form("higgs_temp_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_el) );
	RooDataSet bg_ds_tempel(Form("bg_temp_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton_el) );
	RooDataSet gjets_ds_tempel("gjets_temp_ds", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets_el) );
	RooDataSet qcd_ds_tempel("qcd_temp_ds", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd_el) );
	RooDataSet data_ds_tempel( Form("data_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data_el) );
	bg_ds_tempel.append( qcd_ds_tempel );   //append bg data into one
	bg_ds_tempel.append( gjets_ds_tempel ); //append bg data into one

	RooDataSet higgs_ds_130_tempel( higgs_ds_tempel, Form("higgs_temp_130_13TeV_%s", elRegionName.c_str() ) );
	RooDataSet bg_ds_130_tempel   ( bg_ds_tempel,    Form("bg_temp_130_13TeV_%s",    elRegionName.c_str() ) );
	RooDataSet data_ds_130_tempel ( data_ds_tempel,  Form("data_temp_130_13TeV_%s",  elRegionName.c_str() ) );

	RooDataSet higgs_ds_125_tempel( higgs_ds_tempel, Form("higgs_temp_125_13TeV_%s", elRegionName.c_str() ) );
	RooDataSet bg_ds_125_tempel   ( bg_ds_tempel,    Form("bg_temp_125_13TeV_%s",    elRegionName.c_str() ) );
	RooDataSet data_ds_125_tempel ( data_ds_tempel,  Form("data_temp_125_13TeV_%s",  elRegionName.c_str() ) );

	RooDataSet higgs_ds_120_tempel( higgs_ds_tempel, Form("higgs_temp_120_13TeV_%s", elRegionName.c_str() ) );
	RooDataSet bg_ds_120_tempel   ( bg_ds_tempel,    Form("bg_temp_120_13TeV_%s",    elRegionName.c_str() ) );
	RooDataSet data_ds_120_tempel ( data_ds_tempel,  Form("data_temp_120_13TeV_%s",  elRegionName.c_str() ) );

	if( isHighPt=="pT0"){
	  higgs_ds_130_el_pT0.append( higgs_ds_130_tempel);
	  data_ds_130_el_pT0.append( data_ds_130_tempel);
	  higgs_ds_125_el_pT0.append( higgs_ds_125_tempel);
	  data_ds_125_el_pT0.append( data_ds_125_tempel);
	  higgs_ds_120_el_pT0.append( higgs_ds_120_tempel);
	  data_ds_120_el_pT0.append( data_ds_120_tempel);
	  bg_ds_130_el_pT0.append(bg_ds_130_tempel);
	  bg_ds_125_el_pT0.append(bg_ds_125_tempel);
	  bg_ds_120_el_pT0.append(bg_ds_120_tempel);
	}else if( isHighPt=="pT1"){
	  higgs_ds_130_el_pT1.append( higgs_ds_130_tempel);
	  data_ds_130_el_pT1.append( data_ds_130_tempel);
	  higgs_ds_125_el_pT1.append( higgs_ds_125_tempel);
	  data_ds_125_el_pT1.append( data_ds_125_tempel);
	  higgs_ds_120_el_pT1.append( higgs_ds_120_tempel);
	  data_ds_120_el_pT1.append( data_ds_120_tempel);
	  bg_ds_130_el_pT1.append(bg_ds_130_tempel);
	  bg_ds_125_el_pT1.append(bg_ds_125_tempel);
	  bg_ds_120_el_pT1.append(bg_ds_120_tempel);
	}	

	std::cout << "Appended a region for the di-Lep region" << std::endl;
	std::cout << regionSaveName << "    " << data_ds_125_tempel.sumEntries() << std::endl;

	//SINGLE MUECTRON REGION
	tree_data_mu = tree_data_ini->CopyTree(treeSel_mu.c_str() );
	tree_higgs_mu = tree_higgs_ini->CopyTree(treeSel_mu.c_str() );
	tree_qcd_mu = tree_qcd_ini->CopyTree(treeSel_mu.c_str() );
	tree_diPhoton_mu = tree_diPhoton_ini->CopyTree(treeSel_mu.c_str() );
	tree_gjets_mu = tree_gjets_ini->CopyTree(treeSel_mu.c_str() );

	RooDataSet higgs_ds_tempmu(Form("higgs_temp_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_mu) );
	RooDataSet bg_ds_tempmu(Form("bg_temp_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton_mu) );
	RooDataSet gjets_ds_tempmu("gjets_temp_ds", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets_mu) );
	RooDataSet qcd_ds_tempmu("qcd_temp_ds", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd_mu) );
	RooDataSet data_ds_tempmu( Form("data_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data_mu) );
	bg_ds_tempmu.append( qcd_ds_tempmu );   //append bg data into one
	bg_ds_tempmu.append( gjets_ds_tempmu ); //append bg data into one

	RooDataSet higgs_ds_130_tempmu( higgs_ds_tempmu, Form("higgs_temp_130_13TeV_%s", muRegionName.c_str() ) );
	RooDataSet bg_ds_130_tempmu   ( bg_ds_tempmu,    Form("bg_temp_130_13TeV_%s",    muRegionName.c_str() ) );
	RooDataSet data_ds_130_tempmu ( data_ds_tempmu,  Form("data_temp_130_13TeV_%s",  muRegionName.c_str() ) );

	RooDataSet higgs_ds_125_tempmu( higgs_ds_tempmu, Form("higgs_temp_125_13TeV_%s", muRegionName.c_str() ) );
	RooDataSet bg_ds_125_tempmu   ( bg_ds_tempmu,    Form("bg_temp_125_13TeV_%s",    muRegionName.c_str() ) );
	RooDataSet data_ds_125_tempmu ( data_ds_tempmu,  Form("data_temp_125_13TeV_%s",  muRegionName.c_str() ) );

	RooDataSet higgs_ds_120_tempmu( higgs_ds_tempmu, Form("higgs_temp_120_13TeV_%s", muRegionName.c_str() ) );
	RooDataSet bg_ds_120_tempmu   ( bg_ds_tempmu,    Form("bg_temp_120_13TeV_%s",    muRegionName.c_str() ) );
	RooDataSet data_ds_120_tempmu ( data_ds_tempmu,  Form("data_temp_120_13TeV_%s",  muRegionName.c_str() ) );

	if( isHighPt=="pT0"){
	  higgs_ds_130_mu_pT0.append( higgs_ds_130_tempmu);
	  data_ds_130_mu_pT0.append( data_ds_130_tempmu);
	  higgs_ds_125_mu_pT0.append( higgs_ds_125_tempmu);
	  data_ds_125_mu_pT0.append( data_ds_125_tempmu);
	  higgs_ds_120_mu_pT0.append( higgs_ds_120_tempmu);
	  data_ds_120_mu_pT0.append( data_ds_120_tempmu);
	  bg_ds_130_mu_pT0.append(bg_ds_130_tempmu);
	  bg_ds_125_mu_pT0.append(bg_ds_125_tempmu);
	  bg_ds_120_mu_pT0.append(bg_ds_120_tempmu);
	}else if( isHighPt=="pT1"){
	  higgs_ds_130_mu_pT1.append( higgs_ds_130_tempmu);
	  data_ds_130_mu_pT1.append( data_ds_130_tempmu);
	  higgs_ds_125_mu_pT1.append( higgs_ds_125_tempmu);
	  data_ds_125_mu_pT1.append( data_ds_125_tempmu);
	  higgs_ds_120_mu_pT1.append( higgs_ds_120_tempmu);
	  data_ds_120_mu_pT1.append( data_ds_120_tempmu);
	  bg_ds_130_mu_pT1.append(bg_ds_130_tempmu);
	  bg_ds_125_mu_pT1.append(bg_ds_125_tempmu);
	  bg_ds_120_mu_pT1.append(bg_ds_120_tempmu);
	}

	std::cout << "Appended a region for the di-Lep region" << std::endl;
	std::cout << regionSaveName << "    " << data_ds_125_tempmu.sumEntries() << std::endl;

      }


      if( doBBbin && mt2BinsIt==1 ){
	/////////////////////BB Z mass window //////////////////////////////
	tree_data_bbZ = tree_data_ini->CopyTree(treeSel_bbZ.c_str() );
	tree_higgs_bbZ = tree_higgs_ini->CopyTree(treeSel_bbZ.c_str() );
	tree_qcd_bbZ = tree_qcd_ini->CopyTree(treeSel_bbZ.c_str() );
	tree_diPhoton_bbZ = tree_diPhoton_ini->CopyTree(treeSel_bbZ.c_str() );
	tree_gjets_bbZ = tree_gjets_ini->CopyTree(treeSel_bbZ.c_str() );

	RooDataSet higgs_ds_tempbbZ(Form("higgs_temp_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_bbZ) );
	RooDataSet bg_ds_tempbbZ(Form("bg_temp_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton_bbZ) );
	RooDataSet gjets_ds_tempbbZ("gjets_temp_ds", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets_bbZ) );
	RooDataSet qcd_ds_tempbbZ("qcd_temp_ds", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd_bbZ) );
	RooDataSet data_ds_tempbbZ( Form("data_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data_bbZ) );
	bg_ds_tempbbZ.append( qcd_ds_tempbbZ );   //append bg data into one
	bg_ds_tempbbZ.append( gjets_ds_tempbbZ ); //append bg data into one

	RooDataSet higgs_ds_130_tempbbZ( higgs_ds_tempbbZ, Form("higgs_temp_130_13TeV_%s", bbZRegionName.c_str() ) );
	RooDataSet bg_ds_130_tempbbZ   ( bg_ds_tempbbZ,    Form("bg_temp_130_13TeV_%s",    bbZRegionName.c_str() ) );
	RooDataSet data_ds_130_tempbbZ ( data_ds_tempbbZ,  Form("data_temp_130_13TeV_%s",  bbZRegionName.c_str() ) );

	RooDataSet higgs_ds_125_tempbbZ( higgs_ds_tempbbZ, Form("higgs_temp_125_13TeV_%s", bbZRegionName.c_str() ) );
	RooDataSet bg_ds_125_tempbbZ   ( bg_ds_tempbbZ,    Form("bg_temp_125_13TeV_%s",    bbZRegionName.c_str() ) );
	RooDataSet data_ds_125_tempbbZ ( data_ds_tempbbZ,  Form("data_temp_125_13TeV_%s",  bbZRegionName.c_str() ) );

	RooDataSet higgs_ds_120_tempbbZ( higgs_ds_tempbbZ, Form("higgs_temp_120_13TeV_%s", bbZRegionName.c_str() ) );
	RooDataSet bg_ds_120_tempbbZ   ( bg_ds_tempbbZ,    Form("bg_temp_120_13TeV_%s",    bbZRegionName.c_str() ) );
	RooDataSet data_ds_120_tempbbZ ( data_ds_tempbbZ,  Form("data_temp_120_13TeV_%s",  bbZRegionName.c_str() ) );

	if( isHighPt=="pT0"){
	  higgs_ds_130_bbz_pT0.append( higgs_ds_130_tempbbZ);
	  data_ds_130_bbz_pT0.append( data_ds_130_tempbbZ);
	  higgs_ds_125_bbz_pT0.append( higgs_ds_125_tempbbZ);
	  data_ds_125_bbz_pT0.append( data_ds_125_tempbbZ);
	  higgs_ds_120_bbz_pT0.append( higgs_ds_120_tempbbZ);
	  data_ds_120_bbz_pT0.append( data_ds_120_tempbbZ);
	  bg_ds_130_bbz_pT0.append(bg_ds_130_tempbbZ);
	  bg_ds_125_bbz_pT0.append(bg_ds_125_tempbbZ);
	  bg_ds_120_bbz_pT0.append(bg_ds_120_tempbbZ);
	}else if( isHighPt=="pT1"){
	  higgs_ds_130_bbz_pT1.append( higgs_ds_130_tempbbZ);
	  data_ds_130_bbz_pT1.append( data_ds_130_tempbbZ);
	  higgs_ds_125_bbz_pT1.append( higgs_ds_125_tempbbZ);
	  data_ds_125_bbz_pT1.append( data_ds_125_tempbbZ);
	  higgs_ds_120_bbz_pT1.append( higgs_ds_120_tempbbZ);
	  data_ds_120_bbz_pT1.append( data_ds_120_tempbbZ);
	  bg_ds_130_bbz_pT1.append(bg_ds_130_tempbbZ);
	  bg_ds_125_bbz_pT1.append(bg_ds_125_tempbbZ);
	  bg_ds_120_bbz_pT1.append(bg_ds_120_tempbbZ);
	}	

	std::cout << "Appended a region for the di-B region" << std::endl;
	std::cout << regionSaveName << "    " << data_ds_125_tempbbZ.sumEntries() << std::endl;


	/////////////////////BB H mass window //////////////////////////////
	tree_data_bbH = tree_data_ini->CopyTree(treeSel_bbH.c_str() );
	tree_higgs_bbH = tree_higgs_ini->CopyTree(treeSel_bbH.c_str() );
	tree_qcd_bbH = tree_qcd_ini->CopyTree(treeSel_bbH.c_str() );
	tree_diPhoton_bbH = tree_diPhoton_ini->CopyTree(treeSel_bbH.c_str() );
	tree_gjets_bbH = tree_gjets_ini->CopyTree(treeSel_bbH.c_str() );

	RooDataSet higgs_ds_tempbbH(Form("higgs_temp_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_bbH) );
	RooDataSet bg_ds_tempbbH(Form("bg_temp_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton_bbH) );
	RooDataSet gjets_ds_tempbbH("gjets_temp_ds", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets_bbH) );
	RooDataSet qcd_ds_tempbbH("qcd_temp_ds", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd_bbH) );
	RooDataSet data_ds_tempbbH( Form("data_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data_bbH) );
	bg_ds_tempbbH.append( qcd_ds_tempbbH );   //append bg data into one
	bg_ds_tempbbH.append( gjets_ds_tempbbH ); //append bg data into one

	RooDataSet higgs_ds_130_tempbbH( higgs_ds_tempbbH, Form("higgs_temp_130_13TeV_%s", bbHRegionName.c_str() ) );
	RooDataSet bg_ds_130_tempbbH   ( bg_ds_tempbbH,    Form("bg_temp_130_13TeV_%s",    bbHRegionName.c_str() ) );
	RooDataSet data_ds_130_tempbbH ( data_ds_tempbbH,  Form("data_temp_130_13TeV_%s",  bbHRegionName.c_str() ) );

	RooDataSet higgs_ds_125_tempbbH( higgs_ds_tempbbH, Form("higgs_temp_125_13TeV_%s", bbHRegionName.c_str() ) );
	RooDataSet bg_ds_125_tempbbH   ( bg_ds_tempbbH,    Form("bg_temp_125_13TeV_%s",    bbHRegionName.c_str() ) );
	RooDataSet data_ds_125_tempbbH ( data_ds_tempbbH,  Form("data_temp_125_13TeV_%s",  bbHRegionName.c_str() ) );

	RooDataSet higgs_ds_120_tempbbH( higgs_ds_tempbbH, Form("higgs_temp_120_13TeV_%s", bbHRegionName.c_str() ) );
	RooDataSet bg_ds_120_tempbbH   ( bg_ds_tempbbH,    Form("bg_temp_120_13TeV_%s",    bbHRegionName.c_str() ) );
	RooDataSet data_ds_120_tempbbH ( data_ds_tempbbH,  Form("data_temp_120_13TeV_%s",  bbHRegionName.c_str() ) );

	if( isHighPt=="pT0"){
	  higgs_ds_130_bbh_pT0.append( higgs_ds_130_tempbbH);
	  data_ds_130_bbh_pT0.append( data_ds_130_tempbbH);
	  higgs_ds_125_bbh_pT0.append( higgs_ds_125_tempbbH);
	  data_ds_125_bbh_pT0.append( data_ds_125_tempbbH);
	  higgs_ds_120_bbh_pT0.append( higgs_ds_120_tempbbH);
	  data_ds_120_bbh_pT0.append( data_ds_120_tempbbH);
	  bg_ds_130_bbh_pT0.append(bg_ds_130_tempbbH);
	  bg_ds_125_bbh_pT0.append(bg_ds_125_tempbbH);
	  bg_ds_120_bbh_pT0.append(bg_ds_120_tempbbH);
	}else if( isHighPt=="pT1"){
	  higgs_ds_130_bbh_pT1.append( higgs_ds_130_tempbbH);
	  data_ds_130_bbh_pT1.append( data_ds_130_tempbbH);
	  higgs_ds_125_bbh_pT1.append( higgs_ds_125_tempbbH);
	  data_ds_125_bbh_pT1.append( data_ds_125_tempbbH);
	  higgs_ds_120_bbh_pT1.append( higgs_ds_120_tempbbH);
	  data_ds_120_bbh_pT1.append( data_ds_120_tempbbH);
	  bg_ds_130_bbh_pT1.append(bg_ds_130_tempbbH);
	  bg_ds_125_bbh_pT1.append(bg_ds_125_tempbbH);
	  bg_ds_120_bbh_pT1.append(bg_ds_120_tempbbH);
	}

	std::cout << "Appended a region for the di-B region" << std::endl;
	std::cout << regionSaveName << "    " << data_ds_125_tempbbH.sumEntries() << std::endl;

      }
   
      //Import
      RooDataSet higgs_ds(Form("higgs_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs) );
      RooDataSet bg_ds(Form("bg_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton) ) ;
      RooDataSet gjets_ds("gjets_ds", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets) ) ;
      RooDataSet qcd_ds("qcd_ds", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd) ) ;

      RooDataSet data_ds( Form("data_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data) ) ;
      bg_ds.append( qcd_ds );   //append bg data into one
      bg_ds.append( gjets_ds ); //append bg data into one

      RooDataSet higgs_ds_130(  higgs_ds, Form("higgs_130_13TeV_%s", regionSaveName.c_str() ) );
      RooDataSet bg_ds_130   (  bg_ds,    Form("bg_130_13TeV_%s",    regionSaveName.c_str() ) );
      RooDataSet data_ds_130 (  data_ds,  Form("data_130_13TeV_%s",  regionSaveName.c_str() ) );

      RooDataSet higgs_ds_125(  higgs_ds, Form("higgs_125_13TeV_%s", regionSaveName.c_str() ) );
      RooDataSet bg_ds_125   (  bg_ds,    Form("bg_125_13TeV_%s",    regionSaveName.c_str() ) );
      RooDataSet data_ds_125 (  data_ds,  Form("data_125_13TeV_%s",  regionSaveName.c_str() ) );

      RooDataSet higgs_ds_120(  higgs_ds, Form("higgs_120_13TeV_%s", regionSaveName.c_str() ) );
      RooDataSet bg_ds_120   (  bg_ds,    Form("bg_120_13TeV_%s",    regionSaveName.c_str() ) );
      RooDataSet data_ds_120 (  data_ds,  Form("data_120_13TeV_%s",  regionSaveName.c_str() ) );

      ((RooRealVar*)higgs_ds_130.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
      ((RooRealVar*)data_ds_130.addColumn( hgg_mass_130 ))->setRange(100.,180.); 

      ((RooRealVar*)higgs_ds_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)data_ds_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 

      ((RooRealVar*)higgs_ds_120.addColumn( hgg_mass_120 ))->setRange(100.,180.);
      ((RooRealVar*)data_ds_120.addColumn( hgg_mass_120 ))->setRange(100.,180.);

      ws_data.import( data_ds_125 );
      ws_bg.import( bg_ds_125   );
      ws_sig.import( higgs_ds_125  );

      ws_data.import( data_ds_130 );
      ws_bg.import( bg_ds_130   );
      ws_sig.import( higgs_ds_130  );

      ws_data.import( data_ds_120 );
      ws_bg.import( bg_ds_120   );
      ws_sig.import( higgs_ds_120  );

      //write number of data events per region into txt file
      dataAmount << regionSaveName << "    " << data_ds_125.sumEntries() << std::endl;






      //NOW for the SIGNAL LIST
      for( unsigned int i=0; i<signalList.size(); i++){

	std::cout << "Working on signal nr " << i << std::endl;

	TTree* tree_sig_ini = signals[i]->get(thisRegion)->tree;  
	TTree* tree_sig = tree_sig_ini->CopyTree(treeSel.c_str() );
	RooDataSet sig_ds(Form("%s_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig) ) ;

	if(doHZ){
	  std::cout << "ADDING HH to HZ" << std::endl;

	  TTree* tree_sig_ini_HH_HH0p25 = signals_HH_HH0p25[i]->get(thisRegion)->tree;  
	  TTree* tree_sig_HH_HH0p25 = tree_sig_ini_HH_HH0p25->CopyTree(treeSel.c_str() );
	  RooDataSet sig_ds_HH_HH0p25(Form("%s_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25) ) ;

	  sig_ds.append( sig_ds_HH_HH0p25 );
	}

	if( doISRsyst ){
	  TH1D* h_isr    = new TH1D("h_isr",    "", 100, 0, 2);
	  TH1D* h_isr_UP = new TH1D("h_isr_UP", "", 100, 0, 2);
	  TH1D* h_isr_DN = new TH1D("h_isr_DN", "", 100, 0, 2);
 
	  tree_sig->Draw("weight_isr    >> h_isr");
	  tree_sig->Draw("weight_isr_UP >> h_isr_UP");
	  tree_sig->Draw("weight_isr_DN >> h_isr_DN");

	  float sign = 1;
	  float isr_uncert = max (fabs(h_isr_UP->GetMean() - h_isr->GetMean()), fabs( h_isr->GetMean() - h_isr_DN->GetMean() ) );

	  if( fabs(h_isr_UP->GetMean() - h_isr->GetMean()) >  fabs(  h_isr->GetMean() - h_isr_DN->GetMean() ))
	    sign = ((h_isr_UP->GetMean() - h_isr->GetMean())  > 0 ) - ((h_isr_UP->GetMean() - h_isr->GetMean())  < 0) ;
	  else
	    sign = ((h_isr->GetMean() - h_isr_DN->GetMean())  > 0 ) - ((h_isr->GetMean() - h_isr_DN->GetMean())  < 0) ;
	  
	  isr <<  signalList[i] << " \t " << regionSaveName << " \t " <<    sign* isr_uncert  <<std::endl;

	  int massSbottom = 0;
	  int massChi = 0;

	  std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	  std::string::size_type sz;   // alias of size_t
	  massChi = std::stoi (mLSP,&sz);
	  std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	  massSbottom = std::stoi (mSbottom,&sz);

	  int binBottom = h2D_isr->GetXaxis()->FindBin( (float)massSbottom );
	  int binChi    = h2D_isr->GetYaxis()->FindBin( (float)massChi );

	  h2D_isr->SetBinContent( binBottom, binChi, h_isr->GetMean() );
	  h2D_isr_err->SetBinContent( binBottom, binChi, sign* isr_uncert );

	  delete h_isr; delete h_isr_UP; delete h_isr_DN;

	  if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	    sigAmount <<  signalList[i] << " \t " << regionSaveName << " \t " <<  sig_ds.sumEntries()  << std::endl;
	
	} 

	if( i==0 )
	  sigAmount <<  signalList[i] << " \t " << regionSaveName << " \t " <<  sig_ds.sumEntries()  << std::endl;

	RooDataSet sig_ds_130(  sig_ds, Form("%s_130_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ) );
	RooDataSet sig_ds_125(  sig_ds, Form("%s_125_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ) );
	RooDataSet sig_ds_120(  sig_ds, Form("%s_120_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ) );

	((RooRealVar*)sig_ds_130.addColumn( hgg_mass_130 ))->setRange(100.,180.);  
	((RooRealVar*)sig_ds_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
	((RooRealVar*)sig_ds_120.addColumn( hgg_mass_120 ))->setRange(100.,180.); 

	ws_sig.import( sig_ds_125  );
	ws_sig.import( sig_ds_130  );
	ws_sig.import( sig_ds_120  );

 	///////////////////////////////////////////
	////// DI-LEPTON BIN //////////////////////
	///////////////////////////////////////////
	if( doLLbin && mt2BinsIt==1 ){

	  tree_sig_ll = tree_sig_ini->CopyTree(treeSel_ll.c_str() );

	  RooDataSet sig_ds_templl( Form("%s_templl_125_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ), "Sig", vars, WeightVar(weight), Import(*tree_sig_ll) );

	  if(doHZ){
	    std::cout << "ADDING HH to HZ for ll" << std::endl;

	    TTree* tree_sig_ini_HH_HH0p25_ll = signals_HH_HH0p25[i]->get(thisRegion)->tree;  
	    TTree* tree_sig_HH_HH0p25_ll = tree_sig_ini_HH_HH0p25_ll->CopyTree(treeSel_ll.c_str() );
	    RooDataSet sig_ds_HH_HH0p25_ll(Form("%s_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_ll) ) ;

	    sig_ds_templl.append( sig_ds_HH_HH0p25_ll );
	  }

	  RooDataSet sig_ds_130_templl( sig_ds_templl, Form("%s_temp_130_13TeV_%s", signalList[i].c_str(), llRegionName.c_str() ) );
	  RooDataSet sig_ds_125_templl( sig_ds_templl, Form("%s_temp_125_13TeV_%s", signalList[i].c_str(), llRegionName.c_str() ) );
	  RooDataSet sig_ds_120_templl( sig_ds_templl, Form("%s_temp_120_13TeV_%s", signalList[i].c_str(), llRegionName.c_str() ) );

	  sig_ds_130_ll[i].append( sig_ds_130_templl );
	  sig_ds_125_ll[i].append( sig_ds_125_templl );
	  sig_ds_120_ll[i].append( sig_ds_120_templl );

	  if( doISRsyst ){
	    TH1D* h_isr    = new TH1D("h_isr",    "", 100, 0, 2);
	    TH1D* h_isr_UP = new TH1D("h_isr_UP", "", 100, 0, 2);
	    TH1D* h_isr_DN = new TH1D("h_isr_DN", "", 100, 0, 2);
 
	    tree_sig_ll->Draw("weight_isr    >> h_isr");
	    tree_sig_ll->Draw("weight_isr_UP >> h_isr_UP");
	    tree_sig_ll->Draw("weight_isr_DN >> h_isr_DN");

	    h_isr_ll_vec[i]->Add( h_isr );
	    h_isr_UP_ll_vec[i]->Add( h_isr_UP );
	    h_isr_DN_ll_vec[i]->Add( h_isr_DN );

	    delete h_isr; delete h_isr_UP; delete h_isr_DN;
	  } 
	}

	if( doSingleLbin && mt2BinsIt==1 ){

	  //Single Electron region
	  tree_sig_el = tree_sig_ini->CopyTree(treeSel_el.c_str() );

	  RooDataSet sig_ds_tempel( Form("%s_tempel_125_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ), "Sig", vars, WeightVar(weight), Import(*tree_sig_el) );

	  if(doHZ){
	    std::cout << "ADDING HH to HZ for el" << std::endl;

	    TTree* tree_sig_ini_HH_HH0p25_el = signals_HH_HH0p25[i]->get(thisRegion)->tree;  
	    TTree* tree_sig_HH_HH0p25_el = tree_sig_ini_HH_HH0p25_el->CopyTree(treeSel_el.c_str() );
	    RooDataSet sig_ds_HH_HH0p25_el(Form("%s_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_el) ) ;

	    sig_ds_tempel.append( sig_ds_HH_HH0p25_el );
	  }

	  RooDataSet sig_ds_130_tempel( sig_ds_tempel, Form("%s_temp_130_13TeV_%s", signalList[i].c_str(), elRegionName.c_str() ) );
	  RooDataSet sig_ds_125_tempel( sig_ds_tempel, Form("%s_temp_125_13TeV_%s", signalList[i].c_str(), elRegionName.c_str() ) );
	  RooDataSet sig_ds_120_tempel( sig_ds_tempel, Form("%s_temp_120_13TeV_%s", signalList[i].c_str(), elRegionName.c_str() ) );

	  if( isHighPt=="pT0"){
	    sig_ds_130_el_pT0[i].append( sig_ds_130_tempel );
	    sig_ds_125_el_pT0[i].append( sig_ds_125_tempel );
	    sig_ds_120_el_pT0[i].append( sig_ds_120_tempel );
	  }else if( isHighPt=="pT1"){
	    sig_ds_130_el_pT1[i].append( sig_ds_130_tempel );
	    sig_ds_125_el_pT1[i].append( sig_ds_125_tempel );
	    sig_ds_120_el_pT1[i].append( sig_ds_120_tempel );
	  }	  

	  if( doISRsyst ){
	    TH1D* h_isr    = new TH1D("h_isr",    "", 100, 0, 2);
	    TH1D* h_isr_UP = new TH1D("h_isr_UP", "", 100, 0, 2);
	    TH1D* h_isr_DN = new TH1D("h_isr_DN", "", 100, 0, 2);
 
	    tree_sig_el->Draw("weight_isr    >> h_isr");
	    tree_sig_el->Draw("weight_isr_UP >> h_isr_UP");
	    tree_sig_el->Draw("weight_isr_DN >> h_isr_DN");

	    if( isHighPt=="pT0"){
	      h_isr_el_pT0_vec[i]->Add( h_isr );
	      h_isr_UP_el_pT0_vec[i]->Add( h_isr_UP );
	      h_isr_DN_el_pT0_vec[i]->Add( h_isr_DN );
	    }else if( isHighPt=="pT1"){
	      h_isr_el_pT1_vec[i]->Add( h_isr );
	      h_isr_UP_el_pT1_vec[i]->Add( h_isr_UP );
	      h_isr_DN_el_pT1_vec[i]->Add( h_isr_DN );
	    }

	    delete h_isr; delete h_isr_UP; delete h_isr_DN;
	  } 


	  //Single Muectron region
	  tree_sig_mu = tree_sig_ini->CopyTree(treeSel_mu.c_str() );

	  RooDataSet sig_ds_tempmu( Form("%s_tempmu_125_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ), "Sig", vars, WeightVar(weight), Import(*tree_sig_mu) );

	  if(doHZ){
	    std::cout << "ADDING HH to HZ for mu" << std::endl;

	    TTree* tree_sig_ini_HH_HH0p25_mu = signals_HH_HH0p25[i]->get(thisRegion)->tree;  
	    TTree* tree_sig_HH_HH0p25_mu = tree_sig_ini_HH_HH0p25_mu->CopyTree(treeSel_mu.c_str() );
	    RooDataSet sig_ds_HH_HH0p25_mu(Form("%s_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_mu) ) ;

	    sig_ds_tempmu.append( sig_ds_HH_HH0p25_mu );
	  }

	  RooDataSet sig_ds_130_tempmu( sig_ds_tempmu, Form("%s_temp_130_13TeV_%s", signalList[i].c_str(), muRegionName.c_str() ) );
	  RooDataSet sig_ds_125_tempmu( sig_ds_tempmu, Form("%s_temp_125_13TeV_%s", signalList[i].c_str(), muRegionName.c_str() ) );
	  RooDataSet sig_ds_120_tempmu( sig_ds_tempmu, Form("%s_temp_120_13TeV_%s", signalList[i].c_str(), muRegionName.c_str() ) );

	  if( isHighPt=="pT0"){
	    sig_ds_130_mu_pT0[i].append( sig_ds_130_tempmu );
	    sig_ds_125_mu_pT0[i].append( sig_ds_125_tempmu );
	    sig_ds_120_mu_pT0[i].append( sig_ds_120_tempmu );
	  }else if( isHighPt=="pT1"){
	    sig_ds_130_mu_pT1[i].append( sig_ds_130_tempmu );
	    sig_ds_125_mu_pT1[i].append( sig_ds_125_tempmu );
	    sig_ds_120_mu_pT1[i].append( sig_ds_120_tempmu );
	  }

	  if( doISRsyst ){
	    TH1D* h_isr    = new TH1D("h_isr",    "", 100, 0, 2);
	    TH1D* h_isr_UP = new TH1D("h_isr_UP", "", 100, 0, 2);
	    TH1D* h_isr_DN = new TH1D("h_isr_DN", "", 100, 0, 2);
 
	    tree_sig_mu->Draw("weight_isr    >> h_isr");
	    tree_sig_mu->Draw("weight_isr_UP >> h_isr_UP");
	    tree_sig_mu->Draw("weight_isr_DN >> h_isr_DN");

	    if( isHighPt=="pT0"){
	      h_isr_mu_pT0_vec[i]->Add( h_isr );
	      h_isr_UP_mu_pT0_vec[i]->Add( h_isr_UP );
	      h_isr_DN_mu_pT0_vec[i]->Add( h_isr_DN );
	    }else if( isHighPt=="pT1"){
	      h_isr_mu_pT1_vec[i]->Add( h_isr );
	      h_isr_UP_mu_pT1_vec[i]->Add( h_isr_UP );
	      h_isr_DN_mu_pT1_vec[i]->Add( h_isr_DN );
	    }

	    delete h_isr; delete h_isr_UP; delete h_isr_DN;
	  } 



	}

	///////////////////////////////////////////
	//////   DI-B BIN    //////////////////////
	///////////////////////////////////////////
	if( doBBbin && mt2BinsIt==1 ){

	  /////////////////////BB Z mass window //////////////////////////////
	  tree_sig_bbZ = tree_sig_ini->CopyTree(treeSel_bbZ.c_str() );

	  RooDataSet sig_ds_tempbbZ( Form("%s_tempBBZ_125_13TeV_%s", signalList[i].c_str(), bbZRegionName.c_str() ), "Sig", vars, WeightVar(weight), Import(*tree_sig_bbZ) );

	  if(doHZ){
	    std::cout << "ADDING HH to HZ for bbZ" << std::endl;

	    TTree* tree_sig_ini_HH_HH0p25_bbZ = signals_HH_HH0p25[i]->get(thisRegion)->tree;  
	    TTree* tree_sig_HH_HH0p25_bbZ = tree_sig_ini_HH_HH0p25_bbZ->CopyTree(treeSel_bbZ.c_str() );
	    RooDataSet sig_ds_HH_HH0p25_bbZ(Form("%s_BBZ_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(), bbHRegionName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_bbZ) ) ;

	    sig_ds_tempbbZ.append( sig_ds_HH_HH0p25_bbZ );
	  }

	  RooDataSet sig_ds_130_tempbbZ( sig_ds_tempbbZ, Form("%s_tempBBZ_130_13TeV_%s", signalList[i].c_str(), bbZRegionName.c_str() ) );
	  RooDataSet sig_ds_125_tempbbZ( sig_ds_tempbbZ, Form("%s_tempBBZ_125_13TeV_%s", signalList[i].c_str(), bbZRegionName.c_str() ) );
	  RooDataSet sig_ds_120_tempbbZ( sig_ds_tempbbZ, Form("%s_tempBBZ_120_13TeV_%s", signalList[i].c_str(), bbZRegionName.c_str() ) );

	  if( isHighPt=="pT0"){
	    sig_ds_130_bbz_pT0[i].append( sig_ds_130_tempbbZ );
	    sig_ds_125_bbz_pT0[i].append( sig_ds_125_tempbbZ );
	    sig_ds_120_bbz_pT0[i].append( sig_ds_120_tempbbZ );
	  }else if( isHighPt=="pT1"){
	    sig_ds_130_bbz_pT1[i].append( sig_ds_130_tempbbZ );
	    sig_ds_125_bbz_pT1[i].append( sig_ds_125_tempbbZ );
	    sig_ds_120_bbz_pT1[i].append( sig_ds_120_tempbbZ );
	  }	  


	  if( doISRsyst ){
	    TH1D* h_isr    = new TH1D("h_isr",    "", 100, 0, 2);
	    TH1D* h_isr_UP = new TH1D("h_isr_UP", "", 100, 0, 2);
	    TH1D* h_isr_DN = new TH1D("h_isr_DN", "", 100, 0, 2);
 
	    tree_sig_bbZ->Draw("weight_isr    >> h_isr");
	    tree_sig_bbZ->Draw("weight_isr_UP >> h_isr_UP");
	    tree_sig_bbZ->Draw("weight_isr_DN >> h_isr_DN");

	    if( isHighPt=="pT0"){
	      h_isr_bbz_pT0_vec[i]->Add( h_isr );
	      h_isr_UP_bbz_pT0_vec[i]->Add( h_isr_UP );
	      h_isr_DN_bbz_pT0_vec[i]->Add( h_isr_DN );
	    }else if( isHighPt=="pT1"){
	      h_isr_bbz_pT1_vec[i]->Add( h_isr );
	      h_isr_UP_bbz_pT1_vec[i]->Add( h_isr_UP );
	      h_isr_DN_bbz_pT1_vec[i]->Add( h_isr_DN );
	    }

	    delete h_isr; delete h_isr_UP; delete h_isr_DN;
	  } 

	  /////////////////////BB H mass window //////////////////////////////
	  tree_sig_bbH = tree_sig_ini->CopyTree(treeSel_bbH.c_str() );

	  RooDataSet sig_ds_tempbbH( Form("%s_tempBBH_125_13TeV_%s", signalList[i].c_str(), bbHRegionName.c_str() ), "Sig", vars, WeightVar(weight), Import(*tree_sig_bbH) );

	  if(doHZ){
	    std::cout << "ADDING HH to HZ for bbH" << std::endl;

	    TTree* tree_sig_ini_HH_HH0p25_bbH = signals_HH_HH0p25[i]->get(thisRegion)->tree;  
	    TTree* tree_sig_HH_HH0p25_bbH = tree_sig_ini_HH_HH0p25_bbH->CopyTree(treeSel_bbH.c_str() );
	    RooDataSet sig_ds_HH_HH0p25_bbH(Form("%s_BBH_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(), bbHRegionName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_bbH) ) ;

	    sig_ds_tempbbH.append( sig_ds_HH_HH0p25_bbH );
	    std::cout << "ADDED HH to HZ for bbH" << std::endl;
	  }

	  RooDataSet sig_ds_130_tempbbH( sig_ds_tempbbH, Form("%s_tempBBH_130_13TeV_%s", signalList[i].c_str(), bbHRegionName.c_str() ) );
	  RooDataSet sig_ds_125_tempbbH( sig_ds_tempbbH, Form("%s_tempBBH_125_13TeV_%s", signalList[i].c_str(), bbHRegionName.c_str() ) );
	  RooDataSet sig_ds_120_tempbbH( sig_ds_tempbbH, Form("%s_tempBBH_120_13TeV_%s", signalList[i].c_str(), bbHRegionName.c_str() ) );

	  if( isHighPt=="pT0"){
	    sig_ds_130_bbh_pT0[i].append( sig_ds_130_tempbbH );
	    sig_ds_125_bbh_pT0[i].append( sig_ds_125_tempbbH );
	    sig_ds_120_bbh_pT0[i].append( sig_ds_120_tempbbH );
	  }else if( isHighPt=="pT1"){
	    sig_ds_130_bbh_pT1[i].append( sig_ds_130_tempbbH );
	    sig_ds_125_bbh_pT1[i].append( sig_ds_125_tempbbH );
	    sig_ds_120_bbh_pT1[i].append( sig_ds_120_tempbbH );
	  }

	  if( doISRsyst ){
	    TH1D* h_isr    = new TH1D("h_isr",    "", 100, 0, 2);
	    TH1D* h_isr_UP = new TH1D("h_isr_UP", "", 100, 0, 2);
	    TH1D* h_isr_DN = new TH1D("h_isr_DN", "", 100, 0, 2);
 
	    tree_sig_bbH->Draw("weight_isr    >> h_isr");
	    tree_sig_bbH->Draw("weight_isr_UP >> h_isr_UP");
	    tree_sig_bbH->Draw("weight_isr_DN >> h_isr_DN");

	    if( isHighPt=="pT0"){
	      h_isr_bbh_pT0_vec[i]->Add( h_isr );
	      h_isr_UP_bbh_pT0_vec[i]->Add( h_isr_UP );
	      h_isr_DN_bbh_pT0_vec[i]->Add( h_isr_DN );
	    }else if( isHighPt=="pT1"){
	      h_isr_bbh_pT1_vec[i]->Add( h_isr );
	      h_isr_UP_bbh_pT1_vec[i]->Add( h_isr_UP );
	      h_isr_DN_bbh_pT1_vec[i]->Add( h_isr_DN );
	    }

	    delete h_isr; delete h_isr_UP; delete h_isr_DN;
	  } 

	}//done bb bin

      }//done signals for 1 region
   
      h2D_isr->Write();
      h2D_isr_err->Write();

      delete h2D_isr; delete h2D_isr_err;

      // // // Print unbinned dataset with default frame binning (100 bins)
      // //RooPlot* frame3 = CMS_hgg_mass.frame(Title("H mass")) ;
      // RooPlot* frame3 = h_mass_var.frame(Title("H mass")) ;
  
      // data_ds.plotOn(frame3, MarkerColor(kRed) ) ;
      // higgs_ds_125.plotOn(frame3, MarkerColor(kGreen) ) ;
      // sig_ds_125.plotOn(frame3) ;
      // bg_ds.plotOn(frame3, MarkerColor(kBlue) ) ;

      // TCanvas* c = new TCanvas("rf102_dataimport","rf102_dataimport",800,800) ;
      // frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;

      // c->SaveAs(  Form("%s/CMS_hgg_mass_%s.png", outputdir.c_str(), regionSaveName.c_str() ) );

    } // loop over mt2 bin(s)

  }




  if(doLLbin){
    ((RooRealVar*)higgs_ds_130_ll.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130_ll .addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_125_ll.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125_ll .addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_120_ll.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120_ll .addColumn( hgg_mass_120 ))->setRange(100.,180.);

    ws_sig .import( higgs_ds_130_ll );
    ws_bg  .import( bg_ds_130_ll    );
    ws_data.import( data_ds_130_ll  );
  
    ws_sig .import( higgs_ds_125_ll );
    ws_bg  .import( bg_ds_125_ll    );
    ws_data.import( data_ds_125_ll  );
  
    ws_sig .import( higgs_ds_120_ll );
    ws_bg  .import( bg_ds_120_ll    );
    ws_data.import( data_ds_120_ll  );

    dataAmount << llRegionName << "    " << data_ds_125_ll.sumEntries() << std::endl;
  }

  if(doSingleLbin){
    ((RooRealVar*)higgs_ds_130_el_pT0.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130_el_pT0 .addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_125_el_pT0.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125_el_pT0 .addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_120_el_pT0.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120_el_pT0 .addColumn( hgg_mass_120 ))->setRange(100.,180.);

    ws_sig .import( higgs_ds_130_el_pT0 );
    ws_bg  .import( bg_ds_130_el_pT0    );
    ws_data.import( data_ds_130_el_pT0  );
  
    ws_sig .import( higgs_ds_125_el_pT0 );
    ws_bg  .import( bg_ds_125_el_pT0    );
    ws_data.import( data_ds_125_el_pT0  );
  
    ws_sig .import( higgs_ds_120_el_pT0 );
    ws_bg  .import( bg_ds_120_el_pT0    );
    ws_data.import( data_ds_120_el_pT0  );

    dataAmount << el_pT0RegionName << "    " << data_ds_125_el_pT0.sumEntries() << std::endl;

    ((RooRealVar*)higgs_ds_130_el_pT1.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130_el_pT1 .addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_125_el_pT1.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125_el_pT1 .addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_120_el_pT1.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120_el_pT1 .addColumn( hgg_mass_120 ))->setRange(100.,180.);

    ws_sig .import( higgs_ds_130_el_pT1 );
    ws_bg  .import( bg_ds_130_el_pT1    );
    ws_data.import( data_ds_130_el_pT1  );
  
    ws_sig .import( higgs_ds_125_el_pT1 );
    ws_bg  .import( bg_ds_125_el_pT1    );
    ws_data.import( data_ds_125_el_pT1  );
  
    ws_sig .import( higgs_ds_120_el_pT1 );
    ws_bg  .import( bg_ds_120_el_pT1    );
    ws_data.import( data_ds_120_el_pT1  );

    dataAmount << el_pT1RegionName << "    " << data_ds_125_el_pT1.sumEntries() << std::endl;

    ((RooRealVar*)higgs_ds_130_mu_pT0.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130_mu_pT0 .addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_125_mu_pT0.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125_mu_pT0 .addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_120_mu_pT0.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120_mu_pT0 .addColumn( hgg_mass_120 ))->setRange(100.,180.);

    ws_sig .import( higgs_ds_130_mu_pT0 );
    ws_bg  .import( bg_ds_130_mu_pT0    );
    ws_data.import( data_ds_130_mu_pT0  );
  
    ws_sig .import( higgs_ds_125_mu_pT0 );
    ws_bg  .import( bg_ds_125_mu_pT0    );
    ws_data.import( data_ds_125_mu_pT0  );
  
    ws_sig .import( higgs_ds_120_mu_pT0 );
    ws_bg  .import( bg_ds_120_mu_pT0    );
    ws_data.import( data_ds_120_mu_pT0  );

    dataAmount << mu_pT0RegionName << "    " << data_ds_125_mu_pT0.sumEntries() << std::endl;

    ((RooRealVar*)higgs_ds_130_mu_pT1.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130_mu_pT1 .addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_125_mu_pT1.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125_mu_pT1 .addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_120_mu_pT1.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120_mu_pT1 .addColumn( hgg_mass_120 ))->setRange(100.,180.);

    ws_sig .import( higgs_ds_130_mu_pT1 );
    ws_bg  .import( bg_ds_130_mu_pT1    );
    ws_data.import( data_ds_130_mu_pT1  );
  
    ws_sig .import( higgs_ds_125_mu_pT1 );
    ws_bg  .import( bg_ds_125_mu_pT1    );
    ws_data.import( data_ds_125_mu_pT1  );
  
    ws_sig .import( higgs_ds_120_mu_pT1 );
    ws_bg  .import( bg_ds_120_mu_pT1    );
    ws_data.import( data_ds_120_mu_pT1  );

    dataAmount << mu_pT1RegionName << "    " << data_ds_125_mu_pT1.sumEntries() << std::endl;

  }
  if(doBBbin){
    /////////////////bb Z window///////////////////////
    ((RooRealVar*)higgs_ds_130_bbz_pT0.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130_bbz_pT0 .addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_125_bbz_pT0.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125_bbz_pT0 .addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_120_bbz_pT0.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120_bbz_pT0 .addColumn( hgg_mass_120 ))->setRange(100.,180.);

    ws_sig .import( higgs_ds_130_bbz_pT0 );    ws_bg  .import( bg_ds_130_bbz_pT0    );    ws_data.import( data_ds_130_bbz_pT0  );
    ws_sig .import( higgs_ds_125_bbz_pT0 );    ws_bg  .import( bg_ds_125_bbz_pT0    );    ws_data.import( data_ds_125_bbz_pT0  );
    ws_sig .import( higgs_ds_120_bbz_pT0 );    ws_bg  .import( bg_ds_120_bbz_pT0    );    ws_data.import( data_ds_120_bbz_pT0  );

    dataAmount << bbz_pT0RegionName << "    " << data_ds_125_bbz_pT0.sumEntries() << std::endl;

  ((RooRealVar*)higgs_ds_130_bbz_pT1.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130_bbz_pT1 .addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_125_bbz_pT1.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125_bbz_pT1 .addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_120_bbz_pT1.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120_bbz_pT1 .addColumn( hgg_mass_120 ))->setRange(100.,180.);

    ws_sig .import( higgs_ds_130_bbz_pT1 );    ws_bg  .import( bg_ds_130_bbz_pT1    );    ws_data.import( data_ds_130_bbz_pT1  );
    ws_sig .import( higgs_ds_125_bbz_pT1 );    ws_bg  .import( bg_ds_125_bbz_pT1    );    ws_data.import( data_ds_125_bbz_pT1  );
    ws_sig .import( higgs_ds_120_bbz_pT1 );    ws_bg  .import( bg_ds_120_bbz_pT1    );    ws_data.import( data_ds_120_bbz_pT1  );

    dataAmount << bbz_pT1RegionName << "    " << data_ds_125_bbz_pT1.sumEntries() << std::endl;

    /////////////////bb H window///////////////////////
    ((RooRealVar*)higgs_ds_130_bbh_pT0.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130_bbh_pT0 .addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_125_bbh_pT0.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125_bbh_pT0 .addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_120_bbh_pT0.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120_bbh_pT0 .addColumn( hgg_mass_120 ))->setRange(100.,180.);

    ws_sig .import( higgs_ds_130_bbh_pT0 );    ws_bg  .import( bg_ds_130_bbh_pT0    );    ws_data.import( data_ds_130_bbh_pT0  );
    ws_sig .import( higgs_ds_125_bbh_pT0 );    ws_bg  .import( bg_ds_125_bbh_pT0    );    ws_data.import( data_ds_125_bbh_pT0  );
    ws_sig .import( higgs_ds_120_bbh_pT0 );    ws_bg  .import( bg_ds_120_bbh_pT0    );    ws_data.import( data_ds_120_bbh_pT0  );

    dataAmount << bbh_pT0RegionName << "    " << data_ds_125_bbh_pT0.sumEntries() << std::endl;

    ((RooRealVar*)higgs_ds_130_bbh_pT1.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130_bbh_pT1 .addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_125_bbh_pT1.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125_bbh_pT1 .addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)higgs_ds_120_bbh_pT1.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120_bbh_pT1 .addColumn( hgg_mass_120 ))->setRange(100.,180.);

    ws_sig .import( higgs_ds_130_bbh_pT1 );    ws_bg  .import( bg_ds_130_bbh_pT1    );    ws_data.import( data_ds_130_bbh_pT1  );
    ws_sig .import( higgs_ds_125_bbh_pT1 );    ws_bg  .import( bg_ds_125_bbh_pT1    );    ws_data.import( data_ds_125_bbh_pT1  );
    ws_sig .import( higgs_ds_120_bbh_pT1 );    ws_bg  .import( bg_ds_120_bbh_pT1    );    ws_data.import( data_ds_120_bbh_pT1  );

    dataAmount << bbh_pT1RegionName << "    " << data_ds_125_bbh_pT1.sumEntries() << std::endl;

  }
 
  for( unsigned int i=0; i<signalList.size(); i++){
    if( doBBbin ){
      /////////////////bb Z window///////////////////////
      ((RooRealVar*)sig_ds_130_bbz_pT0[i].addColumn( hgg_mass_130 ))->setRange(100.,180.);  
      ((RooRealVar*)sig_ds_125_bbz_pT0[i].addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)sig_ds_120_bbz_pT0[i].addColumn( hgg_mass_120 ))->setRange(100.,180.); 
      ws_sig.import( sig_ds_130_bbz_pT0[i]);
      ws_sig.import( sig_ds_125_bbz_pT0[i]);
      ws_sig.import( sig_ds_120_bbz_pT0[i]);
      ((RooRealVar*)sig_ds_130_bbz_pT1[i].addColumn( hgg_mass_130 ))->setRange(100.,180.);  
      ((RooRealVar*)sig_ds_125_bbz_pT1[i].addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)sig_ds_120_bbz_pT1[i].addColumn( hgg_mass_120 ))->setRange(100.,180.); 
      ws_sig.import( sig_ds_130_bbz_pT1[i]);
      ws_sig.import( sig_ds_125_bbz_pT1[i]);
      ws_sig.import( sig_ds_120_bbz_pT1[i]);
      /////////////////bb H window///////////////////////
      ((RooRealVar*)sig_ds_130_bbh_pT0[i].addColumn( hgg_mass_130 ))->setRange(100.,180.);  
      ((RooRealVar*)sig_ds_125_bbh_pT0[i].addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)sig_ds_120_bbh_pT0[i].addColumn( hgg_mass_120 ))->setRange(100.,180.); 
      ws_sig.import( sig_ds_130_bbh_pT0[i]);
      ws_sig.import( sig_ds_125_bbh_pT0[i]);
      ws_sig.import( sig_ds_120_bbh_pT0[i]);
      ((RooRealVar*)sig_ds_130_bbh_pT1[i].addColumn( hgg_mass_130 ))->setRange(100.,180.);  
      ((RooRealVar*)sig_ds_125_bbh_pT1[i].addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)sig_ds_120_bbh_pT1[i].addColumn( hgg_mass_120 ))->setRange(100.,180.); 
      ws_sig.import( sig_ds_130_bbh_pT1[i]);
      ws_sig.import( sig_ds_125_bbh_pT1[i]);
      ws_sig.import( sig_ds_120_bbh_pT1[i]);

      if( i==0 ){
	sigAmount <<  signalList[i] << " \t " << bbz_pT0RegionName << " \t " <<  sig_ds_125_bbz_pT0[i].sumEntries()  << std::endl;
	sigAmount <<  signalList[i] << " \t " << bbh_pT0RegionName << " \t " <<  sig_ds_125_bbh_pT0[i].sumEntries()  << std::endl;
	sigAmount <<  signalList[i] << " \t " << bbz_pT1RegionName << " \t " <<  sig_ds_125_bbz_pT1[i].sumEntries()  << std::endl;
	sigAmount <<  signalList[i] << " \t " << bbh_pT1RegionName << " \t " <<  sig_ds_125_bbh_pT1[i].sumEntries()  << std::endl;
      }
    }
    if( doLLbin ){
      ((RooRealVar*)sig_ds_130_ll[i].addColumn( hgg_mass_130 ))->setRange(100.,180.);  
      ((RooRealVar*)sig_ds_125_ll[i].addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)sig_ds_120_ll[i].addColumn( hgg_mass_120 ))->setRange(100.,180.); 
      ws_sig.import( sig_ds_130_ll[i]);
      ws_sig.import( sig_ds_125_ll[i]);
      ws_sig.import( sig_ds_120_ll[i]);
      if( i==0 ){
	sigAmount <<  signalList[i] << " \t " << llRegionName << " \t " <<  sig_ds_125_ll[i].sumEntries()  << std::endl;
      }
    }      

    if( doSingleLbin ){

      ((RooRealVar*)sig_ds_130_el_pT0[i].addColumn( hgg_mass_130 ))->setRange(100.,180.);  
      ((RooRealVar*)sig_ds_125_el_pT0[i].addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)sig_ds_120_el_pT0[i].addColumn( hgg_mass_120 ))->setRange(100.,180.); 
      ws_sig.import( sig_ds_130_el_pT0[i]);
      ws_sig.import( sig_ds_125_el_pT0[i]);
      ws_sig.import( sig_ds_120_el_pT0[i]);
      if( i==0 ){
	sigAmount <<  signalList[i] << " \t " << el_pT0RegionName << " \t " <<  sig_ds_125_el_pT0[i].sumEntries()  << std::endl;
      }

      ((RooRealVar*)sig_ds_130_el_pT1[i].addColumn( hgg_mass_130 ))->setRange(100.,180.);  
      ((RooRealVar*)sig_ds_125_el_pT1[i].addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)sig_ds_120_el_pT1[i].addColumn( hgg_mass_120 ))->setRange(100.,180.); 
      ws_sig.import( sig_ds_130_el_pT1[i]);
      ws_sig.import( sig_ds_125_el_pT1[i]);
      ws_sig.import( sig_ds_120_el_pT1[i]);
      if( i==0 ){
	sigAmount <<  signalList[i] << " \t " << el_pT1RegionName << " \t " <<  sig_ds_125_el_pT1[i].sumEntries()  << std::endl;
      }
 
      ((RooRealVar*)sig_ds_130_mu_pT0[i].addColumn( hgg_mass_130 ))->setRange(100.,180.);  
      ((RooRealVar*)sig_ds_125_mu_pT0[i].addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)sig_ds_120_mu_pT0[i].addColumn( hgg_mass_120 ))->setRange(100.,180.); 
      ws_sig.import( sig_ds_130_mu_pT0[i]);
      ws_sig.import( sig_ds_125_mu_pT0[i]);
      ws_sig.import( sig_ds_120_mu_pT0[i]);
      if( i==0 ){
	sigAmount <<  signalList[i] << " \t " << mu_pT0RegionName << " \t " <<  sig_ds_125_mu_pT0[i].sumEntries()  << std::endl;
      }

      ((RooRealVar*)sig_ds_130_mu_pT1[i].addColumn( hgg_mass_130 ))->setRange(100.,180.);  
      ((RooRealVar*)sig_ds_125_mu_pT1[i].addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)sig_ds_120_mu_pT1[i].addColumn( hgg_mass_120 ))->setRange(100.,180.); 
      ws_sig.import( sig_ds_130_mu_pT1[i]);
      ws_sig.import( sig_ds_125_mu_pT1[i]);
      ws_sig.import( sig_ds_120_mu_pT1[i]);
      if( i==0 ){
	sigAmount <<  signalList[i] << " \t " << mu_pT1RegionName << " \t " <<  sig_ds_125_mu_pT1[i].sumEntries()  << std::endl;
      }

    }


 
    if( doISRsyst ){
      
      if( doLLbin ){
	float sign = 1;
	float isr_uncert = max (fabs(h_isr_UP_ll_vec[i]->GetMean()-h_isr_ll_vec[i]->GetMean()), fabs( h_isr_ll_vec[i]->GetMean()-h_isr_DN_ll_vec[i]->GetMean()));

	if( fabs(h_isr_UP_ll_vec[i]->GetMean() - h_isr_ll_vec[i]->GetMean()) >  fabs(  h_isr_ll_vec[i]->GetMean() - h_isr_DN_ll_vec[i]->GetMean() )){
	  sign = ((h_isr_UP_ll_vec[i]->GetMean() - h_isr_ll_vec[i]->GetMean())  > 0 ) - ((h_isr_UP_ll_vec[i]->GetMean() - h_isr_ll_vec[i]->GetMean())  < 0) ;
	}else{
	  sign = ((h_isr_ll_vec[i]->GetMean() - h_isr_DN_ll_vec[i]->GetMean())  > 0 ) - ((h_isr_ll_vec[i]->GetMean() - h_isr_DN_ll_vec[i]->GetMean())  < 0) ;
	}
	isr <<  signalList[i] << " \t " << llRegionName << " \t " <<    sign* isr_uncert  <<std::endl;

	int massSbottom = 0;    int massChi = 0;

	std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	std::string::size_type sz;   // alias of size_t
	massChi = std::stoi (mLSP,&sz);
	std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	massSbottom = std::stoi (mSbottom,&sz);

	int binBottom = h2D_isr_ll->GetXaxis()->FindBin( (float)massSbottom );
	int binChi    = h2D_isr_ll->GetYaxis()->FindBin( (float)massChi );

	h2D_isr_ll->SetBinContent( binBottom, binChi, h_isr_ll_vec[i]->GetMean() );
	h2D_isr_err_ll->SetBinContent( binBottom, binChi, sign* isr_uncert );

	if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  sigAmount <<  signalList[i] << " \t " << llRegionName << " \t " <<  sig_ds_125_ll[i].sumEntries()  << std::endl;
      }  

      if( doSingleLbin ){
	float sign = 1;
	float isr_uncert = max (fabs(h_isr_UP_el_pT0_vec[i]->GetMean()-h_isr_el_pT0_vec[i]->GetMean()), fabs( h_isr_el_pT0_vec[i]->GetMean()-h_isr_DN_el_pT0_vec[i]->GetMean()));

	if( fabs(h_isr_UP_el_pT0_vec[i]->GetMean() - h_isr_el_pT0_vec[i]->GetMean()) >  fabs(  h_isr_el_pT0_vec[i]->GetMean() - h_isr_DN_el_pT0_vec[i]->GetMean() )){
	  sign = ((h_isr_UP_el_pT0_vec[i]->GetMean() - h_isr_el_pT0_vec[i]->GetMean())  > 0 ) - ((h_isr_UP_el_pT0_vec[i]->GetMean() - h_isr_el_pT0_vec[i]->GetMean())  < 0) ;
	}else{
	  sign = ((h_isr_el_pT0_vec[i]->GetMean() - h_isr_DN_el_pT0_vec[i]->GetMean())  > 0 ) - ((h_isr_el_pT0_vec[i]->GetMean() - h_isr_DN_el_pT0_vec[i]->GetMean())  < 0) ;
	}
	isr <<  signalList[i] << " \t " << el_pT0RegionName << " \t " <<    sign* isr_uncert  <<std::endl;

	int massSbottom = 0;    int massChi = 0;

	std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	std::string::size_type sz;   // alias of size_t
	massChi = std::stoi (mLSP,&sz);
	std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	massSbottom = std::stoi (mSbottom,&sz);

	int binBottom = h2D_isr_el_pT0->GetXaxis()->FindBin( (float)massSbottom );
	int binChi    = h2D_isr_el_pT0->GetYaxis()->FindBin( (float)massChi );

	h2D_isr_el_pT0->SetBinContent( binBottom, binChi, h_isr_el_pT0_vec[i]->GetMean() );
	h2D_isr_err_el_pT0->SetBinContent( binBottom, binChi, sign* isr_uncert );

	if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  sigAmount <<  signalList[i] << " \t " << el_pT0RegionName << " \t " <<  sig_ds_125_el_pT0[i].sumEntries()  << std::endl;
      }

      if( doSingleLbin ){
	float sign = 1;
	float isr_uncert = max (fabs(h_isr_UP_el_pT1_vec[i]->GetMean()-h_isr_el_pT1_vec[i]->GetMean()), fabs( h_isr_el_pT1_vec[i]->GetMean()-h_isr_DN_el_pT1_vec[i]->GetMean()));

	if( fabs(h_isr_UP_el_pT1_vec[i]->GetMean() - h_isr_el_pT1_vec[i]->GetMean()) >  fabs(  h_isr_el_pT1_vec[i]->GetMean() - h_isr_DN_el_pT1_vec[i]->GetMean() )){
	  sign = ((h_isr_UP_el_pT1_vec[i]->GetMean() - h_isr_el_pT1_vec[i]->GetMean())  > 0 ) - ((h_isr_UP_el_pT1_vec[i]->GetMean() - h_isr_el_pT1_vec[i]->GetMean())  < 0) ;
	}else{
	  sign = ((h_isr_el_pT1_vec[i]->GetMean() - h_isr_DN_el_pT1_vec[i]->GetMean())  > 0 ) - ((h_isr_el_pT1_vec[i]->GetMean() - h_isr_DN_el_pT1_vec[i]->GetMean())  < 0) ;
	}
	isr <<  signalList[i] << " \t " << el_pT1RegionName << " \t " <<    sign* isr_uncert  <<std::endl;

	int massSbottom = 0;    int massChi = 0;

	std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	std::string::size_type sz;   // alias of size_t
	massChi = std::stoi (mLSP,&sz);
	std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	massSbottom = std::stoi (mSbottom,&sz);

	int binBottom = h2D_isr_el_pT1->GetXaxis()->FindBin( (float)massSbottom );
	int binChi    = h2D_isr_el_pT1->GetYaxis()->FindBin( (float)massChi );

	h2D_isr_el_pT1->SetBinContent( binBottom, binChi, h_isr_el_pT1_vec[i]->GetMean() );
	h2D_isr_err_el_pT1->SetBinContent( binBottom, binChi, sign* isr_uncert );

	if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  sigAmount <<  signalList[i] << " \t " << el_pT1RegionName << " \t " <<  sig_ds_125_el_pT1[i].sumEntries()  << std::endl;
      }



      if( doSingleLbin ){
	float sign = 1;
	float isr_uncert = max (fabs(h_isr_UP_mu_pT0_vec[i]->GetMean()-h_isr_mu_pT0_vec[i]->GetMean()), fabs( h_isr_mu_pT0_vec[i]->GetMean()-h_isr_DN_mu_pT0_vec[i]->GetMean()));

	if( fabs(h_isr_UP_mu_pT0_vec[i]->GetMean() - h_isr_mu_pT0_vec[i]->GetMean()) >  fabs(  h_isr_mu_pT0_vec[i]->GetMean() - h_isr_DN_mu_pT0_vec[i]->GetMean() )){
	  sign = ((h_isr_UP_mu_pT0_vec[i]->GetMean() - h_isr_mu_pT0_vec[i]->GetMean())  > 0 ) - ((h_isr_UP_mu_pT0_vec[i]->GetMean() - h_isr_mu_pT0_vec[i]->GetMean())  < 0) ;
	}else{
	  sign = ((h_isr_mu_pT0_vec[i]->GetMean() - h_isr_DN_mu_pT0_vec[i]->GetMean())  > 0 ) - ((h_isr_mu_pT0_vec[i]->GetMean() - h_isr_DN_mu_pT0_vec[i]->GetMean())  < 0) ;
	}
	isr <<  signalList[i] << " \t " << mu_pT0RegionName << " \t " <<    sign* isr_uncert  <<std::endl;

	int massSbottom = 0;    int massChi = 0;

	std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	std::string::size_type sz;   // alias of size_t
	massChi = std::stoi (mLSP,&sz);
	std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	massSbottom = std::stoi (mSbottom,&sz);

	int binBottom = h2D_isr_mu_pT0->GetXaxis()->FindBin( (float)massSbottom );
	int binChi    = h2D_isr_mu_pT0->GetYaxis()->FindBin( (float)massChi );

	h2D_isr_mu_pT0->SetBinContent( binBottom, binChi, h_isr_mu_pT0_vec[i]->GetMean() );
	h2D_isr_err_mu_pT0->SetBinContent( binBottom, binChi, sign* isr_uncert );

	if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  sigAmount <<  signalList[i] << " \t " << mu_pT0RegionName << " \t " <<  sig_ds_125_mu_pT0[i].sumEntries()  << std::endl;
      }


      if( doSingleLbin ){
	float sign = 1;
	float isr_uncert = max (fabs(h_isr_UP_mu_pT1_vec[i]->GetMean()-h_isr_mu_pT1_vec[i]->GetMean()), fabs( h_isr_mu_pT1_vec[i]->GetMean()-h_isr_DN_mu_pT1_vec[i]->GetMean()));

	if( fabs(h_isr_UP_mu_pT1_vec[i]->GetMean() - h_isr_mu_pT1_vec[i]->GetMean()) >  fabs(  h_isr_mu_pT1_vec[i]->GetMean() - h_isr_DN_mu_pT1_vec[i]->GetMean() )){
	  sign = ((h_isr_UP_mu_pT1_vec[i]->GetMean() - h_isr_mu_pT1_vec[i]->GetMean())  > 0 ) - ((h_isr_UP_mu_pT1_vec[i]->GetMean() - h_isr_mu_pT1_vec[i]->GetMean())  < 0) ;
	}else{
	  sign = ((h_isr_mu_pT1_vec[i]->GetMean() - h_isr_DN_mu_pT1_vec[i]->GetMean())  > 0 ) - ((h_isr_mu_pT1_vec[i]->GetMean() - h_isr_DN_mu_pT1_vec[i]->GetMean())  < 0) ;
	}
	isr <<  signalList[i] << " \t " << mu_pT1RegionName << " \t " <<    sign* isr_uncert  <<std::endl;

	int massSbottom = 0;    int massChi = 0;

	std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	std::string::size_type sz;   // alias of size_t
	massChi = std::stoi (mLSP,&sz);
	std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	massSbottom = std::stoi (mSbottom,&sz);

	int binBottom = h2D_isr_mu_pT1->GetXaxis()->FindBin( (float)massSbottom );
	int binChi    = h2D_isr_mu_pT1->GetYaxis()->FindBin( (float)massChi );

	h2D_isr_mu_pT1->SetBinContent( binBottom, binChi, h_isr_mu_pT1_vec[i]->GetMean() );
	h2D_isr_err_mu_pT1->SetBinContent( binBottom, binChi, sign* isr_uncert );

	if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  sigAmount <<  signalList[i] << " \t " << mu_pT1RegionName << " \t " <<  sig_ds_125_mu_pT1[i].sumEntries()  << std::endl;
      }

      if( doBBbin ){
	/////////////////bb Z window///////////////////////
	float sign_bbz_pT0 = 1;
	float isr_uncert_bbz_pT0 = max(fabs(h_isr_UP_bbz_pT0_vec[i]->GetMean()-h_isr_bbz_pT0_vec[i]->GetMean()), fabs( h_isr_bbz_pT0_vec[i]->GetMean()-h_isr_DN_bbz_pT0_vec[i]->GetMean()));

	if( fabs(h_isr_UP_bbz_pT0_vec[i]->GetMean() - h_isr_bbz_pT0_vec[i]->GetMean()) >  fabs(  h_isr_bbz_pT0_vec[i]->GetMean() - h_isr_DN_bbz_pT0_vec[i]->GetMean() )){
	  sign_bbz_pT0 = ((h_isr_UP_bbz_pT0_vec[i]->GetMean() - h_isr_bbz_pT0_vec[i]->GetMean())  > 0 ) - ((h_isr_UP_bbz_pT0_vec[i]->GetMean() - h_isr_bbz_pT0_vec[i]->GetMean())  < 0) ;
	}else{
	  sign_bbz_pT0 = ((h_isr_bbz_pT0_vec[i]->GetMean() - h_isr_DN_bbz_pT0_vec[i]->GetMean())  > 0 ) - ((h_isr_bbz_pT0_vec[i]->GetMean() - h_isr_DN_bbz_pT0_vec[i]->GetMean())  < 0) ;
	}
	isr <<  signalList[i] << " \t " << bbz_pT0RegionName << " \t " <<    sign_bbz_pT0 * isr_uncert_bbz_pT0  <<std::endl;

	int massSbottom = 0;    int massChi = 0;

	std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	std::string::size_type sz;   // alias of size_t
	massChi = std::stoi (mLSP,&sz);
	std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	massSbottom = std::stoi (mSbottom,&sz);

	int binBottom = h2D_isr_bbz_pT0->GetXaxis()->FindBin( (float)massSbottom );
	int binChi    = h2D_isr_bbz_pT0->GetYaxis()->FindBin( (float)massChi );

	h2D_isr_bbz_pT0->SetBinContent( binBottom, binChi, h_isr_bbz_pT0_vec[i]->GetMean() );
	h2D_isr_err_bbz_pT0->SetBinContent( binBottom, binChi, sign_bbz_pT0* isr_uncert_bbz_pT0 );

	if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  sigAmount <<  signalList[i] << " \t " << bbz_pT0RegionName << " \t " <<  sig_ds_125_bbz_pT0[i].sumEntries()  << std::endl;
    

	/////////////////bb H window///////////////////////
	float sign_bbz_pT1 = 1;
	float isr_uncert_bbz_pT1 = max(fabs(h_isr_UP_bbz_pT1_vec[i]->GetMean()-h_isr_bbz_pT1_vec[i]->GetMean()), fabs( h_isr_bbz_pT1_vec[i]->GetMean()-h_isr_DN_bbz_pT1_vec[i]->GetMean()));

	if( fabs(h_isr_UP_bbz_pT1_vec[i]->GetMean() - h_isr_bbz_pT1_vec[i]->GetMean()) >  fabs(  h_isr_bbz_pT1_vec[i]->GetMean() - h_isr_DN_bbz_pT1_vec[i]->GetMean() )){
	  sign_bbz_pT1 = ((h_isr_UP_bbz_pT1_vec[i]->GetMean() - h_isr_bbz_pT1_vec[i]->GetMean())  > 0 ) - ((h_isr_UP_bbz_pT1_vec[i]->GetMean() - h_isr_bbz_pT1_vec[i]->GetMean())  < 0) ;
	}else{
	  sign_bbz_pT1 = ((h_isr_bbz_pT1_vec[i]->GetMean() - h_isr_DN_bbz_pT1_vec[i]->GetMean())  > 0 ) - ((h_isr_bbz_pT1_vec[i]->GetMean() - h_isr_DN_bbz_pT1_vec[i]->GetMean())  < 0) ;
	}
	isr <<  signalList[i] << " \t " << bbz_pT1RegionName << " \t " <<    sign_bbz_pT1* isr_uncert_bbz_pT1  <<std::endl;

	h2D_isr_bbz_pT1->SetBinContent( binBottom, binChi, h_isr_bbz_pT1_vec[i]->GetMean() );
	h2D_isr_err_bbz_pT1->SetBinContent( binBottom, binChi, sign_bbz_pT1* isr_uncert_bbz_pT1 );

	if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  sigAmount <<  signalList[i] << " \t " << bbz_pT1RegionName << " \t " <<  sig_ds_125_bbz_pT1[i].sumEntries()  << std::endl;

	/////////////////bb H window///////////////////////
	float sign_bbh_pT0 = 1;
	float isr_uncert_bbh_pT0 = max(fabs(h_isr_UP_bbh_pT0_vec[i]->GetMean()-h_isr_bbh_pT0_vec[i]->GetMean()), fabs( h_isr_bbh_pT0_vec[i]->GetMean()-h_isr_DN_bbh_pT0_vec[i]->GetMean()));

	if( fabs(h_isr_UP_bbh_pT0_vec[i]->GetMean() - h_isr_bbh_pT0_vec[i]->GetMean()) >  fabs(  h_isr_bbh_pT0_vec[i]->GetMean() - h_isr_DN_bbh_pT0_vec[i]->GetMean() )){
	  sign_bbh_pT0 = ((h_isr_UP_bbh_pT0_vec[i]->GetMean() - h_isr_bbh_pT0_vec[i]->GetMean())  > 0 ) - ((h_isr_UP_bbh_pT0_vec[i]->GetMean() - h_isr_bbh_pT0_vec[i]->GetMean())  < 0) ;
	}else{
	  sign_bbh_pT0 = ((h_isr_bbh_pT0_vec[i]->GetMean() - h_isr_DN_bbh_pT0_vec[i]->GetMean())  > 0 ) - ((h_isr_bbh_pT0_vec[i]->GetMean() - h_isr_DN_bbh_pT0_vec[i]->GetMean())  < 0) ;
	}
	isr <<  signalList[i] << " \t " << bbh_pT0RegionName << " \t " <<    sign_bbh_pT0* isr_uncert_bbh_pT0  <<std::endl;

	h2D_isr_bbh_pT0->SetBinContent( binBottom, binChi, h_isr_bbh_pT0_vec[i]->GetMean() );
	h2D_isr_err_bbh_pT0->SetBinContent( binBottom, binChi, sign_bbh_pT0* isr_uncert_bbh_pT0 );

	if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  sigAmount <<  signalList[i] << " \t " << bbh_pT0RegionName << " \t " <<  sig_ds_125_bbh_pT0[i].sumEntries()  << std::endl;

	/////////////////bb H window///////////////////////
	float sign_bbh_pT1 = 1;
	float isr_uncert_bbh_pT1 = max(fabs(h_isr_UP_bbh_pT1_vec[i]->GetMean()-h_isr_bbh_pT1_vec[i]->GetMean()), fabs( h_isr_bbh_pT1_vec[i]->GetMean()-h_isr_DN_bbh_pT1_vec[i]->GetMean()));

	if( fabs(h_isr_UP_bbh_pT1_vec[i]->GetMean() - h_isr_bbh_pT1_vec[i]->GetMean()) >  fabs(  h_isr_bbh_pT1_vec[i]->GetMean() - h_isr_DN_bbh_pT1_vec[i]->GetMean() )){
	  sign_bbh_pT1 = ((h_isr_UP_bbh_pT1_vec[i]->GetMean() - h_isr_bbh_pT1_vec[i]->GetMean())  > 0 ) - ((h_isr_UP_bbh_pT1_vec[i]->GetMean() - h_isr_bbh_pT1_vec[i]->GetMean())  < 0) ;
	}else{
	  sign_bbh_pT1 = ((h_isr_bbh_pT1_vec[i]->GetMean() - h_isr_DN_bbh_pT1_vec[i]->GetMean())  > 0 ) - ((h_isr_bbh_pT1_vec[i]->GetMean() - h_isr_DN_bbh_pT1_vec[i]->GetMean())  < 0) ;
	}
	isr <<  signalList[i] << " \t " << bbh_pT1RegionName << " \t " <<    sign_bbh_pT1* isr_uncert_bbh_pT1  <<std::endl;

	h2D_isr_bbh_pT1->SetBinContent( binBottom, binChi, h_isr_bbh_pT1_vec[i]->GetMean() );
	h2D_isr_err_bbh_pT1->SetBinContent( binBottom, binChi, sign_bbh_pT1* isr_uncert_bbh_pT1 );

	if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  sigAmount <<  signalList[i] << " \t " << bbh_pT1RegionName << " \t " <<  sig_ds_125_bbh_pT1[i].sumEntries()  << std::endl;

      }

    }//done if isr
    
  }//done loop over signals

  //  delete h_isr_ll; delete h_isr_UP_ll; delete h_isr_DN_ll;
  //  delete h_isr_bbz_pT0; delete h_isr_UP_bbz_pT0; delete h_isr_DN_bbz_pT0;



  if( doBBbin && doISRsyst ){
    h2D_isr_bbz_pT0->Write();
    h2D_isr_err_bbz_pT0->Write();
    h2D_isr_bbh_pT0->Write();
    h2D_isr_err_bbh_pT0->Write();

    h2D_isr_bbz_pT1->Write();
    h2D_isr_err_bbz_pT1->Write();
    h2D_isr_bbh_pT1->Write();
    h2D_isr_err_bbh_pT1->Write();
  }



  if( doLLbin && doISRsyst ){

    h2D_isr_ll->Write();
    h2D_isr_err_ll->Write();
  }
  
  if( doSingleLbin && doISRsyst ){
    h2D_isr_el_pT0->Write();
    h2D_isr_err_el_pT0->Write(); 
    h2D_isr_el_pT1->Write();
    h2D_isr_err_el_pT1->Write(); 
 
    h2D_isr_mu_pT0->Write();
    h2D_isr_err_mu_pT0->Write(); 
    h2D_isr_mu_pT1->Write();
    h2D_isr_err_mu_pT1->Write(); 
  }


  dataAmount.close();

  isr.close();
 
  ws_data.Write();
  ws_sig.Write();
  ws_bg.Write();

  dataFile->Close();

  return 0;
}






















    // Float_t h_mass_input;
    // tree_higgs->SetBranchAddress("h_mass", &h_mass_input);
    // tree_sig->SetBranchAddress("h_mass", &h_mass_input);
    // tree_qcd->SetBranchAddress("h_mass", &h_mass_input);
    // tree_diPhoton->SetBranchAddress("h_mass", &h_mass_input);
    // tree_gjets->SetBranchAddress("h_mass", &h_mass_input);

    // tree_higgs->SetAlias( "h_mass_input", "h_mass");
    // tree_sig->SetAlias( "h_mass_input", "h_mass");
    // tree_qcd->SetAlias( "h_mass_input", "h_mass");
    // tree_diPhoton->SetAlias( "h_mass_input", "h_mass");
    // tree_gjets->SetAlias( "h_mass_input", "h_mass");
    // tree_higgs->GetEntry(50);
    // std::cout << h_mass_input << std::endl;
