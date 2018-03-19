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

bool doMT2binning = false;

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


  vector<std::string> signalList_T2bH = { "SMS_T2bH_mSbottom250_mLSP100", "SMS_T2bH_mSbottom250_mLSP1", "SMS_T2bH_mSbottom250_mLSP50", "SMS_T2bH_mSbottom300_mLSP100", "SMS_T2bH_mSbottom300_mLSP150", "SMS_T2bH_mSbottom300_mLSP1", "SMS_T2bH_mSbottom300_mLSP50", "SMS_T2bH_mSbottom350_mLSP100", "SMS_T2bH_mSbottom350_mLSP150", "SMS_T2bH_mSbottom350_mLSP1", "SMS_T2bH_mSbottom350_mLSP200", "SMS_T2bH_mSbottom350_mLSP50", "SMS_T2bH_mSbottom400_mLSP100", "SMS_T2bH_mSbottom400_mLSP150", "SMS_T2bH_mSbottom400_mLSP1", "SMS_T2bH_mSbottom400_mLSP200", "SMS_T2bH_mSbottom400_mLSP250", "SMS_T2bH_mSbottom400_mLSP50", "SMS_T2bH_mSbottom450_mLSP1", "SMS_T2bH_mSbottom450_mLSP100", "SMS_T2bH_mSbottom450_mLSP150", "SMS_T2bH_mSbottom450_mLSP200", "SMS_T2bH_mSbottom450_mLSP250", "SMS_T2bH_mSbottom450_mLSP300", "SMS_T2bH_mSbottom450_mLSP50", "SMS_T2bH_mSbottom500_mLSP100", "SMS_T2bH_mSbottom500_mLSP150", "SMS_T2bH_mSbottom500_mLSP200", "SMS_T2bH_mSbottom500_mLSP250", "SMS_T2bH_mSbottom500_mLSP300", "SMS_T2bH_mSbottom500_mLSP50", "SMS_T2bH_mSbottom600_mLSP1", "SMS_T2bH_mSbottom600_mLSP100", "SMS_T2bH_mSbottom600_mLSP200", "SMS_T2bH_mSbottom600_mLSP300" };
  //  vector<std::string> signalList_T2bH = { "SMS_T2bH_mSbottom250_mLSP100", "SMS_T2bH_mSbottom250_mLSP1", "SMS_T2bH_mSbottom250_mLSP50", "SMS_T2bH_mSbottom300_mLSP100", "SMS_T2bH_mSbottom300_mLSP150", "SMS_T2bH_mSbottom300_mLSP1", "SMS_T2bH_mSbottom300_mLSP50", "SMS_T2bH_mSbottom350_mLSP100", "SMS_T2bH_mSbottom350_mLSP150", "SMS_T2bH_mSbottom350_mLSP1", "SMS_T2bH_mSbottom350_mLSP200", "SMS_T2bH_mSbottom350_mLSP50", "SMS_T2bH_mSbottom400_mLSP100", "SMS_T2bH_mSbottom400_mLSP150", "SMS_T2bH_mSbottom400_mLSP1", "SMS_T2bH_mSbottom400_mLSP200", "SMS_T2bH_mSbottom400_mLSP250", "SMS_T2bH_mSbottom400_mLSP50", "SMS_T2bH_mSbottom450_mLSP1", "SMS_T2bH_mSbottom450_mLSP100", "SMS_T2bH_mSbottom450_mLSP150", "SMS_T2bH_mSbottom450_mLSP200", "SMS_T2bH_mSbottom450_mLSP250", "SMS_T2bH_mSbottom450_mLSP300", "SMS_T2bH_mSbottom450_mLSP50", "SMS_T2bH_mSbottom500_mLSP100", "SMS_T2bH_mSbottom500_mLSP150", "SMS_T2bH_mSbottom500_mLSP1", "SMS_T2bH_mSbottom500_mLSP200", "SMS_T2bH_mSbottom500_mLSP250", "SMS_T2bH_mSbottom500_mLSP300", "SMS_T2bH_mSbottom500_mLSP50", "SMS_T2bH_mSbottom600_mLSP1", "SMS_T2bH_mSbottom600_mLSP100", "SMS_T2bH_mSbottom600_mLSP200", "SMS_T2bH_mSbottom600_mLSP300" };

  vector<std::string> signalList_WH = { "SMS_TChiWH_HToGG_127_1","SMS_TChiWH_HToGG_150_1","SMS_TChiWH_HToGG_150_24"};
  //  vector<std::string> signalList_WH = { "SMS_TChiWH_HToGG_127_1","SMS_TChiWH_HToGG_150_1","SMS_TChiWH_HToGG_150_24","SMS_TChiWH_HToGG_175_1","SMS_TChiWH_HToGG_175_25","SMS_TChiWH_HToGG_175_49","SMS_TChiWH_HToGG_200_1","SMS_TChiWH_HToGG_200_25","SMS_TChiWH_HToGG_200_50","SMS_TChiWH_HToGG_200_74" };

  vector<std::string> signalList_HZ = { "SMS_TChiHZ_HToGG_m127", "SMS_TChiHZ_HToGG_m150" };
  // vector<std::string> signalList_HZ = { "SMS_TChiHZ_HToGG_m127", "SMS_TChiHZ_HToGG_m150", "SMS_TChiHZ_HToGG_m175", "SMS_TChiHZ_HToGG_m200" };
  
  //  vector<std::string> signalList_HZ = { "SMS_TChiHZ_HToGG_m127", "SMS_TChiHZ_HToGG_m150", "SMS_TChiHZ_HToGG_m175", "SMS_TChiHZ_HToGG_m200", "SMS_TChiHZ_HToGG_m225", "SMS_TChiHZ_HToGG_m250", "SMS_TChiHZ_HToGG_m275", "SMS_TChiHZ_HToGG_m300", "SMS_TChiHZ_HToGG_m325", "SMS_TChiHZ_HToGG_m350", "SMS_TChiHZ_HToGG_m375", "SMS_TChiHZ_HToGG_m400", "SMS_TChiHZ_HToGG_m425", "SMS_TChiHZ_HToGG_m450", "SMS_TChiHZ_HToGG_m475", "SMS_TChiHZ_HToGG_m500", "SMS_TChiHZ_HToGG_m525", "SMS_TChiHZ_HToGG_m550", "SMS_TChiHZ_HToGG_m575", "SMS_TChiHZ_HToGG_m600", "SMS_TChiHZ_HToGG_m625", "SMS_TChiHZ_HToGG_m650", "SMS_TChiHZ_HToGG_m675", "SMS_TChiHZ_HToGG_m700", "SMS_TChiHZ_HToGG_m725", "SMS_TChiHZ_HToGG_m750", "SMS_TChiHZ_HToGG_m775", "SMS_TChiHZ_HToGG_m800", "SMS_TChiHZ_HToGG_m825", "SMS_TChiHZ_HToGG_m850", "SMS_TChiHZ_HToGG_m875", "SMS_TChiHZ_HToGG_m900", "SMS_TChiHZ_HToGG_m925", "SMS_TChiHZ_HToGG_m950", "SMS_TChiHZ_HToGG_m975", "SMS_TChiHZ_HToGG_m1000" };

  //vector<std::string> signalList_HH_HH0p25 = { "SMS_TChiHH_HToGG_m127_HH0p25", "SMS_TChiHH_HToGG_m150_HH0p25"};
  vector<std::string> signalList_HH_HH0p25 = { "SMS_TChiHH_HToGG_m127_HH0p25", "SMS_TChiHH_HToGG_m150_HH0p25", "SMS_TChiHH_HToGG_m175_HH0p25", "SMS_TChiHH_HToGG_m200_HH0p25"};

  //  vector<std::string> signalList_HH_HH0p25 = { "SMS_TChiHH_HToGG_m127_HH0p25", "SMS_TChiHH_HToGG_m150_HH0p25", "SMS_TChiHH_HToGG_m175_HH0p25", "SMS_TChiHH_HToGG_m200_HH0p25", "SMS_TChiHH_HToGG_m225_HH0p25", "SMS_TChiHH_HToGG_m250_HH0p25", "SMS_TChiHH_HToGG_m275_HH0p25", "SMS_TChiHH_HToGG_m300_HH0p25", "SMS_TChiHH_HToGG_m325_HH0p25", "SMS_TChiHH_HToGG_m350_HH0p25", "SMS_TChiHH_HToGG_m375_HH0p25", "SMS_TChiHH_HToGG_m400_HH0p25", "SMS_TChiHH_HToGG_m425_HH0p25", "SMS_TChiHH_HToGG_m450_HH0p25", "SMS_TChiHH_HToGG_m475_HH0p25", "SMS_TChiHH_HToGG_m500_HH0p25", "SMS_TChiHH_HToGG_m525_HH0p25", "SMS_TChiHH_HToGG_m550_HH0p25", "SMS_TChiHH_HToGG_m575_HH0p25", "SMS_TChiHH_HToGG_m600_HH0p25", "SMS_TChiHH_HToGG_m625_HH0p25", "SMS_TChiHH_HToGG_m650_HH0p25", "SMS_TChiHH_HToGG_m675_HH0p25", "SMS_TChiHH_HToGG_m700_HH0p25", "SMS_TChiHH_HToGG_m725_HH0p25", "SMS_TChiHH_HToGG_m750_HH0p25", "SMS_TChiHH_HToGG_m775_HH0p25", "SMS_TChiHH_HToGG_m800_HH0p25", "SMS_TChiHH_HToGG_m825_HH0p25", "SMS_TChiHH_HToGG_m850_HH0p25", "SMS_TChiHH_HToGG_m875_HH0p25", "SMS_TChiHH_HToGG_m900_HH0p25", "SMS_TChiHH_HToGG_m925_HH0p25", "SMS_TChiHH_HToGG_m950_HH0p25", "SMS_TChiHH_HToGG_m975_HH0p25", "SMS_TChiHH_HToGG_m1000_HH0p25" };


  //  vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127"};
  vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127", "SMS_TChiHH_HToGG_m150", "SMS_TChiHH_HToGG_m175", "SMS_TChiHH_HToGG_m200"};

  //vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127", "SMS_TChiHH_HToGG_m150", "SMS_TChiHH_HToGG_m175", "SMS_TChiHH_HToGG_m200", "SMS_TChiHH_HToGG_m225", "SMS_TChiHH_HToGG_m250", "SMS_TChiHH_HToGG_m275", "SMS_TChiHH_HToGG_m300", "SMS_TChiHH_HToGG_m325", "SMS_TChiHH_HToGG_m350", "SMS_TChiHH_HToGG_m375", "SMS_TChiHH_HToGG_m400", "SMS_TChiHH_HToGG_m425", "SMS_TChiHH_HToGG_m450", "SMS_TChiHH_HToGG_m475", "SMS_TChiHH_HToGG_m500", "SMS_TChiHH_HToGG_m525", "SMS_TChiHH_HToGG_m550", "SMS_TChiHH_HToGG_m575", "SMS_TChiHH_HToGG_m600", "SMS_TChiHH_HToGG_m625", "SMS_TChiHH_HToGG_m650", "SMS_TChiHH_HToGG_m675", "SMS_TChiHH_HToGG_m700", "SMS_TChiHH_HToGG_m725", "SMS_TChiHH_HToGG_m750", "SMS_TChiHH_HToGG_m775", "SMS_TChiHH_HToGG_m800", "SMS_TChiHH_HToGG_m825", "SMS_TChiHH_HToGG_m850", "SMS_TChiHH_HToGG_m875", "SMS_TChiHH_HToGG_m900", "SMS_TChiHH_HToGG_m925", "SMS_TChiHH_HToGG_m950", "SMS_TChiHH_HToGG_m975", "SMS_TChiHH_HToGG_m1000" };

  std::vector<std::string> regionSelection;

  regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && ((h_pt/h_mass) <  0.8)) " );
  regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && ((h_pt/h_mass) >= 0.8)) " );

  regionSelection.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 0 && ((h_pt/h_mass) <  0.8)) " );
  regionSelection.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 0 && ((h_pt/h_mass) >= 0.8)) " );
  regionSelection.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 1 && ((h_pt/h_mass) <  0.8)) " );
  regionSelection.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 1 && ((h_pt/h_mass) >= 0.8)) " );
  regionSelection.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets >= 2 && ((h_pt/h_mass) <  0.8)) " );
  regionSelection.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets >= 2 && ((h_pt/h_mass) >= 0.8)) " );

  regionSelection.push_back(" (gg_nJets>=4 && gg_nBJets == 0 && ((h_pt/h_mass) <  0.8)) " );
  regionSelection.push_back(" (gg_nJets>=4 && gg_nBJets == 0 && ((h_pt/h_mass) >= 0.8)) " );
  regionSelection.push_back(" (gg_nJets>=4 && gg_nBJets == 1 && ((h_pt/h_mass) <  0.8)) " );
  regionSelection.push_back(" (gg_nJets>=4 && gg_nBJets == 1 && ((h_pt/h_mass) >= 0.8)) " );
  regionSelection.push_back(" (gg_nJets>=4 && gg_nBJets >= 2 && ((h_pt/h_mass) <  0.8)) " );
  regionSelection.push_back(" (gg_nJets>=4 && gg_nBJets >= 2 && ((h_pt/h_mass) >= 0.8)) " );

  for( unsigned int i=0; i < regionSelection.size(); i++){
  //   std::cout << "selection pre " << regionSelection[i] << std::endl;
    regionSelection[i] +=  "  && ( (!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )";
  //   std::cout << "selection post " << regionSelection[i] << std::endl;
  }

  
  regionSelection.push_back(" ( is1El) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  regionSelection.push_back(" ( is1Mu) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(is1El) && !(isDiLepZ) )"  );
  regionSelection.push_back(" ( is1El) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  regionSelection.push_back(" ( is1Mu) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(is1El) && !(isDiLepZ) )"  );

  regionSelection.push_back(" ( isDiLepZ ) && ( (!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) )"     );

  regionSelection.push_back(" ( isDiBH ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  regionSelection.push_back(" ( isDiBZ ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBH) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  regionSelection.push_back(" ( isDiBH ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  regionSelection.push_back(" ( isDiBZ ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBH) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );



  std::vector<std::string> regionNames;
  regionNames.push_back("j0_b0toInf_pT0");
  regionNames.push_back("j0_b0toInf_pT1");
  regionNames.push_back("j1to3_b0_pT0");
  regionNames.push_back("j1to3_b0_pT1");
  regionNames.push_back("j1to3_b1_pT0");
  regionNames.push_back("j1to3_b1_pT1");
  regionNames.push_back("j1to3_b2toInf_pT0");
  regionNames.push_back("j1to3_b2toInf_pT1");
  regionNames.push_back("j4toInf_b0_pT0");
  regionNames.push_back("j4toInf_b0_pT1");
  regionNames.push_back("j4toInf_b1_pT0");
  regionNames.push_back("j4toInf_b1_pT1");
  regionNames.push_back("j4toInf_b2toInf_pT0");
  regionNames.push_back("j4toInf_b2toInf_pT1");
  regionNames.push_back("is1El_pT0");
  regionNames.push_back("is1Mu_pT0");
  regionNames.push_back("is1El_pT1");
  regionNames.push_back("is1Mu_pT1");
  regionNames.push_back("diLepZ");
  regionNames.push_back("diBBH_pT0");
  regionNames.push_back("diBBZ_pT0");
  regionNames.push_back("diBBH_pT1");
  regionNames.push_back("diBBZ_pT1");

  //std::vector<int> mt2_bins = { 0 };
  std::vector<int> mt2_bins = { 0, 30 };
  unsigned int regionLength = regionSelection.size();


  for( unsigned int j=0; j<mt2_bins.size(); j++){ 
    for( unsigned int i=0; i < regionLength; i++){// starts at 2 because 0j bins have no "real" second object to make mt2 with
      if( i>15 ) continue;
      if( j==0  ){
  	if ( i<2 ) continue;

  	regionSelection.push_back( regionSelection[i] );
  	regionSelection[i] +=  "  && ( hgg_mt2<30   )";

  	regionNames.push_back( regionNames[i] );
  	regionNames[i] +=  "_mt2_0";
  	//	std::cout << "selection 0 " << regionSelection[i] << std::endl;
      }else{
  	if ( i>(15-2) ) continue;
  	regionSelection[ regionLength + i ] +=  "  && ( hgg_mt2>=30   )";
  	regionNames[ regionLength + i ] +=  "_mt2_30";
  	//	std::cout << "selection 1 " << regionSelection[regionLength + i ] << std::endl;
      }

    }

  }



  // for( unsigned int i=0; i < regionSelection.size(); i++){
  //   std::cout << regionNames[ i ] << " with selection" << regionSelection[i] << std::endl;
  // }



  //BACKGROUNDS
  MT2Analysis<MT2EstimateTree>* qcd      = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_qcd.root", ggDir.c_str()  ), "qcd");
  MT2Analysis<MT2EstimateTree>* diPhoton = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diPhoton.root", ggDir.c_str() ), "diPhoton");
  MT2Analysis<MT2EstimateTree>* gjets    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_gjets.root", ggDir.c_str() ), "gjets");
  MT2Analysis<MT2EstimateTree>* higgs    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir.c_str() ), "higgs");
  
  //DATA
  MT2Analysis<MT2EstimateTree>* data     = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ggDir.c_str() ), "diPhoton_data");

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
  if( doLLbin )       outFileName += "_ll";
  if( doMT2binning )  outFileName += "_mt2";
  if( doSingleLbin )  outFileName += "_onelPT";
  if( doBBbin )       outFileName += "_bb";
  if( mt2_bins.size() > 1 ) outFileName += "_mt2_30";
  
  TFile* dataFile = new TFile(  Form("%s.root",outFileName.c_str()) , "RECREATE" );
  RooWorkspace ws_sig ("ws_sig");
  RooWorkspace ws_bg  ("ws_bg");
  RooWorkspace ws_data("ws_data");

  std::string dataAmountName( Form("%s/dataAmount",  outputdir.c_str() ));
  if( doLLbin )             dataAmountName += "_ll";
  if( doMT2binning )        dataAmountName += "_mt2";
  if( doSingleLbin )        dataAmountName += "_onelPT";
  if( doBBbin )             dataAmountName += "_bb";
  if( mt2_bins.size() > 1 ) dataAmountName += "_mt2_30";
  std::ofstream dataAmount( Form("%s.txt", dataAmountName.c_str() ) );

  std::string isrName( Form("%s/isr",  outputdir.c_str() ));
  if( doLLbin )    isrName += "_ll";
  if( doMT2binning )    isrName += "_mt2";
  if( doSingleLbin )    isrName += "_onelPT";
  if( doBBbin )    isrName += "_bb";
  if( mt2_bins.size() > 1 ) isrName += "_mt2_30";
  std::ofstream isr( Form("%s.txt", isrName.c_str()) );

  std::string sigAmountName( Form("%s/sigAmount_%s",  outputdir.c_str(), sigName.c_str() ));
  if( doLLbin )    sigAmountName += "_ll";
  if( doMT2binning )    sigAmountName += "_mt2";
  if( doSingleLbin )    sigAmountName += "_onelPT";
  if( doBBbin )    sigAmountName += "_bb";
  if( mt2_bins.size() > 1 ) sigAmountName += "_mt2_30";
  std::ofstream sigAmount( Form("%s.txt", sigAmountName.c_str()) );
        
  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this


  // if( doMT2binning )
  //   std::cout << "doing more than 1 MT2 bins, whooop" << std::endl;

  // std::string lowMT2sel  = " && (mt2 < 30 ) ";
  // std::string highMT2sel = " && (mt2 >= 30 ) ";
  // std::string lowLowPtSel  = " && ( (h_pt/h_mass) <   0.4 ) ";
  // std::string highLowPtSel = " && ( (h_pt/h_mass) >=  0.4 ) ";
  // std::string lowHighPtSel  = " && ( (h_pt/h_mass) <   1. ) ";
  // std::string highHighPtSel = " && ( (h_pt/h_mass) >=  1. ) ";

  // // std::string lowMT2sel  = " && (gg_mt2 < 60 ) ";
  // // std::string highMT2sel = " && (gg_mt2 >= 60 ) ";

  // std::string lowMT2name = "loMT2";
  // std::string highMT2name = "hiMT2";
 
  // std::string lowPtname = "loPt";
  // std::string highPtname = "hiPt";



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


  ///////////////////////////////////////////////////////////////////////////////////
  ///// BEGIN LOOP OVER REGIONS /////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  std::set<MT2Region> MT2Regions = higgs->getRegions();

  std::set<MT2Region>::iterator iMT2 = MT2Regions.begin();
  MT2Region thisRegion( (*iMT2) );


  for( unsigned int i=0; i < regionSelection.size(); i++){
    //  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    // int mt2BinsIt = 0;
    // if( !doMT2binning )
    //   mt2BinsIt = 1;

    // for( ; mt2BinsIt<2; mt2BinsIt++){

    string regionSaveName =  regionNames[i];

    //INITIAL SELECTION 
    string treeSel = "((( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )" ;

    //    string treeSel = "((( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 && ( (h_mass>=105 &&h_mass<120  )||( h_mass>=130 &&h_mass<160 )) )" ;
    //string treeSel = "((( ptGamma0) > 40) && ((ptGamma1 )> 25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 &&( (h_mass>=105 &&h_mass<120  )||( h_mass>=130 &&h_mass<160 ) ) )" ;

    //string treeSel = "((( ptGamma0) > 40) && ((ptGamma1 )> 25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )" ;
    
    


    treeSel += " && " + regionSelection[i];
    //string treeSelWideEta = "((( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=2.4  && fabs(etaGamma1)<=2.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )" ;
      
    // string regionSaveName =  thisRegion.getName();
    // if( doMT2binning && !(regionSaveName=="HT0toInf_j4toInf_b2toInf_pT1")   && !(regionSaveName=="HT0toInf_j4toInf_b2toInf_pT0") && !(regionSaveName=="HT0toInf_j1to3_b2toInf_pT1")   && !(regionSaveName=="HT0toInf_j0_b0toInf_pT0")  && !(regionSaveName=="HT0toInf_j0_b0toInf_pT1") ){
    // 	if( mt2BinsIt==0){
    // 	  treeSel += lowMT2sel;
    // 	  regionSaveName += "_"+lowMT2name;
    // 	}else{
    // 	  treeSel += highMT2sel;
    // 	  regionSaveName += "_"+highMT2name;	
    // 	}
    // }else if(doMT2binning && ((regionSaveName=="HT0toInf_j0_b0toInf_pT0")  || (regionSaveName=="HT0toInf_j0_b0toInf_pT1")) ){
    // 	if( (regionSaveName=="HT0toInf_j0_b0toInf_pT0")){
    // 	  if( mt2BinsIt==0){
    // 	    treeSel += lowLowPtSel;
    // 	    regionSaveName += "_"+lowPtname;
    // 	  }else{
    // 	    treeSel += highLowPtSel;
    // 	    regionSaveName += "_"+highPtname;	
    // 	  }
    // 	}else{
    // 	  if( mt2BinsIt==0){
    // 	    treeSel += lowHighPtSel;
    // 	    regionSaveName += "_"+lowPtname;
    // 	  }else{
    // 	    treeSel += highHighPtSel;
    // 	    regionSaveName += "_"+highPtname;	
    // 	  }
    // 	}
    // }

    // std::string isHighPt = regionSaveName.substr( regionSaveName.find("_pT")+1);
    // std::cout << isHighPt << std::endl;
    // std::cout << (isHighPt=="pT1") << std::endl;
    // if( isHighPt == "pT1_loMT2" || isHighPt == "pT1_loPt" || isHighPt == "pT1_hiMT2" || isHighPt == "pT1_hiPt" )
    // 	isHighPt = "pT1";
    // if( isHighPt == "pT0_loMT2" || isHighPt == "pT0_loPt" || isHighPt == "pT0_hiMT2" || isHighPt == "pT0_hiPt" )
    // 	isHighPt = "pT0";


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

    ////////////////////////////////////////////////////////////////////////
    //NOW for the SIGNAL LIST
    ////////////////////////////////////////////////////////////////////////
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

    //  } // loop over mt2 bin(s)
      

    
  }//done loop over signals
 
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
