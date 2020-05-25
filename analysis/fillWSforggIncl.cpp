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


  bool doSMHDonly = 0;

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
    else if( inputName == "SMHD" )
      doSMHDonly = 1;
    else 
      std::cout << "Unknown new physics you are looking for is not here" << std::endl;
  }


  bool  doISRsyst = false;
  //  if( doT2bH )
  doISRsyst = true; //doing it now for all, cause I can

  float lumi2016=cfg.lumi();
  float lumi2017=0.;

  std::string ggDir_2017;
  std::string ggDir_2017_JECup;
  std::string ggDir_2017_JECdn;
  bool doComb = false;
  if( argc>3 ) {
    std::cout << "Doing 2016 and 2017 combination" << std::endl;
    doComb = true;

    std::string configFileName_2017(argv[3]);
    MT2Config cfg_2017( configFileName_2017);
    ggDir_2017 = cfg_2017.getEventYieldDir() + "/diPhotonControlRegion/";
    ggDir_2017_JECup = cfg_2017.getEventYieldDir() + "_JECup/diPhotonControlRegion/";
    ggDir_2017_JECdn = cfg_2017.getEventYieldDir() + "_JECdn/diPhotonControlRegion/";

    lumi2017=cfg_2017.lumi();
  }
 


  unsigned int sigPerPart = 50;
  unsigned int partNr = 1;

  if( argc>4 ) {
    int x= atoi((argv[4]));
    sigPerPart = 5;
    partNr= x;
  }

  std::cout << "working on part nr " << partNr << std::endl;



  std::string ggDir_JECup = cfg.getEventYieldDir() + "_JECup/diPhotonControlRegion/";
  std::string ggDir_JECdn = cfg.getEventYieldDir() + "_JECdn/diPhotonControlRegion/";

  std::cout << "lumi 2016 = " << lumi2016 << " and lumi2017 = " << lumi2017 << std::endl;



  vector<std::string> signalList_T2bH = {"SMS_T2bH_mSbottom250_mLSP1", "SMS_T2bH_mSbottom250_mLSP50", "SMS_T2bH_mSbottom250_mLSP100", "SMS_T2bH_mSbottom300_mLSP1", "SMS_T2bH_mSbottom300_mLSP50", "SMS_T2bH_mSbottom300_mLSP100", "SMS_T2bH_mSbottom300_mLSP150", "SMS_T2bH_mSbottom350_mLSP1", "SMS_T2bH_mSbottom350_mLSP50", "SMS_T2bH_mSbottom350_mLSP100", "SMS_T2bH_mSbottom350_mLSP150", "SMS_T2bH_mSbottom350_mLSP200", "SMS_T2bH_mSbottom400_mLSP1", "SMS_T2bH_mSbottom400_mLSP50", "SMS_T2bH_mSbottom400_mLSP100", "SMS_T2bH_mSbottom400_mLSP150", "SMS_T2bH_mSbottom400_mLSP200", "SMS_T2bH_mSbottom400_mLSP250", "SMS_T2bH_mSbottom450_mLSP1", "SMS_T2bH_mSbottom450_mLSP50", "SMS_T2bH_mSbottom450_mLSP100", "SMS_T2bH_mSbottom450_mLSP150", "SMS_T2bH_mSbottom450_mLSP200", "SMS_T2bH_mSbottom450_mLSP250", "SMS_T2bH_mSbottom450_mLSP300", "SMS_T2bH_mSbottom500_mLSP1", "SMS_T2bH_mSbottom500_mLSP50", "SMS_T2bH_mSbottom500_mLSP100", "SMS_T2bH_mSbottom500_mLSP150", "SMS_T2bH_mSbottom500_mLSP200", "SMS_T2bH_mSbottom500_mLSP250", "SMS_T2bH_mSbottom500_mLSP300", "SMS_T2bH_mSbottom500_mLSP350", "SMS_T2bH_mSbottom550_mLSP1", "SMS_T2bH_mSbottom550_mLSP50", "SMS_T2bH_mSbottom550_mLSP100", "SMS_T2bH_mSbottom550_mLSP150", "SMS_T2bH_mSbottom550_mLSP200", "SMS_T2bH_mSbottom550_mLSP250", "SMS_T2bH_mSbottom550_mLSP300", "SMS_T2bH_mSbottom550_mLSP350", "SMS_T2bH_mSbottom550_mLSP400", "SMS_T2bH_mSbottom600_mLSP1", "SMS_T2bH_mSbottom600_mLSP50", "SMS_T2bH_mSbottom600_mLSP100", "SMS_T2bH_mSbottom600_mLSP150", "SMS_T2bH_mSbottom600_mLSP200", "SMS_T2bH_mSbottom600_mLSP250", "SMS_T2bH_mSbottom600_mLSP300", "SMS_T2bH_mSbottom600_mLSP350", "SMS_T2bH_mSbottom600_mLSP400", "SMS_T2bH_mSbottom600_mLSP450", "SMS_T2bH_mSbottom650_mLSP1", "SMS_T2bH_mSbottom650_mLSP50", "SMS_T2bH_mSbottom650_mLSP100", "SMS_T2bH_mSbottom650_mLSP150", "SMS_T2bH_mSbottom650_mLSP200", "SMS_T2bH_mSbottom650_mLSP250", "SMS_T2bH_mSbottom650_mLSP300", "SMS_T2bH_mSbottom650_mLSP350", "SMS_T2bH_mSbottom650_mLSP400", "SMS_T2bH_mSbottom650_mLSP450", "SMS_T2bH_mSbottom650_mLSP500", "SMS_T2bH_mSbottom700_mLSP1", "SMS_T2bH_mSbottom700_mLSP50", "SMS_T2bH_mSbottom700_mLSP100", "SMS_T2bH_mSbottom700_mLSP150", "SMS_T2bH_mSbottom700_mLSP200", "SMS_T2bH_mSbottom700_mLSP250", "SMS_T2bH_mSbottom700_mLSP300", "SMS_T2bH_mSbottom700_mLSP350", "SMS_T2bH_mSbottom700_mLSP400", "SMS_T2bH_mSbottom700_mLSP450", "SMS_T2bH_mSbottom700_mLSP500", "SMS_T2bH_mSbottom700_mLSP550", "SMS_T2bH_mSbottom750_mLSP1", "SMS_T2bH_mSbottom750_mLSP50", "SMS_T2bH_mSbottom750_mLSP100", "SMS_T2bH_mSbottom750_mLSP150", "SMS_T2bH_mSbottom750_mLSP200", "SMS_T2bH_mSbottom750_mLSP250", "SMS_T2bH_mSbottom750_mLSP300", "SMS_T2bH_mSbottom750_mLSP350", "SMS_T2bH_mSbottom750_mLSP400", "SMS_T2bH_mSbottom750_mLSP450", "SMS_T2bH_mSbottom750_mLSP500", "SMS_T2bH_mSbottom750_mLSP550", "SMS_T2bH_mSbottom750_mLSP600", "SMS_T2bH_mSbottom800_mLSP1", "SMS_T2bH_mSbottom800_mLSP50", "SMS_T2bH_mSbottom800_mLSP100", "SMS_T2bH_mSbottom800_mLSP150", "SMS_T2bH_mSbottom800_mLSP200", "SMS_T2bH_mSbottom800_mLSP250", "SMS_T2bH_mSbottom800_mLSP300", "SMS_T2bH_mSbottom800_mLSP350", "SMS_T2bH_mSbottom800_mLSP400", "SMS_T2bH_mSbottom800_mLSP450", "SMS_T2bH_mSbottom800_mLSP500", "SMS_T2bH_mSbottom800_mLSP550", "SMS_T2bH_mSbottom800_mLSP600", "SMS_T2bH_mSbottom800_mLSP650", "SMS_T2bH_mSbottom850_mLSP1", "SMS_T2bH_mSbottom850_mLSP50", "SMS_T2bH_mSbottom850_mLSP100", "SMS_T2bH_mSbottom850_mLSP150", "SMS_T2bH_mSbottom850_mLSP200", "SMS_T2bH_mSbottom850_mLSP250", "SMS_T2bH_mSbottom850_mLSP300", "SMS_T2bH_mSbottom850_mLSP350", "SMS_T2bH_mSbottom850_mLSP400", "SMS_T2bH_mSbottom850_mLSP450", "SMS_T2bH_mSbottom850_mLSP500", "SMS_T2bH_mSbottom850_mLSP550", "SMS_T2bH_mSbottom850_mLSP600", "SMS_T2bH_mSbottom850_mLSP650", "SMS_T2bH_mSbottom850_mLSP700", "SMS_T2bH_mSbottom900_mLSP1", "SMS_T2bH_mSbottom900_mLSP50", "SMS_T2bH_mSbottom900_mLSP100", "SMS_T2bH_mSbottom900_mLSP150", "SMS_T2bH_mSbottom900_mLSP200", "SMS_T2bH_mSbottom900_mLSP250", "SMS_T2bH_mSbottom900_mLSP300", "SMS_T2bH_mSbottom900_mLSP350", "SMS_T2bH_mSbottom900_mLSP400", "SMS_T2bH_mSbottom900_mLSP450", "SMS_T2bH_mSbottom900_mLSP500", "SMS_T2bH_mSbottom900_mLSP550", "SMS_T2bH_mSbottom900_mLSP600", "SMS_T2bH_mSbottom900_mLSP650", "SMS_T2bH_mSbottom900_mLSP700", "SMS_T2bH_mSbottom900_mLSP750" };


   vector<std::string> signalList_WH = { "SMS_TChiWH_HToGG_m127_m1","SMS_TChiWH_HToGG_m150_m1","SMS_TChiWH_HToGG_m150_m24","SMS_TChiWH_HToGG_m175_m1","SMS_TChiWH_HToGG_m175_m25","SMS_TChiWH_HToGG_m175_m49","SMS_TChiWH_HToGG_m200_m1","SMS_TChiWH_HToGG_m200_m25","SMS_TChiWH_HToGG_m200_m50","SMS_TChiWH_HToGG_m200_m74","SMS_TChiWH_HToGG_m225_m1","SMS_TChiWH_HToGG_m225_m25","SMS_TChiWH_HToGG_m225_m50","SMS_TChiWH_HToGG_m225_m75","SMS_TChiWH_HToGG_m225_m99","SMS_TChiWH_HToGG_m250_m1","SMS_TChiWH_HToGG_m250_m25","SMS_TChiWH_HToGG_m250_m50","SMS_TChiWH_HToGG_m250_m75","SMS_TChiWH_HToGG_m250_m100","SMS_TChiWH_HToGG_m250_m124","SMS_TChiWH_HToGG_m275_m1","SMS_TChiWH_HToGG_m275_m25","SMS_TChiWH_HToGG_m275_m50","SMS_TChiWH_HToGG_m275_m75","SMS_TChiWH_HToGG_m275_m100","SMS_TChiWH_HToGG_m275_m125","SMS_TChiWH_HToGG_m275_m149","SMS_TChiWH_HToGG_m300_m1","SMS_TChiWH_HToGG_m300_m25","SMS_TChiWH_HToGG_m300_m50","SMS_TChiWH_HToGG_m300_m75","SMS_TChiWH_HToGG_m300_m100","SMS_TChiWH_HToGG_m300_m125","SMS_TChiWH_HToGG_m300_m150","SMS_TChiWH_HToGG_m300_m174","SMS_TChiWH_HToGG_m325_m1","SMS_TChiWH_HToGG_m325_m25","SMS_TChiWH_HToGG_m325_m50","SMS_TChiWH_HToGG_m325_m75","SMS_TChiWH_HToGG_m325_m100","SMS_TChiWH_HToGG_m325_m125","SMS_TChiWH_HToGG_m325_m150","SMS_TChiWH_HToGG_m325_m175","SMS_TChiWH_HToGG_m325_m199","SMS_TChiWH_HToGG_m350_m1","SMS_TChiWH_HToGG_m350_m25","SMS_TChiWH_HToGG_m350_m50","SMS_TChiWH_HToGG_m350_m75","SMS_TChiWH_HToGG_m350_m100","SMS_TChiWH_HToGG_m350_m125","SMS_TChiWH_HToGG_m350_m150","SMS_TChiWH_HToGG_m350_m175","SMS_TChiWH_HToGG_m350_m200","SMS_TChiWH_HToGG_m350_m224","SMS_TChiWH_HToGG_m375_m1","SMS_TChiWH_HToGG_m375_m25","SMS_TChiWH_HToGG_m375_m50","SMS_TChiWH_HToGG_m375_m75","SMS_TChiWH_HToGG_m375_m100","SMS_TChiWH_HToGG_m375_m125","SMS_TChiWH_HToGG_m375_m150","SMS_TChiWH_HToGG_m375_m175","SMS_TChiWH_HToGG_m375_m200","SMS_TChiWH_HToGG_m375_m225","SMS_TChiWH_HToGG_m375_m249","SMS_TChiWH_HToGG_m400_m1","SMS_TChiWH_HToGG_m400_m25","SMS_TChiWH_HToGG_m400_m50","SMS_TChiWH_HToGG_m400_m75","SMS_TChiWH_HToGG_m400_m100","SMS_TChiWH_HToGG_m400_m125","SMS_TChiWH_HToGG_m400_m150","SMS_TChiWH_HToGG_m400_m175","SMS_TChiWH_HToGG_m400_m200","SMS_TChiWH_HToGG_m400_m225","SMS_TChiWH_HToGG_m400_m250","SMS_TChiWH_HToGG_m400_m274"};
//,"SMS_TChiWH_HToGG_425_m1","SMS_TChiWH_HToGG_425_25","SMS_TChiWH_HToGG_425_50","SMS_TChiWH_HToGG_425_75","SMS_TChiWH_HToGG_425_m100","SMS_TChiWH_HToGG_425_m125","SMS_TChiWH_HToGG_425_m150","SMS_TChiWH_HToGG_425_m175","SMS_TChiWH_HToGG_425_200","SMS_TChiWH_HToGG_425_225","SMS_TChiWH_HToGG_425_250","SMS_TChiWH_HToGG_425_275","SMS_TChiWH_HToGG_425_299"};//,"SMS_TChiWH_HToGG_450_m1","SMS_TChiWH_HToGG_450_25","SMS_TChiWH_HToGG_450_50","SMS_TChiWH_HToGG_450_75","SMS_TChiWH_HToGG_450_m100","SMS_TChiWH_HToGG_450_m125","SMS_TChiWH_HToGG_450_m150","SMS_TChiWH_HToGG_450_m175","SMS_TChiWH_HToGG_450_200","SMS_TChiWH_HToGG_450_225","SMS_TChiWH_HToGG_450_250","SMS_TChiWH_HToGG_450_275","SMS_TChiWH_HToGG_450_300","SMS_TChiWH_HToGG_475_m1","SMS_TChiWH_HToGG_475_25","SMS_TChiWH_HToGG_475_50","SMS_TChiWH_HToGG_475_75","SMS_TChiWH_HToGG_475_m100","SMS_TChiWH_HToGG_475_m125","SMS_TChiWH_HToGG_475_m150","SMS_TChiWH_HToGG_475_m175","SMS_TChiWH_HToGG_475_200","SMS_TChiWH_HToGG_475_225","SMS_TChiWH_HToGG_475_250","SMS_TChiWH_HToGG_475_275","SMS_TChiWH_HToGG_475_300","SMS_TChiWH_HToGG_500_m1","SMS_TChiWH_HToGG_500_25","SMS_TChiWH_HToGG_500_50","SMS_TChiWH_HToGG_500_75","SMS_TChiWH_HToGG_500_m100","SMS_TChiWH_HToGG_500_m125","SMS_TChiWH_HToGG_500_m150","SMS_TChiWH_HToGG_500_m175","SMS_TChiWH_HToGG_500_200","SMS_TChiWH_HToGG_500_225","SMS_TChiWH_HToGG_500_250","SMS_TChiWH_HToGG_500_275","SMS_TChiWH_HToGG_500_300","SMS_TChiWH_HToGG_525_m1","SMS_TChiWH_HToGG_525_25","SMS_TChiWH_HToGG_525_50","SMS_TChiWH_HToGG_525_75","SMS_TChiWH_HToGG_525_m100","SMS_TChiWH_HToGG_525_m125","SMS_TChiWH_HToGG_525_m150","SMS_TChiWH_HToGG_525_m175","SMS_TChiWH_HToGG_525_200","SMS_TChiWH_HToGG_525_225","SMS_TChiWH_HToGG_525_250","SMS_TChiWH_HToGG_525_275","SMS_TChiWH_HToGG_525_300","SMS_TChiWH_HToGG_550_m1","SMS_TChiWH_HToGG_550_25","SMS_TChiWH_HToGG_550_50","SMS_TChiWH_HToGG_550_75","SMS_TChiWH_HToGG_550_m100","SMS_TChiWH_HToGG_550_m125","SMS_TChiWH_HToGG_550_m150","SMS_TChiWH_HToGG_550_m175","SMS_TChiWH_HToGG_550_200","SMS_TChiWH_HToGG_550_225","SMS_TChiWH_HToGG_550_250","SMS_TChiWH_HToGG_550_275","SMS_TChiWH_HToGG_550_300","SMS_TChiWH_HToGG_575_m1","SMS_TChiWH_HToGG_575_25","SMS_TChiWH_HToGG_575_50","SMS_TChiWH_HToGG_575_75","SMS_TChiWH_HToGG_575_m100","SMS_TChiWH_HToGG_575_m125","SMS_TChiWH_HToGG_575_m150","SMS_TChiWH_HToGG_575_m175","SMS_TChiWH_HToGG_575_200","SMS_TChiWH_HToGG_575_225","SMS_TChiWH_HToGG_575_250","SMS_TChiWH_HToGG_575_275","SMS_TChiWH_HToGG_575_300","SMS_TChiWH_HToGG_600_m1","SMS_TChiWH_HToGG_600_25","SMS_TChiWH_HToGG_600_50","SMS_TChiWH_HToGG_600_75","SMS_TChiWH_HToGG_600_m100","SMS_TChiWH_HToGG_600_m125","SMS_TChiWH_HToGG_600_m150","SMS_TChiWH_HToGG_600_m175","SMS_TChiWH_HToGG_600_200","SMS_TChiWH_HToGG_600_225","SMS_TChiWH_HToGG_600_250","SMS_TChiWH_HToGG_600_275","SMS_TChiWH_HToGG_600_300","SMS_TChiWH_HToGG_625_m1","SMS_TChiWH_HToGG_625_25","SMS_TChiWH_HToGG_625_50","SMS_TChiWH_HToGG_625_75","SMS_TChiWH_HToGG_625_m100","SMS_TChiWH_HToGG_625_m125","SMS_TChiWH_HToGG_625_m150","SMS_TChiWH_HToGG_625_m175","SMS_TChiWH_HToGG_625_200","SMS_TChiWH_HToGG_625_225","SMS_TChiWH_HToGG_625_250","SMS_TChiWH_HToGG_625_275","SMS_TChiWH_HToGG_625_300","SMS_TChiWH_HToGG_650_m1","SMS_TChiWH_HToGG_650_25","SMS_TChiWH_HToGG_650_50","SMS_TChiWH_HToGG_650_75","SMS_TChiWH_HToGG_650_m100","SMS_TChiWH_HToGG_650_m125","SMS_TChiWH_HToGG_650_m150","SMS_TChiWH_HToGG_650_m175","SMS_TChiWH_HToGG_650_200","SMS_TChiWH_HToGG_650_225","SMS_TChiWH_HToGG_650_250","SMS_TChiWH_HToGG_650_275","SMS_TChiWH_HToGG_650_300","SMS_TChiWH_HToGG_675_m1","SMS_TChiWH_HToGG_675_25","SMS_TChiWH_HToGG_675_50","SMS_TChiWH_HToGG_675_75","SMS_TChiWH_HToGG_675_m100","SMS_TChiWH_HToGG_675_m125","SMS_TChiWH_HToGG_675_m150","SMS_TChiWH_HToGG_675_m175","SMS_TChiWH_HToGG_675_200","SMS_TChiWH_HToGG_675_225","SMS_TChiWH_HToGG_675_250","SMS_TChiWH_HToGG_675_275","SMS_TChiWH_HToGG_675_300","SMS_TChiWH_HToGG_700_m1","SMS_TChiWH_HToGG_700_25","SMS_TChiWH_HToGG_700_50","SMS_TChiWH_HToGG_700_75","SMS_TChiWH_HToGG_700_m100","SMS_TChiWH_HToGG_700_m125","SMS_TChiWH_HToGG_700_m150","SMS_TChiWH_HToGG_700_m175","SMS_TChiWH_HToGG_700_200","SMS_TChiWH_HToGG_700_225","SMS_TChiWH_HToGG_700_250","SMS_TChiWH_HToGG_700_275","SMS_TChiWH_HToGG_700_300"};//

   vector<std::string> signalList_HZ = { "SMS_TChiHZ_HToGG_m127","SMS_TChiHZ_HToGG_m150", "SMS_TChiHZ_HToGG_m175", "SMS_TChiHZ_HToGG_m200", "SMS_TChiHZ_HToGG_m225", "SMS_TChiHZ_HToGG_m250", "SMS_TChiHZ_HToGG_m275", "SMS_TChiHZ_HToGG_m300", "SMS_TChiHZ_HToGG_m325", "SMS_TChiHZ_HToGG_m350", "SMS_TChiHZ_HToGG_m375", "SMS_TChiHZ_HToGG_m400", "SMS_TChiHZ_HToGG_m425", "SMS_TChiHZ_HToGG_m450", "SMS_TChiHZ_HToGG_m475", "SMS_TChiHZ_HToGG_m500", "SMS_TChiHZ_HToGG_m525"};
   //, "SMS_TChiHZ_HToGG_m550", "SMS_TChiHZ_HToGG_m575", "SMS_TChiHZ_HToGG_m600", "SMS_TChiHZ_HToGG_m625", "SMS_TChiHZ_HToGG_m650", "SMS_TChiHZ_HToGG_m675", "SMS_TChiHZ_HToGG_m700", "SMS_TChiHZ_HToGG_m725", "SMS_TChiHZ_HToGG_m750", "SMS_TChiHZ_HToGG_m775", "SMS_TChiHZ_HToGG_m800", "SMS_TChiHZ_HToGG_m825", "SMS_TChiHZ_HToGG_m850", "SMS_TChiHZ_HToGG_m875", "SMS_TChiHZ_HToGG_m900", "SMS_TChiHZ_HToGG_m925", "SMS_TChiHZ_HToGG_m950", "SMS_TChiHZ_HToGG_m975", "SMS_TChiHZ_HToGG_m1000" };

   vector<std::string> signalList_HH_HH0p25 = { "SMS_TChiHZ_HToGG_m127_HH0p25", "SMS_TChiHZ_HToGG_m150_HH0p25", "SMS_TChiHZ_HToGG_m175_HH0p25", "SMS_TChiHZ_HToGG_m200_HH0p25", "SMS_TChiHZ_HToGG_m225_HH0p25", "SMS_TChiHZ_HToGG_m250_HH0p25", "SMS_TChiHZ_HToGG_m275_HH0p25", "SMS_TChiHZ_HToGG_m300_HH0p25", "SMS_TChiHZ_HToGG_m325_HH0p25", "SMS_TChiHZ_HToGG_m350_HH0p25", "SMS_TChiHZ_HToGG_m375_HH0p25", "SMS_TChiHZ_HToGG_m400_HH0p25", "SMS_TChiHZ_HToGG_m425_HH0p25", "SMS_TChiHZ_HToGG_m450_HH0p25", "SMS_TChiHZ_HToGG_m475_HH0p25", "SMS_TChiHZ_HToGG_m500_HH0p25", "SMS_TChiHZ_HToGG_m525_HH0p25"}; 
   //, "SMS_TChiHZ_HToGG_m550_HH0p25", "SMS_TChiHZ_HToGG_m575_HH0p25", "SMS_TChiHZ_HToGG_m600_HH0p25", "SMS_TChiHZ_HToGG_m625_HH0p25", "SMS_TChiHZ_HToGG_m650_HH0p25", "SMS_TChiHZ_HToGG_m675_HH0p25", "SMS_TChiHZ_HToGG_m700_HH0p25", "SMS_TChiHZ_HToGG_m725_HH0p25", "SMS_TChiHZ_HToGG_m750_HH0p25", "SMS_TChiHZ_HToGG_m775_HH0p25", "SMS_TChiHZ_HToGG_m800_HH0p25", "SMS_TChiHZ_HToGG_m825_HH0p25", "SMS_TChiHZ_HToGG_m850_HH0p25", "SMS_TChiHZ_HToGG_m875_HH0p25", "SMS_TChiHZ_HToGG_m900_HH0p25", "SMS_TChiHZ_HToGG_m925_HH0p25", "SMS_TChiHZ_HToGG_m950_HH0p25", "SMS_TChiHZ_HToGG_m975_HH0p25", "SMS_TChiHZ_HToGG_m1000_HH0p25" };


   //   vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127"};
   vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127", "SMS_TChiHH_HToGG_m150", "SMS_TChiHH_HToGG_m175", "SMS_TChiHH_HToGG_m200", "SMS_TChiHH_HToGG_m225", "SMS_TChiHH_HToGG_m250", "SMS_TChiHH_HToGG_m275", "SMS_TChiHH_HToGG_m300", "SMS_TChiHH_HToGG_m325", "SMS_TChiHH_HToGG_m350", "SMS_TChiHH_HToGG_m375", "SMS_TChiHH_HToGG_m400", "SMS_TChiHH_HToGG_m425", "SMS_TChiHH_HToGG_m450", "SMS_TChiHH_HToGG_m475", "SMS_TChiHH_HToGG_m500", "SMS_TChiHH_HToGG_m525"};
   //, "SMS_TChiHH_HToGG_m550", "SMS_TChiHH_HToGG_m575", "SMS_TChiHH_HToGG_m600", "SMS_TChiHH_HToGG_m625", "SMS_TChiHH_HToGG_m650", "SMS_TChiHH_HToGG_m675", "SMS_TChiHH_HToGG_m700", "SMS_TChiHH_HToGG_m725", "SMS_TChiHH_HToGG_m750", "SMS_TChiHH_HToGG_m775", "SMS_TChiHH_HToGG_m800", "SMS_TChiHH_HToGG_m825", "SMS_TChiHH_HToGG_m850", "SMS_TChiHH_HToGG_m875", "SMS_TChiHH_HToGG_m900", "SMS_TChiHH_HToGG_m925", "SMS_TChiHH_HToGG_m950", "SMS_TChiHH_HToGG_m975", "SMS_TChiHH_HToGG_m1000" };



  std::vector<std::string> regionSelection_ph;
  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 0 && ((h_pt/h_mass) <  0.6)) " );
  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 0 && ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0)  ) " );
  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 0 && ((h_pt/h_mass) >= 1.0)) " );
  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 1 && ((h_pt/h_mass) <  0.6)) " );
  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 1 && ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0)  ) " );
  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 1 && ((h_pt/h_mass) >= 1.0)) " );
  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets >= 2 && ((h_pt/h_mass) <  0.6)) " );
  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets >= 2&& ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0)  ) " );
  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets >= 2 && ((h_pt/h_mass) >= 1.0)) " );

  regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 0 && ((h_pt/h_mass) <  0.6)) " );
  regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 0 && ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0)  ) " );
  regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 0 && ((h_pt/h_mass) >= 1.0)) " );
  regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 1 && ((h_pt/h_mass) <  0.6)) " );
  regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 1 && ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0)  ) " );
  regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 1 && ((h_pt/h_mass) >= 1.0)) " );
  regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets >= 2 && ((h_pt/h_mass) <  0.6)) " );
  regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets >= 2 && ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0)  ) " );
  regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets >= 2 && ((h_pt/h_mass) >= 1.0)) " );

  for( unsigned int i=0; i < regionSelection_ph.size(); i++){
    regionSelection_ph[i] +=  "  && ( (!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )";
  }
  
  regionSelection_ph.push_back(" ( is1El) && ((h_pt/h_mass) <  0.6)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  regionSelection_ph.push_back(" ( is1Mu) && ((h_pt/h_mass) <  0.6)   && ( (!isDiBH && !isDiBZ) && !(isDiLepZ) )"  );
  regionSelection_ph.push_back(" ( is1El) && ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0)  && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  regionSelection_ph.push_back(" ( is1Mu) && ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0)  && ( (!isDiBH && !isDiBZ) && !(isDiLepZ) )"  );
  regionSelection_ph.push_back(" ( is1El) && ((h_pt/h_mass) >= 1.0)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  regionSelection_ph.push_back(" ( is1Mu) && ((h_pt/h_mass) >= 1.0)   && ( (!isDiBH && !isDiBZ) && !(isDiLepZ) )"  );
  regionSelection_ph.push_back(" ( isDiBH ) && ((h_pt/h_mass) <  0.6) && ( (!isDiBZ) && !(isDiLepZ) )" );
  regionSelection_ph.push_back(" ( isDiBZ ) && ((h_pt/h_mass) <  0.6) && ( (!isDiBH) && !(isDiLepZ) )" );
  regionSelection_ph.push_back(" ( isDiBH ) && ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0) && ( (!isDiBZ) && !(isDiLepZ) )" );
  regionSelection_ph.push_back(" ( isDiBZ ) && ((h_pt/h_mass) >=  0.6) && ((h_pt/h_mass) <  1.0) && ( (!isDiBH) && !(isDiLepZ) )" );
  regionSelection_ph.push_back(" ( isDiBH ) && ((h_pt/h_mass) >= 1.0) && ( (!isDiBZ) && !(isDiLepZ) )" );
  regionSelection_ph.push_back(" ( isDiBZ ) && ((h_pt/h_mass) >= 1.0) && ( (!isDiBH) && !(isDiLepZ) )" );

  std::vector<std::string> regionNames_ph;
  regionNames_ph.push_back("j1to3_b0_pT0");
  regionNames_ph.push_back("j1to3_b0_pT1");
  regionNames_ph.push_back("j1to3_b0_pT2");
  regionNames_ph.push_back("j1to3_b1_pT0");
  regionNames_ph.push_back("j1to3_b1_pT1");
  regionNames_ph.push_back("j1to3_b1_pT2");
  regionNames_ph.push_back("j1to3_b2toInf_pT0");
  regionNames_ph.push_back("j1to3_b2toInf_pT1");
  regionNames_ph.push_back("j1to3_b2toInf_pT2");
  regionNames_ph.push_back("j4toInf_b0_pT0");
  regionNames_ph.push_back("j4toInf_b0_pT1");
  regionNames_ph.push_back("j4toInf_b0_pT2");
  regionNames_ph.push_back("j4toInf_b1_pT0");
  regionNames_ph.push_back("j4toInf_b1_pT1");
  regionNames_ph.push_back("j4toInf_b1_pT2");
  regionNames_ph.push_back("j4toInf_b2toInf_pT0");
  regionNames_ph.push_back("j4toInf_b2toInf_pT1");
  regionNames_ph.push_back("j4toInf_b2toInf_pT2");
  regionNames_ph.push_back("is1El_pT0");
  regionNames_ph.push_back("is1Mu_pT0");
  regionNames_ph.push_back("is1El_pT1");
  regionNames_ph.push_back("is1Mu_pT1");
  regionNames_ph.push_back("is1El_pT2");
  regionNames_ph.push_back("is1Mu_pT2");
  regionNames_ph.push_back("diBBH_pT0");
  regionNames_ph.push_back("diBBZ_pT0");
  regionNames_ph.push_back("diBBH_pT1");
  regionNames_ph.push_back("diBBZ_pT1");
  regionNames_ph.push_back("diBBH_pT2");
  regionNames_ph.push_back("diBBZ_pT2");
 

  std::vector<std::string> mt2_bins = { "0", "30" };
  unsigned int regionLength = regionSelection_ph.size();

  std::vector<std::string> mt2_sel;
  mt2_sel.push_back( "(hgg_mt2<30)" );
  mt2_sel.push_back( "(hgg_mt2>=30)" );

  std::vector<std::string> mt2_genMET_sel;
  mt2_genMET_sel.push_back( "(hgg_mt2_genMET<30)" );
  mt2_genMET_sel.push_back( "(hgg_mt2_genMET>=30)" );

  std::cout << "defined the selections " << std::endl;

  std::vector<std::string> regionNames;
  std::vector<std::string> regionSelection;
  std::vector<std::string> regionSelection_genMET;

  for( unsigned int j=0; j<mt2_bins.size(); j++){ 
    for( unsigned int i=0; i < regionLength; i++){
      regionSelection.push_back( regionSelection_ph[i] );
      regionSelection[ i + j*regionLength ] +=  "  && " +  mt2_sel[j] ;

      regionSelection_genMET.push_back( regionSelection_ph[i] );
      regionSelection_genMET[ i + j*regionLength ] +=  "  && " +  mt2_genMET_sel[j] ;

      regionNames.push_back( regionNames_ph[i] );
      regionNames[i+j*regionLength ] = ( regionNames_ph[i] + "_mt2_" + mt2_bins[j]);
    }
  }

  
  //and now push back the non mt2 split bins 
  regionNames.push_back("j0_b0toInf_pT0");
  regionNames.push_back("j0_b0toInf_pT1");
  regionNames.push_back("j0_b0toInf_pT2");
  regionNames.push_back("diLepZ");

  regionSelection.push_back("(gg_nJets==0 && gg_nBJets>=0 && ((h_pt/h_mass) <  0.6)) && ((!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ))" );
  regionSelection.push_back("(gg_nJets==0 && gg_nBJets>=0 && ((h_pt/h_mass) >= 0.6) && ((h_pt/h_mass) <  1.0)) && ((!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ))" );
  regionSelection.push_back("(gg_nJets==0 && gg_nBJets>=0 && ((h_pt/h_mass) >= 1.0)) && ((!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ))" );
  regionSelection.push_back("( isDiLepZ ) " );


  regionSelection_genMET.push_back("(gg_nJets==0 && gg_nBJets>=0 && ((h_pt/h_mass) <  0.6)) && ((!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ))" );
  regionSelection_genMET.push_back("(gg_nJets==0 && gg_nBJets>=0 && ((h_pt/h_mass) >= 0.6) && ((h_pt/h_mass) <  1.0)) && ((!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ))" );
  regionSelection_genMET.push_back("(gg_nJets==0 && gg_nBJets>=0 && ((h_pt/h_mass) >= 1.0)) && ((!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ))" );
  regionSelection_genMET.push_back("( isDiLepZ ) " );




  //2016 data and MC
  MT2Analysis<MT2EstimateTree>* data     = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ggDir.c_str() ), "diPhoton_data");
  // MT2Analysis<MT2EstimateTree>* qcd      = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_qcd.root", ggDir.c_str()  ), "qcd");
  // MT2Analysis<MT2EstimateTree>* diPhoton = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diPhoton.root", ggDir.c_str() ), "diPhoton");
  // MT2Analysis<MT2EstimateTree>* gjets    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_gjets.root", ggDir.c_str() ), "gjets");
  MT2Analysis<MT2EstimateTree>* higgs    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir.c_str() ), "higgs");
  // MT2Analysis<MT2EstimateTree>* diH      = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diH.root", ggDir.c_str()  ), "diH");
  MT2Analysis<MT2EstimateTree>* higgs_JECup    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir_JECup.c_str() ), "higgs");
  MT2Analysis<MT2EstimateTree>* higgs_JECdn    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir_JECdn.c_str() ), "higgs");
  


  //2017 data and MC  
  MT2Analysis<MT2EstimateTree>* data_2017;
  // MT2Analysis<MT2EstimateTree>* qcd_2017;     
  // //  MT2Analysis<MT2EstimateTree>* diH_2017;     
  // MT2Analysis<MT2EstimateTree>* diPhoton_2017;
  // MT2Analysis<MT2EstimateTree>* gjets_2017;
  MT2Analysis<MT2EstimateTree>* higgs_2017;  

  MT2Analysis<MT2EstimateTree>* higgs_2017_JECup;  
  MT2Analysis<MT2EstimateTree>* higgs_2017_JECdn;  
  

  if( doComb ) {
    data_2017     = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ggDir_2017.c_str() ), "diPhoton_data");
    // qcd_2017      = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_qcd.root", ggDir_2017.c_str()  ), "qcd");
    // //    diH_2017      = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diH.root", ggDir_2017.c_str()  ), "diH");
    // diPhoton_2017 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diPhoton.root", ggDir_2017.c_str() ), "diPhoton");
    // gjets_2017    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_gjets.root", ggDir_2017.c_str() ), "gjets");
    higgs_2017    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir_2017.c_str() ), "higgs");

    higgs_2017_JECup  = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir_2017_JECup.c_str() ), "higgs");
    higgs_2017_JECdn  = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir_2017_JECdn.c_str() ), "higgs");
  
  }

  //SIGNAL
  string sigName;
  vector<std::string> signalList_tmp;

  if( doT2bH ){        signalList_tmp = signalList_T2bH;     sigName = "SMS_T2bH_mSbottom";
  }else if( doHZ ){    signalList_tmp = signalList_HZ;       sigName = "SMS_TChiHZ_HToGG";
  }else if( doHH ){    signalList_tmp = signalList_HH;       sigName = "SMS_TChiHH_HToGG";     
  }else if( doWH ){    signalList_tmp = signalList_WH;       sigName = "SMS_TChiWH_HToGG";
  }else if( doSMHDonly ){    signalList_tmp = {};       sigName = "SMHDonly";     }


  vector<std::string> signalList;

  vector< MT2Analysis<MT2EstimateTree>* > signals;

  vector< MT2Analysis<MT2EstimateTree>* > signals_JECup;
  vector< MT2Analysis<MT2EstimateTree>* > signals_JECdn;


  for( unsigned int i=0; i<signalList_tmp.size(); i++){

    if( (i>= sigPerPart * partNr) || (i < ((partNr-1)*sigPerPart ))   )
      continue;

    signalList.push_back( signalList_tmp[i] );

    if( doT2bH ){
      MT2Analysis<MT2EstimateTree>* sigTemp = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_T2bH.root", ggDir.c_str() ), signalList_tmp[i]  );
      signals.push_back( sigTemp );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECup = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_T2bH.root", ggDir_JECup.c_str() ), signalList_tmp[i]  );
      signals_JECup.push_back( sigTemp_JECup );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECdn = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_T2bH.root", ggDir_JECdn.c_str() ), signalList_tmp[i]  );
      signals_JECdn.push_back( sigTemp_JECdn );
    }
    else if( doHH ){
      MT2Analysis<MT2EstimateTree>* sigTemp = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HH.root", ggDir.c_str() ), signalList_tmp[i]  );
      signals.push_back( sigTemp );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECup = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HH.root", ggDir_JECup.c_str() ), signalList_tmp[i]  );
      signals_JECup.push_back( sigTemp_JECup );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECdn = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HH.root", ggDir_JECdn.c_str() ), signalList_tmp[i]  );
      signals_JECdn.push_back( sigTemp_JECdn );
    }
    else if( doHZ ){
      MT2Analysis<MT2EstimateTree>* sigTemp = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir.c_str() ), signalList_tmp[i]  );
      signals.push_back( sigTemp );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECup = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_JECup.c_str() ), signalList_tmp[i]  );
      signals_JECup.push_back( sigTemp_JECup );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECdn = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_JECdn.c_str() ), signalList_tmp[i]  );
      signals_JECdn.push_back( sigTemp_JECdn );
    }
    else{
      MT2Analysis<MT2EstimateTree>* sigTemp = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_WH.root", ggDir.c_str() ), signalList_tmp[i]  );
      signals.push_back( sigTemp );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECup = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_WH.root", ggDir_JECup.c_str() ), signalList_tmp[i]  );
      signals_JECup.push_back( sigTemp_JECup );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECdn = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_WH.root", ggDir_JECdn.c_str() ), signalList_tmp[i]  );
      signals_JECdn.push_back( sigTemp_JECdn );
    }

  }
 
  vector< MT2Analysis<MT2EstimateTree>* > signals_HH_HH0p25;
  vector< MT2Analysis<MT2EstimateTree>* > signals_HH_HH0p25_JECup;
  vector< MT2Analysis<MT2EstimateTree>* > signals_HH_HH0p25_JECdn;

  if( doHZ )
    for( unsigned int i=0; i<signalList_HH_HH0p25.size(); i++){

      if( (i>= sigPerPart * partNr) || (i < ((partNr-1)*sigPerPart ))   )
	continue;

      MT2Analysis<MT2EstimateTree>* sigTemp_HH_HH0p25 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir.c_str() ), signalList_HH_HH0p25[i]  );
      signals_HH_HH0p25.push_back( sigTemp_HH_HH0p25 );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECup_HH_HH0p25 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_JECup.c_str() ), signalList_HH_HH0p25[i]  );
      signals_HH_HH0p25_JECup.push_back( sigTemp_JECup_HH_HH0p25 );
      MT2Analysis<MT2EstimateTree>* sigTemp_JECdn_HH_HH0p25 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_JECdn.c_str() ), signalList_HH_HH0p25[i]  );
      signals_HH_HH0p25_JECdn.push_back( sigTemp_JECdn_HH_HH0p25 );
    }

  // 2017 version
  vector< MT2Analysis<MT2EstimateTree>* > signals_2017;
  vector< MT2Analysis<MT2EstimateTree>* > signals_2017_JECup;
  vector< MT2Analysis<MT2EstimateTree>* > signals_2017_JECdn;
  if( doComb ) 

    for( unsigned int i=0; i<signalList_tmp.size(); i++){

      if( (i>= sigPerPart * partNr) || (i < ((partNr-1)*sigPerPart ))   )
	continue;

      if( doT2bH ){
	MT2Analysis<MT2EstimateTree>* sigTemp_2017 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_T2bH.root", ggDir_2017.c_str() ), signalList_tmp[i]  );
	signals_2017.push_back( sigTemp_2017 );
	MT2Analysis<MT2EstimateTree>* sigTemp_2017_JECup = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_T2bH.root", ggDir_2017_JECup.c_str() ), signalList_tmp[i]  );
	signals_2017_JECup.push_back( sigTemp_2017_JECup );
	MT2Analysis<MT2EstimateTree>* sigTemp_2017_JECdn = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_T2bH.root", ggDir_2017_JECdn.c_str() ), signalList_tmp[i]  );
	signals_2017_JECdn.push_back( sigTemp_2017_JECdn );

      }
      else if (doHH){
	MT2Analysis<MT2EstimateTree>* sigTemp_2017 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HH.root", ggDir_2017.c_str() ), signalList_tmp[i]  );
	signals_2017.push_back( sigTemp_2017 );
	MT2Analysis<MT2EstimateTree>* sigTemp_2017_JECup = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HH.root", ggDir_2017_JECup.c_str() ), signalList_tmp[i]  );
	signals_2017_JECup.push_back( sigTemp_2017_JECup );
	MT2Analysis<MT2EstimateTree>* sigTemp_2017_JECdn = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HH.root", ggDir_2017_JECdn.c_str() ), signalList_tmp[i]  );
	signals_2017_JECdn.push_back( sigTemp_2017_JECdn );
      }
      else if (doHZ){
	MT2Analysis<MT2EstimateTree>* sigTemp_2017 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_2017.c_str() ), signalList_tmp[i]  );
	signals_2017.push_back( sigTemp_2017 );
	MT2Analysis<MT2EstimateTree>* sigTemp_2017_JECup = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_2017_JECup.c_str() ), signalList_tmp[i]  );
	signals_2017_JECup.push_back( sigTemp_2017_JECup );
	MT2Analysis<MT2EstimateTree>* sigTemp_2017_JECdn = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_2017_JECdn.c_str() ), signalList_tmp[i]  );
	signals_2017_JECdn.push_back( sigTemp_2017_JECdn );
      }else{ 
	MT2Analysis<MT2EstimateTree>* sigTemp_2017 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_WH.root", ggDir_2017.c_str() ), signalList_tmp[i]  );
	signals_2017.push_back( sigTemp_2017 );
	MT2Analysis<MT2EstimateTree>* sigTemp_2017_JECup = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_WH.root", ggDir_2017_JECup.c_str() ), signalList_tmp[i]  );
	signals_2017_JECup.push_back( sigTemp_2017_JECup );
	MT2Analysis<MT2EstimateTree>* sigTemp_2017_JECdn = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_WH.root", ggDir_2017_JECdn.c_str() ), signalList_tmp[i]  );
	signals_2017_JECdn.push_back( sigTemp_2017_JECdn );
      } 
    }

  vector< MT2Analysis<MT2EstimateTree>* > signals_HH_HH0p25_2017;
  vector< MT2Analysis<MT2EstimateTree>* > signals_HH_HH0p25_2017_JECup;
  vector< MT2Analysis<MT2EstimateTree>* > signals_HH_HH0p25_2017_JECdn;
  if( doComb ) 
    if( doHZ )
      for( unsigned int i=0; i<signalList_HH_HH0p25.size(); i++){

	if( (i>= sigPerPart * partNr) || (i < ((partNr-1)*sigPerPart ))   )
	  continue;

	MT2Analysis<MT2EstimateTree>* sigTemp_HH_HH0p25_2017 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_2017.c_str() ), signalList_HH_HH0p25[i]  );
	signals_HH_HH0p25_2017.push_back( sigTemp_HH_HH0p25_2017 );

	MT2Analysis<MT2EstimateTree>* sigTemp_HH_HH0p25_2017_JECup = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_2017_JECup.c_str() ), signalList_HH_HH0p25[i]  );
	signals_HH_HH0p25_2017_JECup.push_back( sigTemp_HH_HH0p25_2017_JECup );
	MT2Analysis<MT2EstimateTree>* sigTemp_HH_HH0p25_2017_JECdn = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir_2017_JECdn.c_str() ), signalList_HH_HH0p25[i]  );
	signals_HH_HH0p25_2017_JECdn.push_back( sigTemp_HH_HH0p25_2017_JECdn );
      }



  std::string outFileName( Form("%s/ws_%s", outputdir.c_str(), sigName.c_str() ));
  if( doComb )       outFileName += "_comb";
  if( mt2_bins.size() > 1 ) outFileName += "_mt2_30";
   
  TFile* dataFile = new TFile(  Form("%s_part%d.root", outFileName.c_str(), partNr) , "RECREATE" );
  RooWorkspace ws_sig ("ws_sig");
  //  RooWorkspace ws_bg  ("ws_bg");
  RooWorkspace ws_data("ws_data");
  
  std::string dataAmountName( Form("%s/dataAmount",  outputdir.c_str() ));
  if( doComb )       dataAmountName += "_comb";
  if( mt2_bins.size() > 1 ) dataAmountName += "_mt2_30";
  std::ofstream dataAmount( Form("%s_%s.txt" , dataAmountName.c_str() , sigName.c_str() ) );
    
  // std::string bgSumAmountName( Form("%s/bgSumAmount",  outputdir.c_str() ));
  // if( doComb )             bgSumAmountName += "_comb";
  // if( mt2_bins.size() > 1 ) bgSumAmountName += "_mt2_30";
  // std::ofstream bgSumAmount( Form("%s_%s.txt" , bgSumAmountName.c_str() , sigName.c_str() ) );

  std::string isrName( Form("%s/isr",  outputdir.c_str() ));
  if( doComb )    isrName += "_comb";
  if( mt2_bins.size() > 1 ) isrName += "_mt2_30";
  std::ofstream isr( Form("%s_%s.txt" , isrName.c_str(), sigName.c_str() ) );

  std::string lepsf_higgsName( Form("%s/lepsf_higgs",  outputdir.c_str() ));
  if( doComb )    lepsf_higgsName += "_comb";
  if( mt2_bins.size() > 1 ) lepsf_higgsName += "_mt2_30";
  std::ofstream lepsf_higgs( Form("%s_%s.txt" , lepsf_higgsName.c_str(), sigName.c_str() ) );

  std::string gammasf_higgsName( Form("%s/gammasf_higgs",  outputdir.c_str() ));
  if( doComb )    gammasf_higgsName += "_comb";
  if( mt2_bins.size() > 1 ) gammasf_higgsName += "_mt2_30";
  std::ofstream gammasf_higgs( Form("%s_%s.txt" , gammasf_higgsName.c_str(), sigName.c_str() ) );

  std::string btagsf_heavy_higgsName( Form("%s/btagsf_heavy_higgs",  outputdir.c_str() ));
  if( doComb )    btagsf_heavy_higgsName += "_comb";
  if( mt2_bins.size() > 1 ) btagsf_heavy_higgsName += "_mt2_30";
  std::ofstream btagsf_heavy_higgs( Form("%s_%s.txt" , btagsf_heavy_higgsName.c_str(), sigName.c_str() ) );

  std::string btagsf_light_higgsName( Form("%s/btagsf_light_higgs",  outputdir.c_str() ));
  if( doComb )    btagsf_light_higgsName += "_comb";
  if( mt2_bins.size() > 1 ) btagsf_light_higgsName += "_mt2_30";
  std::ofstream btagsf_light_higgs( Form("%s_%s.txt" , btagsf_light_higgsName.c_str(), sigName.c_str() ) );

  std::string lepsf_higgs_2017Name( Form("%s/lepsf_higgs_2017",  outputdir.c_str() ));
  if( doComb )    lepsf_higgs_2017Name += "_comb";
  if( mt2_bins.size() > 1 ) lepsf_higgs_2017Name += "_mt2_30";
  std::ofstream lepsf_higgs_2017( Form("%s_%s.txt" , lepsf_higgs_2017Name.c_str(), sigName.c_str() ) );

  std::string gammasf_higgs_2017Name( Form("%s/gammasf_higgs_2017",  outputdir.c_str() ));
  if( doComb )    gammasf_higgs_2017Name += "_comb";
  if( mt2_bins.size() > 1 ) gammasf_higgs_2017Name += "_mt2_30";
  std::ofstream gammasf_higgs_2017( Form("%s_%s.txt" , gammasf_higgs_2017Name.c_str(), sigName.c_str() ) );

  std::string btagsf_heavy_higgs_2017Name( Form("%s/btagsf_heavy_higgs_2017",  outputdir.c_str() ));
  if( doComb )    btagsf_heavy_higgs_2017Name += "_comb";
  if( mt2_bins.size() > 1 ) btagsf_heavy_higgs_2017Name += "_mt2_30";
  std::ofstream btagsf_heavy_higgs_2017( Form("%s_%s.txt" , btagsf_heavy_higgs_2017Name.c_str(), sigName.c_str() ) );

  std::string btagsf_light_higgs_2017Name( Form("%s/btagsf_light_higgs_2017",  outputdir.c_str() ));
  if( doComb )    btagsf_light_higgs_2017Name += "_comb";
  if( mt2_bins.size() > 1 ) btagsf_light_higgs_2017Name += "_mt2_30";
  std::ofstream btagsf_light_higgs_2017( Form("%s_%s.txt" , btagsf_light_higgs_2017Name.c_str(), sigName.c_str() ) );


  std::string scaleVar_higgsName( Form("%s/scaleVar_higgs",  outputdir.c_str() ));
  if( doComb )    scaleVar_higgsName += "_comb";
  if( mt2_bins.size() > 1 ) scaleVar_higgsName += "_mt2_30";
  std::ofstream scaleVar_higgs( Form("%s_%s.txt" , scaleVar_higgsName.c_str(), sigName.c_str() ) );

  std::string scaleVar_higgs_2017Name( Form("%s/scaleVar_higgs_2017",  outputdir.c_str() ));
  if( doComb )    scaleVar_higgs_2017Name += "_comb";
  if( mt2_bins.size() > 1 ) scaleVar_higgs_2017Name += "_mt2_30";
  std::ofstream scaleVar_higgs_2017( Form("%s_%s.txt" , scaleVar_higgs_2017Name.c_str(), sigName.c_str() ) );

  std::string pdfVar_higgsName( Form("%s/pdfVar_higgs",  outputdir.c_str() ));
  if( doComb )    pdfVar_higgsName += "_comb";
  if( mt2_bins.size() > 1 ) pdfVar_higgsName += "_mt2_30";
  std::ofstream pdfVar_higgs( Form("%s_%s.txt" , pdfVar_higgsName.c_str(), sigName.c_str() ) );

  std::string pdfVar_higgs_2017Name( Form("%s/pdfVar_higgs_2017",  outputdir.c_str() ));
  if( doComb )    pdfVar_higgs_2017Name += "_comb";
  if( mt2_bins.size() > 1 ) pdfVar_higgs_2017Name += "_mt2_30";
  std::ofstream pdfVar_higgs_2017( Form("%s_%s.txt" , pdfVar_higgs_2017Name.c_str(), sigName.c_str() ) );

  std::string sigAmountName( Form("%s/sigAmount_%s",  outputdir.c_str(), sigName.c_str() ));
  if( doComb )    sigAmountName += "_comb";
  if( mt2_bins.size() > 1 ) sigAmountName += "_mt2_30";
  std::ofstream sigAmount( Form("%s_%s.txt", sigAmountName.c_str() , sigName.c_str()) );

  std::string higgsGluGluName( Form("%s/higgsGluGlu",  outputdir.c_str() ));
  if( doComb )    higgsGluGluName += "_comb";
  if( mt2_bins.size() > 1 ) higgsGluGluName += "_mt2_30";
  std::ofstream higgsGluGlu( Form("%s_%s.txt" , higgsGluGluName.c_str(), sigName.c_str() ) );

  std::string higgsVBFName( Form("%s/higgsVBF",  outputdir.c_str() ));
  if( doComb )    higgsVBFName += "_comb";
  if( mt2_bins.size() > 1 ) higgsVBFName += "_mt2_30";
  std::ofstream higgsVBF( Form("%s_%s.txt" , higgsVBFName.c_str(), sigName.c_str() ) );

  std::string higgsVHName( Form("%s/higgsVH",  outputdir.c_str() ));
  if( doComb )    higgsVHName += "_comb";
  if( mt2_bins.size() > 1 ) higgsVHName += "_mt2_30";
  std::ofstream higgsVH( Form("%s_%s.txt" , higgsVHName.c_str(), sigName.c_str() ) );

  std::string higgsttHName( Form("%s/higgsttH",  outputdir.c_str() ));
  if( doComb )    higgsttHName += "_comb";
  if( mt2_bins.size() > 1 ) higgsttHName += "_mt2_30";
  std::ofstream higgsttH( Form("%s_%s.txt" , higgsttHName.c_str(), sigName.c_str() ) );

  std::string higgsbbHName( Form("%s/higgsbbH",  outputdir.c_str() ));
  if( doComb )    higgsbbHName += "_comb";
  if( mt2_bins.size() > 1 ) higgsbbHName += "_mt2_30";
  std::ofstream higgsbbH( Form("%s_%s.txt" , higgsbbHName.c_str(), sigName.c_str() ) );

  // std::string higgsTHXName( Form("%s/higgsTHX",  outputdir.c_str() ));
  // if( doComb )    higgsTHXName += "_comb";
  // if( mt2_bins.size() > 1 ) higgsTHXName += "_mt2_30";
  // std::ofstream higgsTHX( Form("%s_%s.txt" , higgsTHXName.c_str(), sigName.c_str() ) );
  // std::string higgsdiHName( Form("%s/higgsdiH",  outputdir.c_str() ));
  // if( doComb )    higgsdiHName += "_comb";
  // if( mt2_bins.size() > 1 ) higgsdiHName += "_mt2_30";
  // std::ofstream higgsdiH( Form("%s_%s.txt" , higgsdiHName.c_str(), sigName.c_str() ) );

  std::string higgsAmountName( Form("%s/higgsAmount",  outputdir.c_str() ));
  if( doComb )    higgsAmountName += "_comb";
  if( mt2_bins.size() > 1 ) higgsAmountName += "_mt2_30";
  std::ofstream higgsAmount( Form("%s_%s.txt" , higgsAmountName.c_str(), sigName.c_str() ) );

  std::string higgsSyst_2016_JECupName( Form("%s/higgsSyst_2016_JECup",  outputdir.c_str() ));
  if( doComb )    higgsSyst_2016_JECupName += "_comb";
  if( mt2_bins.size() > 1 ) higgsSyst_2016_JECupName += "_mt2_30";
  std::ofstream higgsSyst_2016_JECup( Form("%s_%s.txt" , higgsSyst_2016_JECupName.c_str(), sigName.c_str() ) );

  std::string higgsSyst_2016_JECdnName( Form("%s/higgsSyst_2016_JECdn",  outputdir.c_str() ));
  if( doComb )    higgsSyst_2016_JECdnName += "_comb";
  if( mt2_bins.size() > 1 ) higgsSyst_2016_JECdnName += "_mt2_30";
  std::ofstream higgsSyst_2016_JECdn( Form("%s_%s.txt" , higgsSyst_2016_JECdnName.c_str(), sigName.c_str() ) );
 
  std::string higgsSyst_2017_JECupName( Form("%s/higgsSyst_2017_JECup",  outputdir.c_str() ));
  if( doComb )    higgsSyst_2017_JECupName += "_comb";
  if( mt2_bins.size() > 1 ) higgsSyst_2017_JECupName += "_mt2_30";
  std::ofstream higgsSyst_2017_JECup( Form("%s_%s.txt" , higgsSyst_2017_JECupName.c_str(), sigName.c_str() ) );

  std::string higgsSyst_2017_JECdnName( Form("%s/higgsSyst_2017_JECdn",  outputdir.c_str() ));
  if( doComb )    higgsSyst_2017_JECdnName += "_comb";
  if( mt2_bins.size() > 1 ) higgsSyst_2017_JECdnName += "_mt2_30";
  std::ofstream higgsSyst_2017_JECdn( Form("%s_%s.txt" , higgsSyst_2017_JECdnName.c_str(), sigName.c_str() ) );



        
  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

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

    string regionSaveName =  regionNames[i];

    //INITIAL SELECTION 
    string treeSel = "((( ptGamma0/h_mass) > 1./3) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 ) && passPrefire" ;

    string treeSel_genMET = treeSel;

    treeSel += " && " + regionSelection[i];

    treeSel_genMET += " && " + regionSelection_genMET[i];
    //string treeSelWideEta = "((( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=2.4  && fabs(etaGamma1)<=2.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )" ;
      

    TH2D* h2D_genMET         = new TH2D( Form("h2D_genMET_%s",      regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_genMET_2017    = new TH2D( Form("h2D_genMET_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );

    TH2D* h2D_isr    = new TH2D( Form("h2D_isr_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_isr_err = new TH2D(Form("h2D_isr_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_isr_2017    = new TH2D( Form("h2D_isr_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_isr_2017_err = new TH2D(Form("h2D_isr_2017_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );

    TH2D* h2D_lepsf    = new TH2D( Form("h2D_lepsf_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_lepsf_err = new TH2D(Form("h2D_lepsf_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_lepsf_2017    = new TH2D( Form("h2D_lepsf_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_lepsf_2017_err = new TH2D(Form("h2D_lepsf_2017_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );

    TH2D* h2D_gammasf    = new TH2D( Form("h2D_gammasf_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_gammasf_err = new TH2D(Form("h2D_gammasf_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_gammasf_2017    = new TH2D( Form("h2D_gammasf_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_gammasf_2017_err = new TH2D(Form("h2D_gammasf_2017_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );

    TH2D* h2D_btagsf_heavy    = new TH2D( Form("h2D_btagsf_heavy_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_btagsf_heavy_err = new TH2D(Form("h2D_btagsf_heavy_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_btagsf_heavy_2017    = new TH2D( Form("h2D_btagsf_heavy_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_btagsf_heavy_2017_err = new TH2D(Form("h2D_btagsf_heavy_2017_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );

    TH2D* h2D_btagsf_light    = new TH2D( Form("h2D_btagsf_light_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_btagsf_light_err = new TH2D(Form("h2D_btagsf_light_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_btagsf_light_2017    = new TH2D( Form("h2D_btagsf_light_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_btagsf_light_2017_err = new TH2D(Form("h2D_btagsf_light_2017_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );

    TH2D* h2D_JECup    = new TH2D( Form("h2D_JECup_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_JECup_2017    = new TH2D( Form("h2D_JECup_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_JECdn    = new TH2D( Form("h2D_JECdn_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_JECdn_2017    = new TH2D( Form("h2D_JECdn_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );

    TH2D* h2D_scaleVar    = new TH2D( Form("h2D_scaleVar_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_scaleVar_2017    = new TH2D( Form("h2D_scaleVar_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_pdfVar    = new TH2D( Form("h2D_pdfVar_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );
    TH2D* h2D_pdfVar_2017    = new TH2D( Form("h2D_pdfVar_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 33, -12.5, 812.5 );


    if( doSMHDonly ==1 ){

      TTree* tree_data_ini = data->get(thisRegion)->tree;
      TTree* tree_higgs_ini = higgs->get(thisRegion)->tree;
      // TTree* tree_qcd_ini = qcd->get(thisRegion)->tree;
      // //    TTree* tree_diH_ini = diH->get(thisRegion)->tree;
      // TTree* tree_diPhoton_ini = diPhoton->get(thisRegion)->tree;
      // TTree* tree_gjets_ini = gjets->get(thisRegion)->tree;

      TTree* tree_data = tree_data_ini->CopyTree(treeSel.c_str() );
      TTree* tree_higgs = tree_higgs_ini->CopyTree(treeSel.c_str() );
      // TTree* tree_qcd = tree_qcd_ini->CopyTree(treeSel.c_str() );
      // //    TTree* tree_diH = tree_diH_ini->CopyTree(treeSel.c_str() );
      // TTree* tree_diPhoton = tree_diPhoton_ini->CopyTree(treeSel.c_str() );
      // TTree* tree_gjets = tree_gjets_ini->CopyTree(treeSel.c_str() );

      TTree* tree_higgs_ini_JECup = higgs_JECup->get(thisRegion)->tree;
      TTree* tree_higgs_JECup = tree_higgs_ini_JECup->CopyTree(treeSel.c_str() );
      TTree* tree_higgs_ini_2017_JECup;
      TTree* tree_higgs_2017_JECup;
      TTree* tree_higgs_ini_JECdn = higgs_JECdn->get(thisRegion)->tree;
      TTree* tree_higgs_JECdn = tree_higgs_ini_JECdn->CopyTree(treeSel.c_str() );
      TTree* tree_higgs_ini_2017_JECdn;
      TTree* tree_higgs_2017_JECdn;

      //comb with 2017
      TTree* tree_data_ini_2017;
      TTree* tree_higgs_ini_2017;
      // TTree* tree_qcd_ini_2017;
      // //    TTree* tree_diH_ini_2017;
      // TTree* tree_diPhoton_ini_2017;
      // TTree* tree_gjets_ini_2017;

      TTree* tree_data_2017;
      TTree* tree_higgs_2017;
      // TTree* tree_qcd_2017;
      // //    TTree* tree_diH_2017;
      // TTree* tree_diPhoton_2017;
      // TTree* tree_gjets_2017;

      if( doComb ){
	tree_data_ini_2017 = data_2017->get(thisRegion)->tree;
	tree_higgs_ini_2017 = higgs_2017->get(thisRegion)->tree;
	// tree_qcd_ini_2017 = qcd_2017->get(thisRegion)->tree;
	// //      tree_diH_ini_2017 = diH_2017->get(thisRegion)->tree;
	// tree_diPhoton_ini_2017 = diPhoton_2017->get(thisRegion)->tree;
	// tree_gjets_ini_2017 = gjets_2017->get(thisRegion)->tree;

	tree_data_2017 = tree_data_ini_2017->CopyTree(treeSel.c_str() );
	tree_higgs_2017 = tree_higgs_ini_2017->CopyTree(treeSel.c_str() );
	// tree_qcd_2017 = tree_qcd_ini_2017->CopyTree(treeSel.c_str() );
	// //      tree_diH_2017 = tree_diH_ini_2017->CopyTree(treeSel.c_str() );
	// tree_diPhoton_2017 = tree_diPhoton_ini_2017->CopyTree(treeSel.c_str() );
	// tree_gjets_2017 = tree_gjets_ini_2017->CopyTree(treeSel.c_str() );

	tree_higgs_ini_2017_JECup = higgs_2017_JECup->get(thisRegion)->tree;
	tree_higgs_2017_JECup = tree_higgs_ini_2017_JECup->CopyTree(treeSel.c_str() );
	tree_higgs_ini_2017_JECdn = higgs_2017_JECdn->get(thisRegion)->tree;
	tree_higgs_2017_JECdn = tree_higgs_ini_2017_JECdn->CopyTree(treeSel.c_str() );
      }


      string higgsGluGlu_treeSel = treeSel + " && (id==951) && (h_mass>=122 && h_mass<=129) ";
      TTree* tree_higgsGluGlu_ini = higgs->get(thisRegion)->tree;
      TTree* tree_higgsGluGlu = tree_higgsGluGlu_ini->CopyTree(higgsGluGlu_treeSel.c_str() );
      TTree* tree_higgsGluGlu_ini_2017 = higgs_2017->get(thisRegion)->tree;
      TTree* tree_higgsGluGlu_2017 = tree_higgsGluGlu_ini_2017->CopyTree((higgsGluGlu_treeSel.c_str()) );
      RooDataSet higgsGluGlu_ds(Form("higgsGluGlu_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsGluGlu) );
      RooDataSet higgsGluGlu_ds_2017(Form("higgsGluGlu_2017_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsGluGlu_2017) );

      string higgsVBF_treeSel = treeSel + " && (id==952) && (h_mass>=122 && h_mass<=129) ";
      TTree* tree_higgsVBF_ini = higgs->get(thisRegion)->tree;
      TTree* tree_higgsVBF = tree_higgsVBF_ini->CopyTree(higgsVBF_treeSel.c_str() );
      TTree* tree_higgsVBF_ini_2017 = higgs_2017->get(thisRegion)->tree;
      TTree* tree_higgsVBF_2017 = tree_higgsVBF_ini_2017->CopyTree((higgsVBF_treeSel.c_str()) );
      RooDataSet higgsVBF_ds(Form("higgsVBF_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsVBF) );
      RooDataSet higgsVBF_ds_2017(Form("higgsVBF_2017_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsVBF_2017) );

      string higgsVH_treeSel = treeSel + " && (id==953) && (h_mass>=122 && h_mass<=129) ";
      TTree* tree_higgsVH_ini = higgs->get(thisRegion)->tree;
      TTree* tree_higgsVH = tree_higgsVH_ini->CopyTree(higgsVH_treeSel.c_str() );
      TTree* tree_higgsVH_ini_2017 = higgs_2017->get(thisRegion)->tree;
      TTree* tree_higgsVH_2017 = tree_higgsVH_ini_2017->CopyTree((higgsVH_treeSel.c_str()) );
      RooDataSet higgsVH_ds(Form("higgsVH_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsVH) );
      RooDataSet higgsVH_ds_2017(Form("higgsVH_2017_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsVH_2017) );

      string higgsttH_treeSel = treeSel + " && ((id==954) || (id==955)) && (h_mass>=122 && h_mass<=129) ";
      TTree* tree_higgsttH_ini = higgs->get(thisRegion)->tree;
      TTree* tree_higgsttH = tree_higgsttH_ini->CopyTree(higgsttH_treeSel.c_str() );
      TTree* tree_higgsttH_ini_2017 = higgs_2017->get(thisRegion)->tree;
      TTree* tree_higgsttH_2017 = tree_higgsttH_ini_2017->CopyTree((higgsttH_treeSel.c_str()) );
      RooDataSet higgsttH_ds(Form("higgsttH_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsttH) );
      RooDataSet higgsttH_ds_2017(Form("higgsttH_2017_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsttH_2017) );

      string higgsbbH_treeSel = treeSel + " && ( (id==956) || (id==957)) && (h_mass>=122 && h_mass<=129) ";
      TTree* tree_higgsbbH_ini = higgs->get(thisRegion)->tree;
      TTree* tree_higgsbbH = tree_higgsbbH_ini->CopyTree(higgsbbH_treeSel.c_str() );
      TTree* tree_higgsbbH_ini_2017 = higgs_2017->get(thisRegion)->tree;
      TTree* tree_higgsbbH_2017 = tree_higgsbbH_ini_2017->CopyTree((higgsbbH_treeSel.c_str()) );
      RooDataSet higgsbbH_ds(Form("higgsbbH_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsbbH) );
      RooDataSet higgsbbH_ds_2017(Form("higgsbbH_2017_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsbbH_2017) );

      // string higgsTHX_treeSel = treeSel + " && ( (id==958) || (id==959)) && (h_mass>=122 && h_mass<=129) ";
      // TTree* tree_higgsTHX_ini = higgs->get(thisRegion)->tree;
      // TTree* tree_higgsTHX = tree_higgsTHX_ini->CopyTree(higgsTHX_treeSel.c_str() );
      // TTree* tree_higgsTHX_ini_2017 = higgs_2017->get(thisRegion)->tree;
      // TTree* tree_higgsTHX_2017 = tree_higgsTHX_ini_2017->CopyTree((higgsTHX_treeSel.c_str()) );
      // RooDataSet higgsTHX_ds(Form("higgsTHX_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsTHX) );
      // RooDataSet higgsTHX_ds_2017(Form("higgsTHX_2017_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgsTHX_2017) );


      //Import 2017 
      RooDataSet higgs_ds(Form("higgs_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs) );
      //    RooDataSet diH_ds(Form("diH_2016_125_13TeV_%s", regionSaveName.c_str() ), "DiH", vars, WeightVar(weight), Import(*tree_diH) );
      //bg    RooDataSet bg_ds(Form("bg_2016_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton) ) ;
      // RooDataSet gjets_ds("gjets_2016_ds", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets) ) ;
      // RooDataSet qcd_ds("qcd_2016_ds", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd) ) ;

      RooDataSet data_ds( Form("data_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data) ) ;
   
      RooDataSet higgs_ds_JECup(Form("higgs_JECup_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_JECup) );
      RooDataSet higgs_ds_JECdn(Form("higgs_JECdn_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_JECdn) );

      if(doComb){
	//Import 2017
	RooDataSet higgs_ds_2017(Form("higgs_2017_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_2017) ); 

	//bg      RooDataSet bg_ds_2017(Form("bg_2017_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton_2017) ) ;
	// RooDataSet gjets_ds_2017("gjets_ds_2017", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets_2017) ) ;
	// RooDataSet qcd_ds_2017("qcd_ds_2017", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd_2017) ) ;
   
	RooDataSet data_ds_2017( Form("data_2017_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data_2017) ) ;
   
	//bg      bg_ds_2017.append( qcd_ds_2017 );   //append bg data into one
	//bg      bg_ds_2017.append( gjets_ds_2017 ); //append bg data into one

	RooDataSet higgs_ds_130_2017(  higgs_ds_2017, Form("higgs_2017_130_13TeV_%s", regionSaveName.c_str() ) );
	//bg      RooDataSet bg_ds_130_2017   (  bg_ds_2017,    Form("bg_2017_130_13TeV_%s",    regionSaveName.c_str() ) );
	//    RooDataSet data_ds_130_2017 (  data_ds_2017,  Form("data_2017_130_13TeV_%s",  regionSaveName.c_str() ) );
	RooDataSet higgs_ds_125_2017(  higgs_ds_2017, Form("higgs_2017_125_13TeV_%s", regionSaveName.c_str() ) );
	//bg      RooDataSet bg_ds_125_2017   (  bg_ds_2017,    Form("bg_2017_125_13TeV_%s",    regionSaveName.c_str() ) );
	//    RooDataSet data_ds_125_2017 (  data_ds_2017,  Form("data_2017_125_13TeV_%s",  regionSaveName.c_str() ) );
	RooDataSet higgs_ds_120_2017(  higgs_ds_2017, Form("higgs_2017_120_13TeV_%s", regionSaveName.c_str() ) );
	//bg      RooDataSet bg_ds_120_2017   (  bg_ds_2017,    Form("bg_2017_120_13TeV_%s",    regionSaveName.c_str() ) );
	//    RooDataSet data_ds_120_2017 (  data_ds_2017,  Form("data_2017_120_13TeV_%s",  regionSaveName.c_str() ) );
	((RooRealVar*)higgs_ds_130_2017.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
	//    ((RooRealVar*)data_ds_130_2017.addColumn( hgg_mass_130_2017 ))->setRange(100.,180.); 
	((RooRealVar*)higgs_ds_125_2017.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
	//    ((RooRealVar*)data_ds_125_2017.addColumn( hgg_mass_125_2017 ))->setRange(100.,180.); 
	((RooRealVar*)higgs_ds_120_2017.addColumn( hgg_mass_120 ))->setRange(100.,180.);
	//    ((RooRealVar*)data_ds_120_2017.addColumn( hgg_mass_120_2017 ))->setRange(100.,180.);

	//  ws_data.import( data_ds_125 );
	//bg      ws_bg.import( bg_ds_125_2017   );
	ws_sig.import( higgs_ds_125_2017  );
	//  ws_data.import( data_ds_130 );
	//bg      ws_bg.import( bg_ds_130_2017   );
	ws_sig.import( higgs_ds_130_2017  );
	//  ws_data.import( data_ds_120 );
	//bg      ws_bg.import( bg_ds_120_2017   );
	ws_sig.import( higgs_ds_120_2017  );

	data_ds.append( data_ds_2017 );

	RooDataSet higgs_ds_2017_JECup(Form("higgs_2017_JECup_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_2017_JECup) );
	RooDataSet higgs_ds_2017_JECdn(Form("higgs_2017_JECdn_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_2017_JECdn) );
	RooDataSet higgs_ds_125_2017_JECup(  higgs_ds_2017_JECup, Form("higgs_2017_JECup_125_13TeV_%s", regionSaveName.c_str() ) );
	RooDataSet higgs_ds_125_2017_JECdn(  higgs_ds_2017_JECdn, Form("higgs_2017_JECdn_125_13TeV_%s", regionSaveName.c_str() ) );
	((RooRealVar*)higgs_ds_125_2017_JECup.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
	((RooRealVar*)higgs_ds_125_2017_JECdn.addColumn( hgg_mass_125 ))->setRange(100.,180.); 

	// std::cout  << regionSaveName << "    " <<  higgs_ds_125_2017.sumEntries() << std::endl;
	// std::cout  << regionSaveName << "    " <<  higgs_ds_125_2017_JECup.sumEntries() << std::endl;
	// std::cout  << regionSaveName << "    " <<  ((higgs_ds_125_2017.sumEntries() - higgs_ds_125_2017_JECup.sumEntries()))  << std::endl;

	higgsSyst_2017_JECup << regionSaveName << "    " <<  (( higgs_ds_125_2017_JECup.sumEntries() - higgs_ds_125_2017.sumEntries() )/float(higgs_ds_125_2017.sumEntries()  )) << std::endl;
	higgsSyst_2017_JECdn << regionSaveName << "    " <<  ((higgs_ds_125_2017.sumEntries() - higgs_ds_125_2017_JECdn.sumEntries() )/float(higgs_ds_125_2017.sumEntries() ) ) << std::endl;


	//////////////////////////////////////////////////////////////////////////////
	//////////////////  UNCERTS FOR HIGGS_2017 ////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	/////////////////photons
	TH1D* h_gammasf_higgs_2017    = new TH1D("h_gammasf_higgs_2017",    "", 100, 0, 2);
	TH1D* h_gammasf_higgs_2017_UP = new TH1D("h_gammasf_higgs_2017_UP", "", 100, 0, 2);
	TH1D* h_gammasf_higgs_2017_DN = new TH1D("h_gammasf_higgs_2017_DN", "", 100, 0, 2);
 
	tree_higgs_2017->Draw("weight_gammasf    >> h_gammasf_higgs_2017");
	tree_higgs_2017->Draw("weight_gammasf_UP >> h_gammasf_higgs_2017_UP");
	tree_higgs_2017->Draw("weight_gammasf_DN >> h_gammasf_higgs_2017_DN");

	float sign_gammasf_higgs_2017 = 1;
	float gammasf_higgs_2017_uncert = max (fabs(h_gammasf_higgs_2017_UP->GetMean() - h_gammasf_higgs_2017->GetMean()), fabs( h_gammasf_higgs_2017->GetMean() - h_gammasf_higgs_2017_DN->GetMean() ) );

	if( fabs(h_gammasf_higgs_2017_UP->GetMean() - h_gammasf_higgs_2017->GetMean()) >  fabs(  h_gammasf_higgs_2017->GetMean() - h_gammasf_higgs_2017_DN->GetMean() ))
	  sign_gammasf_higgs_2017 = ((h_gammasf_higgs_2017_UP->GetMean() - h_gammasf_higgs_2017->GetMean())  > 0 ) - ((h_gammasf_higgs_2017_UP->GetMean() - h_gammasf_higgs_2017->GetMean())  < 0) ;
	else
	  sign_gammasf_higgs_2017 = ((h_gammasf_higgs_2017->GetMean() - h_gammasf_higgs_2017_DN->GetMean())  > 0 ) - ((h_gammasf_higgs_2017->GetMean() - h_gammasf_higgs_2017_DN->GetMean())  < 0) ;
	  
	gammasf_higgs_2017 << regionSaveName << " \t " <<    sign_gammasf_higgs_2017* gammasf_higgs_2017_uncert  <<std::endl;

	delete h_gammasf_higgs_2017; delete h_gammasf_higgs_2017_UP; delete h_gammasf_higgs_2017_DN;

	/////////////////leptons
	TH1D* h_lepsf_higgs_2017    = new TH1D("h_lepsf_higgs_2017",    "", 100, 0, 2);
	TH1D* h_lepsf_higgs_2017_UP = new TH1D("h_lepsf_higgs_2017_UP", "", 100, 0, 2);
	TH1D* h_lepsf_higgs_2017_DN = new TH1D("h_lepsf_higgs_2017_DN", "", 100, 0, 2);
 
	tree_higgs_2017->Draw("weight_lepsf    >> h_lepsf_higgs_2017");
	tree_higgs_2017->Draw("weight_lepsf_UP >> h_lepsf_higgs_2017_UP");
	tree_higgs_2017->Draw("weight_lepsf_DN >> h_lepsf_higgs_2017_DN");

	float sign_lepsf_higgs_2017 = 1;
	float lepsf_higgs_2017_uncert = max (fabs(h_lepsf_higgs_2017_UP->GetMean() - h_lepsf_higgs_2017->GetMean()), fabs( h_lepsf_higgs_2017->GetMean() - h_lepsf_higgs_2017_DN->GetMean() ) );

	if( fabs(h_lepsf_higgs_2017_UP->GetMean() - h_lepsf_higgs_2017->GetMean()) >  fabs(  h_lepsf_higgs_2017->GetMean() - h_lepsf_higgs_2017_DN->GetMean() ))
	  sign_lepsf_higgs_2017 = ((h_lepsf_higgs_2017_UP->GetMean() - h_lepsf_higgs_2017->GetMean())  > 0 ) - ((h_lepsf_higgs_2017_UP->GetMean() - h_lepsf_higgs_2017->GetMean())  < 0) ;
	else
	  sign_lepsf_higgs_2017 = ((h_lepsf_higgs_2017->GetMean() - h_lepsf_higgs_2017_DN->GetMean())  > 0 ) - ((h_lepsf_higgs_2017->GetMean() - h_lepsf_higgs_2017_DN->GetMean())  < 0) ;
	  
	lepsf_higgs_2017 << regionSaveName << " \t " <<    sign_lepsf_higgs_2017* lepsf_higgs_2017_uncert  <<std::endl;

	delete h_lepsf_higgs_2017; delete h_lepsf_higgs_2017_UP; delete h_lepsf_higgs_2017_DN;

	//btag
	TH1D* h_btagsf_heavy_higgs_2017    = new TH1D("h_btagsf_heavy_higgs_2017",    "", 100, 0, 2);
	TH1D* h_btagsf_heavy_higgs_2017_UP = new TH1D("h_btagsf_heavy_higgs_2017_UP", "", 100, 0, 2);
	TH1D* h_btagsf_heavy_higgs_2017_DN = new TH1D("h_btagsf_heavy_higgs_2017_DN", "", 100, 0, 2);
 
	tree_higgs_2017->Draw("weight_btagsf    >> h_btagsf_heavy_higgs_2017");
	tree_higgs_2017->Draw("weight_btagsf_heavy_UP >> h_btagsf_heavy_higgs_2017_UP");
	tree_higgs_2017->Draw("weight_btagsf_heavy_DN >> h_btagsf_heavy_higgs_2017_DN");

	float sign_btagsf_heavy_higgs_2017 = 1;
	float btagsf_heavy_higgs_2017_uncert = max (fabs(h_btagsf_heavy_higgs_2017_UP->GetMean() - h_btagsf_heavy_higgs_2017->GetMean()), fabs( h_btagsf_heavy_higgs_2017->GetMean() - h_btagsf_heavy_higgs_2017_DN->GetMean() ) );

	if( fabs(h_btagsf_heavy_higgs_2017_UP->GetMean() - h_btagsf_heavy_higgs_2017->GetMean()) >  fabs(  h_btagsf_heavy_higgs_2017->GetMean() - h_btagsf_heavy_higgs_2017_DN->GetMean() ))
	  sign_btagsf_heavy_higgs_2017 = ((h_btagsf_heavy_higgs_2017_UP->GetMean() - h_btagsf_heavy_higgs_2017->GetMean())  > 0 ) - ((h_btagsf_heavy_higgs_2017_UP->GetMean() - h_btagsf_heavy_higgs_2017->GetMean())  < 0) ;
	else
	  sign_btagsf_heavy_higgs_2017 = ((h_btagsf_heavy_higgs_2017->GetMean() - h_btagsf_heavy_higgs_2017_DN->GetMean())  > 0 ) - ((h_btagsf_heavy_higgs_2017->GetMean() - h_btagsf_heavy_higgs_2017_DN->GetMean())  < 0) ;
	  
	btagsf_heavy_higgs_2017 << regionSaveName << " \t " <<    sign_btagsf_heavy_higgs_2017* btagsf_heavy_higgs_2017_uncert  <<std::endl;

	delete h_btagsf_heavy_higgs_2017; delete h_btagsf_heavy_higgs_2017_UP; delete h_btagsf_heavy_higgs_2017_DN;


	TH1D* h_btagsf_light_higgs_2017    = new TH1D("h_btagsf_light_higgs_2017",    "", 100, 0, 2);
	TH1D* h_btagsf_light_higgs_2017_UP = new TH1D("h_btagsf_light_higgs_2017_UP", "", 100, 0, 2);
	TH1D* h_btagsf_light_higgs_2017_DN = new TH1D("h_btagsf_light_higgs_2017_DN", "", 100, 0, 2);
 
	tree_higgs_2017->Draw("weight_btagsf    >> h_btagsf_light_higgs_2017");
	tree_higgs_2017->Draw("weight_btagsf_light_UP >> h_btagsf_light_higgs_2017_UP");
	tree_higgs_2017->Draw("weight_btagsf_light_DN >> h_btagsf_light_higgs_2017_DN");

	float sign_btagsf_light_higgs_2017 = 1;
	float btagsf_light_higgs_2017_uncert = max (fabs(h_btagsf_light_higgs_2017_UP->GetMean() - h_btagsf_light_higgs_2017->GetMean()), fabs( h_btagsf_light_higgs_2017->GetMean() - h_btagsf_light_higgs_2017_DN->GetMean() ) );

	if( fabs(h_btagsf_light_higgs_2017_UP->GetMean() - h_btagsf_light_higgs_2017->GetMean()) >  fabs(  h_btagsf_light_higgs_2017->GetMean() - h_btagsf_light_higgs_2017_DN->GetMean() ))
	  sign_btagsf_light_higgs_2017 = ((h_btagsf_light_higgs_2017_UP->GetMean() - h_btagsf_light_higgs_2017->GetMean())  > 0 ) - ((h_btagsf_light_higgs_2017_UP->GetMean() - h_btagsf_light_higgs_2017->GetMean())  < 0) ;
	else
	  sign_btagsf_light_higgs_2017 = ((h_btagsf_light_higgs_2017->GetMean() - h_btagsf_light_higgs_2017_DN->GetMean())  > 0 ) - ((h_btagsf_light_higgs_2017->GetMean() - h_btagsf_light_higgs_2017_DN->GetMean())  < 0) ;
	  
	btagsf_light_higgs_2017 << regionSaveName << " \t " <<    sign_btagsf_light_higgs_2017* btagsf_light_higgs_2017_uncert  <<std::endl;

	delete h_btagsf_light_higgs_2017; delete h_btagsf_light_higgs_2017_UP; delete h_btagsf_light_higgs_2017_DN;


	/////////////////scale variation
	vector< float > scaleUncert;
	for( int i=0; i<110; i++){
	  TH1D* h_scalesTemp_higgs_2017    = new TH1D("h_scalesTemp_higgs_2017",    "", 100, 0, 2);
	  tree_higgs_2017->Draw(Form("weight_scales_%d    >> h_scalesTemp_higgs_2017", i));
	  float newValue = h_scalesTemp_higgs_2017->GetMean();
	  if( (newValue > 1.4 ) || (newValue < 0.6 ))
	    newValue = 1.;  // savety first, in case of crazy big diff, it will thus take a smaller diff
	  scaleUncert.push_back( newValue );

	  delete h_scalesTemp_higgs_2017;
	}
	double scaleMax = *max_element(scaleUncert.begin(), scaleUncert.end()-100);
	double scaleMin = *min_element(scaleUncert.begin(), scaleUncert.end()-100);
	double pdfMax   = *max_element(scaleUncert.begin()+10, scaleUncert.end());
	double pdfMin   = *min_element(scaleUncert.begin()+10, scaleUncert.end());
	float scaleVar_higgs_2017_uncert = ((scaleMax -1.)>(1.-scaleMin)) ?  scaleMax : scaleMin;
	float pdfVar_higgs_2017_uncert = ((pdfMax -1.)>(1.-pdfMin)) ?  pdfMax : pdfMin;
	scaleVar_higgs_2017 << regionSaveName << " \t " <<    scaleVar_higgs_2017_uncert  <<std::endl;
	pdfVar_higgs_2017   << regionSaveName << " \t " <<    pdfVar_higgs_2017_uncert  <<std::endl;



      } //end 2017



      //bg    bg_ds.append( qcd_ds );   //append bg data into one
      //bg    bg_ds.append( gjets_ds ); //append bg data into one

      RooDataSet higgs_ds_130(  higgs_ds, Form("higgs_2016_130_13TeV_%s", regionSaveName.c_str() ) );
      //bg    RooDataSet bg_ds_130   (  bg_ds,    Form("bg_2016_130_13TeV_%s",    regionSaveName.c_str() ) );
      RooDataSet data_ds_130 (  data_ds,  Form("data_130_13TeV_%s",  regionSaveName.c_str() ) );

      RooDataSet higgs_ds_125(  higgs_ds, Form("higgs_2016_125_13TeV_%s", regionSaveName.c_str() ) );
      //bg    RooDataSet bg_ds_125   (  bg_ds,    Form("bg_2016_125_13TeV_%s",    regionSaveName.c_str() ) );
      RooDataSet data_ds_125 (  data_ds,  Form("data_125_13TeV_%s",  regionSaveName.c_str() ) );

      RooDataSet higgs_ds_120(  higgs_ds, Form("higgs_2016_120_13TeV_%s", regionSaveName.c_str() ) );
      //bg    RooDataSet bg_ds_120   (  bg_ds,    Form("bg_2016_120_13TeV_%s",    regionSaveName.c_str() ) );
      RooDataSet data_ds_120 (  data_ds,  Form("data_120_13TeV_%s",  regionSaveName.c_str() ) );

      ((RooRealVar*)higgs_ds_130.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
      ((RooRealVar*)data_ds_130.addColumn( hgg_mass_130 ))->setRange(100.,180.); 

      ((RooRealVar*)higgs_ds_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)data_ds_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 

      ((RooRealVar*)higgs_ds_120.addColumn( hgg_mass_120 ))->setRange(100.,180.);
      ((RooRealVar*)data_ds_120.addColumn( hgg_mass_120 ))->setRange(100.,180.);

      //    RooDataSet diH_ds_125(  diH_ds, Form("diH_2016_125_13TeV_%s", regionSaveName.c_str() ) );
      //    ((RooRealVar*)diH_ds_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 

      //    RooDataSet diH_ds_2017(Form("diH_2017_125_13TeV_%s", regionSaveName.c_str() ), "DiH", vars, WeightVar(weight), Import(*tree_diH_2017) );
      //    RooDataSet diH_ds_125_2017(  diH_ds_2017, Form("diH_2017_125_13TeV_%s", regionSaveName.c_str() ) );
      //    ((RooRealVar*)diH_ds_125_2017.addColumn( hgg_mass_125 ))->setRange(100.,180.); 

      RooDataSet higgsGluGlu_ds_125(  higgsGluGlu_ds, Form("higgsGluGlu_2016_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsGluGlu_ds_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      RooDataSet higgsVBF_ds_125(  higgsVBF_ds, Form("higgsVBF_2016_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsVBF_ds_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      RooDataSet higgsVH_ds_125(  higgsVH_ds, Form("higgsVH_2016_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsVH_ds_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      RooDataSet higgsttH_ds_125(  higgsttH_ds, Form("higgsttH_2016_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsttH_ds_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      RooDataSet higgsbbH_ds_125(  higgsbbH_ds, Form("higgsbbH_2016_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsbbH_ds_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      // RooDataSet higgsTHX_ds_125(  higgsTHX_ds, Form("higgsTHX_2016_125_13TeV_%s", regionSaveName.c_str() ) );
      // ((RooRealVar*)higgsTHX_ds_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
    
      RooDataSet higgsGluGlu_ds_2017_125(  higgsGluGlu_ds_2017, Form("higgsGluGlu_2017_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsGluGlu_ds_2017_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      RooDataSet higgsVBF_ds_2017_125(  higgsVBF_ds_2017, Form("higgsVBF_2017_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsVBF_ds_2017_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      RooDataSet higgsVH_ds_2017_125(  higgsVH_ds_2017, Form("higgsVH_2017_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsVH_ds_2017_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      RooDataSet higgsttH_ds_2017_125(  higgsttH_ds_2017, Form("higgsttH_2017_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsttH_ds_2017_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      RooDataSet higgsbbH_ds_2017_125(  higgsbbH_ds_2017, Form("higgsbbH_2017_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgsbbH_ds_2017_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
      // RooDataSet higgsTHX_ds_2017_125(  higgsTHX_ds_2017, Form("higgsTHX_2017_125_13TeV_%s", regionSaveName.c_str() ) );
      // ((RooRealVar*)higgsTHX_ds_2017_125.addColumn( hgg_mass_125 ))->setRange(122.,129.); 
 

      //Import the 2016 versions
      ws_data.import( data_ds_125  );
      //bg    ws_bg  .import( bg_ds_125    );
      ws_sig .import( higgs_ds_125 );

      ws_data.import( data_ds_130  );
      //bg    ws_bg  .import( bg_ds_130    );
      ws_sig .import( higgs_ds_130 );

      ws_data.import( data_ds_120  );
      //bg    ws_bg  .import( bg_ds_120    );
      ws_sig .import( higgs_ds_120 );

      //write number of data events per region into txt file

      RooDataSet higgs_ds_2016_JECup(Form("higgs_2016_JECup_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_JECup) );
      RooDataSet higgs_ds_2016_JECdn(Form("higgs_2016_JECdn_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_JECdn) );
      RooDataSet higgs_ds_125_2016_JECup(  higgs_ds_2016_JECup, Form("higgs_2016_JECup_125_13TeV_%s", regionSaveName.c_str() ) );
      RooDataSet higgs_ds_125_2016_JECdn(  higgs_ds_2016_JECdn, Form("higgs_2016_JECdn_125_13TeV_%s", regionSaveName.c_str() ) );
      ((RooRealVar*)higgs_ds_125_2016_JECup.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      ((RooRealVar*)higgs_ds_125_2016_JECdn.addColumn( hgg_mass_125 ))->setRange(100.,180.); 

      // std::cout  << regionSaveName << "    " <<  higgs_ds_125.sumEntries() << std::endl;
      // std::cout  << regionSaveName << "    " <<  higgs_ds_125_2016_JECup.sumEntries() << std::endl;
      // std::cout  << regionSaveName << "  up  " <<  float(( higgs_ds_125_2016_JECup.sumEntries() - higgs_ds_125.sumEntries() ))/float(higgs_ds_125.sumEntries() )   << std::endl;
      // std::cout  << regionSaveName << "  up  " <<  float(higgs_ds_125.sumEntries() )   << std::endl;
      // std::cout  << regionSaveName << "  up  " <<  float(( higgs_ds_125_2016_JECup.sumEntries() - higgs_ds_125.sumEntries() ))    << std::endl;
      // std::cout  << regionSaveName << "  dn  " <<  ((higgs_ds_125.sumEntries() - higgs_ds_125_2016_JECdn.sumEntries()))  << std::endl;

    
      higgsSyst_2016_JECup << regionSaveName << "    " <<  float(( higgs_ds_125_2016_JECup.sumEntries() - higgs_ds_125.sumEntries() ))/float(higgs_ds_125.sumEntries() )  << std::endl;
      higgsSyst_2016_JECdn << regionSaveName << "    " <<  float((higgs_ds_125.sumEntries() - higgs_ds_125_2016_JECdn.sumEntries() ) )/float(higgs_ds_125.sumEntries() )  << std::endl;


      dataAmount << regionSaveName << "    " << data_ds_125.sumEntries() << std::endl;
      //bg    bgSumAmount << regionSaveName << "    " << 35.9* bg_ds_125.sumEntries() << std::endl; //not ideal way of scaling the bg actually

      higgsAmount << regionSaveName << "    " <<  higgs_ds_125.sumEntries() << std::endl;

      higgsGluGlu << regionSaveName << "    " <<  higgsGluGlu_ds_125.sumEntries()*lumi2016 +  higgsGluGlu_ds_2017_125.sumEntries()*lumi2017 << std::endl;
      higgsVBF    << regionSaveName << "    " <<  higgsVBF_ds_125.sumEntries()*lumi2016    +  higgsVBF_ds_2017_125.sumEntries()*lumi2017 << std::endl;
      higgsVH     << regionSaveName << "    " <<  higgsVH_ds_125.sumEntries()*lumi2016     +  higgsVH_ds_2017_125.sumEntries()*lumi2017 << std::endl;
      higgsttH    << regionSaveName << "    " <<  higgsttH_ds_125.sumEntries()*lumi2016    +  higgsttH_ds_2017_125.sumEntries()*lumi2017 << std::endl;
      higgsbbH    << regionSaveName << "    " <<  higgsbbH_ds_125.sumEntries()*lumi2016    +  higgsbbH_ds_2017_125.sumEntries()*lumi2017 << std::endl;
      //    higgsTHX    << regionSaveName << "    " <<  higgsTHX_ds_125.sumEntries()*lumi2016    +  higgsTHX_ds_2017_125.sumEntries()*lumi2017 << std::endl;

      //    higgsdiH    << regionSaveName << "    " <<  diH_ds_125.sumEntries()*lumi2016    +  diH_ds_125_2017.sumEntries()*lumi2017 << std::endl;

      //////////////////////////////////////////////////////////////////////////////
      //////////////////  UNCERTS FOR HIGGS ////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////

      TH1D* h_gammasf_higgs    = new TH1D("h_gammasf_higgs",    "", 100, 0, 2);
      TH1D* h_gammasf_higgs_UP = new TH1D("h_gammasf_higgs_UP", "", 100, 0, 2);
      TH1D* h_gammasf_higgs_DN = new TH1D("h_gammasf_higgs_DN", "", 100, 0, 2);
 
      tree_higgs->Draw("weight_gammasf    >> h_gammasf_higgs");
      tree_higgs->Draw("weight_gammasf_UP >> h_gammasf_higgs_UP");
      tree_higgs->Draw("weight_gammasf_DN >> h_gammasf_higgs_DN");

      float sign_gammasf_higgs = 1;
      float gammasf_higgs_uncert = max (fabs(h_gammasf_higgs_UP->GetMean() - h_gammasf_higgs->GetMean()), fabs( h_gammasf_higgs->GetMean() - h_gammasf_higgs_DN->GetMean() ) );

      if( fabs(h_gammasf_higgs_UP->GetMean() - h_gammasf_higgs->GetMean()) >  fabs(  h_gammasf_higgs->GetMean() - h_gammasf_higgs_DN->GetMean() ))
	sign_gammasf_higgs = ((h_gammasf_higgs_UP->GetMean() - h_gammasf_higgs->GetMean())  > 0 ) - ((h_gammasf_higgs_UP->GetMean() - h_gammasf_higgs->GetMean())  < 0) ;
      else
	sign_gammasf_higgs = ((h_gammasf_higgs->GetMean() - h_gammasf_higgs_DN->GetMean())  > 0 ) - ((h_gammasf_higgs->GetMean() - h_gammasf_higgs_DN->GetMean())  < 0) ;
	  
      gammasf_higgs << regionSaveName << " \t " <<    sign_gammasf_higgs* gammasf_higgs_uncert  <<std::endl;

      delete h_gammasf_higgs; delete h_gammasf_higgs_UP; delete h_gammasf_higgs_DN;






      //////////////////////leptons //////////////////////////////////
      TH1D* h_lepsf_higgs    = new TH1D("h_lepsf_higgs",    "", 100, 0, 2);
      TH1D* h_lepsf_higgs_UP = new TH1D("h_lepsf_higgs_UP", "", 100, 0, 2);
      TH1D* h_lepsf_higgs_DN = new TH1D("h_lepsf_higgs_DN", "", 100, 0, 2);
 
      tree_higgs->Draw("weight_lepsf    >> h_lepsf_higgs");
      tree_higgs->Draw("weight_lepsf_UP >> h_lepsf_higgs_UP");
      tree_higgs->Draw("weight_lepsf_DN >> h_lepsf_higgs_DN");

      float sign_lepsf_higgs = 1;
      float lepsf_higgs_uncert = max (fabs(h_lepsf_higgs_UP->GetMean() - h_lepsf_higgs->GetMean()), fabs( h_lepsf_higgs->GetMean() - h_lepsf_higgs_DN->GetMean() ) );

      if( fabs(h_lepsf_higgs_UP->GetMean() - h_lepsf_higgs->GetMean()) >  fabs(  h_lepsf_higgs->GetMean() - h_lepsf_higgs_DN->GetMean() ))
	sign_lepsf_higgs = ((h_lepsf_higgs_UP->GetMean() - h_lepsf_higgs->GetMean())  > 0 ) - ((h_lepsf_higgs_UP->GetMean() - h_lepsf_higgs->GetMean())  < 0) ;
      else
	sign_lepsf_higgs = ((h_lepsf_higgs->GetMean() - h_lepsf_higgs_DN->GetMean())  > 0 ) - ((h_lepsf_higgs->GetMean() - h_lepsf_higgs_DN->GetMean())  < 0) ;
	  
      lepsf_higgs << regionSaveName << " \t " <<    sign_lepsf_higgs* lepsf_higgs_uncert  <<std::endl;

      delete h_lepsf_higgs; delete h_lepsf_higgs_UP; delete h_lepsf_higgs_DN;

      
 


      //btag
      TH1D* h_btagsf_heavy_higgs    = new TH1D("h_btagsf_heavy_higgs",    "", 100, 0, 2);
      TH1D* h_btagsf_heavy_higgs_UP = new TH1D("h_btagsf_heavy_higgs_UP", "", 100, 0, 2);
      TH1D* h_btagsf_heavy_higgs_DN = new TH1D("h_btagsf_heavy_higgs_DN", "", 100, 0, 2);
 
      tree_higgs->Draw("weight_btagsf    >> h_btagsf_heavy_higgs");
      tree_higgs->Draw("weight_btagsf_heavy_UP >> h_btagsf_heavy_higgs_UP");
      tree_higgs->Draw("weight_btagsf_heavy_DN >> h_btagsf_heavy_higgs_DN");

      float sign_btagsf_heavy_higgs = 1;
      float btagsf_heavy_higgs_uncert = max (fabs(h_btagsf_heavy_higgs_UP->GetMean() - h_btagsf_heavy_higgs->GetMean()), fabs( h_btagsf_heavy_higgs->GetMean() - h_btagsf_heavy_higgs_DN->GetMean() ) );

      if( fabs(h_btagsf_heavy_higgs_UP->GetMean() - h_btagsf_heavy_higgs->GetMean()) >  fabs(  h_btagsf_heavy_higgs->GetMean() - h_btagsf_heavy_higgs_DN->GetMean() ))
	sign_btagsf_heavy_higgs = ((h_btagsf_heavy_higgs_UP->GetMean() - h_btagsf_heavy_higgs->GetMean())  > 0 ) - ((h_btagsf_heavy_higgs_UP->GetMean() - h_btagsf_heavy_higgs->GetMean())  < 0) ;
      else
	sign_btagsf_heavy_higgs = ((h_btagsf_heavy_higgs->GetMean() - h_btagsf_heavy_higgs_DN->GetMean())  > 0 ) - ((h_btagsf_heavy_higgs->GetMean() - h_btagsf_heavy_higgs_DN->GetMean())  < 0) ;
	  
      btagsf_heavy_higgs << regionSaveName << " \t " <<    sign_btagsf_heavy_higgs* btagsf_heavy_higgs_uncert  <<std::endl;

      delete h_btagsf_heavy_higgs; delete h_btagsf_heavy_higgs_UP; delete h_btagsf_heavy_higgs_DN;


      TH1D* h_btagsf_light_higgs    = new TH1D("h_btagsf_light_higgs",    "", 100, 0, 2);
      TH1D* h_btagsf_light_higgs_UP = new TH1D("h_btagsf_light_higgs_UP", "", 100, 0, 2);
      TH1D* h_btagsf_light_higgs_DN = new TH1D("h_btagsf_light_higgs_DN", "", 100, 0, 2);
 
      tree_higgs->Draw("weight_btagsf    >> h_btagsf_light_higgs");
      tree_higgs->Draw("weight_btagsf_light_UP >> h_btagsf_light_higgs_UP");
      tree_higgs->Draw("weight_btagsf_light_DN >> h_btagsf_light_higgs_DN");

      float sign_btagsf_light_higgs = 1;
      float btagsf_light_higgs_uncert = max (fabs(h_btagsf_light_higgs_UP->GetMean() - h_btagsf_light_higgs->GetMean()), fabs( h_btagsf_light_higgs->GetMean() - h_btagsf_light_higgs_DN->GetMean() ) );

      if( fabs(h_btagsf_light_higgs_UP->GetMean() - h_btagsf_light_higgs->GetMean()) >  fabs(  h_btagsf_light_higgs->GetMean() - h_btagsf_light_higgs_DN->GetMean() ))
	sign_btagsf_light_higgs = ((h_btagsf_light_higgs_UP->GetMean() - h_btagsf_light_higgs->GetMean())  > 0 ) - ((h_btagsf_light_higgs_UP->GetMean() - h_btagsf_light_higgs->GetMean())  < 0) ;
      else
	sign_btagsf_light_higgs = ((h_btagsf_light_higgs->GetMean() - h_btagsf_light_higgs_DN->GetMean())  > 0 ) - ((h_btagsf_light_higgs->GetMean() - h_btagsf_light_higgs_DN->GetMean())  < 0) ;
	  
      btagsf_light_higgs << regionSaveName << " \t " <<    sign_btagsf_light_higgs* btagsf_light_higgs_uncert  <<std::endl;

      delete h_btagsf_light_higgs; delete h_btagsf_light_higgs_UP; delete h_btagsf_light_higgs_DN;


      /////////////////scale variation
      vector< float > scaleUncert;
      for( int i=0; i<110; i++){
	TH1D* h_scalesTemp_higgs    = new TH1D("h_scalesTemp_higgs",    "", 100, 0, 2);
	tree_higgs->Draw(Form("weight_scales_%d    >> h_scalesTemp_higgs", i));
	float newValue = h_scalesTemp_higgs->GetMean();
	if( (newValue > 1.4 ) || (newValue < 0.6 ))
	  newValue = 1.;  // savety first, in case of crazy big diff, it will thus take a diff variation
	scaleUncert.push_back( newValue );
	//scaleUncert.push_back( h_scalesTemp_higgs->GetMean() );
	delete h_scalesTemp_higgs;
      }
      double scaleMax = *max_element(scaleUncert.begin(), scaleUncert.end()-100);
      double scaleMin = *min_element(scaleUncert.begin(), scaleUncert.end()-100);
      double pdfMax   = *max_element(scaleUncert.begin()+10, scaleUncert.end());
      double pdfMin   = *min_element(scaleUncert.begin()+10, scaleUncert.end());
      float scaleVar_higgs_uncert = ((scaleMax -1.)>(1.-scaleMin)) ?  scaleMax : scaleMin;
      float pdfVar_higgs_uncert = ((pdfMax -1.)>(1.-pdfMin)) ?  pdfMax : pdfMin;
      scaleVar_higgs << regionSaveName << " \t " <<    scaleVar_higgs_uncert  <<std::endl;
      pdfVar_higgs   << regionSaveName << " \t " <<    pdfVar_higgs_uncert  <<std::endl;

    }//if SMHDonlu




    if( doSMHDonly == 0  ){

    ////////////////////////////////////////////////////////////////////////
    //NOW for the SIGNAL LIST
    ////////////////////////////////////////////////////////////////////////
      //    for( unsigned int i=0; i<signalList.size(); i++){
    for( unsigned int i=0; i<signals.size(); i++){


      // if( (i>= sigPerPart * partNr) || (i < ((partNr-1)*sigPerPart ))   )
      // 	continue;

      std::cout << "Working on signal nr " << i << std::endl;

      TTree* tree_sig_ini = signals[i]->get(thisRegion)->tree;  
      TTree* tree_sig = tree_sig_ini->CopyTree(treeSel.c_str() );
      RooDataSet sig_ds(Form("%s_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig) ) ;

      TTree* tree_sig_genMET_ini = signals[i]->get(thisRegion)->tree;  
      TTree* tree_sig_genMET = tree_sig_genMET_ini->CopyTree(treeSel_genMET.c_str() );
      RooDataSet sig_genMET_ds(Form("%s_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_genMET) ) ;

      TTree* tree_sig_ini_JECup = signals_JECup[i]->get(thisRegion)->tree;  
      TTree* tree_sig_JECup = tree_sig_ini_JECup->CopyTree(treeSel.c_str() );
      RooDataSet sig_ds_2016_JECup(Form("%s_JECup_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_JECup) ) ;

      TTree* tree_sig_ini_JECdn = signals_JECdn[i]->get(thisRegion)->tree;  
      TTree* tree_sig_JECdn = tree_sig_ini_JECdn->CopyTree(treeSel.c_str() );
      RooDataSet sig_ds_2016_JECdn(Form("%s_JECdn_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_JECdn) ) ;




      if(doHZ){
	std::cout << "ADDING HH to HZ" << std::endl;

	TTree* tree_sig_ini_HH_HH0p25 = signals_HH_HH0p25[i]->get(thisRegion)->tree;  
	TTree* tree_sig_HH_HH0p25 = tree_sig_ini_HH_HH0p25->CopyTree(treeSel.c_str() );
	RooDataSet sig_ds_HH_HH0p25(Form("%s_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25) ) ;
	sig_ds.append( sig_ds_HH_HH0p25 );

	TTree* tree_sig_genMET_ini_HH_HH0p25 = signals_HH_HH0p25[i]->get(thisRegion)->tree;  
	TTree* tree_sig_genMET_HH_HH0p25 = tree_sig_genMET_ini_HH_HH0p25->CopyTree(treeSel_genMET.c_str() );
	RooDataSet sig_genMET_ds_HH_HH0p25(Form("%s_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_genMET_HH_HH0p25) ) ;
	sig_genMET_ds.append( sig_genMET_ds_HH_HH0p25 );

	TTree* tree_sig_ini_HH_HH0p25_JECup = signals_HH_HH0p25_JECup[i]->get(thisRegion)->tree;  
	TTree* tree_sig_HH_HH0p25_JECup = tree_sig_ini_HH_HH0p25_JECup->CopyTree(treeSel.c_str() );
	RooDataSet sig_ds_HH_HH0p25_JECup(Form("%s_JECup_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_JECup) ) ;
	sig_ds_2016_JECup.append( sig_ds_HH_HH0p25_JECup );

	TTree* tree_sig_ini_HH_HH0p25_JECdn = signals_HH_HH0p25_JECdn[i]->get(thisRegion)->tree;  
	TTree* tree_sig_HH_HH0p25_JECdn = tree_sig_ini_HH_HH0p25_JECdn->CopyTree(treeSel.c_str() );
	RooDataSet sig_ds_HH_HH0p25_JECdn(Form("%s_JECdn_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_JECdn) ) ;
	sig_ds_2016_JECdn.append( sig_ds_HH_HH0p25_JECdn );
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

	TString sigSampleName(signalList[i]);

	std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	std::string::size_type sz;   // alias of size_t

	if( !sigSampleName.Contains("TChi")  )
	  massChi = std::stoi (mLSP,&sz);
	else if( sigSampleName.Contains("TChiH")  )
	  massChi = 1;
	else if( sigSampleName.Contains("TChiWH")  ){ 
	  mLSP = signalList[i].substr( signalList[i].find("_HToGG_m")+13);
	  std::cout << "m2 " << mLSP << std::endl;
	  massChi = std::stoi (mLSP,&sz);
	}

	std::cout << "mlsp " << mLSP << " mass   " << massChi <<std::endl;

	std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	if( sigSampleName.Contains("TChiH")  ){
	  std::cout << "getting EW mass names out of signal name" << std::endl;
	  mSbottom = signalList[i].substr( signalList[i].find("_m")+2);
	}
	if( sigSampleName.Contains("TChiWH")  ){

	  std::cout << "getting WH EW mass names out of signal name" << std::endl;
	  mSbottom = signalList[i].substr( signalList[i].find("_HToGG_m")+8);
	  std::cout << " names   " << mSbottom << std::endl;
	}

	massSbottom = std::stoi (mSbottom,&sz);
	std::cout << "bottom " << mSbottom << " mass   " << massSbottom <<std::endl;

	std::cout << "mass x " << massSbottom << " mass y   " << massChi <<std::endl;

	int binBottom = h2D_isr->GetXaxis()->FindBin( (float)massSbottom );
	int binChi    = h2D_isr->GetYaxis()->FindBin( (float)massChi );

	h2D_isr->SetBinContent( binBottom, binChi, h_isr->GetMean() );
	h2D_isr_err->SetBinContent( binBottom, binChi, sign* isr_uncert );

	delete h_isr; delete h_isr_UP; delete h_isr_DN;

	///////////////////////////// GAMMA SF///////////////////////////////////////////////////
	TH1D* h_gammasf    = new TH1D("h_gammasf",    "", 100, 0, 2);
	TH1D* h_gammasf_UP = new TH1D("h_gammasf_UP", "", 100, 0, 2);
	TH1D* h_gammasf_DN = new TH1D("h_gammasf_DN", "", 100, 0, 2);
 
	tree_sig->Draw("weight_gammasf    >> h_gammasf");
	tree_sig->Draw("weight_gammasf_UP >> h_gammasf_UP");
	tree_sig->Draw("weight_gammasf_DN >> h_gammasf_DN");

	float sign_gammasf = 1;
	float gammasf_uncert = max (fabs(h_gammasf_UP->GetMean() - h_gammasf->GetMean()), fabs( h_gammasf->GetMean() - h_gammasf_DN->GetMean() ) );

	if( fabs(h_gammasf_UP->GetMean() - h_gammasf->GetMean()) >  fabs(  h_gammasf->GetMean() - h_gammasf_DN->GetMean() ))
	  sign_gammasf = ((h_gammasf_UP->GetMean() - h_gammasf->GetMean())  > 0 ) - ((h_gammasf_UP->GetMean() - h_gammasf->GetMean())  < 0) ;
	else
	  sign_gammasf = ((h_gammasf->GetMean() - h_gammasf_DN->GetMean())  > 0 ) - ((h_gammasf->GetMean() - h_gammasf_DN->GetMean())  < 0) ;
	  
	//	gammasf <<  signalList[i] << " \t " << regionSaveName << " \t " <<    sign_gammasf* gammasf_uncert  <<std::endl;

	h2D_gammasf->SetBinContent( binBottom, binChi, h_gammasf->GetMean() );
	h2D_gammasf_err->SetBinContent( binBottom, binChi, sign_gammasf* gammasf_uncert );

	delete h_gammasf; delete h_gammasf_UP; delete h_gammasf_DN;



	///////////////////////////// LEPTON SF///////////////////////////////////////////////////
	TH1D* h_lepsf    = new TH1D("h_lepsf",    "", 100, 0, 2);
	TH1D* h_lepsf_UP = new TH1D("h_lepsf_UP", "", 100, 0, 2);
	TH1D* h_lepsf_DN = new TH1D("h_lepsf_DN", "", 100, 0, 2);
 
	tree_sig->Draw("weight_lepsf    >> h_lepsf");
	tree_sig->Draw("weight_lepsf_UP >> h_lepsf_UP");
	tree_sig->Draw("weight_lepsf_DN >> h_lepsf_DN");

	float sign_lepsf = 1;
	float lepsf_uncert = max (fabs(h_lepsf_UP->GetMean() - h_lepsf->GetMean()), fabs( h_lepsf->GetMean() - h_lepsf_DN->GetMean() ) );

	if( fabs(h_lepsf_UP->GetMean() - h_lepsf->GetMean()) >  fabs(  h_lepsf->GetMean() - h_lepsf_DN->GetMean() ))
	  sign_lepsf = ((h_lepsf_UP->GetMean() - h_lepsf->GetMean())  > 0 ) - ((h_lepsf_UP->GetMean() - h_lepsf->GetMean())  < 0) ;
	else
	  sign_lepsf = ((h_lepsf->GetMean() - h_lepsf_DN->GetMean())  > 0 ) - ((h_lepsf->GetMean() - h_lepsf_DN->GetMean())  < 0) ;
	  
	//	lepsf <<  signalList[i] << " \t " << regionSaveName << " \t " <<    sign_lepsf* lepsf_uncert  <<std::endl;

	h2D_lepsf->SetBinContent( binBottom, binChi, h_lepsf->GetMean() );
	h2D_lepsf_err->SetBinContent( binBottom, binChi, sign_lepsf* lepsf_uncert );

	delete h_lepsf; delete h_lepsf_UP; delete h_lepsf_DN;


	///////////////////////////// BTAG SF///////////////////////////////////////////////////
	TH1D* h_btagsf_heavy    = new TH1D("h_btagsf_heavy",    "", 100, 0, 2);
	TH1D* h_btagsf_heavy_UP = new TH1D("h_btagsf_heavy_UP", "", 100, 0, 2);
	TH1D* h_btagsf_heavy_DN = new TH1D("h_btagsf_heavy_DN", "", 100, 0, 2);
 
	tree_sig->Draw("weight_btagsf    >> h_btagsf_heavy");
	tree_sig->Draw("weight_btagsf_heavy_UP >> h_btagsf_heavy_UP");
	tree_sig->Draw("weight_btagsf_heavy_DN >> h_btagsf_heavy_DN");

	float sign_btagsf_heavy = 1;
	float btagsf_heavy_uncert = max (fabs(h_btagsf_heavy_UP->GetMean() - h_btagsf_heavy->GetMean()), fabs( h_btagsf_heavy->GetMean() - h_btagsf_heavy_DN->GetMean() ) );

	if( fabs(h_btagsf_heavy_UP->GetMean() - h_btagsf_heavy->GetMean()) >  fabs(  h_btagsf_heavy->GetMean() - h_btagsf_heavy_DN->GetMean() ))
	  sign_btagsf_heavy = ((h_btagsf_heavy_UP->GetMean() - h_btagsf_heavy->GetMean())  > 0 ) - ((h_btagsf_heavy_UP->GetMean() - h_btagsf_heavy->GetMean())  < 0) ;
	else
	  sign_btagsf_heavy = ((h_btagsf_heavy->GetMean() - h_btagsf_heavy_DN->GetMean())  > 0 ) - ((h_btagsf_heavy->GetMean() - h_btagsf_heavy_DN->GetMean())  < 0) ;
	  
	//	btagsf_heavy <<  signalList[i] << " \t " << regionSaveName << " \t " <<    sign_btagsf_heavy* btagsf_heavy_uncert  <<std::endl;

	h2D_btagsf_heavy->SetBinContent( binBottom, binChi, h_btagsf_heavy->GetMean() );
	h2D_btagsf_heavy_err->SetBinContent( binBottom, binChi, sign_btagsf_heavy* btagsf_heavy_uncert );

	delete h_btagsf_heavy; delete h_btagsf_heavy_UP; delete h_btagsf_heavy_DN;

	///////////////////////////// BTAG SF///////////////////////////////////////////////////
	TH1D* h_btagsf_light    = new TH1D("h_btagsf_light",    "", 100, 0, 2);
	TH1D* h_btagsf_light_UP = new TH1D("h_btagsf_light_UP", "", 100, 0, 2);
	TH1D* h_btagsf_light_DN = new TH1D("h_btagsf_light_DN", "", 100, 0, 2);
 
	tree_sig->Draw("weight_btagsf    >> h_btagsf_light");
	tree_sig->Draw("weight_btagsf_light_UP >> h_btagsf_light_UP");
	tree_sig->Draw("weight_btagsf_light_DN >> h_btagsf_light_DN");

	float sign_btagsf_light = 1;
	float btagsf_light_uncert = max (fabs(h_btagsf_light_UP->GetMean() - h_btagsf_light->GetMean()), fabs( h_btagsf_light->GetMean() - h_btagsf_light_DN->GetMean() ) );

	if( fabs(h_btagsf_light_UP->GetMean() - h_btagsf_light->GetMean()) >  fabs(  h_btagsf_light->GetMean() - h_btagsf_light_DN->GetMean() ))
	  sign_btagsf_light = ((h_btagsf_light_UP->GetMean() - h_btagsf_light->GetMean())  > 0 ) - ((h_btagsf_light_UP->GetMean() - h_btagsf_light->GetMean())  < 0) ;
	else
	  sign_btagsf_light = ((h_btagsf_light->GetMean() - h_btagsf_light_DN->GetMean())  > 0 ) - ((h_btagsf_light->GetMean() - h_btagsf_light_DN->GetMean())  < 0) ;
	//	btagsf_light <<  signalList[i] << " \t " << regionSaveName << " \t " <<    sign_btagsf_light* btagsf_light_uncert  <<std::endl;
	h2D_btagsf_light->SetBinContent( binBottom, binChi, h_btagsf_light->GetMean() );
	h2D_btagsf_light_err->SetBinContent( binBottom, binChi, sign_btagsf_light* btagsf_light_uncert );
	delete h_btagsf_light; delete h_btagsf_light_UP; delete h_btagsf_light_DN;
 


	/////////////////scale variation SIGNAL 2016
	vector< float > scaleUncert;
	for( int i=0; i<110; i++){
	  TH1D* h_scalesTemp    = new TH1D("h_scalesTemp",    "", 100, 0, 2);
	  tree_sig->Draw(Form("weight_scales_%d    >> h_scalesTemp", i));
	  float newValue = h_scalesTemp->GetMean();
	  if( (newValue > 1.4 ) || (newValue < 0.6 ))
	    newValue = 1.;  // savety first, in case of crazy big diff, it will thus take a diff variation
	  scaleUncert.push_back( newValue );
	  // 	  scaleUncert.push_back( h_scalesTemp->GetMean() );
	  delete h_scalesTemp;
	}
	double scaleMax = *max_element(scaleUncert.begin(), scaleUncert.end()-100);
	double scaleMin = *min_element(scaleUncert.begin(), scaleUncert.end()-100);
	double pdfMax   = *max_element(scaleUncert.begin()+10, scaleUncert.end());
	double pdfMin   = *min_element(scaleUncert.begin()+10, scaleUncert.end());
	float scaleVar_uncert = ((scaleMax -1.)>(1.-scaleMin)) ?  scaleMax : scaleMin;
	float pdfVar_uncert = ((pdfMax -1.)>(1.-pdfMin)) ?  pdfMax : pdfMin;

	h2D_scaleVar->SetBinContent( binBottom, binChi, scaleVar_uncert  );
	h2D_pdfVar  ->SetBinContent( binBottom, binChi, pdfVar_uncert  );

	// if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	//   sigAmount <<  signalList[i] << " \t " << regionSaveName << " \t " <<  sig_ds.sumEntries()  << std::endl;



	///////////////////////////// GEN MET correction ///////////////////////////////////////
	/////////// factor of 1/2 already in the weight, so don't divide by it again! //////////
	float genMET_uncert = ( sig_ds.sumEntries() - sig_genMET_ds.sumEntries()  );

	if( ((sig_ds.sumEntries() + sig_genMET_ds.sumEntries()) == 0 ) || i >=( regionSelection.size()-4)  ) 
	  genMET_uncert = 0.;
	else
	  genMET_uncert	 /= ( sig_ds.sumEntries() + sig_genMET_ds.sumEntries() );

	h2D_genMET->SetBinContent( binBottom, binChi, genMET_uncert );
 

	std::cout << "genMET 2016 syst                                " <<  genMET_uncert << std::endl;



	std::cout << "genMET 2016 yield " << sig_ds.sumEntries()  << std::endl;
	std::cout << "genMET 2016 gen yield " << sig_genMET_ds.sumEntries()  << std::endl;
	std::cout << "sig yield summed  " <<  sig_ds.sumEntries() + sig_genMET_ds.sumEntries() << std::endl;	  

      } 


      // 2017 version ///////////////////////////////////////////////////////////////
      if( doComb ) {
	std::cout << "Working on 2017 signal nr " << i << std::endl;

	TTree* tree_sig_ini_2017 = signals_2017[i]->get(thisRegion)->tree;  
	TTree* tree_sig_2017 = tree_sig_ini_2017->CopyTree(treeSel.c_str() );
	RooDataSet sig_ds_2017(Form("%s_2017_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_2017) ) ;

	TTree* tree_sig_genMET_ini_2017 = signals_2017[i]->get(thisRegion)->tree;  
	TTree* tree_sig_genMET_2017 = tree_sig_genMET_ini_2017->CopyTree(treeSel_genMET.c_str() );
	RooDataSet sig_genMET_ds_2017(Form("%s_2017_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_genMET_2017) ) ;

	TTree* tree_sig_ini_JECup_2017 = signals_2017_JECup[i]->get(thisRegion)->tree;  
	TTree* tree_sig_JECup_2017 = tree_sig_ini_JECup_2017->CopyTree(treeSel.c_str() );
	RooDataSet sig_ds_2017_JECup(Form("%s_JECup_2017_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_JECup_2017) ) ;

	TTree* tree_sig_ini_JECdn_2017 = signals_2017_JECdn[i]->get(thisRegion)->tree;  
	TTree* tree_sig_JECdn_2017 = tree_sig_ini_JECdn_2017->CopyTree(treeSel.c_str() );
	RooDataSet sig_ds_2017_JECdn(Form("%s_JECdn_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_JECdn_2017) ) ;

	if(doHZ){
	  std::cout << "ADDING HH to HZ" << std::endl;

	  TTree* tree_sig_ini_HH_HH0p25_2017 = signals_HH_HH0p25_2017[i]->get(thisRegion)->tree;  
	  TTree* tree_sig_HH_HH0p25_2017 = tree_sig_ini_HH_HH0p25_2017->CopyTree(treeSel.c_str() );
	  RooDataSet sig_ds_HH_HH0p25_2017(Form("%s_2017_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_2017) ) ;
	  sig_ds_2017.append( sig_ds_HH_HH0p25_2017 );

	  TTree* tree_sig_genMET_ini_HH_HH0p25_2017 = signals_HH_HH0p25_2017[i]->get(thisRegion)->tree;  
	  TTree* tree_sig_genMET_HH_HH0p25_2017 = tree_sig_genMET_ini_HH_HH0p25_2017->CopyTree(treeSel_genMET.c_str() );
	  RooDataSet sig_genMET_ds_HH_HH0p25_2017(Form("%s_2017_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_genMET_HH_HH0p25_2017) ) ;
	  sig_genMET_ds_2017.append( sig_genMET_ds_HH_HH0p25_2017 );

	  TTree* tree_sig_ini_HH_HH0p25_JECup = signals_HH_HH0p25_JECup[i]->get(thisRegion)->tree;  
	  TTree* tree_sig_HH_HH0p25_JECup = tree_sig_ini_HH_HH0p25_JECup->CopyTree(treeSel.c_str() );
	  RooDataSet sig_ds_HH_HH0p25_JECup(Form("%s_JECup_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_JECup) ) ;
	  sig_ds_2017_JECup.append( sig_ds_HH_HH0p25_JECup );

	  TTree* tree_sig_ini_HH_HH0p25_JECdn = signals_HH_HH0p25_JECdn[i]->get(thisRegion)->tree;  
	  TTree* tree_sig_HH_HH0p25_JECdn = tree_sig_ini_HH_HH0p25_JECdn->CopyTree(treeSel.c_str() );
	  RooDataSet sig_ds_HH_HH0p25_JECdn(Form("%s_JECdn_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_JECdn) ) ;
	  sig_ds_2017_JECdn.append( sig_ds_HH_HH0p25_JECdn );
	}


	if( doISRsyst ){
	  TH1D* h_isr_2017    = new TH1D("h_isr_2017",    "", 100, 0, 2);
	  TH1D* h_isr_2017_UP = new TH1D("h_isr_2017_UP", "", 100, 0, 2);
	  TH1D* h_isr_2017_DN = new TH1D("h_isr_2017_DN", "", 100, 0, 2);
 
	  tree_sig_2017->Draw("weight_isr    >> h_isr_2017");
	  tree_sig_2017->Draw("weight_isr_UP >> h_isr_2017_UP");
	  tree_sig_2017->Draw("weight_isr_DN >> h_isr_2017_DN");

	  float sign = 1;
	  float isr_uncert = max (fabs(h_isr_2017_UP->GetMean() - h_isr_2017->GetMean()), fabs( h_isr_2017->GetMean() - h_isr_2017_DN->GetMean() ) );

	  if( fabs(h_isr_2017_UP->GetMean() - h_isr_2017->GetMean()) >  fabs(  h_isr_2017->GetMean() - h_isr_2017_DN->GetMean() ))
	    sign = ((h_isr_2017_UP->GetMean() - h_isr_2017->GetMean())  > 0 ) - ((h_isr_2017_UP->GetMean() - h_isr_2017->GetMean())  < 0) ;
	  else
	    sign = ((h_isr_2017->GetMean() - h_isr_2017_DN->GetMean())  > 0 ) - ((h_isr_2017->GetMean() - h_isr_2017_DN->GetMean())  < 0) ;
	  
	  isr <<  signalList[i] << "_2017"  << " \t " << regionSaveName << " \t " <<    sign* isr_uncert  <<std::endl;

	  int massSbottom = 0;
	  int massChi = 0;

	  TString sigSampleName(signalList[i]);

	  std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	  std::string::size_type sz;   // alias of size_t

	  if( !sigSampleName.Contains("TChi")  )
	    massChi = std::stoi (mLSP,&sz);
	  else if( sigSampleName.Contains("TChiH")  )
	    massChi = 1;
	  else if( sigSampleName.Contains("TChiWH")  ){ 
	    mLSP = signalList[i].substr( signalList[i].find("_HToGG_")+13);
	    std::cout << "m2 " << mLSP << std::endl;
	    massChi = std::stoi (mLSP,&sz);
	  }

	  std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	  if( sigSampleName.Contains("TChiH")  ){
	    mSbottom = signalList[i].substr( signalList[i].find("_m")+2);
	  }
	  if( sigSampleName.Contains("TChiWH")  ){
	    mSbottom = signalList[i].substr( signalList[i].find("_HToGG_")+8);
	  }

	  massSbottom = std::stoi (mSbottom,&sz);
	
	  int binBottom = h2D_isr_2017->GetXaxis()->FindBin( (float)massSbottom );
	  int binChi    = h2D_isr_2017->GetYaxis()->FindBin( (float)massChi );

	  h2D_isr_2017->SetBinContent( binBottom, binChi, h_isr_2017->GetMean() );
	  h2D_isr_2017_err->SetBinContent( binBottom, binChi, sign* isr_uncert );

	  delete h_isr_2017; delete h_isr_2017_UP; delete h_isr_2017_DN;

	  // if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	  //   sigAmount <<  signalList[i] << "_2017" << " \t " << regionSaveName << " \t " <<  sig_ds.sumEntries()  << std::endl;
	


	  ///////////////////////////// GAMMA SF_2017///////////////////////////////////////////////////
	  TH1D* h_gammasf_2017    = new TH1D("h_gammasf_2017",    "", 100, 0, 2);
	  TH1D* h_gammasf_2017_UP = new TH1D("h_gammasf_2017_UP", "", 100, 0, 2);
	  TH1D* h_gammasf_2017_DN = new TH1D("h_gammasf_2017_DN", "", 100, 0, 2);
 
	  tree_sig_2017->Draw("weight_gammasf    >> h_gammasf_2017");
	  tree_sig_2017->Draw("weight_gammasf_UP >> h_gammasf_2017_UP");
	  tree_sig_2017->Draw("weight_gammasf_DN >> h_gammasf_2017_DN");

	  float sign_gammasf_2017 = 1;
	  float gammasf_2017_uncert = max (fabs(h_gammasf_2017_UP->GetMean() - h_gammasf_2017->GetMean()), fabs( h_gammasf_2017->GetMean() - h_gammasf_2017_DN->GetMean() ) );

	  if( fabs(h_gammasf_2017_UP->GetMean() - h_gammasf_2017->GetMean()) >  fabs(  h_gammasf_2017->GetMean() - h_gammasf_2017_DN->GetMean() ))
	    sign_gammasf_2017 = ((h_gammasf_2017_UP->GetMean() - h_gammasf_2017->GetMean())  > 0 ) - ((h_gammasf_2017_UP->GetMean() - h_gammasf_2017->GetMean())  < 0) ;
	  else
	    sign_gammasf_2017 = ((h_gammasf_2017->GetMean() - h_gammasf_2017_DN->GetMean())  > 0 ) - ((h_gammasf_2017->GetMean() - h_gammasf_2017_DN->GetMean())  < 0) ;
	  
	  //	gammasf_2017 <<  signalList[i] << " \t " << regionSaveName << " \t " <<    sign_gammasf_2017* gammasf_2017_uncert  <<std::endl;

	  h2D_gammasf_2017->SetBinContent( binBottom, binChi, h_gammasf_2017->GetMean() );
	  h2D_gammasf_2017_err->SetBinContent( binBottom, binChi, sign_gammasf_2017* gammasf_2017_uncert );

	  delete h_gammasf_2017; delete h_gammasf_2017_UP; delete h_gammasf_2017_DN;



	  ///////////////////////////// LEPTON SF_2017///////////////////////////////////////////////////
	  TH1D* h_lepsf_2017    = new TH1D("h_lepsf_2017",    "", 100, 0, 2);
	  TH1D* h_lepsf_2017_UP = new TH1D("h_lepsf_2017_UP", "", 100, 0, 2);
	  TH1D* h_lepsf_2017_DN = new TH1D("h_lepsf_2017_DN", "", 100, 0, 2);
 
	  tree_sig_2017->Draw("weight_lepsf    >> h_lepsf_2017");
	  tree_sig_2017->Draw("weight_lepsf_UP >> h_lepsf_2017_UP");
	  tree_sig_2017->Draw("weight_lepsf_DN >> h_lepsf_2017_DN");

	  float sign_lepsf_2017 = 1;
	  float lepsf_2017_uncert = max (fabs(h_lepsf_2017_UP->GetMean() - h_lepsf_2017->GetMean()), fabs( h_lepsf_2017->GetMean() - h_lepsf_2017_DN->GetMean() ) );

	  if( fabs(h_lepsf_2017_UP->GetMean() - h_lepsf_2017->GetMean()) >  fabs(  h_lepsf_2017->GetMean() - h_lepsf_2017_DN->GetMean() ))
	    sign_lepsf_2017 = ((h_lepsf_2017_UP->GetMean() - h_lepsf_2017->GetMean())  > 0 ) - ((h_lepsf_2017_UP->GetMean() - h_lepsf_2017->GetMean())  < 0) ;
	  else
	    sign_lepsf_2017 = ((h_lepsf_2017->GetMean() - h_lepsf_2017_DN->GetMean())  > 0 ) - ((h_lepsf_2017->GetMean() - h_lepsf_2017_DN->GetMean())  < 0) ;
	  
	  //	lepsf_2017 <<  signalList[i] << " \t " << regionSaveName << " \t " <<    sign_lepsf_2017* lepsf_2017_uncert  <<std::endl;

	  h2D_lepsf_2017->SetBinContent( binBottom, binChi, h_lepsf_2017->GetMean() );
	  h2D_lepsf_2017_err->SetBinContent( binBottom, binChi, sign_lepsf_2017* lepsf_2017_uncert );

	  delete h_lepsf_2017; delete h_lepsf_2017_UP; delete h_lepsf_2017_DN;
	

	  ///////////////////////////// BTAG SF_2017///////////////////////////////////////////////////
	  TH1D* h_btagsf_heavy_2017    = new TH1D("h_btagsf_heavy_2017",    "", 100, 0, 2);
	  TH1D* h_btagsf_heavy_2017_UP = new TH1D("h_btagsf_heavy_2017_UP", "", 100, 0, 2);
	  TH1D* h_btagsf_heavy_2017_DN = new TH1D("h_btagsf_heavy_2017_DN", "", 100, 0, 2);
 
	  tree_sig_2017->Draw("weight_btagsf    >> h_btagsf_heavy_2017");
	  tree_sig_2017->Draw("weight_btagsf_heavy_UP >> h_btagsf_heavy_2017_UP");
	  tree_sig_2017->Draw("weight_btagsf_heavy_DN >> h_btagsf_heavy_2017_DN");

	  float sign_btagsf_heavy_2017 = 1;
	  float btagsf_heavy_2017_uncert = max (fabs(h_btagsf_heavy_2017_UP->GetMean() - h_btagsf_heavy_2017->GetMean()), fabs( h_btagsf_heavy_2017->GetMean() - h_btagsf_heavy_2017_DN->GetMean() ) );

	  if( fabs(h_btagsf_heavy_2017_UP->GetMean() - h_btagsf_heavy_2017->GetMean()) >  fabs(  h_btagsf_heavy_2017->GetMean() - h_btagsf_heavy_2017_DN->GetMean() ))
	    sign_btagsf_heavy_2017 = ((h_btagsf_heavy_2017_UP->GetMean() - h_btagsf_heavy_2017->GetMean())  > 0 ) - ((h_btagsf_heavy_2017_UP->GetMean() - h_btagsf_heavy_2017->GetMean())  < 0) ;
	  else
	    sign_btagsf_heavy_2017 = ((h_btagsf_heavy_2017->GetMean() - h_btagsf_heavy_2017_DN->GetMean())  > 0 ) - ((h_btagsf_heavy_2017->GetMean() - h_btagsf_heavy_2017_DN->GetMean())  < 0) ;
	  
	  //	btagsf_heavy_2017 <<  signalList[i] << " \t " << regionSaveName << " \t " <<    sign_btagsf_heavy_2017* btagsf_heavy_2017_uncert  <<std::endl;

	  h2D_btagsf_heavy_2017->SetBinContent( binBottom, binChi, h_btagsf_heavy_2017->GetMean() );
	  h2D_btagsf_heavy_2017_err->SetBinContent( binBottom, binChi, sign_btagsf_heavy_2017* btagsf_heavy_2017_uncert );

	  delete h_btagsf_heavy_2017; delete h_btagsf_heavy_2017_UP; delete h_btagsf_heavy_2017_DN;


	  ///////////////////////////// BTAG SF_2017///////////////////////////////////////////////////
	  TH1D* h_btagsf_light_2017    = new TH1D("h_btagsf_light_2017",    "", 100, 0, 2);
	  TH1D* h_btagsf_light_2017_UP = new TH1D("h_btagsf_light_2017_UP", "", 100, 0, 2);
	  TH1D* h_btagsf_light_2017_DN = new TH1D("h_btagsf_light_2017_DN", "", 100, 0, 2);
 
	  tree_sig_2017->Draw("weight_btagsf    >> h_btagsf_light_2017");
	  tree_sig_2017->Draw("weight_btagsf_light_UP >> h_btagsf_light_2017_UP");
	  tree_sig_2017->Draw("weight_btagsf_light_DN >> h_btagsf_light_2017_DN");

	  float sign_btagsf_light_2017 = 1;
	  float btagsf_light_2017_uncert = max (fabs(h_btagsf_light_2017_UP->GetMean() - h_btagsf_light_2017->GetMean()), fabs( h_btagsf_light_2017->GetMean() - h_btagsf_light_2017_DN->GetMean() ) );

	  if( fabs(h_btagsf_light_2017_UP->GetMean() - h_btagsf_light_2017->GetMean()) >  fabs(  h_btagsf_light_2017->GetMean() - h_btagsf_light_2017_DN->GetMean() ))
	    sign_btagsf_light_2017 = ((h_btagsf_light_2017_UP->GetMean() - h_btagsf_light_2017->GetMean())  > 0 ) - ((h_btagsf_light_2017_UP->GetMean() - h_btagsf_light_2017->GetMean())  < 0) ;
	  else
	    sign_btagsf_light_2017 = ((h_btagsf_light_2017->GetMean() - h_btagsf_light_2017_DN->GetMean())  > 0 ) - ((h_btagsf_light_2017->GetMean() - h_btagsf_light_2017_DN->GetMean())  < 0) ;
	  
	  //	btagsf_light_2017 <<  signalList[i] << " \t " << regionSaveName << " \t " <<    sign_btagsf_light_2017* btagsf_light_2017_uncert  <<std::endl;

	  h2D_btagsf_light_2017->SetBinContent( binBottom, binChi, h_btagsf_light_2017->GetMean() );
	  h2D_btagsf_light_2017_err->SetBinContent( binBottom, binChi, sign_btagsf_light_2017* btagsf_light_2017_uncert );

	  delete h_btagsf_light_2017; delete h_btagsf_light_2017_UP; delete h_btagsf_light_2017_DN;




	  /////////////////scale variation SIGNAL
	  vector< float > scaleUncert;
	  for( int i=0; i<110; i++){
	    TH1D* h_scalesTemp_2017    = new TH1D("h_scalesTemp_2017",    "", 100, 0, 2);
	    tree_sig_2017->Draw(Form(" weight_scales_%d    >> h_scalesTemp_2017", i));
	    float newValue = h_scalesTemp_2017->GetMean();
	    if( (newValue > 1.4 ) || (newValue < 0.6 ))
	      newValue = 1.;  // savety first, in case of crazy big diff, it will thus take a diff variation
	    scaleUncert.push_back( newValue );
	    //	    scaleUncert.push_back( h_scalesTemp_2017->GetMean() );
	    delete h_scalesTemp_2017;
	  }
	  double scaleMax = *max_element(scaleUncert.begin(), scaleUncert.end()-100);
	  double scaleMin = *min_element(scaleUncert.begin(), scaleUncert.end()-100);
	  double pdfMax   = *max_element(scaleUncert.begin()+10, scaleUncert.end());
	  double pdfMin   = *min_element(scaleUncert.begin()+10, scaleUncert.end());
	  float scaleVar_2017_uncert = ((scaleMax -1.)>(1.-scaleMin)) ?  scaleMax : scaleMin;
	  float pdfVar_2017_uncert = ((pdfMax -1.)>(1.-pdfMin)) ?  pdfMax : pdfMin;

	  h2D_scaleVar_2017->SetBinContent( binBottom, binChi, scaleVar_2017_uncert  );
	  h2D_pdfVar_2017  ->SetBinContent( binBottom, binChi, pdfVar_2017_uncert  );




	  RooDataSet sig_ds_125_2017_JECup(  sig_ds_2017_JECup, Form("sig_2017_JECup_125_13TeV_%s", regionSaveName.c_str() ) );
	  RooDataSet sig_ds_125_2017_JECdn(  sig_ds_2017_JECdn, Form("sig_2017_JECdn_125_13TeV_%s", regionSaveName.c_str() ) );
	  ((RooRealVar*)sig_ds_125_2017_JECup.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
	  ((RooRealVar*)sig_ds_125_2017_JECdn.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
	  float var2017_JECup = float(( sig_ds_125_2017_JECup.sumEntries() - sig_ds_2017.sumEntries()*2. )) /(float(sig_ds_2017.sumEntries())*2. ) ; 
	  float var2017_JECdn = float((  sig_ds_2017.sumEntries()*2.  - sig_ds_2017_JECdn.sumEntries() )) /(float(sig_ds_2017.sumEntries())*2. ) ; 

	  if( float(sig_ds_2017.sumEntries() == 0 )){
	    var2017_JECup = 0;
	    var2017_JECdn = 0;
	  }
	  // std::cout<< "central = " <<  sig_ds_2017_125.sumEntries()  << std::endl;
	  // std::cout<< "up      = " <<  sig_ds_125_2017_JECup.sumEntries()  << std::endl;
	  // std::cout<< "dn      = " <<  sig_ds_125_2017_JECdn.sumEntries()  << std::endl;

	  h2D_JECup_2017->SetBinContent( binBottom, binChi, var2017_JECup  );
	  h2D_JECdn_2017->SetBinContent( binBottom, binChi, var2017_JECdn  );
	  std::cout << "var up 2017 " <<  var2017_JECup << std::endl;
	  std::cout << "var dn 2017 " <<  var2017_JECdn << std::endl;

	  RooDataSet sig_ds_125_2016_JECup(  sig_ds_2016_JECup, Form("sig_2016_JECup_125_13TeV_%s", regionSaveName.c_str() ) );
	  RooDataSet sig_ds_125_2016_JECdn(  sig_ds_2016_JECdn, Form("sig_2016_JECdn_125_13TeV_%s", regionSaveName.c_str() ) );
	  ((RooRealVar*)sig_ds_125_2016_JECup.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
	  ((RooRealVar*)sig_ds_125_2016_JECdn.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
	  float var2016_JECup = float(( sig_ds_125_2016_JECup.sumEntries() - sig_ds.sumEntries()*2. ))  /(float(sig_ds.sumEntries())*2. ) ; 
	  float var2016_JECdn = float((sig_ds.sumEntries()*2.  -  sig_ds_125_2016_JECdn.sumEntries() )) /(float(sig_ds.sumEntries())*2. ) ;

	  if( float(sig_ds.sumEntries() == 0) ){
	    var2016_JECup = 0;
	    var2016_JECdn = 0;
	  }

	  // std::cout<< "central = " <<  sig_ds_125.sumEntries()  << std::endl;
	  // std::cout<< "up      = " <<  sig_ds_125_2016_JECup.sumEntries()  << std::endl;
	  // std::cout<< "dn      = " <<  sig_ds_125_2016_JECdn.sumEntries()  << std::endl;

	  h2D_JECup->SetBinContent( binBottom, binChi, var2016_JECup  );
	  h2D_JECdn->SetBinContent( binBottom, binChi, var2016_JECdn  );
	  std::cout << "var up 2016 " <<  var2016_JECup << std::endl;
	  std::cout << "var dn 2016 " <<  var2016_JECdn << std::endl;



	  ///////////////////////////// GEN MET correction ///////////////////////////////////////
	  /////////// factor of 1/2 already in the weight, so don't divide by it again! //////////
	  float genMET_uncert_2017 =  sig_ds_2017.sumEntries() - sig_genMET_ds_2017.sumEntries() ;

	  if( ((sig_ds_2017.sumEntries() + sig_genMET_ds_2017.sumEntries() ) == 0) || (i >=( regionSelection.size()-4) ) )
	    genMET_uncert_2017 = 0;
	  else
	    genMET_uncert_2017	 /= ( sig_ds_2017.sumEntries() + sig_genMET_ds_2017.sumEntries() );
	   
	  h2D_genMET_2017->SetBinContent( binBottom, binChi, genMET_uncert_2017 );

	  std::cout << "genMET 2017 syst                         " <<  genMET_uncert_2017 << std::endl;

	  std::cout << "genMET 2017 yield " << sig_ds_2017.sumEntries()  << std::endl;
	  std::cout << "genMET 2017 gen yield " << sig_genMET_ds_2017.sumEntries()  << std::endl;
	  std::cout << "sig yield 2017 summed " << sig_ds_2017.sumEntries() + sig_genMET_ds_2017.sumEntries()  << std::endl;

	  // std::cout << "at end 2017 syst "  << std::endl;

	}//end of 2017 syst


	sig_ds_2017.append( sig_genMET_ds_2017 );

	RooDataSet sig_ds_2017_130(  sig_ds_2017, Form("%s_2017_130_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ) );
	RooDataSet sig_ds_2017_125(  sig_ds_2017, Form("%s_2017_125_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ) ); 
	RooDataSet sig_ds_2017_120(  sig_ds_2017, Form("%s_2017_120_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ) );

	((RooRealVar*)sig_ds_2017_130.addColumn( hgg_mass_130 ))->setRange(100.,180.);  
	((RooRealVar*)sig_ds_2017_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
	((RooRealVar*)sig_ds_2017_120.addColumn( hgg_mass_120 ))->setRange(100.,180.); 

	ws_sig.import( sig_ds_2017_125  );
	ws_sig.import( sig_ds_2017_130  );
	ws_sig.import( sig_ds_2017_120  );

	if( i==0 )
	  sigAmount <<  signalList[i] << "_2017 \t " << regionSaveName << " \t " <<  sig_ds_2017.sumEntries()  << std::endl;

	//	  std::cout << "done 2017 syst "  << std::endl;

      }

      sig_ds.append( sig_genMET_ds );


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

      //	  std::cout << "done 2017 "  << std::endl;

    }//done signals for 1 region
   

    //    std::cout << "done a region "  << std::endl;


    // h2D_isr->Write();
    // h2D_isr_err->Write();

    // delete h2D_isr; delete h2D_isr_err;

    // h2D_isr_2017->Write();
    // h2D_isr_2017_err->Write();

    // delete h2D_isr_2017; delete h2D_isr_2017_err;


    // h2D_gammasf->Write();
    // h2D_gammasf_err->Write();

    // delete h2D_gammasf; delete h2D_gammasf_err;

    // h2D_gammasf_2017->Write();
    // h2D_gammasf_2017_err->Write();

    // delete h2D_gammasf_2017; delete h2D_gammasf_2017_err;

    // h2D_lepsf->Write();
    // h2D_lepsf_err->Write();

    // delete h2D_lepsf; delete h2D_lepsf_err;

    // h2D_lepsf_2017->Write();
    // h2D_lepsf_2017_err->Write();

    // delete h2D_lepsf_2017; delete h2D_lepsf_2017_err;

    // h2D_btagsf_heavy->Write();
    // h2D_btagsf_heavy_err->Write();

    // delete h2D_btagsf_heavy; delete h2D_btagsf_heavy_err;

    // h2D_btagsf_heavy_2017->Write();
    // h2D_btagsf_heavy_2017_err->Write();

    // delete h2D_btagsf_heavy_2017; delete h2D_btagsf_heavy_2017_err;


    // h2D_btagsf_light->Write();
    // h2D_btagsf_light_err->Write();

    // delete h2D_btagsf_light; delete h2D_btagsf_light_err;

    // h2D_btagsf_light_2017->Write();
    // h2D_btagsf_light_2017_err->Write();

    // delete h2D_btagsf_light_2017; delete h2D_btagsf_light_2017_err;


    // h2D_scaleVar->Write();
    // delete h2D_scaleVar;
    // h2D_scaleVar_2017->Write();
    // delete h2D_scaleVar_2017;

    // h2D_pdfVar->Write();
    // delete h2D_pdfVar;
    // h2D_pdfVar_2017->Write();
    // delete h2D_pdfVar_2017;


    // h2D_JECup->Write();
    // delete h2D_JECup;
    // h2D_JECup_2017->Write();
    // delete h2D_JECup_2017;

    // h2D_JECdn->Write();
    // delete h2D_JECdn;
    // h2D_JECdn_2017->Write();
    // delete h2D_JECdn_2017;


    //  } // loop over mt2 bin(s)
      
    }


    h2D_genMET->Write();
    delete h2D_genMET;

    h2D_genMET_2017->Write();
    delete h2D_genMET_2017;


    h2D_isr->Write();
    h2D_isr_err->Write();

    delete h2D_isr; delete h2D_isr_err;

    h2D_isr_2017->Write();
    h2D_isr_2017_err->Write();

    delete h2D_isr_2017; delete h2D_isr_2017_err;


    h2D_gammasf->Write();
    h2D_gammasf_err->Write();

    delete h2D_gammasf; delete h2D_gammasf_err;

    h2D_gammasf_2017->Write();
    h2D_gammasf_2017_err->Write();

    delete h2D_gammasf_2017; delete h2D_gammasf_2017_err;

    h2D_lepsf->Write();
    h2D_lepsf_err->Write();

    delete h2D_lepsf; delete h2D_lepsf_err;

    h2D_lepsf_2017->Write();
    h2D_lepsf_2017_err->Write();

    delete h2D_lepsf_2017; delete h2D_lepsf_2017_err;

    h2D_btagsf_heavy->Write();
    h2D_btagsf_heavy_err->Write();

    delete h2D_btagsf_heavy; delete h2D_btagsf_heavy_err;

    h2D_btagsf_heavy_2017->Write();
    h2D_btagsf_heavy_2017_err->Write();

    delete h2D_btagsf_heavy_2017; delete h2D_btagsf_heavy_2017_err;


    h2D_btagsf_light->Write();
    h2D_btagsf_light_err->Write();

    delete h2D_btagsf_light; delete h2D_btagsf_light_err;

    h2D_btagsf_light_2017->Write();
    h2D_btagsf_light_2017_err->Write();

    delete h2D_btagsf_light_2017; delete h2D_btagsf_light_2017_err;


    h2D_scaleVar->Write();
    delete h2D_scaleVar;
    h2D_scaleVar_2017->Write();
    delete h2D_scaleVar_2017;

    h2D_pdfVar->Write();
    delete h2D_pdfVar;
    h2D_pdfVar_2017->Write();
    delete h2D_pdfVar_2017;


    h2D_JECup->Write();
    delete h2D_JECup;
    h2D_JECup_2017->Write();
    delete h2D_JECup_2017;

    h2D_JECdn->Write();
    delete h2D_JECdn;
    h2D_JECdn_2017->Write();
    delete h2D_JECdn_2017;
    
  }//done loop over signals


  dataAmount.close();
  //bgSumAmount.close();
  higgsAmount.close();

  higgsGluGlu.close();
  higgsVBF.close();
  higgsVH.close();
  higgsttH.close();
  higgsbbH.close();
  //  higgsTHX.close();
  //  higgsdiH.close();


  lepsf_higgs.close();
  gammasf_higgs.close();
  btagsf_heavy_higgs.close();
  btagsf_light_higgs.close();
  lepsf_higgs_2017.close();
  gammasf_higgs_2017.close();
  btagsf_heavy_higgs_2017.close();
  btagsf_light_higgs_2017.close();
  scaleVar_higgs.close();
  scaleVar_higgs_2017.close();
  pdfVar_higgs.close();
  pdfVar_higgs_2017.close();


  higgsSyst_2016_JECup.close();
  higgsSyst_2016_JECdn.close();

  higgsSyst_2017_JECup.close();
  higgsSyst_2017_JECdn.close();

  isr.close();
 
  ws_data.Write();
  ws_sig.Write();
  //bg  ws_bg.Write();

  dataFile->Close();

  return 0;
}










  // old selection  
  // regionSelection_ph.push_back(" ( is1El) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  // regionSelection_ph.push_back(" ( is1Mu) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(isDiLepZ) )"  );
  // regionSelection_ph.push_back(" ( is1El) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  // regionSelection_ph.push_back(" ( is1Mu) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(isDiLepZ) )"  );
  // regionSelection_ph.push_back(" ( isDiBH ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBZ) && !(isDiLepZ) )" );
  // regionSelection_ph.push_back(" ( isDiBZ ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBH) && !(isDiLepZ) )" );
  // regionSelection_ph.push_back(" ( isDiBH ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBZ) && !(isDiLepZ) )" );
  // regionSelection_ph.push_back(" ( isDiBZ ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBH) && !(isDiLepZ) )" );
  // old selection  
  // regionSelection_ph.push_back(" ( is1El) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  // regionSelection_ph.push_back(" ( is1Mu) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(is1El) && !(isDiLepZ) )"  );
  // regionSelection_ph.push_back(" ( is1El) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  // regionSelection_ph.push_back(" ( is1Mu) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(is1El) && !(isDiLepZ) )"  );
  // regionSelection_ph.push_back(" ( isDiLepZ ) && ( (!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) )"     );
  // regionSelection_ph.push_back(" ( isDiBH ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  // regionSelection_ph.push_back(" ( isDiBZ ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBH) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  // regionSelection_ph.push_back(" ( isDiBH ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  // regionSelection_ph.push_back(" ( isDiBZ ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBH) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );

  // regionNames_ph.push_back("j0_b0toInf_pT2");  // regionNames_ph.push_back("j0_b0toInf_pT3");  // regionNames_ph.push_back("j0_b0toInf_pT4");  // regionNames_ph.push_back("j0_b0toInf_pT5");  // regionNames_ph.push_back("j0_b0toInf_pT6");  // regionNames_ph.push_back("j0_b0toInf_pT7");  // regionNames_ph.push_back("j0_b0toInf_pT8");   // regionNames_ph.push_back("j0_b0toInf_pT9");





    // // // Print unbinned dataset with default frame binning (100 bins)
    // //RooPlot* frame3 = CMS_hgg_mass.frame(Title("H mass")) ;
    // RooPlot* frame3 = h_mass_var.frame(Title("H mass")) ;
  
    // data_ds.plotOn(frame3, MarkerColor(kRed) ) ;
    // higgs_ds_125.plotOn(frame3, MarkerColor(kGreen) ) ;
    // sig_ds_125.plotOn(frame3) ;
    //bg    // bg_ds.plotOn(frame3, MarkerColor(kBlue) ) ;

    // TCanvas* c = new TCanvas("rf102_dataimport","rf102_dataimport",800,800) ;
    // frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;

    // c->SaveAs(  Form("%s/CMS_hgg_mass_%s.png", outputdir.c_str(), regionSaveName.c_str() ) );



  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && ((h_pt/h_mass) <  0.1)) " );
  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.1<=(h_pt/h_mass) && (h_pt/h_mass) < 0.2)) " );
  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.2<=(h_pt/h_mass) && (h_pt/h_mass) < 0.3)) " );
  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.3<=(h_pt/h_mass) && (h_pt/h_mass) < 0.4)) " );
  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.4<=(h_pt/h_mass) && (h_pt/h_mass) < 0.6)) " );
  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.6<=(h_pt/h_mass) && (h_pt/h_mass) < 0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.8<=(h_pt/h_mass) && (h_pt/h_mass) < 1.0)) " );
  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (1.0<=(h_pt/h_mass) && (h_pt/h_mass) < 1.2)) " );
  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (1.2<=(h_pt/h_mass) && (h_pt/h_mass) < 1.4)) " );
  // regionSelection_ph.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && ((h_pt/h_mass) >= 1.4)) " );
  // int n0j = 10;
  // regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 0 && ((h_pt/h_mass) <  0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 0 && ((h_pt/h_mass) >= 0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 1 && ((h_pt/h_mass) <  0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 1 && ((h_pt/h_mass) >= 0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets >= 2 && ((h_pt/h_mass) <  0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets >= 2 && ((h_pt/h_mass) >= 0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 0 && ((h_pt/h_mass) <  0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 0 && ((h_pt/h_mass) >= 0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 1 && ((h_pt/h_mass) <  0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets == 1 && ((h_pt/h_mass) >= 0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets >= 2 && ((h_pt/h_mass) <  0.8)) " );
  // regionSelection_ph.push_back(" (gg_nJets>=4 && gg_nBJets >= 2 && ((h_pt/h_mass) >= 0.8)) " );

















  // for( unsigned int j=0; j<mt2_bins.size(); j++){ 
  //   for( unsigned int i=0; i < regionLength; i++){
  //     //if( i>nMT2split + n0j   ) continue;      // if( j==0 )      //	if ( i<n0j ) continue;      //std::cout 
  //     //std::cout << mt2_bins[j] << " with sel " << mt2_sel[j] << std::endl;      //	std::cout << regionSelection_ph[i] << std::endl;
  //     regionSelection.push_back( regionSelection_ph[i] );
  //     // std::cout << regionSelection[ i + j * regionLength ] << std::endl;      // std::cout << i + j * regionLength << std::endl;
  //     regionSelection[ i + j*regionLength ] +=  "  && " +  mt2_sel[j] ;
  //     regionNames.push_back( regionNames_ph[i] );
  //     regionNames[i+j*regionLength ] = ( regionNames_ph[i] + "_mt2_" + mt2_bins[j]);
  //     //	std::cout << "selection 0 " << regionSelection[i] << std::endl;
  //     // if( j==0  ){      //   if ( i<n0j ) continue;
  //     //   regionSelection.push_back( regionSelection[i] );      //   //	  regionSelection[i] +=  "  && ( hgg_mt2<30   )";
  //     //   regionSelection[i] +=  "  && ( hgg_mt2<30   )";      //   regionNames.push_back( regionNames[i] );
  //     //   regionNames[i] +=  "_mt2_0";      //   //	std::cout << "selection 0 " << regionSelection[i] << std::endl;
  //     // }else{      //   if ( i>(nMT2split )  ) continue;      //   //  	if ( i>(nMT2split - n0j) ) continue;
  //     //   regionSelection[ regionLength + i ] +=  "  && ( hgg_mt2>=30   )";
  //     //   regionNames[ regionLength + i ] +=  "_mt2_30";
  //     //   //	std::cout << "selection 1 " << regionSelection[regionLength + i ] << std::endl;       // }
  //   }
  // }





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
