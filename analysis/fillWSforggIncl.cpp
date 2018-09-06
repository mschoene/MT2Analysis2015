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
  //  if( doT2bH )
  doISRsyst = true; //doing it now for all, cause I can

  std::string ggDir_2017;
  bool doComb = false;
  if( argc>3 ) {
    std::cout << "Doing 2016 and 2017 combination" << std::endl;
    doComb = true;

    std::string configFileName_2017(argv[3]);
    MT2Config cfg_2017( configFileName_2017);
    ggDir_2017 = cfg_2017.getEventYieldDir() + "/diPhotonControlRegion/";
  }
 



  //  vector<std::string> signalList_T2bH = { "SMS_T2bH_mSbottom250_mLSP100"};

  //w/o 500 1  vector<std::string> signalList_T2bH = { "SMS_T2bH_mSbottom250_mLSP100", "SMS_T2bH_mSbottom250_mLSP1", "SMS_T2bH_mSbottom250_mLSP50", "SMS_T2bH_mSbottom300_mLSP100", "SMS_T2bH_mSbottom300_mLSP150", "SMS_T2bH_mSbottom300_mLSP1", "SMS_T2bH_mSbottom300_mLSP50", "SMS_T2bH_mSbottom350_mLSP100", "SMS_T2bH_mSbottom350_mLSP150", "SMS_T2bH_mSbottom350_mLSP1", "SMS_T2bH_mSbottom350_mLSP200", "SMS_T2bH_mSbottom350_mLSP50", "SMS_T2bH_mSbottom400_mLSP100", "SMS_T2bH_mSbottom400_mLSP150", "SMS_T2bH_mSbottom400_mLSP1", "SMS_T2bH_mSbottom400_mLSP200", "SMS_T2bH_mSbottom400_mLSP250", "SMS_T2bH_mSbottom400_mLSP50", "SMS_T2bH_mSbottom450_mLSP1", "SMS_T2bH_mSbottom450_mLSP100", "SMS_T2bH_mSbottom450_mLSP150", "SMS_T2bH_mSbottom450_mLSP200", "SMS_T2bH_mSbottom450_mLSP250", "SMS_T2bH_mSbottom450_mLSP300", "SMS_T2bH_mSbottom450_mLSP50", "SMS_T2bH_mSbottom500_mLSP100", "SMS_T2bH_mSbottom500_mLSP150", "SMS_T2bH_mSbottom500_mLSP200", "SMS_T2bH_mSbottom500_mLSP250", "SMS_T2bH_mSbottom500_mLSP300", "SMS_T2bH_mSbottom500_mLSP50", "SMS_T2bH_mSbottom600_mLSP1", "SMS_T2bH_mSbottom600_mLSP100", "SMS_T2bH_mSbottom600_mLSP200", "SMS_T2bH_mSbottom600_mLSP300" };

  vector<std::string> signalList_T2bH = { "SMS_T2bH_mSbottom250_mLSP100", "SMS_T2bH_mSbottom250_mLSP1", "SMS_T2bH_mSbottom250_mLSP50", "SMS_T2bH_mSbottom300_mLSP100", "SMS_T2bH_mSbottom300_mLSP150", "SMS_T2bH_mSbottom300_mLSP1", "SMS_T2bH_mSbottom300_mLSP50", "SMS_T2bH_mSbottom350_mLSP100", "SMS_T2bH_mSbottom350_mLSP150", "SMS_T2bH_mSbottom350_mLSP1", "SMS_T2bH_mSbottom350_mLSP200", "SMS_T2bH_mSbottom350_mLSP50", "SMS_T2bH_mSbottom400_mLSP100", "SMS_T2bH_mSbottom400_mLSP150", "SMS_T2bH_mSbottom400_mLSP1", "SMS_T2bH_mSbottom400_mLSP200", "SMS_T2bH_mSbottom400_mLSP250", "SMS_T2bH_mSbottom400_mLSP50", "SMS_T2bH_mSbottom450_mLSP1", "SMS_T2bH_mSbottom450_mLSP100", "SMS_T2bH_mSbottom450_mLSP150", "SMS_T2bH_mSbottom450_mLSP200", "SMS_T2bH_mSbottom450_mLSP250", "SMS_T2bH_mSbottom450_mLSP300", "SMS_T2bH_mSbottom450_mLSP50", "SMS_T2bH_mSbottom500_mLSP100", "SMS_T2bH_mSbottom500_mLSP150", "SMS_T2bH_mSbottom500_mLSP1", "SMS_T2bH_mSbottom500_mLSP200", "SMS_T2bH_mSbottom500_mLSP250", "SMS_T2bH_mSbottom500_mLSP300", "SMS_T2bH_mSbottom500_mLSP50", "SMS_T2bH_mSbottom600_mLSP1", "SMS_T2bH_mSbottom600_mLSP100", "SMS_T2bH_mSbottom600_mLSP200", "SMS_T2bH_mSbottom600_mLSP300" };

  //  vector<std::string> signalList_WH = { "SMS_TChiWH_HToGG_127_1","SMS_TChiWH_HToGG_150_1","SMS_TChiWH_HToGG_150_24"};
  //  vector<std::string> signalList_WH = { "SMS_TChiWH_HToGG_127_1","SMS_TChiWH_HToGG_150_1","SMS_TChiWH_HToGG_150_24","SMS_TChiWH_HToGG_175_1","SMS_TChiWH_HToGG_175_25","SMS_TChiWH_HToGG_175_49","SMS_TChiWH_HToGG_200_1","SMS_TChiWH_HToGG_200_25","SMS_TChiWH_HToGG_200_50","SMS_TChiWH_HToGG_200_74" };


  vector<std::string> signalList_WH = {"SMS_TChiWH_HToGG_127_1","SMS_TChiWH_HToGG_150_1","SMS_TChiWH_HToGG_150_24","SMS_TChiWH_HToGG_175_1","SMS_TChiWH_HToGG_175_25","SMS_TChiWH_HToGG_175_49","SMS_TChiWH_HToGG_200_1","SMS_TChiWH_HToGG_200_25","SMS_TChiWH_HToGG_200_50","SMS_TChiWH_HToGG_200_74","SMS_TChiWH_HToGG_225_1","SMS_TChiWH_HToGG_225_25","SMS_TChiWH_HToGG_225_50","SMS_TChiWH_HToGG_225_75","SMS_TChiWH_HToGG_225_99","SMS_TChiWH_HToGG_250_1","SMS_TChiWH_HToGG_250_25","SMS_TChiWH_HToGG_250_50","SMS_TChiWH_HToGG_250_75","SMS_TChiWH_HToGG_250_100","SMS_TChiWH_HToGG_250_124"};//,"SMS_TChiWH_HToGG_275_1","SMS_TChiWH_HToGG_275_25","SMS_TChiWH_HToGG_275_50","SMS_TChiWH_HToGG_275_75","SMS_TChiWH_HToGG_275_100","SMS_TChiWH_HToGG_275_125","SMS_TChiWH_HToGG_275_149","SMS_TChiWH_HToGG_300_1","SMS_TChiWH_HToGG_300_25","SMS_TChiWH_HToGG_300_50","SMS_TChiWH_HToGG_300_75","SMS_TChiWH_HToGG_300_100","SMS_TChiWH_HToGG_300_125","SMS_TChiWH_HToGG_300_150","SMS_TChiWH_HToGG_300_174","SMS_TChiWH_HToGG_325_1","SMS_TChiWH_HToGG_325_25","SMS_TChiWH_HToGG_325_50","SMS_TChiWH_HToGG_325_75","SMS_TChiWH_HToGG_325_100","SMS_TChiWH_HToGG_325_125","SMS_TChiWH_HToGG_325_150","SMS_TChiWH_HToGG_325_175","SMS_TChiWH_HToGG_325_199"};//,"SMS_TChiWH_HToGG_350_1","SMS_TChiWH_HToGG_350_25","SMS_TChiWH_HToGG_350_50","SMS_TChiWH_HToGG_350_75","SMS_TChiWH_HToGG_350_100","SMS_TChiWH_HToGG_350_125","SMS_TChiWH_HToGG_350_150","SMS_TChiWH_HToGG_350_175","SMS_TChiWH_HToGG_350_200","SMS_TChiWH_HToGG_350_224","SMS_TChiWH_HToGG_375_1","SMS_TChiWH_HToGG_375_25","SMS_TChiWH_HToGG_375_50","SMS_TChiWH_HToGG_375_75","SMS_TChiWH_HToGG_375_100","SMS_TChiWH_HToGG_375_125","SMS_TChiWH_HToGG_375_150","SMS_TChiWH_HToGG_375_175","SMS_TChiWH_HToGG_375_200","SMS_TChiWH_HToGG_375_225","SMS_TChiWH_HToGG_375_249","SMS_TChiWH_HToGG_400_1","SMS_TChiWH_HToGG_400_25","SMS_TChiWH_HToGG_400_50","SMS_TChiWH_HToGG_400_75","SMS_TChiWH_HToGG_400_100","SMS_TChiWH_HToGG_400_125","SMS_TChiWH_HToGG_400_150","SMS_TChiWH_HToGG_400_175","SMS_TChiWH_HToGG_400_200","SMS_TChiWH_HToGG_400_225","SMS_TChiWH_HToGG_400_250","SMS_TChiWH_HToGG_400_274","SMS_TChiWH_HToGG_425_1","SMS_TChiWH_HToGG_425_25","SMS_TChiWH_HToGG_425_50","SMS_TChiWH_HToGG_425_75","SMS_TChiWH_HToGG_425_100","SMS_TChiWH_HToGG_425_125","SMS_TChiWH_HToGG_425_150","SMS_TChiWH_HToGG_425_175","SMS_TChiWH_HToGG_425_200","SMS_TChiWH_HToGG_425_225","SMS_TChiWH_HToGG_425_250","SMS_TChiWH_HToGG_425_275","SMS_TChiWH_HToGG_425_299"};//,"SMS_TChiWH_HToGG_450_1","SMS_TChiWH_HToGG_450_25","SMS_TChiWH_HToGG_450_50","SMS_TChiWH_HToGG_450_75","SMS_TChiWH_HToGG_450_100","SMS_TChiWH_HToGG_450_125","SMS_TChiWH_HToGG_450_150","SMS_TChiWH_HToGG_450_175","SMS_TChiWH_HToGG_450_200","SMS_TChiWH_HToGG_450_225","SMS_TChiWH_HToGG_450_250","SMS_TChiWH_HToGG_450_275","SMS_TChiWH_HToGG_450_300","SMS_TChiWH_HToGG_475_1","SMS_TChiWH_HToGG_475_25","SMS_TChiWH_HToGG_475_50","SMS_TChiWH_HToGG_475_75","SMS_TChiWH_HToGG_475_100","SMS_TChiWH_HToGG_475_125","SMS_TChiWH_HToGG_475_150","SMS_TChiWH_HToGG_475_175","SMS_TChiWH_HToGG_475_200","SMS_TChiWH_HToGG_475_225","SMS_TChiWH_HToGG_475_250","SMS_TChiWH_HToGG_475_275","SMS_TChiWH_HToGG_475_300","SMS_TChiWH_HToGG_500_1","SMS_TChiWH_HToGG_500_25","SMS_TChiWH_HToGG_500_50","SMS_TChiWH_HToGG_500_75","SMS_TChiWH_HToGG_500_100","SMS_TChiWH_HToGG_500_125","SMS_TChiWH_HToGG_500_150","SMS_TChiWH_HToGG_500_175","SMS_TChiWH_HToGG_500_200","SMS_TChiWH_HToGG_500_225","SMS_TChiWH_HToGG_500_250","SMS_TChiWH_HToGG_500_275","SMS_TChiWH_HToGG_500_300"};//,"SMS_TChiWH_HToGG_525_1","SMS_TChiWH_HToGG_525_25","SMS_TChiWH_HToGG_525_50","SMS_TChiWH_HToGG_525_75","SMS_TChiWH_HToGG_525_100","SMS_TChiWH_HToGG_525_125","SMS_TChiWH_HToGG_525_150","SMS_TChiWH_HToGG_525_175","SMS_TChiWH_HToGG_525_200","SMS_TChiWH_HToGG_525_225","SMS_TChiWH_HToGG_525_250","SMS_TChiWH_HToGG_525_275","SMS_TChiWH_HToGG_525_300","SMS_TChiWH_HToGG_550_1","SMS_TChiWH_HToGG_550_25","SMS_TChiWH_HToGG_550_50","SMS_TChiWH_HToGG_550_75","SMS_TChiWH_HToGG_550_100","SMS_TChiWH_HToGG_550_125","SMS_TChiWH_HToGG_550_150","SMS_TChiWH_HToGG_550_175","SMS_TChiWH_HToGG_550_200","SMS_TChiWH_HToGG_550_225","SMS_TChiWH_HToGG_550_250","SMS_TChiWH_HToGG_550_275","SMS_TChiWH_HToGG_550_300","SMS_TChiWH_HToGG_575_1","SMS_TChiWH_HToGG_575_25","SMS_TChiWH_HToGG_575_50","SMS_TChiWH_HToGG_575_75","SMS_TChiWH_HToGG_575_100","SMS_TChiWH_HToGG_575_125","SMS_TChiWH_HToGG_575_150","SMS_TChiWH_HToGG_575_175","SMS_TChiWH_HToGG_575_200","SMS_TChiWH_HToGG_575_225","SMS_TChiWH_HToGG_575_250","SMS_TChiWH_HToGG_575_275","SMS_TChiWH_HToGG_575_300","SMS_TChiWH_HToGG_600_1","SMS_TChiWH_HToGG_600_25","SMS_TChiWH_HToGG_600_50","SMS_TChiWH_HToGG_600_75","SMS_TChiWH_HToGG_600_100","SMS_TChiWH_HToGG_600_125","SMS_TChiWH_HToGG_600_150","SMS_TChiWH_HToGG_600_175","SMS_TChiWH_HToGG_600_200","SMS_TChiWH_HToGG_600_225","SMS_TChiWH_HToGG_600_250","SMS_TChiWH_HToGG_600_275","SMS_TChiWH_HToGG_600_300","SMS_TChiWH_HToGG_625_1","SMS_TChiWH_HToGG_625_25","SMS_TChiWH_HToGG_625_50","SMS_TChiWH_HToGG_625_75","SMS_TChiWH_HToGG_625_100","SMS_TChiWH_HToGG_625_125","SMS_TChiWH_HToGG_625_150","SMS_TChiWH_HToGG_625_175","SMS_TChiWH_HToGG_625_200","SMS_TChiWH_HToGG_625_225","SMS_TChiWH_HToGG_625_250","SMS_TChiWH_HToGG_625_275","SMS_TChiWH_HToGG_625_300","SMS_TChiWH_HToGG_650_1","SMS_TChiWH_HToGG_650_25","SMS_TChiWH_HToGG_650_50","SMS_TChiWH_HToGG_650_75","SMS_TChiWH_HToGG_650_100","SMS_TChiWH_HToGG_650_125","SMS_TChiWH_HToGG_650_150","SMS_TChiWH_HToGG_650_175","SMS_TChiWH_HToGG_650_200","SMS_TChiWH_HToGG_650_225","SMS_TChiWH_HToGG_650_250","SMS_TChiWH_HToGG_650_275","SMS_TChiWH_HToGG_650_300","SMS_TChiWH_HToGG_675_1","SMS_TChiWH_HToGG_675_25","SMS_TChiWH_HToGG_675_50","SMS_TChiWH_HToGG_675_75","SMS_TChiWH_HToGG_675_100","SMS_TChiWH_HToGG_675_125","SMS_TChiWH_HToGG_675_150","SMS_TChiWH_HToGG_675_175","SMS_TChiWH_HToGG_675_200","SMS_TChiWH_HToGG_675_225","SMS_TChiWH_HToGG_675_250","SMS_TChiWH_HToGG_675_275","SMS_TChiWH_HToGG_675_300","SMS_TChiWH_HToGG_700_1","SMS_TChiWH_HToGG_700_25","SMS_TChiWH_HToGG_700_50","SMS_TChiWH_HToGG_700_75","SMS_TChiWH_HToGG_700_100","SMS_TChiWH_HToGG_700_125","SMS_TChiWH_HToGG_700_150","SMS_TChiWH_HToGG_700_175","SMS_TChiWH_HToGG_700_200","SMS_TChiWH_HToGG_700_225","SMS_TChiWH_HToGG_700_250","SMS_TChiWH_HToGG_700_275","SMS_TChiWH_HToGG_700_300"};//


  //  vector<std::string> signalList_WH = {"SMS_TChiWH_HToGG_127_1","SMS_TChiWH_HToGG_150_1","SMS_TChiWH_HToGG_150_24","SMS_TChiWH_HToGG_175_1","SMS_TChiWH_HToGG_175_25","SMS_TChiWH_HToGG_175_49","SMS_TChiWH_HToGG_200_1","SMS_TChiWH_HToGG_200_25","SMS_TChiWH_HToGG_200_50","SMS_TChiWH_HToGG_200_74","SMS_TChiWH_HToGG_225_1","SMS_TChiWH_HToGG_225_25","SMS_TChiWH_HToGG_225_50","SMS_TChiWH_HToGG_225_75","SMS_TChiWH_HToGG_225_99","SMS_TChiWH_HToGG_250_1","SMS_TChiWH_HToGG_250_25","SMS_TChiWH_HToGG_250_50","SMS_TChiWH_HToGG_250_75","SMS_TChiWH_HToGG_250_100","SMS_TChiWH_HToGG_250_124","SMS_TChiWH_HToGG_275_1","SMS_TChiWH_HToGG_275_25","SMS_TChiWH_HToGG_275_50","SMS_TChiWH_HToGG_275_75","SMS_TChiWH_HToGG_275_100","SMS_TChiWH_HToGG_275_125","SMS_TChiWH_HToGG_275_149","SMS_TChiWH_HToGG_300_1","SMS_TChiWH_HToGG_300_25","SMS_TChiWH_HToGG_300_50","SMS_TChiWH_HToGG_300_75","SMS_TChiWH_HToGG_300_100","SMS_TChiWH_HToGG_300_125","SMS_TChiWH_HToGG_300_150","SMS_TChiWH_HToGG_300_174","SMS_TChiWH_HToGG_325_1","SMS_TChiWH_HToGG_325_25","SMS_TChiWH_HToGG_325_50","SMS_TChiWH_HToGG_325_75","SMS_TChiWH_HToGG_325_100","SMS_TChiWH_HToGG_325_125","SMS_TChiWH_HToGG_325_150","SMS_TChiWH_HToGG_325_175","SMS_TChiWH_HToGG_325_199","SMS_TChiWH_HToGG_350_1","SMS_TChiWH_HToGG_350_25","SMS_TChiWH_HToGG_350_50","SMS_TChiWH_HToGG_350_75","SMS_TChiWH_HToGG_350_100","SMS_TChiWH_HToGG_350_125","SMS_TChiWH_HToGG_350_150","SMS_TChiWH_HToGG_350_175","SMS_TChiWH_HToGG_350_200","SMS_TChiWH_HToGG_350_224","SMS_TChiWH_HToGG_375_1","SMS_TChiWH_HToGG_375_25","SMS_TChiWH_HToGG_375_50","SMS_TChiWH_HToGG_375_75","SMS_TChiWH_HToGG_375_100","SMS_TChiWH_HToGG_375_125","SMS_TChiWH_HToGG_375_150","SMS_TChiWH_HToGG_375_175","SMS_TChiWH_HToGG_375_200","SMS_TChiWH_HToGG_375_225","SMS_TChiWH_HToGG_375_249","SMS_TChiWH_HToGG_400_1","SMS_TChiWH_HToGG_400_25","SMS_TChiWH_HToGG_400_50","SMS_TChiWH_HToGG_400_75","SMS_TChiWH_HToGG_400_100","SMS_TChiWH_HToGG_400_125","SMS_TChiWH_HToGG_400_150","SMS_TChiWH_HToGG_400_175","SMS_TChiWH_HToGG_400_200","SMS_TChiWH_HToGG_400_225","SMS_TChiWH_HToGG_400_250","SMS_TChiWH_HToGG_400_274","SMS_TChiWH_HToGG_425_1","SMS_TChiWH_HToGG_425_25","SMS_TChiWH_HToGG_425_50","SMS_TChiWH_HToGG_425_75","SMS_TChiWH_HToGG_425_100","SMS_TChiWH_HToGG_425_125","SMS_TChiWH_HToGG_425_150","SMS_TChiWH_HToGG_425_175","SMS_TChiWH_HToGG_425_200","SMS_TChiWH_HToGG_425_225","SMS_TChiWH_HToGG_425_250","SMS_TChiWH_HToGG_425_275","SMS_TChiWH_HToGG_425_299","SMS_TChiWH_HToGG_450_1","SMS_TChiWH_HToGG_450_25","SMS_TChiWH_HToGG_450_50","SMS_TChiWH_HToGG_450_75","SMS_TChiWH_HToGG_450_100","SMS_TChiWH_HToGG_450_125","SMS_TChiWH_HToGG_450_150","SMS_TChiWH_HToGG_450_175","SMS_TChiWH_HToGG_450_200","SMS_TChiWH_HToGG_450_225","SMS_TChiWH_HToGG_450_250","SMS_TChiWH_HToGG_450_275","SMS_TChiWH_HToGG_450_300","SMS_TChiWH_HToGG_475_1","SMS_TChiWH_HToGG_475_25","SMS_TChiWH_HToGG_475_50","SMS_TChiWH_HToGG_475_75","SMS_TChiWH_HToGG_475_100","SMS_TChiWH_HToGG_475_125","SMS_TChiWH_HToGG_475_150","SMS_TChiWH_HToGG_475_175","SMS_TChiWH_HToGG_475_200","SMS_TChiWH_HToGG_475_225","SMS_TChiWH_HToGG_475_250","SMS_TChiWH_HToGG_475_275","SMS_TChiWH_HToGG_475_300","SMS_TChiWH_HToGG_500_1","SMS_TChiWH_HToGG_500_25","SMS_TChiWH_HToGG_500_50","SMS_TChiWH_HToGG_500_75","SMS_TChiWH_HToGG_500_100","SMS_TChiWH_HToGG_500_125","SMS_TChiWH_HToGG_500_150","SMS_TChiWH_HToGG_500_175","SMS_TChiWH_HToGG_500_200","SMS_TChiWH_HToGG_500_225","SMS_TChiWH_HToGG_500_250","SMS_TChiWH_HToGG_500_275","SMS_TChiWH_HToGG_500_300","SMS_TChiWH_HToGG_525_1","SMS_TChiWH_HToGG_525_25","SMS_TChiWH_HToGG_525_50","SMS_TChiWH_HToGG_525_75","SMS_TChiWH_HToGG_525_100","SMS_TChiWH_HToGG_525_125","SMS_TChiWH_HToGG_525_150","SMS_TChiWH_HToGG_525_175","SMS_TChiWH_HToGG_525_200","SMS_TChiWH_HToGG_525_225","SMS_TChiWH_HToGG_525_250","SMS_TChiWH_HToGG_525_275","SMS_TChiWH_HToGG_525_300","SMS_TChiWH_HToGG_550_1","SMS_TChiWH_HToGG_550_25","SMS_TChiWH_HToGG_550_50","SMS_TChiWH_HToGG_550_75","SMS_TChiWH_HToGG_550_100","SMS_TChiWH_HToGG_550_125","SMS_TChiWH_HToGG_550_150","SMS_TChiWH_HToGG_550_175","SMS_TChiWH_HToGG_550_200","SMS_TChiWH_HToGG_550_225","SMS_TChiWH_HToGG_550_250","SMS_TChiWH_HToGG_550_275","SMS_TChiWH_HToGG_550_300","SMS_TChiWH_HToGG_575_1","SMS_TChiWH_HToGG_575_25","SMS_TChiWH_HToGG_575_50","SMS_TChiWH_HToGG_575_75","SMS_TChiWH_HToGG_575_100","SMS_TChiWH_HToGG_575_125","SMS_TChiWH_HToGG_575_150","SMS_TChiWH_HToGG_575_175","SMS_TChiWH_HToGG_575_200","SMS_TChiWH_HToGG_575_225","SMS_TChiWH_HToGG_575_250","SMS_TChiWH_HToGG_575_275","SMS_TChiWH_HToGG_575_300","SMS_TChiWH_HToGG_600_1","SMS_TChiWH_HToGG_600_25","SMS_TChiWH_HToGG_600_50","SMS_TChiWH_HToGG_600_75","SMS_TChiWH_HToGG_600_100","SMS_TChiWH_HToGG_600_125","SMS_TChiWH_HToGG_600_150","SMS_TChiWH_HToGG_600_175","SMS_TChiWH_HToGG_600_200","SMS_TChiWH_HToGG_600_225","SMS_TChiWH_HToGG_600_250","SMS_TChiWH_HToGG_600_275","SMS_TChiWH_HToGG_600_300","SMS_TChiWH_HToGG_625_1","SMS_TChiWH_HToGG_625_25","SMS_TChiWH_HToGG_625_50","SMS_TChiWH_HToGG_625_75","SMS_TChiWH_HToGG_625_100","SMS_TChiWH_HToGG_625_125","SMS_TChiWH_HToGG_625_150","SMS_TChiWH_HToGG_625_175","SMS_TChiWH_HToGG_625_200","SMS_TChiWH_HToGG_625_225","SMS_TChiWH_HToGG_625_250","SMS_TChiWH_HToGG_625_275","SMS_TChiWH_HToGG_625_300","SMS_TChiWH_HToGG_650_1","SMS_TChiWH_HToGG_650_25","SMS_TChiWH_HToGG_650_50","SMS_TChiWH_HToGG_650_75","SMS_TChiWH_HToGG_650_100","SMS_TChiWH_HToGG_650_125","SMS_TChiWH_HToGG_650_150","SMS_TChiWH_HToGG_650_175","SMS_TChiWH_HToGG_650_200","SMS_TChiWH_HToGG_650_225","SMS_TChiWH_HToGG_650_250","SMS_TChiWH_HToGG_650_275","SMS_TChiWH_HToGG_650_300","SMS_TChiWH_HToGG_675_1","SMS_TChiWH_HToGG_675_25","SMS_TChiWH_HToGG_675_50","SMS_TChiWH_HToGG_675_75","SMS_TChiWH_HToGG_675_100","SMS_TChiWH_HToGG_675_125","SMS_TChiWH_HToGG_675_150","SMS_TChiWH_HToGG_675_175","SMS_TChiWH_HToGG_675_200","SMS_TChiWH_HToGG_675_225","SMS_TChiWH_HToGG_675_250","SMS_TChiWH_HToGG_675_275","SMS_TChiWH_HToGG_675_300","SMS_TChiWH_HToGG_700_1","SMS_TChiWH_HToGG_700_25","SMS_TChiWH_HToGG_700_50","SMS_TChiWH_HToGG_700_75","SMS_TChiWH_HToGG_700_100","SMS_TChiWH_HToGG_700_125","SMS_TChiWH_HToGG_700_150","SMS_TChiWH_HToGG_700_175","SMS_TChiWH_HToGG_700_200","SMS_TChiWH_HToGG_700_225","SMS_TChiWH_HToGG_700_250","SMS_TChiWH_HToGG_700_275","SMS_TChiWH_HToGG_700_300"};





  // vector<std::string> signalList_WH = { "SMS_TChiWH_HToGG_127_1","SMS_TChiWH_HToGG_150_1","SMS_TChiWH_HToGG_150_24","SMS_TChiWH_HToGG_175_1","SMS_TChiWH_HToGG_175_25","SMS_TChiWH_HToGG_175_49","SMS_TChiWH_HToGG_200_1","SMS_TChiWH_HToGG_200_25","SMS_TChiWH_HToGG_200_50","SMS_TChiWH_HToGG_200_74" };


  //  vector<std::string> signalList_HZ = { "SMS_TChiHZ_HToGG_m127", "SMS_TChiHZ_HToGG_m150" };
  // vector<std::string> signalList_HZ = { "SMS_TChiHZ_HToGG_m127", "SMS_TChiHZ_HToGG_m150", "SMS_TChiHZ_HToGG_m175", "SMS_TChiHZ_HToGG_m200" };
  
  vector<std::string> signalList_HZ = { "SMS_TChiHZ_HToGG_m127", "SMS_TChiHZ_HToGG_m150", "SMS_TChiHZ_HToGG_m175", "SMS_TChiHZ_HToGG_m200", "SMS_TChiHZ_HToGG_m225", "SMS_TChiHZ_HToGG_m250", "SMS_TChiHZ_HToGG_m275", "SMS_TChiHZ_HToGG_m300", "SMS_TChiHZ_HToGG_m325", "SMS_TChiHZ_HToGG_m350", "SMS_TChiHZ_HToGG_m375", "SMS_TChiHZ_HToGG_m400", "SMS_TChiHZ_HToGG_m425", "SMS_TChiHZ_HToGG_m450", "SMS_TChiHZ_HToGG_m475", "SMS_TChiHZ_HToGG_m500"};//, "SMS_TChiHZ_HToGG_m525", "SMS_TChiHZ_HToGG_m550", "SMS_TChiHZ_HToGG_m575", "SMS_TChiHZ_HToGG_m600", "SMS_TChiHZ_HToGG_m625", "SMS_TChiHZ_HToGG_m650", "SMS_TChiHZ_HToGG_m675", "SMS_TChiHZ_HToGG_m700", "SMS_TChiHZ_HToGG_m725", "SMS_TChiHZ_HToGG_m750", "SMS_TChiHZ_HToGG_m775", "SMS_TChiHZ_HToGG_m800", "SMS_TChiHZ_HToGG_m825", "SMS_TChiHZ_HToGG_m850", "SMS_TChiHZ_HToGG_m875", "SMS_TChiHZ_HToGG_m900", "SMS_TChiHZ_HToGG_m925", "SMS_TChiHZ_HToGG_m950", "SMS_TChiHZ_HToGG_m975", "SMS_TChiHZ_HToGG_m1000" };

  //  vector<std::string> signalList_HH_HH0p25 = { "SMS_TChiHZ_HToGG_m127_HH0p25", "SMS_TChiHZ_HToGG_m150_HH0p25"};
  //  vector<std::string> signalList_HH_HH0p25 = { "SMS_TChiHZ_HToGG_m127_HH0p25", "SMS_TChiHZ_HToGG_m150_HH0p25", "SMS_TChiHZ_HToGG_m175_HH0p25", "SMS_TChiHZ_HToGG_m200_HH0p25"};

  vector<std::string> signalList_HH_HH0p25 = { "SMS_TChiHZ_HToGG_m127_HH0p25", "SMS_TChiHZ_HToGG_m150_HH0p25", "SMS_TChiHZ_HToGG_m175_HH0p25", "SMS_TChiHZ_HToGG_m200_HH0p25", "SMS_TChiHZ_HToGG_m225_HH0p25", "SMS_TChiHZ_HToGG_m250_HH0p25", "SMS_TChiHZ_HToGG_m275_HH0p25", "SMS_TChiHZ_HToGG_m300_HH0p25", "SMS_TChiHZ_HToGG_m325_HH0p25", "SMS_TChiHZ_HToGG_m350_HH0p25", "SMS_TChiHZ_HToGG_m375_HH0p25", "SMS_TChiHZ_HToGG_m400_HH0p25", "SMS_TChiHZ_HToGG_m425_HH0p25", "SMS_TChiHZ_HToGG_m450_HH0p25", "SMS_TChiHZ_HToGG_m475_HH0p25", "SMS_TChiHZ_HToGG_m500_HH0p25"}; //, "SMS_TChiHZ_HToGG_m525_HH0p25", "SMS_TChiHZ_HToGG_m550_HH0p25", "SMS_TChiHZ_HToGG_m575_HH0p25", "SMS_TChiHZ_HToGG_m600_HH0p25", "SMS_TChiHZ_HToGG_m625_HH0p25", "SMS_TChiHZ_HToGG_m650_HH0p25", "SMS_TChiHZ_HToGG_m675_HH0p25", "SMS_TChiHZ_HToGG_m700_HH0p25", "SMS_TChiHZ_HToGG_m725_HH0p25", "SMS_TChiHZ_HToGG_m750_HH0p25", "SMS_TChiHZ_HToGG_m775_HH0p25", "SMS_TChiHZ_HToGG_m800_HH0p25", "SMS_TChiHZ_HToGG_m825_HH0p25", "SMS_TChiHZ_HToGG_m850_HH0p25", "SMS_TChiHZ_HToGG_m875_HH0p25", "SMS_TChiHZ_HToGG_m900_HH0p25", "SMS_TChiHZ_HToGG_m925_HH0p25", "SMS_TChiHZ_HToGG_m950_HH0p25", "SMS_TChiHZ_HToGG_m975_HH0p25", "SMS_TChiHZ_HToGG_m1000_HH0p25" };


  //  vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127"};
  //  vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127", "SMS_TChiHH_HToGG_m150", "SMS_TChiHH_HToGG_m175", "SMS_TChiHH_HToGG_m200"};

  vector<std::string> signalList_HH = { "SMS_TChiHH_HToGG_m127", "SMS_TChiHH_HToGG_m150", "SMS_TChiHH_HToGG_m175", "SMS_TChiHH_HToGG_m200", "SMS_TChiHH_HToGG_m225", "SMS_TChiHH_HToGG_m250", "SMS_TChiHH_HToGG_m275", "SMS_TChiHH_HToGG_m300", "SMS_TChiHH_HToGG_m325", "SMS_TChiHH_HToGG_m350", "SMS_TChiHH_HToGG_m375", "SMS_TChiHH_HToGG_m400", "SMS_TChiHH_HToGG_m425", "SMS_TChiHH_HToGG_m450", "SMS_TChiHH_HToGG_m475", "SMS_TChiHH_HToGG_m500", "SMS_TChiHH_HToGG_m525", "SMS_TChiHH_HToGG_m550", "SMS_TChiHH_HToGG_m575", "SMS_TChiHH_HToGG_m600", "SMS_TChiHH_HToGG_m625", "SMS_TChiHH_HToGG_m650", "SMS_TChiHH_HToGG_m675", "SMS_TChiHH_HToGG_m700", "SMS_TChiHH_HToGG_m725", "SMS_TChiHH_HToGG_m750", "SMS_TChiHH_HToGG_m775", "SMS_TChiHH_HToGG_m800", "SMS_TChiHH_HToGG_m825", "SMS_TChiHH_HToGG_m850", "SMS_TChiHH_HToGG_m875", "SMS_TChiHH_HToGG_m900", "SMS_TChiHH_HToGG_m925", "SMS_TChiHH_HToGG_m950", "SMS_TChiHH_HToGG_m975", "SMS_TChiHH_HToGG_m1000" };

  std::vector<std::string> regionSelection;

  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && ((h_pt/h_mass) <  0.1)) " );
  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.1<=(h_pt/h_mass) && (h_pt/h_mass) < 0.2)) " );

  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.2<=(h_pt/h_mass) && (h_pt/h_mass) < 0.3)) " );
  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.3<=(h_pt/h_mass) && (h_pt/h_mass) < 0.4)) " );

  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.4<=(h_pt/h_mass) && (h_pt/h_mass) < 0.6)) " );
  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.6<=(h_pt/h_mass) && (h_pt/h_mass) < 0.8)) " );

  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (0.8<=(h_pt/h_mass) && (h_pt/h_mass) < 1.0)) " );
  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (1.0<=(h_pt/h_mass) && (h_pt/h_mass) < 1.2)) " );

  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && (1.2<=(h_pt/h_mass) && (h_pt/h_mass) < 1.4)) " );
  // regionSelection.push_back(" (gg_nJets==0 && gg_nBJets>= 0 && ((h_pt/h_mass) >= 1.4)) " );

  // int n0j = 10;

  int n0j = 2;
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
    regionSelection[i] +=  "  && ( (!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )";
  }

  
  regionSelection.push_back(" ( is1El) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  regionSelection.push_back(" ( is1Mu) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(isDiLepZ) )"  );
  regionSelection.push_back(" ( is1El) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  regionSelection.push_back(" ( is1Mu) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(isDiLepZ) )"  );

  regionSelection.push_back(" ( isDiLepZ ) "     );

  regionSelection.push_back(" ( isDiBH ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBZ) && !(isDiLepZ) )" );
  regionSelection.push_back(" ( isDiBZ ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBH) && !(isDiLepZ) )" );
  regionSelection.push_back(" ( isDiBH ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBZ) && !(isDiLepZ) )" );
  regionSelection.push_back(" ( isDiBZ ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBH) && !(isDiLepZ) )" );


  // old selection  
  // regionSelection.push_back(" ( is1El) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  // regionSelection.push_back(" ( is1Mu) && ((h_pt/h_mass) <  0.8)   && ( (!isDiBH && !isDiBZ) && !(is1El) && !(isDiLepZ) )"  );
  // regionSelection.push_back(" ( is1El) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(is1Mu) && !(isDiLepZ) )"  );
  // regionSelection.push_back(" ( is1Mu) && ((h_pt/h_mass) >= 0.8)   && ( (!isDiBH && !isDiBZ) && !(is1El) && !(isDiLepZ) )"  );

  // regionSelection.push_back(" ( isDiLepZ ) && ( (!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) )"     );

  // regionSelection.push_back(" ( isDiBH ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  // regionSelection.push_back(" ( isDiBZ ) && ((h_pt/h_mass) <  0.8) && ( (!isDiBH) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  // regionSelection.push_back(" ( isDiBH ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  // regionSelection.push_back(" ( isDiBZ ) && ((h_pt/h_mass) >= 0.8) && ( (!isDiBH) && !(is1El) && !(is1Mu) && !(isDiLepZ) )" );
  


  std::vector<std::string> regionNames;
  regionNames.push_back("j0_b0toInf_pT0");
  regionNames.push_back("j0_b0toInf_pT1");
 
  // regionNames.push_back("j0_b0toInf_pT2");
  // regionNames.push_back("j0_b0toInf_pT3");
  // regionNames.push_back("j0_b0toInf_pT4");
  // regionNames.push_back("j0_b0toInf_pT5");
  // regionNames.push_back("j0_b0toInf_pT6");
  // regionNames.push_back("j0_b0toInf_pT7");
  // regionNames.push_back("j0_b0toInf_pT8");
  // regionNames.push_back("j0_b0toInf_pT9");

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

  //  std::vector<int> mt2_bins = { 0 };
  
  std::vector<int> mt2_bins = { 0, 30 };
  unsigned int regionLength = regionSelection.size();


  unsigned int nMT2split = 15; // 13 if 1 lep high pt is not split in mt2

    for( unsigned int j=0; j<mt2_bins.size(); j++){ 
      for( unsigned int i=0; i < regionLength; i++){// starts at 2 because 0j bins have no "real" second object to make mt2 with
	if( i>nMT2split + n0j   ) continue;
	if( j==0  ){
	  if ( i<n0j ) continue;

	  regionSelection.push_back( regionSelection[i] );
	  regionSelection[i] +=  "  && ( hgg_mt2<30   )";

	  regionNames.push_back( regionNames[i] );
	  regionNames[i] +=  "_mt2_0";
	  //	std::cout << "selection 0 " << regionSelection[i] << std::endl;
	}else{
	  if ( i>(nMT2split )  ) continue;
	  //  	if ( i>(nMT2split - n0j) ) continue;
	  regionSelection[ regionLength + i ] +=  "  && ( hgg_mt2>=30   )";
	  regionNames[ regionLength + i ] +=  "_mt2_30";
	  //	std::cout << "selection 1 " << regionSelection[regionLength + i ] << std::endl;
	}

      }

    }



  for( unsigned int i=0; i < regionSelection.size(); i++){
    std::cout << regionNames[ i ] << std::endl;
    // std::cout << regionNames[ i ] << " with selection" << regionSelection[i] << std::endl;
  }



  //2016 data and MC
  MT2Analysis<MT2EstimateTree>* data     = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ggDir.c_str() ), "diPhoton_data");
  MT2Analysis<MT2EstimateTree>* qcd      = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_qcd.root", ggDir.c_str()  ), "qcd");
  MT2Analysis<MT2EstimateTree>* diPhoton = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diPhoton.root", ggDir.c_str() ), "diPhoton");
  MT2Analysis<MT2EstimateTree>* gjets    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_gjets.root", ggDir.c_str() ), "gjets");
  MT2Analysis<MT2EstimateTree>* higgs    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir.c_str() ), "higgs");
  


  //2017 data and MC  
  MT2Analysis<MT2EstimateTree>* data_2017;
  MT2Analysis<MT2EstimateTree>* qcd_2017;     
  MT2Analysis<MT2EstimateTree>* diPhoton_2017;
  MT2Analysis<MT2EstimateTree>* gjets_2017;
  MT2Analysis<MT2EstimateTree>* higgs_2017;  
  

  if( doComb ) {
    data_2017     = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ggDir_2017.c_str() ), "diPhoton_data");
    qcd_2017      = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_qcd.root", ggDir_2017.c_str()  ), "qcd");
    diPhoton_2017 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diPhoton.root", ggDir_2017.c_str() ), "diPhoton");
    gjets_2017    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_gjets.root", ggDir_2017.c_str() ), "gjets");
    higgs_2017    = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir_2017.c_str() ), "higgs");
  
  }

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

  // 2017 version
  vector< MT2Analysis<MT2EstimateTree>* > signals_2017;
  if( doComb ) 
    for( unsigned int i=0; i<signalList.size(); i++){
      MT2Analysis<MT2EstimateTree>* sigTemp_2017 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir_2017.c_str() ), signalList[i]  );
      signals_2017.push_back( sigTemp_2017 );
    }

  vector< MT2Analysis<MT2EstimateTree>* > signals_HH_HH0p25_2017;
  if( doComb ) 
    if( doHZ )
      for( unsigned int i=0; i<signalList_HH_HH0p25.size(); i++){
	MT2Analysis<MT2EstimateTree>* sigTemp_HH_HH0p25_2017 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir_2017.c_str() ), signalList_HH_HH0p25[i]  );
	signals_HH_HH0p25_2017.push_back( sigTemp_HH_HH0p25_2017 );
      }


  std::string outFileName( Form("%s/ws_%s", outputdir.c_str(), sigName.c_str() ));
  //  if( doLLbin )       outFileName += "_ll";
  //  if( doMT2binning )  outFileName += "_0jsplit_mt2";
  //  if( doSingleLbin )  outFileName += "_onelPT";
  //  if( doBBbin )       outFileName += "_bb";
  if( doComb )       outFileName += "_comb";
  //  if( mt2_bins.size() > 1 ) outFileName += "_0jsplit_mt2_30";
  if( mt2_bins.size() > 1 ) outFileName += "_mt2_30";
  
  TFile* dataFile = new TFile(  Form("%s.root",outFileName.c_str()) , "RECREATE" );
  RooWorkspace ws_sig ("ws_sig");
  RooWorkspace ws_bg  ("ws_bg");
  RooWorkspace ws_data("ws_data");

  std::string dataAmountName( Form("%s/dataAmount",  outputdir.c_str() ));
  // if( doLLbin )             dataAmountName += "_ll";
  // if( doMT2binning )        dataAmountName += "_mt2";
  // if( doSingleLbin )        dataAmountName += "_onelPT";
  // if( doBBbin )             dataAmountName += "_bb";
  if( doComb )       dataAmountName += "_comb";
  if( mt2_bins.size() > 1 ) dataAmountName += "_mt2_30";
  std::ofstream dataAmount( Form("%s.txt", dataAmountName.c_str() ) );

  std::string bgSumAmountName( Form("%s/bgSumAmount",  outputdir.c_str() ));
  // if( doLLbin )             bgSumAmountName += "_ll";
  // if( doMT2binning )        bgSumAmountName += "_mt2";
  // if( doSingleLbin )        bgSumAmountName += "_onelPT";
  // if( doBBbin )             bgSumAmountName += "_bb";
  if( doComb )             bgSumAmountName += "_comb";
  if( mt2_bins.size() > 1 ) bgSumAmountName += "_mt2_30";
  std::ofstream bgSumAmount( Form("%s.txt", bgSumAmountName.c_str() ) );

  std::string isrName( Form("%s/isr",  outputdir.c_str() ));
  // if( doLLbin )    isrName += "_ll";
  // if( doMT2binning )    isrName += "_mt2";
  // if( doSingleLbin )    isrName += "_onelPT";
  // if( doBBbin )    isrName += "_bb";
  if( doComb )    isrName += "_comb";
  if( mt2_bins.size() > 1 ) isrName += "_mt2_30";
  std::ofstream isr( Form("%s.txt", isrName.c_str()) );

  std::string lepsf_higgsName( Form("%s/lepsf_higgs",  outputdir.c_str() ));
  if( doComb )    lepsf_higgsName += "_comb";
  if( mt2_bins.size() > 1 ) lepsf_higgsName += "_mt2_30";
  std::ofstream lepsf_higgs( Form("%s.txt", lepsf_higgsName.c_str()) );

  std::string btagsf_heavy_higgsName( Form("%s/btagsf_heavy_higgs",  outputdir.c_str() ));
  if( doComb )    btagsf_heavy_higgsName += "_comb";
  if( mt2_bins.size() > 1 ) btagsf_heavy_higgsName += "_mt2_30";
  std::ofstream btagsf_heavy_higgs( Form("%s.txt", btagsf_heavy_higgsName.c_str()) );

  std::string btagsf_light_higgsName( Form("%s/btagsf_light_higgs",  outputdir.c_str() ));
  if( doComb )    btagsf_light_higgsName += "_comb";
  if( mt2_bins.size() > 1 ) btagsf_light_higgsName += "_mt2_30";
  std::ofstream btagsf_light_higgs( Form("%s.txt", btagsf_light_higgsName.c_str()) );

  std::string lepsf_higgs_2017Name( Form("%s/lepsf_higgs_2017",  outputdir.c_str() ));
  if( doComb )    lepsf_higgs_2017Name += "_comb";
  if( mt2_bins.size() > 1 ) lepsf_higgs_2017Name += "_mt2_30";
  std::ofstream lepsf_higgs_2017( Form("%s.txt", lepsf_higgs_2017Name.c_str()) );

  std::string btagsf_heavy_higgs_2017Name( Form("%s/btagsf_heavy_higgs_2017",  outputdir.c_str() ));
  if( doComb )    btagsf_heavy_higgs_2017Name += "_comb";
  if( mt2_bins.size() > 1 ) btagsf_heavy_higgs_2017Name += "_mt2_30";
  std::ofstream btagsf_heavy_higgs_2017( Form("%s.txt", btagsf_heavy_higgs_2017Name.c_str()) );

  std::string btagsf_light_higgs_2017Name( Form("%s/btagsf_light_higgs_2017",  outputdir.c_str() ));
  if( doComb )    btagsf_light_higgs_2017Name += "_comb";
  if( mt2_bins.size() > 1 ) btagsf_light_higgs_2017Name += "_mt2_30";
  std::ofstream btagsf_light_higgs_2017( Form("%s.txt", btagsf_light_higgs_2017Name.c_str()) );


  std::string sigAmountName( Form("%s/sigAmount_%s",  outputdir.c_str(), sigName.c_str() ));
  // if( doLLbin )    sigAmountName += "_ll";
  // if( doMT2binning )    sigAmountName += "_mt2";
  // if( doSingleLbin )    sigAmountName += "_onelPT";
  // if( doBBbin )    sigAmountName += "_bb";
  if( doComb )    sigAmountName += "_comb";
  if( mt2_bins.size() > 1 ) sigAmountName += "_mt2_30";
  std::ofstream sigAmount( Form("%s.txt", sigAmountName.c_str()) );

  std::string higgsAmountName( Form("%s/higgsAmount",  outputdir.c_str() ));
  // if( doLLbin )    higgsAmountName += "_ll";
  // if( doMT2binning )    higgsAmountName += "_mt2";
  // if( doSingleLbin )    higgsAmountName += "_onelPT";
  // if( doBBbin )    higgsAmountName += "_bb";
  if( doComb )    higgsAmountName += "_comb";
  if( mt2_bins.size() > 1 ) higgsAmountName += "_mt2_30";
  std::ofstream higgsAmount( Form("%s.txt", higgsAmountName.c_str()) );
        
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
    //  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    // int mt2BinsIt = 0;
    // if( !doMT2binning )
    //   mt2BinsIt = 1;

    // for( ; mt2BinsIt<2; mt2BinsIt++){

    string regionSaveName =  regionNames[i];

    //INITIAL SELECTION 
    string treeSel = "((( ptGamma0/h_mass) > 1./3) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )" ;

    //    string treeSel = "((( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 && ( (h_mass>=105 &&h_mass<120  )||( h_mass>=130 &&h_mass<160 )) )" ;
    //string treeSel = "((( ptGamma0) > 40) && ((ptGamma1 )> 25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 &&( (h_mass>=105 &&h_mass<120  )||( h_mass>=130 &&h_mass<160 ) ) )" ;

    //string treeSel = "((( ptGamma0) > 40) && ((ptGamma1 )> 25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )" ;
    
    


    treeSel += " && " + regionSelection[i];
    //string treeSelWideEta = "((( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=2.4  && fabs(etaGamma1)<=2.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )" ;
      

    TH2D* h2D_isr    = new TH2D( Form("h2D_isr_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );
    TH2D* h2D_isr_err = new TH2D(Form("h2D_isr_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );

    TH2D* h2D_isr_2017    = new TH2D( Form("h2D_isr_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );
    TH2D* h2D_isr_2017_err = new TH2D(Form("h2D_isr_2017_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );

    TH2D* h2D_lepsf    = new TH2D( Form("h2D_lepsf_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );
    TH2D* h2D_lepsf_err = new TH2D(Form("h2D_lepsf_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );

    TH2D* h2D_lepsf_2017    = new TH2D( Form("h2D_lepsf_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );
    TH2D* h2D_lepsf_2017_err = new TH2D(Form("h2D_lepsf_2017_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );

    TH2D* h2D_btagsf_heavy    = new TH2D( Form("h2D_btagsf_heavy_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );
    TH2D* h2D_btagsf_heavy_err = new TH2D(Form("h2D_btagsf_heavy_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );

    TH2D* h2D_btagsf_heavy_2017    = new TH2D( Form("h2D_btagsf_heavy_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );
    TH2D* h2D_btagsf_heavy_2017_err = new TH2D(Form("h2D_btagsf_heavy_2017_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );

    TH2D* h2D_btagsf_light    = new TH2D( Form("h2D_btagsf_light_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );
    TH2D* h2D_btagsf_light_err = new TH2D(Form("h2D_btagsf_light_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );

    TH2D* h2D_btagsf_light_2017    = new TH2D( Form("h2D_btagsf_light_2017_%s", regionSaveName.c_str() ),    "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );
    TH2D* h2D_btagsf_light_2017_err = new TH2D(Form("h2D_btagsf_light_2017_err_%s", regionSaveName.c_str() ), "", 36, 125-12.5, 1012.5, 13 , -12.5, 312.5 );

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

    //comb with 2017
    TTree* tree_data_ini_2017;
    TTree* tree_higgs_ini_2017;
    TTree* tree_qcd_ini_2017;
    TTree* tree_diPhoton_ini_2017;
    TTree* tree_gjets_ini_2017;

    TTree* tree_data_2017;
    TTree* tree_higgs_2017;
    TTree* tree_qcd_2017;
    TTree* tree_diPhoton_2017;
    TTree* tree_gjets_2017;

    if( doComb ){
      tree_data_ini_2017 = data_2017->get(thisRegion)->tree;
      tree_higgs_ini_2017 = higgs_2017->get(thisRegion)->tree;
      tree_qcd_ini_2017 = qcd_2017->get(thisRegion)->tree;
      tree_diPhoton_ini_2017 = diPhoton_2017->get(thisRegion)->tree;
      tree_gjets_ini_2017 = gjets_2017->get(thisRegion)->tree;

      tree_data_2017 = tree_data_ini_2017->CopyTree(treeSel.c_str() );
      tree_higgs_2017 = tree_higgs_ini_2017->CopyTree(treeSel.c_str() );
      tree_qcd_2017 = tree_qcd_ini_2017->CopyTree(treeSel.c_str() );
      tree_diPhoton_2017 = tree_diPhoton_ini_2017->CopyTree(treeSel.c_str() );
      tree_gjets_2017 = tree_gjets_ini_2017->CopyTree(treeSel.c_str() );
    }


   
    //Import 2017
    RooDataSet higgs_ds(Form("higgs_2016_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs) );
    RooDataSet bg_ds(Form("bg_2016_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton) ) ;
    RooDataSet gjets_ds("gjets_2016_ds", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets) ) ;
    RooDataSet qcd_ds("qcd_2016_ds", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd) ) ;

    RooDataSet data_ds( Form("data_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data) ) ;
   
    if(doComb){
      //Import 2017
      RooDataSet higgs_ds_2017(Form("higgs_2017_125_13TeV_%s", regionSaveName.c_str() ), "Higgs", vars, WeightVar(weight), Import(*tree_higgs_2017) );
      RooDataSet bg_ds_2017(Form("bg_2017_125_13TeV_%s", regionSaveName.c_str() ), "Background",  vars, WeightVar(weight), Import(*tree_diPhoton_2017) ) ;
      RooDataSet gjets_ds_2017("gjets_ds_2017", "gjets_ds",                                       vars, WeightVar(weight), Import(*tree_gjets_2017) ) ;
      RooDataSet qcd_ds_2017("qcd_ds_2017", "qcd_ds",                                             vars, WeightVar(weight), Import(*tree_qcd_2017) ) ;
      RooDataSet data_ds_2017( Form("data_2017_125_13TeV_%s", regionSaveName.c_str() ), "Data",   vars, WeightVar(weight), Import(*tree_data_2017) ) ;
   
      bg_ds_2017.append( qcd_ds_2017 );   //append bg data into one
      bg_ds_2017.append( gjets_ds_2017 ); //append bg data into one

      RooDataSet higgs_ds_130_2017(  higgs_ds_2017, Form("higgs_2017_130_13TeV_%s", regionSaveName.c_str() ) );
      RooDataSet bg_ds_130_2017   (  bg_ds_2017,    Form("bg_2017_130_13TeV_%s",    regionSaveName.c_str() ) );
      //    RooDataSet data_ds_130_2017 (  data_ds_2017,  Form("data_2017_130_13TeV_%s",  regionSaveName.c_str() ) );
      RooDataSet higgs_ds_125_2017(  higgs_ds_2017, Form("higgs_2017_125_13TeV_%s", regionSaveName.c_str() ) );
      RooDataSet bg_ds_125_2017   (  bg_ds_2017,    Form("bg_2017_125_13TeV_%s",    regionSaveName.c_str() ) );
      //    RooDataSet data_ds_125_2017 (  data_ds_2017,  Form("data_2017_125_13TeV_%s",  regionSaveName.c_str() ) );
      RooDataSet higgs_ds_120_2017(  higgs_ds_2017, Form("higgs_2017_120_13TeV_%s", regionSaveName.c_str() ) );
      RooDataSet bg_ds_120_2017   (  bg_ds_2017,    Form("bg_2017_120_13TeV_%s",    regionSaveName.c_str() ) );
      //    RooDataSet data_ds_120_2017 (  data_ds_2017,  Form("data_2017_120_13TeV_%s",  regionSaveName.c_str() ) );
      ((RooRealVar*)higgs_ds_130_2017.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
      //    ((RooRealVar*)data_ds_130_2017.addColumn( hgg_mass_130_2017 ))->setRange(100.,180.); 
      ((RooRealVar*)higgs_ds_125_2017.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
      //    ((RooRealVar*)data_ds_125_2017.addColumn( hgg_mass_125_2017 ))->setRange(100.,180.); 
      ((RooRealVar*)higgs_ds_120_2017.addColumn( hgg_mass_120 ))->setRange(100.,180.);
      //    ((RooRealVar*)data_ds_120_2017.addColumn( hgg_mass_120_2017 ))->setRange(100.,180.);

      //  ws_data.import( data_ds_125 );
      ws_bg.import( bg_ds_125_2017   );
      ws_sig.import( higgs_ds_125_2017  );
      //  ws_data.import( data_ds_130 );
      ws_bg.import( bg_ds_130_2017   );
      ws_sig.import( higgs_ds_130_2017  );
      //  ws_data.import( data_ds_120 );
      ws_bg.import( bg_ds_120_2017   );
      ws_sig.import( higgs_ds_120_2017  );

      data_ds.append( data_ds_2017 );



    //////////////////////////////////////////////////////////////////////////////
    //////////////////  UNCERTS FOR HIGGS_2017 ////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

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



    }

    bg_ds.append( qcd_ds );   //append bg data into one
    bg_ds.append( gjets_ds ); //append bg data into one

    RooDataSet higgs_ds_130(  higgs_ds, Form("higgs_2016_130_13TeV_%s", regionSaveName.c_str() ) );
    RooDataSet bg_ds_130   (  bg_ds,    Form("bg_2016_130_13TeV_%s",    regionSaveName.c_str() ) );
    RooDataSet data_ds_130 (  data_ds,  Form("data_130_13TeV_%s",  regionSaveName.c_str() ) );

    RooDataSet higgs_ds_125(  higgs_ds, Form("higgs_2016_125_13TeV_%s", regionSaveName.c_str() ) );
    RooDataSet bg_ds_125   (  bg_ds,    Form("bg_2016_125_13TeV_%s",    regionSaveName.c_str() ) );
    RooDataSet data_ds_125 (  data_ds,  Form("data_125_13TeV_%s",  regionSaveName.c_str() ) );

    RooDataSet higgs_ds_120(  higgs_ds, Form("higgs_2016_120_13TeV_%s", regionSaveName.c_str() ) );
    RooDataSet bg_ds_120   (  bg_ds,    Form("bg_2016_120_13TeV_%s",    regionSaveName.c_str() ) );
    RooDataSet data_ds_120 (  data_ds,  Form("data_120_13TeV_%s",  regionSaveName.c_str() ) );

    ((RooRealVar*)higgs_ds_130.addColumn( hgg_mass_130 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_130.addColumn( hgg_mass_130 ))->setRange(100.,180.); 

    ((RooRealVar*)higgs_ds_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
    ((RooRealVar*)data_ds_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 

    ((RooRealVar*)higgs_ds_120.addColumn( hgg_mass_120 ))->setRange(100.,180.);
    ((RooRealVar*)data_ds_120.addColumn( hgg_mass_120 ))->setRange(100.,180.);


    //Import the 2016 versions
    ws_data.import( data_ds_125  );
    ws_bg  .import( bg_ds_125    );
    ws_sig .import( higgs_ds_125 );

    ws_data.import( data_ds_130  );
    ws_bg  .import( bg_ds_130    );
    ws_sig .import( higgs_ds_130 );

    ws_data.import( data_ds_120  );
    ws_bg  .import( bg_ds_120    );
    ws_sig .import( higgs_ds_120 );

    //write number of data events per region into txt file
    dataAmount << regionSaveName << "    " << data_ds_125.sumEntries() << std::endl;
    bgSumAmount << regionSaveName << "    " << 36.* bg_ds_125.sumEntries() << std::endl;

    higgsAmount << regionSaveName << "    " <<  higgs_ds_125.sumEntries() << std::endl;

    //////////////////////////////////////////////////////////////////////////////
    //////////////////  UNCERTS FOR HIGGS ////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

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

	TString sigSampleName(signalList[i]);

	std::string mLSP = signalList[i].substr( signalList[i].find("mLSP")+4);
	std::string::size_type sz;   // alias of size_t

	if( !sigSampleName.Contains("TChi")  )
	  massChi = std::stoi (mLSP,&sz);
	else if( sigSampleName.Contains("TChiH")  )
	  massChi = 1;
	else if( sigSampleName.Contains("TChiWH")  ){ 
	  mLSP = signalList[i].substr( signalList[i].find("_HToGG_")+11);
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
	  mSbottom = signalList[i].substr( signalList[i].find("_HToGG_")+7);
	  std::cout << " names   " << mSbottom << std::endl;
	}

	massSbottom = std::stoi (mSbottom,&sz);
	std::cout << "bottom " << mSbottom << " mass   " << massSbottom <<std::endl;

	int binBottom = h2D_isr->GetXaxis()->FindBin( (float)massSbottom );
	int binChi    = h2D_isr->GetYaxis()->FindBin( (float)massChi );

	h2D_isr->SetBinContent( binBottom, binChi, h_isr->GetMean() );
	h2D_isr_err->SetBinContent( binBottom, binChi, sign* isr_uncert );

	delete h_isr; delete h_isr_UP; delete h_isr_DN;

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


      // 2017 version ///////////////////////////////////////////////////////////////
      if( doComb ) {
	std::cout << "Working on 2017 signal nr " << i << std::endl;

	TTree* tree_sig_ini_2017 = signals_2017[i]->get(thisRegion)->tree;  
	TTree* tree_sig_2017 = tree_sig_ini_2017->CopyTree(treeSel.c_str() );
	RooDataSet sig_ds_2017(Form("%s_2017_125_13TeV_%s", signalList[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_2017) ) ;

	if(doHZ){
	  std::cout << "ADDING HH to HZ" << std::endl;

	  TTree* tree_sig_ini_HH_HH0p25_2017 = signals_HH_HH0p25_2017[i]->get(thisRegion)->tree;  
	  TTree* tree_sig_HH_HH0p25_2017 = tree_sig_ini_HH_HH0p25_2017->CopyTree(treeSel.c_str() );
	  RooDataSet sig_ds_HH_HH0p25_2017(Form("%s_2017_125_13TeV_%s", signalList_HH_HH0p25[i].c_str(),regionSaveName.c_str() ), "Signal", vars, WeightVar(weight), Import(*tree_sig_HH_HH0p25_2017) ) ;

	  sig_ds_2017.append( sig_ds_HH_HH0p25_2017 );
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
	    mLSP = signalList[i].substr( signalList[i].find("_HToGG_")+11);
	    std::cout << "m2 " << mLSP << std::endl;
	    massChi = std::stoi (mLSP,&sz);
	  }

	  std::string mSbottom = signalList[i].substr( signalList[i].find("bottom")+6);
	  if( sigSampleName.Contains("TChiH")  ){
	    mSbottom = signalList[i].substr( signalList[i].find("_m")+2);
	  }
	  if( sigSampleName.Contains("TChiWH")  ){
	    mSbottom = signalList[i].substr( signalList[i].find("_HToGG_")+7);
	  }

	  massSbottom = std::stoi (mSbottom,&sz);
	
	  int binBottom = h2D_isr_2017->GetXaxis()->FindBin( (float)massSbottom );
	  int binChi    = h2D_isr_2017->GetYaxis()->FindBin( (float)massChi );

	  h2D_isr_2017->SetBinContent( binBottom, binChi, h_isr_2017->GetMean() );
	  h2D_isr_2017_err->SetBinContent( binBottom, binChi, sign* isr_uncert );

	  delete h_isr_2017; delete h_isr_2017_UP; delete h_isr_2017_DN;

	  if( massSbottom==450 && ( massChi==1 || massChi==300  ) )
	    sigAmount <<  signalList[i] << "_2017" << " \t " << regionSaveName << " \t " <<  sig_ds.sumEntries()  << std::endl;
	


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



	} 

	if( i==0 )
	  sigAmount <<  signalList[i] << " \t " << regionSaveName << " \t " <<  sig_ds.sumEntries()  << std::endl;

	RooDataSet sig_ds_2017_130(  sig_ds_2017, Form("%s_2017_130_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ) );
	RooDataSet sig_ds_2017_125(  sig_ds_2017, Form("%s_2017_125_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ) );
	RooDataSet sig_ds_2017_120(  sig_ds_2017, Form("%s_2017_120_13TeV_%s", signalList[i].c_str(), regionSaveName.c_str() ) );

	((RooRealVar*)sig_ds_2017_130.addColumn( hgg_mass_130 ))->setRange(100.,180.);  
	((RooRealVar*)sig_ds_2017_125.addColumn( hgg_mass_125 ))->setRange(100.,180.); 
	((RooRealVar*)sig_ds_2017_120.addColumn( hgg_mass_120 ))->setRange(100.,180.); 

	ws_sig.import( sig_ds_2017_125  );
	ws_sig.import( sig_ds_2017_130  );
	ws_sig.import( sig_ds_2017_120  );
      }

    }//done signals for 1 region
   
    h2D_isr->Write();
    h2D_isr_err->Write();

    delete h2D_isr; delete h2D_isr_err;

    h2D_isr_2017->Write();
    h2D_isr_2017_err->Write();

    delete h2D_isr_2017; delete h2D_isr_2017_err;


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
  bgSumAmount.close();
  higgsAmount.close();

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
