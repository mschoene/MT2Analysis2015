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

using namespace std;

#define mt2_cxx
#include "interface/mt2.h"


double lumiErr = 0.026;
bool shapeNorm = false;

void drawShape( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, const std::string& saveName, const std::string& varName, const std::string& selection, const std::string& selection2, const std::string& selection3, unsigned int nSel, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="" );

void drawStacks(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  MT2Analysis<MT2EstimateTree>* data,const MT2Region thisRegion, std::string cut, float lumi);

void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="" );

std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units );


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

  //  std::string outputdir = cfg.getEventYieldDir() + "/zllPurity";

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
  

  //  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string ggDir = cfg.getEventYieldDir() + "/diPhotonControlRegion/";
 

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this


  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_qcd.root", ggDir.c_str()  ), "qcd");
  qcd->setColor(kQCD);
  qcd->setFullName("QCD");

  MT2Analysis<MT2EstimateTree>* diPhoton = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diPhoton.root", ggDir.c_str() ), "diPhoton");
  diPhoton->setColor(kPink-8);
  diPhoton->setFullName("#gamma#gamma");
  
  MT2Analysis<MT2EstimateTree>* gjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_gjets.root", ggDir.c_str() ), "gjets");
  gjets->setColor(kPink-4);
  gjets->setFullName("#gamma+jets");

  MT2Analysis<MT2EstimateTree>* higgs = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir.c_str() ), "higgs");
  higgs->setColor(kCyan-6);
  higgs->setFullName("H #gamma#gamma");




  // MT2Analysis<MT2EstimateTree>* sig_2 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_T2bH_mSbottom450_mLSP300");
  // sig_2->setColor(kBlue-3);
  // sig_2->setFullName("T2bH 450,300");


  // MT2Analysis<MT2EstimateTree>* sig_9 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_TChiHZ_HToGG_m127");
  // sig_9->setColor(kMagenta+1);
  // sig_9->setFullName("#chiHZ 127");


  // MT2Analysis<MT2EstimateTree>* sig_99 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_TChiHH_HToGG_m127");
  // sig_99->setColor(kCyan-4);
  // sig_99->setFullName("#chiHH 127");


  // MT2Analysis<MT2EstimateTree>* sig_11 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_TChiWH_HToGG_127_1");
  // sig_11->setColor(kMagenta+1);
  // sig_11->setFullName("#chiWH 127,1");



  // MT2Analysis<MT2EstimateTree>* sig = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_T2bH_mSbottom450_mLSP1");
  // // //  sig->setColor(kCyan-8);
  // sig->setColor(kCyan+3);
  // sig->setFullName("T2bH 450 1");

  // // MT2Analysis<MT2EstimateTree>* sig_3 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_T2bH_mSbottom250_mLSP1");
  // // sig_3->setColor(kGreen+1);
  // // sig_3->setFullName("T2bH 250 1");



  // // // MT2Analysis<MT2EstimateTree>* sig_10 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_TChiHZ_HToGG_m300");
  // // // sig_10->setColor(kGray+1);
  // // // sig_10->setFullName("#chiHZ 300");

  // // // MT2Analysis<MT2EstimateTree>* sig_11 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_TChiHZ_HToGG_m400");
  // // // sig_11->setColor(kCyan-4);
  // // // sig_11->setFullName("#chiHZ 400");


  // // // MT2Analysis<MT2EstimateTree>* sig_10 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_TChiHZ_HToGG_m200");
  // // // sig_10->setColor(kGray+1);
  // // // sig_10->setFullName("#chiHZ 200");


  // // // MT2Analysis<MT2EstimateTree>* sig_9 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_TChiHH_HToGG_m1000");
  // // // sig_9->setColor(kMagenta+1);
  // // // sig_9->setFullName("#chiHH 1000");

  // // // MT2Analysis<MT2EstimateTree>* sig_10 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_TChiHZ_HToGG_m1000");
  // // // sig_10->setColor(kGray+1);
  // // // sig_10->setFullName("#chiHZ 1000");

  // // // MT2Analysis<MT2EstimateTree>* sig_11 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_sig.root", ggDir.c_str() ), "SMS_TChiWH_HToGG_m1000");
  // // // sig_11->setColor(kCyan-4);
  // // // sig_11->setFullName("#chiWH 1000");






  std::string plotsDir = cfg.getEventYieldDir() + "/diPhotonControlRegion/plotsDataMC";
  
  plotsDir = cfg.getEventYieldDir() + "/diPhotonControlRegion/plotsDataMC_sig/";
  // plotsDir = cfg.getEventYieldDir() + "/diPhotonControlRegion/plotsDataMC_EW_masses/";
  //plotsDir = cfg.getEventYieldDir() + "/diPhotonControlRegion/plotsDataMC_EW/";


  if( shapeNorm ) plotsDir += "_shape";

  bool doSignal = false;
  if( doSignal)
    plotsDir += "_sig";

  //if( sig) plotsDir += "_T2bH_450_1";

  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ggDir.c_str()) , "diPhoton_data");
  data->setFullName("Data 2016");

  MT2Analysis<MT2EstimateTree>* data2 = MT2Analysis<MT2EstimateTree>::readFromFile( "EventYields_data_2017_nov13/diPhotonControlRegion/data.root"  , "diPhoton_data");
  data2->setColor( kData );
  data2->setFullName("Data 2017");


  MT2Analysis<MT2EstimateTree>* higgs_2017 = MT2Analysis<MT2EstimateTree>::readFromFile( "EventYields_data_2017_nov13/diPhotonControlRegion/mc_higgs.root" , "higgs");
  higgs_2017->setColor(kCyan-6);
  higgs_2017->setFullName("H #gamma#gamma 2017");

  
  //This looks and works ok for the fractions but not the stack
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
  //  bgYields.push_back( diPhoton );
  //  bgYields.push_back( gjets );
  //  bgYields.push_back( qcd );
  bgYields.push_back( higgs );

  // std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
  // bgYields.push_back( diPhoton );
  // bgYields.push_back( gjets );
  // bgYields.push_back( qcd );
  // bgYields.push_back( higgs );

  // bgYields.push_back( sig_4 );
  // bgYields.push_back( sig_5 );
  // bgYields.push_back( sig_6 );
  // bgYields.push_back( sig_7 );
  // bgYields.push_back( sig_8 );

  // if( doSignal){ 
  //   // bgYields.push_back( sig );
  //   bgYields.push_back( sig_2 );
  //   // // bgYields.push_back( sig_3 );
  //   bgYields.push_back( sig_9 );
  //   bgYields.push_back( sig_99 );
  //   bgYields.push_back( sig_11 );
  // }
 
  // //For the lovely stack
  // std::vector<MT2Analysis<MT2EstimateTree>* > bgYields2; 
  // bgYields2.push_back( diPhoton );
  // bgYields2.push_back( gjets );
  // bgYields2.push_back( qcd );
  // bgYields2.push_back( higgs );

  // bgYields2.push_back( sig_4 );
  // bgYields2.push_back( sig_10 );
  // bgYields2.push_back( sig_9 );
  // bgYields2.push_back( sig );

  std::vector<MT2Analysis<MT2EstimateTree>* > data2017Yields; 
  data2017Yields.push_back( data2 );

  MT2DrawTools dt(plotsDir, cfg.lumi() );
  dt.set_shapeNorm( shapeNorm );
 
  
  if( !doSignal)
    dt.set_data( data );
  
  //  dt.set_data2( data2 );

  //  dt.set_data2( data2 );
  //dt.set_mcSF2( 42.4/35.96 );
  
  dt.set_mc( &bgYields );
  dt.set_lumi( cfg.lumi() );


  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields2; 
  dt.set_data( higgs );
  bgYields2.push_back( higgs_2017 );
  dt.set_mc( &bgYields2 );
  dt.set_lumi( 1. );

  dt.set_addOverflow( 0 );
  //  dt.set_addOverflow( 1 );


  // +++++++++++++++++++++++++
  // +++      w/ monjet   +++
  // +++++++++++++++++++++++++

  float htMin=0, htMax=-1;
  std::string cutsLabel = getCutLabel(htMin, htMax, "H_{T}", "GeV");
  cutsLabel = " 100<m_{#gamma#gamma}<180";

  std::string selection = "";//"(ht>200. && nJets==1 && met>200 && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && abs(Z_mass-91.19)<20 )";








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
  std::vector<std::string> regionNames;
  std::vector<std::string> regionSelection;

  unsigned int nMT2split = regionNames_ph.size(); // 13 if 1 lep high pt is not split in mt2
  for( unsigned int j=0; j<mt2_bins.size(); j++){ 
    for( unsigned int i=0; i < regionLength; i++){
      regionSelection.push_back( regionSelection_ph[i] );
      regionSelection[ i + j*regionLength ] +=  "  && " +  mt2_sel[j] ;
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



  for( unsigned int i=0; i < regionSelection.size(); i++){
    string treeSel = "((( ptGamma0/h_mass) > 1./3) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 ) && passPrefire" ;
    std::string  selHere = treeSel + " && " + regionSelection[i];
    //    selHere = "weight * (" + selHere + ")";



    dt.drawRegionFractions_fromTree(Form( "H_mass_%s", regionNames[i].c_str() )  , "h_mass"   , selHere, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, regionNames[i] );

    //   dt.drawRegionYields_fromTree(Form( "H_mass_%s", regionNames[i].c_str() )  , "h_mass"   , selHere, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, regionNames[i] );

 
    
  }





  //  ((( ptGamma0/h_mass) > 1./3) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 )


  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && h_mass>=100 && h_mass<=180 && (!isDiBH && !isDiBZ) && passPrefire && (id>950 && id<958)";

  dt.drawRegionFractions_fromTree( "bl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");


  dt.drawRegionFractions_fromTree( "bl_gg_nJets" , "gg_nJets" , selection, 10, -0.5, 9.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionFractions_fromTree( "bl_gg_nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");

  dt.drawRegionFractions_fromTree( "bl_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionFractions_fromTree( "bl_std_mt2"   , "std_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "bl_hgg_met"   , "hgg_met"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionFractions_fromTree( "bl_met"   , "met"   , selection, 50, 0., 200., "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");


  dt.drawRegionYields_fromTree( "bl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");

  dt.drawRegionYields_fromTree( "bl_gg_nJets" , "gg_nJets" , selection, 10, -0.5, 9.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionYields_fromTree( "bl_gg_nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");

  dt.drawRegionYields_fromTree( "bl_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionYields_fromTree( "bl_std_mt2"   , "std_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionYields_fromTree( "bl_met"   , "met"   , selection, 50, 0., 200., "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

  dt.drawRegionYields_fromTree( "bl_std_met"   , "std_met"   , selection, 50, 0., 200., "Standard ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");


  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4442 && fabs(etaGamma1)<=1.4442  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && ((h_mass>136. || h_mass<114.) ) && passPrefire";
  //  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.442 && fabs(etaGamma1)<=1.442  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && ((h_mass>136. || h_mass<114.) || weight!=1) ";
  dt.drawRegionFractions_fromTree( "bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionYields_fromTree( "bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442  && r9Gamma0 > 0.5 && r9Gamma1>0.5   && h_mass>=100 && h_mass<=180 && passPrefire";


  // //Data data comparison

  // ////////////////////////////////////////////////////////////
  // /////////  BASELINE  but data against data   ///////////////////////
  // ////////////////////////////////////////////////////////////

  //  dt.set_data( 0 );


  dt.set_data2( 0 );
  dt.set_mc( & data2017Yields );

  dt.set_lumi( 35.922/41.529 );

  dt.drawRegionFractions_fromTree( "dataComp_bl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionYields_fromTree( "dataComp_bl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");

  dt.drawRegionFractions_fromTree( "dataComp_bl_gg_nJets" , "gg_nJets" , selection, 10, -0.5, 9.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionFractions_fromTree( "dataComp_bl_gg_nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");

  dt.drawRegionFractions_fromTree( "dataComp_bl_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionFractions_fromTree( "dataComp_bl_std_mt2"   , "std_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");

  dt.drawRegionYields_fromTree( "dataComp_bl_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionYields_fromTree( "dataComp_bl_gg_nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");

  dt.drawRegionYields_fromTree( "dataComp_bl_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");

  dt.drawRegionYields_fromTree( "dataComp_bl_met"   , "met"   , selection, 50, 0., 200., "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionFractions_fromTree( "dataComp_bl_met"   , "met"   , selection, 50, 0., 200., "ME_{T} ", "GeV", cutsLabel, "#geq0j, #geq0b");

  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4442 && fabs(etaGamma1)<=1.4442  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && ((h_mass>136. || h_mass<114.) ) && passPrefire";
  //  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.442 && fabs(etaGamma1)<=1.442  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && ((h_mass>136. || h_mass<114.) || weight!=1) ";
  dt.drawRegionYields_fromTree( "dataComp_bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
  dt.drawRegionFractions_fromTree( "dataComp_bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && h_mass>=100 && h_mass<=180&& passPrefire ";

  // dt.set_mc( &data2016Yields );
  // dt.set_lumi( 41.37/35.96 );
  // //dt.set_lumi( 1 ); 

  // // dt.set_displaySF( true );
  // dt.set_displaySF( false );

  // //selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4";

  //selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && ((h_mass>135. || h_mass<115.) || weight!=1) ";
  // dt.drawRegionFractions_fromTree( "sig_bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

  // selection = " ( gg_nJets==0 && ( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5   && ( (!isDiBH && !isDiBZ) && !(is1El) && !(is1Mu) && !(isDiLepZ) )";
  
  // dt.drawRegionFractions_fromTree( "0j_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "0j_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");



  // selection = " ( gg_nJets>0 &&  ( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && gg_nJets>=4  ";

  // dt.drawRegionFractions_fromTree( "sig_4j_0btoInf_mt2"    , "mt2"   , selection, 50, 0., 200., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_0btoInf_gg_met" , "gg_met"   , selection, 50, 0, 500, "ME_{T} without #gamma#gamma", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_0btoInf_met"    , "met"   , selection, 50, 0, 500, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_0btoInf_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T} without #gamma#gamma", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_0btoInf_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
   
  // dt.drawRegionYields_fromTree( "sig_4j_0btoInf_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_0btoInf_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");

  // selection = " ( gg_nJets>0 &&  ( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && gg_nJets>=4 && (h_pt/h_mass)<=0.8 ";

  // dt.drawRegionYields_fromTree( "sig_4j_0btoInf_pt0_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_0btoInf_pt0_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");

  // selection = " ( gg_nJets>0 &&  ( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && gg_nJets>=4 && (h_pt/h_mass)>0.8  ";

  // dt.drawRegionYields_fromTree( "sig_4j_0btoInf_pt1_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_0btoInf_pt1_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");


  // selection = " ( gg_nJets>0 &&  ( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && gg_nJets>=4 && gg_nBJets>=2 && (h_pt/h_mass)>0.8 ";

  // dt.drawRegionYields_fromTree( "sig_4j_2b_pt1_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_2b_pt1_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");

  // selection = " ( gg_nJets>0 &&  ( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && gg_nJets>=4 && gg_nBJets==1 && (h_pt/h_mass)>0.8 ";

  // dt.drawRegionYields_fromTree( "sig_4j_1b_pt1_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_1b_pt1_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");



  // selection = " ( gg_nJets>0 &&  ( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && gg_nJets>=4 && gg_nBJets>=2 && (h_pt/h_mass)<=0.8 ";

  // dt.drawRegionYields_fromTree( "sig_4j_2b_pt0_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_2b_pt0_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");

  // selection = " ( gg_nJets>0 &&  ( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && gg_nJets>=4 && gg_nBJets==1 && (h_pt/h_mass)<=0.8 ";

  // dt.drawRegionYields_fromTree( "sig_4j_1b_pt0_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_4j_1b_pt0_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");









  // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && ((h_mass>136. || h_mass<114.) || weight!=1) ";
  // dt.drawRegionYields_fromTree( "sig_bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");



  // selection = " ( gg_nJets>0 &&  ( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 ";
  // //  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 ";

  // dt.drawRegionYields_fromTree( "sig_bl_mt2"   , "mt2"   , selection, 50, 0., 200., "full M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "sig_bl_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 200., "M_{T2} w/o #gamma", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "sig_bl_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionYields_fromTree( "sig_bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "sig_bl_mt2"   , "mt2"   , selection, 50, 0., 200., "full M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 200., "M_{T2} w/o #gamma", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_hgg_mt2"   , "hgg_mt2"   , selection, 50, 0., 200., "M_{T2} with H", "GeV", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionFractions_fromTree( "sig_bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");





  // dt.drawRegionFractions_fromTree( "sig_is1El" , "is1El" , selection, 2, -0.5, 1.5, "Single Electron", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_is1Mu" , "is1Mu" , selection, 2, -0.5, 1.5, "Single Muon", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_isDiLepZ" , "isDiLepZ" , selection, 2, -0.5, 1.5, "Di Lepton", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_diLepId" , "diLepId" , selection, 4, 10.5, 13.5, "Di Lepton Id", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_diLepMass" , "diLepMass" , selection, 40, 70, 110, "Di Lepton Mass", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_nVert" , "nVert" , selection, 60, 0, 60, "nVert", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "sig_isDiBZ" , "isDiBZ" , selection, 2, -0.5, 1.5, "Di B Z window", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_isDiBH" , "isDiBH" , selection, 2, -0.5, 1.5, "Di B H window", "", cutsLabel, "#geq0j, #geq0b");
  // // dt.drawRegionFractions_fromTree( "sig_bl_scProd_prod"   , "(scProd11*scProd12)>0"   , selection, 2,-.5, 1.5, "Same hemisphere?", "", cutsLabel, "#geq0j, #geq0b");
  // // dt.drawRegionFractions_fromTree( "sig_bl_scProd11"   , "scProd11"   , selection, 100, -5000., 50000., "Same hemisphere? 1", "", cutsLabel, "#geq0j, #geq0b");
  // // dt.drawRegionFractions_fromTree( "sig_bl_scProd12"   , "scProd12"   , selection, 100, -5000., 50000., "Same hemisphere? 2", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gg_deltaPhiMin"    , "gg_deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gg_diffMetMht"    , "gg_diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gg_nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionFractions_fromTree( "sig_bl_gg_nBJetsCSV", "gg_nBJetsCSV", selection, 7, -0.5, 6.5, "Number of CSV b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_deltaPhiMin"    , "deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_diffMetMht"    , "diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_nJets" , "nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "sig_bl_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionFractions_fromTree( "sig_bl_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionFractions_fromTree( "sig_bl_H_eta"   , "h_eta"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "sig_bl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_pt0t_oM"   , "ptGamma0/h_mass"   , selection, 50, 0., 2., "Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_pt1t_oM"   , "ptGamma1/h_mass"   , selection, 50, 0., 2., "Sub-Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");


  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_pt0"   , "ptGamma0"   , selection, 40, 0., 700., "Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_pt1"   , "ptGamma1"   , selection, 40, 0., 400., "Sub-Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_eta0"   , "etaGamma0"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_eta1"   , "etaGamma1"   , selection, 30, -3., 3., "Sub-Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_phi0"   , "phiGamma0"   , selection,30 , 0., 3.3, "Leading #gamma #phi", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_phi1"   , "phiGamma1"   , selection, 30, 0., 3.3, "Sub-Leading #gamma #phi", "", cutsLabel, "#geq0j, #geq0b");

 
  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_r90"   , "r9Gamma0"   , selection, 40, 0., 2., "Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_r91"   , "r9Gamma1"   , selection, 40, 0., 2., "Sub-Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_sigmaIetaIeta0"   , "sigmaIetaIetaGamma0"   , selection, 40, 0., 0.04, "Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "sig_bl_gamma_sigmaIetaIeta1"   , "sigmaIetaIetaGamma1"   , selection, 40, 0., 0.04, "Sub-Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");

  // //  dt.drawRegionFractions_fromTree( "sig_bl_gamma_chHadIso0"   , "chHadIsoGamma0"   , selection, 40, 0., 2.5, "Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionFractions_fromTree( "sig_bl_gamma_chHadIso1"   , "chHadIsoGamma1"   , selection, 40, 0., 2.5, "Sub-Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");


  // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && nVert>15 && nVert<=25 ";
  // dt.drawRegionFractions_fromTree( "sig_bl_met_nVert"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");









  // std::string selection = "(ht>200. && nJets==1 && met>200 && mt2>200. && deltaPhiMin>0.3 && diffMetMht<0.5*met && abs(Z_mass-91.19)<10 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)) )";



  // string selection1 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 &&( r9Gamma0<0.8 && r9Gamma1<0.8 ) "; //Very bad
  // string selection2 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 &&( r9Gamma0<0.8 || r9Gamma1<0.8 ) && !(r9Gamma0<0.8 && r9Gamma1<0.8) && !(r9Gamma0>0.8 && r9Gamma1>0.8)  "; // semi bad
  // string selection3 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && r9Gamma0>0.8 && r9Gamma1>0.8  "; //good

//   string selection1 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 &&( r9Gamma0<0.8 && r9Gamma1<0.8 ) "; //Very bad
//   string selection2 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 &&( r9Gamma0<0.8 || r9Gamma1<0.8 ) && !(r9Gamma0<0.8 && r9Gamma1<0.8) && !(r9Gamma0>0.8 && r9Gamma1>0.8)  "; // semi bad
//   string selection3 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0>0.8 && r9Gamma1>0.8  "; //good
//   drawShape(  cfg, sig, "incl_r9at0p8_mgg", "h_mass", selection1, selection2, selection3, 3,  50, 110., 140., "M_{#gamma#gamma}", "GeV" );

//   // selection1 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 &&( r9Gamma0<0.9 || r9Gamma1<0.9 ) ";
//   // selection2 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0>0.9 && r9Gamma1>0.9  ";
//   // drawShape(  cfg, sig, "incl_r9at0p9_mgg", "h_mass", selection1, selection2, 2,  50, 110., 140., "M_{#gamma#gamma}", "GeV" );


//   selection1 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 &&( r9Gamma0<0.94 && r9Gamma1<0.94 ) "; //Very bad
//   selection2 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 &&( r9Gamma0<0.94 || r9Gamma1<0.94 ) && !(r9Gamma0<0.94 && r9Gamma1<0.94) && !(r9Gamma0>0.94 && r9Gamma1>0.94)  "; // semi bad
//   selection3 = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0>0.94 && r9Gamma1>0.94  "; //good

//   drawShape(  cfg, sig, "incl_r9at0p94_mgg", "h_mass", selection3, selection2, selection1, 3,  50, 110., 140., "M_{#gamma#gamma}", "GeV" );



//   selection1 = " fabs(etaGamma0)>1.4  && fabs(etaGamma1)>1.4  "; //Very bad
//   selection2 = "( fabs(etaGamma0)<=1.4 || fabs(etaGamma1)<=1.4 ) && !(fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 ) && !(fabs(etaGamma0)>1.4  && fabs(etaGamma1)>1.4 )"; // semi bad
//   selection3 = " fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  "; //good

//   drawShape(  cfg, sig, "incl_etaCut_mgg", "h_mass", selection3, selection2, selection1, 3,  50, 110., 140., "M_{#gamma#gamma}", "GeV" );

// //void drawShape( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, const std::string& saveName, const std::string& varName, const std::string& selection, const std::string& selection2, int nSel, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="" );





  // selection = "((( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180  && ( isDiBH || isDiBZ ) )" ;
  // //    treeSel += " && (!isDiBH && !isDiBZ) ";


  // dt.drawRegionFractions_fromTree( "diB_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  
  // dt.drawRegionFractions_fromTree( "diB_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");




// ///////////////TESTING///////////////////


//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && (( (scProd11*scProd12)>0) || gg_nJets<1 ) ";


//   dt.drawRegionFractions_fromTree( "sameHemi_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "sameHemi_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "sameHemi_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   //  dt.drawRegionFractions_fromTree( "sameHemi_deltaPhiMin"    , "deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionFractions_fromTree( "sameHemi_diffMetMht"    , "diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "sameHemi_nJets" , "nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "sameHemi_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "sameHemi_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionFractions_fromTree( "sameHemi_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");



//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.8 && r9Gamma1>0.8 &&  (sqrt(2*h_pt*met*(1-(cos(met_phi-h_phi)))) )>80 ";


//   // dt.drawRegionFractions_fromTree( "highMT_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionFractions_fromTree( "highMT_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionFractions_fromTree( "highMT_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

//   // dt.drawRegionFractions_fromTree( "highMT_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionFractions_fromTree( "highMT_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionFractions_fromTree( "highMT_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionFractions_fromTree( "highMT_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionFractions_fromTree( "highMT_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");

//   // dt.drawRegionYields_fromTree( "highMT_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");




// // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 ";

// //  dt.drawRegionFractions_fromTree( "bl_mThiggs"   , "sqrt(2*h_pt*met*(1-(cos(met_phi-h_phi))))"   , selection, 50, 0., 200, "M_{T}(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");

// //  dt.drawRegionFractions_fromTree( "bl_gg_mThiggs"   , "sqrt(2*h_pt*met*(1-(cos(gg_met_phi-h_phi))))"   , selection, 50, 0., 200, "M_{T}(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");



// // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && met>50";

// //  dt.drawRegionFractions_fromTree( "bl_test_mThiggs"   , "sqrt(2*h_pt*met*(1-(cos(met_phi-h_phi))))"   , selection, 50, 0., 200, "M_{T}(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");

// //  dt.drawRegionFractions_fromTree( "bl_test_gg_mThiggs"   , "sqrt(2*h_pt*met*(1-(cos(gg_met_phi-h_phi))))"   , selection, 50, 0., 200, "M_{T}(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");




// // ///////////////TESTING///////////////////
// //   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.8 && r9Gamma1>0.8 &&  ( h_pt/h_mass )>0.8 && met >100 ";

// //   dt.drawRegionYields_fromTree( "highMet_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

// //   dt.drawRegionFractions_fromTree( "highMet_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highMet_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highMet_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

// //   dt.drawRegionFractions_fromTree( "highMet_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highMet_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highMet_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highMet_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highMet_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");

// //   dt.drawRegionFractions_fromTree( "highMet_deltaPhiMetGG_gg"   , "acos(cos(gg_met_phi-h_phi))"   , selection, 50, 0., 3.2, "#Delta#Phi(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highMet_deltaPhiMetGG"   , "acos(cos(met_phi-h_phi))"   , selection, 50, 0., 3.2, "#Delta#Phi(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");






// // ///////////////TESTING///////////////////
// //   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.8 && r9Gamma1>0.8 &&  ( h_pt/h_mass )>0.8 ";


// //   dt.drawRegionFractions_fromTree( "highPt_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highPt_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highPt_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

// //   dt.drawRegionFractions_fromTree( "highPt_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highPt_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highPt_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highPt_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highPt_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");

// //   dt.drawRegionFractions_fromTree( "highPt_deltaPhiMetGG_gg"   , "acos(cos(gg_met_phi-h_phi))"   , selection, 50, 0., 3.2, "#Delta#Phi(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highPt_deltaPhiMetGG"   , "acos(cos(met_phi-h_phi))"   , selection, 50, 0., 3.2, "#Delta#Phi(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");




// // ///////////////TESTING///////////////////
// //   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.8 && r9Gamma1>0.8 ";


// //   dt.drawRegionFractions_fromTree( "highr9_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highr9_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highr9_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

// //   dt.drawRegionFractions_fromTree( "highr9_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highr9_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highr9_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highr9_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionFractions_fromTree( "highr9_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");










// //   // +++++++++++++++++++++++++
// //   // +++     inclusive     +++
// //   // +++++++++++++++++++++++++


//   // selection = "  ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 ))";
//   selection = "";
 
//   dt.drawRegionFractions_fromTree( "incl_gamma_pt0t_oM_noPtCuts"   , "ptGamma0/h_mass"   , selection, 50, 0., 2., "Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gamma_pt1t_oM_noPtCuts"   , "ptGamma1/h_mass"   , selection, 50, 0., 2., "Sub-Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_etaFrac_noPtCuts"   , "(fabs(etaGamma0)>1.4) + (fabs(etaGamma1)>1.4)"   , selection, 3, -0.5, 2.5, "if(isEE) + if(isEE)", "", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_r9Frac0p85_noPtCuts" , "((r9Gamma0)<0.85) + ((r9Gamma1)<0.85)"   , selection, 3, -0.5, 2.5, "if(r9<0.85) + if(r9<0.85)", "", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_htMinus_noPtCuts"    , "ht-ptGamma0-ptGamma1"    , selection, 50, 0., 2000., "H_{T}-#gamma(1) p_{T} -#gamma(2) p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");



//   dt.drawRegionFractions_fromTree( "incl_gamma_r90_noPtCuts"   , "r9Gamma0"   , selection, 40, 0., 2., "Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gamma_r91_noPtCuts"   , "r9Gamma1"   , selection, 40, 0., 2., "Sub-Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");



//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && met>50 ";
//   //dt.drawRegionFractions_fromTree( "incl_angleGammas"   , "angleGammas"   , selection, 50, 0., 3.2, "Angle btwn photons", "", cutsLabel, "#geq0j, #geq0b");
//   //  dt.drawRegionFractions_fromTree( "incl_deltaPhiMetGG"   , "acos(cos(met_phi-h_phi))"   , selection, 50, 0., 3.2, "#Delta#Phi(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");



//   //selection = "  ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) && ptGamma0/h_mass > 0.33 && ptGamma1/h_mass > 0.25 ";
  

//   //for the fraction in eta range
//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) ";
//   dt.drawRegionFractions_fromTree( "incl_etaFrac" , "(fabs(etaGamma0)>1.4) + (fabs(etaGamma1)>1.4)"   , selection, 3, -0.5, 2.5, "if(isEE) + if(isEE)", "", cutsLabel, "#geq0j, #geq0b");



//   ////////////////////////////////////////////////////////////
//   /////////  INCL (ONLY PT and eta cut ///////////////////////
//   ////////////////////////////////////////////////////////////
 
//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_nJets==0 ";
//   dt.drawRegionFractions_fromTree( "incl_gg_deltaPhiMin_set0j" , "gg_deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "0j, #geq0b");

 
//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_nJets==1 ";
//   dt.drawRegionFractions_fromTree( "incl_gg_deltaPhiMin_1j" , "gg_deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "1j, #geq0b");

//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_nJets>=2";
//   dt.drawRegionFractions_fromTree( "incl_gg_deltaPhiMin_gt2j" , "gg_deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq2j, #geq0b");


//   //BLIND the signal region
//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && ((h_mass>135. || h_mass<115.)|| weight!=1 )";
//   dt.drawRegionFractions_fromTree( "incl_H_mass"   , "h_mass"   , selection, 50, 100., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");


//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 ";
//   dt.drawRegionFractions_fromTree( "incl_scProd_prod"   , "(scProd11*scProd12)>0"   , selection, 2,-.5, 1.5, "Same hemisphere?", "", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_scProd11"   , "scProd11"   , selection, 100, -5000., 50000., "Same hemisphere? 1", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_scProd12"   , "scProd12"   , selection, 100, -5000., 50000., "Same hemisphere? 2", "", cutsLabel, "#geq0j, #geq0b");



//   dt.drawRegionFractions_fromTree( "incl_htMinus"    , "ht-ptGamma0-ptGamma1"    , selection, 50, 0., 2000., "H_{T}-#gamma(1) p_{T} -#gamma(2) p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
    
//   dt.drawRegionFractions_fromTree( "incl_r9Frac0p8" , "((r9Gamma0)<0.8) + ((r9Gamma1)<0.8)"   , selection, 3, -0.5, 2.5, "if(r9<0.8) + if(r9<0.8)", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_r9Frac0p95" , "((r9Gamma0)<0.95) + ((r9Gamma1)<0.95)"   , selection, 3, -0.5, 2.5, "if(r9<0.95) + if(r9<0.95)", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_r9Frac0p85" , "((r9Gamma0)<0.85) + ((r9Gamma1)<0.85)"   , selection, 3, -0.5, 2.5, "if(r9<0.85) + if(r9<0.85)", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_r9Frac0p9" , "((r9Gamma0)<0.9) + ((r9Gamma1)<0.9)"   , selection, 3, -0.5, 2.5, "if(r9<0.9) + if(r9<0.9 )", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_r9Frac0p94" , "((r9Gamma0)<0.94) + ((r9Gamma1)<0.94)"   , selection, 3, -0.5, 2.5, "if(r9<0.94) + if(r9<0.94 )", "", cutsLabel, "#geq0j, #geq0b");


//   dt.drawRegionFractions_fromTree( "incl_deltaPhiMetGG_gg"   , "acos(cos(gg_met_phi-h_phi))"   , selection, 50, 0., 3.2, "#Delta#Phi(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_deltaPhiMetGG"   , "acos(cos(met_phi-h_phi))"   , selection, 50, 0., 3.2, "#Delta#Phi(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gg_deltaPhiMin" , "gg_deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gg_diffMetMht"  , "gg_diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gg_nJets" ,"gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)","", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gg_nBJets","gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)","", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_deltaPhiMin"    , "deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_diffMetMht"    , "diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_nJets" , "nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
 
//   dt.drawRegionFractions_fromTree( "incl_H_eta"   , "h_eta"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_gamma_pt0", "ptGamma0"   , selection, 40, 0., 700., "Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gamma_pt1", "ptGamma1"   , selection, 40, 0., 400., "Sub-Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gamma_pt0t_oM"   , "ptGamma0/h_mass"   , selection, 50, 0., 2., "Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gamma_pt1t_oM"   , "ptGamma1/h_mass"   , selection, 50, 0., 2., "Sub-Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");


//   dt.drawRegionFractions_fromTree( "incl_gamma_eta0"   , "etaGamma0"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gamma_eta1"   , "etaGamma1"   , selection, 30, -3., 3., "Sub-Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_gamma_r90"   , "r9Gamma0"   , selection, 40, 0., 2., "Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gamma_r91"   , "r9Gamma1"   , selection, 40, 0., 2., "Sub-Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_gamma_sigmaIetaIeta0"   , "sigmaIetaIetaGamma0"   , selection, 40, 0., 0.04, "Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gamma_sigmaIetaIeta1"   , "sigmaIetaIetaGamma1"   , selection, 40, 0., 0.04, "Sub-Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionFractions_fromTree( "incl_gamma_chHadIso0"   , "chHadIsoGamma0"   , selection, 40, 0., 2.5, "Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionFractions_fromTree( "incl_gamma_chHadIso1"   , "chHadIsoGamma1"   , selection, 40, 0., 2.5, "Sub-Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
//   ////////////////////////////////////////////////////////////
//   ////////////////////////////////////////////////////////////
//   ////////////////////////////////////////////////////////////





















 //  ////////////////////////////////////////////////////////////
 //  ///////// BASELINE                   ///////////////////////
 //  ////////////////////////////////////////////////////////////
 //  //  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 ";
 //  // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 ";

 //  //BLINDING
 //  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && ((h_mass>135. || h_mass<115.)|| weight!=1 )";
 //  dt.drawRegionFractions_fromTree( "bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

 // //SELECTION
 // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0 > 0.5 && r9Gamma1>0.5 ";

 //  dt.drawRegionFractions_fromTree( "bl_scProd_prod"   , "(scProd11*scProd12)>0"   , selection, 2,-.5, 1.5, "Same hemisphere?", "", cutsLabel, "#geq0j, #geq0b");

 //  dt.drawRegionFractions_fromTree( "bl_scProd11"   , "scProd11"   , selection, 100, -5000., 50000., "Same hemisphere? 1", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_scProd12"   , "scProd12"   , selection, 100, -5000., 50000., "Same hemisphere? 2", "", cutsLabel, "#geq0j, #geq0b");


  
 //  dt.drawRegionFractions_fromTree( "bl_htMinus"    , "ht-ptGamma0-ptGamma1"    , selection, 50, 0., 2000., "H_{T}-#gamma(1) p_{T} -#gamma(2) p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
    
 //  dt.drawRegionFractions_fromTree( "bl_r9Frac0p8" , "((r9Gamma0)<0.8) + ((r9Gamma1)<0.8)"   , selection, 3, -0.5, 2.5, "if(r9<0.8) + if(r9<0.8)", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_r9Frac0p95" , "((r9Gamma0)<0.95) + ((r9Gamma1)<0.95)"   , selection, 3, -0.5, 2.5, "if(r9<0.95) + if(r9<0.95)", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_r9Frac0p85" , "((r9Gamma0)<0.85) + ((r9Gamma1)<0.85)"   , selection, 3, -0.5, 2.5, "if(r9<0.85) + if(r9<0.85)", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_r9Frac0p9" , "((r9Gamma0)<0.9) + ((r9Gamma1)<0.9)"   , selection, 3, -0.5, 2.5, "if(r9<0.9) + if(r9<0.9 )", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_r9Frac0p94" , "((r9Gamma0)<0.94) + ((r9Gamma1)<0.94)"   , selection, 3, -0.5, 2.5, "if(r9<0.94) + if(r9<0.94 )", "", cutsLabel, "#geq0j, #geq0b");


 //  dt.drawRegionFractions_fromTree( "bl_deltaPhiMetGG_gg"   , "acos(cos(gg_met_phi-h_phi))"   , selection, 50, 0., 3.2, "#Delta#Phi(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_deltaPhiMetGG"   , "acos(cos(met_phi-h_phi))"   , selection, 50, 0., 3.2, "#Delta#Phi(ME_{T}, H)", "GeV", cutsLabel, "#geq0j, #geq0b");

 //  dt.drawRegionFractions_fromTree( "bl_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gg_deltaPhiMin" , "gg_deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gg_diffMetMht"  , "gg_diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gg_nJets" ,"gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)","", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gg_nBJets","gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)","", cutsLabel, "#geq0j, #geq0b");

  

 //  dt.drawRegionFractions_fromTree( "bl_gg_jet1_pt"   , "gg_jet1_pt"   , selection, 50, -100., 1000., "Pseudo jet1 p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gg_jet2_pt"   , "gg_jet2_pt"   , selection, 50, -100., 600., "Pseudo jet2 p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");


 //  dt.drawRegionFractions_fromTree( "bl_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_deltaPhiMin"    , "deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_diffMetMht"    , "diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_nJets" , "nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  //  dt.drawRegionFractions_fromTree( "bl_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
 
 //  dt.drawRegionFractions_fromTree( "bl_H_eta"   , "h_eta"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");


 //  dt.drawRegionFractions_fromTree( "bl_gamma_pt0", "ptGamma0"   , selection, 40, 0., 700., "Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gamma_pt1", "ptGamma1"   , selection, 40, 0., 400., "Sub-Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

 //  dt.drawRegionFractions_fromTree( "bl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gamma_pt0t_oM"   , "ptGamma0/h_mass"   , selection, 50, 0., 2., "Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gamma_pt1t_oM"   , "ptGamma1/h_mass"   , selection, 50, 0., 2., "Sub-Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");


 //  dt.drawRegionFractions_fromTree( "bl_gamma_eta0"   , "etaGamma0"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gamma_eta1"   , "etaGamma1"   , selection, 30, -3., 3., "Sub-Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

 //  dt.drawRegionFractions_fromTree( "bl_gamma_r90"   , "r9Gamma0"   , selection, 40, 0., 2., "Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gamma_r91"   , "r9Gamma1"   , selection, 40, 0., 2., "Sub-Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");

 //  dt.drawRegionFractions_fromTree( "bl_gamma_sigmaIetaIeta0"   , "sigmaIetaIetaGamma0"   , selection, 40, 0., 0.04, "Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gamma_sigmaIetaIeta1"   , "sigmaIetaIetaGamma1"   , selection, 40, 0., 0.04, "Sub-Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");

 //  dt.drawRegionFractions_fromTree( "bl_gamma_chHadIso0"   , "chHadIsoGamma0"   , selection, 40, 0., 2.5, "Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  dt.drawRegionFractions_fromTree( "bl_gamma_chHadIso1"   , "chHadIsoGamma1"   , selection, 40, 0., 2.5, "Sub-Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
 //  ////////////////////////////////////////////////////////////
 //  ////////////////////////////////////////////////////////////
 //  ////////////////////////////////////////////////////////////









//   ////////////////////////////////////////////////////////////
//   /////////  INCL (ONLY PT and eta cut ///////////////////////
//   ////////////////////////////////////////////////////////////
//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && ((h_mass>135. || h_mass<115.)|| weight!=1 )";
//   dt.set_mc( &bgYields2 );


//   dt.drawRegionYields_fromTree( "incl_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b"); 


//   selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 ";

//   dt.drawRegionYields_fromTree( "incl_scProd_prod"   , "(scProd11*scProd12)>0"   , selection, 2,-0.5, 1.5, "Same hemisphere?", "", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionYields_fromTree( "incl_scProd11"   , "scProd11"   , selection, 100, -5000., 50000., "Same hemisphere? 1", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_scProd12"   , "scProd12"   , selection, 100, -5000., 50000., "Same hemisphere? 2", "", cutsLabel, "#geq0j, #geq0b");


//   dt.drawRegionYields_fromTree( "incl_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gg_deltaPhiMin"    , "gg_deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gg_diffMetMht"    , "gg_diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gg_nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");



//   dt.drawRegionYields_fromTree( "incl_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_deltaPhiMin"    , "deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_diffMetMht"    , "diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_nJets" , "nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");


//   dt.drawRegionYields_fromTree( "incl_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionYields_fromTree( "incl_H_eta"   , "h_eta"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");


//   dt.drawRegionYields_fromTree( "incl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gamma_pt0t_oM"   , "ptGamma0/h_mass"   , selection, 50, 0., 2., "Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gamma_pt1t_oM"   , "ptGamma1/h_mass"   , selection, 50, 0., 2., "Sub-Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");


//   dt.drawRegionYields_fromTree( "incl_gamma_pt0"   , "ptGamma0"   , selection, 40, 0., 700., "Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gamma_pt1"   , "ptGamma1"   , selection, 40, 0., 400., "Sub-Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gamma_eta0"   , "etaGamma0"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gamma_eta1"   , "etaGamma1"   , selection, 30, -3., 3., "Sub-Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

 
//   dt.drawRegionYields_fromTree( "incl_gamma_r90"   , "r9Gamma0"   , selection, 40, 0., 2., "Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gamma_r91"   , "r9Gamma1"   , selection, 40, 0., 2., "Sub-Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionYields_fromTree( "incl_gamma_sigmaIetaIeta0"   , "sigmaIetaIetaGamma0"   , selection, 40, 0., 0.04, "Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gamma_sigmaIetaIeta1"   , "sigmaIetaIetaGamma1"   , selection, 40, 0., 0.04, "Sub-Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");

//   dt.drawRegionYields_fromTree( "incl_gamma_chHadIso0"   , "chHadIsoGamma0"   , selection, 40, 0., 2.5, "Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
//   dt.drawRegionYields_fromTree( "incl_gamma_chHadIso1"   , "chHadIsoGamma1"   , selection, 40, 0., 2.5, "Sub-Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");










// // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.8 && r9Gamma1>0.8 ";

// //  dt.drawRegionYields_fromTree( "blR9at0p8_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   //  dt.drawRegionYields_fromTree( "bl_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionYields_fromTree( "blR9at0p8_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");



// // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && fabs(h_eta)< 1.4 ";

// //  dt.drawRegionYields_fromTree( "blHeta_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   //  dt.drawRegionYields_fromTree( "bl_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
// //   dt.drawRegionYields_fromTree( "blHeta_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");




  // ////////////////////////////////////////////////////////////
  // /////////  BASELINE                  ///////////////////////
  // ////////////////////////////////////////////////////////////
  // //selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4";

  // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && ((h_mass>135. || h_mass<115.) || weight!=1) ";
  // dt.drawRegionYields_fromTree( "bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

  // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 ";



  // // dt.drawRegionYields_fromTree( "bl_scProd_prod"   , "(scProd11*scProd12)>0"   , selection, 2,-.5, 1.5, "Same hemisphere?", "", cutsLabel, "#geq0j, #geq0b");

  // // dt.drawRegionYields_fromTree( "bl_scProd11"   , "scProd11"   , selection, 100, -5000., 50000., "Same hemisphere? 1", "", cutsLabel, "#geq0j, #geq0b");
  // // dt.drawRegionYields_fromTree( "bl_scProd12"   , "scProd12"   , selection, 100, -5000., 50000., "Same hemisphere? 2", "", cutsLabel, "#geq0j, #geq0b");


  // dt.drawRegionYields_fromTree( "bl_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gg_deltaPhiMin"    , "gg_deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gg_diffMetMht"    , "gg_diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gg_nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionYields_fromTree( "bl_gg_nBJetsDeepCSV", "gg_nBJetsDeepCSV", selection, 7, -0.5, 6.5, "Number of DeepCSV b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");


  // dt.drawRegionYields_fromTree( "bl_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_deltaPhiMin"    , "deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_diffMetMht"    , "diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_nJets" , "nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");


  // dt.drawRegionYields_fromTree( "bl_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionYields_fromTree( "bl_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionYields_fromTree( "bl_H_eta"   , "h_eta"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");


  // dt.drawRegionYields_fromTree( "bl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gamma_pt0t_oM"   , "ptGamma0/h_mass"   , selection, 50, 0., 2., "Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gamma_pt1t_oM"   , "ptGamma1/h_mass"   , selection, 50, 0., 2., "Sub-Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");


  // dt.drawRegionYields_fromTree( "bl_gamma_pt0"   , "ptGamma0"   , selection, 40, 0., 700., "Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gamma_pt1"   , "ptGamma1"   , selection, 40, 0., 400., "Sub-Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gamma_eta0"   , "etaGamma0"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gamma_eta1"   , "etaGamma1"   , selection, 30, -3., 3., "Sub-Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

 
  // dt.drawRegionYields_fromTree( "bl_gamma_r90"   , "r9Gamma0"   , selection, 40, 0., 2., "Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gamma_r91"   , "r9Gamma1"   , selection, 40, 0., 2., "Sub-Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionYields_fromTree( "bl_gamma_sigmaIetaIeta0"   , "sigmaIetaIetaGamma0"   , selection, 40, 0., 0.04, "Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gamma_sigmaIetaIeta1"   , "sigmaIetaIetaGamma1"   , selection, 40, 0., 0.04, "Sub-Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionYields_fromTree( "bl_gamma_chHadIso0"   , "chHadIsoGamma0"   , selection, 40, 0., 2.5, "Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "bl_gamma_chHadIso1"   , "chHadIsoGamma1"   , selection, 40, 0., 2.5, "Sub-Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
  // ////////////////////////////////////////////////
  // ////////////////////////////////////////////////
  // ////////////////////////////////////////////////












  // //Data data comparison

  // ////////////////////////////////////////////////////////////
  // /////////  BASELINE  but data against data                ///////////////////////
  // ////////////////////////////////////////////////////////////

  // dt.set_data2( 0 );
  // dt.set_mc( &data2016Yields );
  // //  dt.set_mcSF( 42.4/35.96 );
  // dt.set_lumi( 1 );

  // // dt.set_displaySF( true );
  // dt.set_displaySF( false );

  // //selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4";

  // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5  && ((h_mass>135. || h_mass<115.) || weight!=1) ";
  // dt.drawRegionFractions_fromTree( "dd_bl_H_mass"   , "h_mass"   , selection, 50, 100., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

  // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 ";


  // dt.drawRegionFractions_fromTree( "dd_is1El" , "is1El" , selection, 2, -0.5, 1.5, "Single Electron", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_is1Mu" , "is1Mu" , selection, 2, -0.5, 1.5, "Single Muon", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_isDiLepZ" , "isDiLepZ" , selection, 2, -0.5, 1.5, "Di Lepton", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_diLepId" , "diLepId" , selection, 4, 10.5, 13.5, "Di Lepton Id", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_diLepMass" , "diLepMass" , selection, 40, 70, 110, "Di Lepton Mass", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_nVert" , "nVert" , selection, 60, 0, 60, "nVert", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "dd_isDiBZ" , "isDiBZ" , selection, 2, -0.5, 1.5, "Di B Z window", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_isDiBH" , "isDiBH" , selection, 2, -0.5, 1.5, "Di B H window", "", cutsLabel, "#geq0j, #geq0b");
  // // dt.drawRegionFractions_fromTree( "dd_bl_scProd_prod"   , "(scProd11*scProd12)>0"   , selection, 2,-.5, 1.5, "Same hemisphere?", "", cutsLabel, "#geq0j, #geq0b");
  // // dt.drawRegionFractions_fromTree( "dd_bl_scProd11"   , "scProd11"   , selection, 100, -5000., 50000., "Same hemisphere? 1", "", cutsLabel, "#geq0j, #geq0b");
  // // dt.drawRegionFractions_fromTree( "dd_bl_scProd12"   , "scProd12"   , selection, 100, -5000., 50000., "Same hemisphere? 2", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gg_deltaPhiMin"    , "gg_deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gg_diffMetMht"    , "gg_diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gg_nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionFractions_fromTree( "dd_bl_gg_nBJetsCSV", "gg_nBJetsCSV", selection, 7, -0.5, 6.5, "Number of CSV b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_mt2"   , "mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_met"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_deltaPhiMin"    , "deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_diffMetMht"    , "diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_nJets" , "nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "dd_bl_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionFractions_fromTree( "dd_bl_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionFractions_fromTree( "dd_bl_H_eta"   , "h_eta"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "dd_bl_H_pt_oM"   , "h_pt/h_mass"   , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_pt0t_oM"   , "ptGamma0/h_mass"   , selection, 50, 0., 2., "Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_pt1t_oM"   , "ptGamma1/h_mass"   , selection, 50, 0., 2., "Sub-Leading #gamma p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");


  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_pt0"   , "ptGamma0"   , selection, 40, 0., 700., "Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_pt1"   , "ptGamma1"   , selection, 40, 0., 400., "Sub-Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_eta0"   , "etaGamma0"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_eta1"   , "etaGamma1"   , selection, 30, -3., 3., "Sub-Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_phi0"   , "phiGamma0"   , selection,30 , 0., 3.3, "Leading #gamma #phi", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_phi1"   , "phiGamma1"   , selection, 30, 0., 3.3, "Sub-Leading #gamma #phi", "", cutsLabel, "#geq0j, #geq0b");

 
  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_r90"   , "r9Gamma0"   , selection, 40, 0., 2., "Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_r91"   , "r9Gamma1"   , selection, 40, 0., 2., "Sub-Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_sigmaIetaIeta0"   , "sigmaIetaIetaGamma0"   , selection, 40, 0., 0.04, "Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionFractions_fromTree( "dd_bl_gamma_sigmaIetaIeta1"   , "sigmaIetaIetaGamma1"   , selection, 40, 0., 0.04, "Sub-Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");

  // //  dt.drawRegionFractions_fromTree( "dd_bl_gamma_chHadIso0"   , "chHadIsoGamma0"   , selection, 40, 0., 2.5, "Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
  // //  dt.drawRegionFractions_fromTree( "dd_bl_gamma_chHadIso1"   , "chHadIsoGamma1"   , selection, 40, 0., 2.5, "Sub-Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");


  // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && r9Gamma0 > 0.5 && r9Gamma1>0.5 && nVert>15 && nVert<=25 ";
  // dt.drawRegionFractions_fromTree( "dd_bl_met_nVert"   , "met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");


  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////







  // /// Let's see some distributions for different selections

  // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2. ";
  // dt.drawRegionYields_fromTree( "ggDPhi2_H_mass"   , "h_mass"   , selection, 50, 100., 200., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "ggDPhi2_H_pt_oM"   , "h_pt/h_mass" , selection, 50, 0., 2., "H p_{T}/M_{#gamma#gamma}", "", cutsLabel, "#geq0j, #geq0b");

  // dt.drawRegionYields_fromTree( "ggDPhi2_gg_mt2"   , "gg_mt2"   , selection, 50, 0., 300., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "ggDPhi2_gg_met"   , "gg_met"   , selection, 50, 0, 400, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "ggDPhi2_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "ggDPhi2_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
  // dt.drawRegionYields_fromTree( "ggDPhi2_gg_nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");






//   // //Nominal/
//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 ";


//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && gg_deltaPhiMin<2.5 && r9Gamma0>0.94 && r9Gamma1>0.94  ";
//   // dt.drawRegionYields_fromTree( "highRes_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4  && gg_deltaPhiMin<2.5 && (r9Gamma0<0.94 || r9Gamma1<0.94 ) ";
//   // dt.drawRegionYields_fromTree( "lowRes_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");


//   // //  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && r9Gamma0>0.9 && r9Gamma1>0.9  ";
//   // // dt.drawRegionYields_fromTree( "highRes_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && (r9Gamma0<0.9 || r9Gamma1<0.9 ) ";
//   // // dt.drawRegionYields_fromTree( "lowRes_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");


//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5  && gg_ht>200  ";
//   // dt.drawRegionYields_fromTree( "highggHT_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "highggHT_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "highggHT_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5  && gg_ht<200 ";
//   // dt.drawRegionYields_fromTree( "lowggHT_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "lowggHT_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "lowggHT_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");


//   // // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5  && ht>200  ";
//   // // dt.drawRegionYields_fromTree( "highHT_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // // dt.drawRegionYields_fromTree( "highHT_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // // dt.drawRegionYields_fromTree( "highHT_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   // // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5  && ht<200 ";
//   // // dt.drawRegionYields_fromTree( "lowHT_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // // dt.drawRegionYields_fromTree( "lowHT_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // // dt.drawRegionYields_fromTree( "lowHT_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");


//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && gg_nJets>=4  ";
//   // dt.drawRegionYields_fromTree( "highNJets_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "highNJets_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "highNJets_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && gg_nJets<4  ";
//   // dt.drawRegionYields_fromTree( "lowNJets_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "lowNJets_gg_ht"    , "gg_ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "lowNJets_gg_nJets" , "gg_nJets" , selection, 14, -0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");




//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && ( h_pt/h_mass) >= 0.8 ";
//   // dt.drawRegionYields_fromTree( "highHpt_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && gg_nBJets<2 && ( h_pt/h_mass) < 0.8 ";
//   // dt.drawRegionYields_fromTree( "lowHpt_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");



//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && gg_nBJets>=2  ";
//   // dt.drawRegionYields_fromTree( "highNBJets_H_mass", "h_mass" , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && gg_nBJets<2  ";
//   // dt.drawRegionYields_fromTree( "lowNBJets_H_mass", "h_mass" , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");



//   // //selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && gg_nBJets>=2 && r9Gamma0>0.94 && r9Gamma1>0.94 && gg_ht>200  && gg_nJets>=4 ";
//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && gg_nBJets>=2 && r9Gamma0>0.94 && r9Gamma1>0.94 && gg_nJets>=4 && ( h_pt/h_mass) >= 0.8";
//   // dt.drawRegionYields_fromTree( "highAll_H_mass"   , "h_mass"   , selection, 50, 100., 200., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // //  selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && gg_nBJets<2 && (r9Gamma0<0.94 || r9Gamma1<0.94 ) && ht<200  && gg_nJets<4 ";
//   // // dt.drawRegionYields_fromTree( "lowAll_H_mass"   , "h_mass"   , selection, 50, 100., 200., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");





//   // selection = " (( ptGamma0/h_mass) > 0.33) && ((ptGamma1/h_mass )> 0.25) && fabs(etaGamma0)<=1.4  && fabs(etaGamma1)<=1.4 && gg_deltaPhiMin<2.5 && gg_nBJets>=2 && r9Gamma0>0.94 && r9Gamma1>0.94 && ( h_pt/h_mass) >= 0.8 && gg_nJets>=4 ";
//   // dt.drawRegionYields_fromTree( "ggDPhi2p5_highAll_H_mass"   , "h_mass"   , selection, 50, 100., 200., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");








//   // selection = "fabs(etaGamma0)<=1.4 && fabs(etaGamma1)<=1.4 && r9Gamma0<1.1 && r9Gamma1<1.1 && sigmaIetaIetaGamma0<0.011 && sigmaIetaIetaGamma1<0.011  && sigmaIetaIetaGamma0>0.005 && sigmaIetaIetaGamma1>0.005   && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) ";

//   // dt.drawRegionYields_fromTree( "purer_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

//   // selection = "fabs(etaGamma0)<=1.4 && fabs(etaGamma1)<=1.4     && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) ";

//   // dt.drawRegionYields_fromTree( "EBEB_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

//   // selection = "fabs(etaGamma0)<=1.4 || fabs(etaGamma1)<=1.4  && !( fabs(etaGamma0)<=1.4 && fabs(etaGamma1)<=1.4  )      && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 
//   // dt.drawRegionYields_fromTree( "EBEE_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

//   // selection = "fabs(etaGamma0)>1.4 && fabs(etaGamma1)>1.4   && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 
//   // dt.drawRegionYields_fromTree( "EEEE_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");


//   // selection = "ht>400 && fabs(etaGamma0)<=1.4 && fabs(etaGamma1)<=1.4  && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 

//   // dt.drawRegionYields_fromTree( "ht400_EBEB_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_EBEB_nJets" , "nJets" , selection, 13, 0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_EBEB_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");


//   // selection = "ht>400 && fabs(etaGamma0)<=1.4 && fabs(etaGamma1)<=1.4 && nBJets>=3  && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 
//   // dt.drawRegionYields_fromTree( "ht400_EBEB_nBJets3_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_EBEB_nBJets3_nJets" , "nJets" , selection, 13, 0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");


//   // selection = "ht>400 && fabs(etaGamma0)<=1.4 && fabs(etaGamma1)<=1.4 && nBJets>=2    && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 
//   // dt.drawRegionYields_fromTree( "ht400_EBEB_nBJets2_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_EBEB_nBJets2_nJets" , "nJets" , selection, 13, 0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");


//   // selection = "ht>400 && fabs(etaGamma0)<=1.4 && fabs(etaGamma1)<=1.4 && nBJets>=2 && nJets>=6  && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 
//   // dt.drawRegionYields_fromTree( "ht400_EBEB_nBJets2_nJets6_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");


//   // selection = "ht>400 && fabs(etaGamma0)<=1.4 || fabs(etaGamma1)<=1.4  && !( fabs(etaGamma0)<=1.4 && fabs(etaGamma1)<=1.4  )   && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 

//   // dt.drawRegionYields_fromTree( "ht400_EBEE_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");


//   // selection = "ht>400 && fabs(etaGamma0)>1.4 && fabs(etaGamma1)>1.4     && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 

//   // dt.drawRegionYields_fromTree( "ht400_EEEE_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");




//   // selection = "ht>400 && h_mass>90 && h_mass <180    && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 


//   // dt.drawRegionYields_fromTree( "ht400_mt2"   , "mt2"   , selection, 50, 0., 1000., "M_{T2}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_met"   , "met"   , selection, 50, 0, 1000, "ME_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_ht"    , "ht"    , selection, 50, 0., 2000., "H_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

//   // dt.drawRegionYields_fromTree( "ht400_deltaPhiMin"    , "deltaPhiMin"    , selection, 50, 0., 3.2, "deltaPhiMin", "", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_diffMetMht"    , "diffMetMht"    , selection, 50, 0., 200., "diffMetMht", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_H_pt"   , "h_pt"   , selection, 50, 10., 1010., "H p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");

//   // dt.drawRegionYields_fromTree( "ht400_H_mass"   , "h_mass"   , selection, 50, 80., 180., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");

//   // dt.drawRegionYields_fromTree( "ht400_nJets" , "nJets" , selection, 13, 0.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", cutsLabel, "#geq0j, #geq0b");

//   // dt.drawRegionYields_fromTree( "ht400_gamma_pt0"   , "ptGamma0"   , selection, 40, 0., 700., "Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_gamma_pt1"   , "ptGamma1"   , selection, 40, 0., 400., "Sub-Leading #gamma p_{T}", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_gamma_eta0"   , "etaGamma0"   , selection,30 , -3., 3., "Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_gamma_eta1"   , "etaGamma1"   , selection, 30, -3., 3., "Sub-Leading #gamma #eta", "", cutsLabel, "#geq0j, #geq0b");

 
//   // dt.drawRegionYields_fromTree( "ht400_gamma_r90"   , "r9Gamma0"   , selection, 40, 0., 2., "Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_gamma_r91"   , "r9Gamma1"   , selection, 40, 0., 2., "Sub-Leading #gamma r9", "GeV", cutsLabel, "#geq0j, #geq0b");

//   // dt.drawRegionYields_fromTree( "ht400_gamma_sigmaIetaIeta0"   , "sigmaIetaIetaGamma0"   , selection, 40, 0., 0.04, "Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_gamma_sigmaIetaIeta1"   , "sigmaIetaIetaGamma1"   , selection, 40, 0., 0.04, "Sub-Leading #gamma #sigmaI#etaI#eta", "GeV", cutsLabel, "#geq0j, #geq0b");

//   // dt.drawRegionYields_fromTree( "ht400_gamma_chHadIso0"   , "chHadIsoGamma0"   , selection, 40, 0., 2.5, "Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");
//   // dt.drawRegionYields_fromTree( "ht400_gamma_chHadIso1"   , "chHadIsoGamma1"   , selection, 40, 0., 2.5, "Sub-Leading #gamma chHadIso", "GeV", cutsLabel, "#geq0j, #geq0b");


//   // selection = "ht>400 && (fabs(etaGamma0)<=1.4 ||  fabs(etaGamma1)<=1.4) && (r9Gamma0>0.85 || sigmaIetaIetaGamma0<0.015)  && (r9Gamma1>0.85 ||  sigmaIetaIetaGamma1<0.015 )  && mt2<200 && met<300 && deltaPhiMin<1.5  && ( ((sigmaIetaIetaGamma0<= 0.0103 && fabs(etaGamma0)<=1.4)) || (sigmaIetaIetaGamma0< 0.0277 && fabs(etaGamma0)>1.4)) && ((sigmaIetaIetaGamma1<= 0.0103 && fabs(etaGamma1)<=1.4 ) || (sigmaIetaIetaGamma1< 0.0277 && fabs(etaGamma1)>1.4 )) "; 


//   // //selection = "ht>400 && (fabs(etaGamma0)<=1.4 ||  fabs(etaGamma1)<=1.4) && (r9Gamma0>0.85 || sigmaIetaIetaGamma0<0.015)  && (r9Gamma1>0.85 ||  sigmaIetaIetaGamma1<0.015 )  && mt2<200 && met<300 && nJets>=4 && nBJets>=2 && deltaPhiMin<1.5";

//   // dt.drawRegionYields_fromTree( "purer2_H_mass"   , "h_mass"   , selection, 50, 80., 280., "M_{#gamma#gamma}", "GeV", cutsLabel, "#geq0j, #geq0b");


  return 0;

}





std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units ) {

  std::string cutLabel;
  if( theMax>theMin ) cutLabel = std::string(Form("%.0f < %s < %.0f %s", theMin, name.c_str(), theMax, units.c_str()) );
  else                cutLabel = std::string(Form("%s > %.0f %s", name.c_str(), theMin, units.c_str()) );

  return cutLabel;

}



void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units ) {


  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;

  std::vector<int> colors;
  if( bgYields.size()==3 ) { // estimates
    colors.push_back(402); 
    colors.push_back(430); 
    colors.push_back(418); 
  } else { // mc
    colors.push_back(430); // other=zll
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    //    colors.push_back(419); // z+jets
    colors.push_back(855); // top
    //colors.push_back(); // other
  }

  std::string fullPathPlots = cfg.getEventYieldDir() + "/plotsDataMC";
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


    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<bgYields.size(); ++i ) { 
      TTree* tree_mc = (bgYields[i]->get(thisRegion)->tree);
      std::string thisName = "h1_" + bgYields[i]->getName();
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_mc->Sumw2();
      if( selection!="" )
	//tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%s/puWeight", selection.c_str()) );
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%s", selection.c_str()) );
      else
        tree_mc->Project( thisName.c_str(), varName.c_str(), "" );
      histos_mc.push_back(h1_mc);
    }

    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

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
    if( shapeNorm ) 
      std::cout << "SF: " << scaleFactor << std::endl;

    TH1D* histo_mc;
    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = bgYields.size() - i - 1;
      histos_mc[index]->SetFillColor( colors[index] );
      histos_mc[index]->SetLineColor( kBlack );
      if( shapeNorm ) {
        histos_mc[index]->Scale( scaleFactor );
      }
      else
	histos_mc[index]->Scale( 145.0/106.5 );

      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack.Add(histos_mc[index]);
    }


    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();
    TPad *pad1 = MT2DrawTools::getCanvasMainPad();
    TPad *pad2 = MT2DrawTools::getCanvasRatioPad();
        
    TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600);
    c1_log->cd();
    TPad *pad1_log = MT2DrawTools::getCanvasMainPad( true );
    TPad *pad2_log = MT2DrawTools::getCanvasRatioPad( true );
 
    float yMaxScale = 1.1;
    float yMax1 = h1_data->GetMaximum()*yMaxScale;
    float yMax2 = yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
    float yMax3 = yMaxScale*(bgStack.GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    if( h1_data->GetNbinsX()<2 ) yMax *=3.;

    std::string xAxisTitle;
    if( units!="" ) 
      xAxisTitle = (std::string)(Form("%s [%s]", axisName.c_str(), units.c_str()) );
    else
      xAxisTitle = (std::string)(Form("%s", axisName.c_str()) );

    std::string binWidthText;
    if( binWidth>=1. )         binWidthText = (std::string)Form("%.0f", binWidth);
    else if( binWidth>=0.1 )   binWidthText = (std::string)Form("%.1f", binWidth);
    else if( binWidth>=0.01 )  binWidthText = (std::string)Form("%.2f", binWidth);
    else if( binWidth>=0.001 ) binWidthText = (std::string)Form("%.3f", binWidth);
    else                       binWidthText = (std::string)Form("%.4f", binWidth);

    std::string yAxisTitle;
    if( units!="" ) 
      yAxisTitle = (std::string)(Form("Events / (%s %s)", binWidthText.c_str(), units.c_str()));
    else
      yAxisTitle = (std::string)(Form("Events / (%s)", binWidthText.c_str()));

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
    pad1->Draw();
    pad1->cd();
    h2_axes->Draw();
   
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.1, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();
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
      regionText->AddText( niceNames[i].c_str() );

      pad1->cd();
      // regionText->Draw("same");
      pad1_log->cd();
      //  regionText->Draw("same");
    }
    
    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      normText->AddText( "#splitline{Shape}{Norm.}" );
      pad1->cd();
      //normText->Draw("same");
      pad1_log->cd();
      //normText->Draw("same");
    }

    TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( gr_data, "Data", "P" );
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], bgYields[i]->getFullName().c_str(), "F" );
    }


    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    
    TPaveText* ratioText = new TPaveText( 0.133, -0.051, 0.4, 0.1 , "brNDC" );
    ratioText->SetTextSize(0.04);
    ratioText->SetTextFont(40);
    ratioText->SetTextColor(2);
    ratioText->SetFillColor(0);
    ratioText->SetTextAlign(11);
    ratioText->AddText( Form("Data/MC = %.2f", scaleFactor) );
    //  ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", scaleFactor, error_datamc) );
     

    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    
    TLine* lineSF = new TLine(xMin, scaleFactor, xMax, scaleFactor);
    lineSF->SetLineColor(2);

    float yMinR=0.0;
    float yMaxR=2.0;

    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
 
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);
   
    //    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );
    TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TGraphErrors* SFFitBand = MT2DrawTools::getSFFitBand(fSF, xMin, xMax);
    TPaveText* fitText = MT2DrawTools::getFitText( fSF );


    c1->cd();
    pad1->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
      fitText->Draw("same");
    // ratioText->Draw("same");
  
    gPad->RedrawAxis();

    c1_log->cd();
    pad1_log->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
     fitText->Draw("same");
    //  ratioText->Draw("same");

    gPad->RedrawAxis();

   /*
    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    TLine* lineSF = MT2DrawTools::getSFLine(integral_data, integral_mc, xMin, xMax);
    TGraphErrors* SFband = MT2DrawTools::getSFBand(integral_data, error_data, integral_mc, error_mc, xMin, xMax);
    */

    c1->cd();
    //   TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
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

      SFFitBand->Draw("3,same");
      fSF->Draw("same");
    }

    g_ratio->Draw("PE,same");    
    gPad->RedrawAxis();


    c1_log->cd();
    // TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
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
    /*
    line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same"); */
    g_ratio->Draw("PE,same");
    gPad->RedrawAxis();


    c1->SaveAs( Form("%s/%s_%s.eps", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.png", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.pdf", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );

    c1_log->SaveAs( Form("%s/%s_%s_log.eps", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1_log->SaveAs( Form("%s/%s_%s_log.png", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1_log->SaveAs( Form("%s/%s_%s_log.pdf", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );

    delete c1;
    delete h2_axes;

    delete c1_log;
    delete h2_axes_log;

    delete h2_axes_ratio;

    delete h1_data;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

  }// for MT2 regions

}














































void drawShape( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, const std::string& saveName, const std::string& varName, const std::string& selection, const std::string& selection2, const std::string& selection3, unsigned int nSel, int nBins, float xMin, float xMax, std::string axisName, const std::string& units ) {

  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;

  std::vector<int> colors;

  colors.push_back(430); // other=zll
  colors.push_back(401); // qcd
  // colors.push_back(417); // w+jets
  // colors.push_back(419); // z+jets
  // colors.push_back(855); // top
  colors.push_back(97); // other


  std::string fullPathPlots = cfg.getEventYieldDir() + "/plotShapeComp";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_mc = data->get(thisRegion)->tree;
    // TH1D* h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
    // tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );

    // TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    // gr_data->SetMarkerStyle(20);
    // gr_data->SetMarkerSize(1.2);

    float lumi = cfg.lumi();
    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<nSel; ++i ) { 
      std::string thisName = "h1_" + i;
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_mc->Sumw2();
      if(i==0){
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi, selection.c_str()) );
      }else if(i==1){
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi, selection2.c_str()) );
      }else if(i==2){
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi, selection3.c_str()) );
      }
      h1_mc->SetLineColor( colors[i] );
      h1_mc->SetLineWidth( 2 );

      h1_mc->Scale( 1./ h1_mc->Integral() );
      histos_mc.push_back(h1_mc);
    }

    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();
    // TPad *pad1 = MT2DrawTools::getCanvasMainPad();
    // TPad *pad2 = MT2DrawTools::getCanvasRatioPad(false);
        
    TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600);
    c1_log->cd();
    c1_log->SetLogy();
    // TPad *pad1_log = MT2DrawTools::getCanvasMainPad( true );
    // TPad *pad2_log = MT2DrawTools::getCanvasRatioPad( false );
 
    float yMaxScale = 1.1;
    float yMax1 = histos_mc[0]->GetMaximum()*yMaxScale;
    float yMax2 = histos_mc[1]->GetMaximum()*yMaxScale;
    float yMax3 = histos_mc[2]->GetMaximum()*yMaxScale;

    float yMax = yMax1;
    if( yMax2 > yMax1 ) yMax = yMax2;
    if( yMax3 > yMax ) yMax = yMax3;

    // float yMax2 = yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
    // float yMax3 = yMaxScale*(bgStack.GetMaximum());
    // float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    // if( yMax3 > yMax ) yMax = yMax3;
    // if( h1_data->GetNbinsX()<2 ) yMax *=3.;

    std::string xAxisTitle;
    if( units!="" ) 
      xAxisTitle = (std::string)(Form("%s [%s]", axisName.c_str(), units.c_str()) );
    else
      xAxisTitle = (std::string)(Form("%s", axisName.c_str()) );

    std::string binWidthText;
    if( binWidth>=1. )         binWidthText = (std::string)Form("%.0f", binWidth);
    else if( binWidth>=0.1 )   binWidthText = (std::string)Form("%.1f", binWidth);
    else if( binWidth>=0.01 )  binWidthText = (std::string)Form("%.2f", binWidth);
    else if( binWidth>=0.001 ) binWidthText = (std::string)Form("%.3f", binWidth);
    else                       binWidthText = (std::string)Form("%.4f", binWidth);

    std::string yAxisTitle;
    if( units!="" ) 
      yAxisTitle = (std::string)(Form("Events / (%s %s)", binWidthText.c_str(), units.c_str()));
    else
      yAxisTitle = (std::string)(Form("Events / (%s)", binWidthText.c_str()));

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
    c1->Draw();
    //   c1->cd();
    h2_axes->Draw();
   
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.001, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();
    c1_log->Draw();
    //   c1_log->cd();
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
      regionText->AddText( niceNames[i].c_str() );

      c1->cd();
      // regionText->Draw("same");
      c1_log->cd();
      //  regionText->Draw("same");
    }
    
    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      normText->AddText( "#splitline{Shape}{Norm.}" );
      c1->cd();
      //normText->Draw("same");
      c1_log->cd();
      //normText->Draw("same");
    }

    TLegend* legend = new TLegend( 0.7, 0.9-(histos_mc.size())*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    //    legend->AddEntry( gr_data, "Data", "P" );
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], Form("Bin %d",i), "F" );
      //legend->AddEntry( histos_mc[i], bgYields[i]->getFullName().c_str(), "F" );
    }


    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    
    // TPaveText* ratioText = new TPaveText( 0.133, -0.051, 0.4, 0.1 , "brNDC" );
    // ratioText->SetTextSize(0.04);
    // ratioText->SetTextFont(40);
    // ratioText->SetTextColor(2);
    // ratioText->SetFillColor(0);
    // ratioText->SetTextAlign(11);
    // ratioText->AddText( Form("Data/MC = %.2f", scaleFactor) );
    // //  ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", scaleFactor, error_datamc) );
     

    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    
    // TLine* lineSF = new TLine(xMin, scaleFactor, xMax, scaleFactor);
    // lineSF->SetLineColor(2);

    // float yMinR=0.0;
    // float yMaxR=2.0;

    // TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
    // TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
 
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    //    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);
   
    //    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );
    //    TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    // TGraphErrors* SFFitBand = MT2DrawTools::getSFFitBand(fSF, xMin, xMax);
    // TPaveText* fitText = MT2DrawTools::getFitText( fSF );


    c1->cd();
    c1->cd();
    legend->Draw("same");
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      histos_mc[i]->Draw("histo same");
    }
    // bgStack.Draw("histo same");
    //    gr_data->Draw("p same");
    labelTop->Draw("same");
    // if( !shapeNorm )
    //   fitText->Draw("same");
    // ratioText->Draw("same");
  
    gPad->RedrawAxis();

    c1_log->cd();
    c1_log->cd();
    legend->Draw("same");
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      histos_mc[i]->Draw("histo same");
    }
    // bgStack.Draw("histo same");
    // gr_data->Draw("p same");
    labelTop->Draw("same");
    // if( !shapeNorm )
    //  fitText->Draw("same");
    //  ratioText->Draw("same");

    gPad->RedrawAxis();

   /*
    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    TLine* lineSF = MT2DrawTools::getSFLine(integral_data, integral_mc, xMin, xMax);
    TGraphErrors* SFband = MT2DrawTools::getSFBand(integral_data, error_data, integral_mc, error_mc, xMin, xMax);
    */

    // c1->cd();
    // //   TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
    // pad2->Draw();
    // pad2->cd();

    //    h2_axes_ratio->Draw("");
 
    /*  line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same");
    */
    // lineCentral->Draw("same");
    // if( !shapeNorm ){

    //   systBand->Draw("3,same");
    //   lineCentral->Draw("same");

    //   SFFitBand->Draw("3,same");
    //   fSF->Draw("same");
    //}

    //    g_ratio->Draw("PE,same");    
  // gPad->RedrawAxis();


    // c1_log->cd();
    // // TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
    // pad2_log->Draw();
    // pad2_log->cd();

    // //    h2_axes_ratio->Draw(""); 
     
    // lineCentral->Draw("same");
    // if( !shapeNorm ){

    //   systBand->Draw("3,same");
    //   lineCentral->Draw("same");

    //   SFFitBand->Draw("3,same");
    //   fSF->Draw("same");
    // }
    // /*
    // line->Draw("same");
    // SFband->Draw("3,same");
    // lineSF->Draw("same"); */
    // //    g_ratio->Draw("PE,same");
    // gPad->RedrawAxis();


    c1->SaveAs( Form("%s/%s_%s.eps", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.png", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.pdf", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );

    c1_log->SaveAs( Form("%s/%s_%s_log.eps", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1_log->SaveAs( Form("%s/%s_%s_log.png", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1_log->SaveAs( Form("%s/%s_%s_log.pdf", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );

    delete c1;
    delete h2_axes;

    delete c1_log;
    delete h2_axes_log;

    // delete h2_axes_ratio;

    // delete h1_data;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

  }// for MT2 regions

}
