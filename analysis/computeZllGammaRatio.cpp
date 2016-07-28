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
#include "TGraphErrors.h"

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2Config.h"
#include "interface/MT2DrawTools.h"

#include <iostream>

#define mt2_cxx
#include "interface/mt2.h"


double lumiErr = 0.046;

void drawRatios(std::string fullPath, double *binss, unsigned int size,  std::string selection, MT2Analysis<MT2EstimateSyst>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma_mc, MT2Analysis<MT2EstimateTree>*  gamma_data, MT2Analysis<MT2EstimateSyst>*  purity, MT2Analysis<MT2EstimateTree>*  zll_mc,MT2Analysis<MT2EstimateTree>*  zll_data,   MT2Analysis<MT2EstimateTree>*  top, MT2Analysis<MT2EstimateSyst>*  zll_yield, MT2Analysis<MT2EstimateSyst>*  zllG_data,  MT2Analysis<MT2EstimateSyst>*  zllG_mc, const MT2Region thisRegion, std::string cut,  std::string cut_gamma, std::string cut_data, std::string cut_gamma_data, float lumi, std::string saveName, bool onlyMC, float top_SF, bool fullUncert , std::string topoCuts="" );

void drawRatiosTopHisto(std::string fullPath, double *binss, unsigned int size,  std::string selection, MT2Analysis<MT2EstimateSyst>*  zll_ratio, MT2Analysis<MT2EstimateTree>* gamma_mc, MT2Analysis<MT2EstimateTree>* gamma_data, MT2Analysis<MT2EstimateSyst>* purity, MT2Analysis<MT2EstimateTree>*  zll_mc, MT2Analysis<MT2EstimateTree>* zll_data, MT2Analysis<MT2Estimate>* top, MT2Analysis<MT2EstimateSyst>* zll_yield, MT2Analysis<MT2EstimateSyst>*  zllG_data, MT2Analysis<MT2EstimateSyst>* zllG_mc, const MT2Region thisRegion, std::string cut, std::string cut_gamma, std::string cut_data, std::string cut_gamma_data, float lumi, std::string saveName, bool onlyMC, float scaleFactor, bool fullUncert , std::string topoCuts="" );

int main(int argc, char* argv[]){

  
  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeZllGammaRatio             |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<1) {
    std::cout << "USAGE: ./computeZllGammaRatio [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  bool onlyData = false;
  bool onlyMC   = false;
  if( argc > 2 ) {
    std::string dataMC(argv[2]);
    if( dataMC=="data" ) onlyData = true;
    else if( dataMC=="MC" ) onlyMC = true;
    else {
      std::cout << "-> You passed a second argument that isn't 'data' nor 'MC', so I don't know what to do about it." << std::endl;
    }
  }

  std::string regionsSet = cfg.crRegionsSet();
  std::cout << "-> Using regions: " << regionsSet << std::endl;

  std::string samples = cfg.mcSamples();

  std::string outputdir = cfg.getEventYieldDir() + "/zllGammaRatio";
  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string gammaControlRegionDir = cfg.getEventYieldDir() + "/gammaControlRegion";
  std::string ZllDir = cfg.getEventYieldDir() + "/zllControlRegion";

  gStyle->SetOptTitle(0);
  MT2DrawTools::setStyle();

  // NO TOP SCALING AT THE MOMENT, will be estimated by Opposite Flavor anyway
  // ifstream SF_file;
  // SF_file.open(Form("%s/plotsDataMCscaling/scaleFactorOF.txt", cfg.getEventYieldDir().c_str() ) );
  float scaleFactor;
  // SF_file >> scaleFactor;
  scaleFactor= 1;
  // scaleFactor= 0.65;
  std::cout<< "Scale Factor = "  << scaleFactor << std::endl;
  


  MT2Analysis<MT2EstimateTree>* gamma_mc = MT2Analysis<MT2EstimateTree>::readFromFile( gammaControlRegionDir + "/mc.root", "gammaCRtree" );
  MT2Analysis<MT2EstimateTree>* gamma_data = MT2Analysis<MT2EstimateTree>::readFromFile( gammaControlRegionDir + "/data.root", "gammaCRtree" );
  

  MT2Analysis<MT2EstimateSyst>* purity_central = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_central_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_mono_ht = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_mono_ht_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_incl_ht = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_incl_ht_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_incl_njets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_incl_njets_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_incl_nbjets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_incl_nbjets_data.root", "purity");

  //  MT2Analysis<MT2EstimateSyst>* purity_multijet_njets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_multijet_njets_data.root", "purity");

  if(onlyMC){
    gamma_data = MT2Analysis<MT2EstimateTree>::readFromFile(Form( "%s/mc.root", gammaControlRegionDir.c_str() ), "gammaCRtree");
  }

  if( gamma_mc==0 ) {
    std::cout << "-> Please run gammaControlRegion first. I need to get the gammaCR yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(193);
  }
 
  MT2Analysis<MT2EstimateTree>* zll_mc = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc.root", ZllDir.c_str()) , "zllCR");
  MT2Analysis<MT2EstimateTree>* zll_data = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ZllDir.c_str()) , "data");
  
  if(onlyMC){
    zll_data = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc.root", ZllDir.c_str()) , "zllCR");
  }
  if( zll_data==0 ) {
    std::cout << "-> Please run computeZinvFromZll first. I need to get the Z->ll MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }

  MT2Analysis<MT2Estimate>* top_ht = MT2Analysis<MT2Estimate>::readFromFile(Form("%s/zllTop/top.root", cfg.getEventYieldDir().c_str()) , "ht");
  MT2Analysis<MT2Estimate>* top_mono_ht = MT2Analysis<MT2Estimate>::readFromFile(Form("%s/zllTop/top.root", cfg.getEventYieldDir().c_str()) , "ht_mono");
  MT2Analysis<MT2Estimate>* top_nJets = MT2Analysis<MT2Estimate>::readFromFile(Form("%s/zllTop/top.root", cfg.getEventYieldDir().c_str()) , "nJets");
  MT2Analysis<MT2Estimate>* top_nBJets = MT2Analysis<MT2Estimate>::readFromFile(Form("%s/zllTop/top.root", cfg.getEventYieldDir().c_str()) , "nBJets");
  MT2Analysis<MT2Estimate>* top_central = MT2Analysis<MT2Estimate>::readFromFile(Form("%s/zllTop/top.root", cfg.getEventYieldDir().c_str()) , "central");


  MT2Analysis<MT2EstimateTree>* top = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "Top");
  // MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str()  ), "QCD");
  // MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "WJets");
  // wjets->setFullName("W+jets");
 
  std::cout << "Everything read in fine" << std::endl; 
  //YIELDS
  MT2Analysis<MT2EstimateSyst>* zll_mt2 = new MT2Analysis<MT2EstimateSyst>( "zll_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zll_ht = new MT2Analysis<MT2EstimateSyst>( "zll_ht", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zll_nJets = new MT2Analysis<MT2EstimateSyst>( "zll_nJets",regionsSet.c_str());
  //MT2Analysis<MT2EstimateSyst>* zll_nJets = new MT2Analysis<MT2EstimateSyst>( "zll_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zll_nBJets = new MT2Analysis<MT2EstimateSyst>("zll_nBJets",regionsSet.c_str());
  MT2Analysis<MT2EstimateSyst>* zll_mono_ht = new MT2Analysis<MT2EstimateSyst>("zll_mono_ht",regionsSet.c_str());
  MT2Analysis<MT2EstimateSyst>* zll_central = new MT2Analysis<MT2EstimateSyst>("zll_central",regionsSet.c_str());
 
  //DATA RATIOS
  MT2Analysis<MT2EstimateSyst>* zllG_data_mt2 = new MT2Analysis<MT2EstimateSyst>( "zllG_data_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_data_ht = new MT2Analysis<MT2EstimateSyst>( "zllG_data_ht", regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_data_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_data_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_data_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_data_nBJets",regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zllG_data_mono_ht = new MT2Analysis<MT2EstimateSyst>("zllG_data_mono_ht",regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zllG_data_central = new MT2Analysis<MT2EstimateSyst>("zllG_data_central",regionsSet.c_str() );

  //MC RATIOS
  MT2Analysis<MT2EstimateSyst>* zllG_mc_mt2 = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_mc_ht = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_ht", regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_mc_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_mc_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_mc_nBJets",regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zllG_mc_mono_ht = new MT2Analysis<MT2EstimateSyst>("zllG_mc_mono_ht",regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zllG_mc_central = new MT2Analysis<MT2EstimateSyst>("zllG_mc_central",regionsSet.c_str() );

  //Double RATIOS
  MT2Analysis<MT2EstimateSyst>* zllG_mt2 = new MT2Analysis<MT2EstimateSyst>( "zllG_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_ht = new MT2Analysis<MT2EstimateSyst>( "zllG_ht", regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_nBJets",regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zllG_mono_ht = new MT2Analysis<MT2EstimateSyst>("zllG_mono_ht",regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zllG_central = new MT2Analysis<MT2EstimateSyst>("zllG_central",regionsSet.c_str() );
  
  

  std::set<MT2Region> MT2Regions = zll_mc->getRegions();

 
  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  TH1F::AddDirectory(kTRUE);


  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
    std::vector<std::string> niceNames = thisRegion.getNiceNames();
    float lumi = cfg.lumi();    

    // if( !onlyMC){
    //  bins_mt2 = {200,300,400,500, 600, 1500 };
    // }else{
    //  bins_mt2 = {200,300,400,500, 600, 800, 1000, 1500 };
    //    }
    //    float bins_ht[] =  {200,250,300,350,400,450,500,550,600,700,800,900,1000,1500,2000};

    double bins_nJets[] = {1,2,4,7,12};
    int size_nJets = sizeof(bins_nJets)/sizeof(double)-1;

    double bins_nBJets[] = {0,1,2,3,6}; 
    int size_nBJets = sizeof(bins_nBJets)/sizeof(double)-1;

    double bins_mt2[] ={200,300,400,500, 600, 1500 };
    int size_mt2 = sizeof(bins_mt2)/sizeof(double)-1;

    double bins_ht[] =  {200,450,575,1000,1500,3000};
    int size_ht = sizeof(bins_ht)/sizeof(double)-1;
 
    double bins_mono_ht[] = {200, 250, 350,450, 575, 700, 1000, 1500}; 
    int size_mono_ht = sizeof(bins_mono_ht)/sizeof(double)-1;

    double bins_central[] = {200, 3000}; 
    int size_central = sizeof(bins_central)/sizeof(double)-1;



    
    //WITH THE loose electron ID // HLT and lepton scale factors are handled in the loop
    std::string cut_incl = "weight*(abs(Z_mass-91.19)<10 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && mt2>200 && ht>200 && nJets>0 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)))";
    //    std::string cut_incl = "weight*(abs(Z_mass-91.19)<10 && met>200 && mt2>200 && ht>200 && nJets>0 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)))";
    std::string cut_incl_data = cut_incl;

    std::string cut_multijet = "weight*(abs(Z_mass-91.19)<10 && met>200 && mt2>200 && ht>200 && nJets>1 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)))";
    std::string cut_multijet_data = cut_multijet;
    std::string cut_multijet_gamma = "weight*(prompt==2 && iso<2.5  && met>200&&  ptGamma>180 && nJets>1 && mt2>200 && ht>200 )*1.23";
    std::string cut_multijet_gamma_data = "weight*( iso<2.5 && ptGamma>180  && met>200 && nJets>0 && mt2>200 && ht>200 )";
  

    std::string cut_incl_gamma = "weight*(prompt==2 && iso<2.5  && ((ht>200 && met>200)||(ht>1000. && met>30.)) &&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23";
    std::string cut_incl_gamma_data = "weight*( iso<2.5 && ptGamma>180 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && nJets>0 && mt2>200 && ht>200 )";
    std::string cut_incl_el = "weight*(abs(Z_mass-91.19)<10 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && mt2>200 && ht>200 && nJets>0 && lep_tightId0>0 && lep_tightId1>0 && Z_lepId==11 )";
    std::string cut_incl_data_el = cut_incl_el;

    std::string cut_incl_mu = "weight*(abs(Z_mass-91.19)<10 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && mt2>200 && ht>200 && nJets>0 && Z_lepId==13 )";
    std::string cut_incl_data_mu = cut_incl_mu;
  
    /*
    std::string cut_incl_gamma = "weight*(prompt==2 && iso<2.5  && met>200&&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23";
    std::string cut_incl_gamma_data = "weight*( iso<2.5 && ptGamma>180  && met>200 && nJets>0 && mt2>200 && ht>200 )";
  
    std::string cut_incl_el = "weight*(abs(Z_mass-91.19)<10 && met>200 && mt2>200 && ht>200 && nJets>0 && lep_tightId0>0 && lep_tightId1>0 && Z_lepId==11 )";
    std::string cut_incl_data_el = cut_incl_el;

    std::string cut_incl_mu = "weight*(abs(Z_mass-91.19)<10 && met>200 && mt2>200 && ht>200 && nJets>0 && Z_lepId==13 )";
    std::string cut_incl_data_mu = cut_incl_mu;
    */ 



    /*
    if(onlyMC){
      cut_data = Form("weight*(abs(Z_mass-91.19)<20&& nBJets<2  &&  Z_pt>180 && mt2>200 && ht>200 && nJets>1 ) *%f", lumi); 
      cut_el_data =  Form("weight*(abs(Z_mass-91.19)<20 && nBJets<2 && ht>200 && mt2>200&&  Z_pt>180 && nJets>1&& Z_lepId==11 )*%f", lumi);
      cut_mu_data =  Form("weight*(abs(Z_mass-91.19)<20 && nBJets<2 && ht>200&& mt2>200 &&  Z_pt>180 && nJets>1&& Z_lepId==13 ) *%f", lumi);
      cut_gamma_data =  Form("weight*( prompt==2 && nBJets<2 && ptGamma>180 && nJets>1 && mt2>200&&ht>200 )*1.23*%f", lumi);
      cut_mono_data =  Form("weight*( abs(Z_mass-91.19)<20 && nBJets<2 &&  Z_pt>180 && mt2>200&& ht>200 && nJets==1 )*%f", lumi);
      cut_gamma_mono_data =  Form("  weight*( prompt==2 && nBJets<2 && ptGamma>180 && nJets==1 && mt2>200&&ht>200 )*1.23*%f", lumi);
    }
    */

    
    MT2EstimateSyst::rebinYields( zll_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zll_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zll_mono_ht,  size_mono_ht, bins_mono_ht);
    MT2EstimateSyst::rebinYields( zll_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zll_nBJets,  size_nBJets, bins_nBJets);

    MT2EstimateSyst::rebinYields( zllG_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zllG_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zllG_mono_ht,  size_mono_ht, bins_mono_ht);
    MT2EstimateSyst::rebinYields( zllG_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zllG_nBJets,  size_nBJets, bins_nBJets);
 
    MT2EstimateSyst::rebinYields( zllG_data_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zllG_data_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zllG_data_mono_ht,  size_mono_ht, bins_mono_ht);
    MT2EstimateSyst::rebinYields( zllG_data_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zllG_data_nBJets,  size_nBJets, bins_nBJets);

    MT2EstimateSyst::rebinYields( zllG_mc_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zllG_mc_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zllG_mc_mono_ht,  size_mono_ht, bins_mono_ht);
    MT2EstimateSyst::rebinYields( zllG_mc_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zllG_mc_nBJets,  size_nBJets, bins_nBJets);

    MT2EstimateSyst::rebinYields( zll_central,  size_central, bins_central);
    MT2EstimateSyst::rebinYields( zllG_central,  size_central, bins_central);
    MT2EstimateSyst::rebinYields( zllG_data_central,  size_central, bins_central);
    MT2EstimateSyst::rebinYields( zllG_mc_central,  size_central, bins_central);
    

    std::cout << "Rebinned the yields " << std::endl;

    //draw ratio also fills the ratio and yield estimates
    //outputdir, bins, nbins, var to project, ratio estimate, gamma mc, gamma data, purity gamma, zll mc, zll data, yield estimate, region, cut zll, cut gamma, cut zll data, cut gamma data, lumi, name, flag , topo region);
    
    //  drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data, top, zll_mt2 , zllG_data_mt2, zllG_mc_mt2 , thisRegion, cut_el, cut_gamma, cut_el_data,  cut_gamma_data,  lumi , "mt2_el" , onlyMC, scaleFactor ,"#geq2j, #geq0b" );
    //  drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data, top, zll_mt2, zllG_data_mt2, zllG_mc_mt2 , thisRegion, cut_mu,  cut_gamma, cut_mu_data, cut_gamma_data, lumi , "mt2_mu" , onlyMC, scaleFactor ,"#geq2j, #geq0b");
    // drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data, top, zll_mt2, zllG_data_mt2, zllG_mc_mt2 ,thisRegion, cut, cut_gamma, cut_data,cut_gamma_data, lumi , "mt2" , onlyMC, scaleFactor ,"#geq2j, #geq0b");
   

    //CENTRAL
    drawRatiosTopHisto( outputdir, bins_central, size_central , "ht",   zllG_central,   gamma_mc, gamma_data, purity_central,  zll_mc, zll_data, top_central, zll_central, zllG_data_central, zllG_mc_central , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"central_dd" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    //drawRatiosTopHisto( outputdir, bins_central, size_central , "ht",   zllG_central,   gamma_mc, gamma_data, purity_central,  zll_mc, zll_data, top_central, zll_central, zllG_data_central, zllG_mc_central , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"central_noPFU_dd" , onlyMC, scaleFactor, 0 ,"#geq1j, #geq0b");

    drawRatios( outputdir, bins_central, size_central , "ht",   zllG_central,   gamma_mc, gamma_data, purity_central,  zll_mc, zll_data, top, zll_central, zllG_data_central, zllG_mc_central , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"central" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    drawRatios( outputdir, bins_central, size_central , "ht",   zllG_central,   gamma_mc, gamma_data, purity_central,  zll_mc, zll_data, top, zll_central, zllG_data_central, zllG_mc_central , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"central_noPFU" , onlyMC, scaleFactor, 0 ,"#geq1j, #geq0b");



    //incl_//////////////////////
     drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_incl_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_incl_el, cut_incl_gamma, cut_incl_data_el,  cut_incl_gamma_data, lumi,"ht_incl_el" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
     drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_incl_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_incl_mu, cut_incl_gamma, cut_incl_data_mu,  cut_incl_gamma_data, lumi,"ht_incl_mu" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");

    //DATA DRIVEN TOP ESTIMATE
     // drawRatiosTopHisto( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_incl_ht,  zll_mc, zll_data, top_ht, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"ht_incl_dd" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    //SIMULATION TOP ESTIMATE
    drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_incl_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"ht_incl" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_incl_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"ht_incl_noPFU" , onlyMC, scaleFactor, 0 ,"#geq1j, #geq0b");
    

   //incl_ new HT BINNING FOR MONOJET//////////////////////
      drawRatios( outputdir, bins_mono_ht, size_mono_ht , "ht",   zllG_mono_ht,   gamma_mc, gamma_data, purity_mono_ht,  zll_mc, zll_data, top, zll_mono_ht, zllG_data_mono_ht, zllG_mc_mono_ht , thisRegion, cut_incl_el, cut_incl_gamma, cut_incl_data_el,  cut_incl_gamma_data, lumi,"ht_mono_el" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
      drawRatios( outputdir, bins_mono_ht, size_mono_ht , "ht",   zllG_mono_ht,  gamma_mc, gamma_data, purity_mono_ht,  zll_mc, zll_data, top, zll_mono_ht, zllG_data_mono_ht, zllG_mc_mono_ht , thisRegion, cut_incl_mu, cut_incl_gamma, cut_incl_data_mu,  cut_incl_gamma_data, lumi,"ht_mono_mu" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
  

      //DATA DRIVEN TOP ESTIMATE
      // drawRatiosTopHisto( outputdir, bins_mono_ht, size_mono_ht , "ht",   zllG_mono_ht,   gamma_mc, gamma_data, purity_mono_ht,  zll_mc, zll_data, top_mono_ht, zll_mono_ht, zllG_data_mono_ht, zllG_mc_mono_ht , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"ht_mono_dd" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
      //SIMULATION TOP ESTIMATE
      drawRatios( outputdir, bins_mono_ht, size_mono_ht , "ht",   zllG_mono_ht,   gamma_mc, gamma_data, purity_mono_ht,  zll_mc, zll_data, top, zll_mono_ht, zllG_data_mono_ht, zllG_mc_mono_ht , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"ht_mono" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
      drawRatios( outputdir, bins_mono_ht, size_mono_ht , "ht",   zllG_mono_ht,   gamma_mc, gamma_data, purity_mono_ht,  zll_mc, zll_data, top, zll_mono_ht, zllG_data_mono_ht, zllG_mc_mono_ht , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"ht_mono_noPFU" , onlyMC, scaleFactor, 0 ,"#geq1j, #geq0b");
    




    //incl//////////////////////////////////////
      drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_incl_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion,cut_incl_el, cut_incl_gamma, cut_incl_data_el,  cut_incl_gamma_data, lumi, "nJets_incl_el" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_incl_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion,cut_incl_mu, cut_incl_gamma, cut_incl_data_mu,  cut_incl_gamma_data, lumi, "nJets_incl_mu" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");

    //DATA DRIVEN TOP ESTIMATE
    //drawRatiosTopHisto( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_incl_njets,  zll_mc, zll_data, top_nJets,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion,cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi, "nJets_incl_dd" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    //SIMULATION TOP ESTIMATE
    drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_incl_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion,cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi, "nJets_incl" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_incl_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion,cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi, "nJets_incl_noPFU" , onlyMC, scaleFactor, 0 ,"#geq1j, #geq0b");



    //incl///////////////////
    drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_incl_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_incl_el, cut_incl_gamma ,  cut_incl_data_el,  cut_incl_gamma_data, lumi, "nBJets_incl_el" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_incl_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_incl_mu, cut_incl_gamma ,  cut_incl_data_mu,  cut_incl_gamma_data, lumi, "nBJets_incl_mu" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");

    //DATA DRIVEN TOP ESTIMATE
    //drawRatiosTopHisto( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_incl_nbjets,  zll_mc, zll_data, top_nBJets, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_incl, cut_incl_gamma ,  cut_incl_data,  cut_incl_gamma_data, lumi, "nBJets_incl_dd" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    //SIMULATION TOP ESTIMATE
    drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_incl_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_incl, cut_incl_gamma ,  cut_incl_data,  cut_incl_gamma_data, lumi, "nBJets_incl" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_incl_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_incl, cut_incl_gamma ,  cut_incl_data,  cut_incl_gamma_data, lumi, "nBJets_incl_noPFU" , onlyMC, scaleFactor, 0 ,"#geq1j, #geq0b");



  }
  
 
  std::string outFile = outputdir + "/zll_ratio.root";
  //  zllG_mt2->writeToFile(outFile);
  zllG_ht->writeToFile( outFile );
  zllG_nJets->addToFile( outFile );
  zllG_nBJets->addToFile( outFile );
  zllG_mono_ht->addToFile( outFile );

  std::string outFile_yield = outputdir + "/zll_yield.root";
  //  zll_mt2->writeToFile(outFile_yield);
  zll_ht->writeToFile( outFile_yield );
  zll_nJets->addToFile( outFile_yield );
  zll_nBJets->addToFile( outFile_yield );
  zll_mono_ht->addToFile( outFile_yield );

  std::string outFile_data = outputdir + "/zllG_data_ratio.root";
  //  zllG_data_mt2->writeToFile(outFile_data);
  zllG_data_ht->writeToFile( outFile_data );
  zllG_data_nJets->addToFile( outFile_data );
  zllG_data_nBJets->addToFile( outFile_data );
  zllG_data_mono_ht->addToFile( outFile_data );

  std::string outFile_mc = outputdir + "/zllG_mc_ratio.root";
  //  zllG_mc_mt2->writeToFile(outFile_mc);
  zllG_mc_ht->writeToFile( outFile_mc );
  zllG_mc_nJets->addToFile( outFile_mc );
  zllG_mc_nBJets->addToFile( outFile_mc );
  zllG_mc_mono_ht->addToFile( outFile_mc );


  return 0;

}


void drawRatios(std::string fullPath, double *binss, unsigned int size,  std::string selection, MT2Analysis<MT2EstimateSyst>*  zll_ratio, MT2Analysis<MT2EstimateTree>* gamma_mc, MT2Analysis<MT2EstimateTree>* gamma_data, MT2Analysis<MT2EstimateSyst>* purity, MT2Analysis<MT2EstimateTree>*  zll_mc, MT2Analysis<MT2EstimateTree>* zll_data, MT2Analysis<MT2EstimateTree>* top, MT2Analysis<MT2EstimateSyst>* zll_yield, MT2Analysis<MT2EstimateSyst>*  zllG_data, MT2Analysis<MT2EstimateSyst>* zllG_mc, const MT2Region thisRegion, std::string cut, std::string cut_gamma, std::string cut_data, std::string cut_gamma_data, float lumi, std::string saveName, bool onlyMC, float scaleFactor, bool fullUncert , std::string topoCuts ){
 


  TH1F::AddDirectory(kTRUE);

  double bins[size+1]; 
  for(unsigned int i=0; i<= size ; ++i)  
    bins[i]=binss[i];
  float xMin = binss[0];   float xMax = binss[size];

  //THE TREES
  TTree *zll_data_tree =  zll_data->get(thisRegion)->tree;
  TTree *gamma_data_tree =  gamma_data->get(thisRegion)->tree;
  TTree *zll_mc_tree =  zll_mc->get(thisRegion)->tree;
  TTree *gamma_mc_tree =  gamma_mc->get(thisRegion)->tree;
  TTree *top_tree =  top->get(thisRegion)->tree;
  TGraphAsymmErrors* this_zinv_purity = purity->get(thisRegion)->getGraph();

  /*
  //Getting the weighted center of the bin, somebody at some point thought this was good idea
  float meanX_err[size+1];
  for(unsigned int i = 0; i < size; ++i){
  TH1D* histo_x= new TH1D("histo_x","",100, bins[i], bins[i+1] );
  zll_mc_tree->Project( "histo_x",selection.c_str(), cut.c_str());
  histo_x->GetMean();
  meanX[i] =  histo_x->GetMean();
  meanX_err[i] =  histo_x->GetMeanError();
  delete histo_x;
  }
  */

  TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
  canny->cd();

  TPad* pad1 = 0;
  if( !onlyMC ) {
    pad1 = MT2DrawTools::getCanvasMainPad();
    pad1->Draw();
    pad1->cd();
  }

  TH1D* h_top = new TH1D("h_top","", size , bins); h_top->Sumw2();  
  top_tree ->Project( "h_top" , selection.c_str(), Form("(%s)*%f", cut.c_str(), lumi ) );
  MT2DrawTools::addOverflowSingleHisto(h_top);

  TH1D* h_mt2 = new TH1D("h_mt2","", size , bins); h_mt2->Sumw2();
  zll_data_tree ->Project( "h_mt2" , selection.c_str(), cut_data.c_str() );
  MT2DrawTools::addOverflowSingleHisto(h_mt2);

  TH1D* g_mt2 = new TH1D("g_mt2","", size , bins); g_mt2->Sumw2();
  gamma_data_tree->Project( "g_mt2", selection.c_str(), cut_gamma_data.c_str() );
  MT2DrawTools::addOverflowSingleHisto(g_mt2);


  //GET THE UNCERTAINTIES
  double f = 0.92;
  double f_uncert = 0.1; // 0.08 from the fragmentation and then o+ ~5% for the mc closure
  if( !fullUncert ) f_uncert = 0.0;
  TH1D* g_Up = new TH1D("g_Up","", size , bins);
  TH1D* g_Down = new TH1D("g_Down","", size , bins);

  int nBinss =  h_mt2->GetNbinsX();
  for(int iBin = 1; iBin <= nBinss; iBin++){

    //Zll//////////////
    double value = h_mt2->GetBinContent(iBin);
    double top = h_top->GetBinContent(iBin);

    h_mt2->SetBinContent(iBin, value - top*scaleFactor);
    //add 50% of the top bg as uncertainty to the estimate
    double zll_err =  sqrt(value + top*top*scaleFactor*scaleFactor*0.25);
    //  if( scaleFactor != 1) zll_err =  sqrt(value + top*top*0.05*0.05 );
    h_mt2->SetBinError(iBin,  zll_err  );

    //Photons//////////
    Double_t x_tmp, p, p_errUp, p_errDown;	       
    this_zinv_purity->GetPoint( iBin-1, x_tmp, p);
    p_errUp   = this_zinv_purity->GetErrorYhigh(iBin-1);
    p_errDown = this_zinv_purity->GetErrorYlow(iBin-1);
    if( !fullUncert ){
       p_errUp = 0.;       p_errDown = 0.;
    }


    double value_g = g_mt2->GetBinContent(iBin);
    g_mt2->SetBinContent(iBin, value_g  * p * f);

    std::cout << "purity = " << p << std::endl;
    std::cout << "yield photons = " << value_g << std::endl;
    std::cout << "yield zll = " << value << std::endl;
    std::cout << "yield top = " << top  << std::endl;


    double uncertUp = sqrt ( value_g*f*f*p*p + ( value_g*value_g*f*f *p_errUp*p_errUp) + (p*p*value_g*value_g* f_uncert*f_uncert) );
    double uncertDown = sqrt ( value_g*f*f*p*p + ( value_g*value_g*f*f *p_errDown*p_errDown) + (p*p*value_g*value_g* f_uncert*f_uncert) );
    if( !fullUncert ){
      uncertUp = sqrt ( value_g*f*f*p*p  );
      uncertDown = sqrt ( value_g*f*f*p*p );
    }
    g_Up->SetBinContent(iBin, uncertUp);
    g_Down->SetBinContent(iBin, uncertDown);

  }

  h_mt2->SetMarkerStyle(20);   h_mt2->SetMarkerSize(1.6); h_mt2->SetLineColor(kBlack);
  
 

  //Filling the YIELD
  zll_yield->get(thisRegion);
  int yieldBins = h_mt2->GetNbinsX();

  for(int bi = 1; bi <= yieldBins ; bi++){
    double value = h_mt2->GetBinContent(bi);
    if( value < 0.) value = 0.;
    double err = h_mt2->GetBinError(bi);
    zll_yield->get(thisRegion)->yield->SetBinContent( bi,   value);
    zll_yield->get(thisRegion)->yield_systUp->SetBinContent( bi,   value + err);
    zll_yield->get(thisRegion)->yield_systDown->SetBinContent( bi,   value - err);

    double value_g = g_mt2->GetBinContent(bi);
    double errUp = g_Up->GetBinContent(bi);
    double errDown = g_Down->GetBinContent(bi);
 
    //store the gamma's temporarily in zllG_data
    zllG_data->get(thisRegion)->yield->SetBinContent( bi,   value_g);
    zllG_data->get(thisRegion)->yield_systUp->SetBinContent( bi,   value_g + errUp);
    zllG_data->get(thisRegion)->yield_systDown->SetBinContent( bi,   value_g - errDown );
  }

  *(zllG_data) = *(zll_yield) / *(zllG_data);

  TH1D* h_mt2_mc = new TH1D("h_mt2_mc","", size  , bins);  h_mt2_mc->Sumw2();
  TH1D* g_mt2_mc = new TH1D("g_mt2_mc","", size  , bins);  g_mt2_mc->Sumw2();
  TH1D* h_HLT_weight = new TH1D("h_HLT_weight","", size  , bins);  h_HLT_weight->Sumw2();
  TH1D* h_lepSF = new TH1D("h_lepSF","", size  , bins);  h_lepSF->Sumw2();
  TH1D* h_lepSF_err = new TH1D("h_lepSF_err","", size  , bins);  h_lepSF_err->Sumw2();

  zll_mc_tree ->Project( "h_lepSF", selection.c_str(), Form("(%s)*%f *%s", cut.c_str(), lumi, "(weight_lep0 * weight_lep1)" ) );
  zll_mc_tree ->Project( "h_lepSF_err", selection.c_str(), Form("(%s)*%f *%s", cut.c_str(), lumi, "(weight_lep_err)" ) );
  zll_mc_tree ->Project( "h_HLT_weight", selection.c_str(), Form("(%s)*%f *%s", cut.c_str(), lumi, "(HLT_weight)" ) );
  zll_mc_tree ->Project( "h_mt2_mc", selection.c_str(), Form("(%s)*%f", cut.c_str(), lumi) );
  gamma_mc_tree->Project( "g_mt2_mc", selection.c_str(),  Form("(%s)*%f", cut_gamma.c_str(), lumi)  );

  MT2DrawTools::addOverflowSingleHisto( h_mt2_mc );
  MT2DrawTools::addOverflowSingleHisto( g_mt2_mc );
  MT2DrawTools::addOverflowSingleHisto( h_lepSF );
  MT2DrawTools::addOverflowSingleHisto( h_lepSF_err );
  MT2DrawTools::addOverflowSingleHisto( h_HLT_weight );

 
  //Filling the single ratios
  for(int bi = 1; bi <= yieldBins ; bi++){
    double value = h_mt2_mc->GetBinContent(bi);
    // double value_err = sqrt(value);//h_mt2_mc->GetBinError(bi);
    double value_err = h_mt2_mc->GetBinError(bi);
    double hlt = h_HLT_weight->GetBinContent(bi) /value ;
    double lepSF = h_lepSF->GetBinContent(bi)/value;
    double lepSF_err = h_lepSF_err->GetBinContent(bi)/value;

    lepSF_err -=lepSF;

    std::cout << "MC Zll = " << value  << " +- " << value_err << std::endl;
    std::cout << "MC Gamma= " << g_mt2_mc->GetBinContent(bi) << std::endl;

    // hlt=0.93;
    std::cout << "HTL corr = " << hlt << std::endl;
    std::cout << "id sf = " << lepSF << " +- " << lepSF_err << std::endl;

    //change baaaack
    //h_mt2_mc->SetBinContent(bi, value* hlt ); //LEPTON TRIGGER EFF of 0.92
    //h_mt2_mc->SetBinContent(bi, value * hlt ); //LEPTON TRIGGER EFF of 0.92
    //h_mt2_mc->SetBinError(bi, sqrt( value_err*value_err* 0.92*0.92 +  value*value* 0.05*0.05 ) ); //tempororaa 5% for the trigger eff

    // h_mt2_mc->SetBinError(bi, sqrt( value_err*value_err ) ); //3% for the trigger eff
  
    //h_mt2_mc->SetBinError(bi, sqrt( value_err*value_err* hlt* hlt +  value*value* 0.05*0.05 ) ); //3% for the trigger eff
  
  
    //    h_mt2_mc->SetBinError(bi, sqrt( value_err*value_err* 0.92*0.92 +  value*value* 0.05*0.05 ) ); //tempororaa 5% for the trigger eff
  
    //change back to this:::  
    //h_mt2_mc->SetBinContent(bi, value * hlt );
    //h_mt2_mc->SetBinError(bi, sqrt( value_err*value_err* hlt* hlt +  value*value* 0.03*0.03 ) ); //3% for the trigger eff

    h_mt2_mc->SetBinContent(bi, value * hlt * lepSF );
    h_mt2_mc->SetBinError(bi, sqrt( value_err*value_err* hlt * lepSF* hlt * lepSF  + value*value*lepSF_err*lepSF_err*hlt*hlt +  value*value* lepSF*lepSF * 0.03*0.03 ) ); //3% for the trigger eff
  } 

  h_mt2_mc->Divide(g_mt2_mc);

  //Filling the single ratios
  for(int bi = 1; bi <= yieldBins ; bi++){
    double value_mc = h_mt2_mc->GetBinContent(bi);
    double err_mc = h_mt2_mc->GetBinError(bi);
    zllG_mc->get(thisRegion)->yield->SetBinContent( bi,   value_mc);
    zllG_mc->get(thisRegion)->yield_systUp->SetBinContent( bi, value_mc + err_mc);
    zllG_mc->get(thisRegion)->yield_systDown->SetBinContent( bi, value_mc - err_mc);
  }


  TGraphAsymmErrors* gr_data = zllG_data->get(thisRegion)->getGraph();
  gr_data->SetMarkerSize(1.4);
  gr_data->SetMarkerStyle(20);
  gr_data->SetLineColor(kBlack);
  gr_data->SetLineWidth(2);
  gr_data->SetMarkerColor(kBlack);

  *(zll_ratio) = *(zllG_data) / *(zllG_mc);

  TGraphAsymmErrors* gr_ratio = zll_ratio->get(thisRegion)->getGraph();
  gr_ratio->SetMarkerSize(1.4);
  gr_ratio->SetMarkerStyle(20);
  gr_ratio->SetLineColor(kBlack);
  gr_ratio->SetLineWidth(2);
  gr_ratio->SetMarkerColor(kBlack);

  for(int bi = 1; bi <= yieldBins ; bi++){
    Double_t x_tmp, p, p_errUp, p_errDown;	       
    gr_data->GetPoint( bi-1, x_tmp, p);
    if(p<0.005){
      gr_data->RemovePoint(bi-1);
      gr_ratio->RemovePoint(bi-1);
    }
  }

  h_mt2_mc->SetLineColor(38); h_mt2_mc->SetLineWidth(2);

  float yMax = 0.13;
  if( selection == "nBJets") yMax = 0.18;
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0.0, yMax);
  if(selection == "mt2"){
    h2_axes->SetXTitle("M_{T2} [GeV]");
  }else    if(selection == "ht"){
    h2_axes->SetXTitle("H_{T} [GeV]");
  }else    if(selection == "nJets"){
    h2_axes->SetXTitle("Jet Multiplicity");
  }else{
    h2_axes->SetXTitle("b Jet Multiplicity" );
  }
  h2_axes->SetYTitle("Z(ll) / #gamma Ratio");
  h2_axes->Draw();

  std::vector<std::string> niceNames2 = thisRegion.getNiceNames();

  for( unsigned i=0+1; i< niceNames2.size(); ++i ) {
    float yMaxText = 0.9-(float)i*0.05 +0.05;
    float yMinText = yMaxText - 0.05;
    TPaveText* regionText = new TPaveText( 0.18, yMinText, 0.55, yMaxText, "brNDC" );
    regionText->SetTextSize(0.04);
    //   regionText->SetTextFont(42);
    regionText->SetFillColor(0);
    regionText->SetTextAlign(11);
    regionText->AddText( niceNames2[i].c_str() );
    // regionText->Draw("same");
  }

  h_mt2_mc->SetFillColor(18);
  h_mt2_mc->SetFillStyle(3001);

  h_mt2_mc->Draw("E2 same");
  h_mt2_mc->Draw("line same");
  
  gr_data->Draw("p same");
  //    h_mt2->DrawClone("p same");
  // if(!onlyMC)
  //   gr_ratio->Draw("same P");
  MT2DrawTools::addLabels( (TCanvas*)pad1, lumi, "CMS Preliminary" );
  gPad->RedrawAxis();

  int i= 2 - int(onlyMC);
  TLegend* legend = new TLegend( 0.2, 0.92-(i)*0.06, 0.4, 0.92 );
  legend->SetTextSize(0.04);
  // legend->SetTextFont(42);
  legend->SetFillColor(0);
  if(!onlyMC){ 
    legend->AddEntry( gr_data ,"Data", "PL" );
  }
  legend->AddEntry( h_mt2_mc ,"Simulation", "L" );
  legend->Draw("same");



 
  gPad->RedrawAxis();

  canny->cd();

  TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , 0.3, 1.7 );
 

  if( !onlyMC){
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
    pad2->SetTopMargin(0.10);
    pad2->SetBottomMargin(0.1);
    pad2->Draw();
    pad2->cd();
 

    float grMax = h_mt2->GetBinLowEdge(size+1);
    float grMin = h_mt2->GetBinLowEdge(1);


    /*
    double sf = 0.902279;
    double sfErr = 0.102958;//0.109381;

    TF1* fSF=new TF1("fSF", "[0]", grMin, grMax);
    fSF->SetLineColor(kRed);
    fSF->SetParameter(0, sf );
    fSF->SetParError(0, sfErr );
    //   g_ratio->Fit(fSF, "0");

    
    double x[2]    ={grMin, grMax};
    double y[2]    ={sf, sf};

    double xerr[2] ={0., 0.};
    double yerr[2] ={sfErr, sfErr};
    TGraphErrors* SFFitBand = new TGraphErrors(2, x, y, xerr, yerr);
    SFFitBand->SetLineColor(0);
    SFFitBand->SetFillColor(kRed);
    SFFitBand->SetFillStyle(3244);
  
    TPaveText* fitText = MT2DrawTools::getFitText( fSF );
    */

    
    TF1* fSF = MT2DrawTools::getSFFit(gr_ratio, grMin, grMax  );
    TPaveText* fitText = MT2DrawTools::getFitText( fSF );
    TGraphErrors* SFFitBand =  MT2DrawTools::getSFFitBand(fSF, xMin, xMax);
    

    h2_axes_rat->SetYTitle("Data / MC");
    h2_axes_rat->GetXaxis()->SetTitleSize(0.2);
    h2_axes_rat->GetXaxis()->SetTitleOffset(4.8);
    h2_axes_rat->GetXaxis()->SetLabelSize(0.00);
    h2_axes_rat->GetXaxis()->SetTickLength(0.09);
    h2_axes_rat->GetYaxis()->SetNdivisions(5,5,0);
    h2_axes_rat->GetYaxis()->SetTitleSize(0.20);
    h2_axes_rat->GetYaxis()->SetTitleOffset(0.39);
    h2_axes_rat->GetYaxis()->SetLabelSize(0.14);
  
  
    h2_axes_rat->Draw();
 

    TLine* line = new TLine(xMin,1 , xMax, 1);
    line->SetLineColor(kBlack);
    line->Draw("same");
 
    SFFitBand->Draw("3,same");
    fSF->Draw("same");
    
    gr_ratio->Draw("same P");


    gPad->RedrawAxis();

    canny->cd();
    pad1->cd();
    fitText->Draw("same");
  }

  canny->SaveAs( Form("%s/%s_ratios.eps", fullPath.c_str(), saveName.c_str() ) );
  canny->SaveAs( Form("%s/%s_ratios.png", fullPath.c_str(), saveName.c_str() ) );
  canny->SaveAs( Form("%s/%s_ratios.pdf", fullPath.c_str(), saveName.c_str() ) );

  delete h_mt2;   delete h_mt2_mc;    delete g_mt2; delete g_mt2_mc; 
  delete h_lepSF; delete h_lepSF_err; delete h_HLT_weight;
  delete h2_axes; delete h2_axes_rat;
  delete gr_ratio;
  delete canny;

  delete h_top;
  delete g_Up; delete g_Down;

}



























void drawRatiosTopHisto(std::string fullPath, double *binss, unsigned int size,  std::string selection, MT2Analysis<MT2EstimateSyst>*  zll_ratio, MT2Analysis<MT2EstimateTree>* gamma_mc, MT2Analysis<MT2EstimateTree>* gamma_data, MT2Analysis<MT2EstimateSyst>* purity, MT2Analysis<MT2EstimateTree>*  zll_mc, MT2Analysis<MT2EstimateTree>* zll_data, MT2Analysis<MT2Estimate>* top, MT2Analysis<MT2EstimateSyst>* zll_yield, MT2Analysis<MT2EstimateSyst>*  zllG_data, MT2Analysis<MT2EstimateSyst>* zllG_mc, const MT2Region thisRegion, std::string cut, std::string cut_gamma, std::string cut_data, std::string cut_gamma_data, float lumi, std::string saveName, bool onlyMC, float scaleFactor, bool fullUncert , std::string topoCuts ){
 
  TH1F::AddDirectory(kTRUE);

  double bins[size+1]; 
  for(unsigned int i=0; i<= size ; ++i)  
    bins[i]=binss[i];
  float xMin = binss[0];   float xMax = binss[size];



  //THE TREES
  TTree *zll_data_tree =  zll_data->get(thisRegion)->tree;
  TTree *gamma_data_tree =  gamma_data->get(thisRegion)->tree;
  TTree *zll_mc_tree =  zll_mc->get(thisRegion)->tree;
  TTree *gamma_mc_tree =  gamma_mc->get(thisRegion)->tree;
  //TTree *top_tree =  top->get(thisRegion)->tree;
  TH1D* h_top = top->get(thisRegion)->yield;
  TGraphAsymmErrors* this_zinv_purity = purity->get(thisRegion)->getGraph();

  std::cout << "Read histograms and graphs" << std::endl;

  TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
  canny->cd();

  TPad* pad1 = 0;
  if( !onlyMC ) {
    pad1 = MT2DrawTools::getCanvasMainPad();
    pad1->Draw();
    pad1->cd();
  }

  // TH1D* h_top = new TH1D("h_top","", size , bins); h_top->Sumw2();  
  //  top_tree ->Project( "h_top" , selection.c_str(), Form("(%s)*%f", cut.c_str(), lumi ) );
  MT2DrawTools::addOverflowSingleHisto(h_top);

  TH1D* h_mt2 = new TH1D("h_mt2","", size , bins); h_mt2->Sumw2();
  zll_data_tree ->Project( "h_mt2" , selection.c_str(), cut_data.c_str() );
  MT2DrawTools::addOverflowSingleHisto(h_mt2);

  TH1D* g_mt2 = new TH1D("g_mt2","", size , bins); g_mt2->Sumw2();
  gamma_data_tree->Project( "g_mt2", selection.c_str(), cut_gamma_data.c_str() );
  MT2DrawTools::addOverflowSingleHisto(g_mt2);


  //GET THE UNCERTAINTIES
  double f = 0.92;
  double f_uncert = 0.1; // 0.08 from the fragmentation and then o+ ~5% for the mc closure
  if( !fullUncert ) f_uncert = 0.0;
  TH1D* g_Up = new TH1D("g_Up","", size , bins);
  TH1D* g_Down = new TH1D("g_Down","", size , bins);

  int nBinss =  h_mt2->GetNbinsX();
  for(int iBin = 1; iBin <= nBinss; iBin++){

    //Zll//////////////
    double value = h_mt2->GetBinContent(iBin);
    double top = h_top->GetBinContent(iBin);
    double top_err = h_top->GetBinError(iBin);

    h_mt2->SetBinContent(iBin, value - top*scaleFactor);
    //add 50% of the top bg as uncertainty to the estimate
    double zll_err =  sqrt(value + top_err*top_err);
    h_mt2->SetBinError(iBin,  zll_err  );

    //Photons//////////
    Double_t x_tmp, p, p_errUp, p_errDown;	       
    this_zinv_purity->GetPoint( iBin-1, x_tmp, p);
    p_errUp   = this_zinv_purity->GetErrorYhigh(iBin-1);
    p_errDown = this_zinv_purity->GetErrorYlow(iBin-1);
    if( !fullUncert ){
       p_errUp = 0.;       p_errDown = 0.;
    }

    double value_g = g_mt2->GetBinContent(iBin);
    g_mt2->SetBinContent(iBin, value_g  * p * f);

    double uncertUp = sqrt ( value_g*f*f*p*p + ( value_g*value_g*f*f *p_errUp*p_errUp) + (p*p*value_g*value_g* f_uncert*f_uncert) );
    double uncertDown = sqrt ( value_g*f*f*p*p + ( value_g*value_g*f*f *p_errDown*p_errDown) + (p*p*value_g*value_g* f_uncert*f_uncert) );
    if( !fullUncert ){
      uncertUp = sqrt ( value_g*f*f*p*p  );
      uncertDown = sqrt ( value_g*f*f*p*p );
    }
    g_Up->SetBinContent(iBin, uncertUp);
    g_Down->SetBinContent(iBin, uncertDown);

  }

  h_mt2->SetMarkerStyle(20);   h_mt2->SetMarkerSize(1.6); h_mt2->SetLineColor(kBlack);
  
 

  //Filling the YIELD
  zll_yield->get(thisRegion);
  int yieldBins = h_mt2->GetNbinsX();

  for(int bi = 1; bi <= yieldBins ; bi++){
    double value = h_mt2->GetBinContent(bi);
    if( value < 0.) value = 0.;
    double err = h_mt2->GetBinError(bi);
    zll_yield->get(thisRegion)->yield->SetBinContent( bi,   value);
    zll_yield->get(thisRegion)->yield_systUp->SetBinContent( bi,   value + err);
    zll_yield->get(thisRegion)->yield_systDown->SetBinContent( bi,   value - err);

    double value_g = g_mt2->GetBinContent(bi);
    double errUp = g_Up->GetBinContent(bi);
    double errDown = g_Down->GetBinContent(bi);
 
    //store the gamma's temporarily in zllG_data
    zllG_data->get(thisRegion)->yield->SetBinContent( bi,   value_g);
    zllG_data->get(thisRegion)->yield_systUp->SetBinContent( bi,   value_g + errUp);
    zllG_data->get(thisRegion)->yield_systDown->SetBinContent( bi,   value_g - errDown );
  }

  *(zllG_data) = *(zll_yield) / *(zllG_data);

  TH1D* h_mt2_mc = new TH1D("h_mt2_mc","", size  , bins);  h_mt2_mc->Sumw2();
  TH1D* g_mt2_mc = new TH1D("g_mt2_mc","", size  , bins);  g_mt2_mc->Sumw2();
  TH1D* h_HLT_weight = new TH1D("h_HLT_weight","", size  , bins);  h_HLT_weight->Sumw2();
  TH1D* h_lepSF = new TH1D("h_lepSF","", size  , bins);  h_lepSF->Sumw2();
  TH1D* h_lepSF_err = new TH1D("h_lepSF_err","", size  , bins);  h_lepSF_err->Sumw2();

  zll_mc_tree ->Project( "h_lepSF", selection.c_str(), Form("(%s)*%f *%s", cut.c_str(), lumi, "(weight_lep0 * weight_lep1)" ) );
  zll_mc_tree ->Project( "h_lepSF_err", selection.c_str(), Form("(%s)*%f *%s", cut.c_str(), lumi, "(weight_lep_err)" ) );
  zll_mc_tree ->Project( "h_HLT_weight", selection.c_str(), Form("(%s)*%f *%s", cut.c_str(), lumi, "(HLT_weight)" ) );
  zll_mc_tree ->Project( "h_mt2_mc", selection.c_str(), Form("(%s)*%f", cut.c_str(), lumi) );
  gamma_mc_tree->Project( "g_mt2_mc", selection.c_str(),  Form("(%s)*%f", cut_gamma.c_str(), lumi)  );

  MT2DrawTools::addOverflowSingleHisto( h_mt2_mc );
  MT2DrawTools::addOverflowSingleHisto( g_mt2_mc );
  MT2DrawTools::addOverflowSingleHisto( h_lepSF );
  MT2DrawTools::addOverflowSingleHisto( h_lepSF_err );
  MT2DrawTools::addOverflowSingleHisto( h_HLT_weight );

 
  //Filling the single ratios
  for(int bi = 1; bi <= yieldBins ; bi++){
    double value = h_mt2_mc->GetBinContent(bi);
    double value_err = h_mt2_mc->GetBinError(bi);
    double hlt = h_HLT_weight->GetBinContent(bi) /value ;
    double lepSF = h_lepSF->GetBinContent(bi)/value;
    double lepSF_err = h_lepSF_err->GetBinContent(bi)/value;
 
    // h_mt2_mc->SetBinContent(bi, value * hlt * lepSF );
    // h_mt2_mc->SetBinError(bi, sqrt( value_err*value_err* hlt * lepSF* hlt * lepSF  + value*value*lepSF_err*lepSF_err*hlt*hlt +  value*value* lepSF*lepSF * 0.03*0.03 ) ); //3% for the trigger eff
  } 

  h_mt2_mc->Divide(g_mt2_mc);

  //Filling the single ratios
  for(int bi = 1; bi <= yieldBins ; bi++){
    double value_mc = h_mt2_mc->GetBinContent(bi);
    double err_mc = h_mt2_mc->GetBinError(bi);
    zllG_mc->get(thisRegion)->yield->SetBinContent( bi,   value_mc);
    zllG_mc->get(thisRegion)->yield_systUp->SetBinContent( bi, value_mc + err_mc);
    zllG_mc->get(thisRegion)->yield_systDown->SetBinContent( bi, value_mc - err_mc);
  }


  TGraphAsymmErrors* gr_data = zllG_data->get(thisRegion)->getGraph();
  gr_data->SetMarkerSize(1.4);
  gr_data->SetMarkerStyle(20);
  gr_data->SetLineColor(kBlack);
  gr_data->SetLineWidth(2);
  gr_data->SetMarkerColor(kBlack);

  *(zll_ratio) = *(zllG_data) / *(zllG_mc);

  TGraphAsymmErrors* gr_ratio = zll_ratio->get(thisRegion)->getGraph();
  gr_ratio->SetMarkerSize(1.4);
  gr_ratio->SetMarkerStyle(20);
  gr_ratio->SetLineColor(kBlack);
  gr_ratio->SetLineWidth(2);
  gr_ratio->SetMarkerColor(kBlack);

  for(int bi = 1; bi <= yieldBins ; bi++){
    Double_t x_tmp, p, p_errUp, p_errDown;	       
    gr_data->GetPoint( bi-1, x_tmp, p);
    if(p<0.005){
      gr_data->RemovePoint(bi-1);
      gr_ratio->RemovePoint(bi-1);
    }
  }

  h_mt2_mc->SetLineColor(38); h_mt2_mc->SetLineWidth(2);

  float yMax = 0.13;
  if( selection == "nBJets") yMax = 0.18;
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0.0, yMax);
  if(selection == "mt2"){
    h2_axes->SetXTitle("M_{T2} [GeV]");
  }else    if(selection == "ht"){
    h2_axes->SetXTitle("H_{T} [GeV]");
  }else    if(selection == "nJets"){
    h2_axes->SetXTitle("Jet Multiplicity");
  }else{
    h2_axes->SetXTitle("b Jet Multiplicity" );
  }
  h2_axes->SetYTitle("Z(ll) / #gamma Ratio");
  h2_axes->Draw();

  std::vector<std::string> niceNames2 = thisRegion.getNiceNames();

  for( unsigned i=0+1; i< niceNames2.size(); ++i ) {
    float yMaxText = 0.9-(float)i*0.05 +0.05;
    float yMinText = yMaxText - 0.05;
    TPaveText* regionText = new TPaveText( 0.18, yMinText, 0.55, yMaxText, "brNDC" );
    regionText->SetTextSize(0.04);
    //   regionText->SetTextFont(42);
    regionText->SetFillColor(0);
    regionText->SetTextAlign(11);
    regionText->AddText( niceNames2[i].c_str() );
    // regionText->Draw("same");
  }

  h_mt2_mc->SetFillColor(18);
  h_mt2_mc->SetFillStyle(3001);

  h_mt2_mc->Draw("E2 same");
  h_mt2_mc->Draw("line same");
  
  gr_data->Draw("p same");
  //    h_mt2->DrawClone("p same");
  // if(!onlyMC)
  //   gr_ratio->Draw("same P");
  MT2DrawTools::addLabels( (TCanvas*)pad1, lumi, "CMS" );
  gPad->RedrawAxis();

  int i= 2 - int(onlyMC);
  TLegend* legend = new TLegend( 0.2, 0.92-(i)*0.06, 0.4, 0.92 );
  legend->SetTextSize(0.04);
  // legend->SetTextFont(42);
  legend->SetFillColor(0);
  if(!onlyMC){ 
    legend->AddEntry( gr_data ,"Data", "P" );
  }
  legend->AddEntry( h_mt2_mc ,"Simulation", "L" );
  legend->Draw("same");



 
  gPad->RedrawAxis();

  canny->cd();

  float rangeYLow = 0.3;
  if( selection == "nBJets") rangeYLow = 0;

  TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , rangeYLow, 1.7 );
 

  if( !onlyMC){
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
    pad2->SetTopMargin(0.10);
    pad2->SetBottomMargin(0.1);
    pad2->Draw();
    pad2->cd();
 
    float grMax = h_mt2->GetBinLowEdge(size+1);
    float grMin = h_mt2->GetBinLowEdge(1);
    TF1* fSF = MT2DrawTools::getSFFit(gr_ratio, grMin, grMax  );
  
    TPaveText* fitText = MT2DrawTools::getFitText( fSF );
    TGraphErrors* SFFitBand =  MT2DrawTools::getSFFitBand(fSF, xMin, xMax);
    

    h2_axes_rat->SetYTitle("Data / MC");
    h2_axes_rat->GetXaxis()->SetTitleSize(0.2);
    h2_axes_rat->GetXaxis()->SetTitleOffset(4.8);
    h2_axes_rat->GetXaxis()->SetLabelSize(0.00);
    h2_axes_rat->GetXaxis()->SetTickLength(0.09);
    h2_axes_rat->GetYaxis()->SetNdivisions(5,5,0);
    h2_axes_rat->GetYaxis()->SetTitleSize(0.20);
    h2_axes_rat->GetYaxis()->SetTitleOffset(0.39);
    h2_axes_rat->GetYaxis()->SetLabelSize(0.14);
  
  
    h2_axes_rat->Draw();
 

    TLine* line = new TLine(xMin,1 , xMax, 1);
    line->SetLineColor(kBlack);
    line->Draw("same");
 
    SFFitBand->Draw("3,same");
    fSF->Draw("same");

    gr_ratio->Draw("same P");


    gPad->RedrawAxis();

    canny->cd();
    pad1->cd();
    fitText->Draw("same");
  }

  canny->SaveAs( Form("%s/%s_ratios.eps", fullPath.c_str(), saveName.c_str() ) );
  canny->SaveAs( Form("%s/%s_ratios.png", fullPath.c_str(), saveName.c_str() ) );
  canny->SaveAs( Form("%s/%s_ratios.pdf", fullPath.c_str(), saveName.c_str() ) );

  delete h_mt2;   delete h_mt2_mc;    delete g_mt2; delete g_mt2_mc; 
  delete h_lepSF; delete h_lepSF_err; delete h_HLT_weight;
  delete h2_axes; delete h2_axes_rat;
  delete gr_ratio;
  delete canny;

  delete h_top;
  delete g_Up; delete g_Down;

}

















