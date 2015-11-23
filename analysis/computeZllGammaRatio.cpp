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

bool HFveto = false;

double lumiErr = 0.12;

void drawRatios(std::string fullPath, double *binss, unsigned int size,  std::string zll_sel, MT2Analysis<MT2EstimateSyst>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma_mc, MT2Analysis<MT2EstimateTree>*  gamma_data, MT2Analysis<MT2EstimateSyst>*  purity, MT2Analysis<MT2EstimateTree>*  zll_mc,MT2Analysis<MT2EstimateTree>*  zll_data,   MT2Analysis<MT2EstimateTree>*  top, MT2Analysis<MT2EstimateSyst>*  zll_yield, MT2Analysis<MT2EstimateSyst>*  zllG_data,  MT2Analysis<MT2EstimateSyst>*  zllG_mc, const MT2Region thisRegion, std::string cut,  std::string cut_gamma, std::string cut_data, std::string cut_gamma_data, float lumi, std::string saveName, bool onlyMC, float top_SF, bool fullUncert , std::string topoCuts="" );



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

  std::string regionsSet = cfg.regionsSet();
  std::cout << "-> Using regions: " << regionsSet << std::endl;

  std::string samples = cfg.mcSamples();

  std::string outputdir = cfg.getEventYieldDir() + "/zllGammaRatio";
  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string gammaControlRegionDir = cfg.getEventYieldDir() + "/gammaControlRegion";
  std::string ZllDir = cfg.getEventYieldDir() + "/zllControlRegion";

  gStyle->SetOptTitle(0);
  MT2DrawTools::setStyle();

  
  ifstream SF_file;
  SF_file.open(Form("%s/plotsDataMCscaling/scaleFactorOF.txt", cfg.getEventYieldDir().c_str() ) );
  float scaleFactor;
  SF_file >> scaleFactor;
  std::cout<< "Scale Factor = "  << scaleFactor << std::endl;
  scaleFactor=1;
  std::cout<< "Scale Factor = "  << scaleFactor << std::endl;
  


 
  MT2Analysis<MT2EstimateTree>* gamma_mc = MT2Analysis<MT2EstimateTree>::readFromFile(gammaControlRegionDir + "/mc.root", "gammaCRtree");
  MT2Analysis<MT2EstimateTree>*  gamma_data = MT2Analysis<MT2EstimateTree>    ::readFromFile( gammaControlRegionDir + "/data.root", "gammaCRtree");
 
 MT2Analysis<MT2EstimateSyst>* purity = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_data.root", "purity");


  MT2Analysis<MT2EstimateSyst>* purity_ht = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_ht_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_njets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_njets_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_nbjets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_nbjets_data.root", "purity");
 
 MT2Analysis<MT2EstimateSyst>* purity_mono_nbjets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_mono_nbjets_data.root", "purity");
 


  MT2Analysis<MT2EstimateSyst>* purity_incl_ht = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_incl_ht_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_incl_njets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_incl_njets_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_incl_nbjets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_incl_nbjets_data.root", "purity");


  if(onlyMC){
    gamma_data = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc.root", gammaControlRegionDir.c_str()) , "gammaCRtree");
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

  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str()  ), "QCD");
  MT2Analysis<MT2EstimateTree>* top = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "Top");
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "WJets");

  wjets->setFullName("W+jets");

  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
  bgYields.push_back( qcd );
  bgYields.push_back( wjets );
  bgYields.push_back( top );
  

  MT2Analysis<MT2EstimateSyst>* zll_ratio = new MT2Analysis<MT2EstimateSyst>( "zll_ratio", regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zll_yield = new MT2Analysis<MT2EstimateSyst>( "zll_yield", regionsSet.c_str() );

  //YIELDS
  MT2Analysis<MT2EstimateSyst>* zll_pt = new MT2Analysis<MT2EstimateSyst>( "zll_pt", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zll_mt2 = new MT2Analysis<MT2EstimateSyst>( "zll_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zll_ht = new MT2Analysis<MT2EstimateSyst>( "zll_ht", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zll_nJets = new MT2Analysis<MT2EstimateSyst>( "zll_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zll_nBJets = new MT2Analysis<MT2EstimateSyst>("zll_nBJets",regionsSet.c_str());
  MT2Analysis<MT2EstimateSyst>* zll_mono_nBJets = new MT2Analysis<MT2EstimateSyst>("zll_mono_nBJets",regionsSet.c_str());

  MT2Analysis<MT2EstimateSyst>* zll_incl_nJets = new MT2Analysis<MT2EstimateSyst>("zll_incl_nJets",regionsSet.c_str());

  //DATA RATIOS
  MT2Analysis<MT2EstimateSyst>* zllG_data_mt2 = new MT2Analysis<MT2EstimateSyst>( "zllG_data_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_data_ht = new MT2Analysis<MT2EstimateSyst>( "zllG_data_ht", regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_data_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_data_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_data_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_data_nBJets",regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zllG_data_mono_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_data_mono_nBJets",regionsSet.c_str() );

  MT2Analysis<MT2EstimateSyst>* zllG_data_incl_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_data_incl_nJets",regionsSet.c_str()); 

  //MC RATIOS
  MT2Analysis<MT2EstimateSyst>* zllG_mc_mt2 = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_mc_ht = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_ht", regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_mc_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_mc_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_mc_nBJets",regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zllG_mc_mono_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_mc_mono_nBJets",regionsSet.c_str() );

  MT2Analysis<MT2EstimateSyst>* zllG_mc_incl_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_incl_nJets",regionsSet.c_str()); 

  //RATIOS
  MT2Analysis<MT2EstimateSyst>* zllG_pt = new MT2Analysis<MT2EstimateSyst>( "zllG_pt", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_mt2 = new MT2Analysis<MT2EstimateSyst>( "zllG_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_ht = new MT2Analysis<MT2EstimateSyst>( "zllG_ht", regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_nBJets",regionsSet.c_str() );
  MT2Analysis<MT2EstimateSyst>* zllG_mono_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_mono_nBJets",regionsSet.c_str() );
  
 MT2Analysis<MT2EstimateSyst>* zllG_incl_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_incl_nJets",regionsSet.c_str());

  std::set<MT2Region> MT2Regions = zll_mc->getRegions();

 
  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  TH1F::AddDirectory(kTRUE);


  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
    std::vector<std::string> niceNames = thisRegion.getNiceNames();
    float lumi = cfg.lumi();    

    double bins_nJets[] = {2,4,7,12};
    double bins_nBJets[] = {0,1,2,3,6}; 
    double bins_mt2[] ={200,300,400,500, 600, 1500 };
    // if( !onlyMC){
    //  bins_mt2 = {200,300,400,500, 600, 1500 };
    // }else{
    //  bins_mt2 = {200,300,400,500, 600, 800, 1000, 1500 };
    //    }
    //    float bins_ht[] =  {200,250,300,350,400,450,500,550,600,700,800,900,1000,1500,2000};
    double bins_ht[] =  {200,450,575,1000,1500,3000};
  
    int size_mt2 = sizeof(bins_mt2)/sizeof(double)-1;
    int size_ht = sizeof(bins_ht)/sizeof(double)-1;
    int size_nJets = sizeof(bins_nJets)/sizeof(double)-1;
    int size_nBJets = sizeof(bins_nBJets)/sizeof(double)-1;
  
    double bins_incl_nJets[] = {1,2,4,7,12};
    int size_incl_nJets = sizeof(bins_incl_nJets)/sizeof(double)-1;
 

    std::string cut_incl = "weight*(abs(Z_mass-91.19)<10 && met>200 && mt2>200 && ht>200 && nJets>0 )";
    std::string cut_incl_gamma = "weight*(prompt==2 && iso<2.5  && met>200&&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23";
    std::string cut_incl_data = "weight*(abs(Z_mass-91.19)<10 && met>200 &&  mt2>200 && ht>200 && nJets>0 )";
    std::string cut_incl_gamma_data = "weight*( iso<2.5 && ptGamma>180  && met>200 && nJets>0 && mt2>200 && ht>200 )";
  

    std::string cut = "weight*(abs(Z_mass-91.19)<10 && met>200 && mt2>200 && ht>200 && nJets>1 )";
    std::string cut_el = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200 && mt2>200 && nJets>1 && Z_lepId==11 )";
    std::string cut_mu = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200&& mt2>200 && nJets>1 && Z_lepId==13 )";
    std::string cut_mono = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200 && nJets==1 )";

    std::string cut_mono_el = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200 && nJets==1 && Z_lepId==11 )";
    std::string cut_mono_mu = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200 && nJets==1 && Z_lepId==13 )";

    std::string cut_gamma = "weight*(prompt==2 && iso<2.5  && met>200&&  ptGamma>180 && nJets>1 && mt2>200 && ht>200 )*1.23";
    std::string cut_gamma_mono = "weight*( prompt==2 && iso<2.5 && met>200 && ptGamma>180 && nJets==1 && ht>200 )*1.23";
 
    //f & purity later in the function
    std::string cut_data = "weight*(abs(Z_mass-91.19)<10 && met>200 &&  mt2>200 && ht>200 && nJets>1 )";
    std::string cut_el_data = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200 && mt2>200 && nJets>1&& Z_lepId==11 )";
    std::string cut_mu_data = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200&& mt2>200 && nJets>1&& Z_lepId==13 )";
    std::string cut_mono_data = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200 && nJets==1 )";

    std::string cut_mono_data_el = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200 && nJets==1 && Z_lepId==11)";
    std::string cut_mono_data_mu = "weight*(abs(Z_mass-91.19)<10 && met>200 && ht>200 && nJets==1 && Z_lepId==13)";

    std::string cut_gamma_data = "weight*( iso<2.5 && ptGamma>180  && met>200 && nJets>1 && mt2>200 && ht>200 )";
    std::string cut_gamma_mono_data = "weight*( iso<2.5 &&  ptGamma>180  && met>200&& nJets==1 && ht>200 )";

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
    double bins_mono_nBJets[] = {0,1,2}; 
    int size_mono_nBJets = sizeof(bins_mono_nBJets)/sizeof(double)-1;
    MT2EstimateSyst::rebinYields( zll_incl_nJets,  size_incl_nJets, bins_incl_nJets);
    MT2EstimateSyst::rebinYields( zllG_incl_nJets,  size_incl_nJets, bins_incl_nJets);
    MT2EstimateSyst::rebinYields( zllG_data_incl_nJets,  size_incl_nJets, bins_incl_nJets);
    MT2EstimateSyst::rebinYields( zllG_mc_incl_nJets,  size_incl_nJets, bins_incl_nJets);
  

    MT2EstimateSyst::rebinYields( zll_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zll_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zll_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zll_nBJets,  size_nBJets, bins_nBJets);
    MT2EstimateSyst::rebinYields( zll_mono_nBJets,  size_mono_nBJets, bins_mono_nBJets);
 
    MT2EstimateSyst::rebinYields( zllG_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zllG_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zllG_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zllG_nBJets,  size_nBJets, bins_nBJets);
    MT2EstimateSyst::rebinYields( zllG_mono_nBJets,  size_mono_nBJets, bins_mono_nBJets);
  
    MT2EstimateSyst::rebinYields( zllG_data_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zllG_data_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zllG_data_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zllG_data_nBJets,  size_nBJets, bins_nBJets);
    MT2EstimateSyst::rebinYields( zllG_data_mono_nBJets,  size_mono_nBJets, bins_mono_nBJets);
 
    MT2EstimateSyst::rebinYields( zllG_mc_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zllG_mc_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zllG_mc_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zllG_mc_nBJets,  size_nBJets, bins_nBJets);
    MT2EstimateSyst::rebinYields( zllG_mc_mono_nBJets,  size_mono_nBJets, bins_mono_nBJets);
 

    //draw ratio also fills the ratio and yield estimates
    //outputdir, bins, nbins, var to project, ratio estimate, gamma mc, gamma data, purity gamma, zll mc, zll data, yield estimate, region, cut zll, cut gamma, cut zll data, cut gamma data, lumi, name, flag , topo region);
    
    //  drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data, top, zll_mt2 , zllG_data_mt2, zllG_mc_mt2 , thisRegion, cut_el, cut_gamma, cut_el_data,  cut_gamma_data,  lumi , "mt2_el" , onlyMC, scaleFactor ,"#geq2j, #geq0b" );
    //  drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data, top, zll_mt2, zllG_data_mt2, zllG_mc_mt2 , thisRegion, cut_mu,  cut_gamma, cut_mu_data, cut_gamma_data, lumi , "mt2_mu" , onlyMC, scaleFactor ,"#geq2j, #geq0b");
    // drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data, top, zll_mt2, zllG_data_mt2, zllG_mc_mt2 ,thisRegion, cut, cut_gamma, cut_data,cut_gamma_data, lumi , "mt2" , onlyMC, scaleFactor ,"#geq2j, #geq0b");
   
    //  drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_el, cut_gamma, cut_el_data,  cut_gamma_data, lumi, "ht_el" , onlyMC, scaleFactor ,"#geq2j, #geq0b");
    // drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_mu, cut_gamma, cut_mu_data,  cut_gamma_data, lumi,"ht_mu" , onlyMC, scaleFactor ,"#geq2j, #geq0b"); 
   
    //incl_//////////////////////
    drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_incl_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"ht_incl" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_incl_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi,"ht_incl_noPFU" , onlyMC, scaleFactor, 0 ,"#geq1j, #geq0b");
    

    //  drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut, cut_gamma, cut_data,  cut_gamma_data, lumi,"ht" , onlyMC , scaleFactor,"#geq2j, #geq0b");
    
    //  drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion, cut_el, cut_gamma, cut_el_data,  cut_gamma_data, lumi, "nJets_el" , onlyMC , scaleFactor,"#geq2j, #geq0b");
    //  drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion, cut_mu, cut_gamma, cut_mu_data,  cut_gamma_data, lumi, "nJets_mu" , onlyMC, scaleFactor ,"#geq2j, #geq0b");
    //incl//////////////////////////////////////
    drawRatios( outputdir, bins_incl_nJets, size_incl_nJets , "nJets",   zllG_incl_nJets,   gamma_mc, gamma_data, purity_incl_njets,  zll_mc, zll_data,top,  zll_incl_nJets, zllG_data_incl_nJets, zllG_mc_incl_nJets ,  thisRegion,cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi, "nJets_incl" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    drawRatios( outputdir, bins_incl_nJets, size_incl_nJets , "nJets",   zllG_incl_nJets,   gamma_mc, gamma_data, purity_incl_njets,  zll_mc, zll_data,top,  zll_incl_nJets, zllG_data_incl_nJets, zllG_mc_incl_nJets ,  thisRegion,cut_incl, cut_incl_gamma, cut_incl_data,  cut_incl_gamma_data, lumi, "nJets_incl_noPFU" , onlyMC, scaleFactor, 0 ,"#geq1j, #geq0b");

    //drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion,cut, cut_gamma, cut_data,  cut_gamma_data, lumi, "nJets" , onlyMC , scaleFactor,"#geq2j, #geq0b");

    // drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_el, cut_gamma ,  cut_el_data,  cut_gamma_data, lumi, "nBJets_el" , onlyMC , scaleFactor,"#geq2j, #geq0b");
    //  drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_mu, cut_gamma ,  cut_mu_data,  cut_gamma_data, lumi, "nBJets_mu" , onlyMC , scaleFactor,"#geq2j, #geq0b");
    //incl///////////////////
    drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_incl_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_incl, cut_incl_gamma ,  cut_incl_data,  cut_incl_gamma_data, lumi, "nBJets_incl" , onlyMC, scaleFactor, 1 ,"#geq1j, #geq0b");
    drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_incl_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_incl, cut_incl_gamma ,  cut_incl_data,  cut_incl_gamma_data, lumi, "nBJets_incl_noPFU" , onlyMC, scaleFactor, 0 ,"#geq1j, #geq0b");

    // drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut, cut_gamma ,  cut_data,  cut_gamma_data, lumi, "nBJets" , onlyMC , scaleFactor,"#geq2j, #geq0b");


   //MONOJET
    //   drawRatios( outputdir, bins_mono_nBJets, size_mono_nBJets , "nBJets",  zllG_mono_nBJets,   gamma_mc, gamma_data, purity_mono_nbjets,  zll_mc, zll_data, top, zll_mono_nBJets, zllG_data_mono_nBJets, zllG_mc_mono_nBJets , thisRegion, cut_mono_el, cut_gamma_mono ,  cut_mono_data_el,  cut_gamma_mono_data, lumi, "mono_nBJets_el" , onlyMC , scaleFactor,"#geq1j, #geq0b");
    // drawRatios( outputdir, bins_mono_nBJets, size_mono_nBJets , "nBJets",  zllG_mono_nBJets,   gamma_mc, gamma_data, purity_mono_nbjets,  zll_mc, zll_data, top, zll_mono_nBJets, zllG_data_mono_nBJets, zllG_mc_mono_nBJets , thisRegion, cut_mono_mu, cut_gamma_mono ,  cut_mono_data_mu,  cut_gamma_mono_data, lumi, "mono_nBJets_mu" , onlyMC , scaleFactor,"#geq1j, #geq0b");

    // drawRatios( outputdir, bins_mono_nBJets, size_mono_nBJets , "nBJets",  zllG_mono_nBJets,   gamma_mc, gamma_data, purity_mono_nbjets,  zll_mc, zll_data, top, zll_mono_nBJets, zllG_data_mono_nBJets, zllG_mc_mono_nBJets , thisRegion, cut_mono, cut_gamma_mono ,  cut_mono_data,  cut_gamma_mono_data, lumi, "mono_nBJets" , onlyMC, scaleFactor ,"#geq1j, #geq0b");

    

  }
  
 
  std::string outFile = outputdir + "/zll_ratio.root";
  zllG_mt2->writeToFile(outFile);
  zllG_ht->addToFile( outFile );
  zllG_nJets->addToFile( outFile );
  zllG_nBJets->addToFile( outFile );
  zllG_mono_nBJets->addToFile( outFile );

  std::string outFile_yield = outputdir + "/zll_yield.root";
  zll_mt2->writeToFile(outFile_yield);
  zll_ht->addToFile( outFile_yield );
  zll_nJets->addToFile( outFile_yield );
  zll_nBJets->addToFile( outFile_yield );
  zll_mono_nBJets->addToFile( outFile_yield );

  std::string outFile_data = outputdir + "/zllG_data_ratio.root";
  zllG_data_mt2->writeToFile(outFile_data);
  zllG_data_ht->addToFile( outFile_data );
  zllG_data_nJets->addToFile( outFile_data );
  zllG_data_nBJets->addToFile( outFile_data );
  zllG_data_mono_nBJets->addToFile( outFile_data );

 std::string outFile_mc = outputdir + "/zllG_mc_ratio.root";
  zllG_mc_mt2->writeToFile(outFile_mc);
  zllG_mc_ht->addToFile( outFile_mc );
  zllG_mc_nJets->addToFile( outFile_mc );
  zllG_mc_nBJets->addToFile( outFile_mc );
  zllG_mc_mono_nBJets->addToFile( outFile_mc );


  return 0;

}


void drawRatios(std::string fullPath, double *binss, unsigned int size,  std::string zll_sel, MT2Analysis<MT2EstimateSyst>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma_mc, MT2Analysis<MT2EstimateTree>*  gamma_data, MT2Analysis<MT2EstimateSyst>*  purity, MT2Analysis<MT2EstimateTree>*  zll_mc,MT2Analysis<MT2EstimateTree>*  zll_data, MT2Analysis<MT2EstimateTree>*  top, MT2Analysis<MT2EstimateSyst>*  zll_yield, MT2Analysis<MT2EstimateSyst>*  zllG_data,  MT2Analysis<MT2EstimateSyst>*  zllG_mc, const MT2Region thisRegion, std::string cut, std::string cut_gamma, std::string cut_data, std::string cut_gamma_data, float lumi, std::string saveName, bool onlyMC, float scaleFactor, bool fullUncert , std::string topoCuts ){
 
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
  if(onlyMC)
    labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  TH1F::AddDirectory(kTRUE);

  double bins[size+1]; 
  for(unsigned int i=0; i<= size ; ++i)  
    bins[i]=binss[i];
  float xMin = binss[0];   float xMax = binss[size];

  std::string gamma_sel = zll_sel;
 
  //THE TREES
  TTree *zll_data_tree =  zll_data->get(thisRegion)->tree;
  TTree *gamma_data_tree =  gamma_data->get(thisRegion)->tree;
  TTree *zll_mc_tree =  zll_mc->get(thisRegion)->tree;
  TTree *gamma_mc_tree =  gamma_mc->get(thisRegion)->tree;
  TTree *top_tree =  top->get(thisRegion)->tree;
  TGraphAsymmErrors* this_zinv_purity = purity->get(thisRegion)->getGraph();


  //Getting the weighted center of the bin
  float meanX[size+1];
  float meanX_err[size+1];
  for(unsigned int i = 0; i < size; ++i){
    TH1D* histo_x= new TH1D("histo_x","",100, bins[i], bins[i+1] );
    zll_mc_tree->Project( "histo_x",zll_sel.c_str(), cut.c_str());
    histo_x->GetMean();
    meanX[i] =  histo_x->GetMean();
    meanX_err[i] =  histo_x->GetMeanError();
    delete histo_x;
  }


  TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
  canny->cd();

  TPad* pad1 = 0;
  if( !onlyMC ) {
    pad1 = MT2DrawTools::getCanvasMainPad();
    pad1->Draw();
    pad1->cd();
  }

  TH1D* h_top = new TH1D("h_top","", size , bins); h_top->Sumw2();
  TH1D* h_mt2 = new TH1D("h_mt2","", size , bins); h_mt2->Sumw2();
  TH1D* g_mt2 = new TH1D("g_mt2","", size , bins); g_mt2->Sumw2();
  
  top_tree ->Project( "h_top" , zll_sel.c_str(), cut_data.c_str() );
  zll_data_tree ->Project( "h_mt2" , zll_sel.c_str(), cut_data.c_str() );
  gamma_data_tree->Project( "g_mt2", gamma_sel.c_str(), cut_gamma_data.c_str() );
  h_top->Scale( lumi );
  h_top->Scale( scaleFactor ); //to come soon, hopefully

  h_top->SetBinContent(size, h_top->GetBinContent(size) + h_top->GetBinContent(size+1));//adding overflow
  h_mt2->SetBinContent(size, h_mt2->GetBinContent(size) + h_mt2->GetBinContent(size+1));//adding overflow
  g_mt2->SetBinContent(size, g_mt2->GetBinContent(size) + g_mt2->GetBinContent(size+1));
  h_mt2->SetBinContent(size+1,0);
  h_top->SetBinContent(size+1,0);
  g_mt2->SetBinContent(size+1,0);


  //GET THE UNCERTAINTIES
  double f = 0.92;
  double f_uncert = 0.1; // 0.08 from the fragmentation and then o+ ~5% for the mc closure
  if( !fullUncert ) f_uncert = 0.0;
  TH1D* g_Up = new TH1D("g_Up","", size , bins);
  TH1D* g_Down = new TH1D("g_Down","", size , bins);

  int nBinss =  h_mt2->GetNbinsX();
  for(int binnie = 1; binnie <= nBinss; binnie++){
    //Zll//////////////
    double value = h_mt2->GetBinContent(binnie);
    //value -= 1.;
    double top = h_top->GetBinContent(binnie);

    std::cout << "Zll " << value << "  top " << top << std::endl;
 
    h_mt2->SetBinContent(binnie, value - top);
    //add 50% of the top bg as uncertainty to the estimate
    h_mt2->SetBinError(binnie,  sqrt(value + top*top*0.25)  );

    //Photons//////////
    Double_t x_tmp, p, p_errUp, p_errDown;	       
    this_zinv_purity->GetPoint( binnie-1, x_tmp, p);
    p_errUp   = this_zinv_purity->GetErrorYhigh(binnie-1);
    p_errDown = this_zinv_purity->GetErrorYlow(binnie-1);
    if( !fullUncert ){
       p_errUp = 0.;       p_errDown = 0.;
    }
    std::cout << "Purity = " << p << std::endl;
    std::cout << "Purity up = " << p_errUp << std::endl;
    std::cout << "Purity down = " << p_errDown << std::endl;

    double value_g = g_mt2->GetBinContent(binnie);

    g_mt2->SetBinContent(binnie, value_g  * p * f);

    double uncertUp = sqrt ( value_g*f*f*p*p + ( value_g*value_g*f*f *p_errUp*p_errUp) + (p*p*value_g*value_g* f_uncert*f_uncert) );
    double uncertDown = sqrt ( value_g*f*f*p*p + ( value_g*value_g*f*f *p_errDown*p_errDown) + (p*p*value_g*value_g* f_uncert*f_uncert) );
    if( !fullUncert ){
      uncertUp = sqrt ( value_g*f*f*p*p  );
      uncertDown = sqrt ( value_g*f*f*p*p );
    }
    g_Up->SetBinContent(binnie, uncertUp);
    g_Down->SetBinContent(binnie, uncertDown);

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

    std::cout << "ZLL uncert " << value  << " +- " <<err << std::endl;
    std::cout << "Gamma uncert " << value_g << " + " << errUp <<  " - " << errDown << std::endl;
  }

  *(zllG_data) = *(zll_yield) / *(zllG_data);

   for(int bi = 1; bi <= yieldBins ; bi++){
     std::cout << "Data ratio = " <<   zllG_data->get(thisRegion)->yield->GetBinContent(bi) << " + " << zllG_data->get(thisRegion)->yield_systUp->GetBinContent(bi) << " - " << zllG_data->get(thisRegion)->yield_systDown->GetBinContent(bi) << std::endl;
   }


  TH1D* h_mt2_mc = new TH1D("h_mt2_mc","", size  , bins);
  TH1D* g_mt2_mc = new TH1D("g_mt2_mc","", size  , bins);
  h_mt2_mc->Sumw2(); g_mt2_mc->Sumw2();

  zll_mc_tree ->Project( "h_mt2_mc" , zll_sel.c_str(), Form("(%s)*%f",cut.c_str(), lumi) );
  gamma_mc_tree->Project( "g_mt2_mc", gamma_sel.c_str(),  Form("(%s)*%f",cut_gamma.c_str(), lumi)  );

  h_mt2_mc->SetBinContent(size, h_mt2_mc->GetBinContent(size) +h_mt2_mc->GetBinContent(size+1));
  g_mt2_mc->SetBinContent(size, g_mt2_mc->GetBinContent(size) +g_mt2_mc->GetBinContent(size+1));
  h_mt2_mc->SetBinContent(size+1, 0);
  g_mt2_mc->SetBinContent(size+1, 0);


  TH1D* h_clone = (TH1D*)h_mt2->Clone("h_clone");
  h_clone->Divide(h_mt2_mc);
  TH1D* g_clone = (TH1D*)g_mt2->Clone("g_clone");
  g_clone->Divide(g_mt2_mc);

  TH2D* h2_axes_yield = new TH2D("axes_yield", "", 10, 0, 1500, 10, 0.1, 3);
  h2_axes_yield->SetXTitle("M_{T2} [GeV]"); h2_axes_yield->SetYTitle("Zll Data/MC ratio");h2_axes_yield->Draw();  
  h_clone->Draw("p same");    labelTop->Draw("same");   gPad->RedrawAxis();    //  gPad->SetLogy();
  canny->SaveAs( Form("%s/yield_zll_%s_%s.eps", fullPath.c_str(), saveName.c_str(), thisRegion.getName().c_str() ) );
  h2_axes_yield->SetYTitle("#gamma Data/MC ratio"); h2_axes_yield->Draw();  
  g_clone->Draw("p same");    labelTop->Draw("same");   gPad->RedrawAxis();    //  gPad->SetLogy();
  canny->SaveAs( Form("%s/yield_gamma_%s_%s.eps",  fullPath.c_str(), saveName.c_str(),thisRegion.getName().c_str() ) );
  
  delete    h2_axes_yield;


  for(int bi = 1; bi <= yieldBins ; bi++){
    std::cout << "MC Zll / Gamma = " << h_mt2_mc->GetBinContent(bi) << " / " << g_mt2_mc->GetBinContent(bi)/1.23 << std::endl;
  }
 
  h_mt2_mc->Divide(g_mt2_mc);

  //Filling the single ratios
  for(int bi = 1; bi <= yieldBins ; bi++){
    double value_mc = h_mt2_mc->GetBinContent(bi);
    double err_mc = h_mt2_mc->GetBinError(bi);
    zllG_mc->get(thisRegion)->yield->SetBinContent( bi,   value_mc);
    zllG_mc->get(thisRegion)->yield_systUp->SetBinContent( bi,   value_mc + err_mc);
    zllG_mc->get(thisRegion)->yield_systDown->SetBinContent( bi,   value_mc - err_mc);
  }


  TGraphAsymmErrors* gr_data = zllG_data->get(thisRegion)->getGraph();
  gr_data->SetMarkerSize(1.4);
  gr_data->SetMarkerStyle(20);
  gr_data->SetLineColor(kBlack);
  gr_data->SetLineWidth(2);
  gr_data->SetMarkerColor(kBlack);
  std::cout << std::endl; 

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


  h_mt2_mc->SetLineColor( 38 ); h_mt2_mc->SetLineWidth(2);
  float yMax = 0.13;
  if( zll_sel == "nBJets") yMax = 0.18;
  //  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 0.125 );
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0.0, yMax);
  if(zll_sel == "mt2"){
    h2_axes->SetXTitle("M_{T2} [GeV]");
  }else    if(zll_sel == "ht"){
    h2_axes->SetXTitle("H_{T} [GeV]");
  }else    if(zll_sel == "nJets"){
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

  labelTop->Draw("same");
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

  TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , 0.3, 1.7 );
 

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
    h2_axes_rat->GetXaxis()->SetTitleOffset(5);
    h2_axes_rat->GetXaxis()->SetLabelSize(0.00);
    h2_axes_rat->GetXaxis()->SetTickLength(0.09);
    h2_axes_rat->GetYaxis()->SetNdivisions(5,5,0);
    h2_axes_rat->GetYaxis()->SetTitleSize(0.2);
    h2_axes_rat->GetYaxis()->SetTitleOffset(0.34);
    h2_axes_rat->GetYaxis()->SetLabelSize(0.17);
  
  
    h2_axes_rat->Draw();
    gr_ratio->Draw("same P");


    TLine* line = new TLine(xMin,1 , xMax, 1);
    line->SetLineColor(kBlack);
    line->Draw("same");
 
    SFFitBand->Draw("3,same");
    fSF->Draw("same");


    gPad->RedrawAxis();

    canny->cd();
    pad1->cd();
    fitText->Draw("same");
    
    }

    canny->SaveAs( Form("%s/%s_ratios_%s.eps", fullPath.c_str(), saveName.c_str(),  thisRegion.getName().c_str() ) );
    canny->SaveAs( Form("%s/%s_ratios_%s.png", fullPath.c_str(), saveName.c_str(),thisRegion.getName().c_str() ) );
    canny->SaveAs( Form("%s/%s_ratios_%s.pdf", fullPath.c_str(), saveName.c_str(),thisRegion.getName().c_str() ) );

    delete h_mt2;   delete h_mt2_mc;    delete g_mt2; delete g_mt2_mc;
    delete h2_axes; delete h2_axes_rat;
    delete gr_ratio;
    delete canny;

    delete h_top;
    delete g_Up; delete g_Down;

}

















