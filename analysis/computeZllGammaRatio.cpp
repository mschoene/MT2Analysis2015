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

void drawRatios(std::string fullPath, double *binss, unsigned int size,  std::string zll_sel, MT2Analysis<MT2EstimateSyst>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma_mc, MT2Analysis<MT2EstimateTree>*  gamma_data, MT2Analysis<MT2EstimateSyst>*  purity, MT2Analysis<MT2EstimateTree>*  zll_mc,MT2Analysis<MT2EstimateTree>*  zll_data,   MT2Analysis<MT2EstimateTree>*  top, MT2Analysis<MT2EstimateSyst>*  zll_yield, MT2Analysis<MT2EstimateSyst>*  zllG_data,  MT2Analysis<MT2EstimateSyst>*  zllG_mc, const MT2Region thisRegion, std::string cut,  std::string cut_gamma, std::string cut_data, std::string cut_gamma_data, float lumi, std::string saveName, bool onlyMC, std::string topoCuts="" );

TH1D drawBGsubtraction(  MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* zllMC, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units, float scaleFactor );
 
void drawCorrelation(std::string fullPath, float *binss, unsigned int size,  float *binss2, unsigned int size2,  std::string zll_sel,  std::string zll_sel2, MT2Analysis<MT2EstimateTree>*  Zll,  const MT2Region thisRegion, std::string cut , float lumi);


//void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string name ,  MT2Analysis<MT2Estimate>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma, MT2Analysis<MT2EstimateTree>*  Zll, MT2Analysis<MT2Estimate>*  zll_yield, const MT2Region thisRegion, std::string cut);


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

    
  float scaleFactor=1;



 
  MT2Analysis<MT2EstimateTree>* gamma_mc = MT2Analysis<MT2EstimateTree>::readFromFile(gammaControlRegionDir + "/mc.root", "gammaCRtree");
  MT2Analysis<MT2EstimateTree>*  gamma_data = MT2Analysis<MT2EstimateTree>    ::readFromFile( gammaControlRegionDir + "/data.root", "gammaCRtree");
 
  MT2Analysis<MT2EstimateSyst>* purity = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_data.root", "purity");

  MT2Analysis<MT2EstimateSyst>* purity_ht = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_ht_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_njets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_njets_data.root", "purity");
  MT2Analysis<MT2EstimateSyst>* purity_nbjets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_nbjets_data.root", "purity");
 
  //  MT2Analysis<MT2EstimateSyst>* purity_mono_nbjets = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/PurityFitsRC/purityFit_mono_nbjets_data.root", "purity");
 

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
  //  MT2Analysis<MT2EstimateSyst>* zll_mono_nBJets = new MT2Analysis<MT2EstimateSyst>("zll_mono_nBJets",regionsSet.c_str());

  //DATA RATIOS
  MT2Analysis<MT2EstimateSyst>* zllG_data_mt2 = new MT2Analysis<MT2EstimateSyst>( "zllG_data_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_data_ht = new MT2Analysis<MT2EstimateSyst>( "zllG_data_ht", regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_data_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_data_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_data_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_data_nBJets",regionsSet.c_str() );
  //  MT2Analysis<MT2EstimateSyst>* zllG_data_mono_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_data_mono_nBJets",regionsSet.c_str() );

  //MC RATIOS
  MT2Analysis<MT2EstimateSyst>* zllG_mc_mt2 = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_mc_ht = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_ht", regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_mc_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_mc_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_mc_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_mc_nBJets",regionsSet.c_str() );
  //  MT2Analysis<MT2EstimateSyst>* zllG_mc_mono_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_mc_mono_nBJets",regionsSet.c_str() );


  //RATIOS
  MT2Analysis<MT2EstimateSyst>* zllG_pt = new MT2Analysis<MT2EstimateSyst>( "zllG_pt", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_mt2 = new MT2Analysis<MT2EstimateSyst>( "zllG_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2EstimateSyst>* zllG_ht = new MT2Analysis<MT2EstimateSyst>( "zllG_ht", regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_nJets = new MT2Analysis<MT2EstimateSyst>( "zllG_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2EstimateSyst>* zllG_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_nBJets",regionsSet.c_str() );
  //  MT2Analysis<MT2EstimateSyst>* zllG_mono_nBJets = new MT2Analysis<MT2EstimateSyst>("zllG_mono_nBJets",regionsSet.c_str() );
  

  //std::set<MT2Region> MT2Regions = zll_ratio->getRegions();
  //std::set<MT2Region> Inclusive = zll_mc->getRegions();
  // MT2Region RegIncl( (*Inclusive.begin()) );
  std::set<MT2Region> MT2Regions = zll_mc->getRegions();

 

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  TH1F::AddDirectory(kTRUE);


  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
    std::vector<std::string> niceNames = thisRegion.getNiceNames();
    float lumi = cfg.lumi();    

    double bins_nJets[] = {1,2,4,7,12};
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
  

    std::string cut = "weight*(abs(Z_mass-91.19)<10  && mt2>200 && ht>200 && nJets>0 )";
    std::string cut_el = "weight*(abs(Z_mass-91.19)<10 && ht>200 && mt2>200 && nJets>0 && Z_lepId==11 )";
    std::string cut_mu = "weight*(abs(Z_mass-91.19)<10 && ht>200&& mt2>200 && nJets>0 && Z_lepId==13 )";
    std::string cut_mono = "weight*(abs(Z_mass-91.19)<10 && ht>200 && nJets==1 )";

    std::string cut_mono_el = "weight*(abs(Z_mass-91.19)<10 && ht>200 && nJets==1 && Z_lepId==11 )";
    std::string cut_mono_mu = "weight*(abs(Z_mass-91.19)<10 && ht>200 && nJets==1 && Z_lepId==13 )";

    std::string cut_gamma = "weight*(prompt==2 && iso<2.5 &&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23";
    std::string cut_gamma_mono = "weight*( prompt==2 && iso<2.5 && ptGamma>180 && nJets==1 && ht>200 )*1.23";
 
    //f = 0.92, purity later in the function
    std::string cut_data = "weight*(abs(Z_mass-91.19)<10 &&  mt2>200 && ht>200 && nJets>0 )";
    std::string cut_el_data = "weight*(abs(Z_mass-91.19)<10 && ht>200 && mt2>200 && nJets>0&& Z_lepId==11 )";
    std::string cut_mu_data = "weight*(abs(Z_mass-91.19)<10 && ht>200&& mt2>200 && nJets>0&& Z_lepId==13 )";
    std::string cut_mono_data = "weight*(abs(Z_mass-91.19)<10 && ht>200 && nJets==1 )";

    std::string cut_mono_data_el = "weight*(abs(Z_mass-91.19)<10 && ht>200 && nJets==1 && Z_lepId==11)";
    std::string cut_mono_data_mu = "weight*(abs(Z_mass-91.19)<10 && ht>200 && nJets==1 && Z_lepId==13)";

    std::string cut_gamma_data = "weight*( iso<2.5 && ptGamma>180  && nJets>0 && mt2>200 && ht>200 )";
    std::string cut_gamma_mono_data = "weight*( iso<2.5 &&  ptGamma>180 && nJets==1 && ht>200 )";

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
  

    MT2EstimateSyst::rebinYields( zll_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zll_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zll_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zll_nBJets,  size_nBJets, bins_nBJets);
    //    MT2EstimateSyst::rebinYields( zll_mono_nBJets,  size_mono_nBJets, bins_mono_nBJets);
 
    MT2EstimateSyst::rebinYields( zllG_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zllG_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zllG_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zllG_nBJets,  size_nBJets, bins_nBJets);
    //    MT2EstimateSyst::rebinYields( zllG_mono_nBJets,  size_mono_nBJets, bins_mono_nBJets);
  
    MT2EstimateSyst::rebinYields( zllG_data_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zllG_data_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zllG_data_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zllG_data_nBJets,  size_nBJets, bins_nBJets);
    //    MT2EstimateSyst::rebinYields( zllG_data_mono_nBJets,  size_mono_nBJets, bins_mono_nBJets);
 
    MT2EstimateSyst::rebinYields( zllG_mc_mt2,  size_mt2, bins_mt2);
    MT2EstimateSyst::rebinYields( zllG_mc_ht,  size_ht, bins_ht);
    MT2EstimateSyst::rebinYields( zllG_mc_nJets,  size_nJets, bins_nJets);
    MT2EstimateSyst::rebinYields( zllG_mc_nBJets,  size_nBJets, bins_nBJets);
    //    MT2EstimateSyst::rebinYields( zllG_mc_mono_nBJets,  size_mono_nBJets, bins_mono_nBJets);
 

    //draw ratio also fills the ratio and yield estimates
    //outputdir, bins, nbins, var to project, ratio estimate, gamma mc, gamma data, purity gamma, zll mc, zll data, yield estimate, region, cut zll, cut gamma, cut zll data, cut gamma data, lumi, name, flag , topo region);
    
    drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data, top, zll_mt2 , zllG_data_mt2, zllG_mc_mt2 , thisRegion, cut_el, cut_gamma, cut_el_data,  cut_gamma_data,  lumi , "mt2_el" , onlyMC ,"#geq2j, #geq0b" );
    drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data, top, zll_mt2, zllG_data_mt2, zllG_mc_mt2 , thisRegion, cut_mu,  cut_gamma, cut_mu_data, cut_gamma_data, lumi , "mt2_mu" , onlyMC ,"#geq2j, #geq0b");
    drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data, top, zll_mt2, zllG_data_mt2, zllG_mc_mt2 ,thisRegion, cut, cut_gamma, cut_data,cut_gamma_data, lumi , "mt2" , onlyMC ,"#geq2j, #geq0b");
   
    drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_el, cut_gamma, cut_el_data,  cut_gamma_data, lumi, "ht_el" , onlyMC ,"#geq2j, #geq0b");
    drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut_mu, cut_gamma, cut_mu_data,  cut_gamma_data, lumi,"ht_mu" , onlyMC ,"#geq2j, #geq0b"); 
    drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity_ht,  zll_mc, zll_data, top, zll_ht, zllG_data_ht, zllG_mc_ht , thisRegion, cut, cut_gamma, cut_data,  cut_gamma_data, lumi,"ht" , onlyMC ,"#geq2j, #geq0b");
    
    drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion, cut_el, cut_gamma, cut_el_data,  cut_gamma_data, lumi, "nJets_el" , onlyMC ,"#geq2j, #geq0b");
   drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion, cut_mu, cut_gamma, cut_mu_data,  cut_gamma_data, lumi, "nJets_mu" , onlyMC ,"#geq2j, #geq0b");
   drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity_njets,  zll_mc, zll_data,top,  zll_nJets, zllG_data_nJets, zllG_mc_nJets ,  thisRegion,cut, cut_gamma, cut_data,  cut_gamma_data, lumi, "nJets" , onlyMC ,"#geq2j, #geq0b");

   drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_el, cut_gamma ,  cut_el_data,  cut_gamma_data, lumi, "nBJets_el" , onlyMC ,"#geq2j, #geq0b");
   drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut_mu, cut_gamma ,  cut_mu_data,  cut_gamma_data, lumi, "nBJets_mu" , onlyMC ,"#geq2j, #geq0b");
   drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity_nbjets,  zll_mc, zll_data, top, zll_nBJets,zllG_data_nBJets, zllG_mc_nBJets ,   thisRegion, cut, cut_gamma ,  cut_data,  cut_gamma_data, lumi, "nBJets" , onlyMC ,"#geq2j, #geq0b");


//   //MONOJET
//   drawRatios( outputdir, bins_mono_nBJets, size_mono_nBJets , "nBJets",  zllG_mono_nBJets,   gamma_mc, gamma_data, purity_mono_nbjets,  zll_mc, zll_data, top, zll_mono_nBJets, zllG_data_mono_nBJets, zllG_mc_mono_nBJets , thisRegion, cut_mono_el, cut_gamma_mono ,  cut_mono_data_el,  cut_gamma_mono_data, lumi, "mono_nBJets_el" , onlyMC ,"#geq1j, #geq0b");
//   drawRatios( outputdir, bins_mono_nBJets, size_mono_nBJets , "nBJets",  zllG_mono_nBJets,   gamma_mc, gamma_data, purity_mono_nbjets,  zll_mc, zll_data, top, zll_mono_nBJets, zllG_data_mono_nBJets, zllG_mc_mono_nBJets , thisRegion, cut_mono_mu, cut_gamma_mono ,  cut_mono_data_mu,  cut_gamma_mono_data, lumi, "mono_nBJets_mu" , onlyMC ,"#geq1j, #geq0b");
//
//    drawRatios( outputdir, bins_mono_nBJets, size_mono_nBJets , "nBJets",  zllG_mono_nBJets,   gamma_mc, gamma_data, purity_mono_nbjets,  zll_mc, zll_data, top, zll_mono_nBJets, zllG_data_mono_nBJets, zllG_mc_mono_nBJets , thisRegion, cut_mono, cut_gamma_mono ,  cut_mono_data,  cut_gamma_mono_data, lumi, "mono_nBJets" , onlyMC ,"#geq1j, #geq0b");

    


    //  TH1D data_bgSub =  drawBGsubtraction(  cfg, zll_data, zll_mc,  bgYields , "mt2", "mt2", cut, 20, 200, 800, "M_{T2}", "GeV", scaleFactor );
    //  TH1D data_bgSub_nBJets =  drawBGsubtraction(  cfg, zll_data, zll_mc,  bgYields , "nBJets", "nBJets", cut, 8, 0, 8, "nBJets", "", scaleFactor );
    /*
    drawCorrelation( outputdir , bins_mt22, size_mt2 ,  bins_ht2, size_ht,  "mt2", "ht",  Zll, thisRegion,  cut_corr, lumi );
    drawCorrelation( outputdir , bins_mt22, size_mt2 ,  bins_nJets2, size_nJets,  "mt2", "nJets",  Zll, thisRegion,  cut_corr , lumi);
    drawCorrelation( outputdir , bins_mt22, size_mt2 ,  bins_nBJets2, size_nBJets,  "mt2", "nBJets",  Zll, thisRegion,  cut_corr , lumi);
    drawCorrelation( outputdir , bins_ht2, size_ht  ,  bins_nJets2, size_nJets ,  "ht", "nJets",  Zll, thisRegion,  cut_corr , lumi);
//  drawCorrelation( outputdir , bins_ht,sizeof(bins_ht)/sizeof(float)-1   ,  bins_nJets, sizeof(bins_nJets)/sizeof(float)-1  ,  "ht", "nJets",  Zll, thisRegion,  cut_corr );
    drawCorrelation( outputdir , bins_ht2, size_ht ,  bins_nBJets2, size_nBJets,  "ht", "nBJets",  Zll, thisRegion,  cut_corr, lumi );
    drawCorrelation( outputdir , bins_nJets2, size_nJets ,  bins_nBJets2, size_nBJets,  "nJets", "nBJets",  Zll, thisRegion,  cut_corr, lumi );
    */
  }
  
 
  std::string outFile = outputdir + "/zll_ratio.root";
  zllG_mt2->writeToFile(outFile);
  zllG_ht->addToFile( outFile );
  zllG_nJets->addToFile( outFile );
  zllG_nBJets->addToFile( outFile );
  //  zllG_mono_nBJets->addToFile( outFile );

  std::string outFile_yield = outputdir + "/zll_yield.root";
  zll_mt2->writeToFile(outFile_yield);
  zll_ht->addToFile( outFile_yield );
  zll_nJets->addToFile( outFile_yield );
  zll_nBJets->addToFile( outFile_yield );
  //  zll_mono_nBJets->addToFile( outFile_yield );


  std::string outFile_data = outputdir + "/zllG_data_ratio.root";
  zllG_data_mt2->writeToFile(outFile_data);
  zllG_data_ht->addToFile( outFile_data );
  zllG_data_nJets->addToFile( outFile_data );
  zllG_data_nBJets->addToFile( outFile_data );
  //  zllG_data_mono_nBJets->addToFile( outFile_data );

 std::string outFile_mc = outputdir + "/zllG_mc_ratio.root";
  zllG_mc_mt2->writeToFile(outFile_mc);
  zllG_mc_ht->addToFile( outFile_mc );
  zllG_mc_nJets->addToFile( outFile_mc );
  zllG_mc_nBJets->addToFile( outFile_mc );
  //  zllG_mc_mono_nBJets->addToFile( outFile_mc );


  return 0;

}


void drawRatios(std::string fullPath, double *binss, unsigned int size,  std::string zll_sel, MT2Analysis<MT2EstimateSyst>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma_mc, MT2Analysis<MT2EstimateTree>*  gamma_data, MT2Analysis<MT2EstimateSyst>*  purity, MT2Analysis<MT2EstimateTree>*  zll_mc,MT2Analysis<MT2EstimateTree>*  zll_data, MT2Analysis<MT2EstimateTree>*  top, MT2Analysis<MT2EstimateSyst>*  zll_yield, MT2Analysis<MT2EstimateSyst>*  zllG_data,  MT2Analysis<MT2EstimateSyst>*  zllG_mc, const MT2Region thisRegion, std::string cut, std::string cut_gamma, std::string cut_data, std::string cut_gamma_data, float lumi, std::string saveName, bool onlyMC, std::string topoCuts ){
 
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
  //h_top->Scale( scaleFactor ); //to come soon, hopefully

  h_top->SetBinContent(size, h_top->GetBinContent(size) + h_top->GetBinContent(size+1));//adding overflow
  h_mt2->SetBinContent(size, h_mt2->GetBinContent(size) + h_mt2->GetBinContent(size+1));//adding overflow
  g_mt2->SetBinContent(size, g_mt2->GetBinContent(size) + g_mt2->GetBinContent(size+1));
  h_mt2->SetBinContent(size+1,0);
  h_top->SetBinContent(size+1,0);
  g_mt2->SetBinContent(size+1,0);


  //GET THE UNCERTAINTIES
  float f = 0.92;
  float f_uncert = 0.1; // 0.08 from the fragmentation and then o+ ~5% for the mc closure
  TH1D* g_Up = new TH1D("g_Up","", size , bins);
  TH1D* g_Down = new TH1D("g_Down","", size , bins);

  int nBinss =  h_mt2->GetNbinsX();
  for(int binnie = 1; binnie <= nBinss; binnie++){
    //Zll//////////////
    double value = h_mt2->GetBinContent(binnie);
    double top = h_top->GetBinContent(binnie);
 
    h_mt2->SetBinContent(binnie, value - top);
    //add 50% of the top bg as uncertainty to the estimate
    h_mt2->SetBinError(binnie,  sqrt(value + top*top*0.25)  );

    //Photons//////////
    Double_t x_tmp, p, p_errUp, p_errDown;	       
    this_zinv_purity->GetPoint( binnie-1, x_tmp, p);
    p_errUp   = this_zinv_purity->GetErrorYhigh(binnie-1);
    p_errDown = this_zinv_purity->GetErrorYlow(binnie-1);
 
    std::cout << "Purity = " << p << std::endl;

    double value_g = g_mt2->GetBinContent(binnie);

    g_mt2->SetBinContent(binnie, value_g  * p * f);

    float uncertUp = sqrt ( value_g*f*f*p*p + ( value_g*value_g*f*f *p_errUp*p_errUp) + (p*p*value_g*value_g* f_uncert*f_uncert) );
    float uncertDown = sqrt ( value_g*f*f*p*p + ( value_g*value_g*f*f *p_errDown*p_errDown) + (p*p*value_g*value_g* f_uncert*f_uncert) );
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

  }

  *(zllG_data) = *(zll_yield) / *(zllG_data);

 


  TH1D* h_mt2_mc = new TH1D("h_mt2_mc","", size  , bins);
  TH1D* g_mt2_mc = new TH1D("g_mt2_mc","", size  , bins);
  h_mt2_mc->Sumw2(); g_mt2_mc->Sumw2();

//  zll_mc_tree ->Project( "h_mt2_mc" , zll_sel.c_str(), Form("(%s)*%f",cut.c_str(), lumi) );
//  gamma_mc_tree->Project( "g_mt2_mc", gamma_sel.c_str(),  Form("(%s)*%f",cut_gamma.c_str(), lumi)  );
//  h_mt2_mc->Scale( lumi );
//  g_mt2_mc->Scale( lumi );

  zll_mc_tree ->Project( "h_mt2_mc" , zll_sel.c_str(), Form("(%s)*%f",cut.c_str(), 1.0) );
  gamma_mc_tree->Project( "g_mt2_mc", gamma_sel.c_str(),  Form("(%s)*%f",cut_gamma.c_str(), 1.0)  );
 
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
    std::cout << "MC Zll / Gamma = " << h_mt2_mc->GetBinContent(bi) << " / " << g_mt2_mc->GetBinContent(bi) << std::endl;
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
    if(p<0.01){
      gr_data->RemovePoint(bi-1);
      gr_ratio->RemovePoint(bi-1);
    }
  }


  h_mt2_mc->SetLineColor( 38 ); h_mt2_mc->SetLineWidth(2);
 
  //  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 0.125 );
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0.0, 0.2 );
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

  TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , 0., 2.0 );
 

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




















TH1D drawBGsubtraction(  MT2Config cfg, MT2Analysis<MT2EstimateTree>* data,  MT2Analysis<MT2EstimateTree>* zllMC, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units, float scaleFactor ) {


  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;

  std::vector<int> colors;
  if( bgYields.size()==3 ) { // estimates
    colors.push_back(402); 
    colors.push_back(430); 
    colors.push_back(418); 
  } else { // mc
    //   colors.push_back(430); // other=zll
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    colors.push_back(419); // z+jets
    colors.push_back(855); // top
    //colors.push_back(); // other
  }

  TH1D* h1_data_bgSub = new TH1D("h1_data_bgSub", "", nBins, xMin, xMax);


  std::string fullPathPlots = cfg.getEventYieldDir() + "/plotsBGsubtr";
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

    TTree* tree_zllMC = zllMC->get(thisRegion)->tree;
    TH1D* h1_zllMC = new TH1D("h1_zllMC", "", nBins, xMin, xMax );
    tree_zllMC->Project( "h1_zllMC", varName.c_str(), selection.c_str() );
    h1_zllMC->SetFillColor(430);
    h1_zllMC->SetLineColor( kBlack );

    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<bgYields.size(); ++i ) { 
      TTree* tree_mc = (bgYields[i]->get(thisRegion)->tree);
      std::string thisName = "h1_" + bgYields[i]->getName();
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_mc->Sumw2();
      if( selection!="" )
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%s", selection.c_str()) );
      else
        tree_mc->Project( thisName.c_str(), varName.c_str(), "" );
      if(i==3)
	h1_mc->Scale(scaleFactor);
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

    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

    std::cout << "Integrals: " << h1_data->Integral(0, nBins+1) << "\t" << h1_zllMC->Integral(0, nBins+1) << std::endl;
    float ratio = h1_data->Integral(0, nBins+1)/h1_zllMC->Integral(0, nBins+1);   
    std::cout << "SF: " << ratio << std::endl;

    h1_data_bgSub = (TH1D*)h1_data->Clone();
    h1_data_bgSub->Add( mc_sum, -1);

    TGraphAsymmErrors* gr_data_bgSub = MT2DrawTools::getPoissonGraph(h1_data_bgSub);
    gr_data_bgSub->SetMarkerStyle(21);
    gr_data_bgSub->SetMarkerSize(1.2);

    /*
    TH1D* histo_mc;
    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = bgYields.size() - i - 1;
      histos_mc[index]->SetFillColor( colors[index] );
      histos_mc[index]->SetLineColor( kBlack );

      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack.Add(histos_mc[index]);
    }
    */

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
    float yMax3 = yMaxScale*(h1_zllMC->GetMaximum());
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
      regionText->Draw("same");
  
      pad1_log->cd();
      regionText->Draw("same");
    }
    
    TLegend* legend = new TLegend( 0.7, 0.9-(3)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( gr_data, "Data", "P" );
    legend->AddEntry( gr_data_bgSub, "Data BG Subtracted", "P" );
    legend->AddEntry( h1_zllMC, "Z+jets", "F");
    /*
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], bgYields[i]->getFullName().c_str(), "F" );
    }
    */

    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    
    TPaveText* ratioText = new TPaveText( 0.133, -0.051, 0.4, 0.1 , "brNDC" );
    ratioText->SetTextSize(0.04);
    ratioText->SetTextFont(40);
    ratioText->SetTextColor(2);
    ratioText->SetFillColor(0);
    ratioText->SetTextAlign(11);
    ratioText->AddText( Form("Data/MC = %.2f", ratio) );
    //  ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", scaleFactor, error_datamc) );
     

    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    
    TLine* lineSF = new TLine(xMin, ratio, xMax, ratio);
    lineSF->SetLineColor(2);

    float yMinR=0.0;
    float yMaxR=2.0;

  
    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, h1_zllMC);
 
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
    h1_zllMC->Draw("histo same");
    gr_data->Draw("p same");
    gr_data_bgSub->Draw("p same");
    labelTop->Draw("same");
    ratioText->Draw("same");
  
    gPad->RedrawAxis();

    c1_log->cd();
    pad1_log->cd();
    legend->Draw("same");
    h1_zllMC->Draw("histo same");
    gr_data->Draw("p same");
    gr_data_bgSub->Draw("p same");
    labelTop->Draw("same");
    ratioText->Draw("same");

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


      systBand->Draw("3,same");
      lineCentral->Draw("same");

      SFFitBand->Draw("3,same");
      fSF->Draw("same");
 

    g_ratio->Draw("PE,same");    
    gPad->RedrawAxis();


    c1_log->cd();
    // TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
    pad2_log->Draw();
    pad2_log->cd();

    h2_axes_ratio->Draw(""); 
    
    lineCentral->Draw("same");

      systBand->Draw("3,same");
      lineCentral->Draw("same");

      SFFitBand->Draw("3,same");
      fSF->Draw("same");
 
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

  return *h1_data_bgSub;
}


























void drawCorrelation(std::string fullPath, float *binss, unsigned int size,  float *binss2, unsigned int size2,  std::string zll_sel,  std::string zll_sel2, MT2Analysis<MT2EstimateTree>*  Zll,  const MT2Region thisRegion, std::string cut, float lumi){


  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.12);

  gStyle->SetPalette(51,0);
 
  TH1F::AddDirectory(kTRUE);

  float bins[size+1]; for(unsigned int i=0; i<= size ; i++)      bins[i]=binss[i];
  float bins2[size2+1]; for(unsigned int i=0; i<= size2 ; i++)      bins2[i]=binss2[i];
 
  float xMin = bins[0];
  float xMax = bins[size];
  float yMin = bins2[0];
  float yMax = bins2[size2];

  //THE TREES
  TTree *zllT =  Zll->get(thisRegion)->tree;

  //  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
 
  TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
  canny->cd();

  
  TH2D* histo = new TH2D("histo","", size, bins, size2, bins2);
  // TH2D* histo = new TH2D("histo","", size*50, xMin, xMax, size2*50, yMin, yMax);
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(0.5);
   

  //  zllT ->Draw(  Form("%s:%s>> histo", zll_sel.c_str(), zll_sel.c_str()), cut.c_str()  );
  zllT ->Project( "histo" , Form("%s:%s", zll_sel2.c_str(), zll_sel.c_str())  );
  //   zllT ->Project( "histo" , Form("%s:%s", zll_sel.c_str(), zll_sel2.c_str()) , cut.c_str() );
 
  //histo = (TH2D*)gDirectory->Get("histo");

  double corr =  histo->GetCorrelationFactor();
  std::cout << corr << std::endl;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax);
  if(zll_sel == "zll_mt2"){
    h2_axes->SetXTitle("M_{T2} [GeV]");
  }else    if(zll_sel == "zll_ht"){
    h2_axes->SetXTitle("H_{T} [GeV]");
  }else    if(zll_sel == "nJets"){
    h2_axes->SetXTitle("Jet Multiplicity");
  }else{
    h2_axes->SetXTitle("b Jet Multiplicity" );
  }

  if(zll_sel2 == "zll_mt2"){
    h2_axes->SetYTitle("M_{T2} [GeV]");
  }else    if(zll_sel2 == "zll_ht"){
    h2_axes->SetYTitle("H_{T} [GeV]");
  }else    if(zll_sel2 == "nJets"){
    h2_axes->SetYTitle("Jet Multiplicity");
  }else{
    h2_axes->SetYTitle("b Jet Multiplicity" );
  }

  gPad->SetLogz();


  h2_axes->Draw();
  histo->Draw("colz same");
  //  labelTop->Draw("same");

  TPaveText* regionText = new TPaveText( 0.45, 0.91-0.03, 0.78, 0.91, "brNDC" );
  regionText->SetTextSize(0.04);
  // regionText->SetTextFont(42);
  regionText->SetFillColor(0);
  regionText->SetTextAlign(11);
  regionText->AddText( Form("Linear Correlation = %.2f",corr));
  regionText->Draw("same"); 

  gPad->RedrawAxis();


  /*
    TLegend* legend = new TLegend( 0.7, 0.92-(2)*0.06, 0.9, 0.92 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( h_mt2 ,"Data", "P" );
    legend->AddEntry( h_mt2_mc ,"Simulation", "L" );
    legend->Draw("same");


    std::vector<std::string> niceNames2 = thisRegion.getNiceNames();

    for( unsigned i=0; i< niceNames2.size(); ++i ) {
    float yMaxText = 0.9-(float)i*0.05;
    float yMinText = yMaxText - 0.05;
    TPaveText* regionText = new TPaveText( 0.18, yMinText, 0.55, yMaxText, "brNDC" );
    regionText->SetTextSize(0.035);
    regionText->SetTextFont(42);
    regionText->SetFillColor(0);
    regionText->SetTextAlign(11);
    regionText->AddText( niceNames2[i].c_str() );
    regionText->Draw("same");
    }

  */
 
 
  canny->SaveAs( Form("%s/correlation_%s_vs_%s.eps", fullPath.c_str(), zll_sel.c_str(),  zll_sel2.c_str() ) );
  canny->SaveAs( Form("%s/correlation_%s_vs_%s.png", fullPath.c_str(), zll_sel.c_str(),  zll_sel2.c_str() ) );
  canny->SaveAs( Form("%s/correlation_%s_vs_%s.pdf", fullPath.c_str(), zll_sel.c_str(), zll_sel2.c_str()  ) );



  delete histo;    delete canny; delete h2_axes;

}






















