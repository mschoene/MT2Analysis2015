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

void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string zll_sel, MT2Analysis<MT2Estimate>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma_mc, MT2Analysis<MT2EstimateTree>*  gamma_data, MT2Analysis<MT2EstimateSyst>*  purity, MT2Analysis<MT2EstimateTree>*  zll_mc,MT2Analysis<MT2EstimateTree>*  zll_data, MT2Analysis<MT2Estimate>*  zll_yield, const MT2Region thisRegion, std::string cut, std::string cut_gamma, float lumi, std::string saveName );

TH1D drawBGsubtraction(  MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* zllMC, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units, float scaleFactor );
 
void drawCorrelation(std::string fullPath, float *binss, unsigned int size,  float *binss2, unsigned int size2,  std::string zll_sel,  std::string zll_sel2, MT2Analysis<MT2EstimateTree>*  Zll,  const MT2Region thisRegion, std::string cut , float lumi);


//void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string name ,  MT2Analysis<MT2Estimate>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma, MT2Analysis<MT2EstimateTree>*  Zll, MT2Analysis<MT2Estimate>*  zll_yield, const MT2Region thisRegion, std::string cut);


int main(int argc, char* argv[]){

  
  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeZinvFromGamma             |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc!=2 ) {
    std::cout << "USAGE: ./computeZllGammaRatio [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

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
  SF_file.open("scaleFactorOF.txt");
  float scaleFactor;
  SF_file >> scaleFactor;
  std::cout << scaleFactor << std::endl;


 
  MT2Analysis<MT2EstimateTree>* gamma_mc = MT2Analysis<MT2EstimateTree>::readFromFile(gammaControlRegionDir + "/mc.root", "gammaCRtree");
  MT2Analysis<MT2EstimateTree>*  gamma_data = MT2Analysis<MT2EstimateTree>    ::readFromFile( gammaControlRegionDir + "/data.root", "gammaCRtree");
  MT2Analysis<MT2EstimateSyst>* purity = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaControlRegionDir + "/purityMC.root", "purity");

 
  //  MT2Analysis<MT2EstimateTree>* gamma = MT2Analysis<MT2EstimateTree>::readFromFile(gammaControlRegionDir + "/data.root", "gammaCRtree");
  if( gamma_mc==0 ) {
    std::cout << "-> Please run gammaControlRegion first. I need to get the gammaCR yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(193);
  }
 
  MT2Analysis<MT2EstimateTree>* zll_mc = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc.root", ZllDir.c_str()) , "zllCR");

 MT2Analysis<MT2EstimateTree>* zll_data = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/data.root", ZllDir.c_str()) , "data");
  if( zll_data==0 ) {
    std::cout << "-> Please run computeZinvFromZll first. I need to get the Z->ll MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }
  

  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str()  ), "QCD");
  MT2Analysis<MT2EstimateTree>* top = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "Top");
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "WJets");

  wjets->setFullName("W+jets");
  //  zjets->setFullName("Z#nu#nu+jets");

  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
  bgYields.push_back( qcd );
  bgYields.push_back( wjets );
  bgYields.push_back( top );
  

  MT2Analysis<MT2Estimate>* zll_ratio = new MT2Analysis<MT2Estimate>( "zll_ratio", regionsSet.c_str() );
  MT2Analysis<MT2Estimate>* zll_yield = new MT2Analysis<MT2Estimate>( "zll_yield", regionsSet.c_str() );

  //YIELDS
  MT2Analysis<MT2Estimate>* zllY_pt = new MT2Analysis<MT2Estimate>( "zllY_pt", regionsSet.c_str() ); 
  MT2Analysis<MT2Estimate>* zllY_mt2 = new MT2Analysis<MT2Estimate>( "zllY_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2Estimate>* zllY_ht = new MT2Analysis<MT2Estimate>( "zllY_ht", regionsSet.c_str() ); 
  MT2Analysis<MT2Estimate>* zllY_nJets = new MT2Analysis<MT2Estimate>( "zllY_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2Estimate>* zllY_nBJets = new MT2Analysis<MT2Estimate>("zllY_nBJets",regionsSet.c_str());

  //RATIOS
  MT2Analysis<MT2Estimate>* zllG_pt = new MT2Analysis<MT2Estimate>( "zllG_pt", regionsSet.c_str() ); 
  MT2Analysis<MT2Estimate>* zllG_mt2 = new MT2Analysis<MT2Estimate>( "zllG_mt2", regionsSet.c_str() ); 
  MT2Analysis<MT2Estimate>* zllG_ht = new MT2Analysis<MT2Estimate>( "zllG_ht", regionsSet.c_str()); 
  MT2Analysis<MT2Estimate>* zllG_nJets = new MT2Analysis<MT2Estimate>( "zllG_nJets",regionsSet.c_str()); 
  MT2Analysis<MT2Estimate>* zllG_nBJets = new MT2Analysis<MT2Estimate>("zllG_nBJets",regionsSet.c_str() );
  


  int counter_total_bins = 0;
  int counter_empty_bins = 0;

  //Loop over regions

  //std::set<MT2Region> MT2Regions = zll_ratio->getRegions();
  std::set<MT2Region> MT2Regions = zll_mc->getRegions();
  //  std::set<MT2Region> Inclusive = zll_mc->getRegions();
  // MT2Region RegIncl( (*Inclusive.begin()) );

 
  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  TH1F::AddDirectory(kTRUE);

  //THE TREES
  TTree *gammaT =  gamma_mc->get( *MT2Regions.begin() )->tree;
  TTree *zllT =  zll_mc->get( *MT2Regions.begin() )->tree;

 


  TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi() );
 

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
    // std::set<MT2Region>::iterator iMT2 = MT2Regions.begin();  
    std::vector<std::string> niceNames = thisRegion.getNiceNames();
 
    /*
    float bins_nJets[] = {2, 6, 12};
    float bins_nBJets[] = {0,3, 6}; 
    float bins_mt2[] = {200,600,  1500 };
    float bins_ht[] =  {450,1000, 2000};
    */

    float bins_nJets[] = {2,4,7,12};
    float bins_nBJets[] = {0,1,2,3,6}; 
    // float bins_mt2[] = {200,300,400,500,600 };
  
    //    float bins_mt2[] = {200,300,400, 600 };
      float bins_mt2[] = {200,300,400,500, 600, 800, 1000, 1500 };
    float bins_ht[] =  {450,575,1000,1500,2000};
  

    std::string cut =  "weight*(abs(Z_mass-91.19)<20 &&  Z_pt>180 )";

    std::string cut_el =  "weight*(abs(Z_mass-91.19)<20 &&  Z_pt>180 && Z_lepId==11 )";
    std::string cut_mu =  "weight*(abs(Z_mass-91.19)<20 &&  Z_pt>180 && Z_lepId==13 )";

    std::string cut_gamma =  "weight*( prompt==2 && ptGamma>180 )*1.27";
    //  std::string cut =  "weight*(mt2>200 && abs(Z_mass-91.19)<15 &&  Z_pt>200 && nJetHF30==0 )";
    //    std::string cut_gamma =  "weight*( prompt==2 &&ptGamma>180 && nJetHF30==0  )*1.27";
 

   //   std::string cut =  "weight*(mt2>200 && abs(Z_mass-91.19)<15 && lep_pt0>25 && lep_pt1>20 )";
   //  std::string cut =  "weight*(mt2>200 && abs(Z_mass-91.19)<15 && lep_pt0>25 && lep_pt1>20 &&( HLT_DoubleMu||HLT_DoubleEl) )";
    //std::string cut =  "weight*(ht>450 && abs(Z_mass-91.19)<25 )";
    //    std::string cut =  "weight*(abs(Z_mass-91.19)<20)";
    //    std::string cut =  "weight*(abs(Z_mass-91.19)<10 &&nBJets<2)";
    std::string cut_corr =  "weight*(abs(Z_mass-91.19)<10 )";
   //    std::string cut_gamma =  "weight*( prompt==2)*1.27";

    
    int size_mt2 = sizeof(bins_mt2)/sizeof(float)-1;
    int size_ht = sizeof(bins_ht)/sizeof(float)-1;
    int size_nJets = sizeof(bins_nJets)/sizeof(float)-1;
    int size_nBJets = sizeof(bins_nBJets)/sizeof(float)-1;
  
    /* 
    float bins_nJets2[] = {2,3,4,5,6,7,8,9,10,11,12};
    float bins_nBJets2[] = {0,1,2,3,4,5,6};
    float bins_mt22[] = {200,250,300,350, 400,450,500,550, 600,650,700, 800, 900,1000, 1100, 1200,1300, 1400, 1500, 1600,1700 };
    float bins_ht2[] =  {450,500,550, 600,650,700,750, 800, 900, 1000,1100, 1200, 1300,1400, 1500,1600};
 
    int size_mt2 = sizeof(bins_mt22)/sizeof(float)-1;
    int size_ht = sizeof(bins_ht2)/sizeof(float)-1;
    int size_nJets = sizeof(bins_nJets2)/sizeof(float)-1;
    int size_nBJets = sizeof(bins_nBJets2)/sizeof(float)-1;
    */ 
    float lumi = cfg.lumi();    
    
    //draw ratio also fills the ratio and yield estimates
    drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data,  zllY_mt2, thisRegion, cut, cut_gamma, lumi , "mt2" );
   
    drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data,  zllY_mt2, thisRegion, cut_el, cut_gamma, lumi , "mt2_el" );
   
    drawRatios( outputdir, bins_mt2, size_mt2 , "mt2",  zllG_mt2,   gamma_mc, gamma_data, purity,  zll_mc, zll_data,  zllY_mt2, thisRegion, cut_mu, cut_gamma, lumi , "mt2_mu" );
   
  
    /*
    drawRatios( outputdir, bins_ht, size_ht , "ht",   zllG_ht,   gamma_mc, gamma_data, purity,  zll_mc, zll_data,  zllY_ht, thisRegion, cut, cut_gamma, lumi );
    
    drawRatios( outputdir, bins_nJets, size_nJets , "nJets",   zllG_nJets,   gamma_mc, gamma_data, purity,  zll_mc, zll_data,  zllY_nJets, thisRegion, cut, cut_gamma, lumi );

    drawRatios( outputdir, bins_nBJets, size_nBJets , "nBJets",  zllG_nBJets,   gamma_mc, gamma_data, purity,  zll_mc, zll_data,  zllY_nBJets, thisRegion, cut, cut_gamma , lumi);
    */
    


     TH1D data_bgSub =  drawBGsubtraction(  cfg, zll_data, zll_mc,  bgYields , "mt2", "mt2", cut, 20, 200, 800, "M_{T2}", "GeV", scaleFactor );

     TH1D data_bgSub_nBJets =  drawBGsubtraction(  cfg, zll_data, zll_mc,  bgYields , "nBJets", "nBJets", cut, 8, 0, 8, "nBJets", "", scaleFactor );




 
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
  //  std::string outFile_ht = outputdir + "/zll_ht.root";
  //  std::string outFile_nJets = outputdir + "/zll_nJets.root";
  //  std::string outFile_nBJets = outputdir + "/zll_nBJets.root";
  /*
    zll_ratio->writeToFile( outFile );
    zll_ht->addToFile( outFile );
    zll_nJets->addToFile( outFile );
    zll_nBJets->addToFile( outFile );
    zll_yield->addToFile( outFile );
  */

  //  zll_ratio->writeToFile( outFile );
  /*
  zllY_mt2->writeToFile(outFile);

  zllG_ht->addToFile( outFile );
  zllG_nJets->addToFile( outFile );
  zllG_nBJets->addToFile( outFile );
  zllG_mt2->addToFile( outFile );
  */
  /*
  zll_ht->writeToFile( outFile_ht );
  zll_nJets->writeToFile( outFile_nJets );
  zll_nBJets->writeToFile( outFile_nBJets );
  */

  return 0;

}


void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string zll_sel, MT2Analysis<MT2Estimate>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma_mc, MT2Analysis<MT2EstimateTree>*  gamma_data, MT2Analysis<MT2EstimateSyst>*  purity, MT2Analysis<MT2EstimateTree>*  zll_mc,MT2Analysis<MT2EstimateTree>*  zll_data, MT2Analysis<MT2Estimate>*  zll_yield, const MT2Region thisRegion, std::string cut, std::string cut_gamma, float lumi, std::string saveName){
 
  std::vector<int> colors;
  colors.push_back(430); // other = zll 
  colors.push_back(401); // qcd
  colors.push_back(417); // w+jets
  colors.push_back(419); // z+jets
  colors.push_back(855); // top
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
  TH1F::AddDirectory(kTRUE);

  float bins[size+1]; for(unsigned int i=0; i<= size ; i++)      bins[i]=binss[i];
  float xMin = binss[0];   float xMax = binss[size];

  std::string gamma_sel = zll_sel;

  //THE TREES
  TTree *zll_data_tree =  zll_data->get(thisRegion)->tree;
  TTree *gamma_data_tree =  gamma_data->get(thisRegion)->tree;
  TTree *zll_mc_tree =  zll_mc->get(thisRegion)->tree;
  TTree *gamma_mc_tree =  gamma_mc->get(thisRegion)->tree;
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

  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();


  TH1D* h_mt2 = new TH1D("h_mt2","", size , bins);
  TH1D* g_mt2 = new TH1D("g_mt2","", size , bins);
  h_mt2->Sumw2(); g_mt2->Sumw2();
  
  zll_data_tree ->Project( "h_mt2" , zll_sel.c_str(), cut.c_str() );
  //will have to add the purity here at some point soon (for now just 0.95)
  gamma_data_tree->Project( "g_mt2", gamma_sel.c_str(),  "weight*0.92*(ptGamma>180 ) " );
  // gamma_data_tree->Project( "g_mt2", gamma_sel.c_str(),  "weight*0.92*(ptGamma>180 && nJetHF30==0) " );

  std::cout <<  h_mt2->GetMean() << std::endl;
  std::cout <<  g_mt2->GetMean() << std::endl;

  h_mt2->SetBinContent(size, h_mt2->GetBinContent(size) + h_mt2->GetBinContent(size+1));//adding overflow
  g_mt2->SetBinContent(size, g_mt2->GetBinContent(size) + g_mt2->GetBinContent(size+1));

  // we'll have to ignore the purity for a sec
  // Double_t x_tmp, p, p_errUp, p_errDown;	       
  // this_zinv_purity->GetPoint( binnie-1, x_tmp, p);
 

  int nBinss =  h_mt2->GetNbinsX();
  for(int binnie = 1; binnie <= nBinss; binnie++){

    Double_t x_tmp, p, p_errUp, p_errDown;	       
    this_zinv_purity->GetPoint( binnie-1, x_tmp, p);

    std::cout << "Purity = " << p << std::endl;

    double value = h_mt2->GetBinContent(binnie);
    h_mt2->SetBinError(binnie,sqrt(value));
    double value_g = g_mt2->GetBinContent(binnie) * p;
    g_mt2->SetBinContent(binnie, value_g);
    g_mt2->SetBinError(binnie,sqrt(value_g));
  }
  h_mt2->SetMarkerStyle(20);   h_mt2->SetMarkerSize(1.6); h_mt2->SetLineColor(kBlack);
  
 

  //Filling the YIELDS
  int yieldBins = h_mt2->GetNbinsX();
  for(int bi = 1; bi <= yieldBins ; bi++){
    double value = h_mt2->GetBinContent(bi);
    double err = h_mt2->GetBinError(bi);
    zll_yield->get(thisRegion)->yield->SetBinContent( bi,   value);
    zll_yield->get(thisRegion)->yield->SetBinError( bi,   err);

  }

  /*TH2D* h2_axes_yield = new TH2D("axes_yield", "", 10, 0, 1500, 10, 0.1, 699);
    h2_axes_yield->SetXTitle("M_{T2} [GeV]");
    h2_axes_yield->SetYTitle("Zll Yield");    h2_axes_yield->Draw();   
    h_mt2->Draw("p same");    labelTop->Draw("same");
    gPad->RedrawAxis();    //  gPad->SetLogy();
    //   canny->SaveAs( Form("%s/yield_mt2_%s.eps", outputdir.c_str(), thisRegion.getName().c_str() ) );
    //   canny->SaveAs( Form("%s/yield_mt2_%s.png", outputdir.c_str() , thisRegion.getName().c_str() ) );
    */


    TH1D* h_mt2_mc = new TH1D("h_mt2_mc","", size  , bins);
    TH1D* g_mt2_mc = new TH1D("g_mt2_mc","", size  , bins);
    h_mt2_mc->Sumw2(); g_mt2_mc->Sumw2();

    zll_mc_tree ->Project( "h_mt2_mc" , zll_sel.c_str(), cut.c_str() );
    gamma_mc_tree->Project( "g_mt2_mc", gamma_sel.c_str(), cut_gamma.c_str() );
  
    h_mt2_mc->SetBinContent(size, h_mt2_mc->GetBinContent(size) +h_mt2_mc->GetBinContent(size+1));
    g_mt2_mc->SetBinContent(size, g_mt2_mc->GetBinContent(size) +g_mt2_mc->GetBinContent(size+1));

    h_mt2->Divide(g_mt2);  
    h_mt2_mc->Divide(g_mt2_mc);

    TGraphErrors *gr_ratio = new TGraphErrors(0);
    for( unsigned int k=0; k< size; ++k){
     if(h_mt2->GetBinContent(k+1)<0.01){
      gr_ratio->SetPoint(k, 9000 , 900);
     }else{  
      gr_ratio->SetPoint(k, h_mt2->GetBinCenter(k+1), h_mt2->GetBinContent(k+1)); 
       //gr_ratio->SetPoint(k, meanX[k], h_mt2->GetBinContent(k+1));

      gr_ratio->SetPointError(k, h_mt2->GetBinWidth(k+1)/2. ,  h_mt2->GetBinError(k+1) );
    }
    }
    gr_ratio->SetMarkerSize(1.4);
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->SetLineColor(kBlack);
    gr_ratio->SetLineWidth(2);
    gr_ratio->SetMarkerColor(kBlack);


    //   h_mt2_mc->SetMarkerStyle(20); h_mt2_mc->SetMarkerSize(1.6);
    //   h_mt2_mc->SetMarkerColor(46); h_mt2_mc->SetLineColor(46);
    h_mt2_mc->SetLineColor(kBlue+1); h_mt2_mc->SetLineWidth(2);

    //  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.2 );
     TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 0.2 );
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
      //   if(i==0)
      //	regionText->AddText( "H_{T} > 180 GeV" );
      //  else
	regionText->AddText( niceNames2[i].c_str() );
      regionText->Draw("same");
    }

    h_mt2_mc->Draw("hist same");
 
    //    h_mt2->DrawClone("p same");
    gr_ratio->Draw("same P");
 
    labelTop->Draw("same");
    gPad->RedrawAxis();

    TLegend* legend = new TLegend( 0.7, 0.92-(2)*0.06, 0.9, 0.92 );
    legend->SetTextSize(0.04);
    // legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( h_mt2 ,"Data", "P" );
    legend->AddEntry( h_mt2_mc ,"Simulation", "L" );
    legend->Draw("same");


 

 
    gPad->RedrawAxis();

    canny->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
    pad2->SetTopMargin(0.10);
    pad2->SetBottomMargin(0.1);
    pad2->Draw();
    pad2->cd();
      
    h_mt2->Divide(h_mt2_mc); //the ratio is now in h_mt2

    TGraphErrors *gr_ratioD = new TGraphErrors(0);
    for( unsigned int k=0; k< size; ++k){  
      gr_ratioD->SetPoint(k, h_mt2->GetBinCenter(k+1), h_mt2->GetBinContent(k+1));
      // gr_ratioD->SetPoint(k, meanX[k], h_mt2->GetBinContent(k+1));
      gr_ratioD->SetPointError(k, 0, h_mt2->GetBinError(k+1));
      // gr_ratioD->SetPointError(k, h_mt2->GetBinWidth(k+1)/2., h_mt2->GetBinError(k+1));
    } 
    gr_ratioD->SetMarkerSize(1.4);
    gr_ratioD->SetMarkerStyle(20);
    gr_ratioD->SetLineColor(kBlack);
    gr_ratioD->SetLineWidth(2);
    gr_ratioD->SetMarkerColor(kBlack);


    TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , 0., 2.0 );
    h2_axes_rat->SetYTitle("Data / MC");
 
    h2_axes_rat->GetXaxis()->SetTitleSize(0.2);
    h2_axes_rat->GetXaxis()->SetTitleOffset(5);
    h2_axes_rat->GetXaxis()->SetLabelSize(0.00);
    h2_axes_rat->GetXaxis()->SetTickLength(0.09);
    // h2_axes_rat->GetYaxis()->SetTitle("ratio");
    h2_axes_rat->GetYaxis()->SetNdivisions(5,5,0);
    h2_axes_rat->GetYaxis()->SetTitleSize(0.2);
    h2_axes_rat->GetYaxis()->SetTitleOffset(0.34);
    h2_axes_rat->GetYaxis()->SetLabelSize(0.17);
  
  
    h2_axes_rat->Draw();
    // h_mt2->Draw("p same");
    gr_ratioD->Draw("same P");

    TLine* line = new TLine(xMin,1 , xMax, 1);
    line->SetLineColor(kBlack);
    line->Draw("same");

    gPad->RedrawAxis();

    canny->cd();
    
    canny->SaveAs( Form("%s/%s_ratios_%s.eps", fullPath.c_str(), saveName.c_str(),  thisRegion.getName().c_str() ) );
    canny->SaveAs( Form("%s/%s_ratios_%s.png", fullPath.c_str(), saveName.c_str(),thisRegion.getName().c_str() ) );
    canny->SaveAs( Form("%s/%s_ratios_%s.pdf", fullPath.c_str(), saveName.c_str(),thisRegion.getName().c_str() ) );
    
    /*
    canny->SaveAs( Form("%s/%s_ratios_wHFveto_%s.eps", fullPath.c_str(), zll_sel.c_str(),  thisRegion.getName().c_str() ) );
    canny->SaveAs( Form("%s/%s_ratios_wHFveto_%s.png", fullPath.c_str(), zll_sel.c_str(),thisRegion.getName().c_str() ) );
    canny->SaveAs( Form("%s/%s_ratios_wHFveto_%s.pdf", fullPath.c_str(), zll_sel.c_str(),thisRegion.getName().c_str() ) );
    */

    //Filling the YIELDS
    int yieldBins_mt2 = h_mt2->GetNbinsX();
    for(int bi = 1; bi <= yieldBins_mt2 ; bi++){
      double value = h_mt2->GetBinContent(bi);
      double err = h_mt2->GetBinError(bi);
      zll_ratio->get(thisRegion)->yield->SetBinContent( bi,   value);
      zll_ratio->get(thisRegion)->yield->SetBinError( bi,   err);
    }
 
    delete h_mt2;   delete h_mt2_mc;    delete g_mt2; delete g_mt2_mc;
    delete h2_axes; delete h2_axes_rat;
    delete gr_ratio;
    delete canny;

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

  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
 
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






















