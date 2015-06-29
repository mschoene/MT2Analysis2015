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

//float lumi = 0.1;
float lumi = 4.;


void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string zll_sel,  std::string gamma_sel, MT2Analysis<MT2Estimate>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma, MT2Analysis<MT2EstimateTree>*  Zll, MT2Analysis<MT2Estimate>*  zll_yield, const MT2Region thisRegion, std::string cut, std::string cut_gamma);


//void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string name ,  MT2Analysis<MT2Estimate>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma, MT2Analysis<MT2EstimateTree>*  Zll, MT2Analysis<MT2Estimate>*  zll_yield, const MT2Region thisRegion, std::string cut);


int main(int argc, char* argv[]){

  std::string regionsSet = "zurich";
  if( argc>2 ) {
    regionsSet = std::string(argv[2]);
  }

  std::cout << "-> Using regions: " << regionsSet << std::endl;

  if( argc<2 ) {
    std::cout << "USAGE: ./coputeZllGammaRatio [configFileName] regionSet" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  std::string configFileName(argv[1]);
  MT2Config cfg("cfgs/" + configFileName + ".txt");
  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat"; 
  std::string samples = cfg.mcSamples();

  std::string outputdir = "ZllGamma_Ratio_" + cfg.mcSamples();
  //  std::string outputdir = "ZllGamma_Ratio_" + regionsSet;
  //  std::string outputdir = "Zll_" + configFileName;
  double intpart;
  double fracpart = modf(lumi, &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("_%.0ffb", intpart ) );
  outputdir += suffix;
  
  system(Form("mkdir -p %s", outputdir.c_str()));
 
  std::cout << "-> Using regions: " << regionsSet << std::endl;



  gStyle->SetOptTitle(0);
  MT2DrawTools::setStyle();
  /*
  TGaxis::SetMaxDigits(3);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetPalette(55,0);
  */
 
  std::string gammaControlRegionDir = "GammaControlRegion_" + samples + "_" + regionsSet + suffix;

  // std::string gammaControlRegionDir = "GammaControlRegion_" + samples + "_" + regionsSet + "_10fb";

 
  MT2Analysis<MT2EstimateTree>* gamma = MT2Analysis<MT2EstimateTree>::readFromFile(gammaControlRegionDir + "/data.root", "gammaCRtree");
  if( gamma==0 ) {
    std::cout << "-> Please run gammaControlRegion first. I need to get the gammaCR yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(193);
  }


  std::string ZllDir = "Zll_CR_" + samples + "_" + regionsSet;
  //  std::string ZllDir = "Zll_" + samples + "_" + regionsSet;

  //  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb/Zll_analyses.root", ZllDir.c_str(), lumi), "DYJets");
  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s%s/Zll_analyses.root", ZllDir.c_str(), suffix.c_str()), "DYJets");
  if( Zll==0 ) {
    std::cout << "-> Please run computeZinvFromZll first. I need to get the Z->ll MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }


  // MT2Analysis<MT2Estimate>* zll_zurich = new MT2Analysis<MT2Estimate>("zll_zurich", "zurich");
  MT2Analysis<MT2Estimate>* zll_ratio = new MT2Analysis<MT2Estimate>( "zll_ratio", "inclusive" );
  MT2Analysis<MT2Estimate>* zll_yield = new MT2Analysis<MT2Estimate>( "zll_yield", "inclusive" );
  //  MT2Analysis<MT2Estimate>* zll_ratio = new MT2Analysis<MT2Estimate>( "zll_ratio", "zurich" );

  //YIELDS
  MT2Analysis<MT2Estimate>* zllY_pt = new MT2Analysis<MT2Estimate>( "zllY_pt", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllY_mt2 = new MT2Analysis<MT2Estimate>( "zllY_mt2", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllY_ht = new MT2Analysis<MT2Estimate>( "zllY_ht", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllY_nJets = new MT2Analysis<MT2Estimate>( "zllY_nJets", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllY_nBJets = new MT2Analysis<MT2Estimate>( "zllY_nBJets", "13TeV_inclusive" );

  //RATIOS
  MT2Analysis<MT2Estimate>* zllG_pt = new MT2Analysis<MT2Estimate>( "zllG_pt", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllG_mt2 = new MT2Analysis<MT2Estimate>( "zllG_mt2", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllG_ht = new MT2Analysis<MT2Estimate>( "zllG_ht", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllG_nJets = new MT2Analysis<MT2Estimate>( "zllG_nJets", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllG_nBJets = new MT2Analysis<MT2Estimate>( "zllG_nBJets", "13TeV_inclusive" );
  

  int counter_total_bins = 0;
  int counter_empty_bins = 0;

  //Loop over regions

  std::set<MT2Region> MT2Regions = zll_ratio->getRegions();
  //std::set<MT2Region> MT2Regions = Zll->getRegions();
  std::set<MT2Region> Inclusive = Zll->getRegions();
  MT2Region RegIncl( (*Inclusive.begin()) );

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  TH1F::AddDirectory(kTRUE);

  //THE TREES
  TTree *zllT =  Zll->get(RegIncl)->tree;
  TTree *gammaT =  gamma->get(RegIncl)->tree;

  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
 

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
    // std::set<MT2Region>::iterator iMT2 = MT2Regions.begin();  
    std::vector<std::string> niceNames = thisRegion.getNiceNames();



    float bins_nJets[] = {2,4,7,12};
    float bins_nBJets[] = {0,1,2,3,6};
    //in MT2
    float bins_mt2[] = {200,300,400,500, 600, 800, 1000, 1500 };
    //in HT
    float bins_ht[] =  {450,575,1000,1500,2000};

    //    float bins_mt2[] = {200,300,400,500, 600, 800, 1000, 1500 , 1900};
    std::string cut =  "weight*(abs(Z_mass-91.19)<10 &&nBJets<2)";
     std::string cut_gamma =  "weight*(prompt==2)";
  
     drawRatios( outputdir, bins_mt2, sizeof(bins_mt2)/sizeof(float)-1 , "zll_mt2", "gamma_mt2" ,  zllG_mt2,   gamma,  Zll,  zllY_mt2, thisRegion, cut, cut_gamma );
   
     drawRatios( outputdir, bins_ht, sizeof(bins_ht)/sizeof(float)-1 , "zll_ht", "gamma_ht" ,  zllG_ht,   gamma,  Zll,  zllY_ht, thisRegion, cut, cut_gamma );
    
     drawRatios( outputdir, bins_nJets, sizeof(bins_nJets)/sizeof(float)-1 , "nJets", "gamma_nJets" ,  zllG_nJets,   gamma,  Zll,  zllY_nJets, thisRegion, cut, cut_gamma );
     drawRatios( outputdir, bins_nBJets, sizeof(bins_nBJets)/sizeof(float)-1 , "nBJets", "gamma_nBJets" ,  zllG_nBJets,   gamma,  Zll,  zllY_nBJets, thisRegion, cut, cut_gamma );
 



     /*

    //void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string name, MT2Analysis<MT2EstimateTree>  zll_ratio,  MT2Analysis<MT2EstimateTree>  gamma, MT2Analysis<MT2EstimateTree>  Zll, MT2Analysis<MT2EstimateTree>  zll_yield, const MT2Region thisRegion, std::string cut);
 

    //Filling the estimate
    TH1D* thisZll  = zll_ratio ->get(*iMT2)->yield;
    TH1D* thisGamma = zll_ratio->get(*iMT2)->yield;
    

    TH1D* g_mt2 = new TH1D(*thisGamma); // so that it gets the same binning
    g_mt2->SetName("g_mt2");
    TH1D* h_mt2  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2->SetName("h_mt2");

    zllT ->Project( "h_mt2" , "zll_mt2", "weight*(( abs(Z_mass-91.19)<10 && nBJets<2))");
    gammaT->Project( "g_mt2", "gamma_mt2", "weight*(prompt==2)");
  
 

    int nBins =  thisZll ->GetNbinsX();
    counter_total_bins += nBins;

    for(int binnie = 1; binnie < nBins+1; ++binnie ){
      if( thisZll->GetBinContent(binnie) < 1 ) {
	std::cout << "Yield = "<< thisZll->GetBinContent(binnie) <<" in " << thisRegion.getName() << " at mt2 = " << thisZll->GetBinCenter(binnie)   << std::endl;
	counter_empty_bins = counter_empty_bins + 1;
      }
    }


    int nBinss =  h_mt2->GetNbinsX();
    for(int binnie = 1; binnie <= nBinss; binnie++){
      double value = h_mt2->GetBinContent(binnie);
      h_mt2->SetBinError(binnie,sqrt(value));
      double value_g = g_mt2->GetBinContent(binnie);
      g_mt2->SetBinError(binnie,sqrt(value_g));
    }
    h_mt2->SetMarkerStyle(20);   h_mt2->SetMarkerSize(1.6); h_mt2->SetLineColor(kBlack);
  


    //Filling the YIELDS
    int yieldBins = h_mt2->GetNbinsX();
    for(int bi = 1; bi <= yieldBins+1 ; bi++){
      double value = h_mt2->GetBinContent(bi);
      double err = h_mt2->GetBinError(bi);
      zll_yield->get(*iMT2)->yield->SetBinContent( bi,   value);
      zll_yield->get(*iMT2)->yield->SetBinError( bi,   err);
    }



    TH2D* h2_axes_yield = new TH2D("axes_yield", "", 10, 0, 1500, 10, 0.1, 699);
    h2_axes_yield->SetXTitle("M_{T2} [GeV]");
    h2_axes_yield->SetYTitle("Zll Yield");
    h2_axes_yield->Draw();
    
    h_mt2->Draw("p same");
    labelTop->Draw("same");
    gPad->RedrawAxis();
    //  gPad->SetLogy();
 
    c1->SaveAs( Form("%s/yield_mt2_%s.eps", outputdir.c_str(), thisRegion.getName().c_str() ) );
    c1->SaveAs( Form("%s/yield_mt2_%s.png", outputdir.c_str() , thisRegion.getName().c_str() ) );
    


    //Checking the projection//
    TH2D* h2_axes_yield_b= new TH2D("axes_yield_b", "", 10, 0, 1500, 10, 0.0005,0.9);
    h2_axes_yield_b->SetXTitle("M_{T2} [GeV]");
    h2_axes_yield_b->SetYTitle("Zll Yield");
    h2_axes_yield_b->Draw();

    TH1D* h_mt2_b0  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_b0->SetName("h_mt2_b0");
    zllT ->Project( "h_mt2_b0" , "zll_mt2", "weight*((abs(Z_mass-91.19)<10 && nBJets<2) && nBJets==0)");
    h_mt2_b0->SetMarkerStyle(20);   h_mt2_b0->SetMarkerSize(1.6);
    h_mt2_b0->SetLineColor(46); h_mt2_b0->SetMarkerColor(46);
    h_mt2_b0->Scale(1./h_mt2_b0->Integral());    h_mt2_b0->Draw("p same");

    TH1D* h_mt2_b1  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_b1->SetName("h_mt2_b1");
    zllT ->Project( "h_mt2_b1" , "zll_mt2", "weight*(abs(Z_mass-91.19)<10 && nBJets<2  && nBJets==1)");
    h_mt2_b1->SetMarkerStyle(20);   h_mt2_b1->SetMarkerSize(1.6); 
    h_mt2_b1->SetLineColor(29); h_mt2_b1->SetMarkerColor(29);
    h_mt2_b1->Scale(1./h_mt2_b1->Integral());    h_mt2_b1->Draw("p same");

    TH1D* h_mt2_b2  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_b2->SetName("h_mt2_b2");
    zllT ->Project( "h_mt2_b2" , "zll_mt2", "weight*(abs(Z_mass-91.19)<10 && nBJets<2  && nBJets==2)");
    h_mt2_b2->SetMarkerStyle(20);   h_mt2_b2->SetMarkerSize(1.6);
    h_mt2_b2->SetLineColor(38);  h_mt2_b2->SetMarkerColor(38);
    h_mt2_b2->Scale(1./h_mt2_b2->Integral());    h_mt2_b2->Draw("p same");

    TH1D* h_mt2_b3  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_b3->SetName("h_mt2_b3");
    zllT ->Project( "h_mt2_b3" , "zll_mt2", "weight*(abs(Z_mass-91.19)<10 && nBJets<2  && nBJets>2)");
    h_mt2_b3->SetMarkerStyle(20);   h_mt2_b3->SetMarkerSize(1.6);
    h_mt2_b3->SetLineColor(42);  h_mt2_b3->SetMarkerColor(42);
    h_mt2_b3->Scale(1./h_mt2_b3->Integral());    h_mt2_b3->Draw("p same");

    TLegend* legend_b = new TLegend( 0.69, 0.9-(5)*0.06, 0.9, 0.9 );
    legend_b->SetTextSize(0.038);
    legend_b->SetTextFont(42);
    legend_b->SetFillColor(0);
    legend_b->AddEntry( h_mt2_b0 ,"b=0", "P" );    legend_b->AddEntry( h_mt2_b1 ,"b=1", "P" );
    legend_b->AddEntry( h_mt2_b2 ,"b=2", "P" );    legend_b->AddEntry( h_mt2_b3 ,"b=3+", "P" );
    legend_b->Draw("same");    labelTop->Draw("same");    gPad->RedrawAxis();

    c1->SaveAs( Form("%s/yield_mt2_b_%s.eps", outputdir.c_str(), thisRegion.getName().c_str() ) );
    c1->SaveAs( Form("%s/yield_mt2_b_%s.png", outputdir.c_str() , thisRegion.getName().c_str() ) );
    


    h2_axes_yield_b->Draw();
    TH1D* h_mt2_j0  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_j0->SetName("h_mt2_j0");
    zllT ->Project( "h_mt2_j0" , "zll_mt2", "weight*((abs(Z_mass-91.19)<10 && nBJets<2) && nJets<4)");
    h_mt2_j0->SetMarkerStyle(20);   h_mt2_j0->SetMarkerSize(1.6);
    h_mt2_j0->SetLineColor(46); h_mt2_j0->SetMarkerColor(46);
    h_mt2_j0->Scale(1./h_mt2_j0->Integral());    h_mt2_j0->Draw("p same");

    TH1D* h_mt2_j1  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_j1->SetName("h_mt2_j1");
    zllT ->Project( "h_mt2_j1" , "zll_mt2", "weight*(abs(Z_mass-91.19)<10 && nBJets<2  && nJets>3 && nJets<=6)");
    h_mt2_j1->SetMarkerStyle(20);   h_mt2_j1->SetMarkerSize(1.6); 
    h_mt2_j1->SetLineColor(29); h_mt2_j1->SetMarkerColor(29);
    h_mt2_j1->Scale(1./h_mt2_j1->Integral());    h_mt2_j1->Draw("p same");

    TH1D* h_mt2_j2  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_j2->SetName("h_mt2_j2");
    zllT ->Project( "h_mt2_j2" , "zll_mt2", "weight*(abs(Z_mass-91.19)<10 && nBJets<2  && nJets>=7)");
    h_mt2_j2->SetMarkerStyle(20);   h_mt2_j2->SetMarkerSize(1.6);
    h_mt2_j2->SetLineColor(38);  h_mt2_j2->SetMarkerColor(38);
    h_mt2_j2->Scale(1./h_mt2_j2->Integral());    h_mt2_j2->Draw("p same");


    TLegend* legend_j = new TLegend( 0.69, 0.9-(5)*0.06, 0.9, 0.9 );
    legend_j->SetTextSize(0.038);
    legend_j->SetTextFont(42);
    legend_j->SetFillColor(0);
    legend_j->AddEntry( h_mt2_j0 ,"jets=2,3", "P" );    legend_j->AddEntry( h_mt2_j1 ,"jets=4,5,6", "P" );
    legend_j->AddEntry( h_mt2_j2 ,"jets=7+", "P" );   
    legend_j->Draw("same");    labelTop->Draw("same");    gPad->RedrawAxis();

    c1->SaveAs( Form("%s/yield_mt2_j_%s.eps", outputdir.c_str(), thisRegion.getName().c_str() ) );
    c1->SaveAs( Form("%s/yield_mt2_j_%s.png", outputdir.c_str() , thisRegion.getName().c_str() ) );
    

    //HT
    h2_axes_yield_b->Draw();

    TH1D* h_mt2_ht0  = new TH1D(*thisZll);   h_mt2_ht0->SetName("h_mt2_ht0");
    zllT ->Project( "h_mt2_ht0" , "zll_mt2", "weight*((abs(Z_mass-91.19)<10 && nBJets<2)&& zll_ht>450&&zll_ht<575 )");
    h_mt2_ht0->SetMarkerStyle(20);   h_mt2_ht0->SetMarkerSize(1.6);
    h_mt2_ht0->SetLineColor(46); h_mt2_ht0->SetMarkerColor(46);
    h_mt2_ht0->Scale(1./h_mt2_ht0->Integral());    h_mt2_ht0->Draw("p same");

    TH1D* h_mt2_ht1  = new TH1D(*thisZll);    h_mt2_ht1->SetName("h_mt2_ht1");
    zllT ->Project( "h_mt2_ht1" , "zll_mt2", "weight*(abs(Z_mass-91.19)<10 && nBJets<2 && zll_ht>575&&zll_ht<1000)");
    h_mt2_ht1->SetMarkerStyle(20);   h_mt2_ht1->SetMarkerSize(1.6); 
    h_mt2_ht1->SetLineColor(29); h_mt2_ht1->SetMarkerColor(29);
    h_mt2_ht1->Scale(1./h_mt2_ht1->Integral());    h_mt2_ht1->Draw("p same");

    TH1D* h_mt2_ht2  = new TH1D(*thisZll);  h_mt2_ht2->SetName("h_mt2_ht2");
    zllT ->Project( "h_mt2_ht2" , "zll_mt2", "weight*(abs(Z_mass-91.19)<10 && nBJets<2&& zll_ht>1000&&zll_ht<1500 )");
    h_mt2_ht2->SetMarkerStyle(20);   h_mt2_ht2->SetMarkerSize(1.6);
    h_mt2_ht2->SetLineColor(38);  h_mt2_ht2->SetMarkerColor(38);
    h_mt2_ht2->Scale(1./h_mt2_ht2->Integral());    h_mt2_ht2->Draw("p same");

    TH1D* h_mt2_ht3  = new TH1D(*thisZll);    h_mt2_ht3->SetName("h_mt2_ht3");
    zllT ->Project( "h_mt2_ht3" , "zll_mt2", "weight*(abs(Z_mass-91.19)<10 && nBJets<2&& zll_ht>1500 )");
    h_mt2_ht3->SetMarkerStyle(20);   h_mt2_ht3->SetMarkerSize(1.6);
    h_mt2_ht3->SetLineColor(42);  h_mt2_ht3->SetMarkerColor(42);
    h_mt2_ht3->Scale(1./h_mt2_ht3->Integral());    h_mt2_ht3->Draw("p same");

    TLegend* legend_ht = new TLegend( 0.69, 0.9-(5)*0.06, 0.9, 0.9 );
    legend_ht->SetTextSize(0.038);
    legend_ht->SetTextFont(42);
    legend_ht->SetFillColor(0);
    legend_ht->AddEntry( h_mt2_ht0 ,"Low HT", "P" );    legend_ht->AddEntry( h_mt2_ht1 ,"Medium HT", "P" );
    legend_ht->AddEntry( h_mt2_ht2 ,"High HT", "P" );    legend_ht->AddEntry( h_mt2_ht3 ,"Extreme HT", "P" );
    legend_ht->Draw("same");    labelTop->Draw("same");    gPad->RedrawAxis();

    c1->SaveAs( Form("%s/yield_mt2_ht_%s.eps", outputdir.c_str(), thisRegion.getName().c_str() ) );
    c1->SaveAs( Form("%s/yield_mt2_ht_%s.png", outputdir.c_str() , thisRegion.getName().c_str() ) );







    h_mt2->Divide(g_mt2);


    //Histos of MT2, ZLL
    //   TH1D* h_mt2 = new TH1D("h_mt2","h_mt2",10, 0, 1500);
    //   zllT->Project("h_mt2", "zll_mt2","weight*(( abs(Z_mass-91.19)<10 && nBJets<2))");
    TH1D* h_mt2_450  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_450->SetName("h_mt2_450");
    //   TH1D* h_mt2_450 = new TH1D("h_mt2_450","h_mt2_450",10, 0, 1500); 
    zllT->Project("h_mt2_450", "zll_mt2","weight*(zll_ht>450&&zll_ht<575 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
    TH1D* h_mt2_575  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_575->SetName("h_mt2_575");
    //   TH1D* h_mt2_575 = new TH1D("h_mt2_575","h_mt2_575",10, 0, 1500); 
    zllT->Project("h_mt2_575", "zll_mt2","weight*(zll_ht>575&&zll_ht<1000 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
    //1000 < HT < 1500
    TH1D* h_mt2_1000  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_1000->SetName("h_mt2_1000");
    //    TH1D* h_mt2_1000 = new TH1D("h_mt2_1000","h_mt2_1000",10, 0, 1500); 
    zllT->Project("h_mt2_1000", "zll_mt2","weight*(zll_ht>1000&&zll_ht<1500 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
    //       HT > 1500
    TH1D* h_mt2_1500  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_1500->SetName("h_mt2_1500");
    //  TH1D* h_mt2_1500 = new TH1D("h_mt2_1500","h_mt2_1500",10, 0, 1500); 
    zllT->Project("h_mt2_1500", "zll_mt2","weight*(zll_ht>1500 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
    //    h_mt2_450->Sumw2();h_mt2_575->Sumw2();h_mt2_1000->Sumw2(); h_mt2_1500->Sumw2();
    //  h_mt2->Sumw2();

    //Histos of MT2, gamma
    //Inclusive
    //  TH1D* g_mt2 = new TH1D("g_mt2","",10, 0, 1500); 
    //   gammaT->Project("g_mt2","gamma_mt2","weight*(prompt==2)"); 
  
      // 450 < HT < 575
    TH1D* g_mt2_450  = new TH1D(*thisZll); // so that it gets the same binning
    g_mt2_450->SetName("g_mt2_450");
    //    TH1D* g_mt2_450 = new TH1D("g_mt2_450","",10, 0, 1500); 
    gammaT->Project("g_mt2_450","gamma_mt2","weight*(prompt==2&&gamma_ht>450 && gamma_ht<575)"); 
    // 575 < HT < 1000
    TH1D* g_mt2_575  = new TH1D(*thisZll); // so that it gets the same binning
    g_mt2_575->SetName("g_mt2_575");
    //   TH1D* g_mt2_575 = new TH1D("g_mt2_575","",10, 0, 1500); 
    gammaT->Project("g_mt2_575","gamma_mt2","weight*(prompt==2&&gamma_ht>575 && gamma_ht<1000)"); 
    // 1000 < HT < 1500
    TH1D* g_mt2_1000  = new TH1D(*thisZll); // so that it gets the same binning
    g_mt2_1000->SetName("g_mt2_1000");
    //    TH1D* g_mt2_1000 = new TH1D("g_mt2_1000","",10, 0, 1500); 
    gammaT->Project("g_mt2_1000","gamma_mt2","weight*(prompt==2&&gamma_ht>1000 && gamma_ht<1500)"); 
    // HT > 1500
    TH1D* g_mt2_1500  = new TH1D(*thisZll); // so that it gets the same binning
    g_mt2_1500->SetName("g_mt2_1500");
    //   TH1D* g_mt2_1500 = new TH1D("g_mt2_1500","",10, 0, 1500); 
    gammaT->Project("g_mt2_1500","gamma_mt2","weight*(prompt==2&&gamma_ht>1500)"); 

    //    g_mt2_450->Sumw2();g_mt2_575->Sumw2();g_mt2_1000->Sumw2();g_mt2_1500->Sumw2();
  

    for(int binnie = 1; binnie <= nBinss; binnie++){
      double value_450 = h_mt2_450->GetBinContent(binnie);
      h_mt2_450->SetBinError(binnie,sqrt(value_450));
      double value_g_450 = g_mt2_450->GetBinContent(binnie);
      g_mt2_450->SetBinError(binnie,sqrt(value_g_450));

      double value_575 = h_mt2_575->GetBinContent(binnie);
      h_mt2_575->SetBinError(binnie,sqrt(value_575));
      double value_g_575 = g_mt2_575->GetBinContent(binnie);
      g_mt2_575->SetBinError(binnie,sqrt(value_g_575));

      double value_1000 = h_mt2_1000->GetBinContent(binnie);
      h_mt2_1000->SetBinError(binnie,sqrt(value_1000));
      double value_g_1000 = g_mt2_1000->GetBinContent(binnie);
      g_mt2_1000->SetBinError(binnie,sqrt(value_g_1000));

      double value_1500 = h_mt2_1500->GetBinContent(binnie);
      h_mt2_1500->SetBinError(binnie,sqrt(value_1500));
      double value_g_1500 = g_mt2_1500->GetBinContent(binnie);
      g_mt2_1500->SetBinError(binnie,sqrt(value_g_1500));
    }  

    //   g_mt2->Sumw2();

    //    h_mt2->Divide(g_mt2);  
   
    h_mt2_450->Divide(g_mt2_450);
    h_mt2_575->Divide(g_mt2_575);
    h_mt2_1000->Divide(g_mt2_1000);
    h_mt2_1500->Divide(g_mt2_1500);
    h_mt2_450->SetMarkerStyle(20);   h_mt2_450->SetMarkerSize(1.6); 
    h_mt2_450->SetMarkerColor(46);  h_mt2_450->SetLineColor(46);
    h_mt2_575->SetMarkerStyle(20);   h_mt2_575->SetMarkerSize(1.6); 
    h_mt2_575->SetMarkerColor(29); h_mt2_575->SetLineColor(29); 
    h_mt2_1000->SetMarkerStyle(20);   h_mt2_1000->SetMarkerSize(1.6);
    h_mt2_1000->SetMarkerColor(38); h_mt2_1000->SetLineColor(38); 
    h_mt2_1500->SetMarkerStyle(20);   h_mt2_1500->SetMarkerSize(1.6); 
    h_mt2_1500->SetMarkerColor(42); h_mt2_1500->SetLineColor(42); 
    

    TH2D* h2_axes = new TH2D("axes", "", 10, 0, 1500, 10, 0., 0.5 );
    h2_axes->SetXTitle("M_{T2} [GeV]");
    h2_axes->SetYTitle("Zll / #gamma Ratio");
    h2_axes->Draw();


    TLegend* legend = new TLegend( 0.2, 0.9-(5)*0.06, 0.5, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( h_mt2 ,"HT inclusive", "P" );
    legend->AddEntry( h_mt2_450 ,"450<HT<575", "P" );
    legend->AddEntry( h_mt2_575 ,"575<HT<1000", "P" );
    legend->AddEntry( h_mt2_1000 ,"1000<HT<1500", "P" );
    legend->AddEntry( h_mt2_1500 ,"HT>1500", "P" );
    h_mt2_450->Draw("p same");
    h_mt2_575->Draw("p same");
    h_mt2_1000->Draw("p same");
    h_mt2_1500->Draw("p same");
    
    h_mt2->Draw("p same");

    labelTop->Draw("same");
    legend->Draw("same");

    gPad->RedrawAxis();

 
    c1->SaveAs( Form("%s/mt2_%s.eps", outputdir.c_str(), thisRegion.getName().c_str() ) );
    c1->SaveAs( Form("%s/mt2_%s.png", outputdir.c_str() , thisRegion.getName().c_str() ) );




    TH1D* g_mt2_mc = new TH1D(*thisGamma); // so that it gets the same binning
    g_mt2_mc->SetName("g_mt2_mc");
    TH1D* h_mt2_mc  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_mc->SetName("h_mt2_mc");
    h_mt2_mc->Sumw2(); g_mt2_mc->Sumw2();

    zllT ->Project( "h_mt2_mc" , "zll_mt2", "weight*(( abs(Z_mass-91.19)<10 && nBJets<2))" );
    gammaT->Project( "g_mt2_mc", "gamma_mt2", "weight*(prompt==2)" );
  

    h_mt2_mc->Divide(g_mt2_mc);// h_mt2_mc->SetMarkerStyle(20); h_mt2_mc->SetMarkerSize(1.6);
    h_mt2->Divide(h_mt2_mc); //the ratio is now in h_mt2

    TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, 0, 1500, 10, 0.,2  );
    h2_axes_rat->SetXTitle("M_{T2} [GeV]");
    h2_axes_rat->SetYTitle("(Zll/ #gamma)_{Data} / (Zll/ #gamma)_{MC}");
    h2_axes_rat->Draw();
    labelTop->Draw("same");
    h_mt2->Draw("p same");

    c1->SaveAs( Form("%s/mt2_ratio_%s.eps", outputdir.c_str(), thisRegion.getName().c_str() ) );
    c1->SaveAs( Form("%s/mt2_ratio_%s.png", outputdir.c_str(), thisRegion.getName().c_str() ) );




    //Filling the YIELDS
    int yieldBins_mt2 = h_mt2->GetNbinsX();
    for(int bi = 1; bi <= yieldBins_mt2+1 ; bi++){
      double value = h_mt2->GetBinContent(bi);
      double err = h_mt2->GetBinError(bi);
      zll_ratio->get(*iMT2)->yield->SetBinContent( bi,   value);
      zll_ratio->get(*iMT2)->yield->SetBinError( bi,   err);
    }
 
    delete h_mt2;   delete h_mt2_mc;    delete g_mt2; delete g_mt2_mc;
    delete h2_axes; delete h2_axes_rat;

    //   delete h_mt2_450; delete h_mt2_575;    delete h_mt2_1000;    delete h_mt2_1500;
    //   delete g_mt2_450; delete g_mt2_575;    delete g_mt2_1000;    delete g_mt2_1500;

  }





  //THE INCLUSIVE VARIABLES////


  //FOR THE HT plot
  float bins_ht[] = {450,575,1000,1500,2000};
  int binNum_ht =  sizeof(bins_ht)/sizeof(float) -1;

  TH1D* h_ht2 = new TH1D("h_ht2","h_ht2",binNum_ht, bins_ht);
  zllT->Project("h_ht2", "zll_ht","weight*(abs(Z_mass-91.19)<10 && nBJets<2)");
  TH1D* g_ht2 = new TH1D("g_ht2","g_ht2",binNum_ht, bins_ht);
    
  gammaT->Project("g_ht2", "gamma_ht","weight*(prompt==2)");
    
  int nBins2 =  h_ht2->GetNbinsX();
  for(int binnie = 1; binnie <= nBins2+1; binnie++){
    double value = h_ht2->GetBinContent(binnie);
    h_ht2->SetBinError(binnie,sqrt(value));
    double value_g = g_ht2->GetBinContent(binnie);
    g_ht2->SetBinError(binnie,sqrt(value_g));
  }
  h_ht2->SetMarkerStyle(20);      h_ht2->SetMarkerSize(1.6);
  h_ht2->SetMarkerColor(kBlack);  h_ht2->SetLineColor(kBlack);
 

  TH2D* h2_axes_yield_ht = new TH2D("axes_yield_ht", "", 10, 400, 2000, 10, 1,900);
  h2_axes_yield_ht->SetXTitle("H_{T} [GeV]");
  h2_axes_yield_ht->SetYTitle("Zll Yield");
  h2_axes_yield_ht->Draw();
    
  h_ht2->Draw("p same");
  labelTop->Draw("same");
  gPad->RedrawAxis();
  //  gPad->SetLogy();

  c1->SaveAs( Form("%s/yield_ht_%s.png", outputdir.c_str() , RegIncl.getName().c_str()) );
  c1->SaveAs( Form("%s/yield_ht_%s.eps", outputdir.c_str() , RegIncl.getName().c_str()) );
    



  h_ht2->Divide(g_ht2); 

  //The histo with the wrong errors
  TH1D* h_ht = new TH1D("h_ht","h_ht",binNum_ht, bins_ht);
  TH1D* g_ht = new TH1D("g_ht","g_ht",binNum_ht, bins_ht);  
  zllT->Project("h_ht", "zll_ht","weight*(abs(Z_mass-91.19)<10 && nBJets<2)");
  gammaT->Project("g_ht", "gamma_ht","weight*(prompt==2)");
 
  h_ht->Sumw2(); g_ht->Sumw2();
  h_ht->Divide(g_ht);    h_ht->SetMarkerStyle(20);    h_ht->SetMarkerSize(1.6);
 

  TH1D* h_ht_mc = new TH1D("h_ht_mc","h_ht_mc",binNum_ht, bins_ht);
  TH1D* g_ht_mc = new TH1D("g_ht_mc","g_ht_mc",binNum_ht, bins_ht);
 
  //Double ratioooooo
  //    TH1D* h_ht_mc = new TH1D("h_ht_mc","h_ht_mc",20, 0, 2000); 
  h_ht_mc->Sumw2();  
  zllT->Project("h_ht_mc", "zll_ht","weight*(abs(Z_mass-91.19)<10 && nBJets<2)");
  //    TH1D* g_ht_mc = new TH1D("g_ht_mc","",20, 0, 2000); 
  g_ht_mc->Sumw2();
  gammaT->Project("g_ht_mc", "gamma_ht","weight*(prompt==2)");
  
  h_ht_mc->Divide(g_ht_mc);h_ht_mc->SetMarkerStyle(20); h_ht_mc->SetMarkerSize(1.6);
  h_ht_mc->SetMarkerColor(kRed);

  TH2D* h2_axes_ht = new TH2D("axes_ht", "", 10, 0, 2099, 10, 0., 0.5 );
  h2_axes_ht->SetYTitle("Zll / #gamma Ratio"); 
  h2_axes_ht->SetXTitle("HT [GeV]");
  h2_axes_ht->Draw();
 
  h_ht2->Draw("p same");
  // h_ht2->Draw("p same");
  // h_ht_mc->Draw("p same");
  labelTop->Draw("same");

  //  legend->Clear();
  // legend->AddEntry( h_ht2 ,"HT", "P" );
  //  legend->AddEntry( h_ht2 ,"looped over bins", "P" );
  //   legend->AddEntry( h_ht_mc, "HT MC", "P");
  //   legend->Draw("same");
    
  c1->SaveAs( Form("%s/ht_%s.png", outputdir.c_str() , RegIncl.getName().c_str()) );
  c1->SaveAs( Form("%s/ht_%s.eps", outputdir.c_str() , RegIncl.getName().c_str()) );
    

  TH2D* h2_axes_ht_mc = new TH2D("axes_ht_mc", "", 10, 0, 2099, 10, 0., 2 );
  h2_axes_ht_mc->SetXTitle("HT_MC [GeV]");
  h2_axes_ht_mc->SetYTitle("(Zll/ #gamma)_{Data} / (Zll/ #gamma)_{MC}");
  h2_axes_ht_mc->Draw();
     
  h_ht2->Divide(h_ht_mc);
  h_ht2->Draw("p same");
  labelTop->Draw("same");

  // legend->Clear();
  // legend->AddEntry( h_ht2 ,"HT", "P" );
  // legend->Draw("same");

  c1->SaveAs( Form("%s/ht_double.png", outputdir.c_str()) );
  c1->SaveAs( Form("%s/ht_double.eps", outputdir.c_str()) );
 



  //Filling the YIELDS
  int yieldBins_ht = h_ht->GetNbinsX();
  for(int bi = 1; bi <= yieldBins_ht ; bi++){
    double value = h_ht2->GetBinContent(bi);
    double err = h_ht2->GetBinError(bi);
    zll_ht->get(RegIncl)->yield->SetBinContent( bi,   value);
    zll_ht->get(RegIncl)->yield->SetBinError( bi,   err);
  }
  

  delete h2_axes_ht; 
  delete h_ht; delete g_ht; delete h_ht_mc; delete g_ht_mc;



  /////////////////////////////////////////////////////////////////////////////

    
  
  float bins_nJet40[] = {2,4,7,12};
  int binNum_nJet40 =  sizeof(bins_nJet40)/sizeof(float) -1;

  //Histos of NJET40, ZLL
  TH1D* h_nJet40 = new TH1D("h_nJet40","h_nJet40",binNum_nJet40, bins_nJet40);
  zllT->Project("h_nJet40", "nJets","weight*(( abs(Z_mass-91.19)<10 && nBJets<2))");
  //    h_nJet40->Sumw2();

  //Histos of NJET40, gamma
  TH1D* g_nJet40 = new TH1D("g_nJet40","g_nJet40",binNum_nJet40, bins_nJet40);
  gammaT->Project("g_nJet40","gamma_nJets","weight*(prompt==2)"); 
  //   g_nJet40->Sumw2();

  //MC ratio
  //Histos of NJET40, ZLL
  TH1D* h_nJet40_mc = new TH1D("h_nJet40_mc","h_nJet40_mc",binNum_nJet40, bins_nJet40);
  h_nJet40_mc->Sumw2();
  zllT->Project("h_nJet40_mc", "nJets","weight*(( abs(Z_mass-91.19)<10 && nBJets<2))");

  //Histos of NJET40, gamma
  TH1D* g_nJet40_mc = new TH1D("g_nJet40_mc","g_nJet40_mc",binNum_nJet40, bins_nJet40);
  g_nJet40_mc->Sumw2();
  gammaT->Project("g_nJet40_mc","gamma_nJets","weight*(prompt==2)"); 




  int nBins_nJet40 =  h_nJet40->GetNbinsX();
  //  h_nJet40->SetBinContent(nBins_nJet40, h_nJet40->GetBinContent(nBins_nJet40)+h_nJet40->GetBinContent(nBins_nJet40+1));
  //  g_nJet40->SetBinContent(nBins_nJet40, g_nJet40->GetBinContent(nBins_nJet40)+g_nJet40->GetBinContent(nBins_nJet40+1));
  for(int binnie = 1; binnie <= nBins_nJet40+1; binnie++){
    double value = h_nJet40->GetBinContent(binnie);
    h_nJet40->SetBinError(binnie,sqrt(value));
    double value_g = g_nJet40->GetBinContent(binnie);
    g_nJet40->SetBinError(binnie,sqrt(value_g));
  }


  h_nJet40->SetMarkerStyle(20);   h_nJet40->SetMarkerSize(1.6);

  TH2D* h2_axes_yield_nJet40 = new TH2D("axes_yield_nJet40", "", 10,0, 13, 10, 1,1000);
  h2_axes_yield_nJet40->SetXTitle("nJet");
  h2_axes_yield_nJet40->SetYTitle("Zll Yield");
  h2_axes_yield_nJet40->Draw();
    
  h_nJet40->Draw("p same");
  labelTop->Draw("same");
  gPad->RedrawAxis();
  // gPad->SetLogy();

  c1->SaveAs( Form("%s/yield_nJet40_%s.png", outputdir.c_str() , RegIncl.getName().c_str()) );
  c1->SaveAs( Form("%s/yield_nJet40_%s.eps", outputdir.c_str() , RegIncl.getName().c_str()) );
   

  h_nJet40->Divide(g_nJet40); 


  float err_data = h_nJet40->GetBinError(3);
  float content_data = h_nJet40->GetBinContent(3);
  std::cout << err_data << std::endl;
  std::cout << content_data << std::endl;

    
  TH1D* h_nJet40_450 = new TH1D("h_nJet40_450","h_nJet40_450",binNum_nJet40, bins_nJet40);
  TH1D* h_nJet40_575 = new TH1D("h_nJet40_575","h_nJet40_575",binNum_nJet40, bins_nJet40);
  TH1D* h_nJet40_1000 = new TH1D("h_nJet40_1000","h_nJet40_1000",binNum_nJet40, bins_nJet40);
  TH1D* h_nJet40_1500 = new TH1D("h_nJet40_1500","h_nJet40_1500",binNum_nJet40, bins_nJet40);
  zllT->Project("h_nJet40_450", "nJets","weight*(zll_ht>450&&zll_ht<575 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
  zllT->Project("h_nJet40_575", "nJets","weight*(zll_ht>575&&zll_ht<1000 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
  zllT->Project("h_nJet40_1000", "nJets","weight*(zll_ht>1000&&zll_ht<1500 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
  zllT->Project("h_nJet40_1500", "nJets","weight*(zll_ht>1500 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
  h_nJet40_450->Sumw2();h_nJet40_575->Sumw2();h_nJet40_1000->Sumw2();h_nJet40_1500->Sumw2();
   
  // 450 < HT < 575
  TH1D* g_nJet40_450 = new TH1D("g_nJet40_450","g_nJet40_450",binNum_nJet40, bins_nJet40);
  TH1D* g_nJet40_575 = new TH1D("g_nJet40_575","g_nJet40_575",binNum_nJet40, bins_nJet40);
  TH1D* g_nJet40_1000 = new TH1D("g_nJet40_1000","g_nJet40_1000",binNum_nJet40, bins_nJet40);
  TH1D* g_nJet40_1500 = new TH1D("g_nJet40_1500","g_nJet40_1500",binNum_nJet40, bins_nJet40);
  
  gammaT->Project("g_nJet40_450","gamma_nJets","weight*(prompt==2&&gamma_ht>450 && gamma_ht<575)"); 
  gammaT->Project("g_nJet40_575","gamma_nJets","weight*(prompt==2&&gamma_ht>575 && gamma_ht<1000)"); 
  gammaT->Project("g_nJet40_1000","gamma_nJets","weight*(prompt==2&&gamma_ht>1000 && gamma_ht<1500)"); 
  gammaT->Project("g_nJet40_1500","gamma_nJets","weight*(prompt==2&&gamma_ht>1500)"); 

  g_nJet40_450->Sumw2();g_nJet40_575->Sumw2();g_nJet40_1000->Sumw2();g_nJet40_1500->Sumw2();
 
  h_nJet40_450->Divide(g_nJet40_450);      h_nJet40_575->Divide(g_nJet40_575);
  h_nJet40_1000->Divide(g_nJet40_1000);      h_nJet40_1500->Divide(g_nJet40_1500);
  h_nJet40_450->SetMarkerStyle(20);h_nJet40_450->SetMarkerSize(1.6);
  h_nJet40_450->SetMarkerColor(46);    h_nJet40_450->SetLineColor(46); 
  h_nJet40_575->SetMarkerStyle(20);h_nJet40_575->SetMarkerSize(1.6);
  h_nJet40_575->SetMarkerColor(29);      h_nJet40_575->SetLineColor(29);  
  h_nJet40_1000->SetMarkerStyle(20);h_nJet40_1000->SetMarkerSize(1.6); 
  h_nJet40_1000->SetMarkerColor(38);     h_nJet40_1000->SetLineColor(38); 
  h_nJet40_1500->SetMarkerStyle(20);   h_nJet40_1500->SetMarkerSize(1.6);
  h_nJet40_1500->SetMarkerColor(42);  h_nJet40_1500->SetLineColor(42); 
    

    
  TH2D* h2_axes_jets = new TH2D("axes_jets", "", 10, 0, 13, 10, 0., 0.5 );
  h2_axes_jets->SetXTitle("nJets");
  h2_axes_jets->SetYTitle("Zll / #gamma Ratio");
  h2_axes_jets->Draw();
  TLegend* legend_jets = new TLegend( 0.2, 0.9-(5)*0.06, 0.5, 0.9 );
  legend_jets->SetTextSize(0.038);
  legend_jets->SetTextFont(42);
  legend_jets->SetFillColor(0);
  legend_jets->AddEntry( h_nJet40 ,"HT inclusive", "P" );
  legend_jets->AddEntry( h_nJet40_450 ,"450<HT<575", "P" );
  legend_jets->AddEntry( h_nJet40_575 ,"575<HT<1000", "P" );
  legend_jets->AddEntry( h_nJet40_1000 ,"1000<HT<1500", "P" );
  legend_jets->AddEntry( h_nJet40_1500 ,"HT>1500", "P" );
  h_nJet40_450->Draw("p same");
  h_nJet40_575->Draw("p same");
  h_nJet40_1000->Draw("p same");
  h_nJet40_1500->Draw("p same");
    
  h_nJet40->Draw("p same");
  labelTop->Draw("same");
  legend_jets->Draw("same");
  gPad->RedrawAxis();
   
  c1->SaveAs( Form("%s/nJet40_%s.eps", outputdir.c_str() , RegIncl.getName().c_str()) );
  c1->SaveAs( Form("%s/nJet40_%s.png", outputdir.c_str() , RegIncl.getName().c_str()) );


  h_nJet40_mc->Divide(g_nJet40_mc);

  float err_mc = h_nJet40_mc->GetBinError(3);
  float content_mc = h_nJet40_mc->GetBinContent(3);
  std::cout << err_mc << std::endl;
  std::cout << content_mc << std::endl;

  h_nJet40->Divide(h_nJet40_mc);


  std::cout << "ERROR IN LAST NJETS BIN = " << std::endl;
  std::cout << sqrt( (err_data*err_data )/(content_mc*content_mc)+ (err_mc*err_mc*content_data*content_data)/(content_mc*content_mc*content_mc*content_mc)  ) << std::endl;

  std::cout << h_nJet40->GetBinError(3);

  TH2D* h2_axes_jets_rat = new TH2D("axes_jets_rat", "", 10, 0, 13, 10, 0., 2 );
  h2_axes_jets_rat->SetXTitle("nJets");
  h2_axes_jets_rat->SetYTitle("(Zll / #gamma)_{Data}/(Zll / #gamma)_{MC}");
  h2_axes_jets_rat->Draw();
  h_nJet40->Draw("p same");
  labelTop->Draw("same");
 
  c1->SaveAs( Form("%s/ratio_nJet40_%s.eps", outputdir.c_str() , RegIncl.getName().c_str()) );
  c1->SaveAs( Form("%s/ratio_nJet40_%s.png", outputdir.c_str() , RegIncl.getName().c_str()) );


  //Filling the YIELDS
  int yieldBins_nJets = h_nJet40->GetNbinsX();
  for(int bi = 1; bi <= yieldBins_nJets+1 ; bi++){
    double value = h_nJet40->GetBinContent(bi);
    double err = h_nJet40->GetBinError(bi);
    zll_nJets->get(RegIncl)->yield->SetBinContent( bi,   value);
    zll_nJets->get(RegIncl)->yield->SetBinError( bi,   err);
  }
  
  
  delete h_nJet40; delete h_nJet40_mc;
  delete g_nJet40; delete g_nJet40_mc;
  delete h_nJet40_450; delete h_nJet40_575; delete h_nJet40_1000;    delete h_nJet40_1500;
  delete g_nJet40_450; delete g_nJet40_575; delete g_nJet40_1000;    delete g_nJet40_1500;
  delete h2_axes_jets; delete h2_axes_jets_rat;









  float bins_nBJet20[] = {0,1,2,3,6};
  int binNum_nBJet20 =  sizeof(bins_nBJet20)/sizeof(float) -1;

  TH1D* h_nBJet20 = new TH1D("h_nBJet20","h_nBJet20",binNum_nBJet20, bins_nBJet20);
  TH1D* g_nBJet20 = new TH1D("g_nBJet20","g_nBJet20",binNum_nBJet20, bins_nBJet20);
   
   

  //Histos of NBJET20, ZLL
  zllT->Project("h_nBJet20", "nBJets","weight*(( abs(Z_mass-91.19)<10 && nBJets<2))");
  //   h_nBJet20->Sumw2(); 
  gammaT->Project("g_nBJet20","gamma_nBJets","weight*(prompt==2)");
  //  g_nBJet20->Sumw2();

  int nBins_nBJet20 =  h_nBJet20->GetNbinsX();
  // h_nBJet20->SetBinContent(nBins_nBJet20, h_nBJet20->GetBinContent(nBins_nBJet20)+h_nBJet20->GetBinContent(nBins_nBJet20+1));
  // g_nBJet20->SetBinContent(nBins_nBJet20, g_nBJet20->GetBinContent(nBins_nBJet20)+g_nBJet20->GetBinContent(nBins_nBJet20+1));
  for(int binnie = 1; binnie <= nBins_nBJet20+1; binnie++){
    double value = h_nBJet20->GetBinContent(binnie);
    h_nBJet20->SetBinError(binnie,sqrt(value));
    double value_g = g_nBJet20->GetBinContent(binnie);
    g_nBJet20->SetBinError(binnie,sqrt(value_g));
  }

  h_nBJet20->SetMarkerStyle(20);h_nBJet20->SetMarkerSize(1.6);

  TH2D* h2_axes_yield_nBJet20 = new TH2D("axes_yield_nBJet20", "", 10,0, 6, 10, 1,1050);
  h2_axes_yield_nBJet20->SetXTitle("nBJet");
  h2_axes_yield_nBJet20->SetYTitle("Zll Yield");
  h2_axes_yield_nBJet20->Draw();
    
  h_nBJet20->Draw("p same");
  labelTop->Draw("same");
  gPad->RedrawAxis();
  //  gPad->SetLogy();

  c1->SaveAs( Form("%s/yield_nBJet20_%s.png", outputdir.c_str() , RegIncl.getName().c_str()) );
  c1->SaveAs( Form("%s/yield_nBJet20_%s.eps", outputdir.c_str() , RegIncl.getName().c_str()) );
   

  h_nBJet20->Divide(g_nBJet20);

  //MC ratio
  //Histos of NJET40, ZLL
  TH1D* h_nBJet20_mc = new TH1D("h_nBJet20_mc","h_nBJet20_mc",binNum_nBJet20, bins_nBJet20);
  zllT->Project("h_nBJet20_mc", "nBJets","weight*(( abs(Z_mass-91.19)<10 && nBJets<2))");
  h_nBJet20_mc->Sumw2();
  //Histos of NJET40, gamma
  TH1D* g_nBJet20_mc = new TH1D("g_nBJet20_mc","g_nBJet20_mc",binNum_nBJet20, bins_nBJet20);
  gammaT->Project("g_nBJet20_mc","gamma_nBJets","weight*(prompt==2)"); 
  g_nBJet20_mc->Sumw2();


  TH1D* h_nBJet20_450 = new TH1D("h_nBJet20_450","h_nBJet20_450",binNum_nBJet20, bins_nBJet20);
  TH1D* h_nBJet20_575 = new TH1D("h_nBJet20_575","h_nBJet20_575",binNum_nBJet20, bins_nBJet20);
  TH1D* h_nBJet20_1000 = new TH1D("h_nBJet20_1000","h_nBJet20_1000",binNum_nBJet20, bins_nBJet20);
  TH1D* h_nBJet20_1500 = new TH1D("h_nBJet20_1500","h_nBJet20_1500",binNum_nBJet20, bins_nBJet20);
   
  TH1D* g_nBJet20_450 = new TH1D("g_nBJet20_450","g_nBJet20_450",binNum_nBJet20, bins_nBJet20);
  TH1D* g_nBJet20_575 = new TH1D("g_nBJet20_575","g_nBJet20_575",binNum_nBJet20, bins_nBJet20);
  TH1D* g_nBJet20_1000 = new TH1D("g_nBJet20_1000","g_nBJet20_1000",binNum_nBJet20, bins_nBJet20);
  TH1D* g_nBJet20_1500 = new TH1D("g_nBJet20_1500","g_nBJet20_1500",binNum_nBJet20, bins_nBJet20);
    

    
  zllT->Project("h_nBJet20_450", "nBJets","weight*(zll_ht>450&&zll_ht<575 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
  zllT->Project("h_nBJet20_575", "nBJets","weight*(zll_ht>575&&zll_ht<1000 &&( abs(Z_mass-91.19)<10 && nBJets<2))");
  zllT->Project("h_nBJet20_1000", "nBJets","weight*(zll_ht>1000&&zll_ht<1500 &&( abs(Z_mass-91.19)<10 && nBJets<2)");
  zllT->Project("h_nBJet20_1500", "nBJets","weight*(zll_ht>1500 &&( abs(Z_mass-91.19)<10 && nBJets<2)");

  h_nBJet20_450->Sumw2();h_nBJet20_575->Sumw2(); h_nBJet20_1000->Sumw2();h_nBJet20_1500->Sumw2();
  
  gammaT->Project("g_nBJet20_450","gamma_nBJets","weight*(prompt==2&&gamma_ht>450 && gamma_ht<575)"); 
  gammaT->Project("g_nBJet20_575","gamma_nBJets","weight*(prompt==2&&gamma_ht>575 && gamma_ht<1000)"); 
  gammaT->Project("g_nBJet20_1000","gamma_nBJets","weight*(prompt==2&&gamma_ht>1000 && gamma_ht<1500)");
  gammaT->Project("g_nBJet20_1500","gamma_nBJets","weight*(prompt==2&&gamma_ht>1500)"); 

  g_nBJet20_450->Sumw2();g_nBJet20_575->Sumw2(); 
  g_nBJet20_1000->Sumw2();g_nBJet20_1500->Sumw2();
  
  h_nBJet20_450->Divide(g_nBJet20_450);
  h_nBJet20_575->Divide(g_nBJet20_575);
  h_nBJet20_1000->Divide(g_nBJet20_1000);
  h_nBJet20_1500->Divide(g_nBJet20_1500);

  h_nBJet20_450->SetMarkerStyle(20);   h_nBJet20_450->SetMarkerSize(1.6);
  h_nBJet20_450->SetMarkerColor(46);    h_nBJet20_450->SetLineColor(46); 
  h_nBJet20_575->SetMarkerStyle(20);   h_nBJet20_575->SetMarkerSize(1.6);
  h_nBJet20_575->SetMarkerColor(29);     h_nBJet20_575->SetLineColor(29);
  h_nBJet20_1000->SetMarkerStyle(20);   h_nBJet20_1000->SetMarkerSize(1.6);
  h_nBJet20_1000->SetMarkerColor(38);    h_nBJet20_1000->SetLineColor(38);
  h_nBJet20_1500->SetMarkerStyle(20);   h_nBJet20_1500->SetMarkerSize(1.6);
  h_nBJet20_1500->SetMarkerColor(42);     h_nBJet20_1500->SetLineColor(42);
    
    
    
    TH2D* h2_axes_bjets = new TH2D("axes_bjets", "", 10, 0, 6, 10, 0., 0.5 );
    h2_axes_bjets->SetXTitle("nBJets");
    h2_axes_bjets->SetYTitle("Zll / #gamma Ratio");
    h2_axes_bjets->Draw();


    TLegend* legend_bjets = new TLegend( 0.2, 0.9-(5)*0.06, 0.5, 0.9 );
    legend_bjets->SetTextSize(0.038);
    legend_bjets->SetTextFont(42);
    legend_bjets->SetFillColor(0);
    legend_bjets->AddEntry( h_nBJet20 ,"HT inclusive", "P" );
    
    legend_bjets->AddEntry( h_nBJet20_450 ,"450<HT<575", "P" );
    legend_bjets->AddEntry( h_nBJet20_575 ,"575<HT<1000", "P" );
    legend_bjets->AddEntry( h_nBJet20_1000 ,"1000<HT<1500", "P" );
    legend_bjets->AddEntry( h_nBJet20_1500 ,"HT>1500", "P" );
    h_nBJet20_450->Draw("p same");
    h_nBJet20_575->Draw("p same");
    h_nBJet20_1000->Draw("p same");
    h_nBJet20_1500->Draw("p same");
        
    h_nBJet20->Draw("p same");

    labelTop->Draw("same");
    legend_bjets->Draw("same");

    gPad->RedrawAxis();
   
    c1->SaveAs( Form("%s/nBJet20_%s.eps", outputdir.c_str() , RegIncl.getName().c_str()) );
    c1->SaveAs( Form("%s/nBJet20_%s.png", outputdir.c_str() , RegIncl.getName().c_str()) );


    h_nBJet20_mc->Divide(g_nBJet20_mc);

    h_nBJet20->Divide(h_nBJet20_mc);
    TH2D* h2_axes_bjets_rat = new TH2D("axes_bjets_rat", "", 10, 0, 6, 10, 0., 2 );
    h2_axes_bjets_rat->SetXTitle("nBJets");
    h2_axes_bjets_rat->SetYTitle("(Zll / #gamma)_{Data}/(Zll / #gamma)_{MC}");
    h2_axes_bjets_rat->Draw();

    h_nBJet20->Draw("p same");
    labelTop->Draw("same");

    c1->SaveAs( Form("%s/ratio_nBJet20_%s.eps", outputdir.c_str() , RegIncl.getName().c_str()) );
    c1->SaveAs( Form("%s/ratio_nBJet20_%s.png", outputdir.c_str() , RegIncl.getName().c_str()) );


    //Filling the YIELDS
    int yieldBins_nBJets = h_nBJet20->GetNbinsX();
    for(int bi = 1; bi <= yieldBins_nBJets+1 ; bi++){
      double value = h_nBJet20->GetBinContent(bi);
      double err = h_nBJet20->GetBinError(bi);
      zll_nBJets->get(RegIncl)->yield->SetBinContent( bi,   value);
      zll_nBJets->get(RegIncl)->yield->SetBinError( bi,   err);
    }
  



    delete h_nBJet20;    delete h_nBJet20_mc;
    delete g_nBJet20;    delete g_nBJet20_mc;
    delete h_nBJet20_450; delete h_nBJet20_575; delete h_nBJet20_1000; delete h_nBJet20_1500;
    delete g_nBJet20_450; delete g_nBJet20_575; delete g_nBJet20_1000; delete g_nBJet20_1500;
    delete h2_axes_bjets; delete h2_axes_bjets_rat;
    

 
  */
  }
  
    //FOR THE Z MASS PLOT
    TH1D* h_Zmass = new TH1D("h_Zmass","h_Zmass", 60, 60, 120);  //  h_Zmass->Sumw2();
    zllT->Project("h_Zmass", "Z_mass","weight");
    //  h_Zmass->SetMarkerStyle(20);
    //  h_Zmass->SetMarkerSize(1.6);
    h_Zmass->SetLineWidth(2);
    h_Zmass->SetLineColor(kBlue+1);

    TH2D* h2_axes2 = new TH2D("axes2", "", 10, 60, 120 , 10, 0., 170 );
    h2_axes2->SetXTitle("M_{ll} [GeV]");
    h2_axes2->SetYTitle("Entries");
    h2_axes2->Draw();
    h_Zmass->Draw("same");    labelTop->Draw("same");
 
    c1->SaveAs( Form("%s/Zmass.png", outputdir.c_str()) );
    c1->SaveAs( Form("%s/Zmass.eps", outputdir.c_str()) );
   


  std::cout << counter_total_bins << std::endl;
  std::cout << counter_empty_bins << std::endl;




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
  zllY_mt2->writeToFile(outFile);

  zllG_ht->addToFile( outFile );
  zllG_nJets->addToFile( outFile );
  zllG_nBJets->addToFile( outFile );
  zllG_mt2->addToFile( outFile );

  /*
  zll_ht->writeToFile( outFile_ht );
  zll_nJets->writeToFile( outFile_nJets );
  zll_nBJets->writeToFile( outFile_nBJets );
  */

  return 0;
}


















void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string zll_sel,  std::string gamma_sel, MT2Analysis<MT2Estimate>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma, MT2Analysis<MT2EstimateTree>*  Zll, MT2Analysis<MT2Estimate>*  zll_yield, const MT2Region thisRegion, std::string cut, std::string cut_gamma){
 
  std::vector<int> colors;
  colors.push_back(430); // other = zll 
  colors.push_back(401); // qcd
  colors.push_back(417); // w+jets
  colors.push_back(419); // z+jets
  colors.push_back(855); // top

  TH1F::AddDirectory(kTRUE);

  float bins[size+1]; for(unsigned int i=0; i<= size ; i++)      bins[i]=binss[i];
  float xMin = binss[0];
  float xMax = binss[size];


  //THE TREES
  TTree *zllT =  Zll->get(thisRegion)->tree;
  TTree *gammaT =  gamma->get(thisRegion)->tree;

  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
 


  
  TCanvas* canny = new TCanvas( "canny", "", 600, 700 );
  canny->cd();

  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();


  TH1D* h_mt2 = new TH1D("h_mt2","", size  , bins);
  TH1D* g_mt2 = new TH1D("g_mt2","", size  , bins);
   
  zllT ->Project( "h_mt2" , zll_sel.c_str(), cut.c_str() );
  gammaT->Project( "g_mt2", gamma_sel.c_str(), cut_gamma.c_str() );

  h_mt2->SetBinContent(size, h_mt2->GetBinContent(size) + h_mt2->GetBinContent(size+1));//adding overflow
  g_mt2->SetBinContent(size, g_mt2->GetBinContent(size) + g_mt2->GetBinContent(size+1));

  int nBinss =  h_mt2->GetNbinsX();
  for(int binnie = 1; binnie <= nBinss; binnie++){
    double value = h_mt2->GetBinContent(binnie);
    h_mt2->SetBinError(binnie,sqrt(value));
    double value_g = g_mt2->GetBinContent(binnie);
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

    zllT ->Project( "h_mt2_mc" , zll_sel.c_str(), cut.c_str() );
    gammaT->Project( "g_mt2_mc", gamma_sel.c_str(), cut_gamma.c_str() );
  
    h_mt2_mc->SetBinContent(size, h_mt2_mc->GetBinContent(size) +h_mt2_mc->GetBinContent(size+1));
    g_mt2_mc->SetBinContent(size, g_mt2_mc->GetBinContent(size) +g_mt2_mc->GetBinContent(size+1));



    h_mt2->Divide(g_mt2);  
    h_mt2_mc->Divide(g_mt2_mc);
    //   h_mt2_mc->SetMarkerStyle(20); h_mt2_mc->SetMarkerSize(1.6);
    //   h_mt2_mc->SetMarkerColor(46); h_mt2_mc->SetLineColor(46);
   h_mt2_mc->SetLineColor(kBlue+1); h_mt2_mc->SetLineWidth(2);

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 0.5 );
    //    TH2D* h2_axes = new TH2D("axes", "", 10, 0, 1500, 10, 0., 0.5 );
    if(zll_sel == "zll_mt2"){
      h2_axes->SetXTitle("M_{T2} [GeV]");
    }else    if(zll_sel == "zll_ht"){
      h2_axes->SetXTitle("H_{T} [GeV]");
    }else    if(zll_sel == "nJets"){
      h2_axes->SetXTitle("Jet Multiplicity");
    }else{
      h2_axes->SetXTitle("b Jet Multiplicity" );
    }

    h2_axes->SetYTitle("Zll / #gamma Ratio");
    h2_axes->Draw();

   h_mt2_mc->Draw("L same");
     h_mt2->DrawClone("p same");
     
    labelTop->Draw("same");
    gPad->RedrawAxis();

    TLegend* legend = new TLegend( 0.2, 0.9-(2)*0.06, 0.5, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( h_mt2 ,"Data", "P" );
    legend->AddEntry( h_mt2_mc ,"Simulation", "L" );
    legend->Draw("same");



 
    gPad->RedrawAxis();

    canny->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.1);
    pad2->Draw();
    pad2->cd();
      
    h_mt2->Divide(h_mt2_mc); //the ratio is now in h_mt2


    TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , 0.,2  );
    // TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, 0, 1500, 5 , 0.,2  );
    // h2_axes_rat->SetXTitle("M_{T2} [GeV]");
    h2_axes_rat->SetYTitle("Data / MC");
    // h2_axes_rat->SetYTitle("(Zll/ #gamma)_{Data} / (Zll/ #gamma)_{MC}");



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

    h_mt2->Draw("p same");



    canny->cd();

    canny->SaveAs( Form("%s/%s_ratios_%s.eps", fullPath.c_str(), zll_sel.c_str(),  thisRegion.getName().c_str() ) );
    canny->SaveAs( Form("%s/%s_ratios_%s.png", fullPath.c_str(), zll_sel.c_str(),thisRegion.getName().c_str() ) );




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




    delete canny;







}











