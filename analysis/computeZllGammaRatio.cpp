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



#include "../interface/MT2DrawTools.h"


#include <iostream>

float lumi = 4.;

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

  std::string outputdir = "ZllGamma_Ratio_" + regionsSet;
  //  std::string outputdir = "Zll_" + configFileName;
  double intpart;
  double fracpart = modf(lumi, &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("_%.0ffb", intpart ) );
  // outputdir += suffix;
  
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
 
  std::string gammaControlRegionDir = "GammaControlRegion_" + samples + "_" + regionsSet + "_4fb";

 
  MT2Analysis<MT2EstimateTree>* gamma = MT2Analysis<MT2EstimateTree>::readFromFile(gammaControlRegionDir + "/data.root", "gammaCRtree");
  if( gamma==0 ) {
    std::cout << "-> Please run gammaControlRegion first. I need to get the gammaCR yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(193);
  }


  std::string ZllDir = "Zll_" + samples + "_" + regionsSet;

  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb/Zll_analyses.root", ZllDir.c_str(), lumi), "DYJets");
  if( Zll==0 ) {
    std::cout << "-> Please run computeZinvFromZll first. I need to get the Z->vv MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }

  


  std::set<MT2Region> MT2Regions = Zll->getRegions();
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
    // std::set<MT2Region>::iterator iMT2 = MT2Regions.begin();  
    std::vector<std::string> niceNames = thisRegion.getNiceNames();


    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    TH1F::AddDirectory(kTRUE);

    //THE TREES
    TTree *zllT =  Zll->get(*iMT2)->tree;
    TTree *gammaT =  gamma->get(*iMT2)->tree;



    //Histos of MT2, ZLL
    TH1D* h_mt2 = new TH1D("h_mt2","h_mt2",10, 0, 1500);
    zllT->Project("h_mt2", "zll_mt2","weight*(( abs(Z_mass-91.19)<20))");
    /*
    TH1D* h_mt2_450 = new TH1D("h_mt2_450","h_mt2_450",10, 0, 1500); 
    zllT->Project("h_mt2_450", "zll_mt2","weight*(zll_ht>450&&zll_ht<575 &&( abs(Z_mass-91.19)<20))");
    TH1D* h_mt2_575 = new TH1D("h_mt2_575","h_mt2_575",10, 0, 1500); 
    zllT->Project("h_mt2_575", "zll_mt2","weight*(zll_ht>575&&zll_ht<1000 &&( abs(Z_mass-91.19)<20))");
    //1000 < HT < 1500
    TH1D* h_mt2_1000 = new TH1D("h_mt2_1000","h_mt2_1000",10, 0, 1500); 
    zllT->Project("h_mt2_1000", "zll_mt2","weight*(zll_ht>1000&&zll_ht<1500 &&( abs(Z_mass-91.19)<20))");
    //       HT > 1500
    TH1D* h_mt2_1500 = new TH1D("h_mt2_1500","h_mt2_1500",10, 0, 1500); 
    zllT->Project("h_mt2_1500", "zll_mt2","weight*(zll_ht>1500 &&( abs(Z_mass-91.19)<20))");
    
    h_mt2_450->Sumw2();h_mt2_575->Sumw2();h_mt2_1000->Sumw2(); h_mt2_1500->Sumw2();
    */

    h_mt2->Sumw2();

    //Histos of MT2, gamma
    //Inclusive
    TH1D* g_mt2 = new TH1D("g_mt2","",10, 0, 1500); 
    gammaT->Project("g_mt2","gamma_mt2","weight*(prompt==2)"); 
  
    /*  // 450 < HT < 575
    TH1D* g_mt2_450 = new TH1D("g_mt2_450","",10, 0, 1500); 
    gammaT->Project("g_mt2_450","gamma_mt2","weight*(prompt==2&&gamma_ht>450 && gamma_ht<575)"); 
    // 575 < HT < 1000
    TH1D* g_mt2_575 = new TH1D("g_mt2_575","",10, 0, 1500); 
    gammaT->Project("g_mt2_575","gamma_mt2","weight*(prompt==2&&gamma_ht>575 && gamma_ht<1000)"); 
    // 1000 < HT < 1500
    TH1D* g_mt2_1000 = new TH1D("g_mt2_1000","",10, 0, 1500); 
    gammaT->Project("g_mt2_1000","gamma_mt2","weight*(prompt==2&&gamma_ht>1000 && gamma_ht<1500)"); 
    // HT > 1500
    TH1D* g_mt2_1500 = new TH1D("g_mt2_1500","",10, 0, 1500); 
    gammaT->Project("g_mt2_1500","gamma_mt2","weight*(prompt==2&&gamma_ht>1500)"); 

    g_mt2_450->Sumw2();g_mt2_575->Sumw2();g_mt2_1000->Sumw2();g_mt2_1500->Sumw2();
    */

    g_mt2->Sumw2();

    h_mt2->Divide(g_mt2);    h_mt2->SetMarkerStyle(20);   h_mt2->SetMarkerSize(1.6);
    /*
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
    */

    TH2D* h2_axes = new TH2D("axes", "", 10, 0, 1500, 10, 0., 1 );
    h2_axes->SetXTitle("M_{T2} [GeV]");
    h2_axes->SetYTitle("Zll / #gamma Ratio");
    h2_axes->Draw();


    TLegend* legend = new TLegend( 0.2, 0.9-(5)*0.06, 0.5, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( h_mt2 ,"HT inclusive", "P" );
    /* legend->AddEntry( h_mt2_450 ,"450<HT<575", "P" );
    legend->AddEntry( h_mt2_575 ,"575<HT<1000", "P" );
    legend->AddEntry( h_mt2_1000 ,"1000<HT<1500", "P" );
    legend->AddEntry( h_mt2_1500 ,"HT>1500", "P" );
 
    h_mt2_450->Draw("p same");
    h_mt2_575->Draw("p same");
    h_mt2_1000->Draw("p same");
    h_mt2_1500->Draw("p same");
    */
    h_mt2->Draw("p same");

    TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
    labelTop->Draw("same");
    // legend->Draw("same");

    gPad->RedrawAxis();
   
    c1->SaveAs( Form("%s/mt2_%s.eps", outputdir.c_str(), thisRegion.getName().c_str() ) );
    c1->SaveAs( Form("%s/mt2_%s.png", outputdir.c_str() , thisRegion.getName().c_str() ) );


    //Double ratio
    TH1D* h_mt2 = new TH1D("h_mt2","h_mt2",10, 0, 1500); h_mt2->Sumw2();
    zllT->Project("h_mt2", "zll_mt2","weight*(( abs(Z_mass-91.19)<20))");
 
    //Histos of MT2, gamma
    //Inclusive
    TH1D* g_mt2 = new TH1D("g_mt2","",10, 0, 1500); 
    gammaT->Project("g_mt2","gamma_mt2","weight*(prompt==2)"); 
  


    delete h_mt2;
    delete g_mt2;
    delete h2_axes;

    //   delete h_mt2_450;    delete h_mt2_575;    delete h_mt2_1000;    delete h_mt2_1500;
    //   delete g_mt2_450;    delete g_mt2_575;    delete g_mt2_1000;    delete g_mt2_1500;



    //FOR THE HT plot
    /*
    TH1D* h_ht2 = new TH1D("h_ht2","h_ht2",20, 0, 2000); 
    zllT->Project("h_ht2", "zll_ht","weight*(abs(Z_mass-91.19)<20)");
    TH1D* g_ht2 = new TH1D("g_ht2","",20, 0, 2000);  
    gammaT->Project("g_ht2", "gamma_ht","weight*(prompt==2)");
    
    int nBins =  h_ht2->GetNbinsX();
    for(int binnie = 1; binnie < nBins; binnie++){
      double value = h_ht2->GetBinContent(binnie);
      h_ht2->SetBinError(binnie,sqrt(value));
      double value_g = g_ht2->GetBinContent(binnie);
      g_ht2->SetBinError(binnie,sqrt(value_g));
    }
    h_ht2->Divide(g_ht2); 
    h_ht2->SetMarkerStyle(20);    h_ht2->SetMarkerSize(1.6);
    h_ht2->SetMarkerColor(kRed);  h_ht2->SetLineColor(kRed);
    */
    TH1D* h_ht = new TH1D("h_ht","h_ht",20, 0, 2000); 
    zllT->Project("h_ht", "zll_ht","weight*(abs(Z_mass-91.19)<20)");
    TH1D* g_ht = new TH1D("g_ht","",20, 0, 2000);  
    gammaT->Project("g_ht", "gamma_ht","weight*(prompt==2)");
 
    h_ht->Sumw2(); g_ht->Sumw2();
    h_ht->Divide(g_ht);    h_ht->SetMarkerStyle(20);    h_ht->SetMarkerSize(1.6);
 

    //Double ratioooooo
    TH1D* h_ht_mc = new TH1D("h_ht_mc","h_ht_mc",20, 0, 2000);   h_ht_mc->Sumw2();  
    zllT->Project("h_ht_mc", "zll_ht","weight*(abs(Z_mass-91.19)<20)");
    TH1D* g_ht_mc = new TH1D("g_ht_mc","",20, 0, 2000);  g_ht_mc->Sumw2();
    gammaT->Project("g_ht_mc", "gamma_ht","weight*(prompt==2)");
  
    h_ht_mc->Divide(g_ht_mc);h_ht_mc->SetMarkerStyle(20); h_ht_mc->SetMarkerSize(1.6);
    h_ht_mc->SetMarkerColor(kRed);

    TH2D* h2_axes_ht = new TH2D("axes_ht", "", 10, 0, 2099, 10, 0., 1 );
    h2_axes_ht->SetYTitle("Zll / #gamma Ratio"); 
    h2_axes_ht->SetXTitle("HT [GeV]");
    h2_axes_ht->Draw();
 
    h_ht->Draw("p same");
    // h_ht2->Draw("p same");
    // h_ht_mc->Draw("p same");
    labelTop->Draw("same");

    legend->Clear();
    legend->AddEntry( h_ht ,"HT", "P" );
    //  legend->AddEntry( h_ht2 ,"looped over bins", "P" );
    //   legend->AddEntry( h_ht_mc, "HT MC", "P");
    //   legend->Draw("same");
    
    c1->SaveAs( Form("%s/ht_%s.png", outputdir.c_str() , thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/ht_%s.eps", outputdir.c_str() , thisRegion.getName().c_str()) );
  

    /*

    h_ht->Divide(h_ht_mc);
    h2_axes_ht->SetYTitle("(Zll/ #gamma)_{Data} / (Zll/ #gamma)_{MC}");
    h2_axes_ht->Draw();
    h_ht->Draw("p same");
    labelTop->Draw("same");

    legend->Clear();
    legend->AddEntry( h_ht ,"HT", "P" );
    //   legend->Draw("same");

    c1->SaveAs( Form("%s/ht_double.png", outputdir.c_str()) );
    c1->SaveAs( Form("%s/ht_double.eps", outputdir.c_str()) );
    */

    delete h2_axes_ht; 
    delete h_ht; delete g_ht; delete h_ht_mc; delete g_ht_mc;




   //Histos of NJET40, ZLL
    TH1D* h_nJet40 = new TH1D("h_nJet40","h_nJet40",15, 0, 15); 
    zllT->Project("h_nJet40", "nJets","weight*(( abs(Z_mass-91.19)<20))");
    /*
    TH1D* h_nJet40_450 = new TH1D("h_nJet40_450","h_nJet40_450",15, 0, 15); 
    zllT->Project("h_nJet40_450", "nJets","weight*(zll_ht>450&&zll_ht<575 &&( abs(Z_mass-91.19)<20))");
    TH1D* h_nJet40_575 = new TH1D("h_nJet40_575","h_nJet40_575",15, 0, 15); 
    zllT->Project("h_nJet40_575", "nJets","weight*(zll_ht>575&&zll_ht<1000 &&( abs(Z_mass-91.19)<20))");
    //1000 < HT < 1500
    TH1D* h_nJet40_1000 = new TH1D("h_nJet40_1000","h_nJet40_1000",15, 0, 15); 
    zllT->Project("h_nJet40_1000", "nJets","weight*(zll_ht>1000&&zll_ht<1500 &&( abs(Z_mass-91.19)<20))");
    //       HT > 1500
    TH1D* h_nJet40_1500 = new TH1D("h_nJet40_1500","h_nJet40_1500",15, 0, 15); 
    zllT->Project("h_nJet40_1500", "nJets","weight*(zll_ht>1500 &&( abs(Z_mass-91.19)<20))");

    h_nJet40_450->Sumw2();h_nJet40_575->Sumw2();h_nJet40_1000->Sumw2();h_nJet40_1500->Sumw2();
    */
    
    h_nJet40->Sumw2();

    //Histos of NJET40, gamma
    //Inclusive
    TH1D* g_nJet40 = new TH1D("g_nJet40","",15, 0, 15);
    gammaT->Project("g_nJet40","gamma_nJets","weight*(prompt==2)"); 
    /*
    // 450 < HT < 575
    TH1D* g_nJet40_450 = new TH1D("g_nJet40_450","",15, 0, 15); 
    gammaT->Project("g_nJet40_450","gamma_nJets","weight*(prompt==2&&gamma_ht>450 && gamma_ht<575)"); 
    // 575 < HT < 1000
    TH1D* g_nJet40_575 = new TH1D("g_nJet40_575","",15, 0, 15);
    gammaT->Project("g_nJet40_575","gamma_nJets","weight*(prompt==2&&gamma_ht>575 && gamma_ht<1000)"); 
    // 1000 < HT < 1500
    TH1D* g_nJet40_1000 = new TH1D("g_nJet40_1000","",15, 0, 15); 
    gammaT->Project("g_nJet40_1000","gamma_nJets","weight*(prompt==2&&gamma_ht>1000 && gamma_ht<1500)"); 
    // HT > 1500
    TH1D* g_nJet40_1500 = new TH1D("g_nJet40_1500","",15, 0, 15); 
    gammaT->Project("g_nJet40_1500","gamma_nJets","weight*(prompt==2&&gamma_ht>1500)"); 

    g_nJet40_450->Sumw2();g_nJet40_575->Sumw2();g_nJet40_1000->Sumw2();g_nJet40_1500->Sumw2();
    */
    g_nJet40->Sumw2();

    h_nJet40->Divide(g_nJet40); h_nJet40->SetMarkerStyle(20);   h_nJet40->SetMarkerSize(1.6);
    /*
    h_nJet40_450->Divide(g_nJet40_450);
    h_nJet40_575->Divide(g_nJet40_575);
    h_nJet40_1000->Divide(g_nJet40_1000);
    h_nJet40_1500->Divide(g_nJet40_1500);

  h_nJet40_450->SetMarkerStyle(20);   h_nJet40_450->SetMarkerSize(1.6); h_nJet40_450->SetMarkerColor(46);
    h_nJet40_450->SetLineColor(46); 
    h_nJet40_575->SetMarkerStyle(20);   h_nJet40_575->SetMarkerSize(1.6); h_nJet40_575->SetMarkerColor(29);  
    h_nJet40_575->SetLineColor(29);  
    h_nJet40_1000->SetMarkerStyle(20);   h_nJet40_1000->SetMarkerSize(1.6); h_nJet40_1000->SetMarkerColor(38); 
    h_nJet40_1000->SetLineColor(38); 
    h_nJet40_1500->SetMarkerStyle(20);   h_nJet40_1500->SetMarkerSize(1.6); h_nJet40_1500->SetMarkerColor(42); 
    h_nJet40_1500->SetLineColor(42); 
    */


    TH2D* h2_axes_jets = new TH2D("axes_jets", "", 10, 0, 13, 10, 0., 1 );
    h2_axes_jets->SetXTitle("nJets");
    h2_axes_jets->SetYTitle("Zll / #gamma Ratio");
    h2_axes_jets->Draw();


    TLegend* legend_jets = new TLegend( 0.45, 0.9-(5)*0.06, 0.93, 0.9 );
    legend_jets->SetTextSize(0.038);
    legend_jets->SetTextFont(42);
    legend_jets->SetFillColor(0);
    legend_jets->AddEntry( h_nJet40 ,"HT inclusive", "P" );
    /*    legend_jets->AddEntry( h_nJet40_450 ,"450<HT<575", "P" );
    legend_jets->AddEntry( h_nJet40_575 ,"575<HT<1000", "P" );
    legend_jets->AddEntry( h_nJet40_1000 ,"1000<HT<1500", "P" );
    legend_jets->AddEntry( h_nJet40_1500 ,"HT>1500", "P" );
 
    h_nJet40_450->Draw("p same");
    h_nJet40_575->Draw("p same");
    h_nJet40_1000->Draw("p same");
    h_nJet40_1500->Draw("p same");
    */
    h_nJet40->Draw("p same");

    labelTop->Draw("same");
    //   legend_jets->Draw("same");

    gPad->RedrawAxis();
   
    c1->SaveAs( Form("%s/nJet40_%s.eps", outputdir.c_str() , thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/nJet40_%s.png", outputdir.c_str() , thisRegion.getName().c_str()) );



  
    delete h_nJet40;
    delete g_nJet40;
    //    delete h_nJet40_450; delete h_nJet40_575; delete h_nJet40_1000;    delete h_nJet40_1500;
    //    delete g_nJet40_450; delete g_nJet40_575; delete g_nJet40_1000;    delete g_nJet40_1500;
    delete h2_axes_jets;











   //Histos of NJET40, ZLL
    TH1D* h_nBJet40 = new TH1D("h_nBJet40","h_nBJet40",6, 0, 6); 
    zllT->Project("h_nBJet40", "nBJets","weight*(( abs(Z_mass-91.19)<20))");
    /*
    TH1D* h_nBJet40_450 = new TH1D("h_nBJet40_450","h_nBJet40_450",6, 0, 6); 
    zllT->Project("h_nBJet40_450", "nBJets","weight*(zll_ht>450&&zll_ht<575 &&( abs(Z_mass-91.19)<20))");
    TH1D* h_nBJet40_575 = new TH1D("h_nBJet40_575","h_nBJet40_575",6, 0, 6); 
    zllT->Project("h_nBJet40_575", "nBJets","weight*(zll_ht>575&&zll_ht<1000 &&( abs(Z_mass-91.19)<20))");
    //1000 < HT < 1500
    TH1D* h_nBJet40_1000 = new TH1D("h_nBJet40_1000","h_nBJet40_1000",6, 0, 6); 
    zllT->Project("h_nBJet40_1000", "nBJets","weight*(zll_ht>1000&&zll_ht<1500 &&( abs(Z_mass-91.19)<20))");
    //       HT > 1500
    TH1D* h_nBJet40_1500 = new TH1D("h_nBJet40_1500","h_nBJet40_1500",6, 0, 6); 
    zllT->Project("h_nBJet40_1500", "nBJets","weight*(zll_ht>1500 &&( abs(Z_mass-91.19)<20))");

h_nBJet40_450->Sumw2();h_nBJet40_575->Sumw2(); h_nBJet40_1000->Sumw2();h_nBJet40_1500->Sumw2();
    */

    h_nBJet40->Sumw2(); 
    //Histos of NJET40, gamma
    //Inclusive
    TH1D* g_nBJet40 = new TH1D("g_nBJet40","",6, 0, 6); 
    gammaT->Project("g_nBJet40","gamma_nBJets","weight*(prompt==2)"); 
    /*
    // 450 < HT < 575
    TH1D* g_nBJet40_450 = new TH1D("g_nBJet40_450","",6, 0, 6); 
    gammaT->Project("g_nBJet40_450","gamma_nBJets","weight*(prompt==2&&gamma_ht>450 && gamma_ht<575)"); 
    // 575 < HT < 1000
    TH1D* g_nBJet40_575 = new TH1D("g_nBJet40_575","",6, 0, 6); 
    gammaT->Project("g_nBJet40_575","gamma_nBJets","weight*(prompt==2&&gamma_ht>575 && gamma_ht<1000)"); 
    // 1000 < HT < 1500
    TH1D* g_nBJet40_1000 = new TH1D("g_nBJet40_1000","",6, 0, 6);
    gammaT->Project("g_nBJet40_1000","gamma_nBJets","weight*(prompt==2&&gamma_ht>1000 && gamma_ht<1500)"); 
    // HT > 1500
    TH1D* g_nBJet40_1500 = new TH1D("g_nBJet40_1500","",6, 0, 6); 
    gammaT->Project("g_nBJet40_1500","gamma_nBJets","weight*(prompt==2&&gamma_ht>1500)"); 

g_nBJet40_450->Sumw2();g_nBJet40_575->Sumw2(); 
    g_nBJet40_1000->Sumw2();g_nBJet40_1500->Sumw2();
    */
    g_nBJet40->Sumw2();

    h_nBJet40->Divide(g_nBJet40); h_nBJet40->SetMarkerStyle(20);h_nBJet40->SetMarkerSize(1.6);
    /*
    h_nBJet40_450->Divide(g_nBJet40_450);
    h_nBJet40_575->Divide(g_nBJet40_575);
    h_nBJet40_1000->Divide(g_nBJet40_1000);
    h_nBJet40_1500->Divide(g_nBJet40_1500);

 h_nBJet40_450->SetMarkerStyle(20);   h_nBJet40_450->SetMarkerSize(1.6); h_nBJet40_450->SetMarkerColor(46);
    h_nBJet40_450->SetLineColor(46); 
    h_nBJet40_575->SetMarkerStyle(20);   h_nBJet40_575->SetMarkerSize(1.6); h_nBJet40_575->SetMarkerColor(29); 
    h_nBJet40_575->SetLineColor(29);
    h_nBJet40_1000->SetMarkerStyle(20);   h_nBJet40_1000->SetMarkerSize(1.6); h_nBJet40_1000->SetMarkerColor(38);
    h_nBJet40_1000->SetLineColor(38);
    h_nBJet40_1500->SetMarkerStyle(20);   h_nBJet40_1500->SetMarkerSize(1.6); h_nBJet40_1500->SetMarkerColor(42); 
    h_nBJet40_1500->SetLineColor(42);
    */
    TH2D* h2_axes_bjets = new TH2D("axes_bjets", "", 10, 0, 6, 10, 0., 1 );
    h2_axes_bjets->SetXTitle("nBJets");
    h2_axes_bjets->SetYTitle("Zll / #gamma Ratio");
    h2_axes_bjets->Draw();


    TLegend* legend_bjets = new TLegend( 0.45, 0.9-(5)*0.06, 0.93, 0.9 );
    legend_bjets->SetTextSize(0.038);
    legend_bjets->SetTextFont(42);
    legend_bjets->SetFillColor(0);
    legend_bjets->AddEntry( h_nBJet40 ,"HT inclusive", "P" );
    /* 
       legend_bjets->AddEntry( h_nBJet40_450 ,"450<HT<575", "P" );
    legend_bjets->AddEntry( h_nBJet40_575 ,"575<HT<1000", "P" );
    legend_bjets->AddEntry( h_nBJet40_1000 ,"1000<HT<1500", "P" );
    legend_bjets->AddEntry( h_nBJet40_1500 ,"HT>1500", "P" );
 
    h_nBJet40_450->Draw("p same");
    h_nBJet40_575->Draw("p same");
    h_nBJet40_1000->Draw("p same");
    h_nBJet40_1500->Draw("p same");
    */ 
    h_nBJet40->Draw("p same");

    labelTop->Draw("same");
    //    legend_bjets->Draw("same");

    gPad->RedrawAxis();
   
    c1->SaveAs( Form("%s/nBJet40_%s.eps", outputdir.c_str() , thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/nBJet40_%s.png", outputdir.c_str() , thisRegion.getName().c_str()) );


    delete h_nBJet40;
    delete g_nBJet40;
    //   delete h_nBJet40_450; delete h_nBJet40_575; delete h_nBJet40_1000; delete h_nBJet40_1500;
    //    delete g_nBJet40_450; delete g_nBJet40_575; delete g_nBJet40_1000; delete g_nBJet40_1500;
    delete h2_axes_bjets;


  }
  
 
    /*
    //FOR THE Z MASS PLOT
    TH1D* h_Zmass = new TH1D("h_Zmass","h_Zmass",200, 0, 200);  //  h_Zmass->Sumw2();
    zllT->Project("h_Zmass", "Z_mass","");
    //  h_Zmass->SetMarkerStyle(20);
    //  h_Zmass->SetMarkerSize(1.6);
    h_Zmass->SetLineWidth(2);

    TH2D* h2_axes2 = new TH2D("axes2", "", 10, 0, 200, 10, 0., 95000 );
    h2_axes2->SetXTitle("M_{ll} [GeV]");
    h2_axes2->SetYTitle("Entries");
    h2_axes2->Draw();
    h_Zmass->Draw("same");    labelTop->Draw("same");
 
    c1->SaveAs( Form("%s/Zmass.png", outputdir.c_str()) );
    c1->SaveAs( Form("%s/Zmass.eps", outputdir.c_str()) );
    */

  return 0;
}
