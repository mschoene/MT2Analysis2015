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


  MT2Analysis<MT2Estimate>* zll_ratio = new MT2Analysis<MT2Estimate>( "zll_ratio", regionsSet );
  

  int counter_total_bins = 0;
  int counter_empty_bins = 0;

  //Loop over regions

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



    //Filling the estimate
    TH1D* thisZll  = Zll ->get(*iMT2)->yield;

    int nBins =  thisZll ->GetNbinsX();
    counter_total_bins += nBins;

    for(int binnie = 1; binnie < nBins+1; ++binnie ){
      if( thisZll->GetBinContent(binnie) < 1 ) {
	std::cout << "Nothing in : "<< thisRegion.getName() << " at mt2 = " << thisZll->GetBinCenter(binnie)   << std::endl;
	counter_empty_bins = counter_empty_bins + 1;
      }
    }



    
    TH1D* thisGamma = gamma->get(*iMT2)->yield;

    TH1D* g_mt2 = new TH1D(*thisGamma); // so that it gets the same binning
    g_mt2->SetName("g_mt2");
    TH1D* h_mt2  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2->SetName("h_mt2");

    zllT ->Project( "h_mt2" , "zll_mt2", "weight*(( abs(Z_mass-91.19)<20))");
    gammaT->Project( "g_mt2", "gamma_mt2", "weight*(prompt==2)");
  
 
    //Filling the MT2Estimate with the ratio
    for( int ibin=1; ibin< nBins+1; ++ibin ) {
      zll_ratio->get(*iMT2)->yield->SetBinContent( ibin, h_mt2->GetBinContent(ibin) );
    }

    int nBinss =  h_mt2->GetNbinsX();
    for(int binnie = 1; binnie <= nBinss; binnie++){
      double value = h_mt2->GetBinContent(binnie);
      h_mt2->SetBinError(binnie,sqrt(value));
      double value_g = g_mt2->GetBinContent(binnie);
      g_mt2->SetBinError(binnie,sqrt(value_g));
    }
  
    h_mt2->Divide(g_mt2);


    //Histos of MT2, ZLL
    //   TH1D* h_mt2 = new TH1D("h_mt2","h_mt2",10, 0, 1500);
    //   zllT->Project("h_mt2", "zll_mt2","weight*(( abs(Z_mass-91.19)<20))");
    TH1D* h_mt2_450  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_450->SetName("h_mt2_450");
    //   TH1D* h_mt2_450 = new TH1D("h_mt2_450","h_mt2_450",10, 0, 1500); 
    zllT->Project("h_mt2_450", "zll_mt2","weight*(zll_ht>450&&zll_ht<575 &&( abs(Z_mass-91.19)<20))");
    TH1D* h_mt2_575  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_575->SetName("h_mt2_575");
    //   TH1D* h_mt2_575 = new TH1D("h_mt2_575","h_mt2_575",10, 0, 1500); 
    zllT->Project("h_mt2_575", "zll_mt2","weight*(zll_ht>575&&zll_ht<1000 &&( abs(Z_mass-91.19)<20))");
    //1000 < HT < 1500
    TH1D* h_mt2_1000  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_1000->SetName("h_mt2_1000");
    //    TH1D* h_mt2_1000 = new TH1D("h_mt2_1000","h_mt2_1000",10, 0, 1500); 
    zllT->Project("h_mt2_1000", "zll_mt2","weight*(zll_ht>1000&&zll_ht<1500 &&( abs(Z_mass-91.19)<20))");
    //       HT > 1500
    TH1D* h_mt2_1500  = new TH1D(*thisZll); // so that it gets the same binning
    h_mt2_1500->SetName("h_mt2_1500");
    //  TH1D* h_mt2_1500 = new TH1D("h_mt2_1500","h_mt2_1500",10, 0, 1500); 
    zllT->Project("h_mt2_1500", "zll_mt2","weight*(zll_ht>1500 &&( abs(Z_mass-91.19)<20))");
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
    h_mt2->SetMarkerStyle(20);   h_mt2->SetMarkerSize(1.6); h_mt2->SetLineColor(kBlack);
    
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

    TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
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

    zllT ->Project( "h_mt2_mc" , "zll_mt2", "weight*(( abs(Z_mass-91.19)<20))" );
    gammaT->Project( "g_mt2_mc", "gamma_mt2", "weight*(prompt==2)" );
  

    /*
    //Double ratio
    TH1D* h_mt2_mc = new TH1D("h_mt2_mc","h_mt2_mc",10, 0, 1500); h_mt2_mc->Sumw2();
    zllT->Project("h_mt2_mc", "zll_mt2","weight*(( abs(Z_mass-91.19)<20))");
 
    //Histos of MT2, gamma
    //Inclusive
    TH1D* g_mt2_mc = new TH1D("g_mt2_mc","",10, 0, 1500); g_mt2_mc->Sumw2();
    gammaT->Project("g_mt2_mc","gamma_mt2","weight*(prompt==2)"); 
    */

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






 











    delete h_mt2;   delete h_mt2_mc;    delete g_mt2; delete g_mt2_mc;
    delete h2_axes; delete h2_axes_rat;

    //   delete h_mt2_450; delete h_mt2_575;    delete h_mt2_1000;    delete h_mt2_1500;
    //   delete g_mt2_450; delete g_mt2_575;    delete g_mt2_1000;    delete g_mt2_1500;

    


    //FOR THE HT plot
    float bins_ht[] = {450,575,1000,1500, 2000};
    int binNum_ht =  sizeof(bins_ht)/sizeof(float) -1;

    TH1D* h_ht2 = new TH1D("h_ht2","h_ht2",binNum_ht, bins_ht);
    zllT->Project("h_ht2", "zll_ht","weight*(abs(Z_mass-91.19)<20)");
    TH1D* g_ht2 = new TH1D("g_ht2","g_ht2",binNum_ht, bins_ht);
    
    gammaT->Project("g_ht2", "gamma_ht","weight*(prompt==2)");
    
    int nBins2 =  h_ht2->GetNbinsX();
    for(int binnie = 1; binnie < nBins2+1; binnie++){
      double value = h_ht2->GetBinContent(binnie);
      h_ht2->SetBinError(binnie,sqrt(value));
      double value_g = g_ht2->GetBinContent(binnie);
      g_ht2->SetBinError(binnie,sqrt(value_g));
    }
    h_ht2->Divide(g_ht2); 
    h_ht2->SetMarkerStyle(20);      h_ht2->SetMarkerSize(1.6);
    h_ht2->SetMarkerColor(kBlack);  h_ht2->SetLineColor(kBlack);
 
    //The histo with the wrong errors
    TH1D* h_ht = new TH1D("h_ht","h_ht",binNum_ht, bins_ht);
    TH1D* g_ht = new TH1D("g_ht","g_ht",binNum_ht, bins_ht);  
    zllT->Project("h_ht", "zll_ht","weight*(abs(Z_mass-91.19)<20)");
    gammaT->Project("g_ht", "gamma_ht","weight*(prompt==2)");
 
    h_ht->Sumw2(); g_ht->Sumw2();
    h_ht->Divide(g_ht);    h_ht->SetMarkerStyle(20);    h_ht->SetMarkerSize(1.6);
 

    TH1D* h_ht_mc = new TH1D("h_ht_mc","h_ht_mc",binNum_ht, bins_ht);
    TH1D* g_ht_mc = new TH1D("g_ht_mc","g_ht_mc",binNum_ht, bins_ht);
 
    //Double ratioooooo
    //    TH1D* h_ht_mc = new TH1D("h_ht_mc","h_ht_mc",20, 0, 2000); 
    h_ht_mc->Sumw2();  
    zllT->Project("h_ht_mc", "zll_ht","weight*(abs(Z_mass-91.19)<20)");
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

    legend->Clear();
    legend->AddEntry( h_ht2 ,"HT", "P" );
    //  legend->AddEntry( h_ht2 ,"looped over bins", "P" );
    //   legend->AddEntry( h_ht_mc, "HT MC", "P");
    //   legend->Draw("same");
    
    c1->SaveAs( Form("%s/ht_%s.png", outputdir.c_str() , thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/ht_%s.eps", outputdir.c_str() , thisRegion.getName().c_str()) );
    

    TH2D* h2_axes_ht_mc = new TH2D("axes_ht_mc", "", 10, 0, 2099, 10, 0., 2 );
    h2_axes_ht_mc->SetXTitle("HT_MC [GeV]");
    h2_axes_ht_mc->SetYTitle("(Zll/ #gamma)_{Data} / (Zll/ #gamma)_{MC}");
    h2_axes_ht_mc->Draw();
     
    h_ht2->Divide(h_ht_mc);
    h_ht2->Draw("p same");
    labelTop->Draw("same");

    legend->Clear();
    legend->AddEntry( h_ht2 ,"HT", "P" );
    // legend->Draw("same");

    c1->SaveAs( Form("%s/ht_double.png", outputdir.c_str()) );
    c1->SaveAs( Form("%s/ht_double.eps", outputdir.c_str()) );
  

    delete h2_axes_ht; 
    delete h_ht; delete g_ht; delete h_ht_mc; delete g_ht_mc;

    
  
    float bins_nJet40[] = {2,4,7,12};
    int binNum_nJet40 =  sizeof(bins_nJet40)/sizeof(float) -1;

    TH1D* h_nJet40 = new TH1D("h_nJet40","h_nJet40",binNum_nJet40, bins_nJet40);
   
    //Histos of NJET40, ZLL
    zllT->Project("h_nJet40", "nJets","weight*(( abs(Z_mass-91.19)<20))");
    //    h_nJet40->Sumw2();

    //Histos of NJET40, gamma
    TH1D* g_nJet40 = new TH1D("g_nJet40","g_nJet40",binNum_nJet40, bins_nJet40);
    gammaT->Project("g_nJet40","gamma_nJets","weight*(prompt==2)"); 
    //   g_nJet40->Sumw2();


    int nBins_nJet40 =  h_nJet40->GetNbinsX();
    h_nJet40->SetBinContent(nBins_nJet40, h_nJet40->GetBinContent(nBins_nJet40)+h_nJet40->GetBinContent(nBins_nJet40+1));
    g_nJet40->SetBinContent(nBins_nJet40, g_nJet40->GetBinContent(nBins_nJet40)+g_nJet40->GetBinContent(nBins_nJet40+1));
    for(int binnie = 1; binnie < nBins_nJet40+1; binnie++){
      double value = h_nJet40->GetBinContent(binnie);
      h_nJet40->SetBinError(binnie,sqrt(value));
      double value_g = g_nJet40->GetBinContent(binnie);
      g_nJet40->SetBinError(binnie,sqrt(value_g));
    }

    h_nJet40->Divide(g_nJet40); h_nJet40->SetMarkerStyle(20);   h_nJet40->SetMarkerSize(1.6);
    
    TH1D* h_nJet40_450 = new TH1D("h_nJet40_450","h_nJet40_450",binNum_nJet40, bins_nJet40);
    TH1D* h_nJet40_575 = new TH1D("h_nJet40_575","h_nJet40_575",binNum_nJet40, bins_nJet40);
    TH1D* h_nJet40_1000 = new TH1D("h_nJet40_1000","h_nJet40_1000",binNum_nJet40, bins_nJet40);
    TH1D* h_nJet40_1500 = new TH1D("h_nJet40_1500","h_nJet40_1500",binNum_nJet40, bins_nJet40);
    zllT->Project("h_nJet40_450", "nJets","weight*(zll_ht>450&&zll_ht<575 &&( abs(Z_mass-91.19)<20))");
    zllT->Project("h_nJet40_575", "nJets","weight*(zll_ht>575&&zll_ht<1000 &&( abs(Z_mass-91.19)<20))");
    zllT->Project("h_nJet40_1000", "nJets","weight*(zll_ht>1000&&zll_ht<1500 &&( abs(Z_mass-91.19)<20))");
    zllT->Project("h_nJet40_1500", "nJets","weight*(zll_ht>1500 &&( abs(Z_mass-91.19)<20))");
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
   
    c1->SaveAs( Form("%s/nJet40_%s.eps", outputdir.c_str() , thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/nJet40_%s.png", outputdir.c_str() , thisRegion.getName().c_str()) );



  
    delete h_nJet40;
    delete g_nJet40;
    delete h_nJet40_450; delete h_nJet40_575; delete h_nJet40_1000;    delete h_nJet40_1500;
    delete g_nJet40_450; delete g_nJet40_575; delete g_nJet40_1000;    delete g_nJet40_1500;
    delete h2_axes_jets;









    float bins_nBJet20[] = {0,1,2,3,6};
    int binNum_nBJet20 =  sizeof(bins_nBJet20)/sizeof(float) -1;

    TH1D* h_nBJet20 = new TH1D("h_nBJet20","h_nBJet20",binNum_nBJet20, bins_nBJet20);
    TH1D* g_nBJet20 = new TH1D("g_nBJet20","g_nBJet20",binNum_nBJet20, bins_nBJet20);
   
   

    //Histos of NJET40, ZLL
    zllT->Project("h_nBJet20", "nBJets","weight*(( abs(Z_mass-91.19)<20))");
    //   h_nBJet20->Sumw2(); 
    //Histos of NJET40, gamma
    //Inclusive
    gammaT->Project("g_nBJet20","gamma_nBJets","weight*(prompt==2)");
    //  g_nBJet20->Sumw2();

    int nBins_nBJet20 =  h_nBJet20->GetNbinsX();
    h_nBJet20->SetBinContent(nBins_nBJet20, h_nBJet20->GetBinContent(nBins_nBJet20)+h_nBJet20->GetBinContent(nBins_nBJet20+1));
    g_nBJet20->SetBinContent(nBins_nBJet20, g_nBJet20->GetBinContent(nBins_nBJet20)+g_nBJet20->GetBinContent(nBins_nBJet20+1));
    for(int binnie = 1; binnie < nBins_nBJet20+1; binnie++){
      double value = h_nBJet20->GetBinContent(binnie);
      h_nBJet20->SetBinError(binnie,sqrt(value));
      double value_g = g_nBJet20->GetBinContent(binnie);
      g_nBJet20->SetBinError(binnie,sqrt(value_g));
    }

    h_nBJet20->Divide(g_nBJet20); h_nBJet20->SetMarkerStyle(20);h_nBJet20->SetMarkerSize(1.6);
    

    TH1D* h_nBJet20_450 = new TH1D("h_nBJet20_450","h_nBJet20_450",binNum_nBJet20, bins_nBJet20);
    TH1D* h_nBJet20_575 = new TH1D("h_nBJet20_575","h_nBJet20_575",binNum_nBJet20, bins_nBJet20);
    TH1D* h_nBJet20_1000 = new TH1D("h_nBJet20_1000","h_nBJet20_1000",binNum_nBJet20, bins_nBJet20);
    TH1D* h_nBJet20_1500 = new TH1D("h_nBJet20_1500","h_nBJet20_1500",binNum_nBJet20, bins_nBJet20);
   
    TH1D* g_nBJet20_450 = new TH1D("g_nBJet20_450","g_nBJet20_450",binNum_nBJet20, bins_nBJet20);
    TH1D* g_nBJet20_575 = new TH1D("g_nBJet20_575","g_nBJet20_575",binNum_nBJet20, bins_nBJet20);
    TH1D* g_nBJet20_1000 = new TH1D("g_nBJet20_1000","g_nBJet20_1000",binNum_nBJet20, bins_nBJet20);
    TH1D* g_nBJet20_1500 = new TH1D("g_nBJet20_1500","g_nBJet20_1500",binNum_nBJet20, bins_nBJet20);
    

    
    zllT->Project("h_nBJet20_450", "nBJets","weight*(zll_ht>450&&zll_ht<575 &&( abs(Z_mass-91.19)<20))");
    zllT->Project("h_nBJet20_575", "nBJets","weight*(zll_ht>575&&zll_ht<1000 &&( abs(Z_mass-91.19)<20))");
    zllT->Project("h_nBJet20_1000", "nBJets","weight*(zll_ht>1000&&zll_ht<1500 &&( abs(Z_mass-91.19)<20))");
    zllT->Project("h_nBJet20_1500", "nBJets","weight*(zll_ht>1500 &&( abs(Z_mass-91.19)<20))");

    h_nBJet20_450->Sumw2();h_nBJet20_575->Sumw2(); h_nBJet20_1000->Sumw2();h_nBJet20_1500->Sumw2();
  
    gammaT->Project("g_nBJet20_450","gamma_nBJets","weight*(prompt==2&&gamma_ht>450 && gamma_ht<575)"); 
    gammaT->Project("g_nBJet20_575","gamma_nBJets","weight*(prompt==2&&gamma_ht>575 && gamma_ht<1000)");   gammaT->Project("g_nBJet20_1000","gamma_nBJets","weight*(prompt==2&&gamma_ht>1000 && gamma_ht<1500)"); 
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
   
    c1->SaveAs( Form("%s/nBJet20_%s.eps", outputdir.c_str() , thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/nBJet20_%s.png", outputdir.c_str() , thisRegion.getName().c_str()) );


    delete h_nBJet20;
    delete g_nBJet20;
    delete h_nBJet20_450; delete h_nBJet20_575; delete h_nBJet20_1000; delete h_nBJet20_1500;
    delete g_nBJet20_450; delete g_nBJet20_575; delete g_nBJet20_1000; delete g_nBJet20_1500;
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
  std::cout << counter_total_bins << std::endl;
  std::cout << counter_empty_bins << std::endl;


  return 0;
}
