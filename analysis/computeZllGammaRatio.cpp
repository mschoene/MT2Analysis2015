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

//float lumi = 0.1;
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


 regionsSet = cfg.regionsSet();


  std::string outputdir( Form("ZllGammaRatio_%s_%s_%.0ffb", samples.c_str(), regionsSet.c_str(), lumi ) );

  //std::string outputdir = "ZllGammaRatio_"+ cfg.mcSamples + regionsSet ;
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


  MT2Analysis<MT2EstimateSyst>* zll_ratio = new MT2Analysis<MT2EstimateSyst>( "zll_ratio", regionsSet );
  MT2Analysis<MT2EstimateSyst>* zll_mc = new MT2Analysis<MT2EstimateSyst>( "zll_mc", regionsSet );
  MT2Analysis<MT2EstimateSyst>* gamma_data = new MT2Analysis<MT2EstimateSyst>( "gamma_data", regionsSet );
  MT2Analysis<MT2EstimateSyst>* gamma_mc = new MT2Analysis<MT2EstimateSyst>( "gamma_mc", regionsSet );
  

  int counter_total_bins = 0;
  int counter_empty_bins = 0;

  //Loop over regions

  std::set<MT2Region> MT2Regions = Zll->getRegions();
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
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
	std::cout << "Yield = "<< thisZll->GetBinContent(binnie)<< " in : "<< thisRegion.getName() << " at mt2 = " << thisZll->GetBinLowEdge(binnie) << " to " << thisZll->GetBinLowEdge(binnie)+ thisZll->GetBinWidth(binnie)   << std::endl;
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
  
 


    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h_mt2);     
    TGraphAsymmErrors* gr_gamma_data = MT2DrawTools::getPoissonGraph(g_mt2);     



    int nBinss =  h_mt2->GetNbinsX();
    for(int binnie = 1; binnie <= nBinss; binnie++){

      double value = h_mt2->GetBinContent(binnie);
      //   h_mt2->SetBinError(binnie,sqrt(value));
      double value_g = g_mt2->GetBinContent(binnie);
      //    g_mt2->SetBinError(binnie,sqrt(value_g));
      zll_ratio->get(*iMT2)->yield->SetBinContent( binnie,   value);
      gamma_data->get(*iMT2)->yield->SetBinContent( binnie, value_g);


      zll_mc->get(*iMT2)->yield->SetBinContent( binnie,   value);
      gamma_mc->get(*iMT2)->yield->SetBinContent( binnie, value_g);

  

      double errYup = gr_data->GetErrorYhigh(binnie-1);
      double errYlow = gr_data->GetErrorYlow(binnie-1);
      zll_ratio->get(*iMT2)->yield_systUp->SetBinContent( binnie,  errYup );
      zll_ratio->get(*iMT2)->yield_systDown->SetBinContent( binnie, errYlow );
      zll_mc->get(*iMT2)->yield_systUp->SetBinContent( binnie,  0 );
      zll_mc->get(*iMT2)->yield_systDown->SetBinContent( binnie, 0 );

      double errYup_g = gr_gamma_data->GetErrorYhigh(binnie-1);
      double errYlow_g  = gr_gamma_data->GetErrorYlow(binnie-1);
      gamma_data->get(*iMT2)->yield_systUp->SetBinContent( binnie,  errYup_g  );
      gamma_data->get(*iMT2)->yield_systDown->SetBinContent( binnie, errYlow_g );
      gamma_mc->get(*iMT2)->yield_systUp->SetBinContent( binnie, 0 );
      gamma_mc->get(*iMT2)->yield_systDown->SetBinContent( binnie, 0 );

    }
  


    /*


    //   h_mt2->Sumw2(); g_mt2->Sumw2();
    
    h_mt2->Divide(g_mt2);

    h_mt2->SetMarkerStyle(20);   h_mt2->SetMarkerSize(1.6);
 
    TH2D* h2_axes = new TH2D("axes", "", 10, 0, 1500, 10, 0., 1 );
    h2_axes->SetXTitle("M_{T2} [GeV]");
    h2_axes->SetYTitle("Zll / #gamma Ratio");
    h2_axes->Draw();


    TLegend* legend = new TLegend( 0.2, 0.9-(5)*0.06, 0.5, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( h_mt2 ,"HT inclusive", "P" );

    h_mt2->Draw("p same");

    TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
    labelTop->Draw("same");
    // legend->Draw("same");
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



    */


    /*
    //Filling the MT2Estimate with the ratio
    for( int ibin=1; ibin< nBins+1; ++ibin ) {
      zll_ratio->get(*iMT2)->yield->SetBinContent( ibin, h_mt2->GetBinContent(ibin) );
      zll_ratio->get(*iMT2)->yield->SetBinError( ibin, h_mt2->GetBinError(ibin) );

    }
    */

    delete h_mt2; //  delete h_mt2_mc; 
    delete g_mt2;// delete g_mt2_mc;
    //  delete h2_axes; delete h2_axes_rat;
    delete c1;


    zll_ratio->finalize();

  }
 


  std::cout << counter_total_bins << std::endl;
  std::cout << counter_empty_bins << std::endl;


  std::string outFile = outputdir + "/ZllRatio.root";
  *zll_ratio = *zll_ratio/(*gamma_data);
  *zll_mc = *zll_mc/(*gamma_mc);
  *zll_ratio = *zll_ratio/(*zll_mc);

  zll_ratio->writeToFile( outFile );



  return 0;
}
