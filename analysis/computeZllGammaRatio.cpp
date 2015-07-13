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

//float lumi = 0.1;
//float lumi = 4.;



void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string zll_sel,  std::string gamma_sel, MT2Analysis<MT2Estimate>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma, MT2Analysis<MT2EstimateTree>*  Zll, MT2Analysis<MT2Estimate>*  zll_yield, const MT2Region thisRegion, std::string cut, std::string cut_gamma, float lumi);

void drawCorrelation(std::string fullPath, float *binss, unsigned int size,  float *binss2, unsigned int size2,  std::string zll_sel,  std::string zll_sel2, MT2Analysis<MT2EstimateTree>*  Zll,  const MT2Region thisRegion, std::string cut , float lumi);


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
  MT2Config cfg( configFileName);
  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat"; 
  std::string samples = cfg.mcSamples();

  std::string outputdir = "ZllGamma_Ratio_" + configFileName ;

  double intpart;
  double fracpart = modf(cfg.lumi(), &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("_%.0ffb", intpart ) );
  //  outputdir += suffix;
  
  system(Form("mkdir -p %s", outputdir.c_str()));
 
  std::cout << "-> Using regions: " << regionsSet << std::endl;



  gStyle->SetOptTitle(0);
  MT2DrawTools::setStyle();
 
  std::string gammaControlRegionDir = "EventYields_" + configFileName + "_dummy/gammaControlRegion" ;
  //  std::string gammaControlRegionDir = "GammaControlRegion_" + samples + "_" + regionsSet ;
 
  MT2Analysis<MT2EstimateTree>* gamma = MT2Analysis<MT2EstimateTree>::readFromFile(gammaControlRegionDir + "/data.root", "gammaCRtree");
  if( gamma==0 ) {
    std::cout << "-> Please run gammaControlRegion first. I need to get the gammaCR yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(193);
  }


  std::string ZllDir = "Zll_CR_" + configFileName;

  //  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb/Zll_analyses.root", ZllDir.c_str(), lumi), "DYJets");
  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/Zll_analyses.root", ZllDir.c_str()) , "DYJets");
  if( Zll==0 ) {
    std::cout << "-> Please run computeZinvFromZll first. I need to get the Z->ll MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }


  MT2Analysis<MT2Estimate>* zll_ratio = new MT2Analysis<MT2Estimate>( "zll_ratio", "13TeV_inclusive" );
  MT2Analysis<MT2Estimate>* zll_yield = new MT2Analysis<MT2Estimate>( "zll_yield", "13TeV_inclusive" );

  //YIELDS
  MT2Analysis<MT2Estimate>* zllY_pt = new MT2Analysis<MT2Estimate>( "zllY_pt", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllY_mt2 = new MT2Analysis<MT2Estimate>( "zllY_mt2", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllY_ht = new MT2Analysis<MT2Estimate>( "zllY_ht", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllY_nJets = new MT2Analysis<MT2Estimate>( "zllY_nJets","13TeV_inclusive"); 
  MT2Analysis<MT2Estimate>* zllY_nBJets = new MT2Analysis<MT2Estimate>("zllY_nBJets","13TeV_inclusive");

  //RATIOS
  MT2Analysis<MT2Estimate>* zllG_pt = new MT2Analysis<MT2Estimate>( "zllG_pt", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllG_mt2 = new MT2Analysis<MT2Estimate>( "zllG_mt2", "13TeV_inclusive" ); 
  MT2Analysis<MT2Estimate>* zllG_ht = new MT2Analysis<MT2Estimate>( "zllG_ht", "13TeV_inclusive"); 
  MT2Analysis<MT2Estimate>* zllG_nJets = new MT2Analysis<MT2Estimate>( "zllG_nJets","13TeV_inclusive"); 
  MT2Analysis<MT2Estimate>* zllG_nBJets = new MT2Analysis<MT2Estimate>("zllG_nBJets","13TeV_inclusive");
  

  int counter_total_bins = 0;
  int counter_empty_bins = 0;

  //Loop over regions

  //std::set<MT2Region> MT2Regions = zll_ratio->getRegions();
  std::set<MT2Region> MT2Regions = Zll->getRegions();
  std::set<MT2Region> Inclusive = Zll->getRegions();
  MT2Region RegIncl( (*Inclusive.begin()) );

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  TH1F::AddDirectory(kTRUE);

  //THE TREES
  TTree *zllT =  Zll->get(RegIncl)->tree;
  TTree *gammaT =  gamma->get(RegIncl)->tree;

  TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi() );
 

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
    // std::set<MT2Region>::iterator iMT2 = MT2Regions.begin();  
    std::vector<std::string> niceNames = thisRegion.getNiceNames();



    float bins_nJets[] = {2,4,7,12};
    float bins_nBJets[] = {0,1,2,3,6};
  

    float bins_nJets2[] = {2,3,4,5,6,7,8,9,10,11,12};
    float bins_nBJets2[] = {0,1,2,3,4,5,6};
    float bins_mt22[] = {200,250,300,350, 400,450,500,550, 600,650,700, 800, 900,1000, 1100, 1200,1300, 1400, 1500, 1600,1700 };
    float bins_ht2[] =  { 450,500,550, 600,650,700,750,  800, 900, 1000,1100, 1200, 1300,1400, 1500,1600};


    //in MT2
    // float bins_mt2[] = {200,300,400,500, 600,  1500 };

    //  float bins_mt2[] = {200,300,400,500, 600, 800,  1500 };
    float bins_mt2[] = {200,300,400,500, 600, 800, 1000, 1500 };
    //in HT
    float bins_ht[] =  {450,575,1000,1500,2000};
 
    //    float bins_mt2[] = {200,300,400,500, 600, 800, 1000, 1500 , 1900};
    std::string cut =  "weight*(abs(Z_mass-91.19)<10 &&nBJets<2)";
    std::string cut_corr =  "weight*(abs(Z_mass-91.19)<10 )";
    std::string cut_gamma =  "weight*(prompt==2)";

    int size_mt2 = sizeof(bins_mt22)/sizeof(float)-1;
    int size_ht = sizeof(bins_ht2)/sizeof(float)-1;
    int size_nJets = sizeof(bins_nJets2)/sizeof(float)-1;
    int size_nBJets = sizeof(bins_nBJets2)/sizeof(float)-1;
  
    float lumi = cfg.lumi();    

    //draw ratio also fills the ratio and yield estimates
    drawRatios( outputdir, bins_mt2, sizeof(bins_mt2)/sizeof(float)-1 , "mt2", "mt2" ,  zllG_mt2,   gamma,  Zll,  zllY_mt2, thisRegion, cut, cut_gamma, lumi  );
   
    drawRatios( outputdir, bins_ht, sizeof(bins_ht)/sizeof(float)-1 , "ht", "ht" ,  zllG_ht,   gamma,  Zll,  zllY_ht, thisRegion, cut, cut_gamma, lumi );
    
    drawRatios( outputdir, bins_nJets, sizeof(bins_nJets)/sizeof(float)-1 , "nJets", "nJets" ,  zllG_nJets,   gamma,  Zll,  zllY_nJets, thisRegion, cut, cut_gamma, lumi );
    drawRatios( outputdir, bins_nBJets, sizeof(bins_nBJets)/sizeof(float)-1 , "nBJets", "nBJets" ,  zllG_nBJets,   gamma,  Zll,  zllY_nBJets, thisRegion, cut, cut_gamma , lumi);
   
 

    drawCorrelation( outputdir , bins_mt22, size_mt2 ,  bins_ht2, size_ht,  "mt2", "ht",  Zll, thisRegion,  cut_corr, lumi );
    drawCorrelation( outputdir , bins_mt22, size_mt2 ,  bins_nJets2, size_nJets,  "mt2", "nJets",  Zll, thisRegion,  cut_corr , lumi);
    drawCorrelation( outputdir , bins_mt22, size_mt2 ,  bins_nBJets2, size_nBJets,  "mt2", "nBJets",  Zll, thisRegion,  cut_corr , lumi);

    drawCorrelation( outputdir , bins_ht2, size_ht  ,  bins_nJets2, size_nJets ,  "ht", "nJets",  Zll, thisRegion,  cut_corr , lumi);
//  drawCorrelation( outputdir , bins_ht,sizeof(bins_ht)/sizeof(float)-1   ,  bins_nJets, sizeof(bins_nJets)/sizeof(float)-1  ,  "ht", "nJets",  Zll, thisRegion,  cut_corr );
    drawCorrelation( outputdir , bins_ht2, size_ht ,  bins_nBJets2, size_nBJets,  "ht", "nBJets",  Zll, thisRegion,  cut_corr, lumi );

    drawCorrelation( outputdir , bins_nJets2, size_nJets ,  bins_nBJets2, size_nBJets,  "nJets", "nBJets",  Zll, thisRegion,  cut_corr, lumi );
  
   


  }
  

  //FOR THE Z MASS PLOT
  TH1D* h_Zmass = new TH1D("h_Zmass","h_Zmass", 60, 60, 120);  //  h_Zmass->Sumw2();
  zllT->Project("h_Zmass", "Z_mass","weight");
  //  h_Zmass->SetMarkerStyle(20);
  //  h_Zmass->SetMarkerSize(1.6);
  h_Zmass->SetLineWidth(2);
  h_Zmass->SetLineColor(kBlue+1);

  double yMaximum = h_Zmass->GetMaximum()*1.1;
  TH2D* h2_axes2 = new TH2D("axes2", "", 10, 60, 120 , 10, 0., yMaximum );
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


















void drawRatios(std::string fullPath, float *binss, unsigned int size,  std::string zll_sel,  std::string gamma_sel, MT2Analysis<MT2Estimate>*  zll_ratio,  MT2Analysis<MT2EstimateTree>*  gamma, MT2Analysis<MT2EstimateTree>*  Zll, MT2Analysis<MT2Estimate>*  zll_yield, const MT2Region thisRegion, std::string cut, std::string cut_gamma, float lumi){
 
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

  //THE TREES
  TTree *zllT =  Zll->get(thisRegion)->tree;
  TTree *gammaT =  gamma->get(thisRegion)->tree;

  float meanX[size+1];
  float meanX_err[size+1];

  for(unsigned int i = 0; i < size; ++i){
    TH1D* histo_x= new TH1D("histo_x","",100, bins[i], bins[i+1] );
    zllT->Project( "histo_x",zll_sel.c_str(), cut.c_str());
    histo_x->GetMean();
    meanX[i] =  histo_x->GetMean();
    meanX_err[i] =  histo_x->GetMeanError();
    /*   std::cout << histo_x->GetMean() << std::endl; 
    std::cout << histo_x->GetMeanError() << std::endl; 
    std::cout << histo_x->GetSumOfWeights() << std::endl; 
    std::cout << histo_x->Integral() << std::endl; 
    */    delete histo_x;
  }


  TCanvas* canny = new TCanvas( "canny", "", 600, 700 );
  canny->cd();

  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();


  TH1D* h_mt2 = new TH1D("h_mt2","", size  , bins);
  TH1D* g_mt2 = new TH1D("g_mt2","", size  , bins);

  h_mt2->Sumw2(); g_mt2->Sumw2();
   
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

    TGraphErrors *gr_ratio = new TGraphErrors(0);
    for( unsigned int k=0; k< size; ++k){
     if(h_mt2->GetBinContent(k+1)<0.01){
      gr_ratio->SetPoint(k, 9000 , 900);
      }else{  gr_ratio->SetPoint(k, meanX[k], h_mt2->GetBinContent(k+1));
      gr_ratio->SetPointError(k, meanX_err[k],  h_mt2->GetBinError(k+1) );
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

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 0.3 );
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

    h_mt2_mc->Draw("hist same");
 
    //    h_mt2->DrawClone("p same");
    gr_ratio->Draw("same P");
 
    labelTop->Draw("same");
    gPad->RedrawAxis();

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


 
    gPad->RedrawAxis();

    canny->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.1);
    pad2->Draw();
    pad2->cd();
      
    h_mt2->Divide(h_mt2_mc); //the ratio is now in h_mt2

   TGraphErrors *gr_ratioD = new TGraphErrors(0);
    for( unsigned int k=0; k< size; ++k){
     
      gr_ratioD->SetPoint(k, meanX[k], h_mt2->GetBinContent(k+1));
      gr_ratioD->SetPointError(k, 5000, h_mt2->GetBinError(k+1));
    } 
    gr_ratioD->SetMarkerSize(1.4);
    gr_ratioD->SetMarkerStyle(20);
    gr_ratioD->SetLineColor(kBlack);
    gr_ratioD->SetLineWidth(2);
    gr_ratioD->SetMarkerColor(kBlack);


    TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , 0.,2  );
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

    canny->cd();
    canny->SaveAs( Form("%s/%s_ratios_%s.eps", fullPath.c_str(), zll_sel.c_str(),  thisRegion.getName().c_str() ) );
    canny->SaveAs( Form("%s/%s_ratios_%s.png", fullPath.c_str(), zll_sel.c_str(),thisRegion.getName().c_str() ) );
    canny->SaveAs( Form("%s/%s_ratios_%s.pdf", fullPath.c_str(), zll_sel.c_str(),thisRegion.getName().c_str() ) );


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
  regionText->SetTextSize(0.038);
  regionText->SetTextFont(42);
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
