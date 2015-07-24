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
#include "TRandom3.h"

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2Config.h"


#include "../interface/MT2DrawTools.h"


#include <iostream>
#include "string.h"


#define mt2_cxx
#include "interface/mt2.h"



void drawMll( const std::string& outputdir,  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, MT2Analysis<MT2EstimateTree>* data, bool of, float lumi ); 

void randomizePoisson( TH1* histo );

void drawStacks(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  MT2Analysis<MT2EstimateTree>* data,const MT2Region thisRegion, std::string cut, float lumi);




int main(int argc, char* argv[]){

  std::string regionsSet = "zurich";
  if( argc>2 ) {
    regionsSet = std::string(argv[2]);
  }

  if( argc<2 ) {
    std::cout << "USAGE: ./coputeZllGammaRatio [configFileName] regionSet" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  std::string configFileName(argv[1]);
  MT2Config cfg( configFileName);
  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat"; 
  std::string samples = cfg.mcSamples();

  regionsSet = cfg.regionsSet();


  std::string outputdir( Form("ZllData_%s", configFileName.c_str() ) );
  std::string outputdir_of( Form("ZllData_OF_%s", configFileName.c_str()  ) );


  std::cout << "-> Using regions: " << regionsSet << std::endl;

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  double intpart;
  double fracpart = modf(cfg.lumi(), &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("_%.0ffb", intpart ) );
  //outputdir += suffix;
  // outputdir_of += suffix;
  
  system(Form("mkdir -p %s", outputdir.c_str()));






  std::string ZllDir = "ZllPurity_" + configFileName;
  std::string ZllDir_of = "ZllPurity_OF_" + configFileName;


  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "DYJets");
  if( Zll==0 ) {
    std::cout << "-> Please run zllPurityTrees first. I need to get the yields from there." << std::endl;    std::cout << "-> Thank you for your cooperation." << std::endl;    exit(197);
  }



  
  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str()  ), "QCD");
  MT2Analysis<MT2EstimateTree>* top = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "Top");
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "WJets");
  MT2Analysis<MT2EstimateTree>* zjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees.root", ZllDir.c_str() ), "ZJets");

  
  //OPPOSITE FLAVOR TREES
  MT2Analysis<MT2EstimateTree>* Zll_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "DYJets");

  MT2Analysis<MT2EstimateTree>* qcd_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str()  ), "QCD");
  MT2Analysis<MT2EstimateTree>* top_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "Top");
  MT2Analysis<MT2EstimateTree>* wjets_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "WJets");
  MT2Analysis<MT2EstimateTree>* zjets_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_of.root", ZllDir_of.c_str() ), "ZJets");
  
  MT2Analysis<MT2EstimateTree>* data_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_data_of.root", ZllDir_of.c_str() ) , "data_of");

  Zll_of->setFullName("Z+jets");
  wjets_of->setFullName("W+jets");
  zjets_of->setFullName("Z#nu#nu+jets");

  data_of->setFullName("Data");

  

  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/ZllPurityTrees_data.root", ZllDir.c_str() ) , "data");


  data->setFullName("Data");

  Zll->setFullName("Z+jets");
  wjets->setFullName("W+jets");
  zjets->setFullName("Z#nu#nu+jets");



  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
  bgYields.push_back( Zll );
  // bgYields.push_back( qcd );
  //  bgYields.push_back( wjets );
  //  bgYields.push_back( zjets );
  bgYields.push_back( top );
  
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_of; 
  bgYields_of.push_back( Zll_of );
  //  bgYields_of.push_back( qcd_of );
  //  bgYields_of.push_back( wjets_of );
  //  bgYields_of.push_back( zjets_of );
  bgYields_of.push_back( top_of );
  
  drawMll( outputdir, bgYields,data,  0 , cfg.lumi() );
  //drawMll( outputdir_of, bgYields_of, data_of,  1, cfg.lumi() );

  return 0;
}




















void drawMll( const std::string& outputdir, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields,  MT2Analysis<MT2EstimateTree> * data, bool of , float lumi) {

  MT2DrawTools::setStyle();

  std::vector<int> colors;
  if( bgYields.size()==3 ) { // estimates
    colors.push_back(402); 
    colors.push_back(430); 
    colors.push_back(418); 
  } else { // mc
    colors.push_back(430); // other = zll 
    //  colors.push_back(401); // qcd
    //  colors.push_back(417); // w+jets
    //    colors.push_back(419); // z+jets
    colors.push_back(855); // top
  }

  TH1F::AddDirectory(kTRUE);

  std::string fullPath = outputdir;  std::string fullPathPlots = outputdir + "/plots";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  std::set<MT2Region> MT2Regions = bgYields[0]->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    THStack bgStack("bgStack", "");
    THStack bgStack_scaled("bgStack_scaled", "");

    TH1D* histo_bg = new TH1D("histo_bg", "", 100, 0, 1000);

    histo_bg->Sumw2();

    for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
      int index = bgYields.size() - i - 1;
      //  TH1D* h1_bg = new TH1D("h1_bg","", 20, 0,250);
      //     h1_bg->Sumw2();
      TH1D* h1_bg2 = new TH1D(Form("h1_bg2_%d_%d",of,index ),"", 100 , 0,1000);
      //     if(of==1){ h1_bg->Rebin(4);
      //     }     else{ h1_bg->Rebin(4); }
      TTree *bgTree = bgYields[index]->get(*iMT2)->tree;
      //   bgTree->Project("h1_bg","Z_mass","weight");
      h1_bg2->Sumw2();
      bgTree->Project(Form("h1_bg2_%d_%d", of,  index ),"Z_mass","weight*(Z_pt>180 &&lep_pt0>25 && lep_pt1>20 )");
      //     bgTree->Project("h1_bg2","Z_mass","weight*(Z_mass>80.)");
      //     h1_bg->SetFillColor( colors[index] );
      //     h1_bg->SetLineColor( kBlack );
      histo_bg->Add(h1_bg2);
      //  bgStack.Add(h1_bg);
    }

    TH1D* h_data = new TH1D("h_data","",  200, 50,150);
    if(of==1){ h_data->Rebin(8);  }     else{ h_data->Rebin(4); } 
    TH1D* histo_data = new TH1D("histo_data","", 100, 0, 1000);  

    TTree *data_Tree = data->get(*iMT2)->tree;

    data_Tree->Project("h_data","Z_mass","weight*(Z_pt>180)");
    data_Tree->Project("histo_data","Z_mass","(Z_mass>80. &&(Z_pt>180 && lep_pt0>25 && lep_pt1>20))");
  
    float bg_int=     histo_bg->Integral();
    float data_int=     h_data->Integral();
    //  float scale = data_int/bg_int;
    // float scale =1;

    for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
      int index = bgYields.size() - i - 1;
      TH1D* h1_bg = new TH1D(Form("h1_bg_%d_%d",of, index ),"", 200, 50,150);
      h1_bg->Sumw2();
      if(of==1){ h1_bg->Rebin(8);
      }     else{ h1_bg->Rebin(4); }
      TTree *bgTree = bgYields[index]->get(*iMT2)->tree;
      bgTree->Project(Form("h1_bg_%d_%d", of, index ),"Z_mass", "weight*(Z_pt>180&&  lep_pt0>25 && lep_pt1>20)" );
      //     bgTree->Project(Form("h1_bg_%d_%d", of, index ),"Z_mass",Form("weight*%f", scale) );
      h1_bg->SetFillColor( colors[index] );
      h1_bg->SetLineColor( kBlack );
      bgStack.Add(h1_bg);
    }

    h_data->SetMarkerStyle(20);
    h_data->SetMarkerColor(kBlack);

    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.2);

    int binMin = histo_data->GetXaxis()->FindBin(80);
    int binMax = histo_data->GetXaxis()->FindBin(10000);
    int binMin_bg = histo_bg->GetXaxis()->FindBin(80);
    int binMax_bg = histo_bg->GetXaxis()->FindBin(10000);
    double int_err=0; 
    double int_err_bg =0;
    // double int_data =  histo_data->Integral(binMin, binMax);
    // double int_bg =    histo_bg->Integral(binMin, binMax) ;

    double int_data =  histo_data->IntegralAndError(binMin, binMax, int_err, "");
    double int_bg =    histo_bg->IntegralAndError(binMin_bg, binMax_bg, int_err_bg, "") ;

    std::cout << "Data = " <<  int_data <<  " +- " << int_err << "  in %: " << int_err/int_data*100 << std::endl;
    std::cout << "MC = " << int_bg  << " +- " << int_err_bg << "  in %: " << int_err_bg/int_bg*100 << std::endl;
    std::cout << "R = " <<  int_data/int_bg << " +- " << sqrt( (int_err*int_err/(int_bg*int_bg)) ) << " (DATA) " << " +- " << sqrt( (int_data*int_data*int_err_bg*int_err_bg/(int_bg*int_bg*int_bg*int_bg)) ) << " (MC) " << std::endl;

    std::cout << "Zll events " << h_data->GetEntries() << std::endl;
    std::cout << "Zll mean " << h_data->GetMean() << std::endl;

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    float yMax = 1.3*(bgStack.GetMaximum());
     float yMax_data = 1.3*(h_data->GetMaximum());
     if(yMax_data > yMax) yMax = yMax_data;
    //if(yMax < 0.3)    //   yMax = 20;    //  if(of==1) yMax=25;

    TH2D* h2_axes = new TH2D("axes", "", 10, 50,150, 10, 0., yMax );
    h2_axes->SetXTitle("M_{ll} [GeV]");
    if(of==1)    h2_axes->SetXTitle("M_{e^{#pm}#mu^{#mp}} [GeV]");
    h2_axes->SetYTitle("Entries");
    h2_axes->Draw();
   
    std::vector<std::string> niceNames = thisRegion.getNiceNames();

    for( unsigned i=0; i<niceNames.size(); ++i ) {
      float yMaxText = 0.9-(float)i*0.05;
      float yMinText = yMaxText - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMinText, 0.55, yMaxText, "brNDC" );
      regionText->SetTextSize(0.035);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      //    if(i==0)
      //	regionText->AddText( "H_{T} > 180 GeV" );
      //     else
	regionText->AddText( niceNames[i].c_str() );
      regionText->Draw("same");
    }
    

    TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    //   legend->AddEntry( h1_data, "Zll", "P" );
    legend->AddEntry( gr_data, "Data", "P" );
 
    for( unsigned i=0; i<bgYields.size(); ++i ) {  
      TH1D* h1_bg1 = bgYields[i]->get(thisRegion)->yield;
      legend->AddEntry( h1_bg1, bgYields[i]->getFullName().c_str(), "F" );
    }

    gPad->Update();

    legend->Draw("same");
    bgStack.Draw("histo same");
    //   h_data->Draw("p same");
    gr_data->Draw("p same");

    TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
    labelTop->Draw("same");
    gPad->RedrawAxis();



    if( of == true){
      c1->SaveAs( Form("%s/mll_of_%s.eps", fullPathPlots.c_str(), thisRegion.getName().c_str()));
      c1->SaveAs( Form("%s/mll_of_%s.png", fullPathPlots.c_str(), thisRegion.getName().c_str()));
      c1->SaveAs( Form("%s/mll_of_%s.pdf", fullPathPlots.c_str(), thisRegion.getName().c_str()));
    }else {
      c1->SaveAs( Form("%s/mll_%s.eps", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
      c1->SaveAs( Form("%s/mll_%s.png", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
      c1->SaveAs( Form("%s/mll_%s.pdf", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
    }


    histo_data->SetLineColor(kBlue);
    histo_data->SetLineWidth(2);
    histo_bg->SetLineWidth(2);
    histo_data->Draw("");
    histo_bg->Draw("same");

    // gPad->SetLogy();
    /*
    if( of == true){
      c1->SaveAs( Form("%s/fuck_%s.eps", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
      c1->SaveAs( Form("%s/fuck_%s.png", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
      c1->SaveAs( Form("%s/fuck_%s.pdf", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
    }
    */
    


    float bins_nVert[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20,21,22,23,24 ,25};
    
    float bins_Z_lepId[] = { 10,11,12,13,14,15};
    
    float bins_nJets[] = {2,3,4,5,6,7, 12};
    float bins_nBJets[] = {0,1,2,3,4,5,6};
    //in MT2
    float bins_mt2[] = {0,50, 100,150,200,250, 300, 400,500, 1500};
    //in HT
    float bins_ht[] =  {0,100,200,300,400,500,600,700,800,1000, 1500,2000};

    float bins_met[] = {0,50,100,150,200,250,300,350,400,450, 500, 600,700,800, 1000};
    float bins_Zpt[] = {0,50,100,150, 200,250,300,350,400,450, 500,550,600,650,700, 750,800,850,900,950, 1000};
    //  float bins_ht[] =  {450,500, 600,700, 800, 900,1000,1100, 1200, 1300, 1500,2000};


    /*    float bins_nJets[] = {2,4,7,12};
    float bins_nBJets[] = {0,1,2,3,6};
    //in MT2
    float bins_mt2[] = {200,300,400,500, 600, 800, 1000, 1500 , 1900};
    //in HT
    float bins_ht[] =  {450,575,1000,1500,2000};
    */

    std::string cut;
    std::string cut2;
 
    std::string cut3 ;
    std::string cut4;


    if(of==0){ 
      // cut =  "weight";
      cut2 =  "weight*(abs(Z_mass-91.19)<20)";
      cut = "weight*(abs(Z_mass-91.19)<15 && Z_pt>180)";
 
      cut3 = "weight*(abs(Z_mass-91.19)<20&&nBJets<2)";
      cut4 ="weight*(abs(Z_mass-91.19)<10&&nBJets<2)";
    }else{
      cut =  "weight*(Z_mass>70)";
      cut2 = "weight*(Z_mass>80)";
 
      cut3 = "weight*(Z_mass>70)";
      cut4 ="weight*(Z_mass>80)";
    }

    std::string cut_nJets3 = "weight*(abs(Z_mass-91.19)<10&&nJets>2)";
   
  


    // leave commented //comment and uncomment the other one when more data
    drawStacks( fullPathPlots,  bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields, data , thisRegion, cut , lumi );
       drawStacks( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , data, thisRegion , cut , lumi);
       drawStacks( fullPathPlots,  bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "mt2", bgYields , data, thisRegion  , cut , lumi  );
       drawStacks( fullPathPlots,  bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "ht", bgYields , data , thisRegion, cut  , lumi);
 drawStacks( fullPathPlots,  bins_met,sizeof(bins_met)/sizeof(float)-1,  "met", bgYields , data, thisRegion  , cut , lumi  );     
 drawStacks( fullPathPlots,  bins_Zpt,sizeof(bins_Zpt)/sizeof(float)-1,  "Z_pt", bgYields , data, thisRegion  , cut , lumi  );     

drawStacks( fullPathPlots,  bins_nVert,sizeof(bins_nVert)/sizeof(float)-1,  "nVert", bgYields , data, thisRegion  , cut , lumi  );     

drawStacks( fullPathPlots,  bins_Z_lepId,sizeof(bins_Z_lepId)/sizeof(float)-1,  "Z_lepId", bgYields , data, thisRegion  , cut , lumi  );     

       /*
    drawStacks( fullPathPlots,  bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields , data,thisRegion, cut2, lumi );
    drawStacks( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , data, thisRegion , cut2 , lumi);
    drawStacks( fullPathPlots,  bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "mt2", bgYields , data, thisRegion  , cut2 , lumi );
    drawStacks( fullPathPlots,  bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "ht", bgYields , data, thisRegion, cut2  , lumi );
    
       */

    delete c1;
    delete h2_axes;
    delete histo_bg;
    delete h_data; delete histo_data;

    //   }//end of loop over beees


  }// for MT2 regions

}











void drawStacks(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  MT2Analysis<MT2EstimateTree>* data, const MT2Region thisRegion, std::string cut , float lumi){
 
  std::vector<int> colors;
  colors.push_back(430); // other = zll 
  //  colors.push_back(401); // qcd
  //  colors.push_back(417); // w+jets
  //  colors.push_back(419); // z+jets
  colors.push_back(855); // top

  TH1F::AddDirectory(kTRUE);

  float bins[size+1]; for(unsigned int i=0; i<= size ; i++)      bins[i]=binss[i];
  float xMin = binss[0];
  float xMax = binss[size];

 TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
  canny->cd();
  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();
 
  // std::string cut = "weight*( abs(Z_mass-91.19)<20)";
  TH1D* h_bg = new TH1D("h_bg","",size,bins);

  THStack bgStack("bgStack", "");
  for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
    int index = bgYields.size() - i - 1;
    TH1D* h1_bg = new TH1D("h1_bg","", size  , bins);
    h1_bg->Sumw2();
    TTree *bgTree = bgYields[index]->get(thisRegion)->tree;
    bgTree->Project("h1_bg",Form("%s",name.c_str()) ,Form("%s",cut.c_str()) );
    h1_bg->SetFillColor( colors[index] );
    h1_bg->SetLineColor( kBlack );
    bgStack.Add(h1_bg);
    h_bg->Add(h1_bg);
  }
    
  TH1D* h_data = new TH1D("h_data","", size, bins);
  TTree *data_Tree = data->get(thisRegion)->tree;
  data_Tree->Project("h_data",Form("%s",name.c_str()) ,Form("%s",cut.c_str()) );
    
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerColor(kBlack);

  TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h_data);
  gr_data->SetMarkerStyle(20);
  gr_data->SetMarkerSize(1.2);
 
 
  float yMax = 1.5*(bgStack.GetMaximum());
  float yMax2 = 1.5*(h_data->GetMaximum());
  if(yMax2>yMax) yMax = yMax2;

  TH2D* h2_axes = new TH2D("axes", "", 10,bins[0] ,bins[size], 10, 0., yMax );
  if(name  == "ht")  
    h2_axes->SetXTitle("H_{T} [GeV]");
  else if(name == "mt2")  
    h2_axes->SetXTitle("M_{T2} [GeV]");
  else if(name == "Z_pt")  
    h2_axes->SetXTitle("Boson p_{T} [GeV]");
  else if(name == "met")  
    h2_axes->SetXTitle("ME_{T} [GeV]");
  else if(name == "Z_lepId")  
    h2_axes->SetXTitle("Lepton Id");
  else
    h2_axes->SetXTitle(name.c_str());
  h2_axes->SetYTitle("Entries");
  h2_axes->Draw();


  TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry(gr_data, "Data", "p");
  for( unsigned i=0; i<bgYields.size(); ++i ) {  
    TH1D* h1_bg1 = bgYields[i]->get(thisRegion)->yield;
    legend->AddEntry( h1_bg1, bgYields[i]->getFullName().c_str(), "F" );
  }

 
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);

  labelTop->Draw("same");
  legend->Draw("same");
  bgStack.Draw("histo same");
  // h_data->DrawCopy("P same");
  gr_data->Draw("p same");
  gPad->RedrawAxis();




  canny->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();

  TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , 0.0, 2.0  );
  h2_axes_rat->SetYTitle("Data / MC");

  h2_axes_rat->GetXaxis()->SetTitleSize(0.2);
  h2_axes_rat->GetXaxis()->SetTitleOffset(5);
  h2_axes_rat->GetXaxis()->SetLabelSize(0.00);
  h2_axes_rat->GetXaxis()->SetTickLength(0.09);
  h2_axes_rat->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_rat->GetYaxis()->SetTitleSize(0.2);
  h2_axes_rat->GetYaxis()->SetTitleOffset(0.34);
  h2_axes_rat->GetYaxis()->SetLabelSize(0.17);
  

  TLine *line = new TLine(xMin, 1, xMax, 1);
  line->SetLineColor(kBlack);

  h2_axes_rat->Draw();
  line->Draw("same");


  h_data->Divide(h_bg);


  h_data->Draw("p same");


  canny->cd();

  std::string extension= "";
  /*
  if(cut == "weight*(abs(Z_mass-91.19)<20)") extension ="mass";
  else if(cut == "weight*(abs(Z_mass-91.19)<10)") extension ="mass10";
  else if(cut == "weight*(abs(Z_mass-91.19)<10&&nJets>2)") extension ="nJets3";
  else if(cut == "weight*(Z_mass>70)") extension = "OF_70GeV";
  else if(cut == "weight*(Z_mass>80)") extension = "OF_80GeV";
  else  extension = "massNbJets";
  // "weight*(abs(Z_mass-91.19)<10&&nJets>2)";
  */

 
  canny->SaveAs( Form("%s/%s_%s_%s.eps", fullPath.c_str(), name.c_str(),extension.c_str() , thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_%s_%s.png", fullPath.c_str(), name.c_str(),extension.c_str() , thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_%s_%s.pdf", fullPath.c_str(), name.c_str(),extension.c_str() , thisRegion.getName().c_str()) );


 
  delete h2_axes;
  delete canny;
}















void randomizePoisson( TH1* histo ) {

  TRandom3 rand(11);;
  //  TRandom3 rand(13);


  //  std::set<MT2HTRegion> HTRegions = data->getHTRegions();
  //  std::set<MT2SignalRegion> signalRegions = data->getSignalRegions();

 
  //  for( std::set<MT2HTRegion>::iterator iHT = HTRegions.begin(); iHT!=HTRegions.end(); ++iHT ) {
  //    for( std::set<MT2SignalRegion>::iterator iSR = signalRegions.begin(); iSR!=signalRegions.end(); ++iSR ) {

  for( int ibin=1; ibin<histo->GetXaxis()->GetNbins()+1; ++ibin ) {

    int poisson_data = rand.Poisson(int(histo->GetBinContent(ibin)));
    histo->SetBinContent(ibin, poisson_data);
    histo->SetBinError(ibin, sqrt(poisson_data));
	  
  }  // for bins

  // return histo;
}



