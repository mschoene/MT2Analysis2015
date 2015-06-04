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

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2Config.h"



#include "../interface/MT2DrawTools.h"


#include <iostream>


#define mt2_cxx
#include "interface/mt2.h"


float lumi = 0.1;
//float lumi = 4.;


void drawMll( const std::string& outputdir,  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields );

void drawSFvsOF( const std::string& outputdir,  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_of );

void drawStacks(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, const MT2Region thisRegion, std::string cut);


void drawPurity(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, const MT2Region thisRegion, std::string cut);


void drawPurityComparison(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, const MT2Region thisRegion, std::string cut, std::string cut2);

MT2Analysis<MT2EstimateTree>* computeYield( const MT2Sample& sample, const MT2Config& cfg, float lumi=1. );
MT2Analysis<MT2EstimateTree>* mergeYields( std::vector< MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max=-1, const std::string& legendName="" );


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
  MT2Config cfg("cfgs/" + configFileName + ".txt");
  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat"; 
  std::string samples = cfg.mcSamples();

 regionsSet = cfg.regionsSet();


  std::string outputdir( Form("ZllPurityFigures_%s_%s_%.0ffb", samples.c_str(), regionsSet.c_str(), lumi ) );
  std::string outputdir_of( Form("ZllPurityFigures_OF_%s_%s_%.0ffb", samples.c_str(), regionsSet.c_str(), lumi ) );


  std::cout << "-> Using regions: " << regionsSet << std::endl;

 

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  //std::string outputdir = "ZllGammaRatio_"+ cfg.mcSamples + regionsSet ;
  //  std::string outputdir = "Zll_" + configFileName;
  double intpart;
  double fracpart = modf(lumi, &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("_%.0ffb", intpart ) );
  outputdir += suffix;
  outputdir_of += suffix;
  
  system(Form("mkdir -p %s", outputdir.c_str()));








  std::string ZllDir = "ZllPurity_" + samples + "_" + regionsSet;
  std::string ZllDir_of = "ZllPurity_OF_" + samples + "_" + regionsSet;


  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees.root", ZllDir.c_str(), lumi, suffix.c_str()), "DYJets");

  //  MT2Analysis<MT2EstimateTree>* Zll = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb/Zll_analyses.root", ZllDir.c_str(), lumi), "DYJets");
  if( Zll==0 ) {
    std::cout << "-> Please run zllPurityTrees first. I need to get the yields from there." << std::endl;    std::cout << "-> Thank you for your cooperation." << std::endl;    exit(197);
  }



  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees.root", ZllDir.c_str(),lumi, suffix.c_str() ), "QCD");
  MT2Analysis<MT2EstimateTree>* top = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees.root", ZllDir.c_str(), lumi, suffix.c_str()), "Top");
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees.root", ZllDir.c_str(), lumi, suffix.c_str()), "WJets");
  MT2Analysis<MT2EstimateTree>* zjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees.root", ZllDir.c_str(), lumi, suffix.c_str()), "ZJets");


  //OPPOSITE FLAVOR TREES
  MT2Analysis<MT2EstimateTree>* Zll_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees_of.root", ZllDir_of.c_str(), lumi, suffix.c_str()), "DYJets");

  MT2Analysis<MT2EstimateTree>* qcd_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees_of.root", ZllDir_of.c_str(),lumi, suffix.c_str() ), "QCD");
  MT2Analysis<MT2EstimateTree>* top_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees_of.root", ZllDir_of.c_str(), lumi, suffix.c_str()), "Top");
  MT2Analysis<MT2EstimateTree>* wjets_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees_of.root", ZllDir_of.c_str(), lumi, suffix.c_str()), "WJets");
  MT2Analysis<MT2EstimateTree>* zjets_of = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb%s/ZllPurityTrees_of.root", ZllDir_of.c_str(), lumi, suffix.c_str()), "ZJets");



  /*
  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb/ZllPurityTrees.root", ZllDir.c_str(), lumi), "QCD");
  MT2Analysis<MT2EstimateTree>* top = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb/ZllPurityTrees.root", ZllDir.c_str(), lumi), "Top");
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb/ZllPurityTrees.root", ZllDir.c_str(), lumi), "WJets");
  MT2Analysis<MT2EstimateTree>* zjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s_%0.ffb/ZllPurityTrees.root", ZllDir.c_str(), lumi), "ZJets");
  */

  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
  bgYields.push_back( Zll );
  bgYields.push_back( qcd );
  bgYields.push_back( wjets );
  bgYields.push_back( zjets );
  bgYields.push_back( top );

  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_of; 
  bgYields_of.push_back( Zll_of );
  bgYields_of.push_back( qcd_of );
  bgYields_of.push_back( wjets_of );
  bgYields_of.push_back( zjets_of );
  bgYields_of.push_back( top_of );

  drawMll( outputdir, bgYields );
  drawMll( outputdir_of, bgYields_of );


 

  return 0;
}



























void drawSFvsOF( const std::string& outputdir, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields_of ) {

  MT2DrawTools::setStyle();

  std::vector<int> colors;
  if( bgYields.size()==3 ) { // estimates
    colors.push_back(402); 
    colors.push_back(430); 
    colors.push_back(418); 
  } else { // mc
    colors.push_back(430); // other = zll 
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    colors.push_back(419); // z+jets
    colors.push_back(855); // top
  }

  TH1F::AddDirectory(kTRUE);

  std::string fullPath = outputdir;
  std::string fullPathPlots = outputdir + "/plots";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  std::set<MT2Region> MT2Regions = bgYields[0]->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    //the loop over the B tagged regions
    //  for(int b= 0; b< 3; b++){

    THStack bgStack("bgStack", "");
    THStack bgStack_of("bgStack_of", "");
    for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
      int index = bgYields.size() - i - 1;
      TH1D* h1_bg = new TH1D("h1_bg","", 50, 0,250);
      TTree *bgTree = bgYields[index]->get(*iMT2)->tree;
      bgTree->Project("h1_bg","Z_mass","weight");
      h1_bg->SetFillColor( colors[index] );
      h1_bg->SetLineColor( kBlack );
      bgStack.Add(h1_bg);
    }



  }//end of regions


}//end of function






void drawMll( const std::string& outputdir, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields ) {

  MT2DrawTools::setStyle();

  std::vector<int> colors;
  if( bgYields.size()==3 ) { // estimates
    colors.push_back(402); 
    colors.push_back(430); 
    colors.push_back(418); 
  } else { // mc
    colors.push_back(430); // other = zll 
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    colors.push_back(419); // z+jets
    colors.push_back(855); // top
  }

  TH1F::AddDirectory(kTRUE);


  std::string fullPath = outputdir;
  std::string fullPathPlots = outputdir + "/plots";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  std::set<MT2Region> MT2Regions = bgYields[0]->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    //the loop over the B tagged regions
    //  for(int b= 0; b< 3; b++){

    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
      int index = bgYields.size() - i - 1;
      TH1D* h1_bg = new TH1D("h1_bg","", 50, 0,250);
      TTree *bgTree = bgYields[index]->get(*iMT2)->tree;
      //  bgTree->Project("h1_bg","Z_mass",Form("weight*(nBJets==%d)",b));
      bgTree->Project("h1_bg","Z_mass","weight");
      h1_bg->SetFillColor( colors[index] );
      h1_bg->SetLineColor( kBlack );
      bgStack.Add(h1_bg);
    }


    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();

    float yMax = 1.1*(bgStack.GetMaximum());
    //if(yMax < 0.3) 
yMax = 0.3; 

    TH2D* h2_axes = new TH2D("axes", "", 10, 0,250, 10, 0., yMax );
    h2_axes->SetXTitle("M_{ll} [GeV]");
    h2_axes->SetYTitle("Entries");
    h2_axes->Draw();
   
    std::vector<std::string> niceNames = thisRegion.getNiceNames();

    for( unsigned i=0; i<niceNames.size(); ++i ) {
      //   for( unsigned i=0; i<niceNames.size(); ++i ) {

      float yMaxText = 0.9-(float)i*0.05;
      float yMinText = yMaxText - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMinText, 0.55, yMaxText, "brNDC" );
      regionText->SetTextSize(0.035);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      regionText->AddText( niceNames[i].c_str() );
      /*  if( i==1 ){
	  regionText->AddText(Form("N(j) #geq 2; N(b) = %d", b) );
	  }else{
	  regionText->AddText( niceNames[i].c_str() );
	  }*/
      regionText->Draw("same");
    }
    

    TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    
    //   legend->AddEntry( h1_data, "Zll", "P" );
    //   legend->AddEntry( gr_data, "Zll", "P" );
    //   histoFile->cd();
    for( unsigned i=0; i<bgYields.size(); ++i ) {  
      TH1D* h1_bg1 = bgYields[i]->get(thisRegion)->yield;
      legend->AddEntry( h1_bg1, bgYields[i]->getFullName().c_str(), "F" );
      //     h1_bg->Write();
    }
    
    //   histoFile->Close();

    legend->Draw("same");
    bgStack.Draw("histo same");
    // h1_data->Draw("p same");
    // gr_data->Draw("p same");


    TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
    labelTop->Draw("same");
    gPad->RedrawAxis();

    /*
      c1->SaveAs(Form("%s/mt2_b%d_%s.eps",fullPathPlots.c_str(),b,thisRegion.getName().c_str()) );
      c1->SaveAs(Form("%s/mt2_b%d_%s.png",fullPathPlots.c_str(),b,thisRegion.getName().c_str()) );
      c1->SaveAs(Form("%s/mt2_b%d_%s.pdf",fullPathPlots.c_str(),b,thisRegion.getName().c_str()) );
    */
   
    c1->SaveAs( Form("%s/mll_%s.eps", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/mll_%s.png", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/mll_%s.pdf", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
   



    float bins_nJets[] = {2,4,7,12};
    float bins_nBJets[] = {0,1,2,3,6};
    //in MT2
    float bins_mt2[] = {200,300,400,500, 600, 800, 1000, 1500 , 1900};
    //in HT
    float bins_ht[] =  {450,575,1000,1500,2000};

    std::string cut =  "weight*(abs(Z_mass-91.19)<20)";
    std::string cut2 = "weight*(abs(Z_mass-91.19)<10)";
 
    std::string cut3 = "weight*(abs(Z_mass-91.19)<20&&nBJets<2)";
    std::string cut4 ="weight*(abs(Z_mass-91.19)<10&&nBJets<2)";

    std::string cut_nJets3 = "weight*(abs(Z_mass-91.19)<10&&nJets>2)";
    
    drawStacks( fullPathPlots,  bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields , thisRegion, cut );
    drawStacks( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , thisRegion , cut);
    drawStacks( fullPathPlots,  bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "zll_mt2", bgYields , thisRegion  , cut  );
    drawStacks( fullPathPlots,  bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "zll_ht", bgYields , thisRegion, cut );
    
    drawStacks( fullPathPlots,  bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields , thisRegion, cut2 );
    drawStacks( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , thisRegion , cut2);
    drawStacks( fullPathPlots,  bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "zll_mt2", bgYields , thisRegion  , cut2  );
    drawStacks( fullPathPlots,  bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "zll_ht", bgYields , thisRegion, cut2 );



    drawStacks( fullPathPlots,  bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields , thisRegion, cut4 );
    drawStacks( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , thisRegion , cut4);
    drawStacks( fullPathPlots,  bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "zll_mt2", bgYields , thisRegion  , cut4  );
    drawStacks( fullPathPlots,  bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "zll_ht", bgYields , thisRegion, cut4 );
    

    drawPurity( fullPathPlots,  bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields , thisRegion, cut );
     drawPurity( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , thisRegion , cut);
    drawPurity( fullPathPlots,  bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "zll_mt2", bgYields , thisRegion  , cut  );
    drawPurity( fullPathPlots,  bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "zll_ht", bgYields , thisRegion, cut );


  //NJETSCUT
   drawStacks( fullPathPlots,  bins_nBJets, sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , thisRegion , cut_nJets3);
   drawPurity( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , thisRegion , cut_nJets3);


    drawPurityComparison( fullPathPlots,  bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields , thisRegion, cut, cut2 );
    drawPurityComparison( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , thisRegion , cut, cut2);
    drawPurityComparison( fullPathPlots,  bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "zll_mt2", bgYields , thisRegion  , cut , cut2 );
    drawPurityComparison( fullPathPlots,  bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "zll_ht", bgYields , thisRegion, cut, cut2 );

    drawPurityComparison( fullPathPlots,  bins_nJets,sizeof(bins_nJets)/sizeof(float)-1,  "nJets", bgYields , thisRegion, cut3, cut4 );
    drawPurityComparison( fullPathPlots,  bins_nBJets,sizeof(bins_nBJets)/sizeof(float)-1,  "nBJets", bgYields , thisRegion , cut3, cut4);
    drawPurityComparison( fullPathPlots,  bins_mt2,sizeof(bins_mt2)/sizeof(float)-1,  "zll_mt2", bgYields , thisRegion  , cut3 , cut4 );
    drawPurityComparison( fullPathPlots,  bins_ht,sizeof(bins_ht)/sizeof(float)-1,  "zll_ht", bgYields , thisRegion, cut3, cut4 );




    delete c1;
    delete h2_axes;

    //   }//end of loop over beees


  }// for MT2 regions

}








void drawPurity(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, const MT2Region thisRegion, std::string cut){
 
  std::vector<int> colors;
  colors.push_back(430); // other = zll 
  colors.push_back(401); // qcd
  colors.push_back(417); // w+jets
  colors.push_back(419); // z+jets
  colors.push_back(855); // top

  TH1F::AddDirectory(kTRUE);

  float bins[size+1]; for(unsigned int i=0; i<= size ; i++)      bins[i]=binss[i];
  
  TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
  canny->cd();

  THStack Stack("Stack", "");
  THStack StackBG("StackBG", "");

  TH1D* purity  = new TH1D("purity", "", size, bins);
  TH1D* sNb  = new TH1D("sNb", "", size, bins);

  for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
    int index = bgYields.size() - i - 1;
    TH1D* h1_bg = new TH1D("h1_bg","", size  , bins);
    TTree *bgTree = bgYields[index]->get(thisRegion)->tree;
    bgTree->Project("h1_bg",Form("%s",name.c_str()) ,Form("%s",cut.c_str()) );
    h1_bg->SetFillColor( colors[index] );
    h1_bg->SetLineColor( kBlack );
    //  Stack.Add(h1_bg);
 
    sNb->Add(h1_bg);
    if(index ==0) continue;
    purity->Add(h1_bg);
  }

  purity->Sumw2(); sNb->Sumw2();
  purity->SetLineWidth(2);
  purity->Divide(sNb);
 
  float yMax = 1.1*(purity->GetMaximum());

  TH2D* h2_axes = new TH2D("axes", "", 10,bins[0] ,bins[size], 10, 0., 0.7 );
  // TH2D* h2_axes = new TH2D("axes", "", 10,bins[0] ,bins[size], 10, 0., yMax );
  if(name  == "zll_ht")  
   h2_axes->SetXTitle("H_{T}");
  else if(name == "zll_mt2")  
   h2_axes->SetXTitle("M_{T2}");
  else
    h2_axes->SetXTitle(name.c_str());
  h2_axes->SetYTitle("B/(S+B)");
  h2_axes->Draw();

 
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);

  labelTop->Draw("same");
  purity->Draw("same");

  gPad->RedrawAxis();

std::string extension= "";
  if(cut == "weight*(abs(Z_mass-91.19)<20)") extension ="mass";
  else if(cut == "weight*(abs(Z_mass-91.19)<10)") extension ="mass10";
  else if(cut == "weight*(abs(Z_mass-91.19)<10&&nJets>2)") extension ="nJets3";
  else  extension = "massNbJets";
  // "weight*(abs(Z_mass-91.19)<10&&nJets>2)";
    

 
  canny->SaveAs( Form("%s/%s_purity_%s_%s.eps", fullPath.c_str(), name.c_str() ,extension.c_str(), thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_purity_%s_%s.png", fullPath.c_str(), name.c_str() ,extension.c_str(), thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_purity_%s_%s.pdf", fullPath.c_str(), name.c_str() ,extension.c_str(), thisRegion.getName().c_str()) );

 
  delete h2_axes;
  delete canny;
}


















void drawPurityComparison(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, const MT2Region thisRegion, std::string cut , std::string cut2 ){
 
  std::vector<int> colors;
  colors.push_back(430); // other = zll 
  colors.push_back(401); // qcd
  colors.push_back(417); // w+jets
  colors.push_back(419); // z+jets
  colors.push_back(855); // top

  TH1F::AddDirectory(kTRUE);

  float bins[size+1]; for(unsigned int i=0; i<= size ; i++)      bins[i]=binss[i];
  
  TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
  canny->cd();

  THStack Stack("Stack", "");
  THStack StackBG("StackBG", "");

  TH1D* purity  = new TH1D("purity", "", size, bins);
  TH1D* sNb  = new TH1D("sNb", "", size, bins);

  THStack Stack2("Stack2", "");
  THStack StackBG2("StackBG2", "");

  TH1D* purity2  = new TH1D("purity2", "", size, bins);
  TH1D* sNb2  = new TH1D("sNb2", "", size, bins);

  for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
    int index = bgYields.size() - i - 1;
    TH1D* h1_bg = new TH1D("h1_bg","", size  , bins);
    TH1D* h1_bg2 = new TH1D("h1_bg2","", size  , bins);
    TTree *bgTree = bgYields[index]->get(thisRegion)->tree;
    bgTree->Project("h1_bg",Form("%s",name.c_str()) ,Form("%s",cut.c_str()) );
    bgTree->Project("h1_bg2",Form("%s",name.c_str()) ,Form("%s",cut2.c_str()) );
    sNb->Add(h1_bg);    sNb2->Add(h1_bg2);
    if(index ==0) continue;
    purity->Add(h1_bg);   purity2->Add(h1_bg2);
  }

  purity->Divide(sNb);  purity2->Divide(sNb2);
  purity->SetLineColor(kBlack); purity->SetLineWidth(2);
  purity2->SetLineColor(kRed); purity2->SetLineWidth(2);


  float yMax = 1.1*(purity->GetMaximum());

  TH2D* h2_axes = new TH2D("axes", "", 10,bins[0] ,bins[size], 10, 0., 0.4);
  if(name  == "zll_ht")  
    h2_axes->SetXTitle("H_{T}");
  else if(name == "zll_mt2")  
   h2_axes->SetXTitle("M_{T2}");
  else
    h2_axes->SetXTitle(name.c_str());
  h2_axes->SetYTitle("B/(S+B)");
  h2_axes->Draw();


  TLegend* legend = new TLegend( 0.2, 0.7-2*0.06, 0.4, 0.7 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( purity, "20GeV", "L" );
  legend->AddEntry( purity2, "10GeV", "L" );
  legend->Draw("same");

 
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);

  labelTop->Draw("same");
  purity->Draw("same");
  purity2->Draw("same");

  gPad->RedrawAxis();

  std::string extension= "";
  if(cut == "weight*(abs(Z_mass-91.19)<10)") extension ="mass";
  else extension = "massNbJets";


 

  

  canny->SaveAs( Form("%s/%s_purity2_%s_%s.eps", fullPath.c_str(), name.c_str() ,extension.c_str() , thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_purity2_%s_%s.png", fullPath.c_str(), name.c_str() , extension.c_str() ,thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_purity2_%s_%s.pdf", fullPath.c_str(), name.c_str() , extension.c_str() ,thisRegion.getName().c_str()) );

 
  delete h2_axes;
  delete canny;
}


























void drawStacks(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, const MT2Region thisRegion, std::string cut){
 
  std::vector<int> colors;
  colors.push_back(430); // other = zll 
  colors.push_back(401); // qcd
  colors.push_back(417); // w+jets
  colors.push_back(419); // z+jets
  colors.push_back(855); // top

  TH1F::AddDirectory(kTRUE);

  float bins[size+1]; for(unsigned int i=0; i<= size ; i++)      bins[i]=binss[i];
 
 
  TCanvas* canny = new TCanvas( "canny", "", 600, 600 );
  canny->cd();
 
  // std::string cut = "weight*( abs(Z_mass-91.19)<20)";

  THStack nJetsStack("nJetsStack", "");

  for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
    int index = bgYields.size() - i - 1;
    TH1D* h1_bg = new TH1D("h1_bg","", size  , bins);
    TTree *bgTree = bgYields[index]->get(thisRegion)->tree;
    bgTree->Project("h1_bg",Form("%s",name.c_str()) ,Form("%s",cut.c_str()) );
    h1_bg->SetFillColor( colors[index] );
    h1_bg->SetLineColor( kBlack );
    nJetsStack.Add(h1_bg);
  }
  float yMax_nJets = 1.1*(nJetsStack.GetMaximum());


  TH2D* h2_axes = new TH2D("axes", "", 10,bins[0] ,bins[size], 10, 0., yMax_nJets );
  if(name  == "zll_ht")  
   h2_axes->SetXTitle("H_{T}");
  else if(name == "zll_mt2")  
   h2_axes->SetXTitle("M_{T2}");
  else
    h2_axes->SetXTitle(name.c_str());
  h2_axes->SetYTitle("Entries");
  h2_axes->Draw();



  TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  for( unsigned i=0; i<bgYields.size(); ++i ) {  
    TH1D* h1_bg1 = bgYields[i]->get(thisRegion)->yield;
    legend->AddEntry( h1_bg1, bgYields[i]->getFullName().c_str(), "F" );
  }
 
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);

  labelTop->Draw("same");
  legend->Draw("same");
  nJetsStack.Draw("histo same");
  gPad->RedrawAxis();

  std::string extension= "";
  if(cut == "weight*(abs(Z_mass-91.19)<20)") extension ="mass";
  else if(cut == "weight*(abs(Z_mass-91.19)<10)") extension ="mass10";
  else if(cut == "weight*(abs(Z_mass-91.19)<10&&nJets>2)") extension ="nJets3";
  else  extension = "massNbJets";
  // "weight*(abs(Z_mass-91.19)<10&&nJets>2)";
    


 
  canny->SaveAs( Form("%s/%s_%s_%s.eps", fullPath.c_str(), name.c_str(),extension.c_str() , thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_%s_%s.png", fullPath.c_str(), name.c_str(),extension.c_str() , thisRegion.getName().c_str()) );
  canny->SaveAs( Form("%s/%s_%s_%s.pdf", fullPath.c_str(), name.c_str(),extension.c_str() , thisRegion.getName().c_str()) );


 
  delete h2_axes;
  delete canny;
}
