#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"

#include "TRandom3.h"

#define mt2_cxx
#include "interface/mt2.h"

double lumi=5.; //fb-1
bool doNminusOne=false;

void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, std::vector<MT2Analysis<MT2EstimateTree>* > sigYields, std::string var, int nBins, float xmin, float xmax, std::string label, std::string selection, bool logY );

int main( int argc, char* argv[] ) {

  MT2DrawTools::setStyle();

  if( argc!=2 ) {
    std::cout << "USAGE: ./inclusivePlots [dir]" << std::endl;
    exit(113);
  }


  std::string dir( argv[1] );
  std::string mc_fileName = dir + "/analyses.root";

  std::string outputdir = dir + "/inclusivePlots";

  MT2Analysis<MT2EstimateTree>* data  = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "data" );
  
  MT2Analysis<MT2EstimateTree>* qcd   = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "QCD" );
  MT2Analysis<MT2EstimateTree>* zjets = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "ZJets" );
  zjets->setFullName("Z+jets");
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "WJets");
  wjets->setFullName("W+jets");
  MT2Analysis<MT2EstimateTree>* top   = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "Top");
  
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields;
  bgYields.push_back(qcd);
  bgYields.push_back(wjets);
  bgYields.push_back(zjets);
  bgYields.push_back(top);

  //  std::vector<MT2Analysis<MT2EstimateTree>*> sigYields = MT2Analysis<MT2EstimateTree>::readAllFromFile( mc_fileName, "SMS" );
  std::vector<MT2Analysis<MT2EstimateTree>*> sigYields = MT2Analysis<MT2EstimateTree>::readAllFromFile( mc_fileName, "DarkMatter" );

  drawHisto( outputdir, data, bgYields, sigYields, "mt2", 60, 0., 1500., "M_{T2} [GeV]", "", kTRUE );
  drawHisto( outputdir, data, bgYields, sigYields, "ht", 100, 0., 2500., "H_{T} [GeV]", "", kTRUE );
  drawHisto( outputdir, data, bgYields, sigYields, "nJets", 12, 0, 12, "N(jet)", "", kFALSE );
  drawHisto( outputdir, data, bgYields, sigYields, "nBJets", 6, 0, 6, "N(b-tag)", "", kFALSE );

  return 0;

}// main



void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields, std::vector<MT2Analysis<MT2EstimateTree>* > sigYields, std::string var, int nBins, float xmin, float xmax, std::string label, std::string selection, bool logY ) {

  MT2DrawTools::setStyle();

  std::vector<int> colors; //mc
  colors.push_back(401); // qcd
  colors.push_back(417); // w+jets
  colors.push_back(419); // z+jets
  colors.push_back(855); // top
  //colors.push_back(); // other
  
  std::vector<int> colorsSig; //PHYS14 mc
//  colorsSig.push_back(1); // t1tttt
//  colorsSig.push_back(2); // t1bbbb
//  colorsSig.push_back(6); // t1qqqq
//  colorsSig.push_back(5); // t2tt
//  colorsSig.push_back(7); // t2bb
//  colorsSig.push_back(9); // t2qq
  //////DM
  colorsSig.push_back(1); // 1000
  colorsSig.push_back(1);
  colorsSig.push_back(2); // 100
  colorsSig.push_back(2);
  colorsSig.push_back(6); // 10
  colorsSig.push_back(6); 
  colorsSig.push_back(9); // 1
  
  std::vector<int> styleSig;
  styleSig.push_back(1);
  styleSig.push_back(2);
  styleSig.push_back(1);
  styleSig.push_back(2);
  styleSig.push_back(1);
  styleSig.push_back(2);
  styleSig.push_back(1);

 
  std::string fullPath = outputdir;
  system( Form("mkdir -p %s", fullPath.c_str()) );

  std::set<MT2Region> MT2Regions = data->getRegions();

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
   
    TH1F::AddDirectory(kTRUE);
    
    TH1D* h_data = new TH1D("h_data", "", nBins, xmin, xmax);
   
    // All selections:
    std::string MinusOneCut=selection;
    int nCuts = 8;
    std::string cut[nCuts];
    cut[0] = "((ht>1000. && met>30.) || (ht>450. && met>200.))";
    cut[1] = "nVert>0";
    cut[2] = "nJets>= 2";
    cut[3] = "nElectrons10==0 && nMuons10==0";
    cut[4] = "nPFLep5LowMT==0 && nPFHad10LowMT==0";
    cut[5] = "dPhiMin>0.3";
    cut[6] = "diffMetMht/met<0.5";
    cut[7] = "mt2>200.";

    if( !doNminusOne ){
    
      if( selection == "" ) selection = Form("(mt2>200.)*weight");
      else selection = Form("(mt2>200. &&  %s)*weight", selection.c_str());
    
    }
    else{
      
      selection = "";
      for( int c=0; c < nCuts; ++c ){

	if( cut[c] == MinusOneCut );
	else selection += cut[c];
      
      }
      
      selection = Form("(%s)*weight", selection.c_str());
      
    }
    
    data->get( thisRegion )->tree->Project( h_data->GetName(), var.c_str(), selection.c_str() );
  
    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.6);

    THStack bgStack("bgStack", "");
    TH1D* h_bg[bgYields.size()];
    for( unsigned i=0; i<bgYields.size(); ++i ) {
    //for( unsigned i=1; i<bgYields.size(); ++i ) { //Not drawing QCD

      //int index=i;
      int index = bgYields.size() - i - 1; // reverse ordered stack is prettier with QCD
      
      h_bg[index] = new TH1D(Form("h_%s_%s", var.c_str(), bgYields[i]->getName().c_str()), "", nBins, xmin, xmax);
      bgYields[index]->get(thisRegion)->tree->Project( h_bg[index]->GetName(), var.c_str(), selection.c_str() );
      h_bg[index]->SetFillColor( colors[index] );
      h_bg[index]->SetLineColor( kBlack );
      bgStack.Add(h_bg[index]);
      
    }
    
    TH1D* h_sig[sigYields.size()];    
    for( unsigned i=0; i<sigYields.size(); ++i ) { 
      
      h_sig[i] = new TH1D(Form("h_%s_%s", var.c_str(), sigYields[i]->getName().c_str()), "", nBins, xmin, xmax);
      sigYields[i]->get(thisRegion)->tree->Project( h_sig[i]->GetName(), var.c_str(), selection.c_str() );
      h_sig[i]->SetLineColor( colorsSig[i] );
      //DM
      h_sig[i]->SetLineStyle( styleSig[i] );
      //
      h_sig[i]->SetLineWidth( 2 );
      //h_sig[i]->Scale(50.);

    }
    
    std::vector< std::string > sigNames;
//    sigNames.push_back("T1tttt 1500, 100 x50");
//    sigNames.push_back("T1bbbb 1500, 100 x50");
//    sigNames.push_back("T1qqqq 1400, 100 x50");
//    sigNames.push_back("T2tt 850, 100 x50");
//    sigNames.push_back("T2bb 900, 100 x50");
//    sigNames.push_back("T2qq 1200, 100 x50");
    for (int i=0; i<sigYields.size(); ++i){
      TString thisName=sigYields[i]->getName().c_str();
      thisName.ReplaceAll( "M1", "M=1" );
      thisName.ReplaceAll( "tsg", "" );
      thisName.ReplaceAll( "V", "(V)" );
      thisName.ReplaceAll( "A(V)", "(AV)" );	    
      thisName.ReplaceAll( "_", " " );

      sigNames.push_back((std::string)thisName);
    }

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();

    float yMin = 0.;
    float yMax1 = h_data->GetMaximum()*1.5;
    float yMax2 = 1.5*(h_data->GetMaximum() + sqrt(h_data->GetMaximum()));
    float yMax3 = 1.5*(bgStack.GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    if( h_data->GetNbinsX()<2 ) yMax *=3.;
    
    if(logY) {
      gPad->SetLogy();
      yMin=1e-1;
      yMax*=50.;
    }
    std::string labelY = "Events";
    if( (xmax - xmin)/nBins != 1. ) labelY = (Form("Events/%.0f GeV", (xmax - xmin)/nBins));
    TH2D* h_axes = new TH2D("axes", "", 10, xmin, xmax, 10, yMin, yMax );
    h_axes->GetXaxis()->SetTitle(label.c_str());
    h_axes->GetYaxis()->SetTitle(labelY.c_str());
    h_axes->GetYaxis()->SetTitleOffset(1.5);
    h_axes->GetYaxis()->SetLabelSize(0.04);
    h_axes->Draw();

    std::vector<std::string> niceNames = thisRegion.getNiceNames();
    //niceNames.push_back("M_{T2} > 200 GeV");
    for( unsigned i=0; i<niceNames.size(); ++i ) {

      float yMax = 0.9-(float)i*0.04;
      float yMin = yMax - 0.04;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      //TPaveText* regionText = new TPaveText( 0.18, yMin, 0.5, yMax, "brNDC" );
      //regionText->SetTextSize(0.035);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      regionText->AddText( niceNames[i].c_str() );
      regionText->Draw("same");

    }
 
    TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+sigYields.size()+1)*0.03, 0.93, 0.9 );
    //TLegend* legend = new TLegend( 0.6, 0.9-(sigYields.size()+1)*0.04, 0.93, 0.9 ); //////Signal only
    //TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+sigYields.size()-1)*0.03, 0.93, 0.9 ); //Not drawing QCD and data --> +1 -> -1
    //legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    
    //legend->AddEntry( gr_data, "Data", "P" ); //Not drawing data yet

    //////Signal only  
    for( unsigned i=0; i<bgYields.size(); ++i ) { //Not drawing QCD
      //for( unsigned i=1; i<bgYields.size(); ++i ) {
      legend->AddEntry( h_bg[i], bgYields[i]->getFullName().c_str(), "F" );
    }
    //////
    
    for( unsigned i=0; i<sigNames.size(); ++i ) {
      legend->AddEntry( h_sig[i], sigNames[i].c_str(), "l" );
    }

    legend->Draw("same");
    ////Only Signal
    bgStack.Draw("histo same");
    bgStack.SetMinimum(yMin);
    bgStack.Draw("histo same");
    ////
    
    //gr_data->Draw("p same"); //Not drawing data yet

    for( unsigned i=0; i<sigYields.size(); ++i ) {
      h_sig[i]->Draw("same");
    }
 
    //TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
    TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi); //Not drawing data yet
    labelTop->Draw("same");

    gPad->RedrawAxis();

    c1->SaveAs( Form("%s/%s_%s.eps", fullPath.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.png", fullPath.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.pdf", fullPath.c_str(), var.c_str(), thisRegion.getName().c_str()) );

    delete c1;
    delete h_axes;
    
    delete h_data;

    for(unsigned i=0; i<sigYields.size(); ++i){

      delete h_sig[i];
    
    }
    
    for(unsigned i=0; i<bgYields.size(); ++i){
      //for(unsigned i=1; i<bgYields.size(); ++i){ //Not drawing QCD
    
      delete h_bg[i];
    
    }
    
    TH1F::AddDirectory(kFALSE);

  }// for MT2 regions

}
