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

double lumi=4.; //fb-1
bool doNminusOne=false;

void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  std::string var, std::string label, std::string selection, bool logY );

int main( int argc, char* argv[] ) {

  MT2DrawTools::setStyle();

  if( argc!=2 ) {
    std::cout << "USAGE: ./scalesVariation [dir]" << std::endl;
    exit(113);
  }


  std::string dir( argv[1] );
  std::string mc_fileName = dir + "/analyses.root";

  std::string outputdir = dir + "/scalesVariations";

  MT2Analysis<MT2EstimateTree>* data  = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "data" );
  
  MT2Analysis<MT2EstimateTree>* top   = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "Top");
  
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields;
  bgYields.push_back(top);

  drawHisto( outputdir, data, bgYields,  "mt2", "M_{T2} [GeV]", "", kFALSE );
//  drawHisto( outputdir, data, bgYields, "ht", 100, 0., 2500., "H_{T} [GeV]", "", kTRUE );
//  drawHisto( outputdir, data, bgYields, "nJets", 12, 0, 12, "N(jet)", "", kFALSE );
//  drawHisto( outputdir, data, bgYields, "nBJets", 6, 0, 6, "N(b-tag)", "", kFALSE );

  return 0;

}// main



void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields, std::string var, std::string label, std::string selection, bool logY ) {

  MT2DrawTools::setStyle();

  std::vector<int> colors; //mc
  colors.push_back(855); // top
    
  std::string fullPath = outputdir;
  system( Form("mkdir -p %s", fullPath.c_str()) );

  std::set<MT2Region> MT2Regions = data->getRegions();

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
   
    TH1F::AddDirectory(kTRUE);

    
    int nBins;
    double* bins;
    thisRegion.getBins(nBins, bins);

      
    TH1D* h_data = new TH1D("h_data", "", nBins, bins);
    
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

    std::string scale_=selection;
    if( !doNminusOne ){
    
      if( selection == "" ) selection = Form("(mt2>200.)*weight*%0.f", lumi);
      else selection = Form("(mt2>200. &&  %s)*weight*%0.f", selection.c_str(), lumi);
    
    }
    else{
      
      selection = "";
      for( int c=0; c < nCuts; ++c ){

	if( cut[c] == MinusOneCut );
	else selection += cut[c];
      
      }
      
      selection = Form("(%s)*weight*%0.f", selection.c_str(), lumi);
      
    }
    
    data->get( thisRegion )->tree->Project( h_data->GetName(), var.c_str(), selection.c_str() );
  
    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.6);

    THStack bgStack("bgStack", "");
    TH1D* h_bg[bgYields.size()];
    
    //    int nVariations=9;
    int nVariations=101;
    int firstId=10;
    std::vector<std::string> variations;
//    variations.push_back("#mu_{R}=1, #mu_{F}=1");
//    variations.push_back("#mu_{R}=1, #mu_{F}=2");
//    variations.push_back("#mu_{R}=1, #mu_{F}=0.5");
//    variations.push_back("#mu_{R}=2, #mu_{F}=1");
//    variations.push_back("#mu_{R}=2, #mu_{F}=2");
//    variations.push_back("#mu_{R}=2, #mu_{F}=0.5");
//    variations.push_back("#mu_{R}=0.5, #mu_{F}=1");
//    variations.push_back("#mu_{R}=0.5, #mu_{F}=2");
//    variations.push_back("#mu_{R}=0.5, #mu_{F}=0.5");
    variations.push_back("NNPDF3.0 LO, #alpha_{s}=0.130");
//    variations.push_back("CTEQ6L1");
    
    TH1D* h_bg_v[bgYields.size()][nVariations];
    TH1D* h_bg_r[bgYields.size()][nVariations];

    float normalization;
    float integral0;
    float sigma2;
    
    for( unsigned i=0; i<bgYields.size(); ++i ) {
    //for( unsigned i=1; i<bgYields.size(); ++i ) { //Not drawing QCD

      sigma2=0.0;
      integral0=0.0;
      normalization=1.0;
      int index=i;
      
      //int index = bgYields.size() - i - 1; // reverse ordered stack is prettier with QCD
      
      h_bg[index] = new TH1D(Form("h_%s_%s", var.c_str(), bgYields[i]->getName().c_str()), "", nBins, bins);      
      h_bg[index]->Sumw2();

      bgYields[index]->get(thisRegion)->tree->Project( h_bg[index]->GetName(), var.c_str(), selection.c_str() );

      std::vector<int> colors_v;
//      colors_v.push_back(2);
//      colors_v.push_back(3);
//      colors_v.push_back(4);
//      colors_v.push_back(5);
//      colors_v.push_back(6);
//      colors_v.push_back(7);
//      colors_v.push_back(8);
//      colors_v.push_back(9);
//      colors_v.push_back(30);
      
      for(int v=0; v < nVariations; ++v){
	
	colors_v.push_back(2);

	std::string scale = scale_;
	normalization = h_bg[index]->Integral();

	h_bg_v[index][v] = new TH1D(Form("h_%s_%s_v%d", var.c_str(), bgYields[i]->getName().c_str(), v+firstId), "", nBins, bins);
	h_bg_v[index][v]->Sumw2();
	
	if( scale == "" ) scale = Form("(mt2>200. && LHEweight_id==%d)*LHEweight_wgt/LHE_originalWeight*weight*%0.f", v+firstId, lumi);
	else scale = Form("(mt2>200. &&  %s && LHEweight_id==%d)*LHEweight_wgt/LHE_originalWeight*weight*%0.f", scale.c_str(), v+firstId, lumi);
//	if( scale == "" ) scale = Form("(mt2>200. && LHEweight_id==%d)*LHEweight_wgt/LHEweight_wgt[%d]*weight*%0.f", v+firstId, firstId, lumi);
//	else scale = Form("(mt2>200. &&  %s && LHEweight_id==%d)*LHEweight_wgt/LHEweight_wgt[%d]*weight*%0.f", scale.c_str(), v+firstId, firstId, lumi);
	bgYields[index]->get(thisRegion)->tree->Project( h_bg_v[index][v]->GetName(), var.c_str(), scale.c_str() );
	
	float thisIntegral = h_bg_v[index][v]->Integral();
	normalization = 1.*normalization/thisIntegral;
	
	integral0 = h_bg_v[index][0]->Integral();
	sigma2 += (thisIntegral-integral0)*(thisIntegral-integral0);
	
	h_bg_v[index][v]->Scale(normalization);
	
	h_bg_v[index][v]->SetLineColor(colors_v[v]);
	h_bg_v[index][v]->SetLineWidth(2);
	h_bg_v[index][v]->SetMarkerColor(colors_v[v]);
	//h_bg_v[index][v]->SetFillColor(0);

	h_bg_r[index][v]=(TH1D*) h_bg_v[index][v]->Clone(Form("h_%s_%s_r%d", var.c_str(), bgYields[i]->getName().c_str(), v+firstId));
	h_bg_r[index][v]->Sumw2();

	h_bg_r[index][v]->Divide(h_bg[index]);

	for(int b=1; b<=nBins; ++b)
	  h_bg_r[index][v]->SetBinError(b, 1e-3);

      }

      h_bg[index]->SetFillColor( colors[index] );
      h_bg[index]->SetLineColor( kBlack );
      bgStack.Add(h_bg[index]);
      
    }
    
  
    float sigma=sqrt(1./(nVariations-1-1)*sigma2);
    std::cout << "sigma = " << sigma << ", central value = " << integral0 << std::endl;

    TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
    c1->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
    pad1->SetBottomMargin(0.15);
    pad1->Draw();
    pad1->cd();

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
    TH2D* h_axes = new TH2D("axes", "", 10, 200, 1500, 10, yMin, yMax );
    h_axes->GetXaxis()->SetTitle(label.c_str());
    h_axes->GetYaxis()->SetTitle(labelY.c_str());
    h_axes->GetYaxis()->SetTitleOffset(1.5);
    h_axes->GetYaxis()->SetLabelSize(0.04);
    h_axes->Draw();

    std::vector<std::string> niceNames = thisRegion.getNiceNames();
    //  niceNames.push_back("M_{T2} > 200 GeV");
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
 
    //    TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+nVariations+1)*0.04, 0.93, 0.9 );
    TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+1)*0.04, 0.93, 0.9 );
    //legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    
    //legend->AddEntry( gr_data, "Data", "P" ); //Not drawing data yet
  
    for( unsigned i=0; i<bgYields.size(); ++i ) { //Not drawing QCD
      //for( unsigned i=1; i<bgYields.size(); ++i ) {
      legend->AddEntry( h_bg[i], bgYields[i]->getFullName().c_str(), "F" );
    }
//    for( int v=0; v<nVariations; ++v)
//      legend->AddEntry( h_bg_v[0][v], variations[v].c_str(), "L");
    for( int v=0; v<1; ++v)
      legend->AddEntry( h_bg_v[0][v], variations[v].c_str(), "L");

    legend->Draw("same");
    bgStack.Draw("histo same");
    bgStack.SetMinimum(yMin);
    bgStack.Draw("histo same");
    //gr_data->Draw("p same"); //Not drawing data yet

    for( unsigned i=0; i<bgYields.size(); ++i )
      for(int v=0; v<nVariations; ++v)
	h_bg_v[i][v]->Draw("histo,same");
 
    //TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
    TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi); //Not drawing data yet
    labelTop->Draw("same");

    gPad->RedrawAxis();

    c1->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.1);
    pad2->Draw();
    pad2->cd();
    
    pad2->SetGridy();
    
    for( unsigned i=0; i<bgYields.size(); ++i )
      for(int v=0; v<nVariations; ++v){

	h_bg_r[i][v]->SetStats(0);
	//h_bg_r[i][v]->SetMarkerStyle(21);
	//h_bg_r[i][v]->SetMarkerSize(0.5);
	h_bg_r[i][v]->GetXaxis()->SetLabelSize(0.00);
	h_bg_r[i][v]->GetXaxis()->SetTickLength(0.09);
	h_bg_r[i][v]->GetYaxis()->SetNdivisions(5,5,0);
	//	h_bg_r[i][v]->GetYaxis()->SetRangeUser(0.8,1.2);
	h_bg_r[i][v]->GetYaxis()->SetRangeUser(0.,2.);
	h_bg_r[i][v]->GetYaxis()->SetTitleSize(0.17);
	h_bg_r[i][v]->GetYaxis()->SetTitleOffset(0.4);
	h_bg_r[i][v]->GetYaxis()->SetLabelSize(0.17);
	h_bg_r[i][v]->GetYaxis()->SetTitle("ratio");
       	
      }
    
    h_bg_r[0][0]->Draw("EP");
    
    for( unsigned i=0; i<bgYields.size(); ++i )
      for(int v=0; v<nVariations; ++v)
	h_bg_r[i][v]->Draw("EP,same");

    gPad->RedrawAxis();
      
    c1->cd();
    
    c1->SaveAs( Form("%s/%s_%s.eps", fullPath.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.png", fullPath.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.pdf", fullPath.c_str(), var.c_str(), thisRegion.getName().c_str()) );

    delete c1;
    delete h_axes;
    
    delete h_data;
    
    delete bins;
    
    for(unsigned i=0; i<bgYields.size(); ++i){
      //for(unsigned i=1; i<bgYields.size(); ++i){ //Not drawing QCD
    
      delete h_bg[i];
      
      for(int v=0; v<nVariations; ++v){
	delete h_bg_v[i][v];
	delete h_bg_r[i][v];
      }
    
    }
    
    TH1F::AddDirectory(kFALSE);

  }// for MT2 regions

}
