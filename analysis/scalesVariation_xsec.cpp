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
#include "interface/MT2Config.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"

#include "TRandom3.h"

#define mt2_cxx
#include "interface/mt2.h"

float lumi;

void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  std::string var, std::string label, std::string selection, bool logY );

int main( int argc, char* argv[] ) {



  MT2DrawTools::setStyle();

  if( argc<2 ) {
    std::cout << "USAGE: ./scalesVariation_xsec [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);
  
  lumi=cfg.lumi();
  
  std::string outputdir = cfg.getEventYieldDir() + "/scaleVariations_xsec05";
  system(Form("mkdir -p %s", outputdir.c_str()));


  std::string mcFile = cfg.getEventYieldDir() + "/analyses.root";
  std::string dataFile = cfg.getEventYieldDir() + "/analyses.root";


  MT2Analysis<MT2EstimateTree>* data  = MT2Analysis<MT2EstimateTree>::readFromFile( dataFile, "data" );
  
  MT2Analysis<MT2EstimateTree>* top   = MT2Analysis<MT2EstimateTree>::readFromFile( mcFile, "Top");
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile( mcFile, "WJets");
  
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields;
  bgYields.push_back(top);
  bgYields.push_back(wjets);

  drawHisto( outputdir, data, bgYields,  "mt2", "M_{T2} [GeV]", "", kTRUE );

  return 0;

}// main



void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields, std::string var, std::string label, std::string selection, bool logY ) {

  MT2DrawTools::setStyle();

  std::vector<int> colors; //mc
  colors.push_back(855); // top
  colors.push_back(417); // w+jets
    
  std::set<MT2Region> MT2Regions = data->getRegions();

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
   
    TH1F::AddDirectory(kTRUE);

    int nBins;
    double* bins;
    thisRegion.getBins(nBins, bins);

      
    TH1D* h_data = new TH1D("h_data", "", nBins, bins);

    std::string selection_ = Form("(mt2>200 && ht>200 && deltaPhiMin>0.3 && diffMetMht/met<0.5)*weight");	  
    data->get( thisRegion )->tree->Project( h_data->GetName(), var.c_str(), selection_.c_str() );
  
    THStack bgStack("bgStack", "");
    TH1D* h_bg[bgYields.size()];

    float integral[3];
    integral[2]=0;

    float relativeF[2];
    float relativeV[2];

    for( unsigned i=0; i<bgYields.size(); ++i ) {

      int index=i;

      h_bg[index] = new TH1D(Form("h_%s_%s", var.c_str(), bgYields[i]->getName().c_str()), "", nBins, bins);
      h_bg[index]->Sumw2();

      selection_ = Form("(mt2>200 && ht>200 && deltaPhiMin>0.3 && diffMetMht/met<0.5)*weight");
      bgYields[index]->get(thisRegion)->tree->Project( h_bg[index]->GetName(), var.c_str(), selection_.c_str() );

      integral[index]=h_bg[index]->Integral();
      integral[2]+=integral[index];

    }
    
    relativeF[0]=integral[0]/integral[2];
    relativeF[1]=integral[1]/integral[2];

    std::cout << relativeF[0] << "\t" << relativeF[1] << std::endl;
    
    if( relativeF[0]>0 && relativeF[1]>0 ){

      if( relativeF[0] >= relativeF[1] ){
	
	relativeV[1]=0.5;
	relativeV[0]=-1+(1-relativeF[1]*relativeV[1]/relativeF[0]);
	
      }
      else{
	
	relativeV[0]=0.5;
	relativeV[1]=-1+(1-relativeF[0]*relativeV[0]/relativeF[1]);
	
      }
    }
    else{
      
      relativeV[0]=0.0;
      relativeV[1]=0.0;

    }

    std::cout << relativeV[0] << "\t" << relativeV[1] << std::endl;
    
    TH1D* h_bg_sum;

    int nVariations=2;   
    TH1D* h_bg_v[bgYields.size()][nVariations];
    TH1D* h_bg_v_sum[nVariations];
    TH1D* h_bg_r[nVariations];

    for( unsigned i=0; i<bgYields.size(); ++i ) {

      int index=i;
      
      std::vector<int> colors_v;
      colors_v.push_back(1);
      colors_v.push_back(2);
      
//      for(int v=0; v < nVariations; ++v){
//	
//	h_bg_v[index][v] = new TH1D(Form("h_%s_%s_v%d", var.c_str(), bgYields[i]->getName().c_str(), v), "", nBins, bins);
//	h_bg_v[index][v]->Sumw2();
//	
//	selection_ = Form("(mt2>200 && ht>200 && deltaPhiMin>0.3 && diffMetMht/met<0.5)*weight*(%f)", 1.+relativeV[index]);
//	std::cout << index << "\t" << selection_ << std::endl;
//
//	bgYields[index]->get(thisRegion)->tree->Project( h_bg_v[index][v]->GetName(), var.c_str(), selection_.c_str() );
//
//	h_bg_v[index][v+1] = new TH1D(Form("h_%s_%s_v%d", var.c_str(), bgYields[i]->getName().c_str(), v+1), "", nBins, bins);
//	h_bg_v[index][v+1]->Sumw2();
//	
//	selection_ = Form("(mt2>200 && ht>200 && deltaPhiMin>0.3 && diffMetMht/met<0.5)*weight*(%f)", 1.-relativeV[index]);
//	std::cout << index << "\t" << selection_ << std::endl;
//
//	bgYields[index]->get(thisRegion)->tree->Project( h_bg_v[index][v+1]->GetName(), var.c_str(), selection_.c_str() );
//	
//	v+=1;
//	
//      }

      for(int v=0; v < nVariations; ++v){
	
	h_bg_v[index][v] = (TH1D*) h_bg[index]->Clone((Form("h_%s_%s_v%d", var.c_str(), bgYields[i]->getName().c_str(), v)));
	h_bg_v[index][v]->Sumw2();
	h_bg_v[index][v]->Scale(1.+relativeV[index]);

	std::cout << index << "\t" << 1.+relativeV[index] << std::endl;

	h_bg_v[index][v+1] = (TH1D*) h_bg[index]->Clone((Form("h_%s_%s_v%d", var.c_str(), bgYields[i]->getName().c_str(), v+1)));
	h_bg_v[index][v+1]->Sumw2();
	h_bg_v[index][v+1]->Scale(1.-relativeV[index]);

	std::cout << index << "\t" << 1.-relativeV[index] << std::endl;
						      	
	v+=1;
	
      }


      h_bg[index]->SetFillColor( colors[index] );
      h_bg[index]->SetLineColor( colors[index] );
      bgStack.Add(h_bg[index]);

      if(index==0){
	h_bg_sum = (TH1D*) h_bg[index]->Clone("h_bg_sum");
	h_bg_sum->Sumw2();
      }
      else
	h_bg_sum->Add(h_bg[index]);

    }

    float normalization = 1./(h_bg_sum->Integral());
    std::cout << 1./normalization << std::endl;
    std::cout << h_bg[0]->Integral() << std::endl;
    std::cout << h_bg[1]->Integral() << std::endl;

    h_bg_v_sum[0] = (TH1D*) h_bg_v[0][0]->Clone("h_bg_v_sum_up");
    h_bg_v_sum[1] = (TH1D*) h_bg_v[0][1]->Clone("h_bg_v_sum_dn");

    h_bg_v_sum[0]->Sumw2();
    h_bg_v_sum[1]->Sumw2();
    
    h_bg_v_sum[0]->Add(h_bg_v[1][0]);
    h_bg_v_sum[1]->Add(h_bg_v[1][1]);

    std::cout << h_bg_v_sum[0]->Integral() << std::endl;
    std::cout << h_bg_v[0][0]->Integral() << std::endl;
    std::cout << h_bg_v[1][0]->Integral() << std::endl;
  
    std::cout << h_bg_v_sum[1]->Integral() << std::endl;
    std::cout << h_bg_v[0][1]->Integral() << std::endl;
    std::cout << h_bg_v[1][1]->Integral() << std::endl;

    h_bg_v_sum[0]->SetLineColor(1);
    h_bg_v_sum[1]->SetLineColor(2);

    h_bg_r[0] = (TH1D*) h_bg_v_sum[0]->Clone("h_bg_r_up");
    h_bg_r[1] = (TH1D*) h_bg_v_sum[1]->Clone("h_bg_r_dn");
    
    h_bg_r[0]->Sumw2();
    h_bg_r[1]->Sumw2();
    
    h_bg_r[0]->Divide(h_bg_sum);
    h_bg_r[1]->Divide(h_bg_sum);
    
//    for(int b=1; b<=nBins; ++b){
//      h_bg_r[0]->SetBinError(b, 1e-3);
//      h_bg_r[1]->SetBinError(b, 1e-3);
//    }    

    TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
    c1->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
    pad1->SetBottomMargin(0.15);
    pad1->Draw();
    pad1->cd();

    float yMin = 0.;
    float yMax1 = h_data->GetMaximum()*1.5;
    float yMax2 = 1.5*(bgStack.GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( h_data->GetNbinsX()<2 ) yMax *=3.;
    
    if(logY) {
      gPad->SetLogy();
      yMin=1e-1;
      yMax*=50.;
    }

    label="M_{T2} [GeV]";
    if(thisRegion.nJetsMax()==1)
      label="H_{T} [GeV]";

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
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      regionText->AddText( niceNames[i].c_str() );
      regionText->Draw("same");

    }
 
    TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+nVariations+1)*0.04, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    
    for( unsigned i=0; i<bgYields.size(); ++i )
      legend->AddEntry( h_bg[i], bgYields[i]->getFullName().c_str(), "F" );
    
    legend->AddEntry( h_bg_v_sum[0], "Variation UP", "L");
    legend->AddEntry( h_bg_v_sum[1], "Variation DN", "L");

    legend->Draw("same");
    bgStack.Draw("histo same");
    bgStack.SetMinimum(yMin);
    bgStack.Draw("histo same");

    for(int v=0; v<nVariations; ++v)
	h_bg_v_sum[v]->Draw("histo,same");
 
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
    
    for(int v=0; v<nVariations; ++v){
      
      h_bg_r[v]->SetStats(0);
      h_bg_r[v]->GetXaxis()->SetLabelSize(0.00);
      h_bg_r[v]->GetXaxis()->SetTickLength(0.09);
      h_bg_r[v]->GetYaxis()->SetNdivisions(5,5,0);
      h_bg_r[v]->GetYaxis()->SetRangeUser(0.75,1.25);
      h_bg_r[v]->GetYaxis()->SetTitleSize(0.17);
      h_bg_r[v]->GetYaxis()->SetTitleOffset(0.4);
      h_bg_r[v]->GetYaxis()->SetLabelSize(0.17);
      h_bg_r[v]->GetYaxis()->SetTitle("ratio");
      
    }
    
    h_bg_r[0]->Draw("EP");
    h_bg_r[1]->Draw("EP,same");
    

    gPad->RedrawAxis();
      
    c1->cd();
    
    c1->SaveAs( Form("%s/%s_%s.eps", outputdir.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.png", outputdir.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.pdf", outputdir.c_str(), var.c_str(), thisRegion.getName().c_str()) );

    delete c1;
    delete h_axes;
    
    delete h_data;
    
    delete h_bg_sum;
    delete h_bg_v_sum[0];
    delete h_bg_v_sum[1];
    delete h_bg_r[0];
    delete h_bg_r[1];
    
    delete bins;
    
    for(unsigned i=0; i<bgYields.size(); ++i){

      delete h_bg[i];
      
      for(int v=0; v<nVariations; ++v){
	delete h_bg_v[i][v];
      }
    
    }
    
  }// for MT2 regions
  
}
