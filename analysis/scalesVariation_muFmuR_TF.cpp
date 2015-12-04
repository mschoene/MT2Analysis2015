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

void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, std::vector<MT2Analysis<MT2EstimateTree>* > crYields,  std::string var, std::string label, std::string selection, bool logY );

MT2Analysis<MT2EstimateTree>* mergeYields( std::vector< MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max=-1, const std::string& legendName="" );

int main( int argc, char* argv[] ) {



  MT2DrawTools::setStyle();

  if( argc<2 ) {
    std::cout << "USAGE: ./scalesVariation_muFmuR_TF [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);
  
  lumi=cfg.lumi();
  
  std::string outputdir = cfg.getEventYieldDir() + "/scaleVariations_TF";
  system(Form("mkdir -p %s", outputdir.c_str()));


  TH1F::AddDirectory(kTRUE);

  std::string mcFile = cfg.getEventYieldDir() + "/analyses.root";
  std::string crFile = cfg.getEventYieldDir() + "/llepControlRegion/mc.root";
  std::string dataFile = cfg.getEventYieldDir() + "/analyses.root";


  MT2Analysis<MT2EstimateTree>* data  = MT2Analysis<MT2EstimateTree>::readFromFile( dataFile, "data" );
  
  MT2Analysis<MT2EstimateTree>* top      = MT2Analysis<MT2EstimateTree>::readFromFile( mcFile, "Top + W+jets");
  //  MT2Analysis<MT2EstimateTree>* wjets      = MT2Analysis<MT2EstimateTree>::readFromFile( mcFile, "WJets");
  

//  std::vector< MT2Analysis<MT2EstimateTree>* > llepSR_vec;
//  llepSR_vec.push_back( top );
//  llepSR_vec.push_back( wjets );
//
//  MT2Analysis<MT2EstimateTree>* llepSR = mergeYields( llepSR_vec, cfg.regionsSet(), "Top + W+jets", 100, 999 );
//  llepSR->setName("Top + W+jets");
//  llepSR->setFullName("Top + W+jets");
  
  MT2Analysis<MT2EstimateTree>* llepCR   = MT2Analysis<MT2EstimateTree>::readFromFile( crFile, "llepCR");
  
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields;
  bgYields.push_back(top);
  //bgYields.push_back(llepSR);

  std::vector<MT2Analysis<MT2EstimateTree>* > crYields;
  crYields.push_back(llepCR);

  drawHisto( outputdir, data, bgYields, crYields, "mt2", "M_{T2} [GeV]", "", kFALSE );

  return 0;

}// main



void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields, std::vector< MT2Analysis<MT2EstimateTree> *> crYields, std::string var, std::string label, std::string selection, bool logY ) {

  MT2DrawTools::setStyle();

  std::vector<int> colors; //mc
  colors.push_back(855); // top
    
  std::set<MT2Region> MT2Regions = data->getRegions();

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
   
    TH1F::AddDirectory(kTRUE);

    int nBins;
    double* bins;
    thisRegion.getBins(nBins, bins);

      
    TH1D* h_data = new TH1D("h_data", "", nBins, bins);

    std::string selection_ = Form("(mt2>200 && ht>200 && deltaPhiMin>0.3 && diffMetMht/met<0.5)*weight*%0.f", lumi);	  
    data->get( thisRegion )->tree->Project( h_data->GetName(), var.c_str(), selection_.c_str() );
  
    TH1D* h_bg[bgYields.size()];
    TH1D* h_cr[crYields.size()];
    TH1D* h_TF[bgYields.size()];
    
    int nVariations=9;
    int firstId=1;
    std::vector<std::string> variations;
    variations.push_back("#mu_{R}=1, #mu_{F}=1");
    variations.push_back("#mu_{R}=1, #mu_{F}=2");
    variations.push_back("#mu_{R}=1, #mu_{F}=0.5");
    variations.push_back("#mu_{R}=2, #mu_{F}=1");
    variations.push_back("#mu_{R}=2, #mu_{F}=2");
    variations.push_back("#mu_{R}=2, #mu_{F}=0.5");
    variations.push_back("#mu_{R}=0.5, #mu_{F}=1");
    variations.push_back("#mu_{R}=0.5, #mu_{F}=2");
    variations.push_back("#mu_{R}=0.5, #mu_{F}=0.5");
    
    TH1D* h_bg_v[bgYields.size()][nVariations];
    TH1D* h_cr_v[bgYields.size()][nVariations];
    TH1D* h_TF_v[bgYields.size()][nVariations];
    TH1D* h_TF_r[bgYields.size()][nVariations];

    for( unsigned i=0; i<bgYields.size(); ++i ) {

      int index=i;
      
      h_bg[index] = new TH1D(Form("h_%s_%s", var.c_str(), bgYields[i]->getName().c_str()), "", 1, 0., 1500.);      
      h_bg[index]->Sumw2();

      h_cr[index] = new TH1D(Form("h_%s_%s", var.c_str(), crYields[i]->getName().c_str()), "", 1, 0., 1500.);      
      h_cr[index]->Sumw2();
      
      selection_ = Form("(mt2>200 && ht>200 && deltaPhiMin>0.3 && diffMetMht/met<0.5 && LHEweight_id==%d)*LHEweight_wgt/LHEweight_original*weight*%0.f", firstId, lumi);
      bgYields[index]->get(thisRegion)->tree->Project( h_bg[index]->GetName(), var.c_str(), selection_.c_str() );
      crYields[index]->get(thisRegion)->tree->Project( h_cr[index]->GetName(), var.c_str(), selection_.c_str() );

      std::cout << h_bg[index]->Integral() << std::endl;

      h_TF[index] = (TH1D*) h_bg[index]->Clone(Form("h_%s_%s_r", var.c_str(), bgYields[i]->getName().c_str()));
      h_TF[index]->Divide(h_cr[index]);

      std::vector<int> colors_v;
      colors_v.push_back(2);
      colors_v.push_back(3);
      colors_v.push_back(4);
      colors_v.push_back(5);
      colors_v.push_back(6);
      colors_v.push_back(7);
      colors_v.push_back(8);
      colors_v.push_back(9);
      colors_v.push_back(30);
      
      for(int v=0; v < nVariations; ++v){
	
	h_bg_v[index][v] = new TH1D(Form("h_%s_%s_v%d", var.c_str(), bgYields[i]->getName().c_str(), v+firstId), "", 1, 0., 1500.);
	h_bg_v[index][v]->Sumw2();

	h_cr_v[index][v] = new TH1D(Form("h_%s_%s_v%d", var.c_str(), crYields[i]->getName().c_str(), v+firstId), "", 1, 0., 1500.);
	h_cr_v[index][v]->Sumw2();
	
	selection_ = Form("(mt2>200 && ht>200 && deltaPhiMin>0.3 && diffMetMht/met<0.5 && LHEweight_id==%d)*LHEweight_wgt/LHEweight_original*weight*%0.f", v+firstId, lumi);
	bgYields[index]->get(thisRegion)->tree->Project( h_bg_v[index][v]->GetName(), var.c_str(), selection_.c_str() );
	crYields[index]->get(thisRegion)->tree->Project( h_cr_v[index][v]->GetName(), var.c_str(), selection_.c_str() );
		
	h_bg_v[index][v]->SetLineColor(colors_v[v]);
	h_bg_v[index][v]->SetLineWidth(2);
	h_bg_v[index][v]->SetMarkerColor(colors_v[v]);

	h_TF_v[index][v]=(TH1D*) h_bg_v[index][v]->Clone(Form("h_%s_%s_tf%d", var.c_str(), bgYields[i]->getName().c_str(), v+firstId));
	h_TF_v[index][v]->Divide(h_cr_v[index][v]);

	h_TF_r[index][v]=(TH1D*) h_TF_v[index][v]->Clone(Form("h_%s_%s_tf_r%d", var.c_str(), bgYields[i]->getName().c_str(), v+firstId));
	h_TF_r[index][v]->Divide(h_TF[0]);

	for(int b=1; b<=1; ++b)
	  h_TF_r[index][v]->SetBinError(b, 1e-3);

      }

      h_TF[index]->SetFillColor( colors[index] );
      h_TF[index]->SetLineColor( kBlack );
      
    }
  
    TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
    c1->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
    pad1->SetBottomMargin(0.15);
    pad1->Draw();
    pad1->cd();

    float yMin = 0.;
    float yMax = h_TF[0]->GetMaximum()*5.0;
    
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
      legend->AddEntry( h_TF[i], bgYields[i]->getFullName().c_str(), "F" );
    
    for( int v=0; v<nVariations; ++v)
      legend->AddEntry( h_TF_v[0][v], variations[v].c_str(), "L");

    legend->Draw("same");
    h_TF[0]->Draw("histo,same");

    for( unsigned i=0; i<bgYields.size(); ++i )
      for(int v=0; v<nVariations; ++v)
	h_TF_v[i][v]->Draw("histo,same");
 
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

	h_TF_r[i][v]->SetStats(0);
	h_TF_r[i][v]->GetXaxis()->SetLabelSize(0.00);
	h_TF_r[i][v]->GetXaxis()->SetTickLength(0.09);
	h_TF_r[i][v]->GetYaxis()->SetNdivisions(5,5,0);
	h_TF_r[i][v]->GetYaxis()->SetRangeUser(0.75,1.25);
	h_TF_r[i][v]->GetYaxis()->SetTitleSize(0.17);
	h_TF_r[i][v]->GetYaxis()->SetTitleOffset(0.4);
	h_TF_r[i][v]->GetYaxis()->SetLabelSize(0.17);
	h_TF_r[i][v]->GetYaxis()->SetTitle("ratio");
       	
      }
    
    h_TF_r[0][0]->Draw("EP");
    
    for( unsigned i=0; i<bgYields.size(); ++i )
      for(int v=0; v<nVariations; ++v)
	h_TF_r[i][v]->Draw("EP,same");

    gPad->RedrawAxis();
      
    c1->cd();
    
    c1->SaveAs( Form("%s/%s_%s.eps", outputdir.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.png", outputdir.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.pdf", outputdir.c_str(), var.c_str(), thisRegion.getName().c_str()) );

    delete c1;
    delete h_axes;
    
    delete h_data;
    
    delete bins;
    
    for(unsigned i=0; i<bgYields.size(); ++i){

      delete h_bg[i];
      delete h_cr[i];
      delete h_TF[i];
      
      for(int v=0; v<nVariations; ++v){
	delete h_bg_v[i][v];
	delete h_cr_v[i][v];
	delete h_TF_v[i][v];
	delete h_TF_r[i][v];
      }
    
    }
    
  }// for MT2 regions
  
}


MT2Analysis<MT2EstimateTree>* mergeYields( std::vector<MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max, const std::string& legendName ) {

  if( id_max<0 ) id_max=id_min;

  MT2Analysis<MT2EstimateTree>* return_EventYield = new MT2Analysis<MT2EstimateTree>(name, regionsSet, id_min, legendName);

  for( unsigned i=0; i<EventYield.size(); ++i ) {

    if( EventYield[i]->getId() >= id_min && EventYield[i]->getId() <= id_max ) {

      *(return_EventYield) += *(EventYield[i]);

    }

  } // for EventYield                                                                                                                                                                                                           
  return return_EventYield;

}
