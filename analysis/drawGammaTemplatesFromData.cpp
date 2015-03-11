#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2DrawTools.h"





void drawSinglePlot( const std::string& outputdir, const std::string& name, const MT2Region& region, TH1D* histo1, const std::string& name1, TH1D* histo2=0, const std::string& name2="", TH1D* histo3=0, const std::string& name3="" );


int main( int argc, char* argv[] ) {


  std::string samples = "PHYS14_v2_Zinv";
  if( argc>1 ) {
    std::string samples_tmp(argv[1]); 
    samples = samples_tmp;
  }


  std::string regionsSet = "13TeV_inclusive";
  //std::string regionsSet = "13TeV_CSA14";

  MT2DrawTools::setStyle();


  std::string outputdir = "Plots_GammaTemplatesFromData_" + samples + "_" + regionsSet;
  system( Form("mkdir -p %s", outputdir.c_str()) );

  MT2Analysis<MT2EstimateZinvGamma>* templatesPromptMC  = MT2Analysis<MT2EstimateZinvGamma>::readFromFile("gammaTemplatesMC_" + samples + "_" + regionsSet + ".root", "templatesPrompt");
  MT2Analysis<MT2EstimateZinvGamma>* templatesFakeMC    = MT2Analysis<MT2EstimateZinvGamma>::readFromFile("gammaTemplatesMC_" + samples + "_" + regionsSet + ".root", "templatesFake");

  MT2Analysis<MT2EstimateZinvGamma>* templatesPrompt    = MT2Analysis<MT2EstimateZinvGamma>::readFromFile("gammaTemplatesData_" + samples + "_" + regionsSet + ".root", "templatesPrompt");
  MT2Analysis<MT2EstimateZinvGamma>* templatesPromptRaw = MT2Analysis<MT2EstimateZinvGamma>::readFromFile("gammaTemplatesData_" + samples + "_" + regionsSet + ".root", "templatesPromptRaw");
  MT2Analysis<MT2EstimateZinvGamma>* templatesFake      = MT2Analysis<MT2EstimateZinvGamma>::readFromFile("gammaTemplatesData_" + samples + "_" + regionsSet + ".root", "templatesFake");

  std::set<MT2Region> regions = templatesPrompt->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    drawSinglePlot( outputdir, "Fake", *iR, templatesFake->get(*iR)->iso, "Data Fakes", 0, "", templatesFakeMC->get(*iR)->iso, "MC Fakes" );
    drawSinglePlot( outputdir, "Prompt", *iR, templatesPromptRaw->get(*iR)->iso, "Data (all)", templatesPrompt->get(*iR)->iso, "Data w/o Fakes", templatesPromptMC->get(*iR)->iso, "MC Prompts" );

  } 


  return 0;

}



void drawSinglePlot( const std::string& outputdir, const std::string& name, const MT2Region& region, TH1D* histo1, const std::string& name1, TH1D* histo2, const std::string& name2, TH1D* histo3, const std::string& name3 ) {


  bool rebin = false;


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  TCanvas* c1_log = new TCanvas( "c1_log", "", 600, 600 );
  c1_log->SetLogy();

  float yMinLegend = 0.75;
  float yMaxScale = 1.25;
  if( histo3!=0 ) {
    if( histo2!=0 ) {
      yMinLegend = 0.68;
      float int2 = histo2->Integral();
      float int3 = histo3->Integral();
      histo3->Scale( int2/int3 );
    } else {
      yMaxScale = 1.7;
      if( rebin ) yMaxScale*=2.;
      float int1 = histo1->Integral();
      float int3 = histo3->Integral();
      histo3->Scale( int1/int3 );
    }
  }

  float xMax = histo1->GetXaxis()->GetXmax();
  float yMax1 = histo1->GetMaximum();
  float yMax2 = (histo2!=0) ? histo2->GetMaximum() : 0.;
  float yMax3 = (histo3!=0) ? histo3->GetMaximum() : 0.;
  float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
  if( yMax3>yMax ) yMax = yMax3;
  yMax *= yMaxScale;

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., xMax, 10, 0., yMax );
  h2_axes->SetXTitle( "Photon Charged Isolation [GeV]" );
  h2_axes->SetYTitle( "Entries" );
  c1->cd();
  h2_axes->Draw();

  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, 0., xMax, 10, 1., 10.*yMax );
  h2_axes_log->SetXTitle( "Photon Charged Isolation [GeV]" );
  h2_axes_log->SetYTitle( "Entries" );
  c1_log->cd();
  h2_axes_log->Draw();

  TLegend* legend = new TLegend( 0.55, yMinLegend, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);

  histo1->SetMarkerStyle( 20 );
  histo1->SetMarkerSize( 1.6 );
  histo1->SetMarkerColor( kBlack );
  histo1->SetLineColor( kBlack );
  if( rebin ) histo1->Rebin(2);
  legend->AddEntry( histo1, name1.c_str(), "P" );

  if( histo2!=0 ) {
    histo2->SetMarkerStyle( 24 );
    histo2->SetMarkerSize( 1.6 );
    histo2->SetMarkerColor( kBlue );
    histo2->SetLineColor( kBlack );
    if( rebin ) histo2->Rebin(2);
    c1->cd();
    histo2->Draw("p same");
    c1_log->cd();
    histo2->Draw("p same");
    legend->AddEntry( histo2, name2.c_str(), "P" );
  }

  if( histo3!=0 ) {
    histo3->SetLineColor( 46 );
    histo3->SetLineWidth( 2 );
    if( rebin ) histo3->Rebin(2);
    c1->cd();
    histo3->Draw("L E same");
    c1_log->cd();
    histo3->Draw("L E same");
    legend->AddEntry( histo3, name3.c_str(), "L" );
  }

  TPaveText* labelTop = MT2DrawTools::getLabelTop(4.);

  c1->cd();
  histo1->Draw("p same" );
  legend->Draw("same");
  labelTop->Draw("same");
  gPad->RedrawAxis();


  c1_log->cd();
  histo1->Draw("p same" );
  legend->Draw("same");
  labelTop->Draw("same");
  gPad->RedrawAxis();


  c1->SaveAs(Form("%s/templateClosure%s_%s.eps", outputdir.c_str(), name.c_str(), region.getName().c_str()));
  c1->SaveAs(Form("%s/templateClosure%s_%s.png", outputdir.c_str(), name.c_str(), region.getName().c_str()));
  c1->SaveAs(Form("%s/templateClosure%s_%s.pdf", outputdir.c_str(), name.c_str(), region.getName().c_str()));

  c1_log->SaveAs(Form("%s/templateClosure%s_%s_log.eps", outputdir.c_str(), name.c_str(), region.getName().c_str()));
  c1_log->SaveAs(Form("%s/templateClosure%s_%s_log.png", outputdir.c_str(), name.c_str(), region.getName().c_str()));
  c1_log->SaveAs(Form("%s/templateClosure%s_%s_log.pdf", outputdir.c_str(), name.c_str(), region.getName().c_str()));

  delete c1;
  delete h2_axes;
  delete c1_log;
  delete h2_axes_log;

}
