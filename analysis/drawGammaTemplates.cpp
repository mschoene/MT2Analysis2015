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



bool useMC = true;


void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2EstimateZinvGamma>* templ, bool isPrompt, const std::string& suffix="" );
void compareSigmaIetaIeta( const std::string& outputdir, TH1D* prompt, TH1D* fake );
void drawSingleSigmaPlot( const std::string& outputdir, const std::string& label, float xMin, float xMax, TH1D* prompt, TH1D* fake );
void compareTemplatesVsMt2( const std::string& outputdir, MT2Analysis<MT2EstimateZinvGamma>* templ, bool isPrompt );


int main( int argc, char* argv[] ) {


  std::string samples = "PHYS14_v2_Zinv";
  //std::string samples = "CSA14_Zinv";
  if( argc>1 ) {
    std::string samples_tmp(argv[1]); 
    samples = samples_tmp;
  }


  //std::string regionsSet = "13TeV_inclusive";
  std::string regionsSet = "13TeV_CSA14";

  MT2DrawTools::setStyle();

  std::string mc_or_not = (useMC) ? "MC" : "";

  std::string outputdir = "Plots_GammaTemplates" + mc_or_not + "_" + samples + "_" + regionsSet;
  system( Form("mkdir -p %s", outputdir.c_str()) );

  MT2Analysis<MT2EstimateZinvGamma>* templatesPrompt = MT2Analysis<MT2EstimateZinvGamma>::readFromFile("gammaTemplates" + mc_or_not + "_" + samples + "_" + regionsSet + ".root", "templatesPrompt");
  MT2Analysis<MT2EstimateZinvGamma>* templatesFake   = MT2Analysis<MT2EstimateZinvGamma>::readFromFile("gammaTemplates" + mc_or_not + "_" + samples + "_" + regionsSet + ".root", "templatesFake");


  if( regionsSet=="13TeV_inclusive" ) {
    compareSigmaIetaIeta( outputdir, templatesPrompt->get(*(templatesPrompt->getRegions().begin()))->sietaieta, templatesFake->get(*(templatesFake->getRegions().begin()))->sietaieta );
    compareTemplatesVsMt2( outputdir, templatesPrompt, true );
    compareTemplatesVsMt2( outputdir, templatesFake, false );
  }

  

  std::vector<MT2Region> r_lowHT_vs_njet;
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2,  3, 0, 0 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, -1, 0, 0 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2,  3, 1, 1 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, -1, 1, 1 ) );

  compareRegions( outputdir, r_lowHT_vs_njet, templatesPrompt, true  );
  compareRegions( outputdir, r_lowHT_vs_njet, templatesFake , false );

  std::vector<MT2Region> r_medHT_vs_njet;
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2,  3, 0, 0 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, -1, 0, 0 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2,  3, 1, 1 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, -1, 1, 1 ) );

  compareRegions( outputdir, r_medHT_vs_njet, templatesPrompt, true );
  compareRegions( outputdir, r_medHT_vs_njet, templatesFake, false );

  std::vector<MT2Region> r_highHT_vs_njet;
  r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 2,  3, 0, 0 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 4, -1, 0, 0 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 2,  3, 1, 1 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 4, -1, 1, 1 ) );

  compareRegions( outputdir, r_highHT_vs_njet, templatesPrompt, true );
  compareRegions( outputdir, r_highHT_vs_njet, templatesFake, false );


  //std::vector<MT2Region> r_medHT_vs_nbjet;
  //r_medHT_vs_nbjet.push_back( MT2Region( 575., 1000., 200., 2, 3, 0, 0 ) );
  //r_medHT_vs_nbjet.push_back( MT2Region( 575., 1000., 200., 2, 3, 1, 1 ) );
  //r_medHT_vs_nbjet.push_back( MT2Region( 575., 1000., 200., 2, 3, 2, 2, 200., 200., 400., false ) );

  //compareRegions( outputdir, r_medHT_vs_nbjet, templatesPrompt, templatesFake, "_vsb" );


  //std::vector<MT2Region> r_highHT_vs_nbjet;
  //r_highHT_vs_nbjet.push_back( MT2Region( 1000., -1., 30., 2, 3, 0, 0 ) );
  //r_highHT_vs_nbjet.push_back( MT2Region( 1000., -1., 30., 2, 3, 1, 1 ) );
  //r_highHT_vs_nbjet.push_back( MT2Region( 1000., -1., 30., 2, 3, 2, 2, 200., 200., 400., false ) );

  //compareRegions( outputdir, r_highHT_vs_nbjet, templatesPrompt, templatesFake, "_vsb" );


  std::vector<MT2Region> r_vsHT;
  //r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 30., 2, -1, 0, 0 ) );
  r_vsHT.push_back( MT2Region( 450., 575. , 2, 3, 0, 0 ) );
  r_vsHT.push_back( MT2Region( 575., 1000., 2, 3, 0, 0 ) );
  r_vsHT.push_back( MT2Region( 1000., -1. , 2, 3, 0, 0 ) );

  compareRegions( outputdir, r_vsHT, templatesPrompt, true );
  compareRegions( outputdir, r_vsHT, templatesFake , false );



  std::vector<MT2Region> r_vsHT2;
  //r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 30., 2, -1, 0, 0 ) );
  r_vsHT2.push_back( MT2Region( 450., 575. , 4, -1, 0, 0 ) );
  r_vsHT2.push_back( MT2Region( 575., 1000., 4, -1, 0, 0 ) );
  r_vsHT2.push_back( MT2Region( 1000., -1. , 4, -1, 0, 0 ) );

  compareRegions( outputdir, r_vsHT2, templatesPrompt, true );
  compareRegions( outputdir, r_vsHT2, templatesFake , false );

  return 0;

}




void compareSigmaIetaIeta( const std::string& outputdir, TH1D* prompt, TH1D* fake ) {

  drawSingleSigmaPlot( outputdir, "Barrel", 0.007, 0.0125, (TH1D*)prompt->Clone(), (TH1D*)fake->Clone() );
  drawSingleSigmaPlot( outputdir, "Endcap", 0.0197, 0.035, (TH1D*)prompt->Clone(), (TH1D*)fake->Clone() );

}



void drawSingleSigmaPlot( const std::string& outputdir, const std::string& label, float xMin, float xMax, TH1D* prompt, TH1D* fake ) {

  int binMin = prompt->FindBin( xMin );
  int binMax = prompt->FindBin( xMax );
  int nbins = binMax-binMin;

  TH1D* thisPrompt = new TH1D("thisPrompt", "", nbins, xMin, xMax );
  TH1D* thisFake = new TH1D("thisFake", "", nbins, xMin, xMax );

  for( unsigned i=0; i<nbins; ++i ) {
    thisPrompt->SetBinContent( i+1, prompt->GetBinContent( i+binMin ));
    thisFake->SetBinContent( i+1, fake->GetBinContent( i+binMin ));
  }

  if( label=="Endcap" ) {
    thisPrompt->Rebin(4);
    thisFake  ->Rebin(4);
  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float yMax_prompt = thisPrompt->GetMaximum()/thisPrompt->Integral();
  float yMax_fake = thisFake->GetMaximum()/thisFake->Integral();
  float yMax = (yMax_prompt>yMax_fake) ? yMax_prompt : yMax_fake;
  yMax *= 1.2;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
  h2_axes->SetXTitle( "Sigma Ieta Ieta" );
  h2_axes->SetYTitle( "Normalized to Unity" );
  h2_axes->Draw();

  thisPrompt->SetFillColor(29);
  thisPrompt->SetLineColor(kBlack);

  thisFake->SetLineColor(kRed+2);
  thisFake->SetLineWidth(2);

  thisPrompt->DrawNormalized("histo same");
  thisFake->DrawNormalized("histo same");

  TLegend* legend = new TLegend( 0.65, 0.7, 0.92, 0.92, label.c_str() );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( thisPrompt, "Prompt", "F" );
  legend->AddEntry( thisFake, "Fake", "L" );
  legend->Draw("Same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/sietaieta_%s.eps", outputdir.c_str(), label.c_str()) );
  c1->SaveAs( Form("%s/sietaieta_%s.png", outputdir.c_str(), label.c_str()) );
  c1->SaveAs( Form("%s/sietaieta_%s.pdf", outputdir.c_str(), label.c_str()) );

  c1->Clear();
  c1->SetLogy();

  c1->cd();

  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.0001, 2. );
  h2_axes_log->SetXTitle( "Sigma Ieta Ieta" );
  h2_axes_log->SetYTitle( "Normalized to Unity" );
  h2_axes_log->Draw();

  thisPrompt->DrawNormalized("histo same");
  thisFake->DrawNormalized("histo same");

  legend->Draw("Same");
  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/sietaieta_%s_log.eps", outputdir.c_str(), label.c_str()) );
  c1->SaveAs( Form("%s/sietaieta_%s_log.png", outputdir.c_str(), label.c_str()) );
  c1->SaveAs( Form("%s/sietaieta_%s_log.pdf", outputdir.c_str(), label.c_str()) );


  delete c1;
  delete h2_axes;
  delete h2_axes_log;
  delete thisPrompt;
  delete thisFake;
  
}




void compareTemplatesVsMt2( const std::string& outputdir, MT2Analysis<MT2EstimateZinvGamma>* templ, bool isPrompt ) {

  std::vector<int> colors;
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 38 );
  colors.push_back( 42 );
  colors.push_back( kGray+2 );

  MT2Region region = *(templ->getRegions().begin());
  std::vector<TH1D*> histos = templ->get(region)->iso_bins_hist;
  int nBins;
  double* bins;
  region.getBins(nBins, bins);


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  TCanvas* c1_log = new TCanvas( "c1_log", "", 600, 600 );
  c1_log->SetLogy();

  float xMax = templ->get(region)->iso->GetXaxis()->GetXmax();

  float yMax = (isPrompt) ? 1. : 0.2;
  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax, 10, 0., yMax );
  h2_axes->SetXTitle( "Photon Charged Isolation [GeV]" );
  h2_axes->SetYTitle( "Normalized to Unity" );
  c1->cd();
  h2_axes->Draw("");


  float yMax_log = (isPrompt) ? 5. : 50.;
  float yMin_log = (isPrompt) ? 0.0001 : 0.001;
  
  TH2D* h2_axes_log = new TH2D( "axes_log", "", 10, 0., xMax, 10, yMin_log, yMax_log );
  h2_axes_log->SetXTitle( "Photon Charged Isolation [GeV]" );
  h2_axes_log->SetYTitle( "Normalized to Unity" );
  c1_log->cd();
  h2_axes_log->Draw("");


  int maxHistos = 4;
  TLegend* legend = new TLegend( 0.45, 0.91-0.07*maxHistos, 0.9, 0.91 );
  legend->SetTextSize(0.038);
  legend->SetFillColor(0);
  

  for( unsigned i=0; i<histos.size(); ++i ) {

    if( i>=maxHistos ) break;

    histos[i]->SetLineColor(colors[i]);
    histos[i]->SetLineWidth(2);
    legend->AddEntry( histos[i], Form("%.0f < M_{T2} < %.0f GeV", bins[i], bins[i+1]), "L" );

    c1->cd();
    histos[i]->DrawNormalized("same");

    c1_log->cd();
    histos[i]->DrawNormalized("same");

  }


  TPaveText* labelTop = MT2DrawTools::getLabelTop();

  float xMinLabel, xMaxLabel, yMinLabel, yMaxLabel;
  if( isPrompt ) {
    xMinLabel = 0.7;
    yMinLabel = 0.2;
    xMaxLabel = 0.9;
    yMaxLabel = 0.3;
  } else {
    xMinLabel = 0.2;
    yMinLabel = 0.8;
    xMaxLabel = 0.49;
    yMaxLabel = 0.9;
  }

  TPaveText* labelPrompt = new TPaveText( xMinLabel, yMinLabel, xMaxLabel, yMaxLabel, "brNDC" );
  labelPrompt->SetFillColor(0);
  labelPrompt->SetTextSize(0.035);
  labelPrompt->SetTextAlign(11); // align left
  if( isPrompt )
    labelPrompt->AddText( "Prompt" );
  else
    labelPrompt->AddText( "Fake" );
  labelPrompt->AddText( "Photons" );


  c1->cd();
  labelTop->Draw("same");
  labelPrompt->Draw("same");
  legend->Draw("same");
  gPad->RedrawAxis();

  c1_log->cd();
  labelTop->Draw("same");
  labelPrompt->Draw("same");
  legend->Draw("same");
  gPad->RedrawAxis();


  std::string prefix = (isPrompt) ? "Prompt" : "Fake";

  c1->SaveAs( Form("%s/templ%s_vsMT2.eps", outputdir.c_str(), prefix.c_str()) );
  c1->SaveAs( Form("%s/templ%s_vsMT2.pdf", outputdir.c_str(), prefix.c_str()) );
  c1->SaveAs( Form("%s/templ%s_vsMT2.png", outputdir.c_str(), prefix.c_str()) );

  c1_log->SaveAs( Form("%s/templ%s_vsMT2_log.eps", outputdir.c_str(), prefix.c_str()) );
  c1_log->SaveAs( Form("%s/templ%s_vsMT2_log.pdf", outputdir.c_str(), prefix.c_str()) );
  c1_log->SaveAs( Form("%s/templ%s_vsMT2_log.png", outputdir.c_str(), prefix.c_str()) );


  delete c1;
  delete c1_log;
  delete h2_axes;
  delete h2_axes_log;



}


void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2EstimateZinvGamma>* templ, bool isPrompt, const std::string& suffix ) {
//void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2EstimateZinvGamma>* templ, MT2Analysis<MT2EstimateZinvGamma>* templ_qcd, const std::string& suffix ) {


  bool loopOnHT=true;

  if( regions.size()>1 ) 
    if( (*(regions[0].htRegion())) == (*(regions[1].htRegion())) ) loopOnHT=false;


  std::vector<int> colors;
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 38 );
  colors.push_back( 42 );
  colors.push_back( kRed );
  colors.push_back( kBlack );
  
  std::vector<int> markers;
  markers.push_back( 21 );
  markers.push_back( 20 );
  markers.push_back( 23 );
  markers.push_back( 24 );
  markers.push_back( 25 );
  markers.push_back( 26 );
  


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  TCanvas* c1_log = new TCanvas( "c1_log", "", 600, 600 );
  c1_log->SetLogy();

  float xMax = templ->get(*(templ->getRegions().begin()))->iso->GetXaxis()->GetXmax();

  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 1. );
  h2_axes->SetXTitle( "Photon Charged Isolation [GeV]" );
  h2_axes->SetYTitle( "Normalized to Unity" );
  c1->cd();
  h2_axes->Draw("");


  float yMax_log = (isPrompt) ? 5. : 50.;
  float yMin_log = (isPrompt) ? 0.0001 : 0.001;
  
  TH2D* h2_axes_log = new TH2D( "axes_log", "", 10, 0., xMax, 10, yMin_log, yMax_log );
  h2_axes_log->SetXTitle( "Photon Charged Isolation [GeV]" );
  h2_axes_log->SetYTitle( "Normalized to Unity" );
  c1_log->cd();
  h2_axes_log->Draw("");


  float xMinLabel, xMaxLabel, yMinLabel, yMaxLabel;
  if( isPrompt ) {
    xMinLabel = 0.7;
    yMinLabel = 0.2;
    xMaxLabel = 0.9;
    yMaxLabel = 0.3;
  } else {
    xMinLabel = 0.2;
    yMinLabel = 0.8;
    xMaxLabel = 0.49;
    yMaxLabel = 0.9;
  }

  TPaveText* labelPrompt = new TPaveText( xMinLabel, yMinLabel, xMaxLabel, yMaxLabel, "brNDC" );
  labelPrompt->SetFillColor(0);
  labelPrompt->SetTextSize(0.035);
  labelPrompt->SetTextAlign(11); // align left
  if( isPrompt )
    labelPrompt->AddText( "Prompt" );
  else
    labelPrompt->AddText( "Fake" );
  labelPrompt->AddText( "Photons" );
  c1->cd();
  labelPrompt->Draw("same");


  float xMinLabelLog = xMinLabel;
  float xMaxLabelLog = xMaxLabel;
  if( isPrompt ) {
    xMinLabelLog = 0.2;
    xMaxLabelLog = 0.38;
  }

  TPaveText* labelPrompt_log = new TPaveText( xMinLabelLog, yMinLabel, xMaxLabelLog, yMaxLabel, "brNDC" );
  labelPrompt_log->SetFillColor(0);
  labelPrompt_log->SetTextSize(0.035);
  labelPrompt_log->SetTextAlign(11); // align left
  if( isPrompt )
    labelPrompt_log->AddText( "Prompt" );
  else
    labelPrompt_log->AddText( "Fake" );
  labelPrompt_log->AddText( "Photons" );
  c1_log->cd();
  labelPrompt_log->Draw("same");


  std::string legendTitle = (loopOnHT) ? regions[0].sigRegion()->getNiceName() : regions[0].htRegion()->getNiceName();

  TLegend* legend = new TLegend( 0.45, 0.9-0.08*regions.size(), 0.9, 0.9, legendTitle.c_str() );
  legend->SetTextSize(0.038);
  legend->SetFillColor(0);
  

  for( unsigned i=0; i<regions.size(); ++i ) {

    MT2EstimateZinvGamma* thisEstimate = templ->get( regions[i] );
    if( thisEstimate==0 ) {
      std::cout << "ERROR! Didn't find estimate for region: " << regions[i].getName() << " ! Exiting." << std::endl;
      exit(119);
    }

    TH1D* thisTempl = (TH1D*)(thisEstimate->iso->Clone());
    //if( !isPrompt )
    //  thisTempl->Rebin(2);
    thisTempl->SetLineColor(colors[i]);
    thisTempl->SetLineWidth(2);

    c1->cd();
    thisTempl->DrawNormalized("same" );
    c1_log->cd();
    thisTempl->DrawNormalized("same" );

    if( loopOnHT )
      legend->AddEntry( thisTempl, regions[i].htRegion()->getNiceName().c_str(), "L" );
    else
      legend->AddEntry( thisTempl, regions[i].sigRegion()->getNiceName().c_str(), "L" );

  }



  TPaveText* labelTop = MT2DrawTools::getLabelTop();



  c1->cd();
  labelTop->Draw("same");
  legend->Draw("same");
  gPad->RedrawAxis();

  c1_log->cd();
  labelTop->Draw("same");
  legend->Draw("same");
  gPad->RedrawAxis();


  std::string saveName = (loopOnHT) ? regions[0].sigRegion()->getName() : regions[0].htRegion()->getName();

  std::string prefix = (isPrompt) ? "Prompt" : "Fake";

  c1->SaveAs( Form("%s/templ%s_%s%s.eps", outputdir.c_str(), prefix.c_str(), saveName.c_str(), suffix.c_str()) );
  c1->SaveAs( Form("%s/templ%s_%s%s.pdf", outputdir.c_str(), prefix.c_str(), saveName.c_str(), suffix.c_str()) );
  c1->SaveAs( Form("%s/templ%s_%s%s.png", outputdir.c_str(), prefix.c_str(), saveName.c_str(), suffix.c_str()) );

  c1_log->SaveAs( Form("%s/templ%s_%s%s_log.eps", outputdir.c_str(), prefix.c_str(), saveName.c_str(), suffix.c_str()) );
  c1_log->SaveAs( Form("%s/templ%s_%s%s_log.pdf", outputdir.c_str(), prefix.c_str(), saveName.c_str(), suffix.c_str()) );
  c1_log->SaveAs( Form("%s/templ%s_%s%s_log.png", outputdir.c_str(), prefix.c_str(), saveName.c_str(), suffix.c_str()) );

  delete c1;
  delete c1_log;
  delete h2_axes;
  delete h2_axes_log;

}
