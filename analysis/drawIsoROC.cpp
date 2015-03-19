#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TChain.h"
#include "TF1.h"
#include "TProfile.h"
#include "TBox.h"
#include "THStack.h"

#include "../interface/MT2DrawTools.h"







void drawIsoRandomCone( const std::string& outputdir, TTree* tree_prompt, TTree* tree_nip, TTree* tree_fake, const std::string& eb_ee );
void drawCompare( const std::string& outputdir, const std::string& saveName, const std::string& axisName, std::vector<TH1D*> v, const std::string& eb_ee );
void drawSietaieta( const std::string& outputdir, TTree* tree_prompt, TTree* tree_nip, TTree* tree_fake, const std::string& eb_ee, float ptMin=-1., float ptMax=100000 );
void drawROC( const std::string& outputdir, TTree* tree_prompt, TTree* tree_nip, TTree* tree_fake, int optionNGI );
TGraph* getRoC( TH1D* h1_prompt, TH1D* h1_fake );
TGraph* getWP( TH1D* h1_prompt, TH1D* h1_fake, float thresh );
void drawTemplatesVsMT2( const std::string& outputdir, const std::string& varName, TTree* tree_prompt, TTree* tree_nip, TTree* tree_fake );
void drawVsMT2( const std::string& outputdir, const std::string& varName, const std::string& name, std::vector<TH1D*> histos, std::vector<float> bins );
void drawIsoVsSigma( const std::string& outputdir, TTree* tree_fake, const std::string& iso1, const std::string& iso2 );
std::string getLongName( const std::string& name );
void setBins( TH1D* h1, TH2D* h2 );



int main() {


  MT2DrawTools::setStyle();

  std::string outputdir = "IsoPlots";
  system(Form("mkdir -p %s", outputdir.c_str()));



  std::string fileName = "GenIsoCheck_PHYS14_v2_Zinv_13TeV_inclusive/genIso.root";

  TFile* file = TFile::Open(fileName.c_str());
  TTree* tree_prompt = (TTree*)file->Get("prompt/HT450toInf_j2toInf_b0toInf/tree_prompt_HT450toInf_j2toInf_b0toInf");
  TTree* tree_fake   = (TTree*)file->Get(  "fake/HT450toInf_j2toInf_b0toInf/tree_fake_HT450toInf_j2toInf_b0toInf");
  TTree* tree_nip    = (TTree*)file->Get(   "nip/HT450toInf_j2toInf_b0toInf/tree_nip_HT450toInf_j2toInf_b0toInf");


  std::cout << "-> Got stuff from file: " << fileName << std::endl;

  drawIsoRandomCone( outputdir, tree_prompt, tree_nip, tree_fake, "" );
  drawIsoRandomCone( outputdir, tree_prompt, tree_nip, tree_fake, "Barrel" );
  drawIsoRandomCone( outputdir, tree_prompt, tree_nip, tree_fake, "Endcap" );

  drawSietaieta( outputdir, tree_prompt, tree_nip, tree_fake, "Barrel" );
  drawSietaieta( outputdir, tree_prompt, tree_nip, tree_fake, "Endcap" );
  drawSietaieta( outputdir, tree_prompt, tree_nip, tree_fake, "Barrel", 200., 300. );
  drawSietaieta( outputdir, tree_prompt, tree_nip, tree_fake, "Barrel", 300., 400. );
  drawSietaieta( outputdir, tree_prompt, tree_nip, tree_fake, "Barrel", 400., 600. );
  drawSietaieta( outputdir, tree_prompt, tree_nip, tree_fake, "Barrel", 600., 800. );
  drawSietaieta( outputdir, tree_prompt, tree_nip, tree_fake, "Barrel", 800. );

  drawIsoVsSigma( outputdir, tree_fake, "iso", "isoCP" );


  drawROC( outputdir, tree_prompt, tree_nip, tree_fake, 0 );  // 0: put them in BG
  drawROC( outputdir, tree_prompt, tree_nip, tree_fake, 1 );  // 1: put them in signal
  drawROC( outputdir, tree_prompt, tree_nip, tree_fake, 2 );  // 2: ignore them

  drawTemplatesVsMT2( outputdir, "iso", tree_prompt, tree_nip, tree_fake );
  drawTemplatesVsMT2( outputdir, "isoCP", tree_prompt, tree_nip, tree_fake );


  return 0;

}




void drawIsoRandomCone( const std::string& outputdir, TTree* tree_prompt, TTree* tree_nip, TTree* tree_fake, const std::string& eb_ee ) {

  float etaMin;
  float etaMax;
  if( eb_ee=="Barrel" ) {
    etaMin = 0.;
    etaMax = 1.479;
  } else if( eb_ee=="Endcap" ) {
    etaMin = 1.479;
    etaMax = 2.5;
  } else if( eb_ee=="" ) {
    etaMin = 0.;
    etaMax = 2.5;
  }

  
  TH1D* h1_iso_prompt = new TH1D("iso_prompt", "", 100, 0., 20.);
  h1_iso_prompt->Sumw2();
  TH1D* h1_iso_nip = new TH1D("iso_nip", "", 100, 0., 20.);
  h1_iso_nip->Sumw2();
  TH1D* h1_iso_fake = new TH1D("iso_fake", "", 100, 0., 20.);
  h1_iso_fake->Sumw2();

  TH1D* h1_isoRC_prompt = new TH1D("isoRC_prompt", "", 100, 0., 20.);
  h1_isoRC_prompt->Sumw2();
  TH1D* h1_isoRC_nip = new TH1D("isoRC_nip", "", 100, 0., 20.);
  h1_isoRC_nip->Sumw2();
  TH1D* h1_isoRC_fake = new TH1D("isoRC_fake", "", 100, 0., 20.);
  h1_isoRC_fake->Sumw2();


  std::string cut(Form("weight*( abs(etaGamma)>=%f && abs(etaGamma)<%f)", etaMin, etaMax) );

  tree_prompt->Project( "iso_prompt"  , "iso*ptGamma", cut.c_str() );
  tree_prompt->Project( "isoRC_prompt", "isoRC", cut.c_str() );

  tree_nip->Project( "iso_nip"  , "iso*ptGamma", cut.c_str() );
  tree_nip->Project( "isoRC_nip", "isoRC", cut.c_str() );

  tree_fake->Project( "iso_fake"  , "iso*ptGamma", cut.c_str() );
  tree_fake->Project( "isoRC_fake", "isoRC", cut.c_str() );


  h1_iso_prompt->SetTitle("Prompt");
  h1_iso_nip->SetTitle("NIP");
  h1_iso_fake->SetTitle("Fake");

  h1_isoRC_prompt->SetTitle("Prompt");
  h1_isoRC_nip->SetTitle("NIP");
  h1_isoRC_fake->SetTitle("Fake");


  std::vector<TH1D*> v1;
  v1.push_back( h1_iso_prompt );
  v1.push_back( h1_iso_nip );
  v1.push_back( h1_iso_fake );
  drawCompare( outputdir, "iso", "Charged Isolation [GeV]", v1, eb_ee );

  std::vector<TH1D*> v1rc;
  v1rc.push_back( h1_isoRC_prompt );
  v1rc.push_back( h1_isoRC_nip );
  v1rc.push_back( h1_isoRC_fake );
  drawCompare( outputdir, "isoRC", "Random Cone Isolation [GeV]", v1rc, eb_ee );


  h1_isoRC_prompt->SetTitle("Random Cone");

  TH1D* h1_iso_prompt_plus_nip = new TH1D(*h1_iso_prompt);
  h1_iso_prompt_plus_nip->SetName("iso_prompt_plus_nip");
  h1_iso_prompt_plus_nip->Add(h1_iso_nip);
  h1_iso_prompt_plus_nip->SetTitle("Prompt+NIP");


  std::vector<TH1D*> v2;
  v2.push_back( h1_iso_prompt );
  v2.push_back( h1_iso_prompt_plus_nip );
  v2.push_back( h1_isoRC_prompt );
  drawCompare( outputdir, "isoRC_vs_iso", "Charged Isolation [GeV]", v2, eb_ee );


  delete h1_iso_prompt;
  delete h1_iso_nip;
  delete h1_iso_fake;

  delete h1_isoRC_prompt;
  delete h1_isoRC_nip;
  delete h1_isoRC_fake;


}



void drawCompare( const std::string& outputdir, const std::string& saveName, const std::string& axisName, std::vector<TH1D*> v, const std::string& eb_ee ) {


  std::vector<int> colors;
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 38 );
  colors.push_back( 42 );
  colors.push_back( kGray+2 );


  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();
  TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600 );
  c1_log->SetLogy();
  c1_log->cd();

  float yMax = 0.;
  for( unsigned i=0; i<v.size(); ++i ) {
    float thismax = v[i]->GetMaximum()/v[i]->Integral();
    if( thismax > yMax ) yMax = thismax;
  }
  yMax *= 1.2;

  float xMax = 10.;

  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax, 10, 0., yMax );  
  h2_axes->SetYTitle("Normalized to Unity");
  h2_axes->SetXTitle(axisName.c_str());
  c1->cd();
  h2_axes->Draw();

  TH2D* h2_axes_log = new TH2D( "axes_log", "", 10, 0., xMax, 10, 0.0001, 5.*yMax );  
  h2_axes_log->SetYTitle("Normalized to Unity");
  h2_axes_log->SetXTitle(axisName.c_str());
  c1_log->cd();
  h2_axes_log->Draw();

  float yMax_leg = 0.9;
  float yMin_leg = yMax_leg - (v.size()+1)*0.07;
  TLegend* legend;
  if( eb_ee!="" ) legend = new TLegend( 0.55, yMin_leg, 0.9, yMax_leg, eb_ee.c_str() );
  else            legend = new TLegend( 0.55, yMin_leg, 0.9, yMax_leg );
  legend->SetTextSize(0.038);
  legend->SetFillColor(0);

  for( unsigned i =0; i<v.size(); ++i ) {

    v[i]->SetLineWidth(2);
    v[i]->SetLineColor(colors[i]);

    c1->cd();
    v[i]->DrawNormalized("l same");

    c1_log->cd();
    v[i]->DrawNormalized("l same");

    legend->AddEntry( v[i], v[i]->GetTitle(), "L" );

  }

  TPaveText* labelTop = MT2DrawTools::getLabelTop();

  c1->cd();
  legend->Draw("same");
  labelTop->Draw("same");
  gPad->RedrawAxis();

  c1_log->cd();
  legend->Draw("same");
  labelTop->Draw("same");
  gPad->RedrawAxis();


  std::string eb_ee_canvas = (eb_ee!="") ? "_"+eb_ee : "";
  c1->SaveAs(Form("%s/%s%s.eps", outputdir.c_str(), saveName.c_str(), eb_ee_canvas.c_str()));
  c1->SaveAs(Form("%s/%s%s.png", outputdir.c_str(), saveName.c_str(), eb_ee_canvas.c_str()));
  c1->SaveAs(Form("%s/%s%s.pdf", outputdir.c_str(), saveName.c_str(), eb_ee_canvas.c_str()));

  c1_log->SaveAs(Form("%s/%s%s_log.eps", outputdir.c_str(), saveName.c_str(), eb_ee_canvas.c_str()));
  c1_log->SaveAs(Form("%s/%s%s_log.png", outputdir.c_str(), saveName.c_str(), eb_ee_canvas.c_str()));
  c1_log->SaveAs(Form("%s/%s%s_log.pdf", outputdir.c_str(), saveName.c_str(), eb_ee_canvas.c_str()));

  delete c1;
  delete h2_axes;
  delete c1_log;
  delete h2_axes_log;

}




void drawSietaieta( const std::string& outputdir, TTree* tree_prompt, TTree* tree_nip, TTree* tree_fake, const std::string& eb_ee, float ptMin, float ptMax ) {


  float etaMin;
  float etaMax;
  float xMin;
  float xMax;
  float xCut;
  float xSBmin;
  float xSBmax;
  int nBins;
  float xMinLegend;
  float xMaxLegend;
  if( eb_ee=="Barrel" ) {
    etaMin = 0.;
    etaMax = 1.479;
    xMin = 0.007;
    xMax = 0.02;
    nBins = 65;
    xCut = 0.010;
    xSBmin = 0.011;
    xSBmax = 0.015;
    xMinLegend = 0.68;
    xMaxLegend = 0.91;
  } else if( eb_ee=="Endcap" ) {
    etaMin = 1.479;
    etaMax = 2.5;
    xMin = 0.02;
    xMax = 0.035;
    xCut = 0.030;
    nBins = 50;
    xSBmin = 0.03;
    xSBmax = 0.035;
    xMinLegend = 0.45;
    xMaxLegend = 0.68;
  } else {
    std::cout << "Unkown ECAL region: " << eb_ee << std::endl;
    return;
  }
    

  TH1D* h1_prompt = new TH1D("hprompt", "", nBins, xMin, xMax );
  h1_prompt->Sumw2();
  TH1D* h1_nip = new TH1D("hnip", "", nBins, xMin, xMax );
  h1_nip->Sumw2();
  TH1D* h1_fake = new TH1D("hfake", "", nBins, xMin, xMax );
  h1_fake->Sumw2();


  std::string cut(Form("weight*(ptGamma>%f && ptGamma<%f && abs(etaGamma)>=%f && abs(etaGamma)<%f)", ptMin, ptMax, etaMin, etaMax));

  tree_prompt->Project( "hprompt", "sietaieta", cut.c_str() );
  tree_nip   ->Project( "hnip"   , "sietaieta", cut.c_str() );
  tree_fake  ->Project( "hfake"  , "sietaieta", cut.c_str() );

  h1_fake->Add(h1_nip);


  std::cout << std::endl;
  std::cout << eb_ee << ":" << std::endl;
  int binCut = h1_prompt->FindBin(xCut);
  int binSBmin = h1_prompt->FindBin(xSBmin);
  int binSBmax = h1_prompt->FindBin(xSBmax);
  std::cout << "Signal region: " << h1_prompt->Integral(1, binCut) + h1_fake->Integral(1, binCut) << std::endl;
  std::cout << "Sideband: " << h1_prompt->Integral(binSBmin, binSBmax) + h1_fake->Integral(binSBmin, binSBmax) << std::endl;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float yMax = (h1_prompt->GetMaximum()+h1_fake->GetMaximum())*1.1;
  if( ptMin>0. ) yMax = 750.;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax);
  h2_axes->SetXTitle("Photon #sigma_{i#eta i#eta}");
  h2_axes->SetYTitle("Events");
  h2_axes->Draw();

  h1_prompt->SetFillColor(kOrange+1);
  h1_prompt->SetLineColor(kBlack);

  h1_fake->SetFillColor(29);
  h1_fake->SetLineColor(kBlack);


  TBox* sbBox = new TBox( xSBmin, 0., xSBmax, yMax );
  sbBox->SetFillColor(kGray);
  sbBox->Draw("same");

  THStack* stack = new THStack();
  stack->Add( h1_fake ); 
  stack->Add( h1_prompt ); 
  stack->Draw("histo same");


  TLine* lineCut = new TLine( xCut, 0., xCut, yMax );
  lineCut->SetLineWidth(2);
  lineCut->Draw("same");


  TLine* lineSB1 = new TLine( xSBmin, 0., xSBmin, yMax );
  lineSB1->SetLineStyle(2);
  lineSB1->SetLineWidth(2);
  lineSB1->Draw("same");

  TLine* lineSB2 = new TLine( xSBmax, 0., xSBmax, yMax );
  lineSB2->SetLineStyle(2);
  lineSB2->SetLineWidth(2);
  lineSB2->Draw("same");

  TLegend* legend;
  if( ptMin>0. ) {
    std::string ptString;
    if( ptMax > 10000. ) ptString = std::string(Form("p_{T} > %.0f", ptMin));
    else                 ptString = std::string(Form("%.0f < p_{T} < %.0f", ptMin, ptMax));
    legend = new TLegend( xMinLegend, 0.6, xMaxLegend, 0.88, Form("#splitline{%s}{%s}", eb_ee.c_str(), ptString.c_str()) );
  } else {
    legend = new TLegend( xMinLegend, 0.67, xMaxLegend, 0.88, eb_ee.c_str());
  }
  legend->SetTextSize(0.035);
  legend->SetFillColor(0);
  legend->AddEntry( h1_prompt, "Prompt", "F" );
  legend->AddEntry( h1_fake, "Fake", "F" );
  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop(4.);
  labelTop->Draw("same");

  gPad->RedrawAxis();

  if( ptMin>0. ) {
    c1->SaveAs(Form("%s/sietaieta%s_pt%.0f_%.0f.eps", outputdir.c_str(), eb_ee.c_str(), ptMin, ptMax));
    c1->SaveAs(Form("%s/sietaieta%s_pt%.0f_%.0f.pdf", outputdir.c_str(), eb_ee.c_str(), ptMin, ptMax));
    c1->SaveAs(Form("%s/sietaieta%s_pt%.0f_%.0f.png", outputdir.c_str(), eb_ee.c_str(), ptMin, ptMax));
  } else {
    c1->SaveAs(Form("%s/sietaieta%s.eps", outputdir.c_str(), eb_ee.c_str()));
    c1->SaveAs(Form("%s/sietaieta%s.pdf", outputdir.c_str(), eb_ee.c_str()));
    c1->SaveAs(Form("%s/sietaieta%s.png", outputdir.c_str(), eb_ee.c_str()));
  }

  delete c1;
  delete h2_axes;
  delete h1_prompt;
  delete h1_fake;
  delete h1_nip;

}



void drawROC( const std::string& outputdir, TTree* tree_prompt, TTree* tree_nip, TTree* tree_fake, int optionNGI ) {


  int nbins = 1200;
  float xmin = 0.;
  //float xmax = 1.2; // rel
  float xmax = 200.; // abs

  TH1D* h1_iso_prompt = new TH1D("iso_prompt", "", nbins, xmin, xmax );
  h1_iso_prompt->Sumw2();
  TH1D* h1_iso_fake = new TH1D("iso_fake", "", nbins, xmin, xmax );
  h1_iso_fake->Sumw2();
  TH1D* h1_iso_nip = new TH1D("iso_nip", "", nbins, xmin, xmax );
  h1_iso_nip->Sumw2();

  TH1D* h1_isoCP_prompt = new TH1D("isoCP_prompt", "", nbins, xmin, xmax );
  h1_isoCP_prompt->Sumw2();
  TH1D* h1_isoCP_fake = new TH1D("isoCP_fake", "", nbins, xmin, xmax );
  h1_isoCP_fake->Sumw2();
  TH1D* h1_isoCP_nip = new TH1D("isoCP_nip", "", nbins, xmin, xmax );
  h1_isoCP_nip->Sumw2();

  TH1D* h1_isoCN_prompt = new TH1D("isoCN_prompt", "", nbins, xmin, xmax );
  h1_isoCN_prompt->Sumw2();
  TH1D* h1_isoCN_fake = new TH1D("isoCN_fake", "", nbins, xmin, xmax );
  h1_isoCN_fake->Sumw2();
  TH1D* h1_isoCN_nip = new TH1D("isoCN_nip", "", nbins, xmin, xmax );
  h1_isoCN_nip->Sumw2();

  TH1D* h1_isoCPN_prompt = new TH1D("isoCPN_prompt", "", nbins, xmin, xmax );
  h1_isoCPN_prompt->Sumw2();
  TH1D* h1_isoCPN_fake = new TH1D("isoCPN_fake", "", nbins, xmin, xmax );
  h1_isoCPN_fake->Sumw2();
  TH1D* h1_isoCPN_nip = new TH1D("isoCPN_nip", "", nbins, xmin, xmax );
  h1_isoCPN_nip->Sumw2();


  std::string sietaietaCut = "weight*( (abs(etaGamma)<1.479 && sietaieta<0.01) || (abs(etaGamma)>1.479 && sietaieta<0.03) )";

  tree_prompt->Project( "iso_prompt", "iso*ptGamma", sietaietaCut.c_str() );
  tree_fake  ->Project( "iso_fake"  , "iso*ptGamma", sietaietaCut.c_str() );
  tree_nip   ->Project( "iso_nip", "iso*ptGamma", sietaietaCut.c_str() );
  
  tree_prompt->Project( "isoCP_prompt", "isoCP*ptGamma", sietaietaCut.c_str() );
  tree_fake  ->Project( "isoCP_fake"  , "isoCP*ptGamma", sietaietaCut.c_str() );
  tree_nip   ->Project( "isoCP_nip", "isoCP*ptGamma", sietaietaCut.c_str() );
  
  tree_prompt->Project( "isoCN_prompt", "isoCN*ptGamma", sietaietaCut.c_str() );
  tree_fake  ->Project( "isoCN_fake"  , "isoCN*ptGamma", sietaietaCut.c_str() );
  tree_nip   ->Project( "isoCN_nip", "isoCN*ptGamma", sietaietaCut.c_str() );
  
  tree_prompt->Project( "isoCPN_prompt", "isoCPN*ptGamma", sietaietaCut.c_str() );
  tree_fake  ->Project( "isoCPN_fake"  , "isoCPN*ptGamma", sietaietaCut.c_str() );
  tree_nip   ->Project( "isoCPN_nip", "isoCPN*ptGamma", sietaietaCut.c_str() );

  if( optionNGI==0 ) {

    // putting NIP in BG:
    h1_iso_fake   ->Add(h1_iso_nip );   
    h1_isoCP_fake ->Add(h1_isoCP_nip ); 
    h1_isoCN_fake ->Add(h1_isoCN_nip ); 
    h1_isoCPN_fake->Add(h1_isoCPN_nip );
    

  } else if( optionNGI==1 ) {

    // putting NIP in signal:
    h1_iso_prompt   ->Add(h1_iso_nip );   
    h1_isoCP_prompt ->Add(h1_isoCP_nip ); 
    h1_isoCN_prompt ->Add(h1_isoCN_nip ); 
    h1_isoCPN_prompt->Add(h1_isoCPN_nip );


  } else if( optionNGI==2 ) {

    // ignoring NIP

  }


  TGraph* roc_iso = getRoC( h1_iso_prompt, h1_iso_fake );
  TGraph* roc_isoCP = getRoC( h1_isoCP_prompt, h1_isoCP_fake );
  TGraph* roc_isoCN = getRoC( h1_isoCN_prompt, h1_isoCN_fake );
  TGraph* roc_isoCPN = getRoC( h1_isoCPN_prompt, h1_isoCPN_fake );

  TGraph* wp_iso_loose = getWP( h1_iso_prompt, h1_iso_fake, 20.);
  TGraph* wp_iso_tight = getWP( h1_iso_prompt, h1_iso_fake, 2.5);

  TGraph* wp_isoCP_loose = getWP( h1_isoCP_prompt, h1_isoCP_fake, 60.);
  TGraph* wp_isoCP_tight = getWP( h1_isoCP_prompt, h1_isoCP_fake, 3.);


  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes;
  h2_axes = new TH2D("axes", "", 10, 0.75, 1.0001, 10, 0.9, 1.0001);
  h2_axes->SetXTitle("Fake Photon Rejection");
  h2_axes->SetYTitle("Prompt Photon Efficiency");
  h2_axes->Draw();


  roc_iso->SetMarkerStyle(20);
  roc_iso->SetMarkerSize(1.6);
  roc_iso->SetMarkerColor(kBlack);

  roc_isoCP->SetMarkerStyle(24);
  roc_isoCP->SetMarkerSize(1.6);
  roc_isoCP->SetMarkerColor(kBlack);

  roc_isoCN->SetMarkerStyle(20);
  roc_isoCN->SetMarkerSize(1.6);
  roc_isoCN->SetMarkerColor(46);

  roc_isoCPN->SetMarkerStyle(24);
  roc_isoCPN->SetMarkerSize(1.6);
  roc_isoCPN->SetMarkerColor(46);

  wp_iso_loose->SetMarkerStyle(30);
  wp_iso_tight->SetMarkerStyle(29);
  wp_iso_loose->SetMarkerSize(2.);
  wp_iso_tight->SetMarkerSize(2.);
  wp_iso_loose->SetMarkerColor(kOrange);
  wp_iso_tight->SetMarkerColor(kOrange);

  wp_isoCP_loose->SetMarkerStyle(30);
  wp_isoCP_tight->SetMarkerStyle(29);
  wp_isoCP_loose->SetMarkerSize(2.);
  wp_isoCP_tight->SetMarkerSize(2.);
  wp_isoCP_loose->SetMarkerColor(kBlue);
  wp_isoCP_tight->SetMarkerColor(kBlue);

  TLegend* legend = new TLegend( 0.2, 0.2, 0.5, 0.48 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( roc_iso, "Charged", "P" );
  legend->AddEntry( roc_isoCP, "Charged + Photon", "P" );
  legend->AddEntry( roc_isoCN, "Charged + Neutral", "P" );
  legend->AddEntry( roc_isoCPN, "Charged + Photon + Neutral", "P" );
  legend->Draw("same");

  TLine* diag = new TLine( 0., 1., 1., 0. );
  diag->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");


  roc_iso->Draw("psame");
  roc_isoCP->Draw("psame");
  roc_isoCN->Draw("psame");
  roc_isoCPN->Draw("psame");

  wp_isoCP_tight->Draw("psame");
  wp_iso_tight->Draw("psame");

  gPad->RedrawAxis();

  std::string suffix;
  if( optionNGI==0 ) suffix = "NGIbg";
  else if( optionNGI==1 ) suffix = "NGIsig";
  else if( optionNGI==2 ) suffix = "NGIig";

  c1->SaveAs(Form("%s/isoROC_%s.png", outputdir.c_str(), suffix.c_str()));
  c1->SaveAs(Form("%s/isoROC_%s.pdf", outputdir.c_str(), suffix.c_str()));
  c1->SaveAs(Form("%s/isoROC_%s.eps", outputdir.c_str(), suffix.c_str()));

  delete c1;
  delete h2_axes;

  delete h1_iso_prompt;
  delete h1_iso_fake;
  delete h1_iso_nip;

  delete h1_isoCP_prompt;
  delete h1_isoCP_fake;
  delete h1_isoCP_nip;

  delete h1_isoCN_prompt;
  delete h1_isoCN_fake;
  delete h1_isoCN_nip;

  delete h1_isoCPN_prompt;
  delete h1_isoCPN_fake;
  delete h1_isoCPN_nip;

}
  


void drawTemplatesVsMT2( const std::string& outputdir, const std::string& varName, TTree* tree_prompt, TTree* tree_nip, TTree* tree_fake ) {

  // now templates vs MT2
  std::vector<float> bins;
  bins.push_back(200.);
  bins.push_back(300.);
  bins.push_back(400.);
  bins.push_back(600.);
  bins.push_back(1000.);
  bins.push_back(1500.);

  std::vector<TH1D*> templates_prompt;
  std::vector<TH1D*> templates_fake;
  std::vector<TH1D*> templates_NIP;

  std::vector<TH1D*> templatesAbs_prompt;
  std::vector<TH1D*> templatesAbs_fake;
  std::vector<TH1D*> templatesAbs_NIP;

  float k = (varName=="iso") ? 1. : 2.;

  int nBinsPlusOne = 12;
  Double_t isoBins[nBinsPlusOne];
  isoBins[0]  = k*0.;
  isoBins[1]  = k*0.005;
  isoBins[2]  = k*0.01;
  isoBins[3]  = k*0.02;
  isoBins[4]  = k*0.03;
  isoBins[5]  = k*0.04;
  isoBins[6]  = k*0.05;
  isoBins[7]  = k*0.06;
  isoBins[8]  = k*0.07;
  isoBins[9]  = k*0.08;
  isoBins[10] = k*0.09;
  isoBins[11] = k*0.1;



  for( unsigned i=0; i<bins.size()-1; ++i ) {

    std::string promptName(Form("prompt%d", i));
    std::string fakeName(Form("fake%d", i));
    std::string NIPName(Form("NIP%d", i));

    TH1D* h1_prompt = new TH1D(promptName.c_str(), "", nBinsPlusOne-1, isoBins);
    TH1D* h1_fake   = new TH1D(  fakeName.c_str(), "", nBinsPlusOne-1, isoBins);
    TH1D* h1_NIP   = new TH1D(  NIPName.c_str(), "", nBinsPlusOne-1, isoBins);
    
    h1_prompt->Sumw2();
    h1_fake  ->Sumw2();
    h1_NIP  ->Sumw2();
    

    std::string cut(Form("weight*( mt2>%f && mt2<%f && (abs(etaGamma)<1.479 && sietaieta<0.01) || (abs(etaGamma)>1.479 && sietaieta<0.03) )", bins[i], bins[i+1]) );;

    tree_prompt->Project( promptName.c_str(), varName.c_str(), cut.c_str() );
    tree_nip   ->Project(    NIPName.c_str(), varName.c_str(), cut.c_str() );
    tree_fake  ->Project(   fakeName.c_str(), varName.c_str(), cut.c_str() );

    templates_prompt.push_back( h1_prompt );
    templates_fake  .push_back( h1_fake );
    templates_NIP   .push_back( h1_NIP );



    std::string promptNameAbs(Form("prompt_abs%d", i));
    std::string fakeNameAbs(Form("fake_abs%d", i));
    std::string NIPNameAbs(Form("NIP_abs%d", i));

    TH1D* h1_abs_prompt;
    TH1D* h1_abs_fake  ;
    TH1D* h1_abs_NIP  ;

    if( varName=="iso" ) {

      h1_abs_prompt = new TH1D(promptNameAbs.c_str(), "", 20, 0., 30.);
      h1_abs_fake   = new TH1D(  fakeNameAbs.c_str(), "", 20, 0., 30.);
      h1_abs_NIP   = new TH1D(  NIPNameAbs.c_str(), "", 20, 0., 30.);

    } else {

      h1_abs_prompt = new TH1D(promptNameAbs.c_str(), "", 30, 0., 60.);
      h1_abs_fake   = new TH1D(  fakeNameAbs.c_str(), "", 30, 0., 60.);
      h1_abs_NIP   = new TH1D(  NIPNameAbs.c_str(), "", 30, 0., 60.);

    }

    
    
    h1_abs_prompt->Sumw2();
    h1_abs_fake  ->Sumw2();
    h1_abs_NIP  ->Sumw2();
    
    tree_prompt->Project( promptNameAbs.c_str(), Form("%s*ptGamma", varName.c_str()), cut.c_str() );
    tree_nip   ->Project(    NIPNameAbs.c_str(), Form("%s*ptGamma", varName.c_str()), cut.c_str() );
    tree_fake  ->Project(   fakeNameAbs.c_str(), Form("%s*ptGamma", varName.c_str()), cut.c_str() );

    templatesAbs_prompt.push_back( h1_abs_prompt );
    templatesAbs_fake  .push_back( h1_abs_fake );
    templatesAbs_NIP  .push_back( h1_abs_NIP );

  }


  drawVsMT2( outputdir, varName, "prompt", templates_prompt, bins );
  drawVsMT2( outputdir, varName, "fake"  , templates_fake  , bins );
  drawVsMT2( outputdir, varName, "NIP"   , templates_NIP   , bins );

  drawVsMT2( outputdir, varName, "prompt", templatesAbs_prompt, bins );
  drawVsMT2( outputdir, varName, "fake"  , templatesAbs_fake  , bins );
  drawVsMT2( outputdir, varName, "NIP"   , templatesAbs_NIP   , bins );

  for( unsigned i=0; i<templates_prompt.size(); ++i ) {
    delete templates_prompt[i];
    delete templates_fake[i];
    delete templates_NIP[i];
    delete templatesAbs_prompt[i];
    delete templatesAbs_fake[i];
    delete templatesAbs_NIP[i];
  }

}



TGraph* getRoC( TH1D* h1_prompt, TH1D* h1_fake ) {

  TGraph* gr = new TGraph(0);

  int nBins = h1_prompt->GetNbinsX()+1;

  for( unsigned iBin=1; iBin<nBins; ++iBin ) {

    float eff_prompt = h1_prompt->Integral(1, iBin)/h1_prompt->Integral(1, nBins+1);
    float eff_fake   = h1_fake  ->Integral(1, iBin)/h1_fake  ->Integral(1, nBins+1);

  //  std::cout << "ibin: " << iBin << std::endl;
  //  std::cout <<  "eff_prompt: " <<  eff_prompt << std::endl;
  //  std::cout <<  "eff_fake  : " <<  eff_fake   << std::endl;

    gr->SetPoint( iBin-1, 1.-eff_fake, eff_prompt );

  }

  return gr;

}


void drawVsMT2( const std::string& outputdir, const std::string& varName, const std::string& name, std::vector<TH1D*> histos, std::vector<float> bins ) {

  bool isPrompt = (name!="fake");
  float xMax = histos[0]->GetXaxis()->GetXmax();
  bool isAbs = (xMax>2.);


  std::vector<int> colors;
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 38 );
  colors.push_back( 42 );
  colors.push_back( kGray+2 );
  

  TCanvas* c1 = new TCanvas( "c2", "", 600, 600 );
  TCanvas* c1_log = new TCanvas( "c2_log", "", 600, 600 );
  c1_log->SetLogy();


  std::string isoLongName = getLongName( varName );

  float yMax = (isPrompt) ? 1. : 0.2;
  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax, 10, 0., yMax );
  if( isAbs )
    h2_axes->SetXTitle( Form("%s Isolation [GeV]", isoLongName.c_str()) );
  else
    h2_axes->SetXTitle( Form("Relative %s Isolation", isoLongName.c_str()) );
  h2_axes->SetYTitle( "Normalized to Unity" );
  c1->cd();
  h2_axes->Draw("");


  float yMax_log, yMin_log;
  if( name=="prompt" ) {
    yMin_log = 0.0001;
    yMax_log = 5.;
  } else if( name=="fake" ) {
    yMin_log = 0.001;
    yMax_log = 50.;
  } else if( name=="NIP" ) {
    yMin_log = 0.001;
    yMax_log = 50.;
  }
  
  
  TH2D* h2_axes_log = new TH2D( "axes_log", "", 10, 0., xMax, 10, yMin_log, yMax_log );
  if( isAbs )
    h2_axes_log->SetXTitle( Form("%s Isolation [GeV]", isoLongName.c_str()) );
  else
    h2_axes_log->SetXTitle( Form("Relative %s Isolation", isoLongName.c_str()) );
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

  float xMin_label;
  float yMin_label;
  if( name=="prompt" ) {
    xMin_label = 0.75;
    yMin_label = 0.2;
  } else if( name=="fake" || name=="NIP" ) {
    xMin_label = 0.2;
    yMin_label = 0.8;
  }
  
  TPaveText* labelPrompt     = new TPaveText( 0.22, 0.8, 0.49, 0.9, "brNDC" );
  TPaveText* labelPrompt_log = new TPaveText( xMin_label, yMin_label, xMin_label+0.2, yMin_label+0.1, "brNDC" );
  labelPrompt   ->SetFillColor(0);
  labelPrompt   ->SetTextSize(0.035);
  labelPrompt   ->SetTextAlign(11); // align left
  labelPrompt_log->SetFillColor(0);
  labelPrompt_log->SetTextSize(0.035);
  labelPrompt_log->SetTextAlign(11); // align left
  if( name=="prompt" ) {
    labelPrompt    ->AddText( "Prompt" );
    labelPrompt_log->AddText( "Prompt" );
  } else if( name=="fake" ) {
    labelPrompt    ->AddText( "Fake" );
    labelPrompt_log->AddText( "Fake" );
  } else if( name=="NIP" ) {
    labelPrompt    ->AddText( "Non-GenIso" );
    labelPrompt_log->AddText( "Non-GenIso" );
  }
  labelPrompt    ->AddText( "Photons" );
  labelPrompt_log->AddText( "Photons" );


  c1->cd();
  labelTop->Draw("same");
  labelPrompt->Draw("same");
  legend->Draw("same");
  gPad->RedrawAxis();

  c1_log->cd();
  labelTop->Draw("same");
  labelPrompt_log->Draw("same");
  legend->Draw("same");
  gPad->RedrawAxis();


  std::string suffix = "";
  if( isAbs ) suffix = "Abs";

  c1->SaveAs( Form("%s/templ_%s%s_%s_vsMT2.eps", outputdir.c_str(), varName.c_str(), suffix.c_str(), name.c_str()) );
  c1->SaveAs( Form("%s/templ_%s%s_%s_vsMT2.pdf", outputdir.c_str(), varName.c_str(), suffix.c_str(), name.c_str()) );
  c1->SaveAs( Form("%s/templ_%s%s_%s_vsMT2.png", outputdir.c_str(), varName.c_str(), suffix.c_str(), name.c_str()) );

  c1_log->SaveAs( Form("%s/templ_%s%s_%s_vsMT2_log.eps", outputdir.c_str(), varName.c_str(), suffix.c_str(), name.c_str()) );
  c1_log->SaveAs( Form("%s/templ_%s%s_%s_vsMT2_log.pdf", outputdir.c_str(), varName.c_str(), suffix.c_str(), name.c_str()) );
  c1_log->SaveAs( Form("%s/templ_%s%s_%s_vsMT2_log.png", outputdir.c_str(), varName.c_str(), suffix.c_str(), name.c_str()) );


  delete c1;
  delete c1_log;
  delete h2_axes;
  delete h2_axes_log;


}




TGraph* getWP( TH1D* h1_prompt, TH1D* h1_fake, float thresh ) {

  int iBin = h1_prompt->FindBin(thresh);
  int nbins = h1_prompt->GetNbinsX() + 1;

  float eff      = h1_prompt->Integral( 1, iBin ) / h1_prompt->Integral( 1, nbins+1 );
  float eff_fake = h1_fake  ->Integral( 1, iBin ) / h1_fake  ->Integral( 1, nbins+1 );
  float rej = 1.-eff_fake;

  TGraph* gr = new TGraph(0);
  gr->SetPoint( 0, rej, eff );

  return gr;

}



void drawIsoVsSigma( const std::string& outputdir, TTree* tree_fake, const std::string& iso1, const std::string& iso2 ) {

  float xmin = 0.008;
  float xmax = 0.015;

  float iso1Max = 20.;
  //float iso2Max = 30.;
  float iso2Max = 60.;

  TH2D* h2_iso1_vs_sigma = new TH2D( "iso1_vs_sigma_2D", "", 20, xmin, xmax, 100, 0., iso1Max);
  TH2D* h2_iso2_vs_sigma = new TH2D( "iso2_vs_sigma_2D", "", 20, xmin, xmax, 100, 0., iso2Max);
  h2_iso1_vs_sigma->Sumw2();
  h2_iso2_vs_sigma->Sumw2();


  tree_fake->Project( "iso1_vs_sigma_2D", Form("%s*ptGamma:sietaieta", iso1.c_str()), Form("weight*(%s*ptGamma<%f)", iso1.c_str(), iso1Max) );
  tree_fake->Project( "iso2_vs_sigma_2D", Form("%s*ptGamma:sietaieta", iso2.c_str()), Form("weight*(%s*ptGamma<%f)", iso2.c_str(), iso2Max) );


  TH1D* h1_iso1_vs_sigma = new TH1D( "iso1_vs_sigma", "", 20, xmin, xmax );
  TH1D* h1_iso2_vs_sigma = new TH1D( "iso2_vs_sigma", "", 20, xmin, xmax );

  setBins( h1_iso1_vs_sigma, h2_iso1_vs_sigma );
  setBins( h1_iso2_vs_sigma, h2_iso2_vs_sigma );


  TF1* line1 = new TF1("line1", "[0] + [1]*x", xmin, xmax );
  line1->SetLineColor(kRed);
  line1->SetLineWidth(2);
  h1_iso1_vs_sigma->Fit(line1, "R+");
  line1->Draw("same");

  TF1* line2 = new TF1("line2", "[0] + [1]*x", xmin, xmax );
  line2->SetLineColor(kRed);
  line2->SetLineWidth(2);
  line2->SetLineStyle(2);
  h1_iso2_vs_sigma->Fit(line2, "R+");
  line2->Draw("same");


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  float yMax = 18.;

  TH2D* h2_axes = new TH2D("axes", "", 10, xmin, xmax, 10, 0., yMax );
  h2_axes->SetXTitle( "#sigma_{i#eta i#eta}" );
  h2_axes->SetYTitle( "Isolation [GeV]");
  h2_axes->Draw();


  h1_iso1_vs_sigma->SetMarkerStyle(20);
  h1_iso1_vs_sigma->SetMarkerSize(1.6);
  h1_iso1_vs_sigma->SetMarkerColor(kBlack);
  h1_iso1_vs_sigma->SetLineColor(kBlack);

  h1_iso2_vs_sigma->SetMarkerStyle(24);
  h1_iso2_vs_sigma->SetMarkerSize(1.6);
  h1_iso2_vs_sigma->SetMarkerColor(kBlack);
  h1_iso2_vs_sigma->SetLineColor(kBlack);

  TLine* lineCut = new TLine( 0.01, 0., 0.01, yMax );
  lineCut->SetLineColor(kBlack);
  lineCut->Draw("same");


  std::string longName1 = getLongName(iso1);
  std::string longName2 = getLongName(iso2);

  TLegend* legend = new TLegend( 0.52, 0.48, 0.9, 0.63 );
  legend->SetTextSize( 0.035 );
  legend->SetFillColor( 0 );
  legend->AddEntry( h1_iso1_vs_sigma, longName1.c_str(), "P" );
  legend->AddEntry( h1_iso2_vs_sigma, longName2.c_str(), "P" );
  legend->Draw("same");


  h1_iso1_vs_sigma->Draw("p same");
  h1_iso2_vs_sigma->Draw("p same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs(Form("%s/iso_vs_sigma.eps", outputdir.c_str()));  
  c1->SaveAs(Form("%s/iso_vs_sigma.pdf", outputdir.c_str()));  
  c1->SaveAs(Form("%s/iso_vs_sigma.png", outputdir.c_str()));  

  delete c1;
  delete h2_axes;

}




std::string getLongName( const std::string& name ) {

  std::string longName;
  if( name=="iso") 
    longName = "Charged";
  else if( name=="isoCP" ) 
    longName = "Charged + Photon";
  else if( name=="isoCN" ) 
    longName = "Charged + Neutral";
  else if( name=="isoCPN" ) 
    longName = "Full PF";

  return longName;

}



void setBins( TH1D* h1, TH2D* h2 ) {

  for( unsigned iBin=1; iBin<h2->GetNbinsX()+1; ++iBin ) {

    TH1D* thisProj = h2->ProjectionY(Form("%s_proj%d", h2->GetName(), iBin), iBin, iBin);

    h1->SetBinContent( iBin, thisProj->GetMean() );
    //h1->SetBinError( iBin, thisProj->GetRMS() );
    h1->SetBinError( iBin, thisProj->GetMeanError() );

  }

}
