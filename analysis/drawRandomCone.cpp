#include <iostream>


#include "interface/MT2DrawTools.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"




int main() {

  MT2DrawTools::setStyle();


  TFile* file15 = TFile::Open("RandomConeStudies/mt2_RC15.root");
  TFile* file20 = TFile::Open("RandomConeStudies/mt2_RC20.root");
  TFile* file30 = TFile::Open("RandomConeStudies/mt2_RC30.root");

  TTree* tree15 = (TTree*)file15->Get("tree");
  TTree* tree20 = (TTree*)file20->Get("tree");
  TTree* tree30 = (TTree*)file30->Get("tree");


  float xMax = 20.;

  TH1D* h1_15 = new TH1D("rc15", "", 50, 0., xMax);
  h1_15->Sumw2();
  TH1D* h1_20 = new TH1D("rc20", "", 50, 0., xMax);
  h1_20->Sumw2();
  TH1D* h1_30 = new TH1D("rc30", "", 50, 0., xMax);
  h1_30->Sumw2();


  tree15->Project("rc15", "gamma_chHadIsoRC[0]", "ngamma>0");
  tree20->Project("rc20", "gamma_chHadIsoRC[0]", "ngamma>0");
  tree30->Project("rc30", "gamma_chHadIsoRC[0]", "ngamma>0");

  h1_15->SetLineColor(46);
  h1_15->SetLineWidth(2);

  h1_20->SetLineColor(38);
  h1_20->SetLineWidth(2);

  h1_30->SetLineColor(29);
  h1_30->SetLineWidth(2);


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  TCanvas* c1_log = new TCanvas( "c1_log", "", 600, 600 );
  c1_log->SetLogy();
  

  float yMax = 6000.;
  
  TH2D* h2_axes = new TH2D("axes", "", 10, 0., xMax, 10, 0., yMax );
  h2_axes->SetXTitle( "Charged Isolation (Random Cone)");
  h2_axes->SetYTitle( "Entries" );
  
  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, 0., xMax, 10, 0.1, 2.*yMax );
  h2_axes_log->SetXTitle( "Charged Isolation (Random Cone)");
  h2_axes_log->SetYTitle( "Entries" );


  TLegend* legend = new TLegend( 0.45, 0.9-3*0.07, 0.9, 0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_20, "Jet Veto = 20 GeV", "L" );
  legend->AddEntry( h1_15, "Jet Veto = 15 GeV", "L" );
  legend->AddEntry( h1_30, "Jet Veto = 30 GeV", "L" );

  TPaveText* labelTop = MT2DrawTools::getLabelTop();


  c1->cd();
  h2_axes->Draw();
  legend->Draw("same");
  h1_15->Draw("same");
  h1_20->Draw("same");
  h1_30->Draw("same");
  labelTop->Draw("same");
  gPad->RedrawAxis();
  
  
  c1_log->cd();
  h2_axes_log->Draw();
  legend->Draw("same");
  h1_15->Draw("same");
  h1_20->Draw("same");
  h1_30->Draw("same");
  labelTop->Draw("same");
  gPad->RedrawAxis();

  c1->SaveAs("RandomConeStudies/rcJetThresh.eps");
  c1->SaveAs("RandomConeStudies/rcJetThresh.pdf");

  c1_log->SaveAs("RandomConeStudies/rcJetThresh_log.eps");
  c1_log->SaveAs("RandomConeStudies/rcJetThresh_log.pdf");

  return 0;

}
