#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"

#include "../interface/MT2DrawTools.h"



int main() {

  MT2DrawTools::setStyle();

  TFile* file_old = TFile::Open("QCD_13TeV_MGMLM_Spring15_bestMatching_angles_withNeutrinos.root");
  TFile* file_new = TFile::Open("templates.root");

  std::string outdir("plotsJetResponseTemplates");
  system(Form("mkdir -p %s", outdir.c_str()));

  for( unsigned ipt=0; ipt<20; ++ipt ) {

    for( unsigned ieta=0; ieta<12; ++ieta ) {

      std::string thisName(Form("h_tot_JetAll_ResponsePt_Pt%d_Eta%d", ipt, ieta));
      TH1D* templ_old = (TH1D*)file_old->Get(thisName.c_str());
      TH1D* templ_new = (TH1D*)file_new->Get(thisName.c_str());

      TCanvas* c1 = new TCanvas("c1", "", 600, 600);
      c1->cd();

      templ_old->SetLineColor(kBlue);
      templ_new->SetLineColor(kRed);

      templ_old->SetXTitle("Jet Response");
      templ_old->SetYTitle("Normalized to Unity");
      templ_old->DrawNormalized();
      templ_new->DrawNormalized("same");

      TLegend* legend = new TLegend( 0.7, 0.77, 0.9, 0.9 );
      legend->SetFillColor(0);
      legend->SetTextSize(0.035);
      legend->AddEntry( templ_old, "Old", "L" );
      legend->AddEntry( templ_new, "New", "L" );
      legend->Draw("same");

      MT2DrawTools::addLabels( c1, -1., "CMS Simulation" );
      gPad->RedrawAxis();

      c1->SaveAs(Form("%s/ResponsePt_Pt%d_Eta%d.eps", outdir.c_str(), ipt, ieta));
      c1->SaveAs(Form("%s/ResponsePt_Pt%d_Eta%d.pdf", outdir.c_str(), ipt, ieta));

      c1->Clear();
      c1->SetLogy();

      templ_old->DrawNormalized();
      templ_new->DrawNormalized("same");
      legend->Draw("same");

      MT2DrawTools::addLabels( c1, 0., "CMS Simulation" );
      gPad->RedrawAxis();

      c1->SaveAs(Form("%s/ResponsePt_Pt%d_Eta%d_log.eps", outdir.c_str(), ipt, ieta));
      c1->SaveAs(Form("%s/ResponsePt_Pt%d_Eta%d_log.pdf", outdir.c_str(), ipt, ieta));

      delete c1;

    } // for eta

  } // for pt

  return 0;
  
}
