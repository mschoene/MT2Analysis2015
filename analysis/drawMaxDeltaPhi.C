{

  gStyle->SetOptStat(0);

  TFile* fdata   = TFile::Open("deltaPhiMax_data_fullPreselection.root");
  TFile* ffilter = TFile::Open("deltaPhiMax_data_HBHEfilter_all.root");
  TFile* fbg     = TFile::Open("deltaPhiMax_bg_fullPreselection.root");
  TFile* fsignal = TFile::Open("deltaPhiMax_signal_fullPreselection.root");

  TH1D* h[4];
  
  h[0] = (TH1D*) fdata->Get("hdeltaPhiMax_pass")->Clone("hdata");
  h[0]->SetLineColor(1);
  h[0]->SetMarkerColor(1);
  h[0]->SetMarkerStyle(20);

  h[1] = (TH1D*) fbg->Get("hdeltaPhiMax_pass")->Clone("hbg");
  h[1]->SetLineColor(kViolet);
  h[1]->SetFillColor(kViolet);
  
  h[2] = (TH1D*) fsignal->Get("hdeltaPhiMax_pass")->Clone("hsignal");
  h[2]->SetLineColor(1);
  h[2]->SetLineStyle(2);
  h[2]->SetLineWidth(2);

  h[3] = (TH1D*) ffilter->Get("hdeltaPhiMax_pass")->Clone("hdata");
  h[3]->SetLineColor(1);
  h[3]->SetMarkerColor(1);
  h[3]->SetMarkerStyle(20);

  TCanvas* c = new TCanvas("c", "", 600, 600);
  c->cd();

  h[0]->GetXaxis()->SetTitle("max #Delta#Phi(jets,MET)");
  h[0]->GetYaxis()->SetTitle("Events/0.01 (normalized)");
  h[0]->GetYaxis()->SetTitleOffset(1.4);

  h[3]->Scale(1.0/(h[0]->Integral()));

  for (int i=0; i<3; ++i)
    h[i]->Scale(1.0/(h[i]->Integral()));
    
  h[0]->Add(h[3], -1.0);

  h[0]->GetYaxis()->SetRangeUser(0., 0.25);
  h[0]->GetXaxis()->SetRangeUser(1.5, 3.2);
  
  h[0]->Draw("EP");
  h[1]->Draw("same");
  h[2]->Draw("same");
  h[0]->Draw("EP,same");

  gPad->RedrawAxis();

  float pi = TMath::Pi();
  TLine* l=new TLine(pi-0.015, 0, pi-0.015, 0.25);
  l->SetLineColor(2);
  l->SetLineWidth(2);
  l->Draw("same");

  TLegend* leg = new TLegend(0.15, 0.7, 0.4, 0.85);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(h[0], "Data", "P");
  leg->AddEntry(h[1], "Full SM", "F");
  leg->AddEntry(h[2], "T1qqqq", "l");

  leg->Draw("same");

  c->SaveAs("deltaPhiMax_fullPreselection_minusHBHE.pdf");

}
