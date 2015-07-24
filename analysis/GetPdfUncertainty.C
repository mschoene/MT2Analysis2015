{
  double lumi = 4.; //fb-1
 
  // set the TStyle
  TStyle* style = new TStyle("DrawBaseStyle", "");
  style->SetCanvasColor(0);
  style->SetPadColor(0);
  style->SetFrameFillColor(0);
  style->SetStatColor(0);
  style->SetOptStat(0);
  style->SetTitleFillColor(0);
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);
  style->SetFrameBorderMode(0);
  style->SetPadBottomMargin(0.12);
  style->SetPadLeftMargin(0.12);
  style->cd();

  // For the canvas:
  style->SetCanvasBorderMode(0);
  style->SetCanvasColor(kWhite);
  style->SetCanvasDefH(600); //Height of canvas
  style->SetCanvasDefW(600); //Width of canvas
  style->SetCanvasDefX(0); //POsition on screen
  style->SetCanvasDefY(0);

  // For the Pad:
  style->SetPadBorderMode(0);
  style->SetPadColor(kWhite);
  style->SetPadGridX(false);
  style->SetPadGridY(false);
  style->SetGridColor(0);
  style->SetGridStyle(3);
  style->SetGridWidth(1);
  
  // For the frame:
  style->SetFrameBorderMode(0);
  style->SetFrameBorderSize(1);
  style->SetFrameFillColor(0);
  style->SetFrameFillStyle(0);
  style->SetFrameLineColor(1);
  style->SetFrameLineStyle(1);
  style->SetFrameLineWidth(1);
  
  // Margins:
  style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.15);//0.13);
  style->SetPadLeftMargin(0.15);//0.16);
  style->SetPadRightMargin(0.05);//0.02);

  // For the Global title:
  style->SetOptTitle(0);
  style->SetTitleFont(42);
  style->SetTitleColor(1);
  style->SetTitleTextColor(1);
  style->SetTitleFillColor(10);
  style->SetTitleFontSize(0.05);

  // For the axis titles:
  style->SetTitleColor(1, "XYZ");
  style->SetTitleFont(42, "XYZ");
  style->SetTitleSize(0.05, "XYZ");
  style->SetTitleXOffset(1.15);//0.9);
  style->SetTitleYOffset(1.3); // => 1.15 if exponents

  // For the axis labels:
  style->SetLabelColor(1, "XYZ");
  style->SetLabelFont(42, "XYZ");
  style->SetLabelOffset(0.007, "XYZ");
  style->SetLabelSize(0.045, "XYZ");

  // For the axis:
  style->SetAxisColor(1, "XYZ");
  style->SetStripDecimals(kTRUE);
  style->SetTickLength(0.03, "XYZ");
  style->SetNdivisions(510, "XYZ");
  style->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  style->SetPadTickY(1);

  // for histograms:
  style->SetHistLineColor(1);

  // for the pallete
  Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red  [5] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[5] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue [5] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 100);
  style->SetNumberContours(100);

  style->cd();


  std::string outputdir = "scalesVariations_LO_Plots/";
  bool logY = false;

  std::vector<std::string> pdfSets;
  pdfSets.push_back("NNPDF3.0 LO, #alpha_{s}=0.130");
  pdfSets.push_back("MMHT2014 LO 68% CL");
  pdfSets.push_back("CTEQ61L");

  TFile* fNNPDF = new TFile("EventYields_mc_74X_v1_dummy_4fb_TTJets_LO/scalesVariations_NNPDF_LO_correctEquations/CT10.root");
  TFile* fMMHT  = new TFile("EventYields_mc_74X_v1_dummy_4fb_TTJets_LO/scalesVariations_MMHT_LO_correctEquations/CT10.root");
  TFile* fCT10  = new TFile("EventYields_mc_74X_v1_dummy_4fb_TTJets_LO/scalesVariations_CTEQ61L_correctEquations/CT10.root");

//  std::vector<std::string> pdfSets;
//  pdfSets.push_back("NNPDF3.0 NLO, #alpha_{s}=0.118");
//  pdfSets.push_back("MMHT2014 NLO 68% CL");
//  pdfSets.push_back("CT10 NLO");
//
//  TFile* fNNPDF = new TFile("EventYields_mc_74X_v1_dummy_4fb_TTJets_LO/scalesVariations_NNPDF30_NLO_shape_rescale/CT10.root");
//  TFile* fMMHT  = new TFile("EventYields_mc_74X_v1_dummy_4fb_TTJets_LO/scalesVariations_MMHT2014_NLO_shape_rescale/CT10.root");
//  TFile* fCT10  = new TFile("EventYields_mc_74X_v1_dummy_4fb_TTJets_LO/scalesVariations_CT10NLO_correctEquations/CT10.root");

  int h=0;
  int nHistos=44*4;

  TH1F* hNNPDF[nHistos];
  TIter next(fNNPDF->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    hNNPDF[h] = (TH1F*)key->ReadObj();
    ++h;
  }

  TH1F* hMMHT[nHistos];
  TIter next(fMMHT->GetListOfKeys());
  TKey *key;
  h=0;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    hMMHT[h] = (TH1F*)key->ReadObj();
    ++h;
  }
  
  TH1F* hCT10[nHistos];
  TIter next(fCT10->GetListOfKeys());
  TKey *key;
  h=0;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    hCT10[h] = (TH1F*)key->ReadObj();
    ++h;
  }

  int r=0;
  while( r < nHistos ){
      
    int nBins = hNNPDF[r+3]->GetNbinsX();
    
    hNNPDF[r+3]->SetLineColor(1);
    hNNPDF[r+3]->SetMarkerColor(1);
    
    hNNPDF[r+2]->SetLineColor(4);
    hNNPDF[r+2]->SetMarkerColor(4);

    hCT10[r+2]->SetLineColor(6);
    hCT10[r+2]->SetMarkerColor(6);
    
    hMMHT[r+2]->SetLineColor(2);
    hMMHT[r+2]->SetMarkerColor(2);
    
    for( int b=1; b<=nBins; ++b ){

      hNNPDF[r+2]->SetBinError(b, 0.);
      hCT10[r+2]->SetBinError(b, 0.);
      hMMHT[r+2]->SetBinError(b, 0.);

    }

    TGraphAsymmErrors* gNNPDF = new TGraphAsymmErrors(hNNPDF[r+2]);
    TGraphAsymmErrors* gCT10 =new TGraphAsymmErrors(hCT10[r+2]);
    TGraphAsymmErrors* gMMHT =new TGraphAsymmErrors(hMMHT[r+2]);

    TH1D* hPDF = (TH1D*) hNNPDF[r+3]->Clone("hPDF");
    hPDF->Divide(hNNPDF[r+3]);
    
    float max=-999.;
    float min=1.e+6;
    for( int b=1; b<=nBins; ++b ){
      
      float max=-999;
      float min=1.e+6;

      max = (hNNPDF[r]->GetBinContent(b) > max)  ? hNNPDF[r]->GetBinContent(b) : max;
      min = (hNNPDF[r+1]->GetBinContent(b) < min)  ? hNNPDF[r+1]->GetBinContent(b) : min;
      
      
      gNNPDF->SetPointEYhigh(b-1, (hNNPDF[r]->GetBinContent(b))-(hNNPDF[r+2]->GetBinContent(b)));
      gNNPDF->SetPointEYlow(b-1, -(hNNPDF[r+1]->GetBinContent(b))+(hNNPDF[r+2]->GetBinContent(b)));     
      
      max = (hMMHT[r]->GetBinContent(b) > max)  ? hMMHT[r]->GetBinContent(b) : max;
      min = (hMMHT[r+1]->GetBinContent(b) < min)  ? hMMHT[r+1]->GetBinContent(b) : min;

      gMMHT->SetPointEYhigh(b-1, (hMMHT[r]->GetBinContent(b))-(hMMHT[r+2]->GetBinContent(b)));
      gMMHT->SetPointEYlow(b-1, -(hMMHT[r+1]->GetBinContent(b))+(hMMHT[r+2]->GetBinContent(b)));     
      
      max = (hCT10[r]->GetBinContent(b) > max)  ? hCT10[r]->GetBinContent(b) : max;
      min = (hCT10[r+1]->GetBinContent(b) < min)? hCT10[r+1]->GetBinContent(b) : min;

      gCT10->SetPointEYhigh(b-1, (hCT10[r]->GetBinContent(b))-(hCT10[r+2]->GetBinContent(b)));
      gCT10->SetPointEYlow(b-1, -(hCT10[r+1]->GetBinContent(b))+(hCT10[r+2]->GetBinContent(b)));

      float pdfuncertainty = (max-min)/2.;
      float yield = (max+min)/2.;
      float pdfuncertainty_rel;

      if ( pdfuncertainty < 1e-6 ) pdfuncertainty = 0.;
      if ( pdfuncertainty < 1e-6 || yield < 1e-6 || pdfuncertainty/yield < 1e-6 ) pdfuncertainty_rel = 0.;
      else pdfuncertainty_rel = pdfuncertainty/yield;

      hPDF->SetBinError(b, pdfuncertainty_rel);
      std::cout << "Topo-Region #" << r/4 << ", MT2 bin #" << b <<"\t pdf uncertainty (%) = " << pdfuncertainty_rel*100.  << std::endl;
      std::cout << max << "\t" << min << std::endl;
      
    }
    
    TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
    c1->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
    pad1->SetBottomMargin(0.15);
    pad1->Draw();
    pad1->cd();

    float yMin = 0.;
    float yMax1 = 1.5*(hNNPDF[r+2]->GetMaximum() + sqrt(hNNPDF[r+2]->GetMaximum()));
    float yMax2 = 1.5*(hMMHT[r+2]->GetMaximum() + sqrt(hMMHT[r+2]->GetMaximum()));
    float yMax3 = 1.5*(hCT10[r+2]->GetMaximum() + sqrt(hCT10[r+2]->GetMaximum()));
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    yMax = (yMax>yMax3) ? yMax : yMax3;

    if( hNNPDF[r+2]->GetNbinsX()<2 ) yMax *=3.;

    if(logY) {
      gPad->SetLogy();
      yMin=1e-1;
      yMax*=50.;
    }

    std::string label = "M_{T2} [GeV]";
    std::string labelY = "Events";
    TH2D* h_axes = new TH2D("axes", "", 10, 200, 1500, 10, yMin, yMax );
    h_axes->GetXaxis()->SetTitle(label.c_str());
    h_axes->GetYaxis()->SetTitle(labelY.c_str());
    h_axes->GetYaxis()->SetTitleOffset(1.5);
    h_axes->GetYaxis()->SetLabelSize(0.04);
    h_axes->Draw();
    
    std::vector<std::string> niceNames;
    niceNames.push_back(Form("Topo-Region #%.0f", r/4.));
    
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

    char text[300];
    sprintf( text, "CMS Simulation, %.1f fb^{-1} at #sqrt{s} = 13 TeV", lumi );

//    TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
//    label_top->SetBorderSize(0);
//    label_top->SetFillColor(kWhite);
//    label_top->SetTextSize(0.038);
//    label_top->SetTextAlign(31); // align right
//    label_top->SetTextFont(62);
//    label_top->AddText(text.c_str());
//    label_top->Draw("same");
    
    TLegend* legend = new TLegend( 0.6, 0.9-(pdfSets.size()+1)*0.04, 0.93, 0.9 );
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    
    legend->AddEntry( hNNPDF[r+3], "TT+jets LO", "PL");
    legend->AddEntry( hNNPDF[r+2], pdfSets[0].c_str(), "PL");
    legend->AddEntry( hMMHT[r+2], pdfSets[1].c_str(), "PL");
    legend->AddEntry( hCT10[r+2], pdfSets[2].c_str(), "PL");
    
    legend->Draw("same");
    
    hNNPDF[r+3]->Draw("same");
    gNNPDF->Draw("P,same");
    gMMHT->Draw("P,same");
    gCT10->Draw("P,same");
    
    gPad->RedrawAxis();
    
    
    c1->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.1);
    pad2->Draw();
    pad2->cd();
    
    pad2->SetGridy();
    
    hPDF->SetStats(0);
    hPDF->GetXaxis()->SetLabelSize(0.00);
    hPDF->GetXaxis()->SetTickLength(0.09);
    hPDF->GetYaxis()->SetNdivisions(5,5,0);
    hPDF->GetYaxis()->SetRangeUser(0.5,1.5);
    hPDF->GetYaxis()->SetTitleSize(0.17);
    hPDF->GetYaxis()->SetTitleOffset(0.4);
    hPDF->GetYaxis()->SetLabelSize(0.17);
    hPDF->GetYaxis()->SetTitle("ratio");
     
    hPDF->Draw("P");
    
    gPad->RedrawAxis();
    
    c1->cd();
    
    c1->SaveAs( Form("%s/TopoRegion_%.0f.pdf", outputdir.c_str(), r/4.) );
    c1->SaveAs( Form("%s/TopoRegion_%.0f.eps", outputdir.c_str(), r/4.) );
    c1->SaveAs( Form("%s/TopoRegion_%.0f.png", outputdir.c_str(), r/4.) );
    
    delete c1;
    delete h_axes;
    delete gNNPDF;
    delete gCT10;
    delete gMMHT;
    delete hPDF;

    r+=4;
    
  }
  
}
