#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TApplication.h"



TH1D* getHisto( TTree* tree, const std::string& name, bool draw=false );
void getPfromHisto( TH1D* h1, float& p, float& p_err );
void getPfromFunc( TF1* line, float& p, float& p_err );
TH1D* fillFromTree( const std::string& name, TTree* tree, float p );
float getCorrection( int njets, float p );
void setStyle();



int nBins = 50;
float xMin = 0.;
float xMax = 1500.;



int main( int argc, char* argv[] ) {

  setStyle();


  TFile* file = TFile::Open("toyMC.root", "recreate");
  file->cd();

  TH1D* h1_pull = new TH1D( "pull", "", 100, -3., 3. );
  

  float p = 0.03;


  TTree* tree;

  int nToys = 500;

  for( int iToy=0; iToy<nToys+1; ++iToy ) {

    std::cout << "-> Starting toy n." << iToy << std::endl;

    if( tree!=0 && iToy!=0 )
      delete tree;

    std::string treeName = (iToy==nToys) ? "tree" : "tmp_tree";
    tree = new TTree(treeName.c_str(), "");


    int njets;
    tree->Branch( "nJets", &njets );
    int nbjets;
    tree->Branch( "nBJets", &nbjets );
    float mt2;
    tree->Branch( "mt2", &mt2 );
    


    TRandom3 rand(iToy);

    int nEntries = (iToy==nToys) ? 100000 : 10000;

    for( unsigned i=0; i<nEntries; ++i ) {

      njets = rand.Poisson(6);

      if( njets<2 ) continue;

      nbjets = rand.Binomial( njets, p );

      mt2 = rand.Exp(200.);

      tree->Fill();

    }

    if( iToy!=nToys ) {

      TH1D* h1_p = getHisto( tree, "toy_tmp" );
      float p_est, p_est_err;
      getPfromHisto( h1_p, p_est, p_est_err );

      h1_pull->Fill( (p_est-p)/p_est_err );

      delete h1_p;

    }


  } // for iToy

  std::cout << "   Done." << std::endl << std::endl;

  std::cout << "-> Starting computation..." << std::endl;
  
  TH1D* h1_p = getHisto( tree, "toy", true );
  float p_est, p_est_err;
  getPfromHisto( h1_p, p_est, p_est_err );


  TH1D* h1_mcTruth = new TH1D("mcTruth", "", nBins, xMin, xMax );
  h1_mcTruth->Sumw2();
  tree->Project( "mcTruth", "mt2", "nBJets==2" );
  TH1D* h1_closure   = fillFromTree( "closure", tree, p_est );
  TH1D* h1_closureMC = fillFromTree( "closureMC", tree, p );


  std::cout << "   Done." << std::endl << std::endl;


  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();
  c1->SetLogy();


  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0.1, 500. );
  h2_axes->SetXTitle( "Dummy M_{T2}" );
  h2_axes->SetYTitle( "Events" );
  h2_axes->Draw();

  h1_mcTruth->SetMarkerStyle(20);
  h1_mcTruth->SetMarkerSize(1.1);
  h1_mcTruth->SetMarkerColor(kBlack);
  h1_mcTruth->SetLineColor(kBlack);

  h1_closure->SetMarkerStyle(25);
  h1_closure->SetMarkerSize(1.1);
  h1_closure->SetMarkerColor(46);
  h1_closure->SetLineColor(46);

  h1_closureMC->SetMarkerStyle(24);
  h1_closureMC->SetMarkerSize(1.1);
  h1_closureMC->SetMarkerColor(kBlack);
  h1_closureMC->SetLineColor(kBlack);

  h1_mcTruth->Draw("P same");
  //h1_closure->Draw("P same");
  h1_closureMC->Draw("P same");


  TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top->SetBorderSize(0);
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(62);
  label_top->AddText("Toy Monte Carlo");
  label_top->Draw("same");

  TLegend* legend = new TLegend( 0.53, 0.73, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( h1_mcTruth, "MC Truth b=2", "P" );
  legend->AddEntry( h1_closureMC, "Extrap. from b=0", "P" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs("closureToy.eps");

  TCanvas* c2 = new TCanvas( "c2", "", 600, 600 );
  c2->cd();

  TH2D* h2_axes2 = new TH2D("axes2", "", 10, xMin, xMax, 10, 0., 2 );
  h2_axes2->SetXTitle( "Dummy M_{T2}" );
  h2_axes2->SetYTitle( "Ratio to MC Truth" );
  h2_axes2->Draw();

  TH1D* h1_ratio = (TH1D*)h1_closure->Clone("ratio");
  h1_ratio->Divide( h1_mcTruth );
  TH1D* h1_ratioMC = (TH1D*)h1_closureMC->Clone("ratioMC");
  h1_ratioMC->Divide( h1_mcTruth );

  h1_ratio->SetMarkerStyle(21);
  h1_ratio->SetMarkerSize(1.1);
  h1_ratio->SetMarkerColor(46);
  h1_ratio->SetLineColor(46);

  h1_ratioMC->SetMarkerStyle(20);
  h1_ratioMC->SetMarkerSize(1.1);
  h1_ratioMC->SetMarkerColor(kBlack);
  h1_ratioMC->SetLineColor(kBlack);

  TLine* lineOne = new TLine( xMin, 1., xMax, 1. );
  lineOne->Draw("Same");

  h1_ratio->Draw("P same");
  h1_ratioMC->Draw("P same");

  label_top->Draw("same");

  TLegend* legend2 = new TLegend( 0.15, 0.2, 0.55, 0.3 );
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry( h1_ratioMC, "Extrap. from b=0", "P" );
  //legend2->AddEntry( h1_ratio, "Estimate from b=0 (p_{est})", "P" );
  legend2->Draw("same");


  gPad->RedrawAxis();

  c2->SaveAs("closureToy_ratio.eps");
  


  std::cout << "-> Writing file..." << std::endl;
  file->cd();
  tree->Write();
  h1_p->Write();
  h1_mcTruth->Write();
  h1_closure->Write();
  h1_closureMC->Write();
  h1_pull->Write();
  file->Close();


  return 0;

}








TH1D* getHisto( TTree* tree, const std::string& name, bool draw ) {


  float njetMin = 2;
  float njetMax = 12;
  int nbins_ratio = njetMax-njetMin + 1;
  TH1D* histo = new TH1D(Form("histo_%s", name.c_str()), "", nbins_ratio, njetMin-0.5, njetMax+0.5 );
  

  float yMin = 0.;
  float yMax = 0.4;

  

  for( int njet=(int)njetMin; njet<=(int)njetMax; ++njet ) {

    std::string name_0b(Form("%s_%dj_0b", name.c_str(), njet));
    std::string name_1b(Form("%s_%dj_1b", name.c_str(), njet));

    TH1D* h1_0b = new TH1D( name_0b.c_str(), "", nBins, 0., 20000. );
    TH1D* h1_1b = new TH1D( name_1b.c_str(), "", nBins, 0., 20000. );

    h1_0b->Sumw2();
    h1_1b->Sumw2();

    tree->Project( name_0b.c_str(), "mt2", Form("nJets==%d && nBJets==0", njet) );
    tree->Project( name_1b.c_str(), "mt2", Form("nJets==%d && nBJets==1", njet) );

    Double_t int_0b_err;
    Double_t int_0b = h1_0b->IntegralAndError(1, nBins, int_0b_err);
    Double_t int_1b_err;
    Double_t int_1b = h1_1b->IntegralAndError(1, nBins, int_1b_err);

    if( int_1b==0 ) continue;

    float integralRatio = int_1b/int_0b;
    float integralRatioErr = (1./int_0b)*(1./int_0b)*int_1b_err*int_1b_err + int_1b*int_1b/(int_0b*int_0b*int_0b*int_0b)*int_0b_err*int_0b_err;
    integralRatioErr = sqrt(integralRatioErr);

    integralRatio*=0.8;
    integralRatioErr*=0.8;

    int bin = histo->FindBin(njet);
    float m = integralRatio/(float)njet;
    float m_err = integralRatioErr/(float)njet;
    float p = m/(1.+m);
    float p_err = 1./( (1.+m)*(1.+m) ) * m_err;


    p *= 100.; // in percent
    p_err *= 100.; // in percent

    histo->SetBinContent( bin, integralRatio );
    histo->SetBinError( bin, integralRatioErr );

    if( integralRatio+integralRatioErr>yMax ) yMax = integralRatio+integralRatioErr;

    delete h1_0b;
    delete h1_1b;

  }

  yMax *= 1.1;

  if( yMax > 1. ) yMax = 1.;


  TF1* line = new TF1("line", "[0] + [1]*x", njetMin-0.5, njetMax+0.5 );
  line->SetLineColor(kRed);
  histo->Fit( line, "RQ" );

  float p, p_err;
  getPfromFunc( line, p, p_err );


  if( draw ) {

    TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
    c1->cd();

    TH2D* h2_axes = new TH2D("axes", "", 10, njetMin-0.5, njetMax+0.5, 10, yMin, yMax );
    h2_axes->SetXTitle("Number of Jets");
    h2_axes->SetYTitle("Events(b=1) / Events(b=0)");
    h2_axes->Draw();

    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(1.6);

    histo->Draw("P same");

   
    TPaveText* pLabel = new TPaveText( 0.2, 0.7, 0.55, 0.8, "brNDC" );
    pLabel->SetFillColor(0);
    pLabel->SetTextSize(0.038);
    pLabel->AddText( Form("p = %.2f #pm %.2f %%", p*100. , p_err*100.) );
    pLabel->Draw("same");


    gPad->RedrawAxis();


    c1->SaveAs("p_toy.eps");

    delete c1;
    delete h2_axes;

  } // if draw

  histo->Scale(1./100.);

  return histo;

}


void getPfromHisto( TH1D* h1, float& p, float& p_err ) {

  getPfromFunc( h1->GetFunction("line"), p, p_err );

}
  


void getPfromFunc( TF1* line, float& p, float& p_err ) {

  float m = line->GetParameter(1);
  float m_err = line->GetParError(1);

  p = m/(1.+m);
  p_err = 1./( (1.+m)*(1.+m) ) * m_err;

}
 



TH1D* fillFromTree( const std::string& name, TTree* tree, float p ) {


  TH1D* h1_histo = new TH1D( name.c_str(), "", nBins, xMin, xMax );
  h1_histo->Sumw2();

  int njets;
  tree->SetBranchAddress( "nJets", &njets );
  int nbjets;
  tree->SetBranchAddress( "nBJets", &nbjets );
  float mt2;
  tree->SetBranchAddress( "mt2", &mt2 );

  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( nbjets==0 ) {
      float corr_mc = getCorrection( njets, p );
      h1_histo->Fill( mt2, corr_mc );
    }

  }

  return h1_histo;

}




float getCorrection( int njets, float p ) {

  float corr = TMath::Binomial(njets,2) * p*p / ( (1.-p)*(1.-p) ); //TMath::Power( 1.-p, njets-2 ); 
  //float corr = (float)njets*(njets-1.)*p*p/((1.-p)*(1.-p)); 
  //float corr = 0.5*(float)njets*(njets-1.)*p*p/((1.-p)*(1.-p)); 

  return corr;

}





void setStyle() {

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

}
