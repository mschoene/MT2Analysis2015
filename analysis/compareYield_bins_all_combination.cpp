#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2DrawTools.h"
#include "interface/MT2Config.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>


#include "TMath.h"
#include "TTreeFormula.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"

#include "RooHistError.h"


float lumi2015; //fb-1
float lumi2016; //fb-1 

void drawYields( const std::string& outputdir, std::string dir2015, std::string dir2016 );


int main( int argc, char* argv[] ) {
  
  
  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeLostLepton                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;
  
  
  if( argc!=3 ) {
    std::cout << "USAGE: ./computeLostLepton [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  


  std::string configFileName2015(argv[1]);
  MT2Config cfg2015(configFileName2015);

  std::string configFileName2016(argv[2]);
  MT2Config cfg2016(configFileName2016);

  lumi2015 = cfg2015.lumi();
  lumi2016 = cfg2016.lumi();

  TH1::AddDirectory(kTRUE);

  std::string dir2015 = cfg2015.getEventYieldDir();
  std::string dir2016 = cfg2016.getEventYieldDir();
  std::string outputdir = "./YieldComparison_dataMC_binned_combination2015plus2016/";

  drawYields( outputdir.c_str(), dir2015, dir2016 );

  return 0;

}

void drawYields( const std::string& outputdir, std::string dir2015, std::string dir2016 ) {
  
  MT2DrawTools::setStyle();

  system(Form("mkdir -p %s", outputdir.c_str()));

  std::vector<int> colors;
//  colors.push_back( 402 );
//  colors.push_back( 430 );
//  colors.push_back( 418 );
  
//  colors.push_back( kYellow+1 );
//  colors.push_back( kAzure+5 );
//  colors.push_back( kGreen+3 );
  
  colors.push_back( kYellow+1 );
  colors.push_back( kAzure+4 );
  colors.push_back( kGreen+2 );

//  colors.push_back( kRed+2 );
//  colors.push_back( kAzure+5 );
//  colors.push_back( kGreen-7 );
  
  unsigned int bgSize = 3;
  
  TFile* file2015 = TFile::Open(Form("%s/YieldComparison_dataMC_binned/histograms_ALL.root", dir2015.c_str()));
  TFile* file2016 = TFile::Open(Form("%s/YieldComparison_dataMC_binned/histograms_ALL.root", dir2016.c_str()));
  TH1D* hdata2015 = (TH1D*) file2015->Get("hdata");
  TH1D* hdata2016 = (TH1D*) file2016->Get("hdata");
  TH1D* hdata = (TH1D*) hdata2015->Clone();
  hdata->Add(hdata2016);
  hdata->GetYaxis()->SetTitle("Entries");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( 1 );
  hdata->SetMarkerColor( 1 );

  TH1D* hestimate_all2015 = (TH1D*) file2015->Get("hestimate_all");
  TH1D* hestimate_all2016 = (TH1D*) file2016->Get("hestimate_all");
  TH1D* hestimate_all = (TH1D*) hestimate_all2015->Clone();
  hestimate_all->Add(hestimate_all2016);

  TH1D* hestimate_all_forRatio2015 = (TH1D*) file2015->Get("hestimate_all_forRatio");
  TH1D* hestimate_all_forRatio2016 = (TH1D*) file2016->Get("hestimate_all_forRatio");
  TH1D* hestimate_all_forRatio = (TH1D*) hestimate_all_forRatio2015->Clone();
  hestimate_all_forRatio->Add(hestimate_all_forRatio2016);

  TH1D* hestimate[bgSize];
  TH1D* hestimate2015[bgSize];
  TH1D* hestimate2016[bgSize];
  for(unsigned int b=0; b<bgSize; ++b){

    hestimate2015[b] = (TH1D*) file2015->Get(Form("hestimate_%d", b));
    hestimate2016[b] = (TH1D*) file2016->Get(Form("hestimate_%d", b));

    hestimate[b] = (TH1D*) hestimate2015[b]->Clone();
    hestimate[b]->Add(hestimate2016[b]);

    hestimate[b]->GetYaxis()->SetTitle("Entries");
    hestimate[b]->SetFillColor(colors[b]);
    hestimate[b]->SetLineColor(1);

  }


  THStack bgStack("bgStack", "");

  TH1D* hPull = new TH1D("hPull", "", 101, -5.05, 5.05);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Est. - Data)/#sigma");
  hPull->GetYaxis()->SetTitle("Entries");
  
  std::string fullPath = outputdir;
  
  std::string labelsMono[12]={"[200,250]","[250,350]","[350,450]","[450,575]","[575,700]","[700,1000]",">1000", "[200,250]","[250,350]","[350,450]","[450,575]",">575"};
  

  for(unsigned int b=0; b<bgSize; ++b){

    hestimate[b]->SetLineWidth(0);
    bgStack.Add(hestimate[b]);
    //bgStack.Add(hestimate_forRatio[b]);

  }


  TGraphAsymmErrors* g_Ratio = MT2DrawTools::getRatioGraph(hdata, hestimate_all_forRatio, "binWidth");  
  g_Ratio->SetMarkerStyle(20);
  g_Ratio->SetMarkerSize(1.6);
  g_Ratio->SetMarkerColor( 1 );
  g_Ratio->SetLineColor(1);
  g_Ratio->GetXaxis()->SetLabelSize(0.00);
  g_Ratio->GetXaxis()->SetTickLength(0.09);
  g_Ratio->GetYaxis()->SetNdivisions(5,5,0);
  g_Ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  g_Ratio->GetYaxis()->SetTitleSize(0.17);
  g_Ratio->GetYaxis()->SetTitleOffset(0.4);
  g_Ratio->GetYaxis()->SetLabelSize(0.17);
  g_Ratio->GetYaxis()->SetTitle("Ratio");

  TGraphAsymmErrors* g_Ratio_zero = new TGraphAsymmErrors(*(g_Ratio));
  g_Ratio_zero->SetMarkerSize(0);
  g_Ratio_zero->SetLineColor( 1 );
  g_Ratio_zero->SetMarkerColor( 1 );
  
  TGraphAsymmErrors* gdata = MT2DrawTools::getPoissonGraph(hdata, true, "binWidth");
  gdata->GetYaxis()->SetTitle("Entries");
  gdata->SetMarkerStyle(20); 
  gdata->SetMarkerSize(1.6);
  gdata->SetLineColor( 1 );
  gdata->SetMarkerColor( 1 );

  TGraphAsymmErrors* gdata_zero = new TGraphAsymmErrors(*(gdata));
  gdata_zero->SetMarkerSize(0);
  gdata_zero->SetLineColor( 1 );
  gdata_zero->SetMarkerColor( 1 );

  
  for(int iBin=1; iBin<=hestimate_all->GetNbinsX(); ++iBin){
    
    float thisData    = hdata->GetBinContent(iBin);
    float thisDataErr = gdata->GetErrorY(iBin-1);
    std::cout << "TEST!" << std::endl;
    std::cout << thisData << "\t" << gdata->GetY()[iBin-1] << "\t" << thisDataErr << std::endl;
    
    float thisEst     = hestimate_all->GetBinContent(iBin);
    float thisEstErr  = hestimate_all->GetBinError(iBin);
    
    
    hPull->Fill( (thisEst-thisData)/( TMath::Sqrt( thisDataErr*thisDataErr + thisEstErr*thisEstErr ) ) );
    
  }


  //  gdata->GetXaxis()->SetRangeUser(0, thisBin);

  TCanvas* c2 = new TCanvas("c2", "", 1100, 600);
  c2->cd();
  
  std::string thisName = Form("%s_ratio", hdata->GetName());
  TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
  h_Ratio->Divide(hestimate_all_forRatio);
  h_Ratio->Write();
  h_Ratio->SetStats(0);
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->GetXaxis()->SetLabelSize(0.00);
  h_Ratio->GetXaxis()->SetTickLength(0.09);
  h_Ratio->GetYaxis()->SetNdivisions(5,5,0);
  h_Ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  h_Ratio->GetYaxis()->SetTitleSize(0.17);
  h_Ratio->GetYaxis()->SetTitleOffset(0.4);
  h_Ratio->GetYaxis()->SetLabelSize(0.17);
  h_Ratio->GetYaxis()->SetTitle("Ratio");

  

  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  pad1->SetLogy();
  
  float yMax_1 = hdata->GetMaximum()*1.5;
  float yMax_2 = 1.2*(hdata->GetMaximum() + hdata->GetBinError(hestimate_all->GetMaximumBin()));
  float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
  float yMax_3 = hestimate_all->GetMaximum()*1.5;
  float yMax_4 = 1.2*(hestimate_all->GetMaximum() + hestimate_all->GetBinError(hestimate_all->GetMaximumBin()));
  float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
  float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
  
  float yMin = 1e-3;
  //  yMin=0;
  yMax*=20.;
  
  int thisBin=174;

  hestimate_all->GetXaxis()->SetRangeUser(0, thisBin);
  hdata->GetXaxis()->SetRangeUser(0, thisBin);
  gdata->GetXaxis()->SetRangeUser(0, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(0, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->GetXaxis()->SetLabelFont(62);
  hestimate_all->SetFillStyle(3244);
  hestimate_all->SetFillColor(kGray+2);

//  TGraphAsymmErrors* g_data = new TGraphAsymmErrors(0);
//  for( int iBin=1; iBin<(hdata->GetXaxis()->GetNbins()+1); ++iBin ) {
//
//    double y;
//    double x, xerr;
//
//    x = hdata->GetBinCenter(iBin);
//    xerr = hdata->GetBinWidth(iBin)/2.;
//
//    y = hdata->GetBinContent(iBin);
//    double yerr = hdata->GetBinError(iBin);
//
//    int thisPoint = g_data->GetN();
//    g_data->SetPoint( thisPoint, x, y );
//    g_data->SetPointError( thisPoint, xerr, xerr, yerr, yerr );
//
//  }
  

//  hdata->Draw("pe");
//  bgStack.Draw("histo, same");
//  hestimate_all->Draw("E2,same");
//  hdata->Draw("pe,same");

  gdata->Draw("pe");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");
  
  TH1D* prefit=new TH1D("prefit", "", 1, 0, 1);
  prefit->SetFillColor(0);
  prefit->SetLineColor(0);

  TLegend* legend = new TLegend( 0.8, 0.9-(bgSize+1-1)*0.06-0.06+0.02+0.02, 0.93, 0.9-0.06+0.02+0.02+0.02 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( hdata, "Data", "PL" );
  //  legend->AddEntry( prefit, "Pre-fit SM", "F");
  legend->AddEntry( hestimate[0], "Multijet", "F");
  legend->AddEntry( hestimate[1], "Lost lepton", "F");
  legend->AddEntry( hestimate[2], "Z #rightarrow #nu#bar{#nu}", "F");

  legend->Draw("same");

  //  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi2015,lumi2016);
  labelTop->Draw("same");
  TPaveText* labelCMS = MT2DrawTools::getLabelCMS();
  labelCMS->Draw("same");
  
  int nHTRegions = 6;
  std::vector< std::string > htRegions;
  htRegions.push_back("1 Jet");
  htRegions.push_back("H_{T} [200, 450] GeV");
  htRegions.push_back("H_{T} [450, 575] GeV");
  htRegions.push_back("H_{T} [575, 1000] GeV");
  htRegions.push_back("H_{T} [1000, 1500] GeV");
  htRegions.push_back("H_{T} > 1500 GeV");
//  htRegions.push_back("#it{H}_{T} [200, 450] GeV");
//  htRegions.push_back("#it{H}_{T} [450, 575] GeV");
//  htRegions.push_back("#it{H}_{T} [575, 1000] GeV");
//  htRegions.push_back("#it{H}_{T} [1000, 1500] GeV");
//  htRegions.push_back("#it{H}_{T} > 1500 GeV");
////  htRegions.push_back("very low H_{T}");
////  htRegions.push_back("1 Jet");
////  htRegions.push_back("low H_{T}");
////  htRegions.push_back("medium H_{T}");
////  htRegions.push_back("high H_{T}");
////  htRegions.push_back("extreme H_{T}");
  
  TPaveText* htBox[5];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){
    htBox[iHT] = new TPaveText(0.4, 0.9-0.06+0.02, 0.7, 0.85+0.02, "brNDC");
    htBox[iHT]->AddText( htRegions[iHT].c_str() );
    
    htBox[iHT]->SetBorderSize(0);
    htBox[iHT]->SetFillColor(kWhite);
    htBox[iHT]->SetTextSize(0.038);
    htBox[iHT]->SetTextAlign(21); // align centered
    htBox[iHT]->SetTextFont(62);
    //    htBox[iHT]->Draw("same");

  }
  //  htBox[0]->Draw("same");

  gPad->RedrawAxis();
  
  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();  

  TH2D* h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0., 2.0 );
  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.4);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitle("Ratio");
  
  TLine* LineCentral = new TLine(0, 1.0, thisBin, 1.0);
  LineCentral->SetLineColor(1);

  
  std::string thisName_Band =  Form("%s_band", hestimate_all->GetName());
  TH1D* h_band = (TH1D*)hestimate_all->Clone(thisName_Band.c_str());
  h_band->SetMarkerSize(0);
  h_band->SetFillColor (kGray+2);
  h_band->SetFillStyle (3244);
  for ( int iBin=1; iBin <= hestimate_all->GetNbinsX(); iBin++){
    
    h_band->SetBinContent(iBin,1);

    double error=0;

    if(hestimate_all_forRatio->GetBinContent(iBin)>0)
      error = hestimate_all->GetBinError(iBin)/hestimate_all->GetBinContent(iBin);
    //else error = hestimate_all->GetBinError(iBin);

    h_band->SetBinError(iBin, error);

  }


  h2_axes_ratio->Draw("");
  h_band->Draw("E2same");
  LineCentral->Draw("same");
  h_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2->cd();
//  //  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.pdf", fullPath.c_str()) );
//  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.png", fullPath.c_str()) );
//  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.eps", fullPath.c_str()) );

  TCanvas* c2_0 = new TCanvas("c2_0", "", 1100, 600);
  c2_0->cd();
  
  TPad *pad1_0 = new TPad("pad1_0","pad1_0",0,0.3-0.1,1,1);
  pad1_0->SetBottomMargin(0.15);
  pad1_0->Draw();
  pad1_0->cd();


  int oldBin=0;
  pad1_0->SetLogy();
    
  thisBin=12;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  //  gdata->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);

  hestimate_all->GetYaxis()->SetTitleOffset(0.95);
  hestimate_all->GetYaxis()->SetLabelSize(0.042);

  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //hdata->Draw("pe,same");
  gdata_zero->Draw("pe1,same");
  gdata->Draw("pe1,same");

  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");
  
  htBox[0]->Draw("same");

  gPad->RedrawAxis();
  
  bool doLogRatio=true;

  c2_0->cd();
  TPad *pad2_0 = new TPad("pad2_0","pad2_0",0,0,1,0.21);
  pad2_0->SetTopMargin(0.05);
  pad2_0->SetBottomMargin(0.1);
  pad2_0->Draw();
  pad2_0->cd();

  TH2D* h2_axes_ratio_0;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
    h_Ratio->GetYaxis()->SetRangeUser(0.1, 10.0);
  }
  else
    h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0., 2.0 );



  h2_axes_ratio_0->SetStats(0);
  h2_axes_ratio_0->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_0->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_0->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_0->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_0->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_0->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_0->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_0 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_0->SetLineColor(1);

  h2_axes_ratio_0->Draw("");
  h_band->Draw("E2same");
  LineCentral_0->Draw("same");
  //h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_0->cd();
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.pdf", fullPath.c_str()) );
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.png", fullPath.c_str()) );
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.eps", fullPath.c_str()) );



  TCanvas* c2_1 = new TCanvas("c2_1", "", 1100, 600);
  c2_1->cd();
  
  TPad *pad1_1 = new TPad("pad1_1","pad1_1",0,0.3-0.1,1,1);
  pad1_1->SetBottomMargin(0.15);
  pad1_1->Draw();
  pad1_1->cd();

  pad1_1->SetLogy();
    
  oldBin=thisBin;
  thisBin=36;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //hdata->Draw("pe,same");
  gdata_zero->Draw("pe1,same");
  gdata->Draw("pe1,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[1]->Draw("same");

  gPad->RedrawAxis();
  
  c2_1->cd();
  TPad *pad2_1 = new TPad("pad2_1","pad2_1",0,0,1,0.21);
  pad2_1->SetTopMargin(0.05);
  pad2_1->SetBottomMargin(0.1);
  pad2_1->Draw();
  pad2_1->cd();

  TH2D* h2_axes_ratio_1;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    //    h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0., 7.5 );
    h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  h2_axes_ratio_1->SetStats(0);
  h2_axes_ratio_1->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_1->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_1->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_1->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_1->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_1->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_1->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_1 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_1->SetLineColor(1);

  h2_axes_ratio_1->Draw("");
  h_band->Draw("E2same");
  LineCentral_1->Draw("same");
  //h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");

  gPad->RedrawAxis();

  c2_1->cd();
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.png", fullPath.c_str()) );
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.eps", fullPath.c_str()) );



  TCanvas* c2_2 = new TCanvas("c2_2", "", 1100, 600);
  c2_2->cd();
  
  TPad *pad1_2 = new TPad("pad1_2","pad1_2",0,0.3-0.1,1,1);
  pad1_2->SetBottomMargin(0.15);
  pad1_2->Draw();
  pad1_2->cd();

  pad1_2->SetLogy();
    
  oldBin=thisBin;
  thisBin=67;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //hdata->Draw("pe,same");
  gdata_zero->Draw("pe1,same");
  gdata->Draw("pe1,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");
    
  htBox[2]->Draw("same");

  gPad->RedrawAxis();
  
  c2_2->cd();
  TPad *pad2_2 = new TPad("pad2_2","pad2_2",0,0,1,0.21);
  pad2_2->SetTopMargin(0.05);
  pad2_2->SetBottomMargin(0.1);
  pad2_2->Draw();
  pad2_2->cd();

  TH2D* h2_axes_ratio_2;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  h2_axes_ratio_2->SetStats(0);
  h2_axes_ratio_2->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_2->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_2->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_2->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_2->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_2->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_2->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_2 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_2->SetLineColor(1);

  h2_axes_ratio_2->Draw("");
  h_band->Draw("E2same");
  LineCentral_2->Draw("same");
  //h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_2->cd();
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.png", fullPath.c_str()) );
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.eps", fullPath.c_str()) );




  TCanvas* c2_3 = new TCanvas("c2_3", "", 1100, 600);
  c2_3->cd();
  
  TPad *pad1_3 = new TPad("pad1_3","pad1_3",0,0.3-0.1,1,1);
  pad1_3->SetBottomMargin(0.15);
  pad1_3->Draw();
  pad1_3->cd();

  pad1_3->SetLogy();
    
  oldBin=thisBin;
  thisBin=109;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //hdata->Draw("pe,same");
  gdata_zero->Draw("pe1,same");
  gdata->Draw("pe1,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[3]->Draw("same");

  gPad->RedrawAxis();
  
  c2_3->cd();
  TPad *pad2_3 = new TPad("pad2_3","pad2_3",0,0,1,0.21);
  pad2_3->SetTopMargin(0.05);
  pad2_3->SetBottomMargin(0.1);
  pad2_3->Draw();
  pad2_3->cd();

  TH2D* h2_axes_ratio_3;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    //    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0., 4.0 );
    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  h2_axes_ratio_3->SetStats(0);
  h2_axes_ratio_3->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_3->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_3->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_3->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_3->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_3->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_3->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_3 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_3->SetLineColor(1);

  h2_axes_ratio_3->Draw("");
  h_band->Draw("E2same");
  LineCentral_3->Draw("same");
  //h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_3->cd();
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.png", fullPath.c_str()) );
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.eps", fullPath.c_str()) );




  TCanvas* c2_4 = new TCanvas("c2_4", "", 1100, 600);
  c2_4->cd();
  
  TPad *pad1_4 = new TPad("pad1_4","pad1_4",0,0.3-0.1,1,1);
  pad1_4->SetBottomMargin(0.15);
  pad1_4->Draw();
  pad1_4->cd();

  pad1_4->SetLogy();
    
  oldBin=thisBin;
  thisBin=144;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //hdata->Draw("pe,same");
  gdata_zero->Draw("pe1,same");
  gdata->Draw("pe1,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[4]->Draw("same");

  gPad->RedrawAxis();
  
  c2_4->cd();
  TPad *pad2_4 = new TPad("pad2_4","pad2_4",0,0,1,0.21);
  pad2_4->SetTopMargin(0.05);
  pad2_4->SetBottomMargin(0.1);
  pad2_4->Draw();
  pad2_4->cd();

  TH2D* h2_axes_ratio_4;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    //    h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0., 6.0 );
    h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0., 2.0 );


  h2_axes_ratio_4->SetStats(0);
  h2_axes_ratio_4->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_4->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_4->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_4->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_4->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_4->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_4->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_4 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_4->SetLineColor(1);

  h2_axes_ratio_4->Draw("");
  h_band->Draw("E2same");
  LineCentral_4->Draw("same");
  //h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_4->cd();
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.png", fullPath.c_str()) );
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.eps", fullPath.c_str()) );




  TCanvas* c2_5 = new TCanvas("c2_5", "", 1100, 600);
  c2_5->cd();
  
  TPad *pad1_5 = new TPad("pad1_5","pad1_5",0,0.3-0.1,1,1);
  pad1_5->SetBottomMargin(0.15);
  pad1_5->Draw();
  pad1_5->cd();

  pad1_5->SetLogy();
    
  oldBin=thisBin;
  thisBin=174;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //hdata->Draw("pe,same");
  gdata_zero->Draw("pe1,same");
  gdata->Draw("pe1,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[5]->Draw("same");

  gPad->RedrawAxis();
  
  c2_5->cd();
  TPad *pad2_5 = new TPad("pad2_5","pad2_5",0,0,1,0.21);
  pad2_5->SetTopMargin(0.05);
  pad2_5->SetBottomMargin(0.1);
  pad2_5->Draw();
  pad2_5->cd();

  TH2D* h2_axes_ratio_5;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    //    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0., 6.0 );
    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  h2_axes_ratio_5->SetStats(0);
  h2_axes_ratio_5->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_5->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_5->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_5->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_5->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_5->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_5->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_5 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_5->SetLineColor(1);

  h2_axes_ratio_5->Draw("");
  h_band->Draw("E2same");
  LineCentral_5->Draw("same");
  //  h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_5->cd();
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.C", fullPath.c_str()) );
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.png", fullPath.c_str()) );
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.eps", fullPath.c_str()) );


  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);
  
  TF1* fgauss= new TF1("fgauss", "gaus", -5, 5);
  fgauss->SetLineColor(2);

  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  hPull->SetStats(1110);
  hPull->Draw("hist");
  hPull->Fit("fgauss");
  fgauss->Draw("same");
  c3->SaveAs( Form("%s/PullDistribution.pdf", fullPath.c_str()) );
  c3->SaveAs( Form("%s/PullDistribution.png", fullPath.c_str()) );
  
}

