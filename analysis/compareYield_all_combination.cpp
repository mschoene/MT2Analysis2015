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
#include "TRandom3.h"
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
  std::string outputdir = "./YieldComparison_dataMC_combination2015plus2016/";
  
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
  colors.push_back( kYellow+1 );
  colors.push_back( kAzure+4 );
  colors.push_back( kGreen+2 );
  
  unsigned int bgSize = 3;
  
  TFile* file2015 = TFile::Open(Form("%s/YieldComparison_dataMC/histograms_ALL.root", dir2015.c_str()));
  TFile* file2016 = TFile::Open(Form("%s/YieldComparison_dataMC/histograms_ALL.root", dir2016.c_str()));
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

  TH1D* hPull = new TH1D("hPull", "", 21, -5.25, 5.25);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Est._{H_{T}, N_{j}, N_{b}}^{pre-fit} - Data) / #sqrt{#sigma_{Est.}^{ 2} + #sigma_{Data}^{ 2}}");
  hPull->GetYaxis()->SetTitle("Entries");

  TH1D* hPvalue = new TH1D("hPvalue", "", 14, 0, 1.05);
  hPvalue->Sumw2();
  hPvalue->GetXaxis()->SetTitle("p-value");
  hPvalue->GetYaxis()->SetTitle("Entries");
  
  std::string fullPath = outputdir;
  
  std::string labelsMono[12]={"[200,250]","[250,350]","[350,450]","[450,575]","[575,700]","[700,1000]",">1000", "[200,250]","[250,350]","[350,450]","[450,575]",">575"};

  for(unsigned int b=0; b<bgSize; ++b){

    bgStack.Add(hestimate[b]);
    //bgStack.Add(hestimate_forRatio[b]);

  }

  for(int iBin=1; iBin<=hestimate_all->GetNbinsX(); ++iBin){
    

      float thisData     = hdata->GetBinContent(iBin);

      float thisEst     = hestimate_all->GetBinContent(iBin);
      float thisEstErr  = hestimate_all->GetBinError(iBin);

      int obs = thisData;
      float meanExp = thisEst;
      float meanUnc = thisEstErr;    
      int N=100000;
      
      TRandom3 randGen;
      
      unsigned int counterN(0);
      unsigned int counterD(0);
      
      bool doLeft = obs<=meanExp;
      
      for(int i=0; i<N; i++){
	double meanShift = randGen.Gaus(0.,(double)meanUnc);
	//cout << "meanShift: " << meanShift << endl;                                                                                                                                    
	double rand = randGen.Poisson(double(meanExp)+meanShift); counterD++;
	
	if(doLeft)
	  {
	    if(rand<=obs) counterN++;
	  }
	else{
	  if(rand>=obs) counterN++;
	}
	
	//cout << "rand " << i << " : " << rand.Poisson(meanExp) << endl;                                                                                                                
      }
      
      double prob=1.0*counterN/counterD;
      double significance  = TMath::NormQuantile(1-prob);
      
//      std::cout << "probability: " << prob  << std::endl;
//      std::cout << "significance: " << significance << std::endl;
//      
//      std::cout << "Bin: " << iBin << std::endl;
//      std::cout << "Obs: " << obs << std::endl;
//      std::cout << "exp: " << meanExp << "\t" << meanUnc << std::endl;
      
      hPvalue->Fill(prob);
      
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

  TCanvas* c2 = new TCanvas("c2", "", 1100, 600);
  c2->cd();
 
  //  c2->SetLeftMargin(0.); 

  std::string thisName = Form("%s_ratio", hdata->GetName());
  TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
  h_Ratio->Divide(hestimate_all_forRatio);
  //h_Ratio->Divide(hestimate_all);
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
  
  int thisBin=67;
  
  hestimate_all->GetXaxis()->SetRangeUser(0, thisBin);
  hdata->GetXaxis()->SetRangeUser(0, thisBin);
  gdata->GetXaxis()->SetRangeUser(0, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hdata->GetYaxis()->SetRangeUser(yMin, yMax);
  //  hestimate_all->GetXaxis()->LabelsOption("v");
  //  hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->GetXaxis()->SetLabelSize(0.043);
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
  
  hestimate_all->GetYaxis()->SetTitleOffset(0.95);
  hestimate_all->GetYaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");
  
  TH1D* prefit=new TH1D("prefit", "", 1, 0, 1);
  prefit->SetFillColor(0);
  prefit->SetLineColor(0);

  TLegend* legend = new TLegend( 0.8, 0.9-(bgSize+1-1)*0.06-0.06, 0.93, 0.9-0.06 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( hdata, "Data", "PL" );
  //  legend->AddEntry( prefit, "A-priori background", "F" );
  legend->AddEntry( hestimate[0], "Multijet", "F");
  legend->AddEntry( hestimate[1], "Lost lepton", "F");
  legend->AddEntry( hestimate[2], "Z #rightarrow #nu#bar{#nu}", "F");


  //  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi2015,lumi2016);
  labelTop->Draw("same");

  TPaveText* labelCMS = MT2DrawTools::getLabelCMS();
  labelCMS->Draw("same");

//  TLine* lHT[5];
//  for( int iHT=0; iHT < 5; iHT++ ){
//    lHT[iHT-1] = new TLine(12+11*iHT, 0.0, 12+11*iHT, yMax );
//    lHT[iHT-1]->SetLineColor(kBlack);
//    lHT[iHT-1]->SetLineStyle(3);
//    lHT[iHT-1]->SetLineWidth(2);
//
//    lHT[iHT-1]->Draw("same");
//  }



  int nHTRegions = 6;
  std::vector< std::string > htRegions;
  htRegions.push_back("1 Jet");
////  htRegions.push_back("very low H_{T}");
////  htRegions.push_back("low H_{T}");
////  htRegions.push_back("medium H_{T}");
////  htRegions.push_back("high H_{T}");
////  htRegions.push_back("extreme H_{T}");
//  htRegions.push_back("#it{H}_{T} [200,450] GeV");
//  htRegions.push_back("#it{H}_{T} [450,575] GeV");
//  htRegions.push_back("#it{H}_{T} [575,1000] GeV");
//  htRegions.push_back("#it{H}_{T} [1000,1500] GeV");
//  htRegions.push_back("#it{H}_{T} >1500 GeV");
  htRegions.push_back("H_{T} [200,450] GeV");
  htRegions.push_back("H_{T} [450,575] GeV");
  htRegions.push_back("H_{T} [575,1000] GeV");
  htRegions.push_back("H_{T} [1000,1500] GeV");
  htRegions.push_back("H_{T} > 1500 GeV");
  
  TPaveText* htBox[5];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){

    if (iHT==0) htBox[iHT] = new TPaveText(0.12+0.15*iHT, 0.9-0.06+0.02, 0.34+0.15*iHT, 0.85+0.02, "brNDC");
    else htBox[iHT] = new TPaveText(0.13+0.13*iHT, 0.9-0.06+0.02, 0.34+0.13*iHT, 0.85+0.02, "brNDC");
    htBox[iHT]->AddText( htRegions[iHT].c_str() );

    htBox[iHT]->SetBorderSize(0);
    htBox[iHT]->SetFillColor(kWhite);
    htBox[iHT]->SetTextSize(0.035);
    htBox[iHT]->SetTextAlign(21); // align centered
    htBox[iHT]->SetTextFont(62);
    htBox[iHT]->Draw("same");

  }

  TLine* lHT[5];
  for( int iHT=0; iHT < 5; iHT++ ){
    lHT[iHT-1] = new TLine(12+11*iHT, 0.0, 12+11*iHT, yMax );
    lHT[iHT-1]->SetLineColor(kBlack);
    lHT[iHT-1]->SetLineStyle(3);
    lHT[iHT-1]->SetLineWidth(2);

    lHT[iHT-1]->Draw("same");
  }

  legend->Draw("same");


  gPad->RedrawAxis();
  
  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();  

  bool doLogRatio = true;
  
  TH2D* h2_axes_ratio;

  if(doLogRatio){
    
    gPad->SetLogy();
    h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0.1, 10.0 );
    h_Ratio->GetYaxis()->SetRangeUser(0.1, 10.0);
    g_Ratio->GetYaxis()->SetRangeUser(0.1, 10.0);
    
  }
  else
    h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0., 2.0 );
  
  //  TH2D* h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0., 3.0 );
  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.18);
  //  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.4);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitle("Data/Est.");
  
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

    if(hestimate_all->GetBinContent(iBin)>0)
      error = hestimate_all->GetBinError(iBin)/hestimate_all->GetBinContent(iBin);
    else error = hestimate_all->GetBinError(iBin);

    h_band->SetBinError(iBin, error);

  }


  h2_axes_ratio->Draw("");
  h_band->Draw("E2same");
  LineCentral->Draw("same");
  //  h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  

  TLine* lHT_b[6];
  for( int iHT=1; iHT < 6; iHT++ ){
    //    lHT_b[iHT-1] = new TLine(12+11*(iHT-1), 0, 12+11*(iHT-1), 3.0 );
    if(doLogRatio)
      lHT_b[iHT-1] = new TLine(12+11*(iHT-1), 0.1, 12+11*(iHT-1), 10.0 );
    else
      lHT_b[iHT-1] = new TLine(12+11*(iHT-1), 0, 12+11*(iHT-1), 2.0 );

    lHT_b[iHT-1]->SetLineColor(kBlack);
    lHT_b[iHT-1]->SetLineStyle(3);
    lHT_b[iHT-1]->SetLineWidth(2);

    lHT_b[iHT-1]->Draw("same");
  }

  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.pdf", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.C", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.png", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.eps", fullPath.c_str()) );
  

  //  gStyle->SetOptStat(0110);
  gStyle->SetOptFit(0011);
  
  TF1* fgauss= new TF1("fgauss", "gaus", -5, 5);
  fgauss->SetLineColor(2);

  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  
  hPull->SetStats(1100);
  //  hPull->GetYaxis()->SetRangeUser(0, 15);
  hPull->Draw("hist");
  hPull->Fit("fgauss");
  fgauss->Draw("same");
  
  labelTop = MT2DrawTools::getLabelTop(lumi2015,lumi2016);
  labelTop->Draw("same");

  labelCMS = MT2DrawTools::getLabelCMS("CMS Supplementary");
  labelCMS->Draw("same");
  
  TPaveText *arxiv = new TPaveText(0.2, 0.9-0.05, 0.35, 0.9, "brNDC");
  arxiv->AddText( "arXiv:1603.04053" );
  arxiv->SetBorderSize(0);
  arxiv->SetFillColor(kWhite);
  //  arxiv->SetTextSize(0.035);
  arxiv->SetTextAlign(11); // align centered                                                                                                                                  
  arxiv->SetTextFont(42);
  arxiv->Draw("same");
  
  c3->SaveAs( Form("%s/PullDistribution.pdf", fullPath.c_str()) );
  c3->SaveAs( Form("%s/PullDistribution.png", fullPath.c_str()) );
  c3->SaveAs( Form("%s/PullDistribution.root", fullPath.c_str()) );

  TCanvas* c4 = new TCanvas("c4", "", 600, 600);
  c4->cd();
  hPvalue->SetStats(1110);
  hPvalue->Draw("hist");
  c4->SaveAs( Form("%s/PvalueDistribution.pdf", fullPath.c_str()) );
  c4->SaveAs( Form("%s/PvalueDistribution.png", fullPath.c_str()) );
  
}



