#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2DrawTools.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>


#include "TMath.h"
#include "TTreeFormula.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"

#include "RooHistError.h"

float lumi =3.0; //fb-1 

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::vector<MT2Analysis<MT2Estimate>* > bgYields );

int main() {

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = "./YieldComparison_QCDv5_allBins/";

//  std::string firstInputFile  = "./EventYields_mc_Spring15_v3_dummy/zinvFromGamma.root";
//  std::string secondInputFile = "./EventYields_mc_Spring15_v3_dummy/zinvFromGamma.root";
//  std::string firstInputFile  = "./EventYields_mc_Spring15_QCDSR_onlyHT_dummy/analyses.root";
//  std::string secondInputFile = "/shome/casal/mt2ana/analysis/qcdINCpred.root";
//  std::string firstInputFile  = "./EventYields_mc_Spring15_QCDSR_dummy/analyses.root";
//  std::string secondInputFile = "/shome/casal/mt2ana/analysis/qcdSRpred.root";
  std::string firstInputFile  = "./EventYields_mc_Spring15_QCDSR_allBins_dummy/analyses.root";
  std::string secondInputFile = "/shome/casal/mt2ana/analysis/qcdSRpred.root";
//  std::string firstInputFile  = "/shome/pandolf/qcdEstimate.root";
//  std::string secondInputFile = "/shome/pandolf/qcdEstimate.root";
  
  MT2Analysis<MT2Estimate>* analysisFirst = MT2Analysis<MT2Estimate>::readFromFile( secondInputFile.c_str(), "estimateSR" );
  MT2Analysis<MT2Estimate>* analysisSecond = MT2Analysis<MT2Estimate>::readFromFile( firstInputFile.c_str(), "QCD");
//  MT2Analysis<MT2Estimate>* analysisFirst = MT2Analysis<MT2Estimate>::readFromFile( secondInputFile.c_str(), "qcdEstimate" );
//  MT2Analysis<MT2Estimate>* analysisSecond = MT2Analysis<MT2Estimate>::readFromFile( firstInputFile.c_str(), "mcTruth");

  std::vector< MT2Analysis<MT2Estimate> *> bgEstimate;
  bgEstimate.push_back( analysisSecond );

  drawYields(outputdir.c_str(), analysisFirst, bgEstimate );

  return 0;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>*  analysisFirst, std::vector< MT2Analysis<MT2Estimate> *> bgYields ) {

  MT2DrawTools::setStyle();

  system(Form("mkdir -p %s", outputdir.c_str()));

  std::vector<int> colors;
  for(unsigned int b=0; b<bgYields.size(); ++b)
    //    colors.push_back( 417 );
    colors.push_back( 401 );
  
  std::set<MT2Region> MT2Regions = analysisFirst->getRegions();
  
  TH1D* hdata = new TH1D("hdata", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  hdata->Sumw2();
  hdata->GetYaxis()->SetTitle("Entries");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( kBlack );
  hdata->SetMarkerColor( kBlack );
  
  TH1D* hestimate = new TH1D("hestimate", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  hestimate->Sumw2();
  hestimate->GetYaxis()->SetTitle("Entries");
  //hestimate->SetFillColor(colors[0]);
  hestimate->SetFillColor(0);
  hestimate->SetLineColor(colors[0]);
  hestimate->SetMarkerColor(colors[0]);
  hestimate->SetMarkerStyle(20);
  hestimate->SetMarkerSize(1.6);
  
  TH1D* hPull = new TH1D("hPull", "", 20, -5, 5);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Data Driven - MC)/#sigma");
  hPull->GetYaxis()->SetTitle("Entries");
  
  std::string thisName = Form("%s_ratio", hdata->GetName());
  TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
  h_Ratio->Sumw2();
  h_Ratio->GetYaxis()->SetTitle("(Data Driven - MC)/#sigma");
  
  
  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::string fullPath = outputdir;

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      TH1D* h_first = analysisFirst->get(*iMT2)->yield;
      TGraphAsymmErrors* g_first = new TGraphAsymmErrors(0); // = new TGraphAsymmErrors(h_first);
      //TGraphAsymmErrors* g_first = MT2DrawTools::getPoissonGraph(h_first);      
      
      int nBins = h_first->GetXaxis()->GetNbins()+1;

      for( int iBin=1; iBin<(h_first->GetXaxis()->GetNbins()+1); ++iBin ) {

	//if(h_first->GetBinCenter(iBin)-h_first->GetBinWidth(iBin)/2.<200.) continue;

	double y; // these are data histograms, so y has to be integer                                                                                                                                                                            
	double x, xerr, yerrplus, yerrminus;
	
	
	x = h_first->GetBinCenter(iBin);
	xerr = h_first->GetBinWidth(iBin)/2.;
	
	y = h_first->GetBinContent(iBin);

	double yerr;
	double y_ = h_first->IntegralAndError(iBin, iBin, yerr);
	
	if(y>0)
	  yerr=yerr/y;
	else
	  yerr=0;
	double yerr2 = yerr*yerr;
	
	double eB;
	if(iMT2->nBJetsMin()==0)
	  eB=0.25;
	else if(iMT2->nBJetsMin()==1)
	  eB=0.25;
	else if(iMT2->nBJetsMin()==2)
	  eB=0.50;
	else
	  eB=0.60;

	double eJ;
	if(iMT2->htMin()<575){

	  if(iMT2->nJetsMin()==2)
	    eJ=6.582806e-01;
	  
	  else if(iMT2->nJetsMin()==4)
	    eJ=4.928055e-01;
	  
	  else
	    eJ=1.0;

	}

	else if(iMT2->htMin()<1000){

	  if(iMT2->nJetsMin()==2)
	    eJ=5.974327e-01;
	  
	  else if(iMT2->nJetsMin()==4)
	    eJ=3.365912e-01;
	  
	  else
	    eJ=3.856405e-01;

	}

	else if(iMT2->htMin()<1500){

	  if(iMT2->nJetsMin()==2)
	    eJ=4.174673e-01;
	  
	  else if(iMT2->nJetsMin()==4)
	    eJ=2.556875e-01;
	  
	  else
	    eJ=3.234015e-01;

	}
	
	else{

	  if(iMT2->nJetsMin()==2)
            eJ=4.620245e-01;

          else if(iMT2->nJetsMin()==4)
            eJ=1.607043e-01;

          else
            eJ=1.778792e-01;

	}
	
	
	yerr2 += eB*eB + eJ*eJ;
	yerr = y*TMath::Sqrt(yerr2);

//	double ym, yp;
//	RooHistError::instance().getPoissonInterval(y,ym,yp,1);
//	yerrplus = yp - y;
//	yerrminus = y - ym;

	h_first->SetBinError(iBin, yerr);

	int thisPoint = g_first->GetN();
	g_first->SetPoint( thisPoint, x, y );
	//	g_first->SetPointError( thisPoint, xerr, xerr, yerrminus, yerrplus );
	g_first->SetPointError( thisPoint, xerr, xerr, yerr, yerr );

      }
      
      TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), iMT2->getName().c_str()), "recreate" );
      histoFile->cd();
      h_first->Write();
      
      g_first->SetMarkerStyle(20);
      g_first->SetMarkerSize(1.6);
      g_first->SetLineColor( kBlack );
      g_first->SetMarkerColor( kBlack );

      //      int firstBin = h_first->FindBin(201);
      
      int firstBin=1;
      double err_data;
      double int_data = h_first->IntegralAndError(firstBin, nBins, err_data);
 
      hdata->SetBinContent(iRegion, int_data);
      hdata->SetBinError(iRegion, err_data);
      
      std::cout << niceNames[1].c_str() << "\t" << int_data << "\t" << err_data << std::endl;

      hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
      
      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.15);
      pad1->Draw();
      pad1->cd();

      //THStack bgStack("bgStack", "");
      TH1D* h_second = bgYields[0]->get(*iMT2)->yield;

      //      TGraphAsymmErrors* g_second = bgYields[0]->get(*iMT2)->getGraph();
      //TGraphAsymmErrors* g_second = MT2DrawTools::getPoissonGraph(h_second);
      TGraphAsymmErrors* g_second = new TGraphAsymmErrors(0);
      
      for( int iBin=1; iBin<(h_second->GetXaxis()->GetNbins()+1); ++iBin ) {

	//if(h_second->GetBinCenter(iBin)-h_second->GetBinWidth(iBin)/2.<200.) continue;

        double y;                                                                                                                                                                                                                          
        double x, xerr, yerrplus, yerrminus;

        x = h_second->GetBinCenter(iBin);
        xerr = h_second->GetBinWidth(iBin)/2.;

        y = h_second->GetBinContent(iBin);
	float yerr = h_second->GetBinError(iBin);
	
//        double ym, yp;
//	RooHistError::instance().getPoissonInterval(y,ym,yp,1);
//        yerrplus = yp - y;
//        yerrminus = y - ym;

        int thisPoint = g_second->GetN();
        g_second->SetPoint( thisPoint, x, y );
	//        g_second->SetPointError( thisPoint, xerr, xerr, yerrminus, yerrplus );
	//        g_second->SetPointError( thisPoint, xerr, xerr, 0., 0. );
        g_second->SetPointError( thisPoint, xerr, xerr, yerr, yerr );

      }

      g_second->SetLineColor( colors[0] );
      g_second->SetMarkerColor( colors[0] );
      g_second->SetMarkerStyle(20);
      g_second->SetMarkerSize(1.6);
      g_second->Write();
      
      firstBin=1;
      double err_int;
      double integral = h_second->IntegralAndError(firstBin, nBins, err_int);
      hestimate->SetBinContent(iRegion, integral);
      hestimate->SetBinError(iRegion, err_int);

      std::cout << niceNames[1].c_str() << "\t" << integral << "\t" << err_int << std::endl;
      
      hestimate->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );

      if(int_data>0 && integral>0){
	hPull->Fill((int_data-integral)/TMath::Sqrt(err_data*err_data+err_int*err_int));
	h_Ratio->SetBinContent(iRegion, (int_data-integral)/TMath::Sqrt(err_data*err_data+err_int*err_int));
	h_Ratio->SetBinError(iRegion, TMath::Sqrt(err_data*err_data+err_int*err_int)/TMath::Sqrt(err_data*err_data+err_int*err_int));
      }

//      for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
//        int index = bgYields.size() - i - 1;
//        TH1D* h_second_ = bgYields[index]->get(*iMT2)->yield;
//        h_second_->SetFillColor( colors[index] );
//        h_second_->SetLineColor( kBlack );
//        bgStack.Add(h_second_);
//	if( i>0 )
//	  h_second->Add(h_second_);
//      }

      float xMin = h_first->GetXaxis()->GetXmin();
      float xMax = h_first->GetXaxis()->GetXmax();
      float yMax_1 = h_first->GetMaximum()*1.5;
      float yMax_2 = 1.2*(h_first->GetMaximum() + h_first->GetBinError(h_first->GetMaximumBin()));
      float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
      float yMax_3 = h_second->GetMaximum()*1.5;
      float yMax_4 = 1.2*(h_second->GetMaximum() + h_second->GetBinError(h_second->GetMaximumBin()));
      float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
      float yMax = (yMax1>yMax2) ? yMax1 : yMax2;

      TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
      h2_axes->SetXTitle("M_{T2} [GeV]");
      h2_axes->SetYTitle("Entries");

      h2_axes->Draw();
 
      //std::vector<std::string> niceNames = iMT2->getNiceNames();
      for( unsigned i=0; i<niceNames.size(); ++i ) {

        float yMax = 0.9-(float)i*0.05;
        float yMin = yMax - 0.05;
        TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
        regionText->SetTextSize(0.035);
        regionText->SetTextFont(42);
        regionText->SetFillColor(0);
        regionText->SetTextAlign(11);
        regionText->AddText( niceNames[i].c_str() );
        regionText->Draw("same");
	
      }


      TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
//      legend->AddEntry( g_first, "data-driven", "P" );
//      legend->AddEntry( g_second, "MC Z(#nu#nu)+jets", "P" );
      legend->AddEntry( g_first, "data-driven", "P" );
      legend->AddEntry( g_second, "MC QCD", "P" );

      legend->Draw("same");

      //      bgStack.Draw("histoE, same");
      g_second->Draw("pe, same");
      g_first->Draw("pe, same");

      TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
      labelTop->Draw("same");

      gPad->RedrawAxis();

      c1->cd();
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.1);
      pad2->Draw();
      pad2->cd();

      std::string thisName = Form("%s_ratio", h_first->GetName());
      TH1D* h_ratio = (TH1D*) h_first->Clone(thisName.c_str());
      h_ratio->Divide(h_second);
      h_ratio->Write();
      h_ratio->SetStats(0);	    
      h_ratio->SetMarkerStyle(20);
      h_ratio->SetLineColor(1);
      //      h_ratio->SetMarkerSize(0.02);
      h_ratio->GetXaxis()->SetLabelSize(0.00);
      h_ratio->GetXaxis()->SetTickLength(0.09);
      h_ratio->GetYaxis()->SetNdivisions(5,5,0);
      h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
      h_ratio->GetYaxis()->SetTitleSize(0.17);
      h_ratio->GetYaxis()->SetTitleOffset(0.4);
      h_ratio->GetYaxis()->SetLabelSize(0.17);
      h_ratio->GetYaxis()->SetTitle("Ratio");
      
      

      h_ratio->SetLineWidth(2);
      //h_ratio->Draw("PE");
      
      TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, 0.0, 2.0 );
      
      TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
      lineCentral->SetLineColor(1);
      
//      TGraphAsymmErrors* graph_ = g_first;
//      TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(0); // MT2DrawTools::getRatioGraph( h_first, h_second );
//
//      for( int i=0; i < graph_->GetN(); ++i){
//
//	Double_t x_tmp, y_tmp, errUp, errDown;
//	graph_->GetPoint( i, x_tmp, y_tmp );
//
//	errUp   = graph_->GetErrorYhigh(i);
//	errDown = graph_->GetErrorYlow(i);
//
//	int iBin = h_second->FindBin(x_tmp);
//	float mc = h_second->GetBinContent(iBin);
//	g_ratio->SetPoint(i, x_tmp, y_tmp/mc);
//	g_ratio->SetPointEYhigh(i, errUp/mc);
//	g_ratio->SetPointEYlow(i, errDown/mc);
//
//      }
//
//      g_ratio->SetLineColor(1);
//      g_ratio->SetMarkerColor(1);
//      g_ratio->SetMarkerStyle(20);

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      //      g_ratio->Draw("pe,same");
      h_ratio->Draw("pe,same");

      gPad->RedrawAxis();
      
      c1->cd();

      c1->SaveAs( Form("%s/mt2_%s.pdf", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.png", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;
      delete h2_axes;
      //delete g_ratio;
      delete h_ratio;
      delete h_first;
      delete h_second;
      
      ++iRegion;

  } // for MT2 regions

  
  TCanvas* c2 = new TCanvas("c2", "", 1200, 600);
  c2->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  pad1->SetLogy();
  
  float yMax_1 = hdata->GetMaximum()*1.5;
  float yMax_2 = 1.2*(hdata->GetMaximum() + hdata->GetBinError(hestimate->GetMaximumBin()));
  float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
  float yMax_3 = hestimate->GetMaximum()*1.5;
  float yMax_4 = 1.2*(hestimate->GetMaximum() + hestimate->GetBinError(hestimate->GetMaximumBin()));
  float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
  float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
  
  float yMin = 1e-1;
  yMax*=20.;
  
  hestimate->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  hestimate->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate->GetXaxis()->LabelsOption("v");
  hestimate->Draw("PE");

  TGraphAsymmErrors* g_data = new TGraphAsymmErrors(0);
  for( int iBin=1; iBin<(hdata->GetXaxis()->GetNbins()+1); ++iBin ) {

    double y;
    double x, xerr, yerrplus, yerrminus;

    x = hdata->GetBinCenter(iBin);
    xerr = hdata->GetBinWidth(iBin)/2.;

    y = hdata->GetBinContent(iBin);
    double yerr = hdata->GetBinError(iBin);

//    double ym, yp;
//    RooHistError::instance().getPoissonInterval(y,ym,yp,1);
//    yerrplus = yp - y;
//    yerrminus = y - ym;

    int thisPoint = g_data->GetN();
    g_data->SetPoint( thisPoint, x, y );
    //    g_data->SetPointError( thisPoint, xerr, xerr, yerrminus, yerrplus );
    g_data->SetPointError( thisPoint, xerr, xerr, yerr, yerr );

  }
  
  //g_data->SetMarkerStyle(20);
  //g_data->SetMarkerSize(1.6);
  //g_data->SetLineColor( kBlack );
  //g_data->SetMarkerColor( kBlack );

  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( kBlack );
  hdata->SetMarkerColor( kBlack );

  hdata->Draw("pe,same");
  //g_data->Draw("pe,same");


//  double epart;
//  std::cout<< "estimate: " << hdata->IntegralAndError(1,11, epart) << "\t" << epart << std::endl;
//  std::cout<< "MC: " << hestimate->IntegralAndError(1,11, epart) << "\t" << epart << std::endl;
//
//  std::cout<< "estimate: " << hdata->IntegralAndError(12,22, epart) << "\t" << epart << std::endl;
//  std::cout<< "MC: " << hestimate->IntegralAndError(12,22, epart) << "\t" << epart << std::endl;
//
//  std::cout<< "estimate: " << hdata->IntegralAndError(23,33, epart) << "\t" << epart << std::endl;
//  std::cout<< "MC: " << hestimate->IntegralAndError(23,33, epart) << "\t" << epart << std::endl;
//
//  std::cout<< "estimate: " << hdata->IntegralAndError(34,44, epart) << "\t" << epart << std::endl;
//  std::cout<< "MC: " << hestimate->IntegralAndError(34,44, epart) << "\t" << epart << std::endl;

  TLegend* legend = new TLegend( 0.8, 0.9-(bgYields.size()+1)*0.06-0.06, 0.93, 0.9-0.06 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
//  legend->AddEntry( hdata, "data-driven", "P" );
//  legend->AddEntry( hestimate, "MC Z(#nu#nu)+jets", "F" );
  legend->AddEntry( hdata, "data-driven", "PL" );
  legend->AddEntry( hestimate, "QCD MC", "PL" );

  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  labelTop->Draw("same");
  
  TLine* lHT[3];
  for( int iHT=1; iHT < 4; iHT++ ){
    lHT[iHT-1] = new TLine(11*iHT, 0.0, 11*iHT, yMax );
    lHT[iHT-1]->SetLineColor(kBlack);
    lHT[iHT-1]->SetLineStyle(3);
    lHT[iHT-1]->SetLineWidth(2);

    lHT[iHT-1]->Draw("same");
  }

  int nHTRegions = 4;
  std::vector< std::string > htRegions;
  htRegions.push_back("low H_{T}");
  htRegions.push_back("medium H_{T}");
  htRegions.push_back("high H_{T}");
  htRegions.push_back("extreme H_{T}");
  
  TPaveText* htBox[nHTRegions];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){
    
    htBox[iHT] = new TPaveText(0.16+0.2*iHT, 0.9-0.06, 0.34+0.2*iHT, 0.9, "brNDC");
    htBox[iHT]->AddText( htRegions[iHT].c_str() );
    
    htBox[iHT]->SetBorderSize(0);
    htBox[iHT]->SetFillColor(kWhite);
    htBox[iHT]->SetTextSize(0.038);
    htBox[iHT]->SetTextAlign(21); // align centered
    htBox[iHT]->SetTextFont(62);
    htBox[iHT]->Draw("same");

  }

  gPad->RedrawAxis();
  
  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();

//  std::string thisName = Form("%s_ratio", hdata->GetName());
//  TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
//  h_Ratio->Divide(hestimate);
//  h_Ratio->Write();
//  h_Ratio->SetStats(0);
//  h_Ratio->SetMarkerStyle(20);
//  h_Ratio->SetLineColor(1);
//  h_Ratio->GetXaxis()->SetLabelSize(0.00);
//  h_Ratio->GetXaxis()->SetTickLength(0.09);
//  h_Ratio->GetYaxis()->SetNdivisions(5,5,0);
//  h_Ratio->GetYaxis()->SetRangeUser(0.0,2.0);
//  h_Ratio->GetYaxis()->SetTitleSize(0.17);
//  h_Ratio->GetYaxis()->SetTitleOffset(0.4);
//  h_Ratio->GetYaxis()->SetLabelSize(0.17);
//  h_Ratio->GetYaxis()->SetTitle("Ratio");
  
  h_Ratio->Write();
  h_Ratio->SetStats(0);
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetMarkerSize(1.6);
  h_Ratio->SetLineColor(1);
  h_Ratio->GetXaxis()->SetLabelSize(0.00);
  h_Ratio->GetXaxis()->SetTickLength(0.09);
  h_Ratio->GetYaxis()->SetNdivisions(5,5,0);
  h_Ratio->GetYaxis()->SetRangeUser(-2,2);
  h_Ratio->GetYaxis()->SetTitleSize(0.17);
  h_Ratio->GetYaxis()->SetTitleOffset(0.4);
  h_Ratio->GetYaxis()->SetLabelSize(0.17);
  //  h_Ratio->GetYaxis()->SetTitle("Ratio");

  //  h_Ratio->SetLineWidth(2);
  //h_Ratio->Draw("P,E");
  

//  TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( 0, 44, 0.0, 2.0 );
//  //TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( 0, 4, 0.0, 2.0 );

  TH2D* h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, 44, 10, -2.0, 2.0 );
  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.3);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitle("Pull");
  
//  TLine* LineCentral = new TLine(0, 1.0, 44, 1.0);
//  //TLine* LineCentral = new TLine(0, 1.0, 4, 1.0);
  TLine* LineCentral = new TLine(0, 0., 44, 0.);
  //TLine* LineCentral = new TLine(0, 1.0, 4, 1.0);
  LineCentral->SetLineColor(1);

//  TGraphAsymmErrors* gData_ = g_data;
//  TGraphAsymmErrors* gRatio = new TGraphAsymmErrors(0); // MT2DrawTools::getRatioGraph( h_first, h_second );                                                                                                                          
//
//  for( int i=0; i < gData_->GetN(); ++i){
//
//    Double_t x_tmp, y_tmp, errUp, errDown;
//    gData_->GetPoint( i, x_tmp, y_tmp );
//
//    errUp   = gData_->GetErrorYhigh(i);
//    errDown = gData_->GetErrorYlow(i);
//
//    int iBin = hestimate->FindBin(x_tmp);
//    float mc = hestimate->GetBinContent(iBin);
//    gRatio->SetPoint(i, x_tmp, y_tmp/mc);
//    gRatio->SetPointEYhigh(i, errUp/mc);
//    gRatio->SetPointEYlow(i, errDown/mc);
//
//  }
//
//  gRatio->SetLineColor(1);
//  gRatio->SetMarkerColor(1);
//  gRatio->SetMarkerStyle(20);

  h2_axes_ratio->Draw("");
  LineCentral->Draw("same");
  h_Ratio->Draw("pe,same");
  //  gRatio->Draw("pe,same");

  TLine* lHT_b[3];
  for( int iHT=1; iHT < 4; iHT++ ){
    lHT_b[iHT-1] = new TLine(11*iHT,-2.0, 11*iHT, 2.0 );
    lHT_b[iHT-1]->SetLineColor(kBlack);
    lHT_b[iHT-1]->SetLineStyle(3);
    lHT_b[iHT-1]->SetLineWidth(2);

    lHT_b[iHT-1]->Draw("same");
  }


  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/mt2_fullEstimate.pdf", outputdir.c_str()) );
  c2->SaveAs( Form("%s/mt2_fullEstimate.png", outputdir.c_str()) );


  gStyle->SetOptStat(1110);
  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  hPull->SetStats(1110);
  hPull->Draw("hist");
  std::cout << hPull->GetEntries() << "\t" << hPull->GetMean() << "\t" << hPull->GetRMS() << std::endl; 
  c3->SaveAs( Form("%s/PullDistribution.pdf", outputdir.c_str()) );
  c3->SaveAs( Form("%s/PullDistribution.png", outputdir.c_str()) );
  
}
