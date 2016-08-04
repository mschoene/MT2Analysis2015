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


float lumi =12.9; //fb-1 


void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data_ext, MT2Analysis<MT2Estimate>* data_bin, std::vector<MT2Analysis<MT2Estimate>* > bgYields );

int main() {

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = "./YieldComparison_lostLepton_binned_12p9ifb_ICHEP/";

//  std::string firstInputFile  = "./EventYields_data_Run2015_25nsGolden_2p3ifb/llepEstimate_binbybin.root";
//  std::string secondInputFile = "./EventYields_data_Run2015_25nsGolden_2p3ifb/llepEstimate_extrapolation.root";
  std::string firstInputFile  = "./EventYields_data_Run2016_12p9ifb_ICHEP/llepEstimate_binbybin.root";
  std::string secondInputFile = "./EventYields_data_Run2016_12p9ifb_ICHEP/llepEstimate_extrapolation.root";
  
>>>>>>> 1e26783bc38f51eef1352958015fd7c820c10d66
  MT2Analysis<MT2Estimate>* analysisFirst_ext = MT2Analysis<MT2Estimate>::readFromFile( secondInputFile.c_str(), "llepEstimate" ); 
  analysisFirst_ext->setName("Data-drien (Std.)");
  MT2Analysis<MT2Estimate>* analysisFirst_bin = MT2Analysis<MT2Estimate>::readFromFile( firstInputFile.c_str(), "llepEstimate" );
  analysisFirst_ext->setName("Data-driven (bin)");

  MT2Analysis<MT2Estimate>* analysisSecond = MT2Analysis<MT2Estimate>::readFromFile( firstInputFile.c_str(), "Top + W+jets");
  analysisSecond->setName("MC");
  
  std::vector< MT2Analysis<MT2Estimate> *> bgEstimate;
  bgEstimate.push_back( analysisSecond );

  drawYields(outputdir.c_str(), analysisFirst_ext, analysisFirst_bin, bgEstimate );

  return 0;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>*  analysisFirst_ext, MT2Analysis<MT2Estimate>*  analysisFirst_bin, std::vector< MT2Analysis<MT2Estimate> *> bgYields ) {

  MT2DrawTools::setStyle();

  system(Form("mkdir -p %s", outputdir.c_str()));

  std::vector<int> colors;
  for(unsigned int b=0; b<bgYields.size(); ++b)
    colors.push_back( 418 );
    //colors.push_back( 401 );
  
  std::set<MT2Region> MT2Regions = analysisFirst_ext->getRegions();
  
  TH1D* hdata = new TH1D("hdata", "", 174, 0, 174);
  //  TH1D* hdata = new TH1D("hdata", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  hdata->Sumw2();
  hdata->GetYaxis()->SetTitle("Entries");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.0);
  hdata->SetLineColor( 4 );
  hdata->SetMarkerColor( 4 );

  TH1D* hdata_bin = new TH1D("hdata_bin", "", 174, 0, 174);
  //  TH1D* hdata_bin = new TH1D("hdata_bin", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  hdata_bin->Sumw2();
  hdata_bin->GetYaxis()->SetTitle("Entries");
  hdata_bin->SetMarkerStyle(24);
  hdata_bin->SetMarkerSize(1.0);
  hdata_bin->SetLineColor( 2 );
  hdata_bin->SetMarkerColor( 2 );
  
  TH1D* hestimate = new TH1D("hestimate", "", 174, 0, 174);
  //  TH1D* hestimate = new TH1D("hestimate", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  hestimate->Sumw2();
  hestimate->GetYaxis()->SetTitle("Entries");
  //hestimate->SetFillColor(colors[0]);
  hestimate->SetFillColor(0);
  hestimate->SetLineColor(colors[0]);
  hestimate->SetMarkerColor(colors[0]);
  hestimate->SetMarkerStyle(21);
  hestimate->SetMarkerSize(1.6);
  
  TH1D* hPull = new TH1D("hPull", "", 20, -5, 5);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Data Driven - MC)/#sigma");
  hPull->GetYaxis()->SetTitle("Entries");
  
  //std::string thisName = Form("%s_ratio", hdata->GetName());
  //TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
  //h_Ratio->Sumw2();
  //h_Ratio->GetYaxis()->SetTitle("(Data Driven - MC)/#sigma");
  
  std::string fullPath = outputdir;
  
  std::string labelsMono[12]={"[200,250]","[250,350]","[350,450]","[450,575]","[575,700]","[700,1000]",">1000", "[200,250]","[250,350]","[350,450]","[450,575]",">575"};

  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      TH1D* h_first = analysisFirst_ext->get(*iMT2)->yield;
      TGraphAsymmErrors* g_first = new TGraphAsymmErrors(0); // = new TGraphAsymmErrors(h_first);
      //TGraphAsymmErrors* g_first = MT2DrawTools::getPoissonGraph(h_first);      

      TH1D* h_first_bin = analysisFirst_bin->get(*iMT2)->yield;
      TGraphAsymmErrors* g_first_bin = new TGraphAsymmErrors(0); // = new TGraphAsymmErrors(h_first);
      
      int nBins;
      double *bins;
      iMT2->getBins(nBins, bins);
      //std::cout << niceNames[0] << ", " << niceNames[1] <<": " << nBins << "\t" << bins[0] << "\t" << bins[nBins-1] << std::endl;
      //      int nBins = h_first->GetXaxis()->GetNbins()+1;

      float absoluteErr=0.;
      for( int iBin=1; iBin<(h_first->GetXaxis()->GetNbins()+1); ++iBin ) {

	//if(h_first->GetBinCenter(iBin)-h_first->GetBinWidth(iBin)/2.<200.) continue;

	double y; // these are data histograms, so y has to be integer
	double x, xerr, yerrplus, yerrminus;
	
	x = h_first->GetBinCenter(iBin);
	xerr = h_first->GetBinWidth(iBin)/2.;
	
	y = h_first->GetBinContent(iBin);

	double yerr;
	double y_ = h_first->IntegralAndError(iBin, iBin, yerr);

	double yerr2;
	if(y>0)
	  yerr2=(yerr/y)*(yerr/y);
	else
	  yerr2=1.0;

	//ADDITION
	yerr2=0.;

	int nbins = h_first->GetNbinsX();
	float relativeErr;
	if(y>0)
	  relativeErr = 0.4 / (nbins-1) * (iBin-1);
	else
	  relativeErr = 0.0;
	
	absoluteErr+=relativeErr*y;
	
	if( iBin==1 ) { 
	  h_first->SetBinError(iBin, absoluteErr);
	  yerr=absoluteErr;
	}
	else {
	  h_first->SetBinError(iBin, relativeErr*y);
	  yerr=relativeErr*y;
	}
	
	
//	yerr2+=0.15*0.15;
//	yerr=y*sqrt(yerr2);
//	
	std::cout << iBin << "\t" << yerr/y << std::endl;
//	h_first->SetBinError(iBin, yerr);

	int thisPoint = g_first->GetN();
	g_first->SetPoint( thisPoint, x, y );
	//	g_first->SetPointError( thisPoint, xerr, xerr, yerrminus, yerrplus );
	g_first->SetPointError( thisPoint, xerr, xerr, yerr, yerr );



	double y_bin; // these are data histograms, so y has to be integer
	double x_bin, xerr_bin;

	x_bin = h_first_bin->GetBinCenter(iBin);
	xerr_bin = h_first_bin->GetBinWidth(iBin)/2.;
	
	y_bin = h_first_bin->GetBinContent(iBin);

	double yerr_bin;
	double y_bin_ = h_first_bin->IntegralAndError(iBin, iBin, yerr_bin);

	double yerr2_bin;
	if(y_bin>0)
	  yerr2_bin=(yerr_bin/y_bin)*(yerr_bin/y_bin);
	else
	  yerr2_bin=1.0;

	yerr_bin=y_bin*sqrt(yerr2_bin);
     

	//ADDITION
	yerr2_bin=0.;
	float tot = h_first_bin->Integral();
	float p;
	if(tot>0) p = y_bin/tot;
	else p=0;
	
	yerr2_bin=tot*p*(1-p);
	yerr_bin=sqrt(yerr2_bin);

	
	h_first_bin->SetBinError(iBin, yerr_bin);
	
	
	int thisPoint_bin = g_first_bin->GetN();
	g_first_bin->SetPoint( thisPoint_bin, x_bin, y_bin );
	g_first_bin->SetPointError( thisPoint_bin, xerr_bin, xerr_bin, yerr_bin, yerr_bin );

      }

      if( h_first->GetNbinsX()>1 && absoluteErr>0 )
	h_first->SetBinError(1, absoluteErr);
      else h_first->SetBinError(1, 0.);
      std::cout << "Error on first bin: " <<  absoluteErr/h_first->GetBinContent(1);

      TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), iMT2->getName().c_str()), "recreate" );
      histoFile->cd();
      h_first->Write();
      h_first_bin->Write();
      
      g_first->SetMarkerStyle(20);
      g_first->SetMarkerSize(1.);
      g_first->SetLineColor( 4 );
      g_first->SetMarkerColor( 4 );

      g_first_bin->SetMarkerStyle(24);
      g_first_bin->SetMarkerSize(1.);
      g_first_bin->SetLineColor( 2 );
      g_first_bin->SetMarkerColor( 2 );


      //THStack bgStack("bgStack", "");
      TH1D* h_second = bgYields[0]->get(*iMT2)->yield;

      TGraphAsymmErrors* g_second = new TGraphAsymmErrors(0);
      
      for( int iBin=1; iBin<(h_second->GetXaxis()->GetNbins()+1); ++iBin ) {

	//if(h_second->GetBinCenter(iBin)-h_second->GetBinWidth(iBin)/2.<200.) continue;

        double y;                                                                                                                                                                                                                          
        double x, xerr, yerrplus, yerrminus;

        x = h_second->GetBinCenter(iBin);
        xerr = h_second->GetBinWidth(iBin)/2.;

        y = h_second->GetBinContent(iBin);
	float yerr = h_second->GetBinError(iBin);
	

        int thisPoint = g_second->GetN();
        g_second->SetPoint( thisPoint, x, y );
        g_second->SetPointError( thisPoint, xerr, xerr, yerr, yerr );

      }

      g_second->SetLineColor( colors[0] );
      g_second->SetMarkerColor( colors[0] );
      g_second->SetMarkerStyle(20);
      g_second->SetMarkerSize(1.6);
      g_second->Write();
      

      int firstBin=1;
      double err_data, err_data_bin;
      double int_data, int_data_bin;
      for (int b=1; b<=nBins; ++b){
	
	int_data = h_first->GetBinContent(b);
	err_data = h_first->GetBinError(b);
	
	int_data_bin = h_first_bin->GetBinContent(b);
	err_data_bin = h_first_bin->GetBinError(b);
	
	hdata->SetBinContent(iRegion, int_data);
	hdata->SetBinError(iRegion, err_data);
	
	hdata_bin->SetBinContent(iRegion, int_data_bin);
	hdata_bin->SetBinError(iRegion, err_data_bin);
	
	
	double err_int = h_second->GetBinError(b);
	double integral = h_second->GetBinContent(b);
	hestimate->SetBinContent(iRegion, integral);
	hestimate->SetBinError(iRegion, err_int);
	
//	if( iMT2->nJetsMax()==1 ){
//	  hestimate->GetXaxis()->SetBinLabel( iRegion, labelsMono[iRegion-1].c_str() );
//	  hdata->GetXaxis()->SetBinLabel( iRegion, labelsMono[iRegion-1].c_str() );
//	}
//	else{
//	  std::string thisLabel;
//	  if( b < nBins )
//	    thisLabel=Form("[%.0lf,%.0lf]", bins[b-1], bins[b]);
//	  else
//	    thisLabel=Form(">%.0lf", bins[b-1]);
//	  hestimate->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
//	  hdata->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
//	}

	std::string thisLabel=Form("%s,%d", niceNames[1].c_str(), b); 
	hestimate->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
	hdata->GetXaxis()->SetBinLabel( iRegion,thisLabel.c_str() );
	
	if(int_data>0 && integral>0){
	  hPull->Fill((int_data-integral)/TMath::Sqrt(err_data*err_data+err_int*err_int));
	}
	
	std::cout <<"iRegion " << iRegion << std::endl;
	
	++iRegion;
	
      }
      
      
      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.15);
      pad1->Draw();
      pad1->cd();

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

      if(iMT2->nJetsMax()==1)
	h2_axes->SetXTitle("H_{T} [GeV]");

      h2_axes->Draw();
 
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

      
      TPaveText* bkg = new TPaveText( 0.6, 0.9-0.06, 0.93, 0.9, "brNDC" );
      bkg->SetTextSize(0.038);
      //      bkg->SetTextFont(42);
      bkg->SetFillColor(0);
      bkg->SetTextAlign(11);
      bkg->AddText("Lost lepton background");
      bkg->Draw("same");
      
      TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+1)*0.06-0.06, 0.93, 0.9-0.06 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
//      legend->AddEntry( g_first, "data-driven (standard)", "P" );
//      legend->AddEntry( g_first_bin, "data-driven (bin by bin)", "P" );
      legend->AddEntry( g_first, "Standard", "P" );
      legend->AddEntry( g_first_bin, "Bin by bin", "P" );
      //      legend->AddEntry( g_second, "MC", "P" );

      legend->Draw("same");

      //      g_second->Draw("pe, same");
      g_first->Draw("pe, same");
      g_first_bin->Draw("pe, same");

      //      TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
      TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
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
      h_ratio->SetMarkerColor(4);
      h_ratio->SetLineColor(4);
      h_ratio->GetXaxis()->SetLabelSize(0.00);
      h_ratio->GetXaxis()->SetTickLength(0.09);
      h_ratio->GetYaxis()->SetNdivisions(5,5,0);
      h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
      h_ratio->GetYaxis()->SetTitleSize(0.17);
      h_ratio->GetYaxis()->SetTitleOffset(0.4);
      h_ratio->GetYaxis()->SetLabelSize(0.17);
      h_ratio->GetYaxis()->SetTitle("Ratio");

      std::string thisName_bin = Form("%s_ratio_bin", h_first_bin->GetName());
      TH1D* h_ratio_bin = (TH1D*) h_first_bin->Clone(thisName_bin.c_str());
      h_ratio_bin->Divide(h_second);
      h_ratio_bin->Divide(h_ratio);
      h_ratio_bin->Write();
      h_ratio_bin->SetStats(0);	    
      h_ratio_bin->SetMarkerStyle(24);
      h_ratio_bin->SetMarkerColor(2);
      h_ratio_bin->SetMarkerSize(1.);
      h_ratio_bin->SetLineColor(2);
      h_ratio_bin->GetXaxis()->SetLabelSize(0.00);
      h_ratio_bin->GetXaxis()->SetTickLength(0.09);
      h_ratio_bin->GetYaxis()->SetNdivisions(5,5,0);
      h_ratio_bin->GetYaxis()->SetRangeUser(0.0,2.0);
      h_ratio_bin->GetYaxis()->SetTitleSize(0.17);
      h_ratio_bin->GetYaxis()->SetTitleOffset(0.26);
      h_ratio_bin->GetYaxis()->SetLabelSize(0.17);
      h_ratio_bin->GetYaxis()->SetTitle("Ratio");
      
            
      TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, 0.0, 2.0 );
      h2_axes_ratio->GetYaxis()->SetTitle("Bin / Std.");

      TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
      lineCentral->SetLineColor(1);

      std::string thisName_Band =  Form("%s_band", h_ratio_bin->GetName());
      TH1D* h_band = (TH1D*)h_ratio_bin->Clone(thisName_Band.c_str());
      h_band->SetMarkerSize(0);
      h_band->SetFillColor (kGray+2);
      h_band->SetFillStyle (3244);
      for ( int iBin=1; iBin <= h_ratio_bin->GetNbinsX(); iBin++){

	h_band->SetBinContent(iBin,1);
	double error;
	if(h_first->GetBinContent(iBin)>0)
	  error = h_first->GetBinError(iBin)/h_first->GetBinContent(iBin);
	else
	  error=0.;
	h_band->SetBinError(iBin, error);

      }

      h2_axes_ratio->Draw("");
      h_band->Draw("e2,same");
      lineCentral->Draw("same");
      //      h_ratio->Draw("pe,same");
      h_ratio_bin->Draw("pe,same");

      gPad->RedrawAxis();
      
      c1->cd();

      c1->SaveAs( Form("%s/mt2_%s.pdf", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.png", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.eps", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;
      delete h2_axes;
      delete h_ratio;
      delete h_ratio_bin;
      delete h_first;
      delete h_first_bin;
      delete h_second;
      
      //      ++iRegion;

  } // for MT2 regions

  
  TCanvas* c2 = new TCanvas("c2", "", 1200, 600);
  c2->cd();
  
  TH1D* hHelp = (TH1D*) hdata->Clone("hHelp");
  for (int iBin=1; iBin <= hHelp->GetNbinsX(); ++iBin )
    hHelp->SetBinError(iBin, 0);
  
  std::string thisName = Form("%s_ratio", hdata->GetName());
  TH1D* h_Ratio = (TH1D*) hHelp->Clone(thisName.c_str());
  h_Ratio->Divide(hestimate);
  h_Ratio->Write();
  h_Ratio->SetStats(0);
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(4);
  h_Ratio->GetXaxis()->SetLabelSize(0.00);
  h_Ratio->GetXaxis()->SetTickLength(0.09);
  h_Ratio->GetYaxis()->SetNdivisions(5,5,0);
  h_Ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  h_Ratio->GetYaxis()->SetTitleSize(0.17);
  h_Ratio->GetYaxis()->SetTitleOffset(0.4);
  h_Ratio->GetYaxis()->SetLabelSize(0.17);
  h_Ratio->GetYaxis()->SetTitle("Ratio");

  std::string thisName_bin = Form("%s_ratio_bin", hdata_bin->GetName());
  TH1D* h_Ratio_bin = (TH1D*) hdata_bin->Clone(thisName_bin.c_str());
  h_Ratio_bin->Divide(hestimate);
  h_Ratio_bin->Divide(h_Ratio);
  h_Ratio_bin->Write();
  h_Ratio_bin->SetStats(0);
  h_Ratio_bin->SetMarkerStyle(24);
  h_Ratio_bin->SetMarkerSize(1.0);
  h_Ratio_bin->SetLineColor(2);
  h_Ratio_bin->GetXaxis()->SetLabelSize(0.00);
  h_Ratio_bin->GetXaxis()->SetTickLength(0.09);
  h_Ratio_bin->GetYaxis()->SetNdivisions(5,5,0);
  h_Ratio_bin->GetYaxis()->SetRangeUser(0.0,2.0);
  h_Ratio_bin->GetYaxis()->SetTitleSize(0.17);
  h_Ratio_bin->GetYaxis()->SetTitleOffset(0.4);
  h_Ratio_bin->GetYaxis()->SetLabelSize(0.17);
  h_Ratio_bin->GetYaxis()->SetTitle("Ratio");

  std::string thisName_Band =  Form("%s_Band", hdata_bin->GetName());
  TH1D* h_Band = (TH1D*)hdata_bin->Clone(thisName_Band.c_str());
  h_Band->SetMarkerSize(0);
  h_Band->SetFillColor (kGray+2);
  h_Band->SetFillStyle (3244);
  for ( int iBin=1; iBin <= hdata_bin->GetNbinsX(); iBin++){

    h_Band->SetBinContent(iBin,1);
    double error;
    if(hdata->GetBinContent(iBin)>0)
      error = hdata->GetBinError(iBin)/hdata->GetBinContent(iBin);
    else error=0.;
    h_Band->SetBinError(iBin, error);

  }


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
  
  float yMin = 1e-2;
  yMax*=20.;
  
  int thisBin=174;
  hestimate->GetXaxis()->SetRangeUser(0, thisBin);
  hdata->GetXaxis()->SetRangeUser(0, thisBin);
  hdata_bin->GetXaxis()->SetRangeUser(0, thisBin);
  //  hestimate->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  hestimate->GetYaxis()->SetRangeUser(yMin, yMax);
  hdata->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate->GetXaxis()->LabelsOption("v");
  hdata->GetXaxis()->LabelsOption("v");
  hdata->GetXaxis()->SetLabelSize(0.04);
  hdata->GetXaxis()->SetLabelFont(62);
  hestimate->GetXaxis()->SetLabelSize(0.04);
  //  hestimate->Draw("PE");

  TGraphAsymmErrors* g_data = new TGraphAsymmErrors(0);
  for( int iBin=1; iBin<(hdata->GetXaxis()->GetNbins()+1); ++iBin ) {

    double y;
    double x, xerr, yerrplus, yerrminus;

    x = hdata->GetBinCenter(iBin);
    xerr = hdata->GetBinWidth(iBin)/2.;

    y = hdata->GetBinContent(iBin);
    double yerr = hdata->GetBinError(iBin);

    int thisPoint = g_data->GetN();
    g_data->SetPoint( thisPoint, x, y );
    g_data->SetPointError( thisPoint, xerr, xerr, yerr, yerr );

  }
  

  hdata->Draw("pe");
  //  hdata->Draw("pe,same");
  hdata_bin->Draw("pe,same");

  TPaveText* bkg = new TPaveText( 0.7, 0.9-0.06, 0.93, 0.9, "brNDC" );
  bkg->SetTextSize(0.042);
  //  bkg->SetTextFont(42);
  bkg->SetFillColor(0);
  bkg->SetTextAlign(11);
  bkg->AddText("Lost lepton background");
  bkg->Draw("same");


  TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06-0.06-0.06, 0.93, 0.9-0.06 );
  legend->SetTextSize(0.042);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
//  legend->AddEntry( hdata, "data-driven (standard)", "PL" );
//  legend->AddEntry( hdata_bin, "data-driven (bin by bin)", "PL" );
  legend->AddEntry( hdata, "Standard", "PL" );
  legend->AddEntry( hdata_bin, "Bin by bin", "PL" );
  //  legend->AddEntry( hestimate, "MC", "PL" );

  legend->Draw("same");

  //  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
  labelTop->Draw("same");
  
  //  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Supplementary");
  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Preliminary");
  labelCMS->Draw("same");
  
  int nHTRegions = 6;
  std::vector< std::string > htRegions;
  htRegions.push_back("1 Jet");
  htRegions.push_back("H_{T} [200, 450] GeV");
  htRegions.push_back("H_{T} [450, 575] GeV");
  htRegions.push_back("H_{T} [575, 1000] GeV");
  htRegions.push_back("H_{T} [1000, 1500] GeV");
  htRegions.push_back("H_{T} > 1500 GeV");
//  htRegions.push_back("very low H_{T}");
//  htRegions.push_back("low H_{T}");
//  htRegions.push_back("medium H_{T}");
//  htRegions.push_back("high H_{T}");
//  htRegions.push_back("extreme H_{T}");
  
  TPaveText* htBox[5];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){
    
    htBox[iHT] = new TPaveText(0.4, 0.9-0.03, 0.7, 0.85+0.03, "brNDC");
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

  bool doLogRatio=false;
  TH2D* h2_axes_ratio;
  if(doLogRatio){
    pad2->SetLogy();
    h_Ratio_bin->GetYaxis()->SetRangeUser(0.1, 10.0);
    h_Band->GetYaxis()->SetRangeUser(0.1, 10.0);
    h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0, 2.0 );

  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitle("Bin / Std.");
  
  TLine* LineCentral = new TLine(0, 1.0, thisBin, 1.0);
  LineCentral->SetLineColor(1);

  h2_axes_ratio->Draw("");
  h_Band->Draw("e2,same");
  LineCentral->Draw("same");
  //  h_Ratio->Draw("pe,same");
  h_Ratio_bin->Draw("pe,same");

  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.pdf", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.png", fullPath.c_str()) );



  TCanvas* c2_0 = new TCanvas("c2_0", "", 1100, 600);
  c2_0->cd();
  
  TPad *pad1_0 = new TPad("pad1_0","pad1_0",0,0.3-0.1,1,1);
  pad1_0->SetBottomMargin(0.15);
  pad1_0->Draw();
  pad1_0->cd();

  pad1_0->SetLogy();
  
  yMin=1e-2;
 
  int oldBin=0;
  thisBin=12;
  hestimate->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  //  hestimate->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  hestimate->GetYaxis()->SetRangeUser(yMin, yMax);
  hdata->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate->GetXaxis()->LabelsOption("v");
  hestimate->GetXaxis()->SetLabelSize(0.04);
  hdata->GetYaxis()->SetTitleOffset(0.95);
  hdata->GetYaxis()->SetLabelSize(0.042);
  //  hestimate->Draw("PE");  

  hdata->Draw("pe");
  //  hdata->Draw("pe,same");
  hdata_bin->Draw("pe,same");

  bkg->Draw("same");
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");
  
  htBox[0]->Draw("same");

  gPad->RedrawAxis();
  
  c2_0->cd();
  TPad *pad2_0 = new TPad("pad2_0","pad2_0",0,0,1,0.21);
  pad2_0->SetTopMargin(0.05);
  pad2_0->SetBottomMargin(0.1);
  pad2_0->Draw();
  pad2_0->cd();
  
  TH2D* h2_axes_ratio_0;
    if(doLogRatio){
    pad2_0->SetLogy();
    h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0, 2.0 );

  h2_axes_ratio_0->SetStats(0);
  h2_axes_ratio_0->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_0->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_0->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_0->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_0->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_0->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_0->GetYaxis()->SetTitle("Bin / Std.");
  
  TLine* LineCentral_0 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_0->SetLineColor(1);

  h2_axes_ratio_0->Draw("");
  h_Band->Draw("e2,same");
  LineCentral_0->Draw("same");
  //  h_Ratio->Draw("pe,same");
  h_Ratio_bin->Draw("pe,same");

  gPad->RedrawAxis();

  c2_0->cd();
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.pdf", fullPath.c_str()) );
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.png", fullPath.c_str()) );



  TCanvas* c2_1 = new TCanvas("c2_1", "", 1100, 600);
  c2_1->cd();
  
  TPad *pad1_1 = new TPad("pad1_1","pad1_1",0,0.3-0.1,1,1);
  pad1_1->SetBottomMargin(0.15);
  pad1_1->Draw();
  pad1_1->cd();

  pad1_1->SetLogy();
    
  yMin=1e-2;

  oldBin=thisBin;
  thisBin=36;
  hestimate->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  //  hestimate->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  hestimate->GetYaxis()->SetRangeUser(yMin, yMax);
  hdata->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate->GetXaxis()->LabelsOption("v");
  hestimate->GetXaxis()->SetLabelSize(0.04);
  //  hestimate->Draw("PE");  

  hdata->Draw("pe");
  //  hdata->Draw("pe,same");
  hdata_bin->Draw("pe,same");

  bkg->Draw("same");

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
    pad2_1->SetLogy();
    h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0, 2.0 );

  h2_axes_ratio_1->SetStats(0);
  h2_axes_ratio_1->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_1->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_1->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_1->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_1->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_1->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_1->GetYaxis()->SetTitle("Bin / Std.");
  
  TLine* LineCentral_1 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_1->SetLineColor(1);

  h2_axes_ratio_1->Draw("");
  h_Band->Draw("e2,same");
  LineCentral_1->Draw("same");
  //  h_Ratio->Draw("pe,same");
  h_Ratio_bin->Draw("pe,same");

  gPad->RedrawAxis();

  c2_1->cd();
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.png", fullPath.c_str()) );



  TCanvas* c2_2 = new TCanvas("c2_2", "", 1100, 600);
  c2_2->cd();
  
  TPad *pad1_2 = new TPad("pad1_2","pad1_2",0,0.3-0.1,1,1);
  pad1_2->SetBottomMargin(0.15);
  pad1_2->Draw();
  pad1_2->cd();

  pad1_2->SetLogy();
  
  yMin=0.5e-2;
    
  oldBin=thisBin;
  thisBin=67;
  hestimate->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  //  hestimate->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  hestimate->GetYaxis()->SetRangeUser(yMin, yMax);
  hdata->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate->GetXaxis()->LabelsOption("v");
  hestimate->GetXaxis()->SetLabelSize(0.04);
  //  hestimate->Draw("PE");  

  hdata->Draw("pe");
  //  hdata->Draw("pe,same");
  hdata_bin->Draw("pe,same");

  bkg->Draw("same");

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
    pad2_2->SetLogy();
    h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0, 2.0 );


  h2_axes_ratio_2->SetStats(0);
  h2_axes_ratio_2->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_2->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_2->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_2->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_2->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_2->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_2->GetYaxis()->SetTitle("Bin / Std.");
  
  TLine* LineCentral_2 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_2->SetLineColor(1);

  h2_axes_ratio_2->Draw("");
  h_Band->Draw("e2,same");
  LineCentral_2->Draw("same");
  //  h_Ratio->Draw("pe,same");
  h_Ratio_bin->Draw("pe,same");

  gPad->RedrawAxis();

  c2_2->cd();
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.png", fullPath.c_str()) );




  TCanvas* c2_3 = new TCanvas("c2_3", "", 1100, 600);
  c2_3->cd();
  
  TPad *pad1_3 = new TPad("pad1_3","pad1_3",0,0.3-0.1,1,1);
  pad1_3->SetBottomMargin(0.15);
  pad1_3->Draw();
  pad1_3->cd();

  pad1_3->SetLogy();

  yMin=0.5e-2;    

  oldBin=thisBin;
  thisBin=109;
  hestimate->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  //  hestimate->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  //  hestimate->GetYaxis()->SetRangeUser(yMin, yMax);
  hdata->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate->GetXaxis()->LabelsOption("v");
  hestimate->GetXaxis()->SetLabelSize(0.04);
  //  hestimate->Draw("PE");  

  hdata->Draw("pe");
  //  hdata->Draw("pe,same");
  hdata_bin->Draw("pe,same");

  bkg->Draw("same");

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
    pad2_3->SetLogy();
    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0, 2.0 );

  h2_axes_ratio_3->SetStats(0);
  h2_axes_ratio_3->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_3->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_3->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_3->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_3->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_3->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_3->GetYaxis()->SetTitle("Bin / Std.");
  
  TLine* LineCentral_3 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_3->SetLineColor(1);

  h2_axes_ratio_3->Draw("");
  h_Band->Draw("e2,same");
  LineCentral_3->Draw("same");
  //  h_Ratio->Draw("pe,same");
  h_Ratio_bin->Draw("pe,same");

  gPad->RedrawAxis();

  c2_3->cd();
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.png", fullPath.c_str()) );




  TCanvas* c2_4 = new TCanvas("c2_4", "", 1100, 600);
  c2_4->cd();
  
  TPad *pad1_4 = new TPad("pad1_4","pad1_4",0,0.3-0.1,1,1);
  pad1_4->SetBottomMargin(0.15);
  pad1_4->Draw();
  pad1_4->cd();

  pad1_4->SetLogy();
    
  yMin=1e-2;

  oldBin=thisBin;
  thisBin=144;
  hestimate->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  //  hestimate->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  //  hestimate->GetYaxis()->SetRangeUser(yMin, yMax);
  hdata->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate->GetXaxis()->LabelsOption("v");
  hestimate->GetXaxis()->SetLabelSize(0.04);
  //  hestimate->Draw("PE");  

  hdata->Draw("pe");
  //  hdata->Draw("pe,same");
  hdata_bin->Draw("pe,same");

  bkg->Draw("same");

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
    pad2_4->SetLogy();
    h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0, 2.0 );

  h2_axes_ratio_4->SetStats(0);
  h2_axes_ratio_4->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_4->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_4->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_4->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_4->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_4->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_4->GetYaxis()->SetTitle("Bin / Std.");
  
  TLine* LineCentral_4 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_4->SetLineColor(1);

  h2_axes_ratio_4->Draw("");
  h_Band->Draw("e2,same");
  LineCentral_4->Draw("same");
  //  h_Ratio->Draw("pe,same");
  h_Ratio_bin->Draw("pe,same");

  gPad->RedrawAxis();

  c2_4->cd();
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.png", fullPath.c_str()) );




  TCanvas* c2_5 = new TCanvas("c2_5", "", 1100, 600);
  c2_5->cd();
  
  TPad *pad1_5 = new TPad("pad1_5","pad1_5",0,0.3-0.1,1,1);
  pad1_5->SetBottomMargin(0.15);
  pad1_5->Draw();
  pad1_5->cd();

  pad1_5->SetLogy();

  yMin=0.5e-3;

  oldBin=thisBin;
  thisBin=174;
  hestimate->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio_bin->GetXaxis()->SetRangeUser(oldBin, thisBin);
  //  hestimate->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  //  hestimate->GetYaxis()->SetRangeUser(yMin, yMax);
  hdata->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate->GetXaxis()->LabelsOption("v");
  hestimate->GetXaxis()->SetLabelSize(0.04);
  //  hestimate->Draw("PE");  

  //  hdata->Draw("pe,same");
  hdata->Draw("pe");
  hdata_bin->Draw("pe,same");

  bkg->Draw("same");

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
    pad2_5->SetLogy();
    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0, 2.0 );

  h2_axes_ratio_5->SetStats(0);
  h2_axes_ratio_5->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_5->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_5->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_5->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_5->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_5->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_5->GetYaxis()->SetTitle("Bin / Std.");
  
  TLine* LineCentral_5 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_5->SetLineColor(1);

  h2_axes_ratio_5->Draw("");
  h_Band->Draw("e2,same");
  LineCentral_5->Draw("same");
  //  h_Ratio->Draw("pe,same");
  h_Ratio_bin->Draw("pe,same");

  gPad->RedrawAxis();

  c2_5->cd();
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.png", fullPath.c_str()) );




  gStyle->SetOptStat(1110);
  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  hPull->SetStats(1110);
  hPull->Draw("hist");
  std::cout << hPull->GetEntries() << "\t" << hPull->GetMean() << "\t" << hPull->GetRMS() << std::endl; 
  c3->SaveAs( Form("%s/PullDistribution.pdf", fullPath.c_str()) );
  c3->SaveAs( Form("%s/PullDistribution.png", fullPath.c_str()) );
  
}
