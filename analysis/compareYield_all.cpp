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

float lumi =0.134; //fb-1 

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::vector<MT2Analysis<MT2Estimate>* > bgYields );

int main() {

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = "./YieldComparison_dataMC_binned/";

  
  std::string llepInputFile = "./EventYields_data_Run2015D_25nsGolden_zurichPlus/llepEstimate.root";
  std::string zinvInputFile = "./EventYields_data_Run2015D_25nsGolden_zurichPlus/zinvFromGamma.root";
  std::string qcdInputFile  = "./EventYields_data_Run2015D_25nsGolden_zurichPlus/qcdEstimateData.root";

  std::string mcInputFile = "./EventYields_data_Run2015D_25nsGolden_zurichPlus/analyses.root";

  std::string dataInputFile = "./EventYields_data_Run2015D_25nsGolden_zurichPlus_149ipb/analyses.root";
  
  MT2Analysis<MT2Estimate>* data = MT2Analysis<MT2Estimate>::readFromFile( dataInputFile.c_str(), "data" );
  
  MT2Analysis<MT2Estimate>* llep = MT2Analysis<MT2Estimate>::readFromFile( llepInputFile.c_str(), "llepEstimate" );
  llep->setName("Lost Lepton");
  
  MT2Analysis<MT2Estimate>* zinv = MT2Analysis<MT2Estimate>::readFromFile( zinvInputFile.c_str(), "ZinvEstimate" );
  zinv->setName("Invisible Z");

  MT2Analysis<MT2Estimate>* qcd  = MT2Analysis<MT2Estimate>::readFromFile( qcdInputFile.c_str(), "qcdEstimate"  );
  qcd->setName("Multijet");

  MT2Analysis<MT2Estimate>* qcd_mc  = MT2Analysis<MT2Estimate>::readFromFile( mcInputFile.c_str(), "QCD"  );

//  MT2Analysis<MT2Estimate>* top = MT2Analysis<MT2Estimate>::readFromFile( dataInputFile.c_str(), "Top" );
//  MT2Analysis<MT2Estimate>* wjets = MT2Analysis<MT2Estimate>::readFromFile( dataInputFile.c_str(), "WJets" );
//  MT2Analysis<MT2Estimate>* zinv = MT2Analysis<MT2Estimate>::readFromFile( dataInputFile.c_str(), "ZJets" );
//  MT2Analysis<MT2Estimate>* qcd  = MT2Analysis<MT2Estimate>::readFromFile( dataInputFile.c_str(), "QCD"  );
  
  std::vector< MT2Analysis<MT2Estimate> *> bgEstimate;
//  bgEstimate.push_back( qcd  );
//  bgEstimate.push_back( top );
//  bgEstimate.push_back( wjets );
//  bgEstimate.push_back( zinv );
  bgEstimate.push_back( qcd  );
  bgEstimate.push_back( llep );
  bgEstimate.push_back( zinv );
  bgEstimate.push_back( qcd_mc );

  drawYields(outputdir.c_str(), data, bgEstimate );

  return 0;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::vector< MT2Analysis<MT2Estimate> *> bgYields ) {

  MT2DrawTools::setStyle();

  system(Form("mkdir -p %s", outputdir.c_str()));

  std::vector<int> colors;
//  colors.push_back( 401 );
//  colors.push_back( 855 );
//  colors.push_back( 417 );
//  colors.push_back( 419 );
  colors.push_back( 402 );
  colors.push_back( 430 );
  colors.push_back( 418 );
  
  std::set<MT2Region> MT2Regions = data->getRegions();
  
  TH1D* hdata = new TH1D("hdata", "", 185, 0, 185);
  //  TH1D* hdata = new TH1D("hdata", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  hdata->Sumw2();
  hdata->GetYaxis()->SetTitle("Entries");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( 1 );
  hdata->SetMarkerColor( 1 );
  
  TH1D* hestimate_all;
  TH1D* hestimate[bgYields.size()];
  for(unsigned int b=0; b<bgYields.size(); ++b){
  
    hestimate[b]= new TH1D(Form("hestimate_%d", b), "", 185, 0, 185);
    hestimate[b]->Sumw2();
    hestimate[b]->GetYaxis()->SetTitle("Entries");
    hestimate[b]->SetFillColor(colors[b]);
    hestimate[b]->SetLineColor(1);
    
  }

  THStack bgStack("bgStack", "");

//  TH1D* hPull = new TH1D("hPull", "", 20, -5, 5);
//  hPull->Sumw2();
//  hPull->GetXaxis()->SetTitle("(Est. - Obs.)/#sigma");
//  hPull->GetYaxis()->SetTitle("Entries");
  
  //std::string thisName = Form("%s_ratio", hdata->GetName());
  //TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
  //h_Ratio->Sumw2();
  //h_Ratio->GetYaxis()->SetTitle("(Data Driven - MC)/#sigma");
  
  std::string fullPath = outputdir;
  
  int iRegion = 14;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      TH1D* h_first = data->get(*iMT2)->yield;
      //TGraphAsymmErrors* g_first = new TGraphAsymmErrors(0); // = new TGraphAsymmErrors(h_first);
      TGraphAsymmErrors* g_first = MT2DrawTools::getPoissonGraph(h_first);      
      
      int nBins = h_first->GetXaxis()->GetNbins()+1;
      
      TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), iMT2->getName().c_str()), "recreate" );
      histoFile->cd();
      h_first->Write();
      
      g_first->SetMarkerStyle(20);
      g_first->SetMarkerSize(1.6);
      g_first->SetLineColor( 1 );
      g_first->SetMarkerColor( 1 );

      THStack bgStack_region("bgStack_region", "");
      
      TH1D* h_second_all;
      TH1D* h_second[bgYields.size()-1];

      for(unsigned int b=0; b<bgYields.size()-1; ++b){

	if(b==0 && iMT2->nJetsMin()==1 )
	  h_second[b] = bgYields[bgYields.size()-1]->get(*iMT2)->yield;
	else
	  h_second[b] = bgYields[b]->get(*iMT2)->yield;
	h_second[b]->SetFillColor( colors[b] );
	h_second[b]->SetLineColor( 1 );

	h_second[b]->Scale(lumi/1.25);
	
	std::cout << "Read yield for background " << b << ", with integral: " << h_second[b]->Integral() << std::endl;
	
      	for( int iBin=1; iBin<(h_second[b]->GetXaxis()->GetNbins()+1); ++iBin ) {
	  
	  double y;                                                                                                                                                                                                            
	  double x, xerr, yerrplus, yerrminus;
	  
	  x = h_second[b]->GetBinCenter(iBin);
	  xerr = h_second[b]->GetBinWidth(iBin)/2.;
	  
	  y = fabs(h_second[b]->GetBinContent(iBin));
	  float yerr = h_second[b]->GetBinError(iBin);
	  float yerr2=1.0;
	  if(y>0) yerr2 = yerr*yerr/(y*y);
	  
	  if(b==0 && iMT2->nJetsMin()==1){

	    yerr=y;

	  } else if(b==1){

	    yerr2+=0.05*0.05;
	    if( h_second[b]->GetXaxis()->GetNbins()>1 )
	      yerr2+=0.15*0.15;
 	    yerr2+=0.05*0.05;
	    yerr2+=0.15*0.15;
	   
	  } else if(b==2){

	    yerr2+=0.21*0.21;
	    if( h_second[b]->GetXaxis()->GetNbins()>1 )
	      yerr2+=0.15*0.15;
	    yerr2+=0.05*0.05;
	   
	  }

	  yerr = TMath::Sqrt(yerr2)*y;
	  h_second[b]->SetBinContent(iBin, y);
	  h_second[b]->SetBinError(iBin, yerr);
	  
	}
	
	std::cout << "Adding to stack histogram for background " << b << " with integral: " << h_second[b]->Integral() << std::endl;
	bgStack_region.Add(h_second[b]);
	
	if(b==0) h_second_all = (TH1D*) h_second[b]->Clone("h_second_all");
	else h_second_all->Add(h_second[b]);

	std::cout << "Summed histogram for background " << b << ", with intermediate integral: "  << h_second_all->Integral() << std::endl;

      }
      
      int firstBin=1;
      double err_data;
      double int_data;
      if( iMT2->nJetsMax()==1 ){
	iRegion=iRegion-24;
	for (int iBin=1; iBin<nBins; ++iBin){
	  
	  int_data = h_first->GetBinContent(iBin);
	  err_data = h_first->GetBinError(iBin);
	  
	  hdata->SetBinContent(iRegion, int_data);
	  hdata->SetBinError(iRegion, err_data);
	  
	  std::cout<<"Filled data  with: " << int_data << std::endl;
	  
	  hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
	  
	  for(unsigned int b=0; b<bgYields.size()-1; ++b){
	    
	    double err_int = h_second[b]->GetBinError(iBin);
	    double integral = h_second[b]->GetBinContent(iBin);
	    hestimate[b]->SetBinContent(iRegion, integral);
	    hestimate[b]->SetBinError(iRegion, err_int);
	    
	    std::cout<<"Filled estimate for background " << b << " with: " << integral << std::endl;
	    
	    std::string thisLabel=Form("%s,%d", niceNames[1].c_str(), iBin); 
	    hestimate[b]->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
	    
	  }
	  
	  std::cout <<"iRegion " << iRegion << std::endl;
	  
	  if(iBin<nBins-1)
	    ++iRegion;
	  
	}
	if(iMT2->nBJetsMax()==0)
          iRegion=iRegion+24;
        else
          iRegion=iRegion+11;
      }
      else{
	
	int_data = h_first->IntegralAndError(firstBin, nBins, err_data);

	hdata->SetBinContent(iRegion, int_data);
	hdata->SetBinError(iRegion, err_data);

	std::cout<<"Filled data  with: " << int_data << std::endl;

	hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );

	for(unsigned int b=0; b<bgYields.size()-1; ++b){

	  double err_int;
	  double integral = h_second[b]->IntegralAndError(firstBin, nBins, err_int);
	  hestimate[b]->SetBinContent(iRegion, integral);
	  hestimate[b]->SetBinError(iRegion, err_int);

	  std::cout<<"Filled estimate for background " << b << " with: " << integral << std::endl;

	  std::string thisLabel=Form("%s", niceNames[1].c_str());
	  hestimate[b]->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );

	}
	
	
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
      
      std::cout << "yMax = " << yMax1 << std::endl;	
      
      float yMax_3 = h_second_all->GetMaximum()*1.5;
      float yMax_4 = 1.2*(h_second_all->GetMaximum() + h_second_all->GetBinError(h_second_all->GetMaximumBin()));
      float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
      float yMax = (yMax1>yMax2) ? yMax1 : yMax2;

      std::cout << "yMax = " << yMax << std::endl;
      std::cout << h_second_all->GetMaximum() << "\t" << h_second_all->GetBinError(h_second_all->GetMaximumBin()) << std::endl;
      
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


      TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+1-1)*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
      legend->AddEntry( g_first, "Data", "P" );

      legend->Draw("same");

      bgStack_region.Draw("histo, same");
      h_first->Draw("pe,same");
      
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
      h_ratio->Divide(h_second_all);
      h_ratio->Write();
      h_ratio->SetStats(0);	    
      h_ratio->SetMarkerStyle(20);
      h_ratio->SetMarkerColor(1);
      h_ratio->SetLineColor(1);
      h_ratio->GetXaxis()->SetLabelSize(0.00);
      h_ratio->GetXaxis()->SetTickLength(0.09);
      h_ratio->GetYaxis()->SetNdivisions(5,5,0);
      h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
      h_ratio->GetYaxis()->SetTitleSize(0.17);
      h_ratio->GetYaxis()->SetTitleOffset(0.4);
      h_ratio->GetYaxis()->SetLabelSize(0.17);
      h_ratio->GetYaxis()->SetTitle("Ratio");
            
      TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, 0.0, 2.0 );
      h2_axes_ratio->GetYaxis()->SetTitle("Ratio");
      
      TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
      lineCentral->SetLineColor(1);

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      h_ratio->Draw("pe,same");
      
      gPad->RedrawAxis();
      
      c1->cd();

      c1->SaveAs( Form("%s/mt2_%s.pdf", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.png", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.eps", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;

      delete h2_axes;
      delete h2_axes_ratio;

      delete h_ratio;

      delete h_first;

      delete h_second_all;

      for(unsigned int b=0; b<bgYields.size()-1; ++b)
	delete h_second[b];
      
      ++iRegion;

  } // for MT2 regions


  for(unsigned int b=0; b<bgYields.size()-1; ++b){

    bgStack.Add(hestimate[b]);
    
    if(b==0) hestimate_all = (TH1D*) hestimate[b]->Clone("hestimate_all");
    else hestimate_all->Add(hestimate[b]);

  }

  
  TCanvas* c2 = new TCanvas("c2", "", 1200, 600);
  c2->cd();
  

  std::string thisName = Form("%s_ratio", hdata->GetName());
  TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
  h_Ratio->Divide(hestimate_all);
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
  
  float yMin = 1e-2;
  yMax*=20.;
  
  hestimate_all->GetXaxis()->SetRangeUser(0, 68);
  hdata->GetXaxis()->SetRangeUser(0, 68);
  //  hestimate_all->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());  
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hdata->GetYaxis()->SetRangeUser(yMin, yMax);
  hdata->GetXaxis()->LabelsOption("v");
  //hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->SetFillStyle(3244);
  hestimate_all->SetFillColor(kGray+2);
  //  hestimate_all->Draw("PE");

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
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  hdata->Draw("pe,same");
  
  TLegend* legend = new TLegend( 0.8, 0.85-(bgYields.size()+1-1)*0.06-0.06, 0.93, 0.85-0.06 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( hdata, "Data", "PL" );
  for(unsigned int b=0; b<bgYields.size()-1; ++b)
    legend->AddEntry( hestimate[b], bgYields[b]->getName().c_str(), "F" );

  legend->Draw("same");

  //  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
  labelTop->Draw("same");
  

  TLine* lHT[5];
  for( int iHT=0; iHT < 5; iHT++ ){
    lHT[iHT-1] = new TLine(13+11*iHT, 0.0, 13+11*iHT, yMax );
    lHT[iHT-1]->SetLineColor(kBlack);
    lHT[iHT-1]->SetLineStyle(3);
    lHT[iHT-1]->SetLineWidth(2);

    lHT[iHT-1]->Draw("same");
  }

  int nHTRegions = 6;
  std::vector< std::string > htRegions;
  htRegions.push_back("1 Jet");
  htRegions.push_back("very low H_{T}");
  htRegions.push_back("low H_{T}");
  htRegions.push_back("medium H_{T}");
  htRegions.push_back("high H_{T}");
  htRegions.push_back("extreme H_{T}");
  
  TPaveText* htBox[5];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){

    if (iHT==0) htBox[iHT] = new TPaveText(0.12+0.15*iHT, 0.9-0.06, 0.34+0.15*iHT, 0.85, "brNDC");
    else htBox[iHT] = new TPaveText(0.16+0.13*iHT, 0.9-0.06, 0.34+0.13*iHT, 0.85, "brNDC");
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

  TH2D* h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, 68, 10, 0., 2.0 );
  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.3);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitle("Ratio");
  
  TLine* LineCentral = new TLine(0, 1.0, 68, 1.0);
  LineCentral->SetLineColor(1);

  std::string thisName_Band =  Form("%s_band", hestimate_all->GetName());
  TH1D* h_band = (TH1D*)hestimate_all->Clone(thisName_Band.c_str());
  h_band->SetMarkerSize(0);
  h_band->SetFillColor (kGray+2);
  h_band->SetFillStyle (3244);
  for ( int iBin=1; iBin <= hestimate_all->GetNbinsX(); iBin++){
    
    h_band->SetBinContent(iBin,1);
    double error;
    error = hestimate_all->GetBinError(iBin)/hestimate_all->GetBinContent(iBin);
    h_band->SetBinError(iBin, error);

  }


  h2_axes_ratio->Draw("");
  h_band->Draw("E2same");
  LineCentral->Draw("same");
  h_Ratio->Draw("pe,same");
  

  TLine* lHT_b[6];
  for( int iHT=1; iHT < 6; iHT++ ){
    lHT_b[iHT-1] = new TLine(13+11*(iHT-1), 0, 13+11*(iHT-1), 2.0 );

    lHT_b[iHT-1]->SetLineColor(kBlack);
    lHT_b[iHT-1]->SetLineStyle(3);
    lHT_b[iHT-1]->SetLineWidth(2);

    lHT_b[iHT-1]->Draw("same");
  }
  
  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/mt2_fullEstimate.pdf", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_fullEstimate.png", fullPath.c_str()) );
  
}
