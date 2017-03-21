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


struct BGTable {

  float zinv;
  float zinv_statUp;
  float zinv_statDn;
  float zinv_systUp;
  float zinv_systDn;

  float llep;
  float llep_statUp;
  float llep_statDn;
  float llep_systUp;
  float llep_systDn;

  float qcd;
  float qcd_statUp;
  float qcd_statDn;
  float qcd_systUp;
  float qcd_systDn;

};


float lumi; //fb-1 

BGTable getTable( const std::string& tableFileName );
void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::string dir );


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
  
  
  if( argc!=2 ) {
    std::cout << "USAGE: ./computeLostLepton [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  
  
  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  // lumi = 18.1;
  lumi = cfg.lumi();
  
  TH1::AddDirectory(kTRUE);
  
  std::string dir = cfg.getEventYieldDir();
  std::string outputdir = cfg.getEventYieldDir() + "/YieldComparison_dataMC_binned";
 
 
  MT2Analysis<MT2Estimate>* analysis = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "data" ); // any one is good, just need to know the regions                                                                    

  std::vector < MT2Analysis<MT2Estimate>* > analysesSignal;
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1tttt_1500_100") );
//  analysesSignal[0]->setName("T1tttt 1500,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1tttt_1200_800") );
//  analysesSignal[1]->setName("T1tttt 1200,800");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1bbbb_1500_100") );
//  analysesSignal[2]->setName("T1bbbb 1500,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1bbbb_1000_900") );
//  analysesSignal[3]->setName("T1bbbb 1000,900");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1qqqq_1400_100") );
//  analysesSignal[4]->setName("T1qqqq 1400,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1qqqq_1000_800") );
//  analysesSignal[5]->setName("T1qqqq 1000,800");

  std::set<MT2Region> regions = analysis->getRegions();

  MT2Analysis<MT2Estimate>* data = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "data" );
  
  drawYields( outputdir.c_str(), data, dir );

  return 0;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::string dir ) {

  
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
  
  std::set<MT2Region> MT2Regions = data->getRegions();
  
  TH1D* hdata = new TH1D("hdata", "", 213, 0, 213);
  hdata->Sumw2();
  hdata->GetYaxis()->SetTitle("Entries");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( 1 );
  hdata->SetMarkerColor( 1 );

  TH1D* hestimate_all;
  TH1D* hestimate[bgSize];
  TH1D* hestimate_all_forRatio;
  TH1D* hestimate_forRatio[bgSize];
  
  for(unsigned int b=0; b<bgSize; ++b){
  
    hestimate[b]= new TH1D(Form("hestimate_%d", b), "", 213, 0, 213);
    hestimate[b]->Sumw2();
    hestimate[b]->GetYaxis()->SetTitle("Entries");
    hestimate[b]->SetFillColor(colors[b]);
    hestimate[b]->SetLineColor(1);

    hestimate_forRatio[b]= new TH1D(Form("hestimate_forRatio%d", b), "", 213, 0, 213);
    hestimate_forRatio[b]->Sumw2();
    hestimate_forRatio[b]->GetYaxis()->SetTitle("Entries");
    hestimate_forRatio[b]->SetFillColor(colors[b]);
    hestimate_forRatio[b]->SetLineColor(1);
    
  }

  THStack bgStack("bgStack", "");

  TH1D* hPull = new TH1D("hPull", "", 101, -5.05, 5.05);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Data - Est.)/#sigma");
  hPull->GetYaxis()->SetTitle("Entries");
  
  std::string fullPath = outputdir;
  
  std::string labelsMono[12]={"[250,350]","[350,450]","[450,575]","[575,700]","[700,1000]","[1000,1200]", ">1200","[250,350]","[350,450]","[450,575]","[575,700]", ">700"};
  
  TFile* bigHistoFile = TFile::Open( Form("%s/histograms_ALL.root", fullPath.c_str()), "recreate" );

  int nBins_[63];

  int iRegion = 1;
  int iTR = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      //      int nBins;
      double *bins;
      iMT2->getBins(nBins_[iTR-1], bins);
      int nBins = nBins_[iTR-1];
      std::cout << nBins << std::endl;
      
      //      TH1D* h_first = data->get(*iMT2)->yield;
      TH1D* h_first_forExtreme = data->get(*iMT2)->yield;

      TH1D* h_first;

      if( iMT2->htMin()==1500 && iMT2->nJetsMin()>1 ){
	double *binsExtreme = bins++;
	h_first = new TH1D("h_first", "", nBins-1, binsExtreme);
	
	for( int iBin=0; iBin<nBins; ++iBin )
	  h_first->SetBinContent( iBin, h_first_forExtreme->GetBinContent(iBin+1) );

      }else 
	h_first = (TH1D*)h_first_forExtreme->Clone("h_first");


      TGraphAsymmErrors* g_first = MT2DrawTools::getPoissonGraph(h_first);    

      //if( iMT2->htMin()==1500 && iMT2->nJetsMin()>1 )
      //	g_first->GetXaxis()->SetLimits(400, 2000);  
      
      TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), iMT2->getName().c_str()), "recreate" );
      histoFile->cd();
      h_first->Write();
      
      g_first->SetMarkerStyle(20);
      g_first->SetMarkerSize(1.6);
      g_first->SetLineColor( 1 );
      g_first->SetMarkerColor( 1 );

      THStack bgStack_region("bgStack_region", "");
      
      TH1D* h_second_all;
      TH1D* h_second_forRatio_all;
      TH1D* h_second[bgSize];
      TH1D* h_second_forRatio[bgSize];

      for(unsigned int b=0; b< bgSize; ++b){

	if( iMT2->htMin()==1500 && iMT2->nJetsMin()>1 && b==0 )
	  nBins--;
	
	h_second[b] = new TH1D(Form("h_second_%d", b), "", nBins, bins);
	
	h_second_forRatio[b] = new TH1D(Form("h_second_%d", b), "", nBins, bins);
	
	h_second[b]->SetFillColor( colors[b] );
	h_second[b]->SetLineColor( 1 );
	
      }
      
      for( int iBin=0; iBin<nBins; ++iBin ) {

	if( iMT2->htMin()==1500 && iMT2->nJetsMin()>1 && bins[iBin]==200 ) continue;

	std::string tableName;
	if(iMT2->nJetsMax()==1){
          tableName = std::string(Form("%s/datacard_templates/table_%s_m0toInf.txt", dir.c_str(), iMT2->getName().c_str() ));
	}
        else{
          if( iBin < nBins-1 ){
	    tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lfto%.0lf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin], bins[iBin+1]) );
          }else
            tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lftoInf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin] ));
        }
//	if( iBin < nBins-1 )
//	  tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lfto%.0lf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin], bins[iBin+1]) );
//	else
//	  tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lftoInf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin] ));

	BGTable thisTable = getTable(tableName);
	
	float totalErr=0.;
	float statErr=0.;
	float systErr=0.;

	//QCD
	statErr = (thisTable.qcd_statUp > thisTable.qcd_statDn) ? thisTable.qcd_statUp : thisTable.qcd_statDn;
	systErr = (thisTable.qcd_systUp > thisTable.qcd_systDn) ? thisTable.qcd_systUp : thisTable.qcd_systDn;
	totalErr = TMath::Sqrt( statErr*statErr + systErr*systErr );
	h_second[0]->SetBinContent(iBin+1, thisTable.qcd);
	h_second[0]->SetBinError(iBin+1, totalErr);
	
	h_second_forRatio[0]->SetBinContent(iBin+1, thisTable.qcd);
	h_second_forRatio[0]->SetBinError(iBin+1, 0);
	
	//Lost Lepton
	statErr = (thisTable.llep_statUp > thisTable.llep_statDn) ? thisTable.llep_statUp : thisTable.llep_statDn;
	systErr = (thisTable.llep_systUp > thisTable.llep_systDn) ? thisTable.llep_systUp : thisTable.llep_systDn;
	totalErr = TMath::Sqrt( statErr*statErr + systErr*systErr );
	h_second[1]->SetBinContent(iBin+1, thisTable.llep);
	h_second[1]->SetBinError(iBin+1, totalErr);

	h_second_forRatio[1]->SetBinContent(iBin+1, thisTable.llep);
	h_second_forRatio[1]->SetBinError(iBin+1, 0);

	//Invisible Z
	statErr = (thisTable.zinv_statUp > thisTable.zinv_statDn) ? thisTable.zinv_statUp : thisTable.zinv_statDn;
	systErr = (thisTable.zinv_systUp > thisTable.zinv_systDn) ? thisTable.zinv_systUp : thisTable.zinv_systDn;
	totalErr = TMath::Sqrt( statErr*statErr + systErr*systErr );
	h_second[2]->SetBinContent(iBin+1, thisTable.zinv);
	h_second[2]->SetBinError(iBin+1, totalErr);

	h_second_forRatio[2]->SetBinContent(iBin+1, thisTable.zinv);
	h_second_forRatio[2]->SetBinError(iBin+1, 0);

      }	
	
      for(unsigned int b=0; b<bgSize; ++b){
      
	bgStack_region.Add(h_second[b]);
	
	if(b==0) h_second_all = (TH1D*) h_second[b]->Clone("h_second_all");
	else h_second_all->Add(h_second[b]);
	
	if(b==0) h_second_forRatio_all = (TH1D*) h_second_forRatio[b]->Clone("h_second_forRatio_all");
	else h_second_forRatio_all->Add(h_second_forRatio[b]);
      
      }
      
//      for(int iBin=1; iBin<=nBins; ++iBin){
//
//	float thisData    = h_first->GetBinContent(iBin);
//	float thisDataErr = h_first->GetBinError(iBin);
//
//	float thisEst     = h_second_all->GetBinContent(iBin);
//	float thisEstErr  = h_second_all->GetBinError(iBin);
//	
//
//	hPull->Fill( (thisEst-thisData)/( TMath::Sqrt( thisDataErr*thisDataErr + thisEstErr*thisEstErr ) ) );
//	
//      }

      
      double err_data;
      double int_data;
      for (int iBin=1; iBin<=nBins; ++iBin){
	
	if( iMT2->htMin()==1500 && iMT2->nJetsMin()>1 && bins[iBin]==200 ) continue;
 
	int_data = h_first->GetBinContent(iBin);
	err_data = h_first->GetBinError(iBin);
	
	hdata->SetBinContent(iRegion, int_data);
	//	hdata->SetBinError(iRegion, err_data);

	hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
	
	for(unsigned int b=0; b<bgSize; ++b){

	  double err_int = fabs(h_second[b]->GetBinError(iBin));
	  double integral = fabs(h_second[b]->GetBinContent(iBin));
	  hestimate[b]->SetBinContent(iRegion, integral);
	  hestimate[b]->SetBinError(iRegion, err_int);

	  double integral_forRatio = fabs(h_second_forRatio[b]->GetBinContent(iBin));
	  hestimate_forRatio[b]->SetBinContent(iRegion, integral_forRatio);
	  hestimate_forRatio[b]->SetBinError(iRegion, 0);
	  
	  //	  std::string thisLabel=Form("%s,%d", niceNames[1].c_str(), iBin); 

	  if( iMT2->nJetsMax()==1 )
            hestimate[b]->GetXaxis()->SetBinLabel( iRegion, labelsMono[iRegion-1].c_str() );
	  else{
	    std::string thisLabel;
	    if( iBin < nBins ){
	      if( iMT2->htMin()==1500 && iMT2->nJetsMin()>1 && bins[iBin]==200 ) continue;

	      thisLabel=Form("[%.0lf,%.0lf]", bins[iBin-1], bins[iBin]);
	    } else
	      thisLabel=Form(">%.0lf", bins[iBin-1]);
	    hestimate[b]->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
	  }
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
      
      float yMax_3 = h_second_all->GetMaximum()*1.5;
      float yMax_4 = 1.2*(h_second_all->GetMaximum() + h_second_all->GetBinError(h_second_all->GetMaximumBin()));
      float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
      float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
      
      TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
      h2_axes->SetXTitle("M_{T2} [GeV]");
      h2_axes->SetYTitle("Entries");

      if(iMT2->nJetsMax()==1)
	h2_axes->SetXTitle("p_{T}^{jet_{1}} [GeV]");

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


      TLegend* legend = new TLegend( 0.6, 0.9-(bgSize)*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
      legend->AddEntry( g_first, "Data", "P" );

      legend->Draw("same");

      bgStack_region.Draw("histo, same");
      g_first->Draw("pe0,same");
      
      //TPaveText* labelTop = MT2DrawTools::getLabelTop("18.1 fb^{-1} (13TeV)" );
      TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
      labelTop->Draw("same");
      TPaveText* labelCMS = MT2DrawTools::getLabelCMS();
      labelCMS->Draw("same");

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
      h_ratio->Draw("pe0,same");
      
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

      for(unsigned int b=0; b<bgSize; ++b)
	delete h_second[b];
      
      //      ++iRegion;

      ++iTR;

  } // for MT2 regions


  for(unsigned int b=0; b<bgSize; ++b){

    hestimate[b]->SetLineWidth(0);
    bgStack.Add(hestimate[b]);
    //bgStack.Add(hestimate_forRatio[b]);
    
    if(b==0) hestimate_all = (TH1D*) hestimate[b]->Clone("hestimate_all");
    else hestimate_all->Add(hestimate[b]);

    if(b==0) hestimate_all_forRatio = (TH1D*) hestimate_forRatio[b]->Clone("hestimate_all_forRatio");
    else hestimate_all_forRatio->Add(hestimate_forRatio[b]);

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
    
    
    hPull->Fill( (-thisEst+thisData)/( TMath::Sqrt( thisDataErr*thisDataErr + thisEstErr*thisEstErr ) ) );
    
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
  pad1->SetBottomMargin(0.18);
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
  
  float yMin = 1e-1; // float yMin = 1e-3;
  //  yMin=0;
  yMax*=20.;
  
  int thisBin=213;

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
//  hdata->Draw("pe0,same");

  gdata->Draw("pe0");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  gdata_zero->Draw("pe0,same");
  gdata->Draw("pe0,same");
  
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
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
  labelTop->Draw("same");
  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Preliminary");
  labelCMS->Draw("same");
  
  int nHTRegions = 6;
  std::vector< std::string > htRegions;
  htRegions.push_back("Monojet Region");
  htRegions.push_back("H_{T} [250, 450] GeV");
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
  pad2->SetTopMargin(0.06);
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
  h_Ratio->Draw("pe0,same");
  
  gPad->RedrawAxis();

  c2->cd();
  //  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.pdf", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.png", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.eps", fullPath.c_str()) );

  bigHistoFile->cd();
  h_band->Write();
  hestimate_all->Write();
  hestimate_all_forRatio->Write();
  hdata->Write();
  for(unsigned int b=0; b<bgSize; ++b)
    hestimate[b]->Write();


  TCanvas* c2_0 = new TCanvas("c2_0", "", 1100, 600);
  c2_0->cd();
  
  TPad *pad1_0 = new TPad("pad1_0","pad1_0",0,0.3-0.1,1,1);
  pad1_0->SetBottomMargin(0.18);
  pad1_0->Draw();
  pad1_0->cd();

  yMin = 1e-1;

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
  //hdata->Draw("pe0,same");
  gdata_zero->Draw("pe0,same");
  gdata->Draw("pe0,same");

  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");
  
  htBox[0]->Draw("same");
  
  float left = pad1_0->GetLeftMargin();
  float right = pad1_0->GetRightMargin();
  float bot = pad1_0->GetBottomMargin();
  float top = pad1_0->GetTopMargin();
  float binWidth = (1.0-right-left)/thisBin;
  
  TLatex* text = new TLatex();
  text->SetNDC(1);
  
  // draw the "Pre-fit background" text
  text->SetTextAlign(13);
  text->SetTextFont(42);
  text->SetTextAngle(0);
  text->SetTextSize(0.05);
  text->DrawLatex(left+0.04,1-top-0.01, "Pre-fit background");


  float ibin = 0;
  int monoBin[2]={7,5};
  TString monoJ[2] = {"1j", "1j"};
  TString monoB[2] = {"0b", "#geq 1b"};
  float xcenter;
  for(int nR=0; nR<2; nR++){ 
    
    xcenter = left+binWidth*(ibin+(monoBin[nR]-1)*0.5);
    text->SetTextAlign(23);
    text->SetTextFont(62);
    text->SetTextSize(0.030);
    
    float y=bot+(1-top-bot)*0.85;
    if (xcenter>1-right-0.19)
      y=0.67;
    
    text->DrawLatex(xcenter, y, monoJ[nR]);
    text->DrawLatex(xcenter,y-text->GetTextSize()-0.001,monoB[nR]);

    ibin+=monoBin[nR];

  }
  
  TLine* line = new TLine();
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  float x=left+monoBin[0]*binWidth;
  line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.85);
  std::cout << bot << std::endl;
  
  gPad->RedrawAxis();
  
  bool doLogRatio=false;

  c2_0->cd();
  TPad *pad2_0 = new TPad("pad2_0","pad2_0",0,0,1,0.21);
  pad2_0->SetTopMargin(0.06);
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
    h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0.0, 2 );



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
  //h_Ratio->Draw("pe0,same");
  g_Ratio->Draw("pe0,same");
  
  line->DrawLine(monoBin[0],0.0,monoBin[0],2.0);

  gPad->RedrawAxis();

  c2_0->cd();
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.pdf", fullPath.c_str()) );
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.png", fullPath.c_str()) );
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.eps", fullPath.c_str()) );



  TCanvas* c2_1 = new TCanvas("c2_1", "", 1100, 600);
  c2_1->cd();
  
  TPad *pad1_1 = new TPad("pad1_1","pad1_1",0,0.3-0.1,1,1);
  pad1_1->SetBottomMargin(0.18);
  pad1_1->Draw();
  pad1_1->cd();

  pad1_1->SetLogy();
    
  oldBin=thisBin;
  thisBin=12+21;
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
  //hdata->Draw("pe0,same");
  gdata_zero->Draw("pe0,same");
  gdata->Draw("pe0,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[1]->Draw("same");

  left = pad1_1->GetLeftMargin();
  right = pad1_1->GetRightMargin();
  bot = pad1_1->GetBottomMargin();
  top = pad1_1->GetTopMargin();
  binWidth = (1.0-right-left)/(thisBin-oldBin);
  
  text->SetTextAlign(13);
  text->SetTextFont(42);
  text->SetTextAngle(0);
  text->SetTextSize(0.05);
  text->DrawLatex(left+0.04,1-top-0.01, "Pre-fit background");


  ibin = 0;
  TString vlJ[7] = {"2-3j","2-3j","2-3j","#geq4j","#geq4j","#geq4j", "#geq2j"};
  TString vlB[7] = {"0b", "1b", "2b", "0b", "1b", "2b", "#geq3b"};

  for(int nR=0; nR<7; nR++){ 
    
    xcenter = left+binWidth*(ibin+(nBins_[oldBin+nR])*0.5);
    text->SetTextAlign(23);
    text->SetTextFont(62);
    text->SetTextSize(0.030);
    
    float y=bot+(1-top-bot)*0.85;
    if (xcenter>1-right-0.19)
      y=0.67;
    
    text->DrawLatex(xcenter, y, vlJ[nR]);
    text->DrawLatex(xcenter,y-text->GetTextSize()-0.001,vlB[nR]);

    ibin+=nBins_[12+nR];

  }
  
  line = new TLine();
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  
  ibin=0;

  for(int nR=0; nR<7; nR++){

    ibin+=nBins_[oldBin+nR];
    x = left+ibin*binWidth;
    
    if(left+binWidth*(ibin+(nBins_[oldBin+nR])*0.5)>1-right-0.19)
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.65);
    else
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.85);
  
    
  }

  gPad->RedrawAxis();
  
  c2_1->cd();
  TPad *pad2_1 = new TPad("pad2_1","pad2_1",0,0,1,0.21);
  pad2_1->SetTopMargin(0.06);
  pad2_1->SetBottomMargin(0.1);
  pad2_1->Draw();
  pad2_1->cd();

  TH2D* h2_axes_ratio_1;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0.0, 2 );
  //h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0., 2.0 );

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
  //h_Ratio->Draw("pe0,same");
  g_Ratio->Draw("pe0,same");
  
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  ibin = oldBin;
  for(int nR=0; nR<7; nR++){
    ibin += nBins_[oldBin+nR];
    line->DrawLine(ibin,0.0,ibin,2.0);
  }
  
  gPad->RedrawAxis();

  c2_1->cd();
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.png", fullPath.c_str()) );
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.eps", fullPath.c_str()) );



  TCanvas* c2_2 = new TCanvas("c2_2", "", 1100, 600);
  c2_2->cd();
  
  TPad *pad1_2 = new TPad("pad1_2","pad1_2",0,0.3-0.1,1,1);
  pad1_2->SetBottomMargin(0.18);
  pad1_2->Draw();
  pad1_2->cd();

  pad1_2->SetLogy();
    
  yMin= 1e-2;

  oldBin=thisBin;
  thisBin=12+21+40;
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
  //hdata->Draw("pe0,same");
  gdata_zero->Draw("pe0,same");
  gdata->Draw("pe0,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");
    
  htBox[2]->Draw("same");

  left = pad1_2->GetLeftMargin();
  right = pad1_2->GetRightMargin();
  bot = pad1_2->GetBottomMargin();
  top = pad1_2->GetTopMargin();
  binWidth = (1.0-right-left)/(thisBin-oldBin);

  text->SetTextAlign(13);
  text->SetTextFont(42);
  text->SetTextAngle(0);
  text->SetTextSize(0.05);
  text->DrawLatex(left+0.04,1-top-0.01, "Pre-fit background");


  ibin = 0;
  TString Jlab[11] = {"2-3j","2-3j","2-3j","4-6j","4-6j","4-6j","#geq7j","#geq7j","#geq7j","2-6j", "#geq7j"};
  TString Blab[11] = {"0b", "1b", "2b", "0b", "1b", "2b", "0b", "1b", "2b","#geq3b","#geq3b"};

  for(int nR=0; nR<11; nR++){

    xcenter = left+binWidth*(ibin+(nBins_[12+7+nR])*0.5);
    text->SetTextAlign(23);
    text->SetTextFont(62);
    text->SetTextSize(0.030);

    float y=bot+(1-top-bot)*0.85;
    if (xcenter>1-right-0.19)
      y=0.67;

    text->DrawLatex(xcenter, y, Jlab[nR]);
    text->DrawLatex(xcenter,y-text->GetTextSize()-0.001,Blab[nR]);

    ibin+=nBins_[12+7+nR];

  }

  line = new TLine();
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);

  ibin=0;

  for(int nR=0; nR<11; nR++){

    ibin+=nBins_[12+7+nR];
    x = left+ibin*binWidth;

    if(left+binWidth*(ibin+(nBins_[12+7+nR])*0.5)>1-right-0.19)
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.65);
    else
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.85);


  }

  gPad->RedrawAxis();
  
  c2_2->cd();
  TPad *pad2_2 = new TPad("pad2_2","pad2_2",0,0,1,0.21);
  pad2_2->SetTopMargin(0.06);
  pad2_2->SetBottomMargin(0.1);
  pad2_2->Draw();
  pad2_2->cd();

  TH2D* h2_axes_ratio_2;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0., 3.5 );

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
  //h_Ratio->Draw("pe0,same");
  g_Ratio->Draw("pe0,same");
  
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  ibin = oldBin;
  for(int nR=0; nR<11; nR++){
    ibin += nBins_[12+7+nR];
    line->DrawLine(ibin,0,ibin,3.5);
  }

  gPad->RedrawAxis();

  c2_2->cd();
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.png", fullPath.c_str()) );
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.eps", fullPath.c_str()) );




  TCanvas* c2_3 = new TCanvas("c2_3", "", 1100, 600);
  c2_3->cd();
  
  TPad *pad1_3 = new TPad("pad1_3","pad1_3",0,0.3-0.1,1,1);
  pad1_3->SetBottomMargin(0.18);
  pad1_3->Draw();
  pad1_3->cd();

  pad1_3->SetLogy();
    
  yMin = 1e-2;
  oldBin=thisBin;
  thisBin=12+21+40+51;
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
  //hdata->Draw("pe0,same");
  gdata_zero->Draw("pe0,same");
  gdata->Draw("pe0,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[3]->Draw("same");

  left = pad1_3->GetLeftMargin();
  right = pad1_3->GetRightMargin();
  bot = pad1_3->GetBottomMargin();
  top = pad1_3->GetTopMargin();
  binWidth = (1.0-right-left)/(thisBin-oldBin);

  text->SetTextAlign(13);
  text->SetTextFont(42);
  text->SetTextAngle(0);
  text->SetTextSize(0.05);
  text->DrawLatex(left+0.04,1-top-0.01, "Pre-fit background");


  ibin = 0;
  for(int nR=0; nR<11; nR++){

    xcenter = left+binWidth*(ibin+(nBins_[12+7+11+nR])*0.5);
    text->SetTextAlign(23);
    text->SetTextFont(62);
    text->SetTextSize(0.030);

    float y=bot+(1-top-bot)*0.85;
    if (xcenter>1-right-0.19)
      y=0.67;

    text->DrawLatex(xcenter, y, Jlab[nR]);
    text->DrawLatex(xcenter,y-text->GetTextSize()-0.001,Blab[nR]);

    ibin+=nBins_[12+7+11+nR];

  }

  line = new TLine();
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);

  ibin=0;

  for(int nR=0; nR<11; nR++){

    ibin+=nBins_[12+7+11+nR];
    x = left+ibin*binWidth;

    if(left+binWidth*(ibin+(nBins_[12+7+11+nR])*0.5)>1-right-0.19)
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.65);
    else
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.85);


  }

  gPad->RedrawAxis();
  
  c2_3->cd();
  TPad *pad2_3 = new TPad("pad2_3","pad2_3",0,0,1,0.21);
  pad2_3->SetTopMargin(0.06);
  pad2_3->SetBottomMargin(0.1);
  pad2_3->Draw();
  pad2_3->cd();

  TH2D* h2_axes_ratio_3;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0., 2.0 );
  //h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0., 2.0 );

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
  //h_Ratio->Draw("pe0,same");
  g_Ratio->Draw("pe0,same");

  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  ibin = oldBin;
  for(int nR=0; nR<11; nR++){
    ibin += nBins_[12+7+11+nR];
    line->DrawLine(ibin,0.,ibin,2.0);
  }

  gPad->RedrawAxis();

  c2_3->cd();
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.png", fullPath.c_str()) );
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.eps", fullPath.c_str()) );




  TCanvas* c2_4 = new TCanvas("c2_4", "", 1100, 600);
  c2_4->cd();
  
  TPad *pad1_4 = new TPad("pad1_4","pad1_4",0,0.3-0.1,1,1);
  pad1_4->SetBottomMargin(0.18);
  pad1_4->Draw();
  pad1_4->cd();

  pad1_4->SetLogy();
    
  oldBin=thisBin;
  thisBin=12+21+40+51+53;
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
  //hdata->Draw("pe0,same");
  gdata_zero->Draw("pe0,same");
  gdata->Draw("pe0,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[4]->Draw("same");

  left = pad1_4->GetLeftMargin();
  right = pad1_4->GetRightMargin();
  bot = pad1_4->GetBottomMargin();
  top = pad1_4->GetTopMargin();
  binWidth = (1.0-right-left)/(thisBin-oldBin);

  text->SetTextAlign(13);
  text->SetTextFont(42);
  text->SetTextAngle(0);
  text->SetTextSize(0.05);
  text->DrawLatex(left+0.04,1-top-0.01, "Pre-fit background");


  ibin = 0;
  for(int nR=0; nR<11; nR++){

    xcenter = left+binWidth*(ibin+(nBins_[12+7+11*2+nR])*0.5);
    text->SetTextAlign(23);
    text->SetTextFont(62);
    text->SetTextSize(0.030);

    float y=bot+(1-top-bot)*0.85;
    if (xcenter>1-right-0.19)
      y=0.67;

    text->DrawLatex(xcenter, y, Jlab[nR]);
    text->DrawLatex(xcenter,y-text->GetTextSize()-0.001,Blab[nR]);

    ibin+=nBins_[12+7+11*2+nR];

  }

  line = new TLine();
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);

  ibin=0;

  for(int nR=0; nR<11; nR++){

    ibin+=nBins_[12+7+11*2+nR];
    x = left+ibin*binWidth;

    if(left+binWidth*(ibin+(nBins_[12+7+11*2+nR])*0.5)>1-right-0.19)
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.65);
    else
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.85);


  }

  gPad->RedrawAxis();
  
  c2_4->cd();
  TPad *pad2_4 = new TPad("pad2_4","pad2_4",0,0,1,0.21);
  pad2_4->SetTopMargin(0.06);
  pad2_4->SetBottomMargin(0.1);
  pad2_4->Draw();
  pad2_4->cd();

  TH2D* h2_axes_ratio_4;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0., 5.0 );
  //h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0., 2.0 );


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
  //h_Ratio->Draw("pe0,same");
  g_Ratio->Draw("pe0,same");
  
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  ibin = oldBin;
  for(int nR=0; nR<11; nR++){
    ibin += nBins_[12+7+11*2+nR];
    line->DrawLine(ibin,0,ibin,5.0);
  }

  gPad->RedrawAxis();

  c2_4->cd();
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.png", fullPath.c_str()) );
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.eps", fullPath.c_str()) );




  TCanvas* c2_5 = new TCanvas("c2_5", "", 1100, 600);
  c2_5->cd();
  
  TPad *pad1_5 = new TPad("pad1_5","pad1_5",0,0.3-0.1,1,1);
  pad1_5->SetBottomMargin(0.18);
  pad1_5->Draw();
  pad1_5->cd();

  pad1_5->SetLogy();
    
  yMax  /= 10;
  oldBin=thisBin;
  thisBin=213;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);

  TH1D* haxes= new TH1D("haxes", "", thisBin+1-oldBin, oldBin, thisBin+1);
  for (int b=1; b<=thisBin-oldBin; ++b)
    haxes->GetXaxis()->SetBinLabel(b, hestimate_all->GetXaxis()->GetBinLabel(b+oldBin));

  haxes->GetYaxis()->SetRangeUser(yMin, yMax);
  haxes->GetXaxis()->LabelsOption("v");
  haxes->GetXaxis()->SetLabelSize(0.042);
  haxes->GetYaxis()->SetTitle("Entries");
  haxes->GetYaxis()->SetTitleOffset(0.95);
  haxes->GetYaxis()->SetLabelSize(0.042);

  haxes->Draw("");
  hestimate_all->Draw("same");  

  //  ((TH1*)(bgStack.GetStack()->Last()))->SetBinContent(214,0);

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //hdata->Draw("pe0,same");
  gdata_zero->Draw("pe0,same");
  gdata->Draw("pe0,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[5]->Draw("same");

  left  = pad1_5->GetLeftMargin();
  right = pad1_5->GetRightMargin();
  bot   = pad1_5->GetBottomMargin();
  top   = pad1_5->GetTopMargin();
  binWidth = (1.0-right-left)/(thisBin+1-oldBin);

  text->SetTextAlign(13);
  text->SetTextFont(42);
  text->SetTextAngle(0);
  text->SetTextSize(0.05);
  text->DrawLatex(left+0.04,1-top-0.01, "Pre-fit background");


  ibin = 0;
  for(int nR=0; nR<11; nR++){

    xcenter = left+binWidth*(ibin+(nBins_[12+7+11*3+nR]-1)*0.5);
    text->SetTextAlign(23);
    text->SetTextFont(62);
    text->SetTextSize(0.030);

    float y=bot+(1-top-bot)*0.85;
    if (xcenter>1-right-0.19)
      y=0.67;

    text->DrawLatex(xcenter, y, Jlab[nR]);
    text->DrawLatex(xcenter,y-text->GetTextSize()-0.001,Blab[nR]);

    ibin+=nBins_[12+7+11*3+nR]-1;

  }

  line = new TLine();
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);

  ibin=0;

  for(int nR=0; nR<11-1; nR++){

    ibin+=(nBins_[12+7+11*3+nR]-1);
    x = left+ibin*binWidth;

    if(left+binWidth*(ibin+(nBins_[12+7+11*3+nR]-1)*0.5)>1-right-0.19)
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.65);
    else
      line->DrawLineNDC(x, bot, x, bot+(1-top-bot)*0.85);


  }

  gPad->RedrawAxis();
  
  c2_5->cd();
  TPad *pad2_5 = new TPad("pad2_5","pad2_5",0,0,1,0.21);
  pad2_5->SetTopMargin(0.06);
  pad2_5->SetBottomMargin(0.1);
  pad2_5->Draw();
  pad2_5->cd();

  TH2D* h2_axes_ratio_5;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin+1, 10, 0.1, 10.0 );
  }
  else
    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin+1, 10, 0., 3.0 );
  //h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  h2_axes_ratio_5->SetStats(0);
  h2_axes_ratio_5->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_5->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_5->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_5->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_5->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_5->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_5->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_5 = new TLine(oldBin, 1.0, thisBin+1, 1.0);
  LineCentral_5->SetLineColor(1);

  h2_axes_ratio_5->Draw("");
  h_band->Draw("E2same");
  LineCentral_5->Draw("same");
  //  h_Ratio->Draw("pe0,same");
  g_Ratio->Draw("pe0,same");
  
  line->SetNDC(1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  ibin = oldBin;
  for(int nR=0; nR<11-1; nR++){
    ibin += (nBins_[12+7+11*3+nR]-1);
    line->DrawLine(ibin,0.,ibin,3.0);
  }

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



BGTable getTable( const std::string& tableFileName ) {

  std::ifstream ifs( tableFileName.c_str() );

  BGTable table;

  while( ifs.good() ) {


    char thisLine[256];
    ifs.getline( thisLine, 256 );
    if( thisLine[0]=='#' ) continue;

    std::istringstream thisLine_iss(thisLine);

    std::string name;
    float yield, statUp, statDn, systUp, systDn;
    thisLine_iss >> name >> yield >> statUp >> statDn >> systUp >> systDn;

    if( name=="zinv" ) {
      table.zinv = yield;
      table.zinv_statUp = statUp;
      table.zinv_statDn = statDn;
      table.zinv_systUp = systUp;
      table.zinv_systDn = systDn;
    } else if( name=="llep" ) {
      table.llep = yield;
      table.llep_statUp = statUp;
      table.llep_statDn = statDn;
      table.llep_systUp = systUp;
      table.llep_systDn = systDn;
    } else if( name=="qcd" ) {
      table.qcd = yield;
      table.qcd_statUp = statUp;
      table.qcd_statDn = statDn;
      table.qcd_systUp = systUp;
      table.qcd_systDn = systDn;
    } else {
      continue;
    }

  }

  return table;

}
