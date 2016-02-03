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

  lumi = cfg.lumi();
  
  TH1::AddDirectory(kTRUE);
  
  std::string dir = cfg.getEventYieldDir();
  std::string outputdir = cfg.getEventYieldDir() + "/YieldComparison_dataMC_binned";
 
 
  MT2Analysis<MT2Estimate>* analysis = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "data" ); // any one is good, just need to know the regions                                                                    

  std::vector < MT2Analysis<MT2Estimate>* > analysesSignal;
  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1tttt_mGluino1500_mLSP100") );
  analysesSignal[0]->setName("T1tttt 1500,100");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1tttt_mGluino1200_mLSP800") );
  analysesSignal[1]->setName("T1tttt 1200,800");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1bbbb_mGluino1500_mLSP100") );
  analysesSignal[2]->setName("T1bbbb 1500,100");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1bbbb_mGluino1000_mLSP900") );
  analysesSignal[3]->setName("T1bbbb 1000,900");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1qqqq_mGluino1400_mLSP100") );
  analysesSignal[4]->setName("T1qqqq 1400,100");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1qqqq_mGluino1000_mLSP800") );
  analysesSignal[5]->setName("T1qqqq 1000,800");

  std::set<MT2Region> regions = analysis->getRegions();

  MT2Analysis<MT2Estimate>* data = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "data" );
  
  drawYields( outputdir.c_str(), data, dir );

  return 0;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::string dir ) {

  
  MT2DrawTools::setStyle();

  system(Form("mkdir -p %s", outputdir.c_str()));

  std::vector<int> colors;
  colors.push_back( 402 );
  colors.push_back( 430 );
  colors.push_back( 418 );
  
  unsigned int bgSize = 3;
  
  std::set<MT2Region> MT2Regions = data->getRegions();
  
  TH1D* hdata = new TH1D("hdata", "", 173, 0, 173);
  hdata->Sumw2();
  hdata->GetYaxis()->SetTitle("Entries");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( 1 );
  hdata->SetMarkerColor( 1 );
  
  TH1D* hestimate_all = new TH1D(Form("hestimate_all"), "", 173, 0, 173);
  hestimate_all->Sumw2();

  TH1D* hestimate[bgSize];

  TH1D* hestimate_all_forRatio = new TH1D(Form("hestimate_all_forRatio"), "", 173, 0, 173);
  hestimate_all_forRatio->Sumw2();

  TH1D* hestimate_forRatio[bgSize];
  
  for(unsigned int b=0; b<bgSize; ++b){
  
    hestimate[b]= new TH1D(Form("hestimate_%d", b), "", 173, 0, 173);
    hestimate[b]->Sumw2();
    hestimate[b]->GetYaxis()->SetTitle("Entries");
    hestimate[b]->SetFillColor(colors[b]);
    hestimate[b]->SetLineColor(1);

    hestimate_forRatio[b]= new TH1D(Form("hestimate_forRatio%d", b), "", 173, 0, 173);
    hestimate_forRatio[b]->Sumw2();
    hestimate_forRatio[b]->GetYaxis()->SetTitle("Entries");
    hestimate_forRatio[b]->SetFillColor(colors[b]);
    hestimate_forRatio[b]->SetLineColor(1);
    
  }

  THStack bgStack("bgStack", "");

  TH1D* hPull = new TH1D("hPull", "", 101, -5.05, 5.05);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Est. - Obs.)/#sigma");
  hPull->GetYaxis()->SetTitle("Entries");
  
  std::string fullPath = outputdir;

  TFile* fmono=TFile::Open("mlfit_monojet.root");
  TFile* fvlht=TFile::Open("mlfit_veryLowHT.root");
  TFile* flht =TFile::Open("mlfit_LowHT.root");
  TFile* fmht =TFile::Open("mlfit_MediumHT.root");
  TFile* fhht =TFile::Open("mlfit_HighHT.root");
  TFile* feht =TFile::Open("mlfit_ExtremeHT.root");
  
  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      int nBins;
      double *bins;
      iMT2->getBins(nBins, bins);
      
      TH1D* h_first = data->get(*iMT2)->yield;
      TGraphAsymmErrors* g_first = MT2DrawTools::getPoissonGraph(h_first);      
      
      TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), iMT2->getName().c_str()), "recreate" );
      histoFile->cd();
      //      h_first->Write();
      
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
	
	h_second[b] = new TH1D(Form("h_second_%d", b), "", nBins, bins);
	
	h_second_forRatio[b] = new TH1D(Form("h_second_%d", b), "", nBins, bins);
	
	h_second[b]->SetFillColor( colors[b] );
	h_second[b]->SetLineColor( 1 );
	
      }

      h_second_all = new TH1D(Form("h_second_all"), "", nBins, bins);
      h_second_all->Sumw2();

      h_second_forRatio_all = new TH1D(Form("h_second_forRatio_all"), "", nBins, bins);
      h_second_forRatio_all->Sumw2();
      
      for( int iBin=0; iBin<nBins; ++iBin ) {

	std::string tableName;
	if( iBin < nBins-1 )
	  tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lfto%.0lf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin], bins[iBin+1]) );
	else
	  tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lftoInf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin] ));

	BGTable thisTable = getTable(tableName);
	

	float totalPost_llep;
	float totalPost_zinv;
	float totalPost_qcd;
	float totalPost_Err_llep;
	float totalPost_Err_zinv;
	float totalPost_Err_qcd;

	float totalPost;
	float totalPost_Err;

	if(iRegion <=24){

	  int ch=iRegion+iBin;
	  fvlht->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
	  
	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;

	  totalPost_Err_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinError(1) : 0;
	  totalPost_Err_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinError(1) : 0;
	  totalPost_Err_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd"))  ? thisqcd ->GetBinError(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);

	  gDirectory->cd("..");
	      
	}
	else if(iRegion >=25 && iRegion <= 37){

	  int ch=iRegion-24+iBin;
	  
	  fmono->cd();
	  gDirectory->cd("shapes_fit_b");
	    
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
	  
	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
	  gDirectory->cd("..");

	}
	else if(iRegion >=38 && iRegion <= 68){
	  
	  int ch=iRegion-37+iBin;
	  
	  flht->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
	  
	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
	  gDirectory->cd("..");
	
	}
	else if(iRegion >=69 && iRegion <= 110){

	  int ch=iRegion-68+iBin;

	  fmht->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");

	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
	  gDirectory->cd("..");
	  
	}
	else if(iRegion >=111 && iRegion <= 145){

	  int ch=iRegion-110+iBin;
	  
	  fhht->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");

	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);	
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//  	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);

	  gDirectory->cd("..");

	}
	else if(iRegion >=146){

	  int ch=iRegion-145+iBin;

	  feht->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//  	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
	  gDirectory->cd("..");

	}
	
	gDirectory->cd();
	
	std::cout << "Filling histograms" << std::endl;
	
	//QCD
	h_second[0]->SetBinContent(iBin+1, totalPost_qcd);
	h_second[0]->SetBinError(iBin+1, 0.);
	
	h_second_forRatio[0]->SetBinContent(iBin+1, totalPost_qcd);
	h_second_forRatio[0]->SetBinError(iBin+1, 0);
	
	//Lost Lepton
	h_second[1]->SetBinContent(iBin+1, totalPost_llep);
	h_second[1]->SetBinError(iBin+1, 0.);

	h_second_forRatio[1]->SetBinContent(iBin+1, totalPost_llep);
	h_second_forRatio[1]->SetBinError(iBin+1, 0.);

	//Invisible Z
	h_second[2]->SetBinContent(iBin+1, totalPost_zinv);
	h_second[2]->SetBinError(iBin+1, 0.);

	h_second_forRatio[2]->SetBinContent(iBin+1, totalPost_zinv);
	h_second_forRatio[2]->SetBinError(iBin+1, 0);

	h_second_all->SetBinContent(iBin+1, totalPost);
	h_second_all->SetBinError(iBin+1, totalPost_Err);

	h_second_forRatio_all->SetBinContent(iBin+1, totalPost);
	h_second_forRatio_all->SetBinError(iBin+1, 0);

      }	
	
      for(unsigned int b=0; b<bgSize; ++b){
      
	bgStack_region.Add(h_second[b]);
	      
      }
      
      for(int iBin=1; iBin<=nBins; ++iBin){

	float thisData    = h_first->GetBinContent(iBin);
	float thisDataErr = h_first->GetBinError(iBin);

	float thisEst     = h_second_all->GetBinContent(iBin);
	float thisEstErr  = h_second_all->GetBinError(iBin);
	

	hPull->Fill( (thisEst-thisData)/( TMath::Sqrt( thisDataErr*thisDataErr + thisEstErr*thisEstErr ) ) );
	
      }

      
      double err_data;
      double int_data;
      for (int iBin=1; iBin<=nBins; ++iBin){
	
	int_data = h_first->GetBinContent(iBin);
	err_data = h_first->GetBinError(iBin);
	
	hdata->SetBinContent(iRegion, int_data);
	hdata->SetBinError(iRegion, err_data);

	hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
	
	for(unsigned int b=0; b<bgSize; ++b){

	  double err_int = fabs(h_second[b]->GetBinError(iBin));
	  double integral = fabs(h_second[b]->GetBinContent(iBin));
	  hestimate[b]->SetBinContent(iRegion, integral);
	  hestimate[b]->SetBinError(iRegion, err_int);

	  double integral_forRatio = fabs(h_second_forRatio[b]->GetBinContent(iBin));
	  hestimate_forRatio[b]->SetBinContent(iRegion, integral_forRatio);
	  hestimate_forRatio[b]->SetBinError(iRegion, 0);
	  
	  std::string thisLabel=Form("%s,%d", niceNames[1].c_str(), iBin); 
	  hestimate[b]->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
	  
	}

	hestimate_all->SetBinContent(iRegion, h_second_all->GetBinContent(iBin));
	hestimate_all->SetBinError(iRegion, h_second_all->GetBinError(iBin));
	std::string thisLabel=Form("%s,%d", niceNames[1].c_str(), iBin); 
	hestimate_all->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
		
	hestimate_all_forRatio->SetBinContent(iRegion, h_second_all->GetBinContent(iBin));
	hestimate_all_forRatio->SetBinError(iRegion, 0);
	hestimate_all_forRatio->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
		
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


      TLegend* legend = new TLegend( 0.6, 0.9-(bgSize)*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
      legend->AddEntry( g_first, "Data", "P" );

      legend->Draw("same");

      bgStack_region.Draw("histo, same");
      g_first->Draw("pe,same");
      
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
      //      h_ratio->Write();
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

      for(unsigned int b=0; b<bgSize; ++b)
	delete h_second[b];
      
      //      ++iRegion;

  } // for MT2 regions


  for(unsigned int b=0; b<bgSize; ++b){

    hestimate[b]->SetLineWidth(0);
    bgStack.Add(hestimate[b]);
    //bgStack.Add(hestimate_forRatio[b]);
    
//    if(b==0) hestimate_all = (TH1D*) hestimate[b]->Clone("hestimate_all");
//    else hestimate_all->Add(hestimate[b]);
//
//    if(b==0) hestimate_all_forRatio = (TH1D*) hestimate_forRatio[b]->Clone("hestimate_all_forRatio");
//    else hestimate_all_forRatio->Add(hestimate_forRatio[b]);

  }

  
  TCanvas* c2 = new TCanvas("c2", "", 1200, 600);
  c2->cd();
  

  std::string thisName = Form("%s_ratio", hdata->GetName());
  TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
  h_Ratio->Divide(hestimate_all_forRatio);
  //  h_Ratio->Write();
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
  
  int thisBin=173;
  
  hestimate_all->GetXaxis()->SetRangeUser(0, thisBin);
  hdata->GetXaxis()->SetRangeUser(0, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->SetFillStyle(3244);
  hestimate_all->SetFillColor(kGray+2);

  TGraphAsymmErrors* g_data = new TGraphAsymmErrors(0);
  for( int iBin=1; iBin<(hdata->GetXaxis()->GetNbins()+1); ++iBin ) {

    double y;
    double x, xerr;

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
  
  TLegend* legend = new TLegend( 0.8, 0.9-(bgSize+1-1)*0.06-0.06, 0.93, 0.9-0.06 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( hdata, "Data", "PL" );
  legend->AddEntry( hestimate[0], "Multijet", "F");
  legend->AddEntry( hestimate[1], "Lost Lepton", "F");
  legend->AddEntry( hestimate[2], "Invisible Z", "F");

  legend->Draw("same");

  //  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
  labelTop->Draw("same");
  
  int nHTRegions = 6;
  std::vector< std::string > htRegions;
  htRegions.push_back("H_{T} [200, 450] GeV");
  htRegions.push_back("1 Jet");
  htRegions.push_back("H_{T} [450, 575] GeV");
  htRegions.push_back("H_{T} [575, 1000] GeV");
  htRegions.push_back("H_{T} [1000, 1500] GeV");
  htRegions.push_back("H_{T} > 1500 GeV");
//  htRegions.push_back("very low H_{T}");
//  htRegions.push_back("1 Jet");
//  htRegions.push_back("low H_{T}");
//  htRegions.push_back("medium H_{T}");
//  htRegions.push_back("high H_{T}");
//  htRegions.push_back("extreme H_{T}");
  
  TPaveText* htBox[5];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){
    
    htBox[iHT] = new TPaveText(0.4, 0.9-0.06, 0.7, 0.85, "brNDC");
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

  TH2D* h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0., 3.0 );
  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.3);
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
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.pdf", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.png", fullPath.c_str()) );



  TCanvas* c2_0 = new TCanvas("c2_0", "", 1200, 600);
  c2_0->cd();
  
  TPad *pad1_0 = new TPad("pad1_0","pad1_0",0,0.3-0.1,1,1);
  pad1_0->SetBottomMargin(0.15);
  pad1_0->Draw();
  pad1_0->cd();


  int oldBin=0;
  pad1_0->SetLogy();
    
  thisBin=24;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  hdata->Draw("pe,same");

  legend->Draw("same");

  labelTop->Draw("same");
  
  htBox[0]->Draw("same");

  gPad->RedrawAxis();
  
  c2_0->cd();
  TPad *pad2_0 = new TPad("pad2_0","pad2_0",0,0,1,0.21);
  pad2_0->SetTopMargin(0.05);
  pad2_0->SetBottomMargin(0.1);
  pad2_0->Draw();
  pad2_0->cd();

  TH2D* h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_0->SetStats(0);
  h2_axes_ratio_0->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_0->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_0->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_0->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio_0->GetYaxis()->SetTitleOffset(0.3);
  h2_axes_ratio_0->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_0->GetYaxis()->SetTitle("Ratio");
  
  TLine* LineCentral_0 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_0->SetLineColor(1);

  h2_axes_ratio_0->Draw("");
  h_band->Draw("E2same");
  LineCentral_0->Draw("same");
  h_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_0->cd();
  c2_0->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_0->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.png", fullPath.c_str()) );



  TCanvas* c2_1 = new TCanvas("c2_1", "", 1200, 600);
  c2_1->cd();
  
  TPad *pad1_1 = new TPad("pad1_1","pad1_1",0,0.3-0.1,1,1);
  pad1_1->SetBottomMargin(0.15);
  pad1_1->Draw();
  pad1_1->cd();

  pad1_1->SetLogy();
    
  oldBin=thisBin;
  thisBin=37;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  hdata->Draw("pe,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  
  htBox[1]->Draw("same");

  gPad->RedrawAxis();
  
  c2_1->cd();
  TPad *pad2_1 = new TPad("pad2_1","pad2_1",0,0,1,0.21);
  pad2_1->SetTopMargin(0.05);
  pad2_1->SetBottomMargin(0.1);
  pad2_1->Draw();
  pad2_1->cd();

  TH2D* h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_1->SetStats(0);
  h2_axes_ratio_1->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_1->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_1->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_1->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio_1->GetYaxis()->SetTitleOffset(0.3);
  h2_axes_ratio_1->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_1->GetYaxis()->SetTitle("Ratio");
  
  TLine* LineCentral_1 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_1->SetLineColor(1);

  h2_axes_ratio_1->Draw("");
  h_band->Draw("E2same");
  LineCentral_1->Draw("same");
  h_Ratio->Draw("pe,same");

  gPad->RedrawAxis();

  c2_1->cd();
  c2_1->SaveAs( Form("%s/mt2_monojet_fullEstimate.pdf", fullPath.c_str()) );
  c2_1->SaveAs( Form("%s/mt2_monojet_fullEstimate.png", fullPath.c_str()) );



  TCanvas* c2_2 = new TCanvas("c2_2", "", 1200, 600);
  c2_2->cd();
  
  TPad *pad1_2 = new TPad("pad1_2","pad1_2",0,0.3-0.1,1,1);
  pad1_2->SetBottomMargin(0.15);
  pad1_2->Draw();
  pad1_2->cd();

  pad1_2->SetLogy();
    
  oldBin=thisBin;
  thisBin=68;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  hdata->Draw("pe,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  
  htBox[2]->Draw("same");

  gPad->RedrawAxis();
  
  c2_2->cd();
  TPad *pad2_2 = new TPad("pad2_2","pad2_2",0,0,1,0.21);
  pad2_2->SetTopMargin(0.05);
  pad2_2->SetBottomMargin(0.1);
  pad2_2->Draw();
  pad2_2->cd();

  TH2D* h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_2->SetStats(0);
  h2_axes_ratio_2->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_2->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_2->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_2->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio_2->GetYaxis()->SetTitleOffset(0.3);
  h2_axes_ratio_2->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_2->GetYaxis()->SetTitle("Ratio");
  
  TLine* LineCentral_2 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_2->SetLineColor(1);

  h2_axes_ratio_2->Draw("");
  h_band->Draw("E2same");
  LineCentral_2->Draw("same");
  h_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_2->cd();
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.png", fullPath.c_str()) );




  TCanvas* c2_3 = new TCanvas("c2_3", "", 1200, 600);
  c2_3->cd();
  
  TPad *pad1_3 = new TPad("pad1_3","pad1_3",0,0.3-0.1,1,1);
  pad1_3->SetBottomMargin(0.15);
  pad1_3->Draw();
  pad1_3->cd();

  pad1_3->SetLogy();
    
  oldBin=thisBin;
  thisBin=110;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  hdata->Draw("pe,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  
  htBox[3]->Draw("same");

  gPad->RedrawAxis();
  
  c2_3->cd();
  TPad *pad2_3 = new TPad("pad2_3","pad2_3",0,0,1,0.21);
  pad2_3->SetTopMargin(0.05);
  pad2_3->SetBottomMargin(0.1);
  pad2_3->Draw();
  pad2_3->cd();

  TH2D* h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_3->SetStats(0);
  h2_axes_ratio_3->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_3->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_3->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_3->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio_3->GetYaxis()->SetTitleOffset(0.3);
  h2_axes_ratio_3->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_3->GetYaxis()->SetTitle("Ratio");
  
  TLine* LineCentral_3 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_3->SetLineColor(1);

  h2_axes_ratio_3->Draw("");
  h_band->Draw("E2same");
  LineCentral_3->Draw("same");
  h_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_3->cd();
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.png", fullPath.c_str()) );




  TCanvas* c2_4 = new TCanvas("c2_4", "", 1200, 600);
  c2_4->cd();
  
  TPad *pad1_4 = new TPad("pad1_4","pad1_4",0,0.3-0.1,1,1);
  pad1_4->SetBottomMargin(0.15);
  pad1_4->Draw();
  pad1_4->cd();

  pad1_4->SetLogy();
    
  oldBin=thisBin;
  thisBin=145;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  hdata->Draw("pe,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  
  htBox[4]->Draw("same");

  gPad->RedrawAxis();
  
  c2_4->cd();
  TPad *pad2_4 = new TPad("pad2_4","pad2_4",0,0,1,0.21);
  pad2_4->SetTopMargin(0.05);
  pad2_4->SetBottomMargin(0.1);
  pad2_4->Draw();
  pad2_4->cd();

  TH2D* h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_4->SetStats(0);
  h2_axes_ratio_4->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_4->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_4->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_4->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio_4->GetYaxis()->SetTitleOffset(0.3);
  h2_axes_ratio_4->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_4->GetYaxis()->SetTitle("Ratio");
  
  TLine* LineCentral_4 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_4->SetLineColor(1);

  h2_axes_ratio_4->Draw("");
  h_band->Draw("E2same");
  LineCentral_4->Draw("same");
  h_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_4->cd();
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.png", fullPath.c_str()) );




  TCanvas* c2_5 = new TCanvas("c2_5", "", 1200, 600);
  c2_5->cd();
  
  TPad *pad1_5 = new TPad("pad1_5","pad1_5",0,0.3-0.1,1,1);
  pad1_5->SetBottomMargin(0.15);
  pad1_5->Draw();
  pad1_5->cd();

  pad1_5->SetLogy();
    
  oldBin=thisBin;
  thisBin=173;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  hdata->Draw("pe,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  
  htBox[5]->Draw("same");

  gPad->RedrawAxis();
  
  c2_5->cd();
  TPad *pad2_5 = new TPad("pad2_5","pad2_5",0,0,1,0.21);
  pad2_5->SetTopMargin(0.05);
  pad2_5->SetBottomMargin(0.1);
  pad2_5->Draw();
  pad2_5->cd();

  TH2D* h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_5->SetStats(0);
  h2_axes_ratio_5->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_5->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_5->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_5->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio_5->GetYaxis()->SetTitleOffset(0.3);
  h2_axes_ratio_5->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_5->GetYaxis()->SetTitle("Ratio");
  
  TLine* LineCentral_5 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_5->SetLineColor(1);

  h2_axes_ratio_5->Draw("");
  h_band->Draw("E2same");
  LineCentral_5->Draw("same");
  h_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_5->cd();
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.png", fullPath.c_str()) );
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.C", fullPath.c_str()) );


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
