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
  std::string outputdir = cfg.getEventYieldDir() + "/QCDForMonoJetValidation_dataMC";
 
 
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

  MT2Analysis<MT2Estimate>* data  = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "data" );
  MT2Analysis<MT2Estimate>* qcdMC = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "QCD" );
  
  drawYields( outputdir.c_str(), qcdMC, dir );

  return 0;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::string dir ) {

  
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
  
  std::set<MT2Region> MT2Regions = data->getRegions();
  
  TH1D* hdata = new TH1D("hdata", "", 63, 0, 63);
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

  TH1D* hmc_qcd = new TH1D("hmc_qcd", "", 63, 0, 63);
  hmc_qcd->Sumw2();

  for(unsigned int b=0; b<bgSize; ++b){
  
    hestimate[b]= new TH1D(Form("hestimate_%d", b), "", 63, 0, 63);
    hestimate[b]->Sumw2();
    hestimate[b]->GetYaxis()->SetTitle("Entries");
    hestimate[b]->SetFillColor(colors[b]);
    hestimate[b]->SetLineColor(1);

    hestimate_forRatio[b]= new TH1D(Form("hestimate_forRatio%d", b), "", 63, 0, 63);
    hestimate_forRatio[b]->Sumw2();
    hestimate_forRatio[b]->GetYaxis()->SetTitle("Entries");
    hestimate_forRatio[b]->SetFillColor(colors[b]);
    hestimate_forRatio[b]->SetLineColor(1);
    
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
  
  std::string labelsMono[12]={"[250,350]","[350,450]","[450,575]","[575,700]","[700,1000]","[1000,1200]", ">1200","[250,350]","[350,450]","[450,575]","[575,700]", ">700"};

  TFile* bigHistoFile = TFile::Open( Form("%s/histograms_ALL.root", fullPath.c_str()), "recreate" );

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

      TH1D* h_qcdMC = data->get(*iMT2)->yield;
      h_qcdMC->Scale(lumi);
      //      h_qcdMC->Scale(18.12);
      
      for(unsigned int b=0; b< bgSize; ++b){
	
	h_second[b] = new TH1D(Form("h_second_%d", b), "", nBins, bins);
	
	h_second_forRatio[b] = new TH1D(Form("h_second_%d", b), "", nBins, bins);
	
	h_second[b]->SetFillColor( colors[b] );
	h_second[b]->SetLineColor( 1 );
	
      }
      
      for( int iBin=0; iBin<nBins; ++iBin ) {

	std::string tableName;
	if(iMT2->nJetsMax()==1){
	  tableName = std::string(Form("%s/datacard_templates/table_%s_m0toInf.txt", dir.c_str(), iMT2->getName().c_str() ));
	}
	else{
	  if( iBin < nBins-1 )
	    tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lfto%.0lf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin], bins[iBin+1]) );
	  else
	    tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lftoInf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin] ));
	}

	std::cout << tableName << std::endl;

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
	
      for(unsigned int b=1; b<bgSize; ++b){ // All but QCD
      
	bgStack_region.Add(h_second[b]);
	
	if(b==1) h_second_all = (TH1D*) h_second[b]->Clone("h_second_all");
	else h_second_all->Add(h_second[b]);
	
	if(b==1) h_second_forRatio_all = (TH1D*) h_second_forRatio[b]->Clone("h_second_forRatio_all");
	else h_second_forRatio_all->Add(h_second_forRatio[b]);
      
      }
      
      h_second_all->Write();

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

      int firstBin=1;
      double err_data;
      double int_data;
     

//      if( iMT2->nJetsMax()==1 ){
//        iRegion=iRegion-24;
//      
//	for (int iBin=1; iBin<=nBins; ++iBin){
//	  
//	  int_data = h_first->GetBinContent(iBin);
//	  err_data = h_first->GetBinError(iBin);
//	  
//	  hdata->SetBinContent(iRegion, int_data);
//	  hdata->SetBinError(iRegion, err_data);
//	  
//	  hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
//	  
//	  for(unsigned int b=0; b<bgSize; ++b){
//	    
//	    double err_int = fabs(h_second[b]->GetBinError(iBin));
//	    double integral = fabs(h_second[b]->GetBinContent(iBin));
//	    hestimate[b]->SetBinContent(iRegion, integral);
//	    hestimate[b]->SetBinError(iRegion, err_int);
//	    
//	    double integral_forRatio = fabs(h_second_forRatio[b]->GetBinContent(iBin));
//	    hestimate_forRatio[b]->SetBinContent(iRegion, integral_forRatio);
//	    hestimate_forRatio[b]->SetBinError(iRegion, 0);
//	    
//	    std::string thisLabel=Form("%s,%d", niceNames[1].c_str(), iBin); 
//	    hestimate[b]->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
//	    
//	  }
//	  
//	  std::cout <<"iRegion " << iRegion << std::endl;
//	  if(iBin<nBins)
//	    ++iRegion;
//	  
//	  
//	}
//	if(iMT2->nBJetsMax()==0)
//	  iRegion=iRegion+24;
//	else
//	  iRegion=iRegion+11;
//      }
//      else{

        int_data = h_first->IntegralAndError(firstBin, nBins+1, err_data);

        hdata->SetBinContent(iRegion, int_data);
        hdata->SetBinError(iRegion, err_data);

	std::cout<<"Filled data  with: " << int_data << std::endl;

	if( iMT2->nJetsMax()==1 )
	  hdata->GetXaxis()->SetBinLabel( iRegion, labelsMono[iRegion-1].c_str() );
	else
	  hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
	
        for(unsigned int b=0; b<bgSize; ++b){

          double err_int;
          double integral = h_second[b]->IntegralAndError(firstBin, nBins+1, err_int);
          hestimate[b]->SetBinContent(iRegion, integral);
          hestimate[b]->SetBinError(iRegion, err_int);

	  double err_int_forRatio;
	  double integral_forRatio = h_second[b]->IntegralAndError(firstBin, nBins+1, err_int_forRatio);;
	  hestimate_forRatio[b]->SetBinContent(iRegion, integral_forRatio);
	  hestimate_forRatio[b]->SetBinError(iRegion, 0);
	  
	  std::cout<<"Filled estimate for background " << b << " with: " << integral << std::endl;

	  std::string thisLabel=Form("%s", niceNames[1].c_str());
	  if( iMT2->nJetsMax()==1 )
	    hestimate[b]->GetXaxis()->SetBinLabel( iRegion, labelsMono[iRegion-1].c_str() );
	  else
	    hestimate[b]->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );

        }
	
	double err_mcQCD;
	double integral_mcQCD = h_qcdMC->IntegralAndError(firstBin, nBins+1, err_mcQCD);
	hmc_qcd->SetBinContent(iRegion, integral_mcQCD);
	hmc_qcd->SetBinError(iRegion, err_mcQCD);
	std::string thisLabel=Form("%s", niceNames[1].c_str());
	if( iMT2->nJetsMax()==1 )
	  hmc_qcd->GetXaxis()->SetBinLabel( iRegion, labelsMono[iRegion-1].c_str() );
	else
	  hmc_qcd->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );

	//      }

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
      
      //      lumi = 18.12;

      TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
      labelTop->Draw("same");

      TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Preliminary");
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

      for(unsigned int b=0; b<bgSize; ++b)
	delete h_second_forRatio[b];
      
      ++iRegion;

  } // for MT2 regions


  for(unsigned int b=0; b<bgSize; ++b){

    if(b>0){
      hestimate[b]->SetLineWidth(0);
      bgStack.Add(hestimate[b]);
    }
    //bgStack.Add(hestimate_forRatio[b]);
    
    if(b==0) continue;

    if(b==1) hestimate_all = (TH1D*) hestimate[b]->Clone("hestimate_all");
    else hestimate_all->Add(hestimate[b]);

    if(b==1) hestimate_all_forRatio = (TH1D*) hestimate_forRatio[b]->Clone("hestimate_all_forRatio");
    else hestimate_all_forRatio->Add(hestimate_forRatio[b]);

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

  std::string thisName = Form("%s_ratio", hestimate[0]->GetName());
  TH1D* h_Ratio = (TH1D*) hestimate[0]->Clone(thisName.c_str());
  h_Ratio->Divide(hestimate_all_forRatio);
  h_Ratio->Write();
  h_Ratio->SetStats(0);
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->GetXaxis()->SetLabelSize(0.042);
  h_Ratio->GetXaxis()->SetLabelFont(62);
  h_Ratio->GetYaxis()->SetNdivisions(5,5,0);
  h_Ratio->GetYaxis()->SetRangeUser(0.0,1.0);
  h_Ratio->GetYaxis()->SetTitle("Multijet / Non-Multijet");

  std::string thisNameMC = Form("%s_ratioMC", hmc_qcd->GetName());
  TH1D* h_RatioMC = (TH1D*) hmc_qcd->Clone(thisName.c_str());
  h_RatioMC->Divide(hestimate_all_forRatio);
  h_RatioMC->Write();
  h_RatioMC->SetStats(0);
  h_RatioMC->SetMarkerStyle(4);
  h_RatioMC->SetMarkerColor(kYellow+1);
  h_RatioMC->SetLineColor(kYellow+1);
  h_RatioMC->GetXaxis()->SetLabelSize(0.042);
  h_RatioMC->GetXaxis()->SetLabelFont(62);
  h_RatioMC->GetYaxis()->SetNdivisions(5,5,0);
  h_RatioMC->GetYaxis()->SetRangeUser(0.0,1.0);
  h_RatioMC->GetYaxis()->SetTitle("Multijet / Non-Multijet");

  int thisBin=12;
  h_Ratio->GetXaxis()->SetRangeUser(0,thisBin);
  h_RatioMC->GetXaxis()->SetRangeUser(0,thisBin);

  TLegend* legend = new TLegend( 0.2, 0.7, 0.4, 0.9 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( h_Ratio, "Multijet (data-driven)", "PL" );
  legend->AddEntry( h_RatioMC, "Multijet (simulation)", "PL" );
//  legend->AddEntry( hestimate[0], "Multijet", "F");
//  legend->AddEntry( hestimate[1], "Lost lepton", "F");
//  legend->AddEntry( hestimate[2], "Z #rightarrow #nu#bar{#nu}", "F");


  bool doLogRatio = false;
  
  TH2D* h2_axes_ratio;

  if(doLogRatio){
    
    gPad->SetLogy();
    h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0.1, 10.0 );
    h_Ratio->GetYaxis()->SetRangeUser(0.1, 10.0);
    g_Ratio->GetYaxis()->SetRangeUser(0.1, 10.0);
    
  }
  else
    h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0., 0.15 );
  // h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0., 3.5 );
  
  h_Ratio  ->GetYaxis()->SetRangeUser(0.000,0.2);
  h_RatioMC->GetYaxis()->SetRangeUser(0.000,0.2);

  //  TH2D* h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0., 3.0 );
  h2_axes_ratio->SetStats(0);
  //  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  //  h2_axes_ratio->GetYaxis()->SetTitleSize(0.18);
  //  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.4);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.042);
  h2_axes_ratio->GetXaxis()->SetLabelFont(62);
  //  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.26);
  //  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitle("Multijet/Non-Multijet");


  
  for (int i=1; i<=12; ++i){

    std::cout << h_Ratio->GetBinContent(i) << " +- " << h_Ratio->GetBinError(i) << std::endl;
    std::cout << h_RatioMC->GetBinContent(i) << " +- " << h_RatioMC->GetBinError(i) << std::endl;
    
  }

//  TLine* LineCentral = new TLine(0, 1.0, thisBin, 1.0);
//  LineCentral->SetLineColor(1);

  
//  std::string thisName_Band =  Form("%s_band", hestimate_all->GetName());
//  TH1D* h_band = (TH1D*)hestimate_all->Clone(thisName_Band.c_str());
//  h_band->SetMarkerSize(0);
//  h_band->SetFillColor (kGray+2);
//  h_band->SetFillStyle (3244);
//  for ( int iBin=1; iBin <= hestimate_all->GetNbinsX(); iBin++){
//    
//    h_band->SetBinContent(iBin,1);
//
//    double error=0;
//
//    if(hestimate_all->GetBinContent(iBin)>0)
//      error = hestimate_all->GetBinError(iBin)/hestimate_all->GetBinContent(iBin);
//    else error = hestimate_all->GetBinError(iBin);
//
//    h_band->SetBinError(iBin, error);
//
//  }


//  h2_axes_ratio->Draw("");
  //  h_band->Draw("E2same");
  //  LineCentral->Draw("same");
  //  h_Ratio->Draw("pe,same");
  
  //  gPad->SetLogy();

  h_Ratio->Draw("pe");
  h_RatioMC->Draw("pe,same");
  
  legend->Draw("same");

  TLine* l_b = new TLine(7, 0, 7, 0.2);
  l_b->SetLineColor(kBlack);
  l_b->SetLineStyle(3);
  l_b->SetLineWidth(2);

  l_b->Draw("same");


  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
  labelTop->Draw("same");

  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Preliminary");
  labelCMS->Draw("same");


  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/QCDForMonojet_validation.pdf", fullPath.c_str()) );
  c2->SaveAs( Form("%s/QCDForMonojet_validation.C", fullPath.c_str()) );
  c2->SaveAs( Form("%s/QCDForMonojet_validation.png", fullPath.c_str()) );
  c2->SaveAs( Form("%s/QCDForMonojet_validation.eps", fullPath.c_str()) );
  
//  bigHistoFile->cd();
//  h_band->Write();
//  hestimate_all->Write();
//  hestimate_all_forRatio->Write();
//  hdata->Write();
//  for(unsigned int b=0; b<bgSize; ++b)
//    hestimate[b]->Write();
//
//  //  gStyle->SetOptStat(0110);
//  gStyle->SetOptFit(0011);
//  
//  TF1* fgauss= new TF1("fgauss", "gaus", -5, 5);
//  fgauss->SetLineColor(2);
//
//  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
//  c3->cd();
//  
//  hPull->SetStats(1100);
//  //  hPull->GetYaxis()->SetRangeUser(0, 15);
//  hPull->Draw("hist");
//  hPull->Fit("fgauss");
//  fgauss->Draw("same");
//  
//  labelTop = MT2DrawTools::getLabelTop(lumi);
//  labelTop->Draw("same");
//
//  labelCMS = MT2DrawTools::getLabelCMS("CMS Supplementary");
//  labelCMS->Draw("same");
//  
//  TPaveText *arxiv = new TPaveText(0.2, 0.9-0.05, 0.35, 0.9, "brNDC");
//  arxiv->AddText( "arXiv:1603.04053" );
//  arxiv->SetBorderSize(0);
//  arxiv->SetFillColor(kWhite);
//  //  arxiv->SetTextSize(0.035);
//  arxiv->SetTextAlign(11); // align centered                                                                                                                                  
//  arxiv->SetTextFont(42);
//  arxiv->Draw("same");
//  
//  c3->SaveAs( Form("%s/PullDistribution.pdf", fullPath.c_str()) );
//  c3->SaveAs( Form("%s/PullDistribution.png", fullPath.c_str()) );
//  c3->SaveAs( Form("%s/PullDistribution.root", fullPath.c_str()) );
//
//  TCanvas* c4 = new TCanvas("c4", "", 600, 600);
//  c4->cd();
//  hPvalue->SetStats(1110);
//  hPvalue->Draw("hist");
//  c4->SaveAs( Form("%s/PvalueDistribution.pdf", fullPath.c_str()) );
//  c4->SaveAs( Form("%s/PvalueDistribution.png", fullPath.c_str()) );
  
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
