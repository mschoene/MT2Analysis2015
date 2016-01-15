#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"

#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2DrawTools.h"





TH1D* drawQCDfrac( const std::string& htName, const std::vector<double> bins, MT2Analysis<MT2Estimate>* qcd, MT2Analysis<MT2Estimate>* zinv, MT2Analysis<MT2Estimate>* llep );





int main() {

  MT2DrawTools::setStyle();

  std::string path = "/scratch/mmasciov/analysisCode_forMerge/analysis/EventYields_data_Run2015D_25nsGolden_miniAODv2_fullStat_original/";

  MT2Analysis<MT2Estimate> *qcd  = MT2Analysis<MT2Estimate>::readFromFile( path + "/qcdEstimateData.root", "qcdEstimate" );
  MT2Analysis<MT2Estimate> *zinv = MT2Analysis<MT2Estimate>::readFromFile( path + "/zinvFromGamma.root", "ZinvEstimate" );
  MT2Analysis<MT2Estimate> *llep = MT2Analysis<MT2Estimate>::readFromFile( path + "/llepEstimate.root", "llepEstimate" );


  std::vector<double> bins_vlowHT;
  bins_vlowHT.push_back( 200. );
  bins_vlowHT.push_back( 300. );
  bins_vlowHT.push_back( 400. );
  bins_vlowHT.push_back( 2500. );

  std::vector<double> bins_lowHT;
  bins_lowHT.push_back( 200. );
  bins_lowHT.push_back( 300. );
  bins_lowHT.push_back( 400. );
  bins_lowHT.push_back( 500. );
  bins_lowHT.push_back( 2500. );

  std::vector<double> bins_medHT;
  bins_medHT.push_back( 200. );
  bins_medHT.push_back( 300. );
  bins_medHT.push_back( 400. );
  bins_medHT.push_back( 600. );
  bins_medHT.push_back( 800. );
  bins_medHT.push_back( 2500. );

  std::vector<double> bins_highHT;
  bins_highHT.push_back( 200. );
  bins_highHT.push_back( 400. );
  bins_highHT.push_back( 600. );
  bins_highHT.push_back( 800. );
  bins_highHT.push_back( 1000. );
  bins_highHT.push_back( 2500. );

  std::vector<double> bins_extHT;
  bins_extHT.push_back( 200. );
  bins_extHT.push_back( 400. );
  bins_extHT.push_back( 600. );
  bins_extHT.push_back( 800. );
  bins_extHT.push_back( 1000. );
  bins_extHT.push_back( 2500. );


  
  TH1D* qcdFrac_HT200to450   = drawQCDfrac( "HT200to450"  , bins_vlowHT, qcd, zinv, llep );
  TH1D* qcdFrac_HT450to575   = drawQCDfrac( "HT450to575"  , bins_lowHT , qcd, zinv, llep );
  TH1D* qcdFrac_HT575to1000  = drawQCDfrac( "HT575to1000" , bins_medHT , qcd, zinv, llep );
  TH1D* qcdFrac_HT1000to1500 = drawQCDfrac( "HT1000to1500", bins_highHT, qcd, zinv, llep );
  TH1D* qcdFrac_HT1500toInf  = drawQCDfrac( "HT1500toInf" , bins_extHT , qcd, zinv, llep );


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  c1->SetLogy();

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 2500., 10, 0.0005, 1. );
  h2_axes->SetYTitle( "QCD / ( Total BG )" );
  h2_axes->SetXTitle( "M_{T2} [GeV]" );
  h2_axes->Draw();

  qcdFrac_HT200to450   ->SetMarkerStyle(20);
  qcdFrac_HT450to575   ->SetMarkerStyle(24);
  qcdFrac_HT575to1000  ->SetMarkerStyle(21);
  qcdFrac_HT1000to1500 ->SetMarkerStyle(25);
  qcdFrac_HT1500toInf  ->SetMarkerStyle(26);

  qcdFrac_HT200to450   ->SetMarkerColor(kBlack);
  qcdFrac_HT450to575   ->SetMarkerColor(46);
  qcdFrac_HT575to1000  ->SetMarkerColor(38);
  qcdFrac_HT1000to1500 ->SetMarkerColor(29);
  qcdFrac_HT1500toInf  ->SetMarkerColor(kBlack);

  qcdFrac_HT200to450   ->SetLineColor(kBlack);
  qcdFrac_HT450to575   ->SetLineColor(46);
  qcdFrac_HT575to1000  ->SetLineColor(38);
  qcdFrac_HT1000to1500 ->SetLineColor(29);
  qcdFrac_HT1500toInf  ->SetLineColor(kBlack);

  qcdFrac_HT200to450   ->SetMarkerSize(1.5);
  qcdFrac_HT450to575   ->SetMarkerSize(1.5);
  qcdFrac_HT575to1000  ->SetMarkerSize(1.5);
  qcdFrac_HT1000to1500 ->SetMarkerSize(1.5);
  qcdFrac_HT1500toInf  ->SetMarkerSize(1.5);

  qcdFrac_HT200to450   ->Draw("psame");
  qcdFrac_HT450to575   ->Draw("psame");
  qcdFrac_HT575to1000  ->Draw("psame");
  qcdFrac_HT1000to1500 ->Draw("psame");
  qcdFrac_HT1500toInf  ->Draw("psame");

  TLegend* legend = new TLegend( 0.5, 0.2, 0.9, 0.5 );
  legend->SetTextSize( 0.033 );
  legend->SetFillColor( 0 );
  legend->AddEntry( qcdFrac_HT200to450,   "200 < H_{T} < 450 GeV", "P" );
  legend->AddEntry( qcdFrac_HT450to575,   "450 < H_{T} < 575 GeV", "P" );
  legend->AddEntry( qcdFrac_HT575to1000,  "575 < H_{T} < 1000 GeV", "P" );
  legend->AddEntry( qcdFrac_HT1000to1500, "1000 < H_{T} < 1500 GeV", "P" );
  legend->AddEntry( qcdFrac_HT1500toInf,  "H_{T} > 1500 GeV", "P" );
  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop( 2.1 );
  labelTop->Draw("same");

  c1->SaveAs("qcdFrac.eps"); 
  c1->SaveAs("qcdFrac.pdf"); 


  return 0;

}



TH1D* drawQCDfrac( const std::string& htName, const std::vector<double> bins, MT2Analysis<MT2Estimate>* qcd, MT2Analysis<MT2Estimate>* zinv, MT2Analysis<MT2Estimate>* llep ) {

  Double_t bins_d[bins.size()];
  for( int i=0;i<bins.size(); ++i ) {
    bins_d[i] = bins[i];
  }

  TH1D* h1_qcd = new TH1D( Form("qcd_%s", htName.c_str()), "", bins.size()-1, bins_d );
  h1_qcd->Sumw2();
  TH1D* h1_all = new TH1D( Form("all_%s", htName.c_str()), "", bins.size()-1, bins_d );
  h1_all->Sumw2();

  std::set<MT2Region> regions = qcd->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    TString regionName(iR->getName());
    if( !(regionName.Contains(htName)) ) continue;

    MT2Estimate* thisQCD  = qcd ->get( *iR );
    MT2Estimate* thisZinv = zinv->get( *iR );
    MT2Estimate* thisLlep = llep->get( *iR );

    if( thisQCD==0 || thisZinv==0 || thisLlep==0 ) {
      std::cout << "THERE MUST BE A PROBLEM" << std::endl;
      std::cout << "qcd: "  << thisQCD << std::endl;
      std::cout << "llep: " << thisLlep << std::endl;
      std::cout << "zinv: " << thisZinv << std::endl;
      exit(101);
    }


    for( int iBin=1; iBin<thisQCD->yield->GetXaxis()->GetNbins()+1; ++iBin ) {

      int ibin_low = h1_qcd->FindBin( thisQCD->yield->GetXaxis()->GetBinLowEdge(iBin) );
      int ibin_hi  = h1_qcd->FindBin( thisQCD->yield->GetXaxis()->GetBinLowEdge(iBin+1)*0.999 );

      if( ibin_low != ibin_hi ) continue; // the estimate bin needs to be included in the broad bin

      float yieldQCD  = thisQCD ->yield->GetBinContent(iBin);
      float yieldZinv = thisZinv->yield->GetBinContent(iBin);
      float yieldLlep = thisLlep->yield->GetBinContent(iBin);
      float yieldAll  = yieldQCD + yieldZinv + yieldLlep;
      
      float xref = thisQCD->yield->GetXaxis()->GetBinLowEdge(iBin);
      h1_qcd->Fill( xref, yieldQCD );
      h1_all->Fill( xref, yieldAll );

    } // for bins

  } // for regions

  TH1D* h1_qcdFrac = (TH1D*)(h1_qcd->Clone());
  h1_qcdFrac->SetName( Form("qcdFrac_%s", htName.c_str()) );
  h1_qcdFrac->Divide( h1_all );


  return h1_qcdFrac; 

}
