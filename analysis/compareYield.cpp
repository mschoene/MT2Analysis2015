#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2DrawTools.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>


#include "TTreeFormula.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"

float lumi =4.0; //fb-1 

void drawYields( const std::string& outputdir, MT2Analysis<MT2EstimateSyst>* data, std::vector<MT2Analysis<MT2EstimateSyst>* > bgYields );

int main() {

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = "./YieldComparison_30GeV_T2qq/";

  //std::string firstInputFile = "/shome/mmasciov/analyses_fp.root";
  std::string firstInputFile  = "./EventYields_mc_PHYS14_v6_dummy_4fb/analyses.root";
  std::string secondInputFile = "./EventYields_mc_PHYS14_v5_dummy_4fb_v0/analyses.root";

  MT2Analysis<MT2EstimateSyst>* analysisFirst = MT2Analysis<MT2EstimateSyst>::readFromFile(firstInputFile.c_str(), "SMS_T2qq_2J_mStop600_mLSP550");
 
  MT2Analysis<MT2EstimateSyst>* analysisSecond = MT2Analysis<MT2EstimateSyst>::readFromFile( secondInputFile.c_str(), "SMS_T2qq_2J_mStop600_mLSP550" );

  std::vector< MT2Analysis<MT2EstimateSyst> *> bgEstimate;
  bgEstimate.push_back( analysisSecond );

  drawYields(outputdir.c_str(), analysisFirst, bgEstimate );

  return 0;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2EstimateSyst>*  analysisFirst, std::vector< MT2Analysis<MT2EstimateSyst> *> bgYields ) {

  MT2DrawTools::setStyle();

  system(Form("mkdir -p %s", outputdir.c_str()));

  std::vector<int> colors;
//  if( bgYields.size()==3 ) { // estimates
//    colors.push_back(402);
//    colors.push_back(430);
//    colors.push_back(418);
//  } else { // mc
//    colors.push_back(401); // qcd
//    colors.push_back(417); // w+jets
//    colors.push_back(419); // z+jets
//    colors.push_back(855); // top
//    //colors.push_back(); // other
//  }
  colors.push_back(855);
  

  std::set<MT2Region> MT2Regions = analysisFirst->getRegions();

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::string fullPath = outputdir;

      TH1D* h_first = analysisFirst->get(*iMT2)->yield;

      TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), iMT2->getName().c_str()), "recreate" );
      histoFile->cd();
      h_first->Write();
      h_first->SetMarkerStyle(20);
      h_first->SetMarkerSize(1.6);
      h_first->SetLineColor( kBlack );
      h_first->SetMarkerColor( kBlack );
      
      h_first->Sumw2();
      
      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();
//      TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
//      pad1->SetBottomMargin(0);
//      pad1->Draw();
//      pad1->cd();
      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.15);
      pad1->Draw();
      pad1->cd();

      THStack bgStack("bgStack", "");
      TH1D* h_second = bgYields[0]->get(*iMT2)->yield;
      h_second->Sumw2();
      for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
        int index = bgYields.size() - i - 1;
        TH1D* h_second_ = bgYields[index]->get(*iMT2)->yield;
        h_second_->SetFillColor( colors[index] );
        h_second_->SetLineColor( kBlack );
	h_second_->Sumw2();
        bgStack.Add(h_second_);
	if( i>0 )
	  h_second->Add(h_second_);
      }



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



      std::vector<std::string> niceNames = iMT2->getNiceNames();

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
      //legend->AddEntry( h_first, "MC only", "P" );
      legend->AddEntry( h_first, "p_{T}(jet) > 30 GeV", "P" );
      legend->AddEntry( h_second, "p_{T}(jet) > 40 GeV", "F" );
//      histoFile->cd();
//      for( unsigned i=0; i<bgYields.size(); ++i ) {
//        TH1D* h_second = bgYields[i]->get(*iMT2)->yield;
//        legend->AddEntry( h_second, bgYields[i]->getName().c_str(), "F" );
//        h_second->Write();
//      }
//
//      histoFile->Close();

      legend->Draw("same");
      bgStack.Draw("histo, same");
      h_first->Draw("p, same");

      TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
      labelTop->Draw("same");

      gPad->RedrawAxis();

      c1->cd();
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.1);
      pad2->Draw();
      pad2->cd();
      //TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
      //pad2->SetTopMargin(0);
      //pad2->Draw();
      //pad2->cd();
      
      TH1D* h_ratio = (TH1D*) h_first->Clone("h_ratio");
      h_ratio->SetStats(0);
      h_ratio->Divide(h_second);
      h_ratio->SetMarkerStyle(21);
      h_ratio->SetMarkerSize(0.02);
      h_ratio->GetXaxis()->SetLabelSize(0.00);
      h_ratio->GetXaxis()->SetTickLength(0.09);
      h_ratio->GetYaxis()->SetNdivisions(5,5,0);
      h_ratio->GetYaxis()->SetRangeUser(0.,2.);
      h_ratio->GetYaxis()->SetTitleSize(0.17);
      h_ratio->GetYaxis()->SetTitleOffset(0.4);
      h_ratio->GetYaxis()->SetLabelSize(0.17);
      h_ratio->GetYaxis()->SetTitle("ratio");
      
//      h_ratio->GetXaxis()->SetTitle("M_{T2} [GeV]");
//      h_ratio->GetXaxis()->SetTitleSize(0.1);
//      h_ratio->GetXaxis()->SetTitleOffset(0.5);
//      h_ratio->GetXaxis()->SetLabelSize(0.05);
//      h_ratio->GetYaxis()->SetTitle("ratio");
//      h_ratio->GetYaxis()->SetTitleSize(0.1);
//      h_ratio->GetYaxis()->SetTitleOffset(0.5);
//      h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);

      h_ratio->Draw("ep");

      gPad->RedrawAxis();
      
      c1->cd();

      c1->SaveAs( Form("%s/mt2_%s.pdf", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.png", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;
      delete h2_axes;

  } // for MT2 regions

}
