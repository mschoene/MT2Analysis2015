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

  std::string outputdir = "/scratch/mmasciov/CMSSW_7_2_3_GammaFunctions/src/MT2Analysis2015/analysis/YieldComparison/";

  //std::string firstInputFile = "/shome/mmasciov/analyses_fp.root";
  std::string firstInputFile  = "/scratch/mmasciov/CMSSW_7_2_3_GammaFunctions/src/MT2Analysis2015/analysis/llep_PHYS14_v4_skimprune_zurich_4fb.root";
  std::string secondInputFile = "/scratch/mmasciov/CMSSW_7_2_3_GammaFunctions/src/MT2Analysis2015/analysis/llep_PHYS14_Zurich_MT2final_13TeV_PHYS14_loJet_hiHT_noMT_4fb.root";

  MT2Analysis<MT2EstimateSyst>* analysisFirst = MT2Analysis<MT2EstimateSyst>::readFromFile(firstInputFile.c_str(), "llep");
 
  MT2Analysis<MT2EstimateSyst>* analysisSecond = MT2Analysis<MT2EstimateSyst>::readFromFile( secondInputFile.c_str(), "llep" );

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
  colors.push_back(430);
  

  std::set<MT2Region> MT2Regions = analysisFirst->getRegions();

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::string fullPath = outputdir;

      TH1D* h_first = analysisFirst->get(*iMT2)->yield;

      TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), iMT2->getName().c_str()), "recreate" );
      histoFile->cd();
      h_first->Write();
      h_first->SetMarkerStyle(20);
      h_first->SetMarkerSize(1.6);

      TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
      c1->cd();

      float xMin = h_first->GetXaxis()->GetXmin();
      float xMax = h_first->GetXaxis()->GetXmax();
      float yMax1 = h_first->GetMaximum()*1.5;
      float yMax2 = 1.2*(h_first->GetMaximum() + h_first->GetBinError(h_first->GetMaximumBin()));
      float yMax = (yMax1>yMax2) ? yMax1 : yMax2;

      TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
      h2_axes->SetXTitle("M_{T2} [GeV]");
      h2_axes->SetYTitle("Entries");

      h2_axes->Draw();

      THStack bgStack("bgStack", "");
      for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
        int index = bgYields.size() - i - 1;
        TH1D* h_second = bgYields[index]->get(*iMT2)->yield;
        h_second->SetFillColor( colors[index] );
        h_second->SetLineColor( kBlack );
        bgStack.Add(h_second);
      }

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


      TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
      legend->AddEntry( h_first, "MC only", "P" );
      histoFile->cd();
      for( unsigned i=0; i<bgYields.size(); ++i ) {
        TH1D* h_second = bgYields[i]->get(*iMT2)->yield;
        legend->AddEntry( h_second, bgYields[i]->getName().c_str(), "F" );
        h_second->Write();
      }

      histoFile->Close();

      legend->Draw("same");
      bgStack.Draw("histo same");
      h_first->Draw("p, same");

      TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
      labelTop->Draw("same");

      gPad->RedrawAxis();

      c1->SaveAs( Form("%s/mt2_%s.eps", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.png", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;
      delete h2_axes;

  } // for MT2 regions

}
