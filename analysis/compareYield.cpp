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
#include "TLatex.h"
#include "TGraphAsymmErrors.h"

float lumi =20.0; //fb-1 

void drawYields( const std::string& outputdir, MT2Analysis<MT2EstimateSyst>* data, std::vector<MT2Analysis<MT2EstimateSyst>* > bgYields );

int main() {

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = "./YieldComparison_InvisibleZ/";

  std::string firstInputFile  = "./EventYields_mc_PHYS14_v6_2_dummy_binbybin/zinvFromGamma.root";
  std::string secondInputFile = "./EventYields_mc_PHYS14_v6_2_dummy_integral/zinvFromGamma.root";

  MT2Analysis<MT2EstimateSyst>* analysisFirst = MT2Analysis<MT2EstimateSyst>::readFromFile( firstInputFile.c_str(), "ZinvEstimate");
  MT2Analysis<MT2EstimateSyst>* analysisSecond = MT2Analysis<MT2EstimateSyst>::readFromFile( secondInputFile.c_str(), "ZinvEstimate" );

  std::vector< MT2Analysis<MT2EstimateSyst> *> bgEstimate;
  bgEstimate.push_back( analysisSecond );

  drawYields(outputdir.c_str(), analysisFirst, bgEstimate );

  return 0;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2EstimateSyst>*  analysisFirst, std::vector< MT2Analysis<MT2EstimateSyst> *> bgYields ) {

  MT2DrawTools::setStyle();

  system(Form("mkdir -p %s", outputdir.c_str()));

  std::vector<int> colors;
  for(unsigned int b=0; b<bgYields.size(); ++b)
    colors.push_back( bgYields[b]->getColor() );
  
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
  hestimate->SetFillColor(colors[0]);
  hestimate->SetLineColor(kBlack);

  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::string fullPath = outputdir;

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      TGraphAsymmErrors* g_first = analysisFirst->get(*iMT2)->getGraph();
      TH1D* h_first = analysisFirst->get(*iMT2)->yield;
      TH1D* h_firstUp = analysisFirst->get(*iMT2)->yield_systUp;
      TH1D* h_firstDn = analysisFirst->get(*iMT2)->yield_systDown;
     
      TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), iMT2->getName().c_str()), "recreate" );
      histoFile->cd();
      h_first->Write();
      
      g_first->SetMarkerStyle(20);
      g_first->SetMarkerSize(1.6);
      g_first->SetLineColor( kBlack );
      g_first->SetMarkerColor( kBlack );
 
      hdata->SetBinContent(iRegion, h_first->Integral());
      hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
      
      float thisError_first = ( fabs( h_first->Integral() - h_firstUp->Integral() ) >  fabs( h_first->Integral() - h_firstDn->Integral() ) ) ? fabs( h_first->Integral() - h_firstUp->Integral() ) : fabs( h_first->Integral() - h_firstDn->Integral() );
      hdata->SetBinError(iRegion, thisError_first);
      
      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.15);
      pad1->Draw();
      pad1->cd();

      //THStack bgStack("bgStack", "");
      TH1D* h_second = bgYields[0]->get(*iMT2)->yield;
      TH1D* h_secondUp = bgYields[0]->get(*iMT2)->yield_systUp;
      TH1D* h_secondDn = bgYields[0]->get(*iMT2)->yield_systDown;

      TGraphAsymmErrors* g_second = bgYields[0]->get(*iMT2)->getGraph();
      g_second->SetLineColor( colors[0] );
      g_second->SetMarkerColor( colors[0] );
      g_second->SetMarkerStyle(20);
      g_second->SetMarkerSize(1.6);
      g_second->Write();
      
      hestimate->SetBinContent(iRegion, h_first->Integral());
      hestimate->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );

      float thisError_second = ( fabs( h_second->Integral() - h_secondUp->Integral() ) >  fabs( h_second->Integral() - h_secondDn->Integral() ) ) ? fabs( h_second->Integral() - h_secondUp->Integral() ) : fabs( h_second->Integral() - h_secondDn->Integral() );
      hestimate->SetBinError(iRegion, thisError_second);

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
      legend->AddEntry( g_first, "Bin by Bin", "P" );
      legend->AddEntry( g_second, "Integration", "P" );

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
      h_ratio->Draw("histo");

      gPad->RedrawAxis();
      
      c1->cd();

      c1->SaveAs( Form("%s/mt2_%s.pdf", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.png", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;
      delete h2_axes;
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
  hestimate->Draw("hist");
  hdata->Draw("pe,same");

  TLegend* legend = new TLegend( 0.8, 0.9-(bgYields.size()+1)*0.06-0.06, 0.93, 0.9-0.06 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( hdata, "Bin by Bin", "P" );
  legend->AddEntry( hestimate, "Integration", "F" );

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

  std::string thisName = Form("%s_ratio", hdata->GetName());
  TH1D* h_ratio = (TH1D*) hdata->Clone(thisName.c_str());
  h_ratio->Divide(hestimate);
  h_ratio->Write();
  h_ratio->SetStats(0);
  h_ratio->SetMarkerStyle(20);
  h_ratio->SetLineColor(1);
  h_ratio->GetXaxis()->SetLabelSize(0.00);
  h_ratio->GetXaxis()->SetTickLength(0.09);
  h_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  h_ratio->GetYaxis()->SetTitleSize(0.17);
  h_ratio->GetYaxis()->SetTitleOffset(0.4);
  h_ratio->GetYaxis()->SetLabelSize(0.17);
  h_ratio->GetYaxis()->SetTitle("Ratio");

  h_ratio->SetLineWidth(2);
  h_ratio->Draw("histo");
  
  //  TLine* lHT[3];
  for( int iHT=1; iHT < 4; iHT++ ){
//    lHT[iHT-1] = new TLine(11*iHT,0.0, 11*iHT, 2.0 );
//    lHT[iHT-1]->SetLineColor(kBlack);
//    lHT[iHT-1]->SetLineStyle(3);
//    lHT[iHT-1]->SetLineWidth(2);

    lHT[iHT-1]->Draw("same");
  }


  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/mt2_fullEstimate.pdf", outputdir.c_str()) );
  
}
