#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"



float lumi = 4.;

int type = 1;

// 0: use bare MC for ratio (pure GJet)
// 1: use GJet+QCD and multiply by fitted purity





void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2Estimate>* analysis, const std::string& suffix="" );
void drawComparison( const std::string& outputdir, std::string plotName, const std::string& axisName, const std::string& units, MT2Region region, MT2Analysis<MT2Estimate>* ana1, const std::string& legendName1, MT2Analysis<MT2Estimate>* ana2=0, const std::string& legendName2="" );
void drawHistos( const std::string& outputdir, std::string plotName, const std::string& axisName, const std::string& units, const MT2Region& region, TH1D* h1_ana1, const std::string& legendName1, TH1D* h1_ana2=0, const std::string& legendName2=""  );
void drawFromTree( const std::string& outputdir, const std::string& varName, int nBins, float xMin, float xMax, const std::string& axisName, const std::string& units, MT2Region region, MT2Analysis<MT2EstimateTree>* ana1, const std::string& legendName1, MT2Analysis<MT2EstimateTree>* ana2=0, const std::string& legendName2="" );


int main( int argc, char* argv[] ) {


  std::string samples = "PHYS14_v2_Zinv";

  std::string regionsSet = "zurich";
  if( argc>1 ) {
    regionsSet = std::string(argv[1]); 
  }



  std::string dir_estimate( Form("ZinvEstimateFromGamma_%s_%s_%.0ffb_type%d", samples.c_str(), regionsSet.c_str(), lumi, type) );
  //std::string dir_gammaCR = "GammaControlRegion_" + samples + "_" + regionsSet;


  MT2DrawTools::setStyle();

  std::string outputdir = dir_estimate + "/plots";
  system( Form("mkdir -p %s", outputdir.c_str()) );

  MT2Analysis<MT2Estimate>* zgammaRatio = MT2Analysis<MT2Estimate>::readFromFile(dir_estimate + "/MT2ZinvEstimate.root", "ZgammaRatio");
  MT2Analysis<MT2EstimateTree>* zinv = MT2Analysis<MT2EstimateTree>::readFromFile(dir_estimate + "/MT2ZinvEstimate.root", "Zinv");
  MT2Analysis<MT2Estimate>* zinvEstimate = MT2Analysis<MT2Estimate>::readFromFile(dir_estimate + "/MT2ZinvEstimate.root", "ZinvEstimate" );
  //MT2Analysis<MT2EstimateTree>* gammaJet = MT2Analysis<MT2EstimateTree>::readFromFile(dir_gammaCR + "/mc.root", "gammaJet");


  std::vector<MT2Region> r_lowHT_vs_njet;
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2, 3, 0, 0 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, 6, 0, 0 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2, 3, 1, 1 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, 6, 1, 1 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 7,-1, 0, 0 ) );

  compareRegions( outputdir, r_lowHT_vs_njet, zgammaRatio );



  std::vector<MT2Region> r_medHT_vs_njet;
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2, 3, 0, 0 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, 6, 0, 0 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2, 3, 1, 1 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, 6, 1, 1 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 7,-1, 0, 0 ) );

  compareRegions( outputdir, r_medHT_vs_njet, zgammaRatio );



  std::vector<MT2Region> r_highHT_vs_njet;
  r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 2, 3, 0, 0 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 4, 6, 0, 0 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 2, 3, 1, 1 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 4, 6, 1, 1 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 7,-1, 0, 0 ) );

  compareRegions( outputdir, r_highHT_vs_njet, zgammaRatio );



  std::vector<MT2Region> r_veryhighHT_vs_njet;
  r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 2, 3, 0, 0 ) );
  r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 4, 6, 0, 0 ) );
  r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 2, 3, 1, 1 ) );
  r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 4, 6, 1, 1 ) );
  r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 7,-1, 0, 0 ) );

  compareRegions( outputdir, r_veryhighHT_vs_njet, zgammaRatio );



  std::vector<MT2Region> r_vsHT;
  r_vsHT.push_back( MT2Region( 450., 575. , 2, 3, 0, 0 ) );
  r_vsHT.push_back( MT2Region( 575., 1000., 2, 3, 0, 0 ) );
  r_vsHT.push_back( MT2Region( 1000.,1500., 2, 3, 0, 0 ) );
  r_vsHT.push_back( MT2Region( 1500., -1. , 2, 3, 0, 0 ) );

  compareRegions( outputdir, r_vsHT, zgammaRatio );




  std::vector<MT2Region> r_vsHT_j46;
  r_vsHT_j46.push_back( MT2Region( 450., 575. , 4, 6, 0, 0 ) );
  r_vsHT_j46.push_back( MT2Region( 575., 1000., 4, 6, 0, 0 ) );
  r_vsHT_j46.push_back( MT2Region( 1000.,1500., 4, 6, 0, 0 ) );
  r_vsHT_j46.push_back( MT2Region( 1500., -1. , 4, 6, 0, 0 ) );

  compareRegions( outputdir, r_vsHT_j46, zgammaRatio );



  std::vector<MT2Region> r_vsHT_b1;
  r_vsHT_b1.push_back( MT2Region( 450., 575. , 2, 3, 1, 1 ) );
  r_vsHT_b1.push_back( MT2Region( 575., 1000., 2, 3, 1, 1 ) );
  r_vsHT_b1.push_back( MT2Region( 1000.,1500., 2, 3, 1, 1 ) );
  r_vsHT_b1.push_back( MT2Region( 1500., -1. , 2, 3, 1, 1 ) );

  compareRegions( outputdir, r_vsHT_b1, zgammaRatio );






  std::set<MT2Region> regions = zinv->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
    drawComparison( outputdir, "closure"             , "M_{T2}", "GeV", *iR, (MT2Analysis<MT2Estimate>*)zinv, "Z#rightarrow#nu#nu (MC)", zinvEstimate, "Z Invisible Estimate" );
    //drawComparison( outputdir, "Z_vs_gamma"          , "M_{T2}", "GeV", *iR, (MT2Analysis<MT2Estimate>*)zinv, "Z#rightarrow#nu#nu", (MT2Analysis<MT2Estimate>*)gammaJet, "#gamma+jet" );
    //drawComparison( outputdir, "ZgammaRatio"         , "M_{T2}", "GeV", *iR, zgammaRatio, "Z#rightarrow#nu#nu / #gamma+jet" );
    //drawFromTree  ( outputdir, "ptV", 100, 0., 1500., "p_{T}(V)", "GeV", *iR, zinv, "Z#rightarrow#nu#nu", gammaJet, "#gamma+jet" );
  }





  return 0;

}



void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2Estimate>* analysis, const std::string& suffix ) {

  bool loopOnHT=true;

  if( regions.size()>1 ) 
    if( (*(regions[0].htRegion())) == (*(regions[1].htRegion())) ) loopOnHT=false;



  std::vector<int> colors;
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 38 );
  colors.push_back( 42 );
  colors.push_back( kRed );
  colors.push_back( kBlack );
  
  std::vector<int> markers;
  markers.push_back( 21 );
  markers.push_back( 20 );
  markers.push_back( 25 );
  markers.push_back( 24 );
  markers.push_back( 26 );
  
  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D( "axes", "", 10, 200., 1500., 10, 0., 1.5 );
  h2_axes->SetXTitle( "M_{T2} [GeV]" );
  h2_axes->SetYTitle( "Z/#gamma Ratio" );
  h2_axes->Draw("");
  

  std::string legendTitle = (loopOnHT) ? regions[0].sigRegion()->getNiceName() : regions[0].htRegion()->getNiceName();

  TLegend* legend = new TLegend( 0.2, 0.9-0.07*regions.size(), 0.55, 0.9, legendTitle.c_str() );
  legend->SetTextSize(0.038);
  legend->SetFillColor(0);
  

  for( unsigned i=0; i<regions.size(); ++i ) {

    MT2Estimate* thisEstimate = analysis->get( regions[i] );
    if( thisEstimate==0 ) {
      std::cout << "-> Didn't find estimate for region: " << regions[i].getName() << " ! Skipping." << std::endl;
      delete c1;
      delete h2_axes;
      return;
    }

    TH1D* thisRatio = thisEstimate->yield;

    thisRatio->SetMarkerSize(1.);
    thisRatio->SetMarkerStyle(markers[i]);
    thisRatio->SetMarkerColor(colors[i]);
    //thisRatio->SetLineColor(kBlack);

    thisRatio->SetLineWidth(2);
    thisRatio->SetLineColor(colors[i]);

    thisRatio->Draw("same" );

    if( loopOnHT )
      legend->AddEntry( thisRatio, regions[i].htRegion()->getNiceName().c_str(), "PL" );
    else
      legend->AddEntry( thisRatio, regions[i].sigRegion()->getNiceName().c_str(), "PL" );

  }


  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  std::string saveName = (loopOnHT) ? regions[0].sigRegion()->getName() : regions[0].htRegion()->getName();

  c1->SaveAs( Form("%s/ZgammaRatio_%s%s.eps", outputdir.c_str(), saveName.c_str(), suffix.c_str()) );
  c1->SaveAs( Form("%s/ZgammaRatio_%s%s.pdf", outputdir.c_str(), saveName.c_str(), suffix.c_str()) );
  c1->SaveAs( Form("%s/ZgammaRatio_%s%s.png", outputdir.c_str(), saveName.c_str(), suffix.c_str()) );

  delete c1;
  delete h2_axes;

}




void drawComparison( const std::string& outputdir, std::string plotName, const std::string& axisName, const std::string& units, MT2Region region, MT2Analysis<MT2Estimate>* ana1, const std::string& legendName1, MT2Analysis<MT2Estimate>* ana2, const std::string& legendName2 ) {


  TH1D* h1_ana1 = ana1->get( region )->yield;
  TH1D* h1_ana2 = 0;
  if( ana2!=0 ) h1_ana2 = ana2->get( region )->yield;


  drawHistos( outputdir, plotName, axisName, units, region, h1_ana1, legendName1, h1_ana2, legendName2 );

}




void drawHistos( const std::string& outputdir, std::string plotName, const std::string& axisName, const std::string& units, const MT2Region& region, TH1D* h1_ana1, const std::string& legendName1, TH1D* h1_ana2, const std::string& legendName2  ) {


  bool use_lumi = true;

  TString ana1_name(h1_ana1->GetName());
  TString ana2_name;
  if( h1_ana2!=0 ) ana2_name = TString(h1_ana2->GetName());
  bool isClosurePlot = h1_ana2!=0 && ana2_name.Contains("ZinvEstimate");
  bool isRatioPlot = ana1_name.Contains("Ratio");

  if( h1_ana2!=0 ) {

    if( isClosurePlot ) {

      h1_ana1->SetLineColor( 46 );
      h1_ana1->SetLineWidth( 2 );

      h1_ana2->SetLineColor( kBlack );
      h1_ana2->SetMarkerColor( kBlack );
      h1_ana2->SetMarkerStyle( 20 );
      h1_ana2->SetMarkerSize( 1.3 );

    } else {

      h1_ana1->SetLineColor( kBlack );
      h1_ana1->SetMarkerColor( kBlack );
      h1_ana1->SetMarkerStyle( 21 );
      h1_ana1->SetMarkerSize( 1.2 );

      h1_ana2->SetLineColor( 46 );
      h1_ana2->SetMarkerColor( 46 );
      h1_ana2->SetMarkerStyle( 20 );
      h1_ana2->SetMarkerSize( 1.3 );

    }

  } else {

    h1_ana1->SetLineColor( 46 );
    h1_ana1->SetMarkerColor( 46 );
    h1_ana1->SetMarkerStyle( 20 );
    h1_ana1->SetMarkerSize( 1.3 );

  }


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  TCanvas* c1_log = new TCanvas( "c1_log", "", 600, 600 );
  c1_log->SetLogy();

  float xMin = h1_ana1->GetXaxis()->GetXmin();
  float xMax = h1_ana1->GetXaxis()->GetXmax();
  float yMax1 = h1_ana1->GetMaximum()*1.5;
  float yMax2 = (h1_ana2!=0) ? h1_ana2->GetMaximum()*1.5 : 0.;
  float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
  if( isRatioPlot ) {
    yMax=2.;
    plotName = "ratio_" + plotName;
    use_lumi = false;
  }
  
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
  if( units!="" )
    h2_axes->SetXTitle(Form("%s [%s]", axisName.c_str(), units.c_str()) );
  else
    h2_axes->SetXTitle(axisName.c_str());
  if( h1_ana2!=0 )
    h2_axes->SetYTitle("Entries");
  else
    h2_axes->SetYTitle(legendName1.c_str());

  TH2D* h2_axes_log = new TH2D( "axes_log", "", 10, xMin, xMax, 10, 0.01, 50.*yMax );
  if( units!="" )
    h2_axes_log->SetXTitle(Form("%s [%s]", axisName.c_str(), units.c_str()) );
  else
    h2_axes_log->SetXTitle(axisName.c_str());
  if( h1_ana2!=0 )
    h2_axes_log->SetYTitle("Entries");
  else
    h2_axes_log->SetYTitle(legendName1.c_str());

  c1->cd();
  h2_axes->Draw();

  c1_log->cd();
  h2_axes_log->Draw("");

  bool jetCuts = region.sigRegion()->getName() != "j2toInf_b0toInf";

  int nEntries = 2;
  if( h1_ana2!=0 ) nEntries++;
  if( jetCuts!=0 ) nEntries++;


  TLegend* legend;
  if( h1_ana2!=0 ) {
    if( jetCuts ) legend = new TLegend( 0.35, 0.9-(float)nEntries*0.06, 0.93, 0.9, Form( "#splitline{%s}{%s}", region.htRegion()->getNiceName().c_str(), region.sigRegion()->getNiceName().c_str() ) );
    else          legend = new TLegend( 0.35, 0.9-(float)nEntries*0.06, 0.93, 0.9, Form( "%s", region.htRegion()->getNiceName().c_str()) );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    if( isClosurePlot ) legend->AddEntry( h1_ana1, legendName1.c_str(), "L" );
    else                legend->AddEntry( h1_ana1, legendName1.c_str(), "P" );
    legend->AddEntry( h1_ana2, legendName2.c_str(), "P" );
  }


  TPaveText* labelTop;
  if( use_lumi )
    labelTop = MT2DrawTools::getLabelTop(lumi);
  else
    labelTop = MT2DrawTools::getLabelTop();

  TLine* line_one = new TLine( xMin, 1., xMax, 1. );

  c1->cd();
  if( isClosurePlot ) { 
    h1_ana1->Draw("histo same");
  } else {
    h1_ana1->Draw("p same");
  }
  if( isRatioPlot )
    line_one->Draw("same");
  if( h1_ana2!=0 ) {
    h1_ana2->Draw("p same");
    legend->Draw("same");
  }
  labelTop->Draw("same");
  gPad->RedrawAxis();

  c1_log->cd();
  if( isClosurePlot ) { 
    h1_ana1->Draw("histo same");
  } else {
    h1_ana1->Draw("p same");
  }
  if( isRatioPlot )
    line_one->Draw("same");
  if( h1_ana2!=0 ) {
    h1_ana2->Draw("p same");
    legend->Draw("same");
  }
  labelTop->Draw("same");
  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s_%s.eps", outputdir.c_str(), plotName.c_str(), region.getName().c_str()) );
  c1->SaveAs( Form("%s/%s_%s.png", outputdir.c_str(), plotName.c_str(), region.getName().c_str()) );
  c1->SaveAs( Form("%s/%s_%s.pdf", outputdir.c_str(), plotName.c_str(), region.getName().c_str()) );

  c1_log->SaveAs( Form("%s/%s_%s_log.eps", outputdir.c_str(), plotName.c_str(), region.getName().c_str()) );
  c1_log->SaveAs( Form("%s/%s_%s_log.png", outputdir.c_str(), plotName.c_str(), region.getName().c_str()) );
  c1_log->SaveAs( Form("%s/%s_%s_log.pdf", outputdir.c_str(), plotName.c_str(), region.getName().c_str()) );


  delete c1;
  delete c1_log;
  delete h2_axes;
  delete h2_axes_log;

  if( isClosurePlot ) {

    TH1D* h1_ratio = new TH1D( *h1_ana1 );
    h1_ratio->SetName(Form("h1_Ratio_vs_%s", plotName.c_str()));
    h1_ratio->Divide( h1_ana2 );

    drawHistos( outputdir, plotName, axisName, units, region, h1_ratio, std::string( legendName1 + " / " + legendName2 ) ); 

  }
    
    
}





void drawFromTree( const std::string& outputdir, const std::string& varName, int nBins, float xMin, float xMax, const std::string& axisName, const std::string& units, MT2Region region, MT2Analysis<MT2EstimateTree>* ana1, const std::string& legendName1, MT2Analysis<MT2EstimateTree>* ana2, const std::string& legendName2 )  {

  // need this = true here. love these root features.
  TH1F::AddDirectory(kTRUE);

  TH1D* h1_ana1 = new TH1D( Form("h1_%s", ana1->getName().c_str()), "", nBins, xMin, xMax );
  h1_ana1->Sumw2();
  TH1D* h1_ana2 = new TH1D( Form("h1_%s", ana2->getName().c_str()), "", nBins, xMin, xMax );
  h1_ana2->Sumw2();

  ana1->get( region )->tree->Project( h1_ana1->GetName(), varName.c_str(), "weight" );
  ana2->get( region )->tree->Project( h1_ana2->GetName(), varName.c_str(), "weight" );

  drawHistos( outputdir, varName, axisName, units, region, h1_ana1, legendName1, h1_ana2, legendName2 );

  TH1D* h1_ratio = new TH1D( *h1_ana1 );
  h1_ratio->SetName(Form("h1_Ratio_vs_%s", varName.c_str()));
  h1_ratio->Divide( h1_ana2 );
  drawHistos( outputdir, varName, axisName, units, region, h1_ratio, "Z#rightarrow#nu#nu / #gamma+jet" );


  delete h1_ana1;
  delete h1_ana2;

  // setting back to false
  TH1F::AddDirectory(kFALSE);

}
