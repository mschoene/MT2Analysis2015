#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "THStack.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"



bool shapeNorm = true;


void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="" );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./drawDataMC [configFileName] [lumi/shape]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  MT2DrawTools::setStyle();

  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  if( argc>2 ) {
    std::string normType(argv[2]);
    if( normType=="lumi" ) shapeNorm=false;
    else if( normType=="shape" ) shapeNorm=true;
    else {
      std::cout << "-> Only 'lumi' and 'shape' are supported normTypes." << std::endl;
      exit(17);
    }
  }


  if( shapeNorm )
    std::cout << "-> Using shape normalization." << std::endl;
  else
    std::cout << "-> Using lumi normalization." << std::endl;


  std::string mcFile = cfg.getEventYieldDir() + "/analyses.root";
  std::string dataFile = cfg.getEventYieldDir() + "/analyses.root";

  MT2Analysis<MT2EstimateTree>* zjets = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "ZJets");
  zjets->setFullName("Z+Jets");
  MT2Analysis<MT2EstimateTree>* wjets = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "WJets");
  wjets->setFullName("W+Jets");
  MT2Analysis<MT2EstimateTree>* top   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "Top");
  top->setFullName("Top");
  MT2Analysis<MT2EstimateTree>* qcd   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "QCD");
  qcd->setFullName("QCD");

  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(dataFile, "data");
  data->setFullName("Data");

  std::vector< MT2Analysis<MT2EstimateTree>* > mc;
  mc.push_back(qcd);
  mc.push_back(wjets);
  mc.push_back(zjets);
  mc.push_back(top);


  std::string selection = "weight*(ht>900. && nJets>1 && met>30.)/puWeight";

  //drawYields( cfg, data, mc, "nVert_noW" , "nVert" , selection, 50, 0.5, 50.5, "Number of Vertices", "" );

  selection = "weight*(ht>450. && ht<900. && met>200. && nJets>1 && mt2>10. && id!=107 && deltaPhiMin>0.3 && diffMetMht<0.5*met)";

  drawYields( cfg, data, mc, "lowHT_nVert" , "nVert" , selection, 25, 0.5, 50.5, "Number of Vertices", "" );
  drawYields( cfg, data, mc, "lowHT_mt2"   , "mt2"   , selection, 18, 10., 910., "M_{T2}", "GeV" );
  drawYields( cfg, data, mc, "lowHT_met"   , "met"   , selection, 20, 200., 1200., "Missing E_{T}", "GeV" );
  drawYields( cfg, data, mc, "lowHT_ht"    , "ht"    , selection, 25, 450., 2950., "H_{T}", "GeV" );
  drawYields( cfg, data, mc, "lowHT_nJets" , "nJets" , selection, 12, 1.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "" );
  drawYields( cfg, data, mc, "lowHT_nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "" );
  drawYields( cfg, data, mc, "lowHT_deltaPhiMin", "deltaPhiMin", selection, 32, 0.0, 3.2, "min #Delta#phi(jets, ME_{T})", "" );
  drawYields( cfg, data, mc, "lowHT_diffMetMht_overMet", "diffMetMht/met", selection, 30, 0.0, 3.0, "|ME_{T}-MH_{T}|/ME_{T}", "" );
  drawYields( cfg, data, mc, "lowHT_diffMetMht", "diffMetMht", selection, 24, 0., 1200., "|ME_{T}-MH_{T}|", "GeV" );


  selection = "weight*(ht>900. && nJets>1 && met>30. && mt2>10. && deltaPhiMin>0.3 && diffMetMht<0.5*met)";

  drawYields( cfg, data, mc, "nVert" , "nVert" , selection, 50, 0.5, 50.5, "Number of Vertices", "" );
  drawYields( cfg, data, mc, "mt2"   , "mt2"   , selection, 60, 10., 310., "M_{T2}", "GeV" );
  drawYields( cfg, data, mc, "met"   , "met"   , selection, 80, 30., 430., "Missing E_{T}", "GeV" );
  drawYields( cfg, data, mc, "ht"    , "ht"    , selection, 50, 900., 3400., "H_{T}", "GeV" );
  drawYields( cfg, data, mc, "nJets" , "nJets" , selection, 12, 1.5, 13.5, "Number of Jets (p_{T} > 30 GeV)", "" );
  drawYields( cfg, data, mc, "nBJets", "nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "" );
  drawYields( cfg, data, mc, "deltaPhiMin", "deltaPhiMin", selection, 32, 0.0, 3.2, "min #Delta#phi(jets, ME_{T})", "" );
  drawYields( cfg, data, mc, "diffMetMht_overMet", "diffMetMht/met", selection, 30, 0., 3.0, "|ME_{T}-MH_{T}|/ME_{T}", "" );
  drawYields( cfg, data, mc, "diffMetMht", "diffMetMht", selection, 100, 0., 500., "|ME_{T}-MH_{T}|", "GeV" );

  drawYields( cfg, data, mc, "mt2_tail"   , "mt2"   , selection, 100, 0., 5000., "M_{T2}", "GeV" );
  drawYields( cfg, data, mc, "met_tail"   , "met"   , selection, 100,  0., 5000., "Missing E_{T}", "GeV" );
  drawYields( cfg, data, mc, "ht_tail"    , "ht"    , selection, 120, 900., 12900., "H_{T}", "GeV" );

  
  return 0;

}



void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units ) {


  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;



  std::vector<int> colors;
  if( bgYields.size()==3 ) { // estimates
    colors.push_back(402); 
    colors.push_back(430); 
    colors.push_back(418); 
  } else { // mc
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    colors.push_back(419); // z+jets
    colors.push_back(855); // top
    //colors.push_back(); // other
  }

  std::string fullPathPlots = cfg.getEventYieldDir() + "/plotsDataMC";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );


    TTree* tree_data = data->get(thisRegion)->tree;
    TH1D* h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
    tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );

    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.2);


    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<bgYields.size(); ++i ) { 
      TTree* tree_mc = (bgYields[i]->get(thisRegion)->tree);
      std::string thisName = "h1_" + bgYields[i]->getName();
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_mc->Sumw2();
      if( selection!="" )
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%s/puWeight", selection.c_str()) );
      //tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%s", selection.c_str()) );
      else
        tree_mc->Project( thisName.c_str(), varName.c_str(), "" );
      histos_mc.push_back(h1_mc);
    }

    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

    TH1D* mc_sum;
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      if( i==0 ) {
        mc_sum = new TH1D( *histos_mc[i] );
        mc_sum->SetName("mc_sum");
      } else {
        mc_sum->Add( histos_mc[i] );
      }
    }

    std::cout << "Integrals: " << h1_data->Integral(0, nBins+1) << "\t" << mc_sum->Integral(0, nBins+1) << std::endl;
    float scaleFactor = h1_data->Integral(0, nBins+1)/mc_sum->Integral(0, nBins+1);
    if( shapeNorm )
      std::cout << "SF: " << scaleFactor << std::endl;

    TH1D* histo_mc;
    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = bgYields.size() - i - 1;
      histos_mc[index]->SetFillColor( colors[index] );
      histos_mc[index]->SetLineColor( kBlack );
      if( shapeNorm )
        histos_mc[index]->Scale( scaleFactor );

      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack.Add(histos_mc[index]);
    }


    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();
    
    TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600);

    float yMaxScale = 1.1;
    float yMax1 = h1_data->GetMaximum()*yMaxScale;
    float yMax2 = yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
    float yMax3 = yMaxScale*(bgStack.GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    if( h1_data->GetNbinsX()<2 ) yMax *=3.;

    std::string xAxisTitle;
    if( units!="" ) 
      xAxisTitle = (std::string)(Form("%s [%s]", axisName.c_str(), units.c_str()) );
    else
      xAxisTitle = (std::string)(Form("%s", axisName.c_str()) );

    std::string yAxisTitle;
    if(binWidth>0.99){
      if( units!="" ) 
	yAxisTitle = (std::string)(Form("Events / (%.0f %s)", binWidth, units.c_str()));
      else
	yAxisTitle = (std::string)(Form("Events / (%.0f)", binWidth));
    }
    else{
      if( units!="" ) 
	yAxisTitle = (std::string)(Form("Events / (%.2f %s)", binWidth, units.c_str()));
      else
	yAxisTitle = (std::string)(Form("Events / (%.2f)", binWidth));
    }

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
    TPad* pad1 = MT2DrawTools::getCanvasMainPad();
    pad1->Draw();
    pad1->cd();
    h2_axes->Draw();

    
   
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.1, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();
    TPad* pad1_log = MT2DrawTools::getCanvasMainPad( true );
    pad1_log->Draw();
    pad1_log->cd();
    h2_axes_log->Draw();
   


    std::vector<std::string> niceNames = thisRegion.getNiceNames();

    for( unsigned i=0; i<niceNames.size(); ++i ) {

      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.035);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      regionText->AddText( niceNames[i].c_str() );

      //c1->cd();
      //regionText->Draw("same");
  
      //c1_log->cd();
      //regionText->Draw("same");
  
    }
    

    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.35, 0.8, 0.75, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      normText->AddText( "#splitline{Shape}{Norm.}" );
      pad1->cd();
      normText->Draw("same");
      pad1_log->cd();
      normText->Draw("same");
    }


    TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( gr_data, "Data", "P" );
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], bgYields[i]->getFullName().c_str(), "F" );
    }


    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    
    TPaveText* ratioText = new TPaveText( 0.133, -0.05, 0.4, 0.1 , "brNDC" );
    ratioText->SetTextSize(0.035);
    ratioText->SetTextFont(40);
    ratioText->SetTextColor(2);
    ratioText->SetFillColor(0);
    ratioText->SetTextAlign(11);
    ratioText->AddText( Form("Data/MC = %.2f", scaleFactor) );
    
    c1->cd();
    pad1->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
      ratioText->Draw("same");

    gPad->RedrawAxis();

    c1_log->cd();
    pad1_log->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
      ratioText->Draw("same");

    gPad->RedrawAxis();

    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    
    TLine* lineSF = new TLine(xMin, scaleFactor, xMax, scaleFactor);
    lineSF->SetLineColor(2);

    float yMinR=0.0;
    float yMaxR=2.0;

    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );

    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
    g_ratio->SetMarkerStyle(20);
    g_ratio->SetLineColor(1);
    g_ratio->SetLineWidth(2);

    c1->cd();
    TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
    pad2->Draw();
    pad2->cd();
    
    h2_axes_ratio->Draw("");
    g_ratio->Draw("PE,same");
    line->Draw("same");
    if( !shapeNorm )
      lineSF->Draw("same");
    
    gPad->RedrawAxis();


    c1_log->cd();
    TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
    pad2_log->Draw();
    pad2_log->cd();
    
    h2_axes_ratio->Draw("");
    g_ratio->Draw("PE,same");
    line->Draw("same");
    if( !shapeNorm )
      lineSF->Draw("same");

    gPad->RedrawAxis();


    c1->SaveAs( Form("%s/%s_%s.eps", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.png", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.pdf", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );

    c1_log->SaveAs( Form("%s/%s_%s_log.eps", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1_log->SaveAs( Form("%s/%s_%s_log.png", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );
    c1_log->SaveAs( Form("%s/%s_%s_log.pdf", fullPathPlots.c_str(), saveName.c_str(), thisRegion.getName().c_str()) );

    delete c1;
    delete h2_axes;

    delete c1_log;
    delete h2_axes_log;
    
    delete h2_axes_ratio;
    
    delete h1_data;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

  }// for MT2 regions

}
