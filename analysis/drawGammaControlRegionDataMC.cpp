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


double lumiErr = 0.12;
bool shapeNorm = false;


void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>* mc, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="", float htMin=450., float htMax=-1);


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




  std::string mcFile = cfg.getEventYieldDir() + "/gammaControlRegion/mc.root";
  std::string dataFile = cfg.getEventYieldDir() + "/gammaControlRegion/data.root";

  MT2Analysis<MT2EstimateTree>* mc   = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "gammaCRtree_loose");
  MT2Analysis<MT2EstimateTree>* data = MT2Analysis<MT2EstimateTree>::readFromFile(dataFile, "gammaCRtree_loose");

  float htMin=450, htMax=-1;

  std::string selection = "ptGamma>180. && ht>450. && nJets>1 && iso<2.5 && deltaPhiMin>0.3 && diffMetMht<0.5*met";

  drawYields( cfg, data, mc, "nVert"     , "nVert"    , selection, 25, 0.5, 50.5, "Number of Vertices", "", htMin, htMax );
  drawYields( cfg, data, mc, "mt2"       , "mt2"      , selection, 12, 0., 600., "M_{T2} (Photon Removed)", "GeV", htMin, htMax );
  drawYields( cfg, data, mc, "met"       , "met"      , selection, 18, 0., 900., "Missing E_{T}", "GeV", htMin, htMax );
  drawYields( cfg, data, mc, "ht"        , "ht"       , selection, 22, 450., 1550., "H_{T}", "GeV", htMin, htMax );
  drawYields( cfg, data, mc, "nJets"     , "nJets"    , selection, 10, 1.5, 11.5, "Number of Jets (p_{T} > 30 GeV)", "", htMin, htMax );
  drawYields( cfg, data, mc, "nBJets"    , "nBJets"   , selection, 6, -0.5, 5.5, "Number of b-Jets (p_{T} > 20 GeV)", "", htMin, htMax );
  drawYields( cfg, data, mc, "nJetHF30"  , "nJetHF30" , selection, 6, -0.5, 5.5, "N(jets, p_{T}>30 GeV & |#eta|>3.0)", "", htMin, htMax );
  drawYields( cfg, data, mc, "ptGamma"   , "ptGamma"  , selection, 18, 180., 1080, "Photon p_{T}", "GeV", htMin, htMax );
  drawYields( cfg, data, mc, "etaGamma"  , "etaGamma" , selection, 10, -2.5, 2.5, "Photon #eta", "", htMin, htMax );
  drawYields( cfg, data, mc, "sietaieta" , "sietaieta", selection, 12, 0.0075, 0.0111, "Photon #sigma_{i#eta i#eta}", "", htMin, htMax );
  drawYields( cfg, data, mc, "gamma_chHadIsoRC" , "gamma_chHadIsoRC[0]", selection, 20, 0., 10., "Charged Isolation RC", "GeV", htMin, htMax );
  
  //drawYields( cfg, data, mc, "iso"       , "iso"      , selection, 20, 0., 10., "Photon Charged Isolation", "GeV", htMin, htMax );

  selection = "ptGamma>180. && ht>450. && nJets>1 && deltaPhiMin>0.3 && diffMetMht<0.5*met";

  //drawYields( cfg, data, mc, "loose_nVert"     , "nVert"    , selection, 25, 0.5, 50.5, "Number of Vertices", "", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_mt2"       , "mt2"      , selection, 15, 0., 750., "M_{T2}", "GeV", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_met"       , "met"      , selection, 18, 0., 900., "Missing E_{T}", "GeV", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_ht"        , "ht"       , selection, 22, 450., 1550., "H_{T}", "GeV", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_nJets"     , "nJets"    , selection, 10, 1.5, 11.5, "Number of Jets (p_{T} > 30 GeV)", "", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_nBJets"    , "nBJets"   , selection, 6, -0.5, 5.5, "Number of b-Jets (p_{T} > 20 GeV)", "", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_jet1_pt"   , "jet1_pt"  , selection, 20, 30., 1030., "Leading Jet p_{T}", "GeV", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_jet2_pt"   , "jet2_pt"  , selection, 13, 30., 680., "Subleading Jet p_{T}", "GeV", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_ptGamma"   , "ptGamma"  , selection, 18, 180., 1080, "Photon p_{T}", "GeV", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_etaGamma"  , "etaGamma" , selection, 10, -2.5, 2.5, "Photon #eta", "", htMin, htMax );
  //drawYields( cfg, data, mc, "loose_sietaieta" , "sietaieta", selection, 12, 0.0075, 0.0111, "Photon #sigma_{i#eta i#eta}", "", htMin, htMax );
  drawYields( cfg, data, mc, "loose_iso"       , "iso"      , selection, 20, 0., 10., "Photon Charged Isolation", "GeV", htMin, htMax );

  
  return 0;

}



void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, MT2Analysis<MT2EstimateTree>*  mc, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units, float htMin, float htMax ) {


  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;



  std::vector<int> colors;
  colors.push_back(18); 
  colors.push_back(38); 
  colors.push_back(46); 


  std::string fullPathPlots = cfg.getEventYieldDir() + "/gammaControlRegion/plotsDataMC";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );


    TTree* tree_data = data->get(thisRegion)->tree;
    TH1D* h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
    tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );

    MT2DrawTools::addOverflowSingleHisto(h1_data);
    
    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.2);



    TH1D* h1_mc_prompt = new TH1D( Form("%s_prompt", mc->getName().c_str()), "", nBins, xMin, xMax );
    h1_mc_prompt->Sumw2();
    h1_mc_prompt->SetTitle("Prompt");
    TH1D* h1_mc_nip = new TH1D( Form("%s_nip", mc->getName().c_str()), "", nBins, xMin, xMax );
    h1_mc_nip->Sumw2();
    h1_mc_nip->SetTitle("Fragm.");
    TH1D* h1_mc_fake = new TH1D( Form("%s_fake", mc->getName().c_str()), "", nBins, xMin, xMax );
    h1_mc_fake->Sumw2();
    h1_mc_fake->SetTitle("Fakes");
    
    TTree* tree_mc = mc->get(thisRegion)->tree;

    if( selection!="" ) {
      tree_mc->Project( h1_mc_prompt->GetName(), varName.c_str(), Form("weight*(%s && prompt==2)", selection.c_str()) );
      tree_mc->Project( h1_mc_nip   ->GetName(), varName.c_str(), Form("weight*(%s && prompt==1)", selection.c_str()) );
      tree_mc->Project( h1_mc_fake  ->GetName(), varName.c_str(), Form("weight*(%s && prompt==0)", selection.c_str()) );
    } else {
      tree_mc->Project( h1_mc_prompt->GetName(), varName.c_str(), "weight*(prompt==2)" );
      tree_mc->Project( h1_mc_nip   ->GetName(), varName.c_str(), "weight*(prompt==1)" );
      tree_mc->Project( h1_mc_fake  ->GetName(), varName.c_str(), "weight*(prompt==0)" );
    }

    std::cout << h1_mc_prompt->GetBinContent(nBins) << std::endl;
    
    MT2DrawTools::addOverflowSingleHisto(h1_mc_prompt);
    MT2DrawTools::addOverflowSingleHisto(h1_mc_nip);
    MT2DrawTools::addOverflowSingleHisto(h1_mc_fake);

    std::cout << h1_mc_prompt->GetBinContent(nBins) << std::endl;

    std::vector< TH1D* > histos_mc;
    histos_mc.push_back( h1_mc_prompt );
    histos_mc.push_back( h1_mc_nip );
    histos_mc.push_back( h1_mc_fake );
    

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
    std::cout << "SF: " << scaleFactor << std::endl;

    double error_data;
    double integral_data = h1_data->IntegralAndError(0, nBins+1, error_data);

    double error_mc;
    double integral_mc = mc_sum->IntegralAndError(0, nBins+1, error_mc);

    double error_datamc = MT2DrawTools::getSFError(integral_data, error_data, integral_mc, error_mc);

    TH1D* histo_mc;
    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = histos_mc.size() - i - 1;
      histos_mc[index]->SetFillColor( colors[index] );
      histos_mc[index]->SetLineColor( kBlack );
      if( shapeNorm ) {
        histos_mc[index]->Scale( scaleFactor );
      } 
      else
	histos_mc[index]->Scale(16.1/20.38);

      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);

      bgStack.Add(histos_mc[index]);
    }

    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);

    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);

    TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TGraphErrors* SFFitBand = MT2DrawTools::getSFFitBand(fSF, xMin, xMax);

    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );


    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    TPad *pad1 = MT2DrawTools::getCanvasMainPad();
    
    TCanvas* c1_log = new TCanvas( "c1_log", "", 600, 600 );
    c1_log->cd();
    TPad *pad1_log = MT2DrawTools::getCanvasMainPad( true );
   
    float yMaxScale = 1.25;
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


    std::string binWidthText;
    if( binWidth>=1. )         binWidthText = (std::string)Form("%.0f", binWidth);
    else if( binWidth>=0.1 )   binWidthText = (std::string)Form("%.1f", binWidth);
    else if( binWidth>=0.01 )  binWidthText = (std::string)Form("%.2f", binWidth);
    else if( binWidth>=0.001 ) binWidthText = (std::string)Form("%.3f", binWidth);
    else                      binWidthText = (std::string)Form("%.4f", binWidth);

    std::string yAxisTitle;
    if( units!="" ) 
      yAxisTitle = (std::string)(Form("Events / (%s %s)", binWidthText.c_str(), units.c_str()));
    else
      yAxisTitle = (std::string)(Form("Events / (%s)", binWidthText.c_str()));

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
    pad1->Draw();
    pad1->cd();
    h2_axes->Draw();


   
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.1, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();
    pad1_log->Draw();
    pad1_log->cd();
    h2_axes_log->Draw();
   


    std::vector<std::string> niceNames = thisRegion.getNiceNames();

//    for( unsigned i=0; i<niceNames.size(); ++i ) {
//
//      float yMax = 0.9-(float)i*0.05;
//      float yMin = yMax - 0.05;
//      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
//      regionText->SetTextSize(0.040);
//      //      regionText->SetTextFont(42);
//      regionText->SetFillColor(0);
//      regionText->SetTextAlign(11);
//      regionText->AddText( niceNames[i].c_str() );
//
//      pad1->cd();
//      regionText->Draw("same");
//  
//      pad1_log->cd();
//      regionText->Draw("same");
//  
//    }

    for( unsigned i=0; i<niceNames.size(); ++i ) {

      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.030);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);

      if(i==0)
        if( htMax > 0. )
          regionText->AddText( Form("%.0f < H_{T} < %.0f GeV", htMin, htMax) );
        else
          regionText->AddText( Form("H_{T} > %.0f GeV", htMin) );
      else
        regionText->AddText( niceNames[i].c_str() );

      pad1->cd();
      regionText->Draw("same");

      pad1_log->cd();
      regionText->Draw("same");

    }


    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      normText->AddText( "#splitline{Shape}{Norm.}" );
      pad1->cd();
      normText->Draw("same");
      pad1_log->cd();
      normText->Draw("same");
    }

    

    TLegend* legend = new TLegend( 0.7, 0.9-(histos_mc.size()+1)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.040);
    //    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( gr_data, "Data", "P" );
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], histos_mc[i]->GetTitle(), "F" );
    }


    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());

    TPaveText* fitText = MT2DrawTools::getFitText( fSF );
    
    //TPaveText* ratioText = new TPaveText( 0.133, -0.05, 0.4, 0.1 , "brNDC" );
    //ratioText->SetTextSize(0.035);
    //ratioText->SetTextFont(40);
    //ratioText->SetTextColor(2);
    //ratioText->SetFillColor(0);
    //ratioText->SetTextAlign(11);
    //ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", scaleFactor, error_datamc) );
    
    float yMinR=0.0;
    float yMaxR=2.0;

    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );

    c1->cd();
    pad1->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    mcBand->Draw("E2 same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
      fitText->Draw("same");
    //    ratioText->Draw("same");

    gPad->RedrawAxis();

    c1_log->cd();
    pad1_log->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    mcBand->Draw("E2 same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
      fitText->Draw("same");
    //    ratioText->Draw("same");

    gPad->RedrawAxis();
    
    c1->cd();
    TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
    pad2->Draw();
    pad2->cd();

    h2_axes_ratio->Draw("");
    lineCentral->Draw("same");
    if( !shapeNorm ){

      systBand->Draw("3,same");
      lineCentral->Draw("same");

      SFFitBand->Draw("3,same");
      fSF->Draw("same");

      //      SFband->Draw("3,same");
      //      lineSF->Draw("same");

    }
    g_ratio->Draw("PE,same");
    gPad->RedrawAxis();
    
    c1_log->cd();
    TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
    pad2_log->Draw();
    pad2_log->cd();

    h2_axes_ratio->Draw("");
    lineCentral->Draw("same");
    if( !shapeNorm ){

      systBand->Draw("3,same");
      lineCentral->Draw("same");

      SFFitBand->Draw("3,same");
      fSF->Draw("same");

      //      SFband->Draw("3,same");
      //      lineSF->Draw("same");

    }
    g_ratio->Draw("PE,same");
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

    delete h1_data;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

  }// for MT2 regions

}




