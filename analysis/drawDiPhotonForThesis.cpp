
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGaxis.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"

#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"


#include <iostream>
#include "string.h"

using namespace std;

#define mt2_cxx
#include "interface/mt2.h"



double lumiErr = 0.026;
bool shapeNorm = false;

void drawShape( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, const std::string& saveName, const std::string& varName, const std::string& selection, const std::string& selection2, const std::string& selection3, unsigned int nSel, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="" );

void drawStacks(std::string fullPath, float *binss, unsigned int size,  std::string name, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  MT2Analysis<MT2EstimateTree>* data,const MT2Region thisRegion, std::string cut, float lumi);

void drawYields( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="" );

std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units );


int main(int argc, char* argv[]){

  std::string regionsSet = "zurich";


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

  if( argc<2 ) {
    std::cout << "USAGE: ./drawZllControlRegion [configFileName] (lumi/shape)" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  std::string configFileName(argv[1]);
  MT2Config cfg( configFileName);
  regionsSet = cfg.regionsSet();

  //  std::string outputdir = cfg.getEventYieldDir() + "/zllPurity";

  std::cout << "-> Using regions: " << regionsSet << std::endl;
  MT2DrawTools::setStyle();
 
  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  double intpart;
  double fracpart = modf(cfg.lumi(), &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("_%.0ffb", intpart ) );
  

  //  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string ggDir = cfg.getEventYieldDir() + "/diPhotonControlRegion/";
 

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this


  MT2Analysis<MT2EstimateTree>* qcd = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_qcd.root", ggDir.c_str()  ), "qcd");
  qcd->setColor(kQCD);
  qcd->setFullName("QCD");

  MT2Analysis<MT2EstimateTree>* diPhoton = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_diPhoton.root", ggDir.c_str() ), "diPhoton");
  diPhoton->setColor(kPink-8);
  diPhoton->setFullName("#gamma#gamma");
  
  MT2Analysis<MT2EstimateTree>* gjets = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_gjets.root", ggDir.c_str() ), "gjets");
  gjets->setColor(kPink-4);
  gjets->setFullName("#gamma+jets");

  MT2Analysis<MT2EstimateTree>* higgs = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_higgs.root", ggDir.c_str() ), "higgs");
  higgs->setColor(kCyan-6);
  higgs->setFullName("SM H");





  MT2Analysis<MT2EstimateTree>* sig = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_T2bH.root", ggDir.c_str() ), "SMS_T2bH_mSbottom450_mLSP1");
  // //  sig->setColor(kCyan-8);
  sig->setColor(kCyan+2);
  sig->setFullName("#tilde{b}#rightarrow bH#tilde{#chi_{0}};450,1");

  MT2Analysis<MT2EstimateTree>* sig_2 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_T2bH.root", ggDir.c_str() ), "SMS_T2bH_mSbottom450_mLSP300");
  sig_2->setColor(kBlue-4);
  sig_2->setFullName("#tilde{b}#rightarrow bH#tilde{#chi_{0}};450,300");

  MT2Analysis<MT2EstimateTree>* sig_9 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HZ.root", ggDir.c_str() ), "SMS_TChiHZ_HToGG_m175");
  sig_9->setColor(kMagenta+1);
  sig_9->setFullName("#tilde{#chi^{x}}#rightarrow H/Z#tilde{G};175");

  MT2Analysis<MT2EstimateTree>* sig_99 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_HH.root", ggDir.c_str() ), "SMS_TChiHH_HToGG_m175");
  sig_99->setColor(kCyan-4);
  sig_99->setFullName("#tilde{#chi^{x}}#rightarrow H#tilde{G};175");

  MT2Analysis<MT2EstimateTree>* sig_11 = MT2Analysis<MT2EstimateTree>::readFromFile(Form("%s/mc_WH.root", ggDir.c_str() ), "SMS_TChiWH_HToGG_m200_m1");
  sig_11->setColor(kPink-5);
  sig_11->setFullName("#tilde{#chi^{x}}#rightarrow H/W#tilde{#chi^{0}};200,1");


  std::string plotsDir = cfg.getEventYieldDir() + "/diPhotonControlRegion/plotsDataMC_sel";
  plotsDir = cfg.getEventYieldDir() + "/diPhotonControlRegion/plotsDataMC_sig_sel";

  if( shapeNorm ) plotsDir += "_shape";
  
  //This looks and works ok for the fractions but not the stack
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 

  bgYields.push_back( diPhoton );
  bgYields.push_back( gjets );
  bgYields.push_back( qcd );
  bgYields.push_back( higgs );

  bgYields.push_back( sig );
  bgYields.push_back( sig_2 );
  bgYields.push_back( sig_99 );
  bgYields.push_back( sig_9 );
  bgYields.push_back( sig_11 );

  MT2DrawTools dt(plotsDir, cfg.lumi() );
  
  dt.set_mc( &bgYields );
  dt.set_lumi( cfg.lumi() );

  //dt.set_addOverflow( 0 );
  
  std::string selection = "((( ptGamma0/h_mass) > 1./3) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 > 0.5 && r9Gamma1>0.5 && h_mass>=100 && h_mass<=180 ) && passPrefire" ;

  dt.set_data( 0 );
  dt.set_data2( 0 );

  dt.drawRegionYields_fromTree( "hgg_mt2"   , "hgg_mt2"   , selection, 20, 0., 100., "M_{T2}", "GeV", " ", " ");
  //  dt.drawRegionYields_fromTree( "pt_higgs"   , "h_pt"   , selection, 50, 0., 200., "p_{T}(H)", "GeV", " ", " ");
  dt.drawRegionYields_fromTree( "ptOM"   , "h_pt/h_mass"   , selection, 20, 0., 2., "p_{T}(H)/M_{H}", "", " ", " ");
  dt.drawRegionYields_fromTree(  "nJets" , "gg_nJets" , selection, 10, -0.5, 9.5, "Number of Jets (p_{T} > 30 GeV)", ""," ", " "); 
dt.drawRegionYields_fromTree( "nBJets", "gg_nBJets", selection, 7, -0.5, 6.5, "Number of b-Jets (p_{T} > 20 GeV)", "", " ", " ");



  //std::vector<TCanvas*> MT2DrawTools::drawRegionYields_fromTree( const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units, const std::string& kinCuts, const std::string& topoCuts ) {
  //  regionSelection_ph.push_back(" (gg_nJets>=1 && gg_nJets<4 && gg_nBJets == 0 && ((h_pt/h_mass) <  0.6)) " );

  //  std::vector<std::string> sel_name; 

  





  return 0;

}
























std::string getCutLabel( float theMin, float theMax, const std::string& name, const std::string& units ) {

  std::string cutLabel;
  if( theMax>theMin ) cutLabel = std::string(Form("%.0f < %s < %.0f %s", theMin, name.c_str(), theMax, units.c_str()) );
  else                cutLabel = std::string(Form("%s > %.0f %s", name.c_str(), theMin, units.c_str()) );

  return cutLabel;

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
    colors.push_back(430); // other=zll
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    //    colors.push_back(419); // z+jets
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
	//tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%s/puWeight", selection.c_str()) );
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%s", selection.c_str()) );
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
      if( shapeNorm ) {
        histos_mc[index]->Scale( scaleFactor );
      }
      else
	histos_mc[index]->Scale( 145.0/106.5 );

      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack.Add(histos_mc[index]);
    }


    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();
    TPad *pad1 = MT2DrawTools::getCanvasMainPad();
    TPad *pad2 = MT2DrawTools::getCanvasRatioPad();
        
    TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600);
    c1_log->cd();
    TPad *pad1_log = MT2DrawTools::getCanvasMainPad( true );
    TPad *pad2_log = MT2DrawTools::getCanvasRatioPad( true );
 
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

    std::string binWidthText;
    if( binWidth>=1. )         binWidthText = (std::string)Form("%.0f", binWidth);
    else if( binWidth>=0.1 )   binWidthText = (std::string)Form("%.1f", binWidth);
    else if( binWidth>=0.01 )  binWidthText = (std::string)Form("%.2f", binWidth);
    else if( binWidth>=0.001 ) binWidthText = (std::string)Form("%.3f", binWidth);
    else                       binWidthText = (std::string)Form("%.4f", binWidth);

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

    for( unsigned i=0; i<niceNames.size(); ++i ) {
      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.04);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      regionText->AddText( niceNames[i].c_str() );

      pad1->cd();
      // regionText->Draw("same");
      pad1_log->cd();
      //  regionText->Draw("same");
    }
    
    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      normText->AddText( "#splitline{Shape}{Norm.}" );
      pad1->cd();
      //normText->Draw("same");
      pad1_log->cd();
      //normText->Draw("same");
    }

    TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( gr_data, "Data", "P" );
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], bgYields[i]->getFullName().c_str(), "F" );
    }


    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    
    TPaveText* ratioText = new TPaveText( 0.133, -0.051, 0.4, 0.1 , "brNDC" );
    ratioText->SetTextSize(0.04);
    ratioText->SetTextFont(40);
    ratioText->SetTextColor(2);
    ratioText->SetFillColor(0);
    ratioText->SetTextAlign(11);
    ratioText->AddText( Form("Data/MC = %.2f", scaleFactor) );
    //  ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", scaleFactor, error_datamc) );
     

    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    
    TLine* lineSF = new TLine(xMin, scaleFactor, xMax, scaleFactor);
    lineSF->SetLineColor(2);

    float yMinR=0.0;
    float yMaxR=2.0;

    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
    TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
 
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);
   
    //    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );
    TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    TGraphErrors* SFFitBand = MT2DrawTools::getSFFitBand(fSF, xMin, xMax);
    TPaveText* fitText = MT2DrawTools::getFitText( fSF );


    c1->cd();
    pad1->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
      fitText->Draw("same");
    // ratioText->Draw("same");
  
    gPad->RedrawAxis();

    c1_log->cd();
    pad1_log->cd();
    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");
    labelTop->Draw("same");
    if( !shapeNorm )
     fitText->Draw("same");
    //  ratioText->Draw("same");

    gPad->RedrawAxis();

   /*
    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    TLine* lineSF = MT2DrawTools::getSFLine(integral_data, integral_mc, xMin, xMax);
    TGraphErrors* SFband = MT2DrawTools::getSFBand(integral_data, error_data, integral_mc, error_mc, xMin, xMax);
    */

    c1->cd();
    //   TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
    pad2->Draw();
    pad2->cd();

    h2_axes_ratio->Draw("");
 
    /*  line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same");
    */
    lineCentral->Draw("same");
    if( !shapeNorm ){

      systBand->Draw("3,same");
      lineCentral->Draw("same");

      SFFitBand->Draw("3,same");
      fSF->Draw("same");
    }

    g_ratio->Draw("PE,same");    
    gPad->RedrawAxis();


    c1_log->cd();
    // TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
    pad2_log->Draw();
    pad2_log->cd();

    h2_axes_ratio->Draw(""); 
    
    lineCentral->Draw("same");
    if( !shapeNorm ){

      systBand->Draw("3,same");
      lineCentral->Draw("same");

      SFFitBand->Draw("3,same");
      fSF->Draw("same");
    }
    /*
    line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same"); */
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

    delete h2_axes_ratio;

    delete h1_data;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

  }// for MT2 regions

}














































void drawShape( MT2Config cfg, MT2Analysis<MT2EstimateTree>* data, const std::string& saveName, const std::string& varName, const std::string& selection, const std::string& selection2, const std::string& selection3, unsigned int nSel, int nBins, float xMin, float xMax, std::string axisName, const std::string& units ) {

  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;

  std::vector<int> colors;

  colors.push_back(430); // other=zll
  colors.push_back(401); // qcd
  // colors.push_back(417); // w+jets
  // colors.push_back(419); // z+jets
  // colors.push_back(855); // top
  colors.push_back(97); // other


  std::string fullPathPlots = cfg.getEventYieldDir() + "/plotShapeComp";
  if( shapeNorm ) fullPathPlots += "_shape";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_mc = data->get(thisRegion)->tree;
    // TH1D* h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
    // tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );

    // TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    // gr_data->SetMarkerStyle(20);
    // gr_data->SetMarkerSize(1.2);

    float lumi = cfg.lumi();
    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<nSel; ++i ) { 
      std::string thisName = "h1_" + i;
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_mc->Sumw2();
      if(i==0){
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi, selection.c_str()) );
      }else if(i==1){
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi, selection2.c_str()) );
      }else if(i==2){
	tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi, selection3.c_str()) );
      }
      h1_mc->SetLineColor( colors[i] );
      h1_mc->SetLineWidth( 2 );

      h1_mc->Scale( 1./ h1_mc->Integral() );
      histos_mc.push_back(h1_mc);
    }

    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();
    // TPad *pad1 = MT2DrawTools::getCanvasMainPad();
    // TPad *pad2 = MT2DrawTools::getCanvasRatioPad(false);
        
    TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600);
    c1_log->cd();
    c1_log->SetLogy();
    // TPad *pad1_log = MT2DrawTools::getCanvasMainPad( true );
    // TPad *pad2_log = MT2DrawTools::getCanvasRatioPad( false );
 
    float yMaxScale = 1.1;
    float yMax1 = histos_mc[0]->GetMaximum()*yMaxScale;
    float yMax2 = histos_mc[1]->GetMaximum()*yMaxScale;
    float yMax3 = histos_mc[2]->GetMaximum()*yMaxScale;

    float yMax = yMax1;
    if( yMax2 > yMax1 ) yMax = yMax2;
    if( yMax3 > yMax ) yMax = yMax3;

    // float yMax2 = yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
    // float yMax3 = yMaxScale*(bgStack.GetMaximum());
    // float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    // if( yMax3 > yMax ) yMax = yMax3;
    // if( h1_data->GetNbinsX()<2 ) yMax *=3.;

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
    else                       binWidthText = (std::string)Form("%.4f", binWidth);

    std::string yAxisTitle;
    if( units!="" ) 
      yAxisTitle = (std::string)(Form("Events / (%s %s)", binWidthText.c_str(), units.c_str()));
    else
      yAxisTitle = (std::string)(Form("Events / (%s)", binWidthText.c_str()));

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
    c1->Draw();
    //   c1->cd();
    h2_axes->Draw();
   
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.001, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();
    c1_log->Draw();
    //   c1_log->cd();
    h2_axes_log->Draw();
   
    std::vector<std::string> niceNames = thisRegion.getNiceNames();

    for( unsigned i=0; i<niceNames.size(); ++i ) {
      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.04);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      regionText->AddText( niceNames[i].c_str() );

      c1->cd();
      // regionText->Draw("same");
      c1_log->cd();
      //  regionText->Draw("same");
    }
    
    if( shapeNorm ) {
      TPaveText* normText = new TPaveText( 0.45, 0.8, 0.68, 0.9, "brNDC" );
      normText->SetFillColor(0);
      normText->SetTextSize(0.035);
      normText->AddText( "#splitline{Shape}{Norm.}" );
      c1->cd();
      //normText->Draw("same");
      c1_log->cd();
      //normText->Draw("same");
    }

    TLegend* legend = new TLegend( 0.7, 0.9-(histos_mc.size())*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    //    legend->AddEntry( gr_data, "Data", "P" );
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], Form("Bin %d",i), "F" );
      //legend->AddEntry( histos_mc[i], bgYields[i]->getFullName().c_str(), "F" );
    }


    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    
    // TPaveText* ratioText = new TPaveText( 0.133, -0.051, 0.4, 0.1 , "brNDC" );
    // ratioText->SetTextSize(0.04);
    // ratioText->SetTextFont(40);
    // ratioText->SetTextColor(2);
    // ratioText->SetFillColor(0);
    // ratioText->SetTextAlign(11);
    // ratioText->AddText( Form("Data/MC = %.2f", scaleFactor) );
    // //  ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", scaleFactor, error_datamc) );
     

    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    
    // TLine* lineSF = new TLine(xMin, scaleFactor, xMax, scaleFactor);
    // lineSF->SetLineColor(2);

    // float yMinR=0.0;
    // float yMaxR=2.0;

    // TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
    // TGraphAsymmErrors* g_ratio = MT2DrawTools::getRatioGraph(h1_data, histo_mc);
 
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    //    TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr);
   
    //    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr );
    //    TF1* fSF = MT2DrawTools::getSFFit(g_ratio, xMin, xMax);
    // TGraphErrors* SFFitBand = MT2DrawTools::getSFFitBand(fSF, xMin, xMax);
    // TPaveText* fitText = MT2DrawTools::getFitText( fSF );


    c1->cd();
    c1->cd();
    legend->Draw("same");
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      histos_mc[i]->Draw("histo same");
    }
    // bgStack.Draw("histo same");
    //    gr_data->Draw("p same");
    labelTop->Draw("same");
    // if( !shapeNorm )
    //   fitText->Draw("same");
    // ratioText->Draw("same");
  
    gPad->RedrawAxis();

    c1_log->cd();
    c1_log->cd();
    legend->Draw("same");
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      histos_mc[i]->Draw("histo same");
    }
    // bgStack.Draw("histo same");
    // gr_data->Draw("p same");
    labelTop->Draw("same");
    // if( !shapeNorm )
    //  fitText->Draw("same");
    //  ratioText->Draw("same");

    gPad->RedrawAxis();

   /*
    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(1);
    TLine* lineSF = MT2DrawTools::getSFLine(integral_data, integral_mc, xMin, xMax);
    TGraphErrors* SFband = MT2DrawTools::getSFBand(integral_data, error_data, integral_mc, error_mc, xMin, xMax);
    */

    // c1->cd();
    // //   TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
    // pad2->Draw();
    // pad2->cd();

    //    h2_axes_ratio->Draw("");
 
    /*  line->Draw("same");
    SFband->Draw("3,same");
    lineSF->Draw("same");
    */
    // lineCentral->Draw("same");
    // if( !shapeNorm ){

    //   systBand->Draw("3,same");
    //   lineCentral->Draw("same");

    //   SFFitBand->Draw("3,same");
    //   fSF->Draw("same");
    //}

    //    g_ratio->Draw("PE,same");    
  // gPad->RedrawAxis();


    // c1_log->cd();
    // // TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
    // pad2_log->Draw();
    // pad2_log->cd();

    // //    h2_axes_ratio->Draw(""); 
     
    // lineCentral->Draw("same");
    // if( !shapeNorm ){

    //   systBand->Draw("3,same");
    //   lineCentral->Draw("same");

    //   SFFitBand->Draw("3,same");
    //   fSF->Draw("same");
    // }
    // /*
    // line->Draw("same");
    // SFband->Draw("3,same");
    // lineSF->Draw("same"); */
    // //    g_ratio->Draw("PE,same");
    // gPad->RedrawAxis();


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

    // delete h2_axes_ratio;

    // delete h1_data;
  
    for( unsigned i=0; i<histos_mc.size(); ++i )
      delete histos_mc[i];

  }// for MT2 regions

}
