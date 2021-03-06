#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateSyst.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"



bool closureTest = false;

bool onlyUseUpToRunG = true; // only use up to run G for rphi, fj, rb
float lumiRatioGtoH = 27.70/36.46;
int lastRunG = 280385;

//float prescales[3] = {7000., 180.0, 60.0}; // 2015
float prescales[3] = {7900., 440.6, 110.2}; // up to run G 27.7 fb-1
//float prescales[3] = {9200., 460.6, 115.2}; // full 2016


void compareFractions( const MT2Config& cfg, const std::string& outputdir, const std::string& dataFile, const std::string& analysisName, const std::string& xaxisName, const std::string& yaxisName, bool logPlot, const std::string& postfix="" );
void drawClosure( const std::string& outputdir, MT2Analysis<MT2Estimate>* estimate, MT2Analysis<MT2Estimate>* mcTruth, float scaleEst, float lumi, MT2Analysis<MT2Estimate>* nonQCD=NULL );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./drawQCDFromDeltaPhi [configFileName] [closureTest=true]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  MT2DrawTools::setStyle();


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);



  if( argc>2 ) {
    std::string closurestr(argv[2]);
    if( closurestr=="closureTest" || closurestr=="true" ) {
      closureTest=true;
      std::cout << "-> Running closure test in validation region" << std::endl;
    }
  } 


  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion/";
  std::string qcdESTdir = qcdCRdir;
  if( closureTest ) qcdESTdir += "/test";
  std::string outputdir = qcdCRdir;
  if( closureTest ) outputdir += "/test";

  std::string mcFile   = qcdESTdir + "/qcdEstimateMC.root";
  std::string dataFile = qcdESTdir + "/qcdEstimateData.root";


  compareFractions( cfg, outputdir, dataFile, "f_jets", "N_{j}"  , "f_{j}", false           );
  compareFractions( cfg, outputdir, dataFile, "f_jets", "N_{j}"  , "f_{j}", false , "_noPS" );
  compareFractions( cfg, outputdir, dataFile, "r_hat" , "N_{b}"  , "r_{b}", true            );


  MT2Analysis<MT2EstimateTree>* qcdTree_mc   = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/mc.root",   "qcdCRtree" );
  MT2Analysis<MT2EstimateTree>* qcdTree_sr   = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/mc.root",   "qcdCRtree" );
  //  MT2Analysis<MT2EstimateTree>* qcdTree_sr   = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/../analyses.root",   "QCD" );
  MT2Analysis<MT2EstimateTree>* qcdTree_data = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/data.root", "qcdCRtree" );

  std::string runRange = "";
  if (onlyUseUpToRunG)
    runRange = Form("&&(run<=%d||ht>1000)",lastRunG);

  MT2Analysis<MT2EstimateTree>* mcTruth;
  MT2Analysis<MT2EstimateTree>* data  ;
  MT2Analysis<MT2EstimateTree>* nonQCD;
  if( closureTest ) { // in validation region we also use id==15
    //    mcTruth = MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "mcTruth", cfg.regionsSet(), qcdTree_mc  , 1, 100, 200 , "(id>150&&(id>151||ht<450)&&(id>152||ht<575)&&(id>153||ht<1000)&&(id>154||ht<1500)) && id<200 && mt2>100 && mt2<200. && deltaPhiMin>0.3"); // signal region for mcTruth
    mcTruth = MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "mcTruth", cfg.regionsSet(), qcdTree_mc  , 1, 100, 200 , "id>150 && id<200 && mt2>100 && mt2<200. && deltaPhiMin>0.3"); // signal region for mcTruth    
    data    = MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "data"   , cfg.regionsSet(), qcdTree_data, 1, 100, 200 , "((id&1)==1"+runRange+") && mt2>100. && mt2<200. && deltaPhiMin>0.3"         ); // signal region for data
    //// nonQCD  = MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "nonQCD" , cfg.regionsSet(), qcdTree_mc  , 1, 100, 200 , "(id>=300 || id<=0) && mt2>100. && mt2<200. && deltaPhiMin>0.3"         ); // signal region for nonQCD mcTruth
    nonQCD  = MT2EstimateTree::makeRebinnedAnalysisFromInclusiveTree( "nonQCD" , cfg.regionsSet(), qcdTree_mc  , 1, 100, 200 , "id>=300 && mt2>100. && mt2<200. && deltaPhiMin>0.3"         ); // signal region for nonQCD mcTruth

  }
  else
    //mcTruth = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mcTruth", cfg.regionsSet(), qcdTree_mc  , "(id>150&&(id>152||ht<450)&&(id>153||ht<575)&&(id>154||ht<1000)&&(id>155||ht<1500)) && id<200 && mt2>200. && deltaPhiMin>0.3 && (met>250 || ht>1000)" ); // signal region for mcTruth
    mcTruth = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mcTruth", cfg.regionsSet(), qcdTree_mc  , "id> 10 && id<200 && mt2>200. && deltaPhiMin>0.3 && (met>250 || ht>1000)" ); // signal region for mcTruth  
    ////mcTruth = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mcTruth", cfg.regionsSet(), qcdTree_mc  , "(id>150&&(id>151||ht<450)&&(id>152||ht<575)&&(id>153||ht<1000)&&(id>154||ht<1500)) && id<200 && mt2>200. && deltaPhiMin>0.3 && (met>250 || ht>1000)" ); // signal region for mcTruth

  MT2Analysis<MT2Estimate>* estimateMC     = MT2Analysis<MT2Estimate>::readFromFile( mcFile  , "qcdEstimate" );
  MT2Analysis<MT2Estimate>* estimateData   = MT2Analysis<MT2Estimate>::readFromFile( dataFile, "qcdEstimate" );

  mcTruth     ->setColor(kQCD  );
  estimateMC  ->setColor(kBlack);
  estimateData->setColor(kBlack);

  std::string plotsDirMC = qcdESTdir + "/plotsMC";
  std::string plotsDirData = qcdESTdir + "/plotsData";
  if ( closureTest ) {
    estimateData->setColor(kQCD       );
    data        ->setColor(kBlack     );
    nonQCD      ->setColor(kLostLepton);
    drawClosure( plotsDirMC  , estimateMC  , (MT2Analysis<MT2Estimate>*)mcTruth, cfg.lumi(), cfg.lumi() );
    drawClosure( plotsDirData, estimateData, (MT2Analysis<MT2Estimate>*)data   ,   1.0     , cfg.lumi(), (MT2Analysis<MT2Estimate>*) nonQCD);
  } 
  else {
    drawClosure( plotsDirMC  , estimateMC  , (MT2Analysis<MT2Estimate>*)mcTruth, cfg.lumi(), cfg.lumi() );
    drawClosure( plotsDirData, estimateData, (MT2Analysis<MT2Estimate>*)mcTruth,   1.0     , cfg.lumi() );
  }


  return 0;

} 

           

 


void compareFractions( const MT2Config& cfg, const std::string& outputdir, const std::string& dataFile, const std::string& analysisName, const std::string& xaxisName, const std::string& yaxisName, bool logPlot, const std::string& postfix ) {

  MT2Analysis<MT2Estimate>* fraction_mc   = MT2Analysis<MT2Estimate>::readFromFile(dataFile, analysisName+"_mc");
  MT2Analysis<MT2Estimate>* fraction_data = MT2Analysis<MT2Estimate>::readFromFile(dataFile, analysisName+"_data"+postfix.c_str());


  std::set<MT2Region> regions = fraction_data->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* thisFraction_mc   = fraction_mc  ->get(*iR);
    MT2Estimate* thisFraction_data = fraction_data->get(*iR);

    TH1D* h1_mc   = thisFraction_mc  ->yield;
    TH1D* h1_data = thisFraction_data->yield;


    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
    if( logPlot ) c1->SetLogy();

    float yMin = (logPlot) ? h1_mc->GetMinimum()/7. : 0.;
    float yMax = 1.2;
    TH2D* h2_axes = new TH2D( "axes", "", 10, h1_mc->GetXaxis()->GetXmin(), h1_mc->GetXaxis()->GetXmax(), 10, yMin, yMax );
    h2_axes->GetXaxis()->SetTitleSize(0.06);
    h2_axes->SetXTitle( xaxisName.c_str() );
    h2_axes->SetYTitle( "" );
    h2_axes->Draw();
    TLatex *ytitle = new TLatex(0.02,0.92,yaxisName.c_str());
    ytitle->SetNDC();
    ytitle->SetTextFont(42);
    ytitle->SetTextSize(0.065);
    if (yaxisName!="f_{j}"){
      ytitle->Draw();
    }
    else{
      ytitle->SetTextFont(12);
      ytitle->DrawLatexNDC(0.026,0.92, "f");
      ytitle->SetTextFont(42);
      ytitle->DrawLatexNDC(0.04,0.92, "_{j}");
    }

    h1_mc->SetLineColor(kRed);
    h1_mc->SetLineWidth(2);

    TH1D* h1_mcBand = MT2DrawTools::getMCBandHisto( h1_mc );
    h1_mcBand->SetLineColor(kRed);
    h1_mcBand->SetLineWidth(2);

    h1_data->SetMarkerStyle(20);
    h1_data->SetMarkerSize(1.6);
    h1_data->SetLineWidth(2);

    h1_mcBand->Draw("e2 same");
    h1_mc  ->Draw("hist norm same");
    h1_data->Draw("p same norm");


    TPaveText* labelTop = MT2DrawTools::getLabelTop( cfg.lumi()*(onlyUseUpToRunG ? lumiRatioGtoH : 1.0) );
    labelTop->Draw("same");
    // //TPaveText* labelCMS = MT2DrawTools::getLabelCMS();
    // //TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Unpublished");
    // TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Supplementary");
    // labelCMS->Draw("same");
    TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Preliminary");
    labelCMS->Draw("same");

    // TLatex *arxiv = new TLatex(); arxiv->SetNDC();
    // arxiv->SetTextSize(0.035); arxiv->SetTextFont(42);
    // if (analysisName=="f_jets")
    //   arxiv->DrawLatex(0.18,0.88,"arXiv:1603.04053");
    // else
    //   arxiv->DrawLatex(0.18,0.2,"arXiv:1603.04053");

    // // TPaveText* labelTop = MT2DrawTools::getLabelTop( "CMS Preliminary, 2.2 fb^{-1} at #sqrt{s} = 13 TeV" );
    // // labelTop->Draw("same");

    std::vector<std::string> regionNiceNames = iR->getNiceNames();
    std::string niceJetName;
    if( iR->nJetsMax()<0 ) 
      niceJetName = std::string( Form( "N_{j} #geq %d", iR->nJetsMin() ) );
    else 
      niceJetName = std::string( Form( "%d #leq N_{j} #leq %d", iR->nJetsMin(), iR->nJetsMax() ) );

    std::string legendTitle = (analysisName=="f_jets") ? regionNiceNames[0].c_str() : niceJetName;
    TLegend* legend = new TLegend( 0.55, 0.72, 0.9, 0.92, legendTitle.c_str() );
    legend->SetFillColor(0);
    legend->SetTextSize(0.038);
    legend->AddEntry( h1_data, "Data", "PL");
    legend->AddEntry( h1_mcBand, "Simulation", "FL");
    legend->Draw("same");

    gPad->RedrawAxis();

    c1->SaveAs( Form("%s/%s_%s%s.eps", outputdir.c_str(), analysisName.c_str(), iR->getName().c_str(), postfix.c_str()) );
    c1->SaveAs( Form("%s/%s_%s%s.pdf", outputdir.c_str(), analysisName.c_str(), iR->getName().c_str(), postfix.c_str()) );
    c1->SaveAs( Form("%s/%s_%s%s.png", outputdir.c_str(), analysisName.c_str(), iR->getName().c_str(), postfix.c_str()) );
    c1->SaveAs( Form("%s/%s_%s%s.C"  , outputdir.c_str(), analysisName.c_str(), iR->getName().c_str(), postfix.c_str()) );

    delete c1;
    delete h2_axes;
  
  }


}





void drawClosure( const std::string& outputdir, MT2Analysis<MT2Estimate>* estimate, MT2Analysis<MT2Estimate>* mcTruth , float scaleEst, float lumi, MT2Analysis<MT2Estimate>* nonQCD) {

  bool doClosureTestData = false;
  if ( nonQCD != NULL ) doClosureTestData = true;

  system(Form("mkdir -p %s/pdf/" , outputdir.c_str()));
  system(Form("mkdir -p %s/png/" , outputdir.c_str()));
  system(Form("mkdir -p %s/eps/" , outputdir.c_str()));
  system(Form("mkdir -p %s/C/"   , outputdir.c_str()));
  system(Form("mkdir -p %s/root/", outputdir.c_str()));

  
  std::set<MT2Region> MT2Regions = estimate->getRegions();
  
  TH1D* h_estimate_tot = new TH1D("h_estimate_tot", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  h_estimate_tot->Sumw2();
  h_estimate_tot->GetYaxis()->SetTitle("Events");
  h_estimate_tot->SetMarkerStyle(20);
  h_estimate_tot->SetMarkerSize(1.6);
  h_estimate_tot->SetLineColor( estimate->getColor() );
  h_estimate_tot->SetMarkerColor( estimate->getColor() );
  
  TH1D* h_nonQCD_tot = new TH1D("h_nonQCD_tot", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  h_nonQCD_tot->Sumw2();
  h_nonQCD_tot->GetYaxis()->SetTitle("Events");
  
  if ( doClosureTestData ){
    h_estimate_tot->SetMarkerStyle(1);
    h_estimate_tot->SetFillColor  (estimate->getColor());
    h_nonQCD_tot  ->SetLineColor  (nonQCD  ->getColor());
    h_nonQCD_tot  ->SetFillColor  (nonQCD  ->getColor());
    h_nonQCD_tot  ->SetMarkerSize (0);
  }
  THStack* stack_tot = new THStack("stack_tot","");

  TH1D* h_mcTruth_tot = new TH1D("h_mcTruth_tot", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  h_mcTruth_tot->Sumw2();
  h_mcTruth_tot->GetYaxis()->SetTitle("Events");
  h_mcTruth_tot->SetFillColor(0);
  h_mcTruth_tot->SetLineColor( mcTruth->getColor() );
  h_mcTruth_tot->SetMarkerColor( mcTruth->getColor() );
  h_mcTruth_tot->SetMarkerStyle(20);
  h_mcTruth_tot->SetMarkerSize(1.6);
  
  TH1D* hPull = new TH1D("hPull", "", 20, -5, 5);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle(!doClosureTestData ? "(Data Driven - MC)/#sigma" : "(Data - Estimate)/#sigma");
  hPull->GetYaxis()->SetTitle("Events");
  
  TH1D* hPull_int = new TH1D("hPull_int", "", 20, -5, 5);
  hPull_int->Sumw2();
  hPull_int->GetXaxis()->SetTitle(!doClosureTestData ? "(Data Driven - MC)/#sigma" : "(Data - Estimate)/#sigma");
  hPull_int->GetYaxis()->SetTitle("Events");
  
  
  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {


      estimate->get(*iMT2)->yield->Scale(scaleEst);
      if ( !doClosureTestData ) // mcTruth is data in this case
	mcTruth ->get(*iMT2)->yield->Scale(lumi);

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      TH1D* h_estimate = estimate->get(*iMT2)->yield;
      h_estimate->SetMarkerStyle(20);
      h_estimate->SetMarkerSize(1.6);
      h_estimate->SetLineColor( estimate->getColor() );
      h_estimate->SetMarkerColor( estimate->getColor() );
//      if(iMT2->htMin()>=1400) h_estimate->SetBinContent(1, 0.0); 
//      if(iMT2->htMin()>=1400) h_estimate->SetBinError(1, 0.0); 

      THStack* stack = new THStack("stack","");

      TH1D* h_nonQCD;
      if ( doClosureTestData ){ 
	float ps = 1./(onlyUseUpToRunG ? lumiRatioGtoH : 1.);
	if     ( iMT2->htMin() < 300. ) ps = prescales[0];
	else if( iMT2->htMin() < 500. ) ps = prescales[1];
	else if( iMT2->htMin() < 600. ) ps = prescales[2];

	// apply scale factors
	float sf = 1.;  // 2015 numbers from lost lepton from Dominick. Confirmed in zll from Myriam
	// if     ( iMT2->htMin() < 500.  ) sf = 0.97; // prescale
	// else if( iMT2->htMin() < 1000. ) sf = 0.85;
	// else if( iMT2->htMin() < 1500. ) sf = 0.82;
	// else                             sf = 0.57;

	h_nonQCD = nonQCD->get(*iMT2)->yield;
	h_nonQCD->Scale(lumi/ps*sf);
	//h_estimate->Add(nonQCD->get(*iMT2)->yield); // add non-QCD mc to estimate. todo: make stack

	h_estimate->SetMarkerStyle(1);
	h_estimate->SetFillColor  (estimate->getColor());
	h_nonQCD  ->SetLineColor  (nonQCD  ->getColor());
	h_nonQCD  ->SetFillColor  (nonQCD  ->getColor());
	h_nonQCD  ->SetMarkerSize (0);
	stack->Add(h_nonQCD);
	stack->Add(h_estimate);
      }

      int nBins = h_estimate->GetXaxis()->GetNbins();
      if ( closureTest )  nBins -= 1;  // remove overflow for validation in 100<mt2<200 (for both mc and data)

      double err_estimate;
      double int_estimate = h_estimate->IntegralAndError(1, nBins+1, err_estimate);
      //      if(iMT2->htMin()>=1400.) h_estimate->IntegralAndError(2, nBins+1, err_estimate);
      
      h_estimate_tot->SetBinContent(iRegion, int_estimate);
      h_estimate_tot->SetBinError  (iRegion, err_estimate);
      
      h_estimate_tot->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );

      double int_nonQCD=0., err_nonQCD=0.;
      if ( doClosureTestData ) {
	int_nonQCD = h_nonQCD->IntegralAndError(1, nBins+1, err_nonQCD);

	// add 50% uncertainty on MC
	err_nonQCD += 0.5 * int_nonQCD; 
 
	h_nonQCD_tot->SetBinContent(iRegion, int_nonQCD);
	h_nonQCD_tot->SetBinError  (iRegion, err_nonQCD);

	h_nonQCD_tot->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
      }

      
      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.15);
      pad1->Draw();
      pad1->cd();

      TH1D* h_mcTruth = mcTruth->get(*iMT2)->yield;

      h_mcTruth->SetLineColor( mcTruth->getColor() );
      h_mcTruth->SetMarkerColor( mcTruth->getColor() );
      h_mcTruth->SetMarkerStyle(20);
      h_mcTruth->SetMarkerSize(1.6);
      
//      if( iMT2->htMin()>1400 )
//	h_mcTruth->SetBinContent(1, 0);
      
      double err_int;
      double int_mcTruth = h_mcTruth->IntegralAndError(1, nBins+1, err_int);
      h_mcTruth_tot->SetBinContent(iRegion, h_mcTruth->Integral());
      h_mcTruth_tot->SetBinError(iRegion, err_int);

//      if(iMT2->htMin()>1400){
//	
//	int_mcTruth = h_mcTruth->IntegralAndError(2, nBins+1, err_int);
//	h_mcTruth_tot->SetBinContent(iRegion, int_mcTruth);
//	h_mcTruth_tot->SetBinError(iRegion, err_int);
//	
//      }
      
      h_mcTruth_tot->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );

      if ( doClosureTestData ) {
	int_estimate += int_nonQCD;
	err_estimate = sqrt(err_estimate*err_estimate + err_nonQCD*err_nonQCD);
      }
      if(int_mcTruth>0)
        hPull_int->Fill((int_estimate-int_mcTruth)/sqrt(err_estimate*err_estimate+err_int*err_int));

      for (int iBin=1; iBin<=h_mcTruth->GetNbinsX(); iBin++){
	double vEst = (!doClosureTestData) ? h_estimate->GetBinContent(iBin) : ((TH1F*)stack->GetStack()->Last())->GetBinContent(iBin);
	double eEst = (!doClosureTestData) ? h_estimate->GetBinError  (iBin) : ((TH1F*)stack->GetStack()->Last())->GetBinContent(iBin);
	double vPre = h_mcTruth ->GetBinContent(iBin);
	double ePre = h_mcTruth ->GetBinError  (iBin);
	
	if ( vPre>0 && eEst+ePre != 0 )
	  hPull->Fill( ( doClosureTestData ? (vPre-vEst) : (vEst-vPre) )/sqrt(eEst*eEst + ePre*ePre) );
      }

      float xMin = h_estimate->GetXaxis()->GetXmin();
      float xMax = h_estimate->GetXaxis()->GetXmax();
      float yMax_1 = h_estimate->GetMaximum()*1.5;
      float yMax_2 = 1.2*(h_estimate->GetMaximum() + h_estimate->GetBinError(h_estimate->GetMaximumBin()));
      float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
      float yMax_3 = h_mcTruth->GetMaximum()*1.5;
      float yMax_4 = 1.2*(h_mcTruth->GetMaximum() + h_mcTruth->GetBinError(h_mcTruth->GetMaximumBin()));
      float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
      float yMax = (yMax1>yMax2) ? yMax1 : yMax2;

      TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
      h2_axes->SetXTitle("M_{T2} [GeV]");
      h2_axes->SetYTitle("Events");

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


      TLegend* legend = new TLegend( 0.6, 0.9-2.*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
      if ( doClosureTestData ) {
	legend->AddEntry( h_nonQCD  , " non-QCD"   , "F" );
	legend->AddEntry( h_estimate, "data-driven", "F" );
	legend->AddEntry( h_mcTruth , "data"       , "PL" );
      }
      else{
	legend->AddEntry( h_estimate, "data-driven", "PL" );
	legend->AddEntry( h_mcTruth , "MC QCD"     , "PL" );
      }

      legend->Draw("same");

      if  (doClosureTestData )
	stack->Draw("histe1 same");
      else
	h_estimate->Draw("Pe same");
      h_mcTruth->Draw("Pe same");


      TPaveText* labelTop;
      if ( scaleEst==1.0 )
	labelTop = MT2DrawTools::getLabelTop( lumi );
      else
	labelTop = MT2DrawTools::getLabelTopSimulation();
      labelTop->Draw("same");
      TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Preliminary");
      labelCMS->Draw("same");

      //TPaveText* labelCMS = MT2DrawTools::getLabelCMS();
      //labelCMS->Draw("same");

      gPad->RedrawAxis();

      c1->cd();
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.1);
      pad2->Draw();
      pad2->cd();

      std::string thisName = Form("%s_ratio", h_estimate->GetName());
      TH1D* h_ratio;
      if ( doClosureTestData ) {
	h_ratio = (TH1D*) ( h_mcTruth->Clone(thisName.c_str()) );
	h_ratio->Divide((TH1D*)stack->GetStack()->Last());
      }
      else {
	h_ratio = (TH1D*) ( h_estimate->Clone(thisName.c_str()) );
	h_ratio->Divide(h_mcTruth);
      }
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
      //h_ratio->Draw("PE");
      
      TH1D* h_band = MT2DrawTools::getBandAtOne( doClosureTestData ? ((TH1D*)stack->GetStack()->Last()) : h_mcTruth );

      TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, 0.0, 2.0 );
      
      TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
      lineCentral->SetLineColor(1);
      

      h2_axes_ratio->Draw(""       );
      h_band       ->Draw("e2same" );
      lineCentral  ->Draw("same"   );
      h_ratio      ->Draw("pe,same");

      gPad->RedrawAxis();
      
      c1->cd();

      TString filename = Form("closure%s_%s_%s", closureTest==false ? "SR" : "VR", iMT2->getName().c_str(), scaleEst==1.0 ? "data" : "mc");
      c1->SaveAs( Form("%s/eps/%s.eps"  , outputdir.c_str(), filename.Data()) );
      c1->SaveAs( Form("%s/pdf/%s.pdf"  , outputdir.c_str(), filename.Data()) );
      c1->SaveAs( Form("%s/png/%s.png"  , outputdir.c_str(), filename.Data()) );
      c1->SaveAs( Form("%s/C/%s.C"      , outputdir.c_str(), filename.Data()) );
      c1->SaveAs( Form("%s/root/%s.root", outputdir.c_str(), filename.Data()) );

      delete c1;
      delete h2_axes;
      delete h2_axes_ratio;
      //delete h_ratio;
      //delete h_estimate;
      //delete h_mcTruth;
      
      ++iRegion;

  } // for MT2 regions


  if ( doClosureTestData ) {
    stack_tot->Add( h_nonQCD_tot   );
    stack_tot->Add( h_estimate_tot );
  }
  
  TCanvas* c2 = new TCanvas("c2", "", 1100, 600);
  c2->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  pad1->SetLogy();
  
  float yMax_1 = h_estimate_tot->GetMaximum();
  float yMax_2 = h_estimate_tot->GetMaximum() + h_estimate_tot->GetBinError(h_mcTruth_tot->GetMaximumBin());
  float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
  float yMax_3 = h_mcTruth_tot->GetMaximum();
  float yMax_4 = h_mcTruth_tot->GetMaximum() + h_mcTruth_tot->GetBinError(h_mcTruth_tot->GetMaximumBin());
  float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
  float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
  yMax*=20.;
  
  float yMin = 1e-1;
  for( int iBin=1; iBin<h_estimate_tot->GetXaxis()->GetNbins()+1; ++iBin ) {
    if( h_estimate_tot    ->GetBinContent(iBin)>0. && h_estimate_tot    ->GetBinContent(iBin)<yMin ) yMin = h_estimate_tot    ->GetBinContent(iBin);
    if( h_mcTruth_tot->GetBinContent(iBin)>0. && h_mcTruth_tot->GetBinContent(iBin)<yMin ) yMin = h_mcTruth_tot->GetBinContent(iBin);
  }
  yMin /= 3.;
  yMin = TMath::Max(yMin, (float)1e-3); // i don't care about anything below 1e-3 and it only makes it ugly
  
  h_mcTruth_tot->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  h_mcTruth_tot->GetYaxis()->SetRangeUser(yMin, yMax);
  h_mcTruth_tot->GetXaxis()->LabelsOption("v");

  h_mcTruth_tot->GetYaxis()->SetTitleOffset(0.95);   h_mcTruth_tot->GetYaxis()->SetLabelSize(0.042); 

  h_mcTruth_tot->Draw("PE");


  if ( doClosureTestData ) {
    stack_tot ->Draw("histe1 same");
    h_mcTruth_tot->Draw("PEsame");
  }
  else
    h_estimate_tot->Draw( "pe,same" );

  TLegend* legend = new TLegend( 0.163, 0.67, 0.3, 0.79 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  if ( doClosureTestData ) {
    legend->AddEntry( h_nonQCD_tot  , "Non-multijet simulation" , "F" );
    //legend->AddEntry( h_estimate_tot, "Data-driven", "F" );
    legend->AddEntry( h_estimate_tot, "Multijet prediction", "F" );
    legend->AddEntry( h_mcTruth_tot , "Data"       , "PL" );
  }
  else{
    legend->AddEntry( h_estimate_tot, "Data-driven", "PL" );
    legend->AddEntry( h_mcTruth_tot , "Multijet simulation"     , "PL" );
  }

  legend->Draw("same");

  
  TPaveText* labelTop;
  if ( scaleEst==1.0 )
    labelTop = MT2DrawTools::getLabelTop( lumi );
  else
    labelTop = MT2DrawTools::getLabelTopSimulation();
  labelTop->Draw("same");
  // //TPaveText* labelCMS = MT2DrawTools::getLabelCMS();
  // //TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Unpublished");
  // TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Supplementary");
  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Preliminary");
  labelCMS->Draw("same");

  // TLatex *arxiv = new TLatex(); arxiv->SetNDC();
  // arxiv->SetTextSize(0.035); arxiv->SetTextFont(42);
  // arxiv->DrawLatex(0.17,0.81,"arXiv:1603.04053");
  
  TLine* lHT[4];
  for( int iHT=1; iHT < 5; iHT++ ){
    lHT[iHT-1] = new TLine(7*iHT+4*(iHT-1), -3, 7*iHT+4*(iHT-1), yMax );
    lHT[iHT-1]->SetLineColor(kBlack);
    lHT[iHT-1]->SetLineStyle(3);
    lHT[iHT-1]->SetLineWidth(2);

    lHT[iHT-1]->Draw("same");
  }


  int nHTRegions = 5;
  std::vector< std::string > htRegions;
  htRegions.push_back("H_{T} [250,450]");
  htRegions.push_back("H_{T} [450,575]");
  htRegions.push_back("H_{T} [575,1000]");
  htRegions.push_back("H_{T} [1000,1500]");
  htRegions.push_back("H_{T} > 1500 GeV");

  TPaveText* htBox[nHTRegions];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){

    if (iHT==0) htBox[iHT] = new TPaveText(0.16, 0.9-0.06+0.02, 0.26, 0.85+0.02, "brNDC");
    else htBox[iHT] = new TPaveText(0.26+0.17*(iHT-1), 0.9-0.06+0.02, 0.26+0.17+0.17*(iHT-1), 0.85+0.02, "brNDC");
    htBox[iHT]->AddText( htRegions[iHT].c_str() );

    htBox[iHT]->SetBorderSize(0);
    htBox[iHT]->SetFillColor(kWhite);
    htBox[iHT]->SetTextSize(0.035);
    htBox[iHT]->SetTextAlign(21); // align centered                                                                                                                                                                        
    htBox[iHT]->SetTextFont(62);
    htBox[iHT]->Draw("same");

  }

//  htRegions.push_back("H_{T} [200,450] GeV");
//  htRegions.push_back("H_{T} [450,575] GeV");
//  htRegions.push_back("H_{T} [575,100] GeV");
//  htRegions.push_back("H_{T} [1000,1500] GeV");
//  htRegions.push_back("H_{T} > 1500 GeV");
//  
//  TPaveText* htBox[nHTRegions];
//  for( int iHT = 0; iHT < nHTRegions; ++iHT){
//    
//    htBox[iHT] = new TPaveText(0.17+0.16*iHT, 0.9-0.06, 0.3+0.16*iHT, 0.9, "brNDC");
//
//    htBox[iHT]->AddText( htRegions[iHT].c_str() );
//    
//    htBox[iHT]->SetBorderSize(0);
//    htBox[iHT]->SetFillColor(kWhite);
//    htBox[iHT]->SetTextSize(0.038);
//    htBox[iHT]->SetTextAlign(21); // align centered
//    htBox[iHT]->SetTextFont(62);
//    htBox[iHT]->Draw("same");
//
//  }

  gPad->RedrawAxis();
  
  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();

  std::string thisName = Form("%s_ratio", h_estimate_tot->GetName());
  TH1D* h_Ratio = (TH1D*) h_estimate_tot->Clone(thisName.c_str());
  for( int iBin=1; iBin<h_Ratio->GetXaxis()->GetNbins()+1; ++iBin ) {
    float mc = doClosureTestData ? ((TH1F*)stack_tot->GetStack()->Last())->GetBinContent(iBin) : h_estimate_tot->GetBinContent(iBin);
    float mc_err = doClosureTestData ?  ((TH1F*)stack_tot->GetStack()->Last())->GetBinError(iBin) : h_estimate_tot->GetBinError(iBin);
    float est = h_mcTruth_tot->GetBinContent(iBin);
    float est_err = h_mcTruth_tot->GetBinError(iBin);
    float denom = sqrt( mc_err*mc_err + est_err*est_err );
    if( denom!=0. && mc>0. ) {
      h_Ratio->SetBinContent( iBin, (mc-est)/denom );
      h_Ratio->SetBinError( iBin, 1. );
    }
  }
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->SetLineWidth(2);
  

  TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( 0, MT2Regions.size(), -3., 3.);
  h2_axes_ratio->SetYTitle("Pull");

  TLine* LineCentral = new TLine(0, 0, MT2Regions.size(), 0);
  LineCentral->SetLineColor(1);

  h2_axes_ratio->Draw("");
  LineCentral->Draw("same");
  h_Ratio->Draw("pe,same");

  for( int iHT=1; iHT < 5; iHT++ ){
    lHT[iHT-1]->Draw("same");
  }


  gPad->RedrawAxis();

  c2->cd();
  TString filename = Form("closure%s_allRegions_pull_%s", closureTest==false ? "SR" : "VR", scaleEst==1.0 ? "data" : "mc");
  c2->SaveAs( Form("%s/eps/%s.eps"  , outputdir.c_str(), filename.Data()) );
  c2->SaveAs( Form("%s/pdf/%s.pdf"  , outputdir.c_str(), filename.Data()) );
  c2->SaveAs( Form("%s/png/%s.png"  , outputdir.c_str(), filename.Data()) );
  c2->SaveAs( Form("%s/C/%s.C"      , outputdir.c_str(), filename.Data()) );
  c2->SaveAs( Form("%s/root/%s.root", outputdir.c_str(), filename.Data()) );

  pad2->cd();
  pad2->Clear();

  float yMax_ratio = 2.0;
  float yMin_ratio = 0.;

  delete h2_axes_ratio;
  h2_axes_ratio = MT2DrawTools::getRatioAxes( 0, MT2Regions.size(),  yMin_ratio, yMax_ratio);
  h2_axes_ratio->SetYTitle(!doClosureTestData ? "Pred. / MC" : "Data / Pred.");
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->Draw("");

  TH1D *h_Band = MT2DrawTools::getBandAtOne( doClosureTestData ? ((TH1D*)stack_tot->GetStack()->Last()) : h_mcTruth_tot );

  h_Band->Draw("e2same");

  TLine* lineOne = new TLine(0, 1., MT2Regions.size(), 1.);
  lineOne->SetLineColor(1);
  lineOne->Draw("same");

  delete h_Ratio;
  if ( doClosureTestData ){
    h_Ratio = (TH1D*) ( h_mcTruth_tot->Clone(thisName.c_str()) );
    h_Ratio->Divide( (TH1D*)stack_tot->GetStack()->Last() );
  }
  else {
    h_Ratio = (TH1D*) ( h_estimate_tot->Clone(thisName.c_str()) );
    h_Ratio->Divide( h_mcTruth_tot );
  }
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->SetLineWidth(2);

  h_Ratio->Draw("pe,same");


  for( int iHT=1; iHT < 5; iHT++ ){
    lHT[iHT-1]->Draw("same");
  }

  gPad->RedrawAxis();

  filename = Form("closure%s_allRegions_ratio_%s", closureTest==false ? "SR" : "VR", scaleEst==1.0 ? "data" : "mc");
  c2->SaveAs( Form("%s/eps/%s.eps"  , outputdir.c_str(), filename.Data()) );
  c2->SaveAs( Form("%s/pdf/%s.pdf"  , outputdir.c_str(), filename.Data()) );
  c2->SaveAs( Form("%s/png/%s.png"  , outputdir.c_str(), filename.Data()) );
  c2->SaveAs( Form("%s/C/%s.C"      , outputdir.c_str(), filename.Data()) );
  c2->SaveAs( Form("%s/root/%s.root", outputdir.c_str(), filename.Data()) );


  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  hPull_int->SetStats(1110);
  TF1* f1_gaus = new TF1("f1_pull", "gaus", -2., 2.);
  f1_gaus->SetLineColor(kRed);
  hPull_int->Fit( f1_gaus, "QRL" );
  TPaveText* fitPars = new TPaveText( 0.2, 0.7, 0.5, 0.9, "brNDC" );
  fitPars->SetTextSize(0.03);
  fitPars->SetTextAlign(11);
  fitPars->SetFillColor(0);
  fitPars->AddText("Gaussian Fit:");
  fitPars->AddText(Form("Mean : %.2f +/- %.2f", f1_gaus->GetParameter(1), f1_gaus->GetParError(1) ));
  fitPars->AddText(Form("Sigma: %.2f +/- %.2f", f1_gaus->GetParameter(2), f1_gaus->GetParError(2) ));
  fitPars->Draw("same");
  hPull_int->Draw("hist same");
  f1_gaus->Draw("l same");

  filename = Form("closure%s_pull_int_%s", closureTest==false ? "SR" : "VR", scaleEst==1.0 ? "data" : "mc");
  c3->SaveAs( Form("%s/eps/%s.eps"  , outputdir.c_str(), filename.Data()) );
  c3->SaveAs( Form("%s/pdf/%s.pdf"  , outputdir.c_str(), filename.Data()) );
  c3->SaveAs( Form("%s/png/%s.png"  , outputdir.c_str(), filename.Data()) );
  c3->SaveAs( Form("%s/C/%s.C"      , outputdir.c_str(), filename.Data()) );
  c3->SaveAs( Form("%s/root/%s.root", outputdir.c_str(), filename.Data()) );

  TCanvas* c4 = new TCanvas("c4", "", 600, 600);
  c4->cd();
  hPull->SetStats(1110);
  f1_gaus->SetLineColor(kRed);
  hPull->Fit( f1_gaus, "QRL" );
  TPaveText* fitPars2 = new TPaveText( 0.2, 0.7, 0.5, 0.9, "brNDC" );
  fitPars2->SetTextSize(0.03);
  fitPars2->SetTextAlign(11);
  fitPars2->SetFillColor(0);
  fitPars2->AddText("Gaussian Fit:");
  fitPars2->AddText(Form("Mean : %.2f +/- %.2f", f1_gaus->GetParameter(1), f1_gaus->GetParError(1) ));
  fitPars2->AddText(Form("Sigma: %.2f +/- %.2f", f1_gaus->GetParameter(2), f1_gaus->GetParError(2) ));
  fitPars2->Draw("same");
  hPull->Draw("hist same");
  f1_gaus->Draw("l same");

  filename = Form("closure%s_pull_%s", closureTest==false ? "SR" : "VR", scaleEst==1.0 ? "data" : "mc");
  c4->SaveAs( Form("%s/eps/%s.eps"  , outputdir.c_str(), filename.Data()) );
  c4->SaveAs( Form("%s/pdf/%s.pdf"  , outputdir.c_str(), filename.Data()) );
  c4->SaveAs( Form("%s/png/%s.png"  , outputdir.c_str(), filename.Data()) );
  c4->SaveAs( Form("%s/C/%s.C"      , outputdir.c_str(), filename.Data()) );
  c4->SaveAs( Form("%s/root/%s.root", outputdir.c_str(), filename.Data()) );


  delete c2;
  delete c3;
  delete c4;
  delete h2_axes_ratio;

  delete h_estimate_tot;
  delete h_mcTruth_tot;
  delete stack_tot;
  delete hPull_int;
  
}


