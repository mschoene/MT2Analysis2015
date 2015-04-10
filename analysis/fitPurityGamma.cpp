#include <cmath>
#include <iostream>

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"


#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateSyst.h"
#include "../interface/MT2EstimateZinvGamma.h"

using namespace RooFit;



struct Purity {

  float purity;
  float purityErrUp;
  float purityErrDown;

};



void fitSinglePurity( const std::string& outputdir, Purity& loose, Purity& tight, RooRealVar* x, RooDataSet* data, TH1D* h1_templPrompt, TH1D* h1_templFake );
void fitPurity( const std::string& outputdir, MT2EstimateSyst* purityLoose, MT2EstimateSyst* purityTight, RooRealVar* x, std::vector<RooDataSet*> data, TH1D* templPrompt, TH1D* templFake );
void checkBoundaries( Purity& p );





int main( int argc, char* argv[] ) {


  std::string samples = "PHYS14_v4_skimprune";


  std::string mc_or_data = "MC";
  if( argc>1 ) {
    mc_or_data = std::string(argv[1]);
    if( mc_or_data=="data" ) mc_or_data="DataRC";
    if( mc_or_data=="mc" ) mc_or_data="MC";
  }



  std::string regionsSet = "zurich";
  if( argc>2 ) {
    regionsSet = std::string(argv[2]);
  }



  TH1::AddDirectory(kFALSE);


  std::string gammaCRdir = "GammaControlRegion_" + samples + "_" + regionsSet;
  MT2Analysis<MT2EstimateZinvGamma>* gammaJet_data = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( gammaCRdir + "/data.root", "gammaCR_loose" );

  MT2Analysis<MT2EstimateZinvGamma>* templates_prompt = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( "gammaTemplates" + mc_or_data + "_" + samples + "_13TeV_inclusive.root", "templatesPrompt" );
  MT2Analysis<MT2EstimateZinvGamma>* templates_fake   = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( "gammaTemplates" + mc_or_data + "_" + samples + "_13TeV_inclusive.root", "templatesFake" );


  std::string outputdir = gammaCRdir + "/PurityFits" + mc_or_data;
  system( Form( "mkdir -p %s/singleFits", outputdir.c_str()) );

  std::set<MT2Region> regions = gammaJet_data->getRegions();

  MT2Analysis<MT2EstimateSyst>* purityLoose = new MT2Analysis<MT2EstimateSyst>( "purityLoose", regionsSet );
  MT2Analysis<MT2EstimateSyst>* purityTight = new MT2Analysis<MT2EstimateSyst>( "purity"     , regionsSet );


  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()>1 ) continue;

    MT2EstimateZinvGamma* thisEstimate = gammaJet_data->get( *iR );

    MT2EstimateZinvGamma* templatePrompt = templates_prompt->get( *(templates_prompt->matchRegion( *iR )) );
    MT2EstimateZinvGamma* templateFake   = templates_fake  ->get( *(templates_fake  ->matchRegion( *iR )) );

    MT2EstimateSyst* thisLoosePurity = purityLoose->get( *iR );
    std::string nameLoose = thisLoosePurity->yield->GetName();

    MT2EstimateSyst* thisTightPurity = purityTight->get( *iR );
    std::string nameTight = thisTightPurity->yield->GetName();

    fitPurity( outputdir, thisLoosePurity, thisTightPurity, thisEstimate->x_, thisEstimate->iso_bins, templatePrompt->iso, templateFake->iso);

    thisLoosePurity->yield->SetName( nameLoose.c_str() );
    thisTightPurity->yield->SetName( nameTight.c_str() );

  }
    

  purityLoose->writeToFile( outputdir + "/purityFit.root" );
  purityTight->addToFile( outputdir + "/purityFit.root" );

  return 0;

}





void fitPurity( const std::string& outputdir, MT2EstimateSyst* purityLoose, MT2EstimateSyst* purityTight, RooRealVar* x, std::vector<RooDataSet*> data, TH1D* templPrompt, TH1D* templFake ) {


  for( unsigned i=0; i<data.size(); ++i ) {

    int ibin = i+1;

    Purity loose, tight;
    fitSinglePurity( outputdir, loose, tight, x, data[i], templPrompt, templFake );

    purityLoose->yield         ->SetBinContent( ibin, loose.purity );
    purityLoose->yield_systUp  ->SetBinContent( ibin, loose.purity + loose.purityErrUp );
    purityLoose->yield_systDown->SetBinContent( ibin, loose.purity - loose.purityErrDown );

    purityTight->yield         ->SetBinContent( ibin, tight.purity );
    purityTight->yield_systUp  ->SetBinContent( ibin, tight.purity + tight.purityErrUp );
    purityTight->yield_systDown->SetBinContent( ibin, tight.purity - tight.purityErrDown );

  }


  return;

}





void fitSinglePurity( const std::string& outputdir, Purity& loose, Purity& tight, RooRealVar* x, RooDataSet* data, TH1D* h1_templPrompt, TH1D* h1_templFake ) {


  float dataIntegral = data->sumEntries();

  if( dataIntegral == 0. ) {
    loose.purity=-1;
    loose.purityErrUp=0.;
    loose.purityErrDown=0.;
    tight.purity=-1;
    tight.purityErrUp=0.;
    tight.purityErrDown=0.;
    return;
  }


  RooDataHist templPrompt("templPrompt", "", *x, h1_templPrompt);
  RooDataHist templFake  ("templFake"  , "", *x, h1_templFake  );

  RooHistPdf pdfPrompt("pdfPrompt", "", *x, templPrompt );
  RooHistPdf pdfFake  ("pdfFake"  , "", *x, templFake   );

  RooRealVar sigFrac("promptFrac","fraction of prompt",0.9,0.,1.) ;
  RooAddPdf  model("model","", RooArgList(pdfPrompt,pdfFake), sigFrac) ;

  int nBins = h1_templPrompt->GetNbinsX();
  float xMin = h1_templPrompt->GetXaxis()->GetXmin();
  float xMax = h1_templPrompt->GetXaxis()->GetXmax();

  float xMaxFit = 0.999*xMax;
  x->setRange( "fittingRange", 0., xMaxFit );
  model.fitTo(*data, SumW2Error(kTRUE), Minos(kTRUE), Range("fittingRange")); 

  loose.purity = sigFrac.getVal();
  loose.purityErrUp = sigFrac.getErrorHi();
  loose.purityErrDown = -sigFrac.getErrorLo();

  float thresh = 2.5;
  int cutBin = h1_templPrompt->FindBin(thresh) - 1;

  float sigEff = h1_templPrompt->Integral(1, cutBin)/h1_templPrompt->Integral(1,nBins);
  float bgEff = h1_templFake->Integral(1, cutBin)/h1_templFake->Integral(1,nBins);
  float sigPassCut = sigFrac.getVal()*sigEff;
  float bgPassCut = (1.-sigFrac.getVal())*bgEff;
  tight.purity = sigPassCut / (sigPassCut+bgPassCut);
  tight.purityErrUp = loose.purityErrUp; // is it ok to assign the same error also to the tight purity?
  tight.purityErrDown = loose.purityErrDown;
  //float factor = tight.purity/loose.purity;
  //tight.purityErrUp = loose.purityErrUp*factor;
  //tight.purityErrDown = loose.purityErrDown*factor;

  checkBoundaries( loose );
  checkBoundaries( tight );

  RooPlot* xframe = x->frame();
  data->plotOn(xframe, Binning(nBins, xMin, xMax));
  model.plotOn(xframe);

  // Overlay the background component of model with a dashed line
  model.plotOn(xframe,Components(pdfFake),LineStyle(kDashed)) ;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  gPad->SetLeftMargin(0.15);
  TH2D* h2_axes = new TH2D("axes", "", 10, 0., xMaxFit, 10, 0., xframe->GetMaximum()*1.1 );
  h2_axes->Draw();
  xframe->GetYaxis()->SetTitleOffset(1.4); 
  xframe->Draw("same");

  c1->SaveAs(Form("%s/singleFits/purityFit_%s.eps", outputdir.c_str(), data->GetName()));
  c1->SaveAs(Form("%s/singleFits/purityFit_%s.png", outputdir.c_str(), data->GetName()));
  c1->SaveAs(Form("%s/singleFits/purityFit_%s.pdf", outputdir.c_str(), data->GetName()));

  delete c1;
  delete xframe;

  return;

}



void checkBoundaries( Purity& p ) {

  if( p.purity - p.purityErrDown < 0. ) p.purityErrDown = p.purity;
  if( p.purity + p.purityErrUp   > 1. ) p.purityErrUp   = 1. - p.purity;

}
