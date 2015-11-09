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


#include "../interface/MT2DrawTools.h"
#include "../interface/MT2Config.h"
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




void fitSinglePurity( const MT2Config& cfg, Purity& loose, Purity& tight, RooRealVar* x, RooDataSet* data, TH1D* h1_templPrompt, TH1D* h1_templFake );
void fitPurity( const MT2Config& cfg, MT2EstimateSyst* purityLoose, MT2EstimateSyst* purityTight, RooRealVar* x, std::vector<RooDataSet*> data, TH1D* templPrompt, TH1D* templFake );
void checkBoundaries( Purity& p );
void makePurity( const MT2Config& cfg, std::string outputdir, MT2Analysis<MT2EstimateZinvGamma>* data, MT2Analysis<MT2EstimateZinvGamma>* temp_prompt,MT2Analysis<MT2EstimateZinvGamma>* temp_fake, bool useMC,  const std::string var="" );


int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|              Running fitPurityGamma                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./fitPurityGamma [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool useMC = true;

  if( argc>2 ) {
    std::string data_or_mc = std::string(argv[2]); 
    if( data_or_mc=="data" || data_or_mc=="Data" || data_or_mc=="DATA" ) useMC=false;
    else if( data_or_mc!="mc" && data_or_mc!="MC" ) {
      std::cout << std::endl;
      std::cout << "-> WARNING! Second argument should be 'data' or 'MC'." << std::endl;
      std::cout << "Exiting." << std::endl;
      std::cout << std::endl;
      exit(817);
    }
  } 

  std::string templateType = cfg.gammaTemplateType();
  if( argc>3 ) {
    templateType = std::string(argv[3]); 
    std::cout << std::endl;
    std::cout << "-> Will disobey the cfg and use templateType = " << argv[2] << std::endl;
    std::cout << std::endl;
  } 
  cfg.set_gammaTemplateType(templateType);


  MT2DrawTools::setStyle();
  TH1::AddDirectory(kFALSE);


  std::string gammaCRdir = cfg.getEventYieldDir() + "/gammaControlRegion"; 
  MT2Analysis<MT2EstimateZinvGamma>* gammaJet_data = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( gammaCRdir + "/data.root", "gammaCR_loose" );
  
 MT2Analysis<MT2EstimateZinvGamma>* gammaJet_data_ht = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( gammaCRdir + "/iso_ht.root", "iso_ht" );

MT2Analysis<MT2EstimateZinvGamma>* gammaJet_data_njets = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( gammaCRdir + "/iso_nJets.root", "iso_njets" );

MT2Analysis<MT2EstimateZinvGamma>* gammaJet_data_nbjets = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( gammaCRdir + "/iso_nBJets.root", "iso_nbjets" );

  std::string templateFileName = gammaCRdir + "/gammaTemplates" + cfg.gammaTemplateType();
  if( useMC ) templateFileName = templateFileName + "_MC";
  else        templateFileName = templateFileName + "_data";
  templateFileName = templateFileName + ".root";


  std::string outputdir = gammaCRdir + "/PurityFits" + cfg.gammaTemplateType();
  system( Form( "mkdir -p %s/singleFits", outputdir.c_str()) );


  MT2Analysis<MT2EstimateZinvGamma>* templates_prompt = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( templateFileName, "templatesPrompt" );
  MT2Analysis<MT2EstimateZinvGamma>* templates_fake   = MT2Analysis<MT2EstimateZinvGamma>::readFromFile( templateFileName, "templatesFake" );


  makePurity( cfg, outputdir,  gammaJet_data, templates_prompt, templates_fake, useMC, "");
  makePurity( cfg, outputdir,  gammaJet_data_ht, templates_prompt, templates_fake, useMC, "ht_");
  makePurity( cfg, outputdir,  gammaJet_data_njets, templates_prompt, templates_fake, useMC, "njets_");
  makePurity( cfg, outputdir,  gammaJet_data_nbjets, templates_prompt, templates_fake, useMC, "nbjets_");


  return 0;

}
















void makePurity( const MT2Config& cfg, std::string outputdir,MT2Analysis<MT2EstimateZinvGamma>* data, MT2Analysis<MT2EstimateZinvGamma>* temp_prompt,MT2Analysis<MT2EstimateZinvGamma>* temp_fake, bool useMC,  const std::string var ){


  MT2EstimateZinvGamma* templatePrompt, *templateFake;
   if( cfg.gammaTemplateRegions()=="13TeV_inclusive" ) { // just get them once
    templatePrompt = temp_prompt->get( MT2Region("HT200toInf_j1toInf_b0toInf") );
    templateFake   = temp_fake  ->get( MT2Region("HT200toInf_j1toInf_b0toInf") );
  }

  MT2Analysis<MT2EstimateSyst>* purityLoose = new MT2Analysis<MT2EstimateSyst>( "purityLoose", cfg.regionsSet() );
  MT2Analysis<MT2EstimateSyst>* purityTight = new MT2Analysis<MT2EstimateSyst>( "purity"     , cfg.regionsSet() );


  std::set<MT2Region> regions = data->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()>2 ) continue;

    MT2EstimateZinvGamma* thisEstimate = data->get( *iR );

    if( cfg.gammaTemplateRegions()!="13TeV_inclusive" ) {
      templatePrompt = temp_prompt->get( *(temp_prompt->matchRegion( *iR )) );
      templateFake   = temp_fake  ->get( *(temp_fake  ->matchRegion( *iR )) );
    }

  
    MT2EstimateSyst* thisLoosePurity = purityLoose->get( *iR );
    std::string nameLoose = thisLoosePurity->yield->GetName();

    MT2EstimateSyst* thisTightPurity = purityTight->get( *iR );
    std::string nameTight = thisTightPurity->yield->GetName();

    int nBins;
    double* bins;
    thisEstimate->MT2Estimate::getYieldBins(nBins, bins);
    thisLoosePurity->rebinYields( (MT2Analysis<MT2Estimate>*)purityLoose,nBins,bins);
    thisTightPurity->rebinYields( (MT2Analysis<MT2Estimate>*)purityTight,nBins,bins);

    fitPurity( cfg, thisLoosePurity, thisTightPurity, thisEstimate->x_, thisEstimate->iso_bins, templatePrompt->iso, templateFake->iso);

    thisLoosePurity->yield->SetName( nameLoose.c_str() );
    thisTightPurity->yield->SetName( nameTight.c_str() );
  }

  if(useMC) {
    purityLoose->writeToFile( outputdir + "/purityFit_"+var +"MC.root" );
    purityTight->addToFile( outputdir + "/purityFit_"+var +"MC.root" );
 } else {
    purityLoose->writeToFile( outputdir + "/purityFit_"+var +"data.root" );
    purityTight->addToFile( outputdir + "/purityFit_"+var +"data.root" );
  }

  return;
}







void fitPurity( const MT2Config& cfg, MT2EstimateSyst* purityLoose, MT2EstimateSyst* purityTight, RooRealVar* x, std::vector<RooDataSet*> data, TH1D* templPrompt, TH1D* templFake ) {

  for( unsigned i=0; i<data.size(); ++i ) {

    int ibin = i+1;

    Purity loose, tight;
    fitSinglePurity( cfg, loose, tight, x, data[i], templPrompt, templFake );

    if( loose.purity>=0. ){
    
      purityLoose->yield         ->SetBinContent( ibin, loose.purity );
      purityLoose->yield_systUp  ->SetBinContent( ibin, loose.purity + loose.purityErrUp );
      purityLoose->yield_systDown->SetBinContent( ibin, loose.purity - loose.purityErrDown );
      
    }
    else{

      purityLoose->yield         ->SetBinContent( ibin, 1.0 );
      purityLoose->yield_systUp  ->SetBinContent( ibin, 1.0 );
      purityLoose->yield_systDown->SetBinContent( ibin, 0.0 );

    }

    if( tight.purity>=0 ){
      
      purityTight->yield         ->SetBinContent( ibin, tight.purity );
      purityTight->yield_systUp  ->SetBinContent( ibin, tight.purity + tight.purityErrUp );
      purityTight->yield_systDown->SetBinContent( ibin, tight.purity - tight.purityErrDown );
    
    }
    else{
      
      purityTight->yield         ->SetBinContent( ibin, 1.0 );
      purityTight->yield_systUp  ->SetBinContent( ibin, 1.0 );
      purityTight->yield_systDown->SetBinContent( ibin, 0.0 );  

    }

  }

  return;
}





void fitSinglePurity( const MT2Config& cfg, Purity& loose, Purity& tight, RooRealVar* x, RooDataSet* data, TH1D* h1_templPrompt, TH1D* h1_templFake ) {


  float dataIntegral = data->sumEntries();
  float thresh = cfg.gammaIsoCut();
  //float data_pass = data->sumEntries(Form("x<%f", thresh));

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

  RooRealVar sigFrac ("promptFrac" ,"fraction of prompt",0.9,0.,1.);
  RooAddPdf  model ("model" ,"", RooArgList(pdfPrompt,pdfFake), sigFrac);

  int nBins = h1_templPrompt->GetNbinsX();
  float xMin = h1_templPrompt->GetXaxis()->GetXmin();
  float xMax = h1_templPrompt->GetXaxis()->GetXmax();

  float xMaxFit = 0.9999*xMax;
  x->setRange( "fittingRange", 0., xMaxFit );
  model .fitTo(*data, SumW2Error(kTRUE), Minos(kTRUE), Range("fittingRange")); 

  loose.purity = sigFrac.getVal();
  loose.purityErrUp = sigFrac.getErrorHi();
  loose.purityErrDown = -sigFrac.getErrorLo();

  float sig = sigFrac.getVal()*dataIntegral;
  float bg  = (1.-sigFrac.getVal())*dataIntegral;

  int cutBin = h1_templPrompt->FindBin(thresh) - 1;
  float sigEff = h1_templPrompt->Integral(1, cutBin)/h1_templPrompt->Integral(1,nBins);
  float bgEff  = h1_templFake  ->Integral(1, cutBin)/h1_templFake  ->Integral(1,nBins);

  float sig_pass = sigEff*sig;
  float bg_pass  = bgEff*bg;

  float purity_tight = sig_pass/(sig_pass + bg_pass);

  tight.purity        = purity_tight;
  tight.purityErrUp   = loose.purityErrUp   * tight.purity / loose.purity;
  tight.purityErrDown = loose.purityErrDown * tight.purity / loose.purity;


  //float sigPassCut = sigFrac.getVal()*sigEff;
  //float bgEff = h1_templFake->Integral(1, cutBin)/h1_templFake->Integral(1,nBins);
  //float bgPassCut = (1.-sigFrac.getVal())*bgEff;
  //tight.purity = sigPassCut / (sigPassCut+bgPassCut);
  //tight.purityErrUp = loose.purityErrUp; // is it ok to assign the same error also to the tight purity?
  //tight.purityErrDown = loose.purityErrDown;
  ////float factor = tight.purity/loose.purity;
  ////tight.purityErrUp = loose.purityErrUp*factor;
  ////tight.purityErrDown = loose.purityErrDown*factor;

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
  h2_axes->SetXTitle("Photon Charged Isolation [GeV]");
  h2_axes->SetYTitle("Events");
  h2_axes->Draw();
  xframe->GetYaxis()->SetTitleOffset(1.4); 
  xframe->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
  labelTop->Draw("same");

  std::string outputdir = cfg.getEventYieldDir() + "/gammaControlRegion/PurityFits" + cfg.gammaTemplateType();
  c1->SaveAs(Form("%s/singleFits/purityFit_%s.eps", outputdir.c_str(), data->GetName()));
  c1->SaveAs(Form("%s/singleFits/purityFit_%s.png", outputdir.c_str(), data->GetName()));
  c1->SaveAs(Form("%s/singleFits/purityFit_%s.pdf", outputdir.c_str(), data->GetName()));

  delete c1;
  delete xframe;

  return;

}



void checkBoundaries( Purity& p ) {

  if( p.purity > 1. ) p.purity = 1.;
  if( p.purity < 0. ) p.purity = 0.;
  if( p.purity - p.purityErrDown < 0. ) p.purityErrDown = p.purity;
  if( p.purity + p.purityErrUp   > 1. ) p.purityErrUp   = 1. - p.purity;

}
