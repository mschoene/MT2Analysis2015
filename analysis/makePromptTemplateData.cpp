#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"



using namespace RooFit;


float lumi = 5.; //fb-1





MT2Analysis<MT2EstimateZinvGamma>* subtractFakes( const std::string& outputdir, MT2Analysis<MT2EstimateZinvGamma>* templatesPromptRaw, MT2Analysis<MT2EstimateZinvGamma>* templatesFake );
void removeNegatives( MT2Analysis<MT2EstimateZinvGamma>* data );
void removeNegativesSingleHisto( TH1D* h1 );




int main( int argc, char* argv[] ) {

  //if( argc==1 ) {
  //  std::cout << "-> You need to pass me the regions set name. Here are some suggestions: " << std::endl;
  //  std::cout << "  13TeV_CSA14" << std::endl;
  //  std::cout << "  13TeV_onlyHT" << std::endl;
  //  std::cout << "  13TeV_onlyJets" << std::endl;
  //  std::cout << "  13TeV_inclusive" << std::endl;
  //  exit(101);
  //}

  //std::string regionsSet = "13TeV_onlyHT";
  std::string regionsSet = "13TeV_inclusive";
  //std::string regionsSet = "13TeV_inclusive";
  //std::string regionsSet = "13TeV_ZinvGammaPurity";
  if( argc>1 ) {
    std::string regionsSet_tmp(argv[1]); 
    regionsSet = regionsSet_tmp;
  }


  std::string samplesName = "PHYS14_v2_Zinv";

  std::string fileName = "gammaTemplatesData_" + samplesName + "_" + regionsSet + ".root";
  MT2Analysis<MT2EstimateZinvGamma>* templatesPromptRaw = MT2Analysis<MT2EstimateZinvGamma>::readFromFile(fileName, "templatesPromptRaw");
  MT2Analysis<MT2EstimateZinvGamma>* templatesFake      = MT2Analysis<MT2EstimateZinvGamma>::readFromFile(fileName, "templatesFake");


  std::string outputdir = "FakeFitsData_" + samplesName + "_" + regionsSet;
  system( Form("mkdir -p %s", outputdir.c_str()) );

  MT2Analysis<MT2EstimateZinvGamma>* templatesPrompt= subtractFakes( outputdir, templatesPromptRaw, templatesFake );

  templatesPrompt->addToFile(fileName, true);


  return 0;

}






MT2Analysis<MT2EstimateZinvGamma>* subtractFakes( const std::string& outputdir, MT2Analysis<MT2EstimateZinvGamma>* templatesPromptRaw, MT2Analysis<MT2EstimateZinvGamma>* templatesFake ) {


  std::set<MT2Region> regions = templatesPromptRaw->getRegions();


  std::set<MT2EstimateZinvGamma*> promptTemplates;


  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateZinvGamma* thisFake = templatesFake->get( *iR );
    MT2EstimateZinvGamma* thisPromptRaw = templatesPromptRaw->get( *iR );

    TH1D* thisFakeTempl = thisFake->iso;
    TH1D* thisPromptRawTempl = thisPromptRaw->iso;

    int binMin = thisPromptRaw->iso->FindBin(8.);
    int binMax = thisPromptRaw->iso->FindBin(20.) -1;
    float intData = thisPromptRawTempl->Integral( binMin, binMax );
    float intMC   = thisFakeTempl     ->Integral( binMin, binMax );
    std::cout <<  binMin << std::endl;
    std::cout <<  binMax << std::endl;
    std::cout << intData << std::endl;
    std::cout << intMC   << std::endl;
    float sf = intData/intMC;

    MT2EstimateZinvGamma* thisPrompt = new MT2EstimateZinvGamma("templatesPrompt", *iR );
    std::string oldName(thisPrompt->iso->GetName());
    thisPrompt->iso = (TH1D*)(thisPromptRawTempl->Clone());
    thisPrompt->iso->Add(thisFakeTempl, -sf);
    thisPrompt->iso->SetName(oldName.c_str());

    promptTemplates.insert( thisPrompt );

  }


  MT2Analysis<MT2EstimateZinvGamma>* templates = new MT2Analysis<MT2EstimateZinvGamma>("templatesPrompt", promptTemplates);

  removeNegatives( templates );

  return templates;

}




void removeNegatives( MT2Analysis<MT2EstimateZinvGamma>* data ) {

  std::set<MT2Region> MT2Regions = data->getRegions();

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
      
    removeNegativesSingleHisto( data->get(thisRegion)->yield );
    removeNegativesSingleHisto( data->get(thisRegion)->iso );
    removeNegativesSingleHisto( data->get(thisRegion)->sietaieta );

  } // for regions

}


void removeNegativesSingleHisto( TH1D* h1 ) {

  for( unsigned ibin=1; ibin<h1->GetXaxis()->GetNbins()+1; ++ibin ) {

    if( h1->GetBinContent(ibin)<0. ) h1->SetBinContent( ibin, 0. );

  }

}

