#include <iostream>

#include "TCanvas.h"
#include "TH1D.h"

#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2DrawTools.h"




int main( int argc, char* argv[] ) {


  MT2DrawTools::setStyle();


  std::string samples = "PHYS14_v4_skimprune";
  if( argc>1 ) {
    std::string samples_tmp(argv[1]); 
    samples = samples_tmp;
  }


  std::string regionsSet = "zurich";
  float lumi = 4.;

  
  std::string gammaControlRegionDir(Form("GammaControlRegion_%s_%s_%.0ffb", samples.c_str(), regionsSet.c_str(), lumi));

  MT2Analysis<MT2Estimate>* purityData = MT2Analysis<MT2Estimate>::readFromFile( gammaControlRegionDir + "/PurityFitsDataRC/purityFit.root", "purity" );
  MT2Analysis<MT2Estimate>* purityMC   = MT2Analysis<MT2Estimate>::readFromFile( gammaControlRegionDir + "/PurityFitsMC/purityFit.root", "purity" );

  if( purityData==0 || purityMC==0 ) {
    std::cout << "ERROR! You need to compute the purities for data and MC first!" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(1);
  }

  std::set<MT2Region> regions = purityData->getRegions();
  TH1D* h1_syst = new TH1D("syst", "", 100., -10., 10.);
  h1_syst->SetXTitle("(Data-MC)/Data [%]");
  h1_syst->SetYTitle("Entries");


  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()>1 ) continue;

    MT2Estimate* thisPurityData = purityData->get(*iR);
    MT2Estimate* thisPurityMC   = purityMC  ->get(*iR);

    int nbins = thisPurityData->yield->GetNbinsX();

    for( int ibin=1; ibin<nbins+1; ++ibin ) {

      float pData = thisPurityData->yield->GetBinContent(ibin);
      float pMC   = thisPurityMC  ->yield->GetBinContent(ibin);

      float relDiff = (pData-pMC)/pData;

      h1_syst->Fill( 100.*relDiff );

    } // for mt2 bins

  } // for regions


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1_syst->Draw();

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  c1->SaveAs("prova.eps");

  c1->Clear();
  c1->SetLogy();
  
  h1_syst->Draw();
  labelTop->Draw("same");

  c1->SaveAs("prova_log.eps");


  return 0;

}
