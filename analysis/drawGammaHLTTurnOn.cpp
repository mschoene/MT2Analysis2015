#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "THStack.h"
#include "TEfficiency.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Sample.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"

#define mt2_cxx
#include "interface/mt2.h"



void drawTurnOn( const MT2Config& cfg, MT2Sample sample );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./drawGammaHLTTurnOn [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  MT2DrawTools::setStyle();

  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";
  std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "SinglePhoton"); 
  if( samples_data.size()==0 ) {
    std::cout << "There must be an error: samples_data is empty!" << std::endl;
    exit(1209);
  }


  drawTurnOn( cfg, samples_data[0] );

  return 0;

}




void drawTurnOn( const MT2Config& cfg, MT2Sample sample ) {


  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");

  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;


  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);


  int nBins = 10;
  float xMin = 120.;
  float xMax = 220.;

  TH1D* h1_num = new TH1D( "num", "", nBins, xMin, xMax );
  TH1D* h1_denom = new TH1D( "denom", "", nBins, xMin, xMax );

  TH1D* h1_num_eb = new TH1D( "num_eb", "", nBins, xMin, xMax );
  TH1D* h1_denom_eb = new TH1D( "denom_eb", "", nBins, xMin, xMax );


  int nentries = tree->GetEntries();


  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( myTree.ngamma==0 ) continue;
    //if( fabs(myTree.gamma_eta[0])>1.444 ) continue;
    if( myTree.gamma_idCutBased[0]==0 ) continue;

    float iso = myTree.gamma_chHadIso[0];
    if( iso>10. ) continue; // loose iso

    if( fabs(myTree.gamma_eta[0])<1.479 && myTree.gamma_sigmaIetaIeta[0]>0.0106 ) continue;
    if( fabs(myTree.gamma_eta[0])>1.479 && myTree.gamma_sigmaIetaIeta[0]>0.0266 ) continue;

    //TLorentzVector gamma;
    //gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );


    if( myTree.HLT_photon120 ) {

      h1_denom->Fill( myTree.gamma_pt[0] );
      if( myTree.HLT_photon165_HE10 ) h1_num->Fill( myTree.gamma_pt[0] );

      if( fabs(myTree.gamma_eta[0])<1.479 ) {
        h1_denom_eb->Fill( myTree.gamma_pt[0] );
        if( myTree.HLT_photon165_HE10 ) h1_num_eb->Fill( myTree.gamma_pt[0] );
      }

    }

  }  // for entries


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0., 1.1 );
  h2_axes->SetXTitle("Photon p_{T} [GeV]");
  h2_axes->SetYTitle("Efficiency of HLT_Photon165_HE10");

  h2_axes->Draw();

  TLine* lineOne = new TLine( xMin, 1., xMax, 1. );
  lineOne->Draw("same");

  TEfficiency* eff = new TEfficiency( *h1_num, *h1_denom );
  eff->SetMarkerStyle(24);
  eff->SetMarkerSize(1.6);
  eff->SetMarkerColor(46);
  eff->SetLineColor(46);

  TEfficiency* eff_eb = new TEfficiency( *h1_num_eb, *h1_denom_eb );
  eff_eb->SetMarkerStyle(20);
  eff_eb->SetMarkerSize(1.6);
  eff_eb->SetMarkerColor(46);
  eff_eb->SetLineColor(46);

  eff->Draw("p same");
  eff_eb->Draw("p same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop(0.0073);
  labelTop->Draw("same");

  TLegend* legend = new TLegend( 0.2, 0.65, 0.45, 0.8 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( eff, "All Photons", "P" );
  legend->AddEntry( eff_eb, "EB Only", "P" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/gammaControlRegion/hltTurnOn.eps", cfg.getEventYieldDir().c_str()) );
  c1->SaveAs( Form("%s/gammaControlRegion/hltTurnOn.pdf", cfg.getEventYieldDir().c_str()) );

}
