#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2DrawTools.h"

#define mt2_cxx
#include "interface/mt2.h"

float lumi = 5.; //fb-1

MT2Analysis<MT2EstimateSyst> computeYield( const MT2Sample& sample, const std::string& regionsSet, float lumi=1. );

int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./computeLostLepton [samplesFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  
  std::string sampleName(argv[1]);
  
  std::string samplesFileName = "../samples/samples_" + sampleName + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName);
  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }

  std::string outputdir = "EventYields_" + sampleName;
  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string regionsSet = "13TeV_CSA14"; 

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  
  MT2Analysis<MT2EstimateSyst>* lostLeptonEstimate = new MT2Analysis<MT2EstimateSyst> ( "LostLepton", regionsSet );  
  for( unsigned i=0; i < fSamples.size(); ++i )
    (*lostLeptonEstimate) += ( computeYield( fSamples[i], regionsSet, lumi ) );
  
  //lostLeptonEstimate->writeToFile("MT2LostLeptonEstimate_v1.root");
  lostLeptonEstimate->writeToFile("TTJets_CR.root");
  //lostLeptonEstimate->writeToFile("MT2QCDEstimate_v1.root.root");

  return 0;
  
}




MT2Analysis<MT2EstimateSyst> computeYield( const MT2Sample& sample, const std::string& regionsSet, float lumi ) {

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  MT2Tree myTree;
  myTree.Init(tree);

  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  MT2Analysis<MT2EstimateSyst> analysis( sample.sname, regionsSet, sample.id );

  int nentries = tree->GetEntries();
    
  //ofstream ofs("events.log");

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    ////if( myTree.nTaus20>0 ) continue;
    //if( myTree.nMuons10>0 ) continue;
    //if( myTree.nElectrons10>0 ) continue;
    ////// ^^ the above should be changed to the following
    if( myTree.nMuons10 > 0) continue;
    if( myTree.nElectrons10 > 0 ) continue;
    if( myTree.nPFLep5LowMT > 0) continue;
    if( myTree.nPFHad10LowMT > 0) continue;

    if( myTree.nVert==0 ) continue;
    if( myTree.nJet40<2 ) continue;
    if( myTree.jet_pt[1]<100. ) continue;
    if( myTree.deltaPhiMin<0.3 ) continue;
    if( myTree.diffMetMht>0.5*myTree.met_pt ) continue;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = myTree.mt2;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet40;
    int nbjets = myTree.nBJet20; 

    ////////To remove for next iteration
    //int nPFLep5LowMT = 0, nPFHad10LowMT = 0;
    //
    //for( int l=0; l < myTree.nisoTrack; ++l ){
    //  
    //  if( sqrt(2*(myTree.isoTrack_pt[l])*met*( 1 + TMath::Cos( myTree.isoTrack_phi[l] - myTree.met_phi ) ) ) > 100. )
    //    continue;
    //  
    //  if ( ( fabs(myTree.isoTrack_pdgId[l]) == 11 || fabs(myTree.isoTrack_pdgId[l]) == 13 ) && myTree.isoTrack_pt[l] > 5 && myTree.isoTrack_absIso[l]/myTree.isoTrack_pt[l] < 0.2)
    //	++nPFLep5LowMT;
    //  
    //  else if ( fabs(myTree.isoTrack_pdgId[l]) == 211 && myTree.isoTrack_pt[l] > 10 && myTree.isoTrack_absIso[l]/myTree.isoTrack_pt[l] < 0.1)
    //	++nPFHad10LowMT;
    //   
    //}
    //
    //if( nPFLep5LowMT > 0 || nPFHad10LowMT > 0 ) continue;
    ////////
    //
    //////// minMTBmet - also to remove at next iteration
    //float minMTBmet = 999999.;
    //for (int j=0; j< myTree.njet; ++j){
    //  
    //  if(myTree.jet_btagCSV[j] < 0.679) continue;
    //  if(myTree.jet_pt[j] < 40 || fabs(myTree.jet_eta[j]) > 2.5) continue;
    //  
    //  float thisMTBmet = sqrt(2*(myTree.jet_pt[j])*met*( 1 + TMath::Cos( myTree.jet_phi[j] - myTree.met_phi ) ) );
    //  if( thisMTBmet < minMTBmet ) minMTBmet = thisMTBmet;
    // 
    //}
    ////////

    Double_t weight = myTree.evt_scale1fb*lumi;

    float fullweight_btagUp = weight;
    float fullweight_btagDown = weight;

    MT2EstimateSyst* thisEstimate = analysis.get( ht, njets, nbjets, met, minMTBmet, mt2 );
    if( thisEstimate==0 ) continue;

    thisEstimate->yield         ->Fill(mt2, weight );
    thisEstimate->yield_btagUp  ->Fill(mt2, fullweight_btagUp );
    thisEstimate->yield_btagDown->Fill(mt2, fullweight_btagDown );

    //ofs << "entry " << iEntry <<  "\tmet " << met << "\tmt2 " << mt2 << "\tminMTBmet " << minMTBmet << std::endl;
    
  } // for entries

  //ofs.close();

  analysis.finalize();

  delete tree;

  file->Close();
  delete file;

  return analysis;

}
