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
#include "TVector2.h"

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2DrawTools.h"

#define mt2_cxx
#include "interface/mt2.h"

float lumi = 4.; //fb-1

struct lepcand {
  
  float pt;
  float eta;
  float phi;
  int pdgId;
  float mt;
  bool isPFCand;

};

MT2Analysis<MT2EstimateSyst> computeYield( const MT2Sample& sample, const std::string& regionsSet, float lumi=1. );

float DeltaR(float eta1, float eta2, float phi1, float phi2);

float DeltaPhi(float phi1, float phi2);

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

  //std::string outputdir = "EventYields_" + sampleName;
  //system(Form("mkdir -p %s", outputdir.c_str()));

  //std::string regionsSet = "13TeV_PHYS14";
  //std::string regionsSet = "13TeV_PHYS14_hiHT"; 
  //std::string regionsSet = "13TeV_PHYS14_hiJet_mergeHT";
  std::string regionsSet = "13TeV_PHYS14_loJet_hiHT";
  //std::string regionsSet = "13TeV_PHYS14_hiJet_extremeHT";
  //std::string regionsSet = "13TeV_CSA14";

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  
  MT2Analysis<MT2EstimateSyst>* lostLeptonEstimate = new MT2Analysis<MT2EstimateSyst> ( "llep", regionsSet );  
  for( unsigned i=0; i < fSamples.size(); ++i )
    (*lostLeptonEstimate) += ( computeYield( fSamples[i], regionsSet, lumi ) );
  
  lostLeptonEstimate->writeToFile(Form("llep_PHYS14_v3_%s.root", regionsSet.c_str()));

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

    if( myTree.nMuons10 > 0 || myTree.nElectrons10 > 0 || myTree.nPFLep5LowMT > 0 || myTree.nPFHad10LowMT > 0) ; //CR
    else continue;

    if( myTree.nVert==0 ) continue;
    if( myTree.nJet40<2 ) continue;
    //if( myTree.jet_pt[1]<100. ) continue;
    if( myTree.deltaPhiMin<0.3 ) continue;
    if( myTree.diffMetMht>0.5*myTree.met_pt ) continue;

    float jetCentral_pt[2];
    int njetsCentral = 0;
    for(int j=0; j<myTree.njet; ++j){
      if( fabs( myTree.jet_eta[j] ) < 2.5 ) {
        jetCentral_pt[njetsCentral] = myTree.jet_pt[j];
        ++njetsCentral;
      }
      if( njetsCentral >= 2 ) break;
    }
    if (jetCentral_pt[1] < 100. ) continue;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = myTree.mt2;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet40;
    int nbjets = myTree.nBJet40; 
    
    int nMuons10 = myTree.nMuons10;
    int nElectrons10 = myTree.nElectrons10;
    int nPFLep5LowMT = myTree.nPFLep5LowMT;
    int nPFHad10LowMT = myTree.nPFHad10LowMT;

    int nlepveto = myTree.nMuons10 + myTree.nElectrons10 + myTree.nPFLep5LowMT + myTree.nPFHad10LowMT;
    int nlep_unique = nlepveto;
    
    // do lepton overlap removal and 1L CR selections
    if (nlepveto >= 1) {
      std::vector<lepcand> all_cands;
      std::vector<lepcand> unique_cands;
      // check reco leptons - apply MT cut later
      // do overlap with PFcands below
      if ( nMuons10 > 0 || nElectrons10 > 0) {
	for (int ilep = 0; ilep < myTree.nlep; ++ilep) {
	  lepcand cand;
	  cand.pt = myTree.lep_pt[ilep];
	  cand.phi = myTree.lep_phi[ilep];
	  cand.mt = sqrt( 2 * myTree.met_pt * cand.pt * ( 1 - cos( myTree.met_phi - cand.phi) ) );
	  cand.eta = myTree.lep_eta[ilep];
	  cand.pdgId = myTree.lep_pdgId[ilep];
	  cand.isPFCand = false;

	  // add cand to vector
	  all_cands.push_back(cand);
	} // loop over reco leps
      }
      // pf leptons: need to find cands passing selection. 
      else if (nPFLep5LowMT > 0) {
	for (int itrk = 0; itrk < myTree.nisoTrack; ++itrk) {
	  lepcand cand;
	  cand.pt = myTree.isoTrack_pt[itrk];
	  cand.phi = myTree.isoTrack_phi[itrk];
	  cand.pdgId = myTree.isoTrack_pdgId[itrk];
	  if (cand.pt < 5.) continue;
	  if (abs(cand.pdgId) != 11 && abs(cand.pdgId) != 13) continue;
	  float absiso = myTree.isoTrack_absIso[itrk];
	  if (absiso/cand.pt > 0.2) continue;
	  cand.mt = sqrt( 2 * myTree.met_pt * cand.pt * ( 1 - cos( myTree.met_phi - cand.phi) ) );
	  cand.eta = myTree.isoTrack_eta[itrk];
	  cand.isPFCand = true;

	  // cand passes cuts: add to vector
	  if (cand.mt > 100.) continue;
	  all_cands.push_back(cand);
	} // loop on isoTracks
      }
      // pf hadrons: need to find cands passing selection. 
      else if (myTree.nPFHad10LowMT > 0) {
	for (int itrk = 0; itrk < myTree.nisoTrack; ++itrk) {
	  lepcand cand;
	  cand.pt = myTree.isoTrack_pt[itrk];
	  cand.phi = myTree.isoTrack_phi[itrk];
	  cand.pdgId = myTree.isoTrack_pdgId[itrk];
	  if (cand.pt < 10.) continue;
	  if (abs(cand.pdgId) != 211) continue;
	  float absiso = myTree.isoTrack_absIso[itrk];
	  if (absiso/cand.pt > 0.1) continue;
	  cand.mt = sqrt( 2 * myTree.met_pt * cand.pt * ( 1 - cos( myTree.met_phi - cand.phi) ) );
	  cand.eta = myTree.isoTrack_eta[itrk];
	  cand.isPFCand = true;

	  // cand passes cuts: add to vector
	  if (cand.mt > 100.) continue;
	  all_cands.push_back(cand);
	} // loop on isoTracks
      }

      // check all_cands for overlaps
      for (unsigned int icand = 0; icand < all_cands.size(); ++icand) {
	bool keep = true;
	for (unsigned int jcand = 0; jcand < all_cands.size(); ++jcand) {
	  float dr = DeltaR(all_cands.at(icand).eta, all_cands.at(jcand).eta, all_cands.at(icand).phi, all_cands.at(jcand).phi);
	  if (dr < 0.1) {
	    // if overlap, check whether the cands have the same pdgId
	    // keep the reco lepton in case of overlap with PF lepton
	    if (all_cands.at(icand).pdgId == all_cands.at(jcand).pdgId && 
		all_cands.at(icand).isPFCand && !all_cands.at(jcand).isPFCand) 
	      keep = false;
	  }
	}
	if (keep) unique_cands.push_back(all_cands.at(icand));
      }
      
      nlep_unique = unique_cands.size() ; // useful counter

      // check size of unique cands. if size == 1 and MT < 100, fill 1L CR plots
      if (unique_cands.size() == 1 && unique_cands.at(0).mt < 100);
      else continue;
      
    } // for 1L control region

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

float DeltaR(float eta1, float eta2, float phi1, float phi2){
  float dEta = eta1 - eta2;
  float dPhi = DeltaPhi(phi1, phi2);
  return TMath::Sqrt(dEta*dEta + dPhi*dPhi);
}

float DeltaPhi(float phi1, float phi2){
  float dPhi = phi1 - phi2;
  while (dPhi  >  TMath::Pi()) dPhi -= 2*TMath::Pi();
  while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
  return fabs(dPhi);
}
