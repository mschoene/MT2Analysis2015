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

#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"

#define mt2_cxx
#include "interface/mt2.h"


struct lepcand {
  
  float pt;
  float eta;
  float phi;
  int pdgId;
  float mt;
  bool isPFCand;

};
void computeYield( const MT2Sample& sample, const MT2Config& cfg, MT2Analysis<MT2EstimateTree>* anaTree );
//MT2Analysis<MT2EstimateTree> computeYield( const MT2Sample& sample, const MT2Config& cfg );

float DeltaR(float eta1, float eta2, float phi1, float phi2);

float DeltaPhi(float phi1, float phi2);

int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeLostLepton_CR             |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc!=2 ) {
    std::cout << "USAGE: ./computeLostLepton_CR [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  
  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 300, 599); // only top (tt, t, ttW, ttZ) and WJets
  //std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 100, 999); // only top (tt, t, ttW, ttZ) and WJets
  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  
  std::string regionsSet;
  if( cfg.regionsSet()=="zurich" )
    regionsSet="zurich_llep";
  else if( cfg.regionsSet()=="zurich_monojet" )
    regionsSet="zurich_monojet_llep";
  else
    regionsSet=cfg.regionsSet();

  std::cout << "using region set " << regionsSet << std::endl;
  
  //  std::string regionsSet_="13TeV_inclusive";

  MT2Analysis<MT2EstimateTree>* lostLeptonEstimate = new MT2Analysis<MT2EstimateTree> ( "llep", regionsSet );  
  //MT2Analysis<MT2EstimateTree>* lostLeptonEstimate = new MT2Analysis<MT2EstimateTree> ( "llep", regionsSet_ );  
//  for( unsigned i=0; i < fSamples.size(); ++i )
//    (*lostLeptonEstimate) += ( computeYield( fSamples[i], cfg ) );
  for( unsigned i=0; i < fSamples.size(); ++i )
    computeYield( fSamples[i], cfg, lostLeptonEstimate );

  TH1D::AddDirectory(kTRUE);
  
//  std::cout << "-> Making MT2EstimateTrees from inclusive tree (might take a sec)...";
//  std::string mcTruthSelection = "id>=151 && id<=157 && mt2>200.";
//  MT2Analysis<MT2EstimateTree>* QCD = MT2EstimateTree::makeAnalysisFromInclusiveTree( "QCD", regionsSet, lostLeptonEstimate, mcTruthSelection+"&&deltaPhiMin>0.3&&diffMetMht/met<0.5" ); 
//
//  mcTruthSelection = "id>=301 && id<=499 && mt2>200.";
//  MT2Analysis<MT2EstimateTree>* Top = MT2EstimateTree::makeAnalysisFromInclusiveTree( "Top", regionsSet, lostLeptonEstimate, mcTruthSelection+"&&deltaPhiMin>0.3&&diffMetMht/met<0.5" ); 
//
//  mcTruthSelection = "id>=501 && id<=599 && mt2>200.";
//  MT2Analysis<MT2EstimateTree>* WJets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "WJets", regionsSet, lostLeptonEstimate, mcTruthSelection+"&&deltaPhiMin>0.3&&diffMetMht/met<0.5" ); 
//  
//  mcTruthSelection = "id>=601 && id<=699 && mt2>200.";
//  MT2Analysis<MT2EstimateTree>* ZJets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "ZJets", regionsSet, lostLeptonEstimate, mcTruthSelection+"&&deltaPhiMin>0.3&&diffMetMht/met<0.5" ); 
//
//  mcTruthSelection = "id>=701 && id<=799 && mt2>200.";
//  MT2Analysis<MT2EstimateTree>* DYJets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "DYJets", regionsSet, lostLeptonEstimate, mcTruthSelection+"&&deltaPhiMin>0.3&&diffMetMht/met<0.5" ); 
//
//  mcTruthSelection = "id>=201 && id<=299 && mt2>200.";
//  MT2Analysis<MT2EstimateTree>* GJets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "GJets", regionsSet, lostLeptonEstimate, mcTruthSelection+"&&deltaPhiMin>0.3&&diffMetMht/met<0.5" ); 

  TH1D::AddDirectory(kFALSE);

  std::string outfile = cfg.getEventYieldDir() + "/llepEstimate.root";

  lostLeptonEstimate->writeToFile( outfile ); //Form("llep_%s_%s_%.0ffb.root", sampleName.c_str(), regionsSet.c_str(), lumi));
//  QCD->addToFile( outfile );
//  Top->addToFile( outfile );
//  WJets->addToFile( outfile );
//  ZJets->addToFile( outfile );
//  DYJets->addToFile( outfile );
//  GJets->addToFile( outfile );

  return 0;
  
}



void computeYield( const MT2Sample& sample, const MT2Config& cfg, MT2Analysis<MT2EstimateTree>* anaTree ){

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  MT2Tree myTree;
  myTree.Init(tree);

  std::string regionsSet;
  if( cfg.regionsSet()=="zurich" )
    regionsSet="zurich_llep";
  else if( cfg.regionsSet()=="zurich_monojet" )
    regionsSet="zurich_monojet_llep";
  else
    regionsSet=cfg.regionsSet();
  

//  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
//  MT2Analysis<MT2EstimateTree> analysis( sample.sname, regionsSet, sample.id );

  int nentries = tree->GetEntries();
    
  //ofstream ofs("events.log");

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( !myTree.passBaseline() ) continue;
    //if( myTree.passLeptonVeto() && myTree.passIsoTrackVeto() ) continue; // OLD lost lepton CR
    if( myTree.nLepLowMT==1 ) ; // New lost lepton CR
    else continue;

    if( myTree.nJet30==1 && (myTree.jet_id[0]<3 || myTree.jet_chHEF[0]<0.05 || myTree.jet_neHEF[0]>0.8 || myTree.jet_phEF[0]>0.7) ) continue;

    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20; 
    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = (njets>1) ? myTree.mt2 : ht;
    float minMTBmet = myTree.minMTBMet;
    
//    //////QCD
//    if( ht > 575. && ht < 1000. && mt2 < 300. ) continue;
//    else if( ht > 1000 && ht < 1500. && mt2 < 300. ) continue;
//    else if( ht > 1500. && mt2 < 400. ) continue;
//    //////

    int nMuons10 = myTree.nMuons10;
    int nElectrons10 = myTree.nElectrons10;
    int nPFLep5LowMT = myTree.nPFLep5LowMT;
    int nPFHad10LowMT = myTree.nPFHad10LowMT;

    int nlepveto = myTree.nMuons10 + myTree.nElectrons10 + myTree.nPFLep5LowMT + myTree.nPFHad10LowMT;
    int nlep_unique = nlepveto;
    
    //// do lepton overlap removal and 1L CR selections
    //if (nlepveto >= 1) {
    //  std::vector<lepcand> all_cands;
    //  std::vector<lepcand> unique_cands;
    //  // check reco leptons - apply MT cut later
    //  // do overlap with PFcands below
    //  if ( nMuons10 > 0 || nElectrons10 > 0) {
    //	for (int ilep = 0; ilep < myTree.nlep; ++ilep) {
    //	  lepcand cand;
    //	  cand.pt = myTree.lep_pt[ilep];
    //	  cand.phi = myTree.lep_phi[ilep];
    //	  cand.mt = sqrt( 2 * myTree.met_pt * cand.pt * ( 1 - cos( myTree.met_phi - cand.phi) ) );
    //	  cand.eta = myTree.lep_eta[ilep];
    //	  cand.pdgId = myTree.lep_pdgId[ilep];
    //	  cand.isPFCand = false;
    //
    //	  // add cand to vector
    //	  all_cands.push_back(cand);
    //	} // loop over reco leps
    //  }
    //  // pf leptons: need to find cands passing selection. 
    //  else if (nPFLep5LowMT > 0) {
    //	for (int itrk = 0; itrk < myTree.nisoTrack; ++itrk) {
    //	  lepcand cand;
    //	  cand.pt = myTree.isoTrack_pt[itrk];
    //	  cand.phi = myTree.isoTrack_phi[itrk];
    //	  cand.pdgId = myTree.isoTrack_pdgId[itrk];
    //	  if (cand.pt < 5.) continue;
    //	  if (abs(cand.pdgId) != 11 && abs(cand.pdgId) != 13) continue;
    //	  float absiso = myTree.isoTrack_absIso[itrk];
    //	  if (absiso/cand.pt > 0.2) continue;
    //	  cand.mt = sqrt( 2 * myTree.met_pt * cand.pt * ( 1 - cos( myTree.met_phi - cand.phi) ) );
    //	  cand.eta = myTree.isoTrack_eta[itrk];
    //	  cand.isPFCand = true;
    //
    //	  // cand passes cuts: add to vector
    //	  if (cand.mt > 100.) continue;
    //	  all_cands.push_back(cand);
    //	} // loop on isoTracks
    //  }
    //  // pf hadrons: need to find cands passing selection. 
    //  else if (myTree.nPFHad10LowMT > 0) {
    //	for (int itrk = 0; itrk < myTree.nisoTrack; ++itrk) {
    //	  lepcand cand;
    //	  cand.pt = myTree.isoTrack_pt[itrk];
    //	  cand.phi = myTree.isoTrack_phi[itrk];
    //	  cand.pdgId = myTree.isoTrack_pdgId[itrk];
    //	  if (cand.pt < 10.) continue;
    //	  if (abs(cand.pdgId) != 211) continue;
    //	  float absiso = myTree.isoTrack_absIso[itrk];
    //	  if (absiso/cand.pt > 0.1) continue;
    //	  cand.mt = sqrt( 2 * myTree.met_pt * cand.pt * ( 1 - cos( myTree.met_phi - cand.phi) ) );
    //	  cand.eta = myTree.isoTrack_eta[itrk];
    //	  cand.isPFCand = true;
    //
    //	  // cand passes cuts: add to vector
    //	  if (cand.mt > 100.) continue;
    //	  all_cands.push_back(cand);
    //	} // loop on isoTracks
    //  }
    //
    //  // check all_cands for overlaps
    //  for (unsigned int icand = 0; icand < all_cands.size(); ++icand) {
    //	bool keep = true;
    //	for (unsigned int jcand = 0; jcand < all_cands.size(); ++jcand) {
    //	  float dr = DeltaR(all_cands.at(icand).eta, all_cands.at(jcand).eta, all_cands.at(icand).phi, all_cands.at(jcand).phi);
    //	  if (dr < 0.1) {
    //	    // if overlap, check whether the cands have the same pdgId
    //	    // keep the reco lepton in case of overlap with PF lepton
    //	    if (all_cands.at(icand).pdgId == all_cands.at(jcand).pdgId && 
    //		all_cands.at(icand).isPFCand && !all_cands.at(jcand).isPFCand) 
    //	      keep = false;
    //	  }
    //	}
    //	if (keep) unique_cands.push_back(all_cands.at(icand));
    //  }
    //  
    //  nlep_unique = unique_cands.size() ; // useful counter
    //
    //  // check size of unique cands. if size == 1 and MT < 100, fill 1L CR plots
    //  if (unique_cands.size() == 1 && unique_cands.at(0).mt < 100);
    //  else continue;
    //  
    //} // for 1L control region

    Double_t weight = myTree.evt_scale1fb*cfg.lumi();

    //float fullweight_btagUp = weight;
    //float fullweight_btagDown = weight;

    MT2EstimateTree* thisEstimate = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisEstimate==0 ) continue;
    
    thisEstimate->assignTree( myTree, weight );
    thisEstimate->tree->Fill();
    thisEstimate->yield->Fill(mt2, weight );
    //thisEstimate->yield_btagUp  ->Fill(mt2, fullweight_btagUp );
    //thisEstimate->yield_btagDown->Fill(mt2, fullweight_btagDown );

    //ofs << "entry " << iEntry <<  "\tmet " << met << "\tmt2 " << mt2 << "\tminMTBmet " << minMTBmet << std::endl;
    
  } // for entries

  //ofs.close();

  anaTree->finalize();

  delete tree;

  file->Close();
  delete file;

  //  return analysis;

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
