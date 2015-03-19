//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 12 23:39:08 2014 by ROOT version 5.34/10
// from TTree mt2/treeProducerSusyFullHad
// found on file: /scratch/mmasciov/12December/POSTPROCESS/TTJets_MSDecaysCKM_central_PU_S14_POSTLS170_babytree_post.root
//////////////////////////////////////////////////////////

#ifndef mt2_h
#define mt2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MT2Tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Int_t           isData;
   Float_t         xsec;
   Int_t           HLT_htXprescale;
   Int_t           HLT_DoubleEl;
   Int_t           HLT_Photons;
   Int_t           HLT_SingleMu;
   Int_t           HLT_ht350met100;
   Int_t           HLT_MET150;
   Int_t           HLT_MuEG;
   Int_t           HLT_HT650;
   Int_t           HLT_DoubleMu;
   Float_t         puWeight;
   Int_t           nTrueInt;
   Float_t         genWeight;
   Float_t         rho;
   Int_t           nVert;
   Int_t           nJet25;
   Int_t           nBJetLoose25;
   Int_t           nBJetMedium25;
   Int_t           nBJetTight25;
   Int_t           nJet40;
   Int_t           nJet40a;
   Int_t           nBJetLoose40;
   Int_t           nBJetMedium40;
   Int_t           nBJetTight40;
   Int_t           nLepGood20;
   Int_t           nLepGood15;
   Int_t           nLepGood10;
   Int_t           GenHeaviestQCDFlavour;
   Float_t         LepEff_1lep;
   Float_t         LepEff_2lep;
   Int_t           GenSusyMScan1;
   Int_t           GenSusyMScan2;
   Int_t           GenSusyMScan3;
   Int_t           GenSusyMScan4;
   Int_t           GenSusyMGluino;
   Int_t           GenSusyMGravitino;
   Int_t           GenSusyMStop;
   Int_t           GenSusyMSbottom;
   Int_t           GenSusyMStop2;
   Int_t           GenSusyMSbottom2;
   Int_t           GenSusyMSquark;
   Int_t           GenSusyMNeutralino;
   Int_t           GenSusyMNeutralino2;
   Int_t           GenSusyMNeutralino3;
   Int_t           GenSusyMNeutralino4;
   Int_t           GenSusyMChargino;
   Int_t           GenSusyMChargino2;
   Float_t         genQScale;
   Float_t         ht;
   Float_t         mht_pt;
   Float_t         mht_phi;
   Float_t         diffMetMht;
   Float_t         deltaPhiMin;
   Float_t         ht_had;
   Float_t         mht_had_pt;
   Float_t         mht_had_phi;
   Float_t         diffMetMht_had;
   Float_t         deltaPhiMin_had;
   Int_t           nBJet40;
   Int_t           nMuons10;
   Int_t           nElectrons10;
   Int_t           nTaus20;
   Int_t           nGammas20;
   Float_t         minMTBMet;
   Int_t           nPFLep5LowMT;
   Int_t           nPFHad10LowMT;
   Float_t         mt2_had;
   Float_t         mt2ViaKt_had;
   Float_t         mt2_bb;
   Float_t         mt2_gen;
   Float_t         mt2;
   Float_t         gamma_mt2;
   Float_t         zll_mt2;
   Int_t           gamma_nJet40;
   Int_t           gamma_nBJet40;
   Float_t         gamma_ht;
   Float_t         gamma_deltaPhiMin;
   Float_t         gamma_diffMetMht;
   Float_t         gamma_mht_pt;
   Float_t         gamma_mht_phi;
   Float_t         gamma_minMTBMet;
   Float_t         zll_deltaPhiMin;
   Float_t         zll_diffMetMht;
   Float_t         zll_mht_pt;
   Float_t         zll_mht_phi;
   Float_t         zll_met_pt;
   Float_t         zll_met_phi;
   Float_t         zll_ht;
   Float_t         zll_pt;
   Float_t         zll_eta;
   Float_t         zll_phi;
   Float_t         zll_mass;
   Float_t         pseudoJet2_pt;
   Float_t         pseudoJet2_eta;
   Float_t         pseudoJet2_phi;
   Float_t         pseudoJet2_mass;
   Float_t         pseudoJet1_pt;
   Float_t         pseudoJet1_eta;
   Float_t         pseudoJet1_phi;
   Float_t         pseudoJet1_mass;
   Float_t         pseudoViaKtJet1_had_pt;
   Float_t         pseudoViaKtJet1_had_eta;
   Float_t         pseudoViaKtJet1_had_phi;
   Float_t         pseudoViaKtJet1_had_mass;
   Float_t         gamma_met_pt;
   Float_t         gamma_met_eta;
   Float_t         gamma_met_phi;
   Float_t         gamma_met_mass;
   Float_t         pseudoViaKtJet2_had_pt;
   Float_t         pseudoViaKtJet2_had_eta;
   Float_t         pseudoViaKtJet2_had_phi;
   Float_t         pseudoViaKtJet2_had_mass;
   Float_t         met_pt;
   Float_t         met_eta;
   Float_t         met_phi;
   Float_t         met_mass;
   Float_t         met_sumEt;
   Float_t         met_genPt;
   Float_t         met_genPhi;
   Float_t         met_genEta;
   Float_t         pseudoJet1_had_pt;
   Float_t         pseudoJet1_had_eta;
   Float_t         pseudoJet1_had_phi;
   Float_t         pseudoJet1_had_mass;
   Float_t         gamma_metNoPU_pt;
   Float_t         gamma_metNoPU_eta;
   Float_t         gamma_metNoPU_phi;
   Float_t         gamma_metNoPU_mass;
   Float_t         pseudoJet2_had_pt;
   Float_t         pseudoJet2_had_eta;
   Float_t         pseudoJet2_had_phi;
   Float_t         pseudoJet2_had_mass;
   Float_t         gamma_pseudoJet1_pt;
   Float_t         gamma_pseudoJet1_eta;
   Float_t         gamma_pseudoJet1_phi;
   Float_t         gamma_pseudoJet1_mass;
   Float_t         metNoPU_pt;
   Float_t         metNoPU_eta;
   Float_t         metNoPU_phi;
   Float_t         metNoPU_mass;
   Float_t         gamma_pseudoJet2_pt;
   Float_t         gamma_pseudoJet2_eta;
   Float_t         gamma_pseudoJet2_phi;
   Float_t         gamma_pseudoJet2_mass;
   Int_t           nfatJet;
   Float_t         fatJet_pt[9];   //[nfatJet]
   Float_t         fatJet_eta[9];   //[nfatJet]
   Float_t         fatJet_phi[9];   //[nfatJet]
   Float_t         fatJet_mass[9];   //[nfatJet]
   Float_t         fatJet_btagCSV[9];   //[nfatJet]
   Float_t         fatJet_rawPt[9];   //[nfatJet]
   Float_t         fatJet_mcPt[9];   //[nfatJet]
   Int_t           fatJet_mcFlavour[9];   //[nfatJet]
   Float_t         fatJet_prunedMass[9];   //[nfatJet]
   Float_t         fatJet_trimmedMass[9];   //[nfatJet]
   Float_t         fatJet_tau1[9];   //[nfatJet]
   Float_t         fatJet_tau2[9];   //[nfatJet]
   Float_t         fatJet_tau3[9];   //[nfatJet]
   Int_t           ngenPart;
   Float_t         genPart_pt[200];   //[ngenPart]
   Float_t         genPart_eta[200];   //[ngenPart]
   Float_t         genPart_phi[200];   //[ngenPart]
   Float_t         genPart_mass[200];   //[ngenPart]
   Int_t           genPart_pdgId[200];   //[ngenPart]
   Float_t         genPart_charge[200];   //[ngenPart]
   Int_t           genPart_status[200];   //[ngenPart]
   Int_t           genPart_motherId[200];   //[ngenPart]
   Int_t           genPart_motherIndex[200];   //[ngenPart]
   Int_t           genPart_grandmaId[200];   //[ngenPart]
   Int_t           ngenTau;
   Float_t         genTau_pt[2];   //[ngenTau]
   Float_t         genTau_eta[2];   //[ngenTau]
   Float_t         genTau_phi[2];   //[ngenTau]
   Float_t         genTau_mass[2];   //[ngenTau]
   Int_t           genTau_pdgId[2];   //[ngenTau]
   Float_t         genTau_charge[2];   //[ngenTau]
   Int_t           genTau_status[2];   //[ngenTau]
   Int_t           genTau_sourceId[2];   //[ngenTau]
   Int_t           nisoTrack;
   Float_t         isoTrack_pt[8];   //[nisoTrack]
   Float_t         isoTrack_eta[8];   //[nisoTrack]
   Float_t         isoTrack_phi[8];   //[nisoTrack]
   Float_t         isoTrack_mass[8];   //[nisoTrack]
   Int_t           isoTrack_pdgId[8];   //[nisoTrack]
   Int_t           isoTrack_charge[8];   //[nisoTrack]
   Float_t         isoTrack_dz[8];   //[nisoTrack]
   Float_t         isoTrack_absIso[8];   //[nisoTrack]
   Int_t           isoTrack_mcMatchId[8];   //[nisoTrack]
   Int_t           njet;
   Float_t         jet_pt[50];   //[njet]
   Float_t         jet_eta[50];   //[njet]
   Float_t         jet_phi[50];   //[njet]
   Float_t         jet_mass[50];   //[njet]
   Float_t         jet_btagCSV[50];   //[njet]
   Float_t         jet_rawPt[50];   //[njet]
   Float_t         jet_mcPt[50];   //[njet]
   Int_t           jet_mcFlavour[50];   //[njet]
   Int_t           jet_mcMatchId[50];   //[njet]
   Int_t           jet_mcMatchFlav[50];   //[njet]
   Float_t         jet_area[50];   //[njet]
   Int_t           jet_puId[50];   //[njet]
   Int_t           jet_id[50];   //[njet]
   Float_t         jet_ptd[50];   //[njet]
   Float_t         jet_axis2[50];   //[njet]
   Int_t           jet_mult[50];   //[njet]
   Float_t         jet_qgl[50];   //[njet]
   Int_t           jet_partonId[50];   //[njet]
   Int_t           jet_partonMotherId[50];   //[njet]
   Int_t           ngenNu;
   Float_t         genNu_pt[9];   //[ngenNu]
   Float_t         genNu_eta[9];   //[ngenNu]
   Float_t         genNu_phi[9];   //[ngenNu]
   Float_t         genNu_mass[9];   //[ngenNu]
   Int_t           genNu_pdgId[9];   //[ngenNu]
   Float_t         genNu_charge[9];   //[ngenNu]
   Int_t           genNu_status[9];   //[ngenNu]
   Int_t           genNu_sourceId[9];   //[ngenNu]
   Int_t           nlep;
   Float_t         lep_pt[5];   //[nlep]
   Float_t         lep_eta[5];   //[nlep]
   Float_t         lep_phi[5];   //[nlep]
   Float_t         lep_mass[5];   //[nlep]
   Int_t           lep_pdgId[5];   //[nlep]
   Int_t           lep_charge[5];   //[nlep]
   Float_t         lep_dxy[5];   //[nlep]
   Float_t         lep_dz[5];   //[nlep]
   Float_t         lep_edxy[5];   //[nlep]
   Float_t         lep_edz[5];   //[nlep]
   Float_t         lep_ip3d[5];   //[nlep]
   Float_t         lep_sip3d[5];   //[nlep]
   Int_t           lep_tightId[5];   //[nlep]
   Int_t           lep_convVeto[5];   //[nlep]
   Int_t           lep_lostHits[5];   //[nlep]
   Int_t           lep_looseIdSusy[5];   //[nlep]
   Float_t         lep_relIso03[5];   //[nlep]
   Float_t         lep_relIso04[5];   //[nlep]
   Float_t         lep_chargedHadRelIso03[5];   //[nlep]
   Float_t         lep_chargedHadRelIso04[5];   //[nlep]
   Int_t           lep_convVetoFull[5];   //[nlep]
   Int_t           lep_eleCutIdCSA14_25ns_v1[5];   //[nlep]
   Int_t           lep_eleCutIdCSA14_50ns_v1[5];   //[nlep]
   Int_t           lep_eleMVAId[5];   //[nlep]
   Int_t           lep_tightCharge[5];   //[nlep]
   Float_t         lep_mvaId[5];   //[nlep]
   Float_t         lep_mvaIdTrig[5];   //[nlep]
   Float_t         lep_nStations[5];   //[nlep]
   Float_t         lep_trkKink[5];   //[nlep]
   Float_t         lep_caloCompatibility[5];   //[nlep]
   Float_t         lep_globalTrackChi2[5];   //[nlep]
   Int_t           lep_trackerLayers[5];   //[nlep]
   Int_t           lep_pixelLayers[5];   //[nlep]
   Float_t         lep_mvaTTH[5];   //[nlep]
   Float_t         lep_mvaSusy[5];   //[nlep]
   Float_t         lep_jetPtRatio[5];   //[nlep]
   Float_t         lep_jetBTagCSV[5];   //[nlep]
   Float_t         lep_jetDR[5];   //[nlep]
   Int_t           lep_mcMatchId[5];   //[nlep]
   Int_t           lep_mcMatchAny[5];   //[nlep]
   Int_t           lep_mcMatchTau[5];   //[nlep]
   Int_t           ngenLepFromTau;
   Float_t         genLepFromTau_pt[5];   //[ngenLepFromTau]
   Float_t         genLepFromTau_eta[5];   //[ngenLepFromTau]
   Float_t         genLepFromTau_phi[5];   //[ngenLepFromTau]
   Float_t         genLepFromTau_mass[5];   //[ngenLepFromTau]
   Int_t           genLepFromTau_pdgId[5];   //[ngenLepFromTau]
   Float_t         genLepFromTau_charge[5];   //[ngenLepFromTau]
   Int_t           genLepFromTau_status[5];   //[ngenLepFromTau]
   Int_t           genLepFromTau_sourceId[5];   //[ngenLepFromTau]
   Int_t           ngamma;
   Float_t         gamma_pt[20];   //[ngamma]
   Float_t         gamma_eta[20];   //[ngamma]
   Float_t         gamma_phi[20];   //[ngamma]
   Float_t         gamma_mass[20];   //[ngamma]
   Int_t           gamma_pdgId[20];   //[ngamma]
   Int_t           gamma_idCutBased[20];   //[ngamma]
   Float_t         gamma_hOverE[20];   //[ngamma]
   Float_t         gamma_r9[20];   //[ngamma]
   Float_t         gamma_sigmaIetaIeta[20];   //[ngamma]
   Float_t         gamma_chHadIso[20];   //[ngamma]
   Float_t         gamma_chHadIso04[20];   //[ngamma]
   Float_t         gamma_chHadIsoRC[20];   //[ngamma]
   Float_t         gamma_chHadIsoRC04[20];   //[ngamma]
   Float_t         gamma_neuHadIso[20];   //[ngamma]
   Float_t         gamma_phIso[20];   //[ngamma]
   Int_t           gamma_mcMatchId[20];   //[ngamma]
   Float_t         gamma_genIso[20];   //[ngamma]
   Float_t         gamma_genIso04[20];   //[ngamma]
   Int_t           nGenP6StatusThree;
   Float_t         GenP6StatusThree_pt[54];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_eta[54];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_phi[54];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_mass[54];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_pdgId[54];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_charge[54];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_status[54];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_motherId[54];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_motherIndex[54];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_grandmaId[54];   //[nGenP6StatusThree]
   Int_t           ngenLep;
   Float_t         genLep_pt[2];   //[ngenLep]
   Float_t         genLep_eta[2];   //[ngenLep]
   Float_t         genLep_phi[2];   //[ngenLep]
   Float_t         genLep_mass[2];   //[ngenLep]
   Int_t           genLep_pdgId[2];   //[ngenLep]
   Float_t         genLep_charge[2];   //[ngenLep]
   Int_t           genLep_status[2];   //[ngenLep]
   Int_t           genLep_sourceId[2];   //[ngenLep]
   Int_t           ntau;
   Float_t         tau_pt[6];   //[ntau]
   Float_t         tau_eta[6];   //[ntau]
   Float_t         tau_phi[6];   //[ntau]
   Float_t         tau_mass[6];   //[ntau]
   Int_t           tau_pdgId[6];   //[ntau]
   Int_t           tau_charge[6];   //[ntau]
   Float_t         tau_dxy[6];   //[ntau]
   Float_t         tau_dz[6];   //[ntau]
   Int_t           tau_idCI3hit[6];   //[ntau]
   Float_t         tau_isoCI3hit[6];   //[ntau]
   Int_t           tau_mcMatchId[6];   //[ntau]
   Float_t         evt_scale1fb;
   Float_t         evt_xsec;
   Float_t         evt_kfactor;
   Float_t         evt_filter;
   Int_t           evt_nEvts;
   Int_t           evt_id;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_HLT_htXprescale;   //!
   TBranch        *b_HLT_DoubleEl;   //!
   TBranch        *b_HLT_Photons;   //!
   TBranch        *b_HLT_SingleMu;   //!
   TBranch        *b_HLT_ht350met100;   //!
   TBranch        *b_HLT_MET150;   //!
   TBranch        *b_HLT_MuEG;   //!
   TBranch        *b_HLT_HT650;   //!
   TBranch        *b_HLT_DoubleMu;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_nJet25;   //!
   TBranch        *b_nBJetLoose25;   //!
   TBranch        *b_nBJetMedium25;   //!
   TBranch        *b_nBJetTight25;   //!
   TBranch        *b_nJet40;   //!
   TBranch        *b_nJet40a;   //!
   TBranch        *b_nBJetLoose40;   //!
   TBranch        *b_nBJetMedium40;   //!
   TBranch        *b_nBJetTight40;   //!
   TBranch        *b_nLepGood20;   //!
   TBranch        *b_nLepGood15;   //!
   TBranch        *b_nLepGood10;   //!
   TBranch        *b_GenHeaviestQCDFlavour;   //!
   TBranch        *b_LepEff_1lep;   //!
   TBranch        *b_LepEff_2lep;   //!
   TBranch        *b_GenSusyMScan1;   //!
   TBranch        *b_GenSusyMScan2;   //!
   TBranch        *b_GenSusyMScan3;   //!
   TBranch        *b_GenSusyMScan4;   //!
   TBranch        *b_GenSusyMGluino;   //!
   TBranch        *b_GenSusyMGravitino;   //!
   TBranch        *b_GenSusyMStop;   //!
   TBranch        *b_GenSusyMSbottom;   //!
   TBranch        *b_GenSusyMStop2;   //!
   TBranch        *b_GenSusyMSbottom2;   //!
   TBranch        *b_GenSusyMSquark;   //!
   TBranch        *b_GenSusyMNeutralino;   //!
   TBranch        *b_GenSusyMNeutralino2;   //!
   TBranch        *b_GenSusyMNeutralino3;   //!
   TBranch        *b_GenSusyMNeutralino4;   //!
   TBranch        *b_GenSusyMChargino;   //!
   TBranch        *b_GenSusyMChargino2;   //!
   TBranch        *b_genQScale;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_mht_pt;   //!
   TBranch        *b_mht_phi;   //!
   TBranch        *b_diffMetMht;   //!
   TBranch        *b_deltaPhiMin;   //!
   TBranch        *b_ht_had;   //!
   TBranch        *b_mht_had_pt;   //!
   TBranch        *b_mht_had_phi;   //!
   TBranch        *b_diffMetMht_had;   //!
   TBranch        *b_deltaPhiMin_had;   //!
   TBranch        *b_nBJet40;   //!
   TBranch        *b_nMuons10;   //!
   TBranch        *b_nElectrons10;   //!
   TBranch        *b_nTaus20;   //!
   TBranch        *b_nGammas20;   //!
   TBranch        *b_minMTBMet;   //!
   TBranch        *b_nPFLep5LowMT;   //!
   TBranch        *b_nPFHad10LowMT;   //!
   TBranch        *b_mt2_had;   //!
   TBranch        *b_mt2ViaKt_had;   //!
   TBranch        *b_mt2_bb;   //!
   TBranch        *b_mt2_gen;   //!
   TBranch        *b_mt2;   //!
   TBranch        *b_gamma_mt2;   //!
   TBranch        *b_zll_mt2;   //!
   TBranch        *b_gamma_nJet40;   //!
   TBranch        *b_gamma_nBJet40;   //!
   TBranch        *b_gamma_ht;   //!
   TBranch        *b_gamma_deltaPhiMin;   //!
   TBranch        *b_gamma_diffMetMht;   //!
   TBranch        *b_gamma_mht_pt;   //!
   TBranch        *b_gamma_mht_phi;   //!
   TBranch        *b_gamma_minMTBMet;   //!
   TBranch        *b_zll_deltaPhiMin;   //!
   TBranch        *b_zll_diffMetMht;   //!
   TBranch        *b_zll_mht_pt;   //!
   TBranch        *b_zll_mht_phi;   //!
   TBranch        *b_zll_met_pt;   //!
   TBranch        *b_zll_met_phi;   //!
   TBranch        *b_zll_ht;   //!
   TBranch        *b_zll_pt;   //!
   TBranch        *b_zll_eta;   //!
   TBranch        *b_zll_phi;   //!
   TBranch        *b_zll_mass;   //!
   TBranch        *b_pseudoJet2_pt;   //!
   TBranch        *b_pseudoJet2_eta;   //!
   TBranch        *b_pseudoJet2_phi;   //!
   TBranch        *b_pseudoJet2_mass;   //!
   TBranch        *b_pseudoJet1_pt;   //!
   TBranch        *b_pseudoJet1_eta;   //!
   TBranch        *b_pseudoJet1_phi;   //!
   TBranch        *b_pseudoJet1_mass;   //!
   TBranch        *b_pseudoViaKtJet1_had_pt;   //!
   TBranch        *b_pseudoViaKtJet1_had_eta;   //!
   TBranch        *b_pseudoViaKtJet1_had_phi;   //!
   TBranch        *b_pseudoViaKtJet1_had_mass;   //!
   TBranch        *b_gamma_met_pt;   //!
   TBranch        *b_gamma_met_eta;   //!
   TBranch        *b_gamma_met_phi;   //!
   TBranch        *b_gamma_met_mass;   //!
   TBranch        *b_pseudoViaKtJet2_had_pt;   //!
   TBranch        *b_pseudoViaKtJet2_had_eta;   //!
   TBranch        *b_pseudoViaKtJet2_had_phi;   //!
   TBranch        *b_pseudoViaKtJet2_had_mass;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_pseudoJet1_had_pt;   //!
   TBranch        *b_pseudoJet1_had_eta;   //!
   TBranch        *b_pseudoJet1_had_phi;   //!
   TBranch        *b_pseudoJet1_had_mass;   //!
   TBranch        *b_gamma_metNoPU_pt;   //!
   TBranch        *b_gamma_metNoPU_eta;   //!
   TBranch        *b_gamma_metNoPU_phi;   //!
   TBranch        *b_gamma_metNoPU_mass;   //!
   TBranch        *b_pseudoJet2_had_pt;   //!
   TBranch        *b_pseudoJet2_had_eta;   //!
   TBranch        *b_pseudoJet2_had_phi;   //!
   TBranch        *b_pseudoJet2_had_mass;   //!
   TBranch        *b_gamma_pseudoJet1_pt;   //!
   TBranch        *b_gamma_pseudoJet1_eta;   //!
   TBranch        *b_gamma_pseudoJet1_phi;   //!
   TBranch        *b_gamma_pseudoJet1_mass;   //!
   TBranch        *b_metNoPU_pt;   //!
   TBranch        *b_metNoPU_eta;   //!
   TBranch        *b_metNoPU_phi;   //!
   TBranch        *b_metNoPU_mass;   //!
   TBranch        *b_gamma_pseudoJet2_pt;   //!
   TBranch        *b_gamma_pseudoJet2_eta;   //!
   TBranch        *b_gamma_pseudoJet2_phi;   //!
   TBranch        *b_gamma_pseudoJet2_mass;   //!
   TBranch        *b_nfatJet;   //!
   TBranch        *b_fatJet_pt;   //!
   TBranch        *b_fatJet_eta;   //!
   TBranch        *b_fatJet_phi;   //!
   TBranch        *b_fatJet_mass;   //!
   TBranch        *b_fatJet_btagCSV;   //!
   TBranch        *b_fatJet_rawPt;   //!
   TBranch        *b_fatJet_mcPt;   //!
   TBranch        *b_fatJet_mcFlavour;   //!
   TBranch        *b_fatJet_prunedMass;   //!
   TBranch        *b_fatJet_trimmedMass;   //!
   TBranch        *b_fatJet_tau1;   //!
   TBranch        *b_fatJet_tau2;   //!
   TBranch        *b_fatJet_tau3;   //!
   TBranch        *b_ngenPart;   //!
   TBranch        *b_genPart_pt;   //!
   TBranch        *b_genPart_eta;   //!
   TBranch        *b_genPart_phi;   //!
   TBranch        *b_genPart_mass;   //!
   TBranch        *b_genPart_pdgId;   //!
   TBranch        *b_genPart_charge;   //!
   TBranch        *b_genPart_status;   //!
   TBranch        *b_genPart_motherId;   //!
   TBranch        *b_genPart_motherIndex;   //!
   TBranch        *b_genPart_grandmaId;   //!
   TBranch        *b_ngenTau;   //!
   TBranch        *b_genTau_pt;   //!
   TBranch        *b_genTau_eta;   //!
   TBranch        *b_genTau_phi;   //!
   TBranch        *b_genTau_mass;   //!
   TBranch        *b_genTau_pdgId;   //!
   TBranch        *b_genTau_charge;   //!
   TBranch        *b_genTau_status;   //!
   TBranch        *b_genTau_sourceId;   //!
   TBranch        *b_nisoTrack;   //!
   TBranch        *b_isoTrack_pt;   //!
   TBranch        *b_isoTrack_eta;   //!
   TBranch        *b_isoTrack_phi;   //!
   TBranch        *b_isoTrack_mass;   //!
   TBranch        *b_isoTrack_pdgId;   //!
   TBranch        *b_isoTrack_charge;   //!
   TBranch        *b_isoTrack_dz;   //!
   TBranch        *b_isoTrack_absIso;   //!
   TBranch        *b_isoTrack_mcMatchId;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_btagCSV;   //!
   TBranch        *b_jet_rawPt;   //!
   TBranch        *b_jet_mcPt;   //!
   TBranch        *b_jet_mcFlavour;   //!
   TBranch        *b_jet_mcMatchId;   //!
   TBranch        *b_jet_mcMatchFlav;   //!
   TBranch        *b_jet_area;   //!
   TBranch        *b_jet_puId;   //!
   TBranch        *b_jet_id;   //!
   TBranch        *b_jet_ptd;   //!
   TBranch        *b_jet_axis2;   //!
   TBranch        *b_jet_mult;   //!
   TBranch        *b_jet_qgl;   //!
   TBranch        *b_jet_partonId;   //!
   TBranch        *b_jet_partonMotherId;   //!
   TBranch        *b_ngenNu;   //!
   TBranch        *b_genNu_pt;   //!
   TBranch        *b_genNu_eta;   //!
   TBranch        *b_genNu_phi;   //!
   TBranch        *b_genNu_mass;   //!
   TBranch        *b_genNu_pdgId;   //!
   TBranch        *b_genNu_charge;   //!
   TBranch        *b_genNu_status;   //!
   TBranch        *b_genNu_sourceId;   //!
   TBranch        *b_nlep;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_mass;   //!
   TBranch        *b_lep_pdgId;   //!
   TBranch        *b_lep_charge;   //!
   TBranch        *b_lep_dxy;   //!
   TBranch        *b_lep_dz;   //!
   TBranch        *b_lep_edxy;   //!
   TBranch        *b_lep_edz;   //!
   TBranch        *b_lep_ip3d;   //!
   TBranch        *b_lep_sip3d;   //!
   TBranch        *b_lep_tightId;   //!
   TBranch        *b_lep_convVeto;   //!
   TBranch        *b_lep_lostHits;   //!
   TBranch        *b_lep_looseIdSusy;   //!
   TBranch        *b_lep_relIso03;   //!
   TBranch        *b_lep_relIso04;   //!
   TBranch        *b_lep_chargedHadRelIso03;   //!
   TBranch        *b_lep_chargedHadRelIso04;   //!
   TBranch        *b_lep_convVetoFull;   //!
   TBranch        *b_lep_eleCutIdCSA14_25ns_v1;   //!
   TBranch        *b_lep_eleCutIdCSA14_50ns_v1;   //!
   TBranch        *b_lep_eleMVAId;   //!
   TBranch        *b_lep_tightCharge;   //!
   TBranch        *b_lep_mvaId;   //!
   TBranch        *b_lep_mvaIdTrig;   //!
   TBranch        *b_lep_nStations;   //!
   TBranch        *b_lep_trkKink;   //!
   TBranch        *b_lep_caloCompatibility;   //!
   TBranch        *b_lep_globalTrackChi2;   //!
   TBranch        *b_lep_trackerLayers;   //!
   TBranch        *b_lep_pixelLayers;   //!
   TBranch        *b_lep_mvaTTH;   //!
   TBranch        *b_lep_mvaSusy;   //!
   TBranch        *b_lep_jetPtRatio;   //!
   TBranch        *b_lep_jetBTagCSV;   //!
   TBranch        *b_lep_jetDR;   //!
   TBranch        *b_lep_mcMatchId;   //!
   TBranch        *b_lep_mcMatchAny;   //!
   TBranch        *b_lep_mcMatchTau;   //!
   TBranch        *b_ngenLepFromTau;   //!
   TBranch        *b_genLepFromTau_pt;   //!
   TBranch        *b_genLepFromTau_eta;   //!
   TBranch        *b_genLepFromTau_phi;   //!
   TBranch        *b_genLepFromTau_mass;   //!
   TBranch        *b_genLepFromTau_pdgId;   //!
   TBranch        *b_genLepFromTau_charge;   //!
   TBranch        *b_genLepFromTau_status;   //!
   TBranch        *b_genLepFromTau_sourceId;   //!
   TBranch        *b_ngamma;   //!
   TBranch        *b_gamma_pt;   //!
   TBranch        *b_gamma_eta;   //!
   TBranch        *b_gamma_phi;   //!
   TBranch        *b_gamma_mass;   //!
   TBranch        *b_gamma_pdgId;   //!
   TBranch        *b_gamma_idCutBased;   //!
   TBranch        *b_gamma_hOverE;   //!
   TBranch        *b_gamma_r9;   //!
   TBranch        *b_gamma_sigmaIetaIeta;   //!
   TBranch        *b_gamma_chHadIso;   //!
   TBranch        *b_gamma_chHadIso04;   //!
   TBranch        *b_gamma_chHadIsoRC;   //!
   TBranch        *b_gamma_chHadIsoRC04;   //!
   TBranch        *b_gamma_neuHadIso;   //!
   TBranch        *b_gamma_phIso;   //!
   TBranch        *b_gamma_mcMatchId;   //!
   TBranch        *b_gamma_genIso;   //!
   TBranch        *b_gamma_genIso04;   //!
   TBranch        *b_nGenP6StatusThree;   //!
   TBranch        *b_GenP6StatusThree_pt;   //!
   TBranch        *b_GenP6StatusThree_eta;   //!
   TBranch        *b_GenP6StatusThree_phi;   //!
   TBranch        *b_GenP6StatusThree_mass;   //!
   TBranch        *b_GenP6StatusThree_pdgId;   //!
   TBranch        *b_GenP6StatusThree_charge;   //!
   TBranch        *b_GenP6StatusThree_status;   //!
   TBranch        *b_GenP6StatusThree_motherId;   //!
   TBranch        *b_GenP6StatusThree_motherIndex;   //!
   TBranch        *b_GenP6StatusThree_grandmaId;   //!
   TBranch        *b_ngenLep;   //!
   TBranch        *b_genLep_pt;   //!
   TBranch        *b_genLep_eta;   //!
   TBranch        *b_genLep_phi;   //!
   TBranch        *b_genLep_mass;   //!
   TBranch        *b_genLep_pdgId;   //!
   TBranch        *b_genLep_charge;   //!
   TBranch        *b_genLep_status;   //!
   TBranch        *b_genLep_sourceId;   //!
   TBranch        *b_ntau;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_mass;   //!
   TBranch        *b_tau_pdgId;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_dxy;   //!
   TBranch        *b_tau_dz;   //!
   TBranch        *b_tau_idCI3hit;   //!
   TBranch        *b_tau_isoCI3hit;   //!
   TBranch        *b_tau_mcMatchId;   //!
   TBranch        *b_evt_scale1fb;   //!
   TBranch        *b_evt_xsec;   //!
   TBranch        *b_evt_kfactor;   //!
   TBranch        *b_evt_filter;   //!
   TBranch        *b_evt_nEvts;   //!
   TBranch        *b_evt_id;   //!

   bool loadGenStuff;

   MT2Tree(TTree *tree=0);
   virtual ~MT2Tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mt2_cxx
MT2Tree::MT2Tree(TTree *tree) : fChain(0) 
{

  loadGenStuff = true;

}

MT2Tree::~MT2Tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MT2Tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MT2Tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MT2Tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
   fChain->SetBranchAddress("HLT_htXprescale", &HLT_htXprescale, &b_HLT_htXprescale);
   fChain->SetBranchAddress("HLT_DoubleEl", &HLT_DoubleEl, &b_HLT_DoubleEl);
   fChain->SetBranchAddress("HLT_Photons", &HLT_Photons, &b_HLT_Photons);
   fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu, &b_HLT_SingleMu);
   fChain->SetBranchAddress("HLT_ht350met100", &HLT_ht350met100, &b_HLT_ht350met100);
   fChain->SetBranchAddress("HLT_MET150", &HLT_MET150, &b_HLT_MET150);
   fChain->SetBranchAddress("HLT_MuEG", &HLT_MuEG, &b_HLT_MuEG);
   fChain->SetBranchAddress("HLT_HT650", &HLT_HT650, &b_HLT_HT650);
   fChain->SetBranchAddress("HLT_DoubleMu", &HLT_DoubleMu, &b_HLT_DoubleMu);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("nJet25", &nJet25, &b_nJet25);
   fChain->SetBranchAddress("nBJetLoose25", &nBJetLoose25, &b_nBJetLoose25);
   fChain->SetBranchAddress("nBJetMedium25", &nBJetMedium25, &b_nBJetMedium25);
   fChain->SetBranchAddress("nBJetTight25", &nBJetTight25, &b_nBJetTight25);
   fChain->SetBranchAddress("nJet40", &nJet40, &b_nJet40);
   fChain->SetBranchAddress("nJet40a", &nJet40a, &b_nJet40a);
   fChain->SetBranchAddress("nBJetLoose40", &nBJetLoose40, &b_nBJetLoose40);
   fChain->SetBranchAddress("nBJetMedium40", &nBJetMedium40, &b_nBJetMedium40);
   fChain->SetBranchAddress("nBJetTight40", &nBJetTight40, &b_nBJetTight40);
   fChain->SetBranchAddress("nLepGood20", &nLepGood20, &b_nLepGood20);
   fChain->SetBranchAddress("nLepGood15", &nLepGood15, &b_nLepGood15);
   fChain->SetBranchAddress("nLepGood10", &nLepGood10, &b_nLepGood10);
   fChain->SetBranchAddress("GenHeaviestQCDFlavour", &GenHeaviestQCDFlavour, &b_GenHeaviestQCDFlavour);
   fChain->SetBranchAddress("LepEff_1lep", &LepEff_1lep, &b_LepEff_1lep);
   fChain->SetBranchAddress("LepEff_2lep", &LepEff_2lep, &b_LepEff_2lep);
   fChain->SetBranchAddress("GenSusyMScan1", &GenSusyMScan1, &b_GenSusyMScan1);
   fChain->SetBranchAddress("GenSusyMScan2", &GenSusyMScan2, &b_GenSusyMScan2);
   fChain->SetBranchAddress("GenSusyMScan3", &GenSusyMScan3, &b_GenSusyMScan3);
   fChain->SetBranchAddress("GenSusyMScan4", &GenSusyMScan4, &b_GenSusyMScan4);
   fChain->SetBranchAddress("GenSusyMGluino", &GenSusyMGluino, &b_GenSusyMGluino);
   fChain->SetBranchAddress("GenSusyMGravitino", &GenSusyMGravitino, &b_GenSusyMGravitino);
   fChain->SetBranchAddress("GenSusyMStop", &GenSusyMStop, &b_GenSusyMStop);
   fChain->SetBranchAddress("GenSusyMSbottom", &GenSusyMSbottom, &b_GenSusyMSbottom);
   fChain->SetBranchAddress("GenSusyMStop2", &GenSusyMStop2, &b_GenSusyMStop2);
   fChain->SetBranchAddress("GenSusyMSbottom2", &GenSusyMSbottom2, &b_GenSusyMSbottom2);
   fChain->SetBranchAddress("GenSusyMSquark", &GenSusyMSquark, &b_GenSusyMSquark);
   fChain->SetBranchAddress("GenSusyMNeutralino", &GenSusyMNeutralino, &b_GenSusyMNeutralino);
   fChain->SetBranchAddress("GenSusyMNeutralino2", &GenSusyMNeutralino2, &b_GenSusyMNeutralino2);
   fChain->SetBranchAddress("GenSusyMNeutralino3", &GenSusyMNeutralino3, &b_GenSusyMNeutralino3);
   fChain->SetBranchAddress("GenSusyMNeutralino4", &GenSusyMNeutralino4, &b_GenSusyMNeutralino4);
   fChain->SetBranchAddress("GenSusyMChargino", &GenSusyMChargino, &b_GenSusyMChargino);
   fChain->SetBranchAddress("GenSusyMChargino2", &GenSusyMChargino2, &b_GenSusyMChargino2);
   fChain->SetBranchAddress("genQScale", &genQScale, &b_genQScale);
   fChain->SetBranchAddress("ht", &ht, &b_ht);
   fChain->SetBranchAddress("mht_pt", &mht_pt, &b_mht_pt);
   fChain->SetBranchAddress("mht_phi", &mht_phi, &b_mht_phi);
   fChain->SetBranchAddress("diffMetMht", &diffMetMht, &b_diffMetMht);
   fChain->SetBranchAddress("deltaPhiMin", &deltaPhiMin, &b_deltaPhiMin);
   fChain->SetBranchAddress("ht_had", &ht_had, &b_ht_had);
   fChain->SetBranchAddress("mht_had_pt", &mht_had_pt, &b_mht_had_pt);
   fChain->SetBranchAddress("mht_had_phi", &mht_had_phi, &b_mht_had_phi);
   fChain->SetBranchAddress("diffMetMht_had", &diffMetMht_had, &b_diffMetMht_had);
   fChain->SetBranchAddress("deltaPhiMin_had", &deltaPhiMin_had, &b_deltaPhiMin_had);
   fChain->SetBranchAddress("nBJet40", &nBJet40, &b_nBJet40);
   fChain->SetBranchAddress("nMuons10", &nMuons10, &b_nMuons10);
   fChain->SetBranchAddress("nElectrons10", &nElectrons10, &b_nElectrons10);
   fChain->SetBranchAddress("nTaus20", &nTaus20, &b_nTaus20);
   fChain->SetBranchAddress("nGammas20", &nGammas20, &b_nGammas20);
   fChain->SetBranchAddress("minMTBMet", &minMTBMet, &b_minMTBMet);
   fChain->SetBranchAddress("nPFLep5LowMT", &nPFLep5LowMT, &b_nPFLep5LowMT);
   fChain->SetBranchAddress("nPFHad10LowMT", &nPFHad10LowMT, &b_nPFHad10LowMT);
   fChain->SetBranchAddress("mt2_had", &mt2_had, &b_mt2_had);
   fChain->SetBranchAddress("mt2ViaKt_had", &mt2ViaKt_had, &b_mt2ViaKt_had);
   fChain->SetBranchAddress("mt2_bb", &mt2_bb, &b_mt2_bb);
   fChain->SetBranchAddress("mt2_gen", &mt2_gen, &b_mt2_gen);
   fChain->SetBranchAddress("mt2", &mt2, &b_mt2);
   fChain->SetBranchAddress("gamma_mt2", &gamma_mt2, &b_gamma_mt2);
   fChain->SetBranchAddress("zll_mt2", &zll_mt2, &b_zll_mt2);
   fChain->SetBranchAddress("gamma_nJet40", &gamma_nJet40, &b_gamma_nJet40);
   fChain->SetBranchAddress("gamma_nBJet40", &gamma_nBJet40, &b_gamma_nBJet40);
   fChain->SetBranchAddress("gamma_ht", &gamma_ht, &b_gamma_ht);
   fChain->SetBranchAddress("gamma_deltaPhiMin", &gamma_deltaPhiMin, &b_gamma_deltaPhiMin);
   fChain->SetBranchAddress("gamma_diffMetMht", &gamma_diffMetMht, &b_gamma_diffMetMht);
   fChain->SetBranchAddress("gamma_mht_pt", &gamma_mht_pt, &b_gamma_mht_pt);
   fChain->SetBranchAddress("gamma_mht_phi", &gamma_mht_phi, &b_gamma_mht_phi);
   fChain->SetBranchAddress("gamma_minMTBMet", &gamma_minMTBMet, &b_gamma_minMTBMet);
   fChain->SetBranchAddress("zll_deltaPhiMin", &zll_deltaPhiMin, &b_zll_deltaPhiMin);
   fChain->SetBranchAddress("zll_diffMetMht", &zll_diffMetMht, &b_zll_diffMetMht);
   fChain->SetBranchAddress("zll_mht_pt", &zll_mht_pt, &b_zll_mht_pt);
   fChain->SetBranchAddress("zll_mht_phi", &zll_mht_phi, &b_zll_mht_phi);
   fChain->SetBranchAddress("zll_met_pt", &zll_met_pt, &b_zll_met_pt);
   fChain->SetBranchAddress("zll_met_phi", &zll_met_phi, &b_zll_met_phi);
   fChain->SetBranchAddress("zll_ht", &zll_ht, &b_zll_ht);
   fChain->SetBranchAddress("zll_pt", &zll_pt, &b_zll_pt);
   fChain->SetBranchAddress("zll_eta", &zll_eta, &b_zll_eta);
   fChain->SetBranchAddress("zll_phi", &zll_phi, &b_zll_phi);
   fChain->SetBranchAddress("zll_mass", &zll_mass, &b_zll_mass);
   fChain->SetBranchAddress("pseudoJet2_pt", &pseudoJet2_pt, &b_pseudoJet2_pt);
   fChain->SetBranchAddress("pseudoJet2_eta", &pseudoJet2_eta, &b_pseudoJet2_eta);
   fChain->SetBranchAddress("pseudoJet2_phi", &pseudoJet2_phi, &b_pseudoJet2_phi);
   fChain->SetBranchAddress("pseudoJet2_mass", &pseudoJet2_mass, &b_pseudoJet2_mass);
   fChain->SetBranchAddress("pseudoJet1_pt", &pseudoJet1_pt, &b_pseudoJet1_pt);
   fChain->SetBranchAddress("pseudoJet1_eta", &pseudoJet1_eta, &b_pseudoJet1_eta);
   fChain->SetBranchAddress("pseudoJet1_phi", &pseudoJet1_phi, &b_pseudoJet1_phi);
   fChain->SetBranchAddress("pseudoJet1_mass", &pseudoJet1_mass, &b_pseudoJet1_mass);
   fChain->SetBranchAddress("pseudoViaKtJet1_had_pt", &pseudoViaKtJet1_had_pt, &b_pseudoViaKtJet1_had_pt);
   fChain->SetBranchAddress("pseudoViaKtJet1_had_eta", &pseudoViaKtJet1_had_eta, &b_pseudoViaKtJet1_had_eta);
   fChain->SetBranchAddress("pseudoViaKtJet1_had_phi", &pseudoViaKtJet1_had_phi, &b_pseudoViaKtJet1_had_phi);
   fChain->SetBranchAddress("pseudoViaKtJet1_had_mass", &pseudoViaKtJet1_had_mass, &b_pseudoViaKtJet1_had_mass);
   fChain->SetBranchAddress("gamma_met_pt", &gamma_met_pt, &b_gamma_met_pt);
   fChain->SetBranchAddress("gamma_met_eta", &gamma_met_eta, &b_gamma_met_eta);
   fChain->SetBranchAddress("gamma_met_phi", &gamma_met_phi, &b_gamma_met_phi);
   fChain->SetBranchAddress("gamma_met_mass", &gamma_met_mass, &b_gamma_met_mass);
   fChain->SetBranchAddress("pseudoViaKtJet2_had_pt", &pseudoViaKtJet2_had_pt, &b_pseudoViaKtJet2_had_pt);
   fChain->SetBranchAddress("pseudoViaKtJet2_had_eta", &pseudoViaKtJet2_had_eta, &b_pseudoViaKtJet2_had_eta);
   fChain->SetBranchAddress("pseudoViaKtJet2_had_phi", &pseudoViaKtJet2_had_phi, &b_pseudoViaKtJet2_had_phi);
   fChain->SetBranchAddress("pseudoViaKtJet2_had_mass", &pseudoViaKtJet2_had_mass, &b_pseudoViaKtJet2_had_mass);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
   fChain->SetBranchAddress("pseudoJet1_had_pt", &pseudoJet1_had_pt, &b_pseudoJet1_had_pt);
   fChain->SetBranchAddress("pseudoJet1_had_eta", &pseudoJet1_had_eta, &b_pseudoJet1_had_eta);
   fChain->SetBranchAddress("pseudoJet1_had_phi", &pseudoJet1_had_phi, &b_pseudoJet1_had_phi);
   fChain->SetBranchAddress("pseudoJet1_had_mass", &pseudoJet1_had_mass, &b_pseudoJet1_had_mass);
   fChain->SetBranchAddress("gamma_metNoPU_pt", &gamma_metNoPU_pt, &b_gamma_metNoPU_pt);
   fChain->SetBranchAddress("gamma_metNoPU_eta", &gamma_metNoPU_eta, &b_gamma_metNoPU_eta);
   fChain->SetBranchAddress("gamma_metNoPU_phi", &gamma_metNoPU_phi, &b_gamma_metNoPU_phi);
   fChain->SetBranchAddress("gamma_metNoPU_mass", &gamma_metNoPU_mass, &b_gamma_metNoPU_mass);
   fChain->SetBranchAddress("pseudoJet2_had_pt", &pseudoJet2_had_pt, &b_pseudoJet2_had_pt);
   fChain->SetBranchAddress("pseudoJet2_had_eta", &pseudoJet2_had_eta, &b_pseudoJet2_had_eta);
   fChain->SetBranchAddress("pseudoJet2_had_phi", &pseudoJet2_had_phi, &b_pseudoJet2_had_phi);
   fChain->SetBranchAddress("pseudoJet2_had_mass", &pseudoJet2_had_mass, &b_pseudoJet2_had_mass);
   fChain->SetBranchAddress("gamma_pseudoJet1_pt", &gamma_pseudoJet1_pt, &b_gamma_pseudoJet1_pt);
   fChain->SetBranchAddress("gamma_pseudoJet1_eta", &gamma_pseudoJet1_eta, &b_gamma_pseudoJet1_eta);
   fChain->SetBranchAddress("gamma_pseudoJet1_phi", &gamma_pseudoJet1_phi, &b_gamma_pseudoJet1_phi);
   fChain->SetBranchAddress("gamma_pseudoJet1_mass", &gamma_pseudoJet1_mass, &b_gamma_pseudoJet1_mass);
   fChain->SetBranchAddress("metNoPU_pt", &metNoPU_pt, &b_metNoPU_pt);
   fChain->SetBranchAddress("metNoPU_eta", &metNoPU_eta, &b_metNoPU_eta);
   fChain->SetBranchAddress("metNoPU_phi", &metNoPU_phi, &b_metNoPU_phi);
   fChain->SetBranchAddress("metNoPU_mass", &metNoPU_mass, &b_metNoPU_mass);
   fChain->SetBranchAddress("gamma_pseudoJet2_pt", &gamma_pseudoJet2_pt, &b_gamma_pseudoJet2_pt);
   fChain->SetBranchAddress("gamma_pseudoJet2_eta", &gamma_pseudoJet2_eta, &b_gamma_pseudoJet2_eta);
   fChain->SetBranchAddress("gamma_pseudoJet2_phi", &gamma_pseudoJet2_phi, &b_gamma_pseudoJet2_phi);
   fChain->SetBranchAddress("gamma_pseudoJet2_mass", &gamma_pseudoJet2_mass, &b_gamma_pseudoJet2_mass);
   fChain->SetBranchAddress("nfatJet", &nfatJet, &b_nfatJet);
   fChain->SetBranchAddress("fatJet_pt", fatJet_pt, &b_fatJet_pt);
   fChain->SetBranchAddress("fatJet_eta", fatJet_eta, &b_fatJet_eta);
   fChain->SetBranchAddress("fatJet_phi", fatJet_phi, &b_fatJet_phi);
   fChain->SetBranchAddress("fatJet_mass", fatJet_mass, &b_fatJet_mass);
   fChain->SetBranchAddress("fatJet_btagCSV", fatJet_btagCSV, &b_fatJet_btagCSV);
   fChain->SetBranchAddress("fatJet_rawPt", fatJet_rawPt, &b_fatJet_rawPt);
   fChain->SetBranchAddress("fatJet_mcPt", fatJet_mcPt, &b_fatJet_mcPt);
   fChain->SetBranchAddress("fatJet_mcFlavour", fatJet_mcFlavour, &b_fatJet_mcFlavour);
   fChain->SetBranchAddress("fatJet_prunedMass", fatJet_prunedMass, &b_fatJet_prunedMass);
   fChain->SetBranchAddress("fatJet_trimmedMass", fatJet_trimmedMass, &b_fatJet_trimmedMass);
   fChain->SetBranchAddress("fatJet_tau1", fatJet_tau1, &b_fatJet_tau1);
   fChain->SetBranchAddress("fatJet_tau2", fatJet_tau2, &b_fatJet_tau2);
   fChain->SetBranchAddress("fatJet_tau3", fatJet_tau3, &b_fatJet_tau3);
   if( loadGenStuff ) {
     fChain->SetBranchAddress("ngenPart", &ngenPart, &b_ngenPart);
     fChain->SetBranchAddress("genPart_pt", genPart_pt, &b_genPart_pt);
     fChain->SetBranchAddress("genPart_eta", genPart_eta, &b_genPart_eta);
     fChain->SetBranchAddress("genPart_phi", genPart_phi, &b_genPart_phi);
     fChain->SetBranchAddress("genPart_mass", genPart_mass, &b_genPart_mass);
     fChain->SetBranchAddress("genPart_pdgId", genPart_pdgId, &b_genPart_pdgId);
     fChain->SetBranchAddress("genPart_charge", genPart_charge, &b_genPart_charge);
     fChain->SetBranchAddress("genPart_status", genPart_status, &b_genPart_status);
     fChain->SetBranchAddress("genPart_motherId", genPart_motherId, &b_genPart_motherId);
     fChain->SetBranchAddress("genPart_motherIndex", genPart_motherIndex, &b_genPart_motherIndex);
     fChain->SetBranchAddress("genPart_grandmaId", genPart_grandmaId, &b_genPart_grandmaId);
     fChain->SetBranchAddress("ngenTau", &ngenTau, &b_ngenTau);
     fChain->SetBranchAddress("genTau_pt", genTau_pt, &b_genTau_pt);
     fChain->SetBranchAddress("genTau_eta", genTau_eta, &b_genTau_eta);
     fChain->SetBranchAddress("genTau_phi", genTau_phi, &b_genTau_phi);
     fChain->SetBranchAddress("genTau_mass", genTau_mass, &b_genTau_mass);
     fChain->SetBranchAddress("genTau_pdgId", genTau_pdgId, &b_genTau_pdgId);
     fChain->SetBranchAddress("genTau_charge", genTau_charge, &b_genTau_charge);
     fChain->SetBranchAddress("genTau_status", genTau_status, &b_genTau_status);
     fChain->SetBranchAddress("genTau_sourceId", genTau_sourceId, &b_genTau_sourceId);
     fChain->SetBranchAddress("ngenNu", &ngenNu, &b_ngenNu);
     fChain->SetBranchAddress("genNu_pt", genNu_pt, &b_genNu_pt);
     fChain->SetBranchAddress("genNu_eta", genNu_eta, &b_genNu_eta);
     fChain->SetBranchAddress("genNu_phi", genNu_phi, &b_genNu_phi);
     fChain->SetBranchAddress("genNu_mass", genNu_mass, &b_genNu_mass);
     fChain->SetBranchAddress("genNu_pdgId", genNu_pdgId, &b_genNu_pdgId);
     fChain->SetBranchAddress("genNu_charge", genNu_charge, &b_genNu_charge);
     fChain->SetBranchAddress("genNu_status", genNu_status, &b_genNu_status);
     fChain->SetBranchAddress("genNu_sourceId", genNu_sourceId, &b_genNu_sourceId);
     fChain->SetBranchAddress("ngenLepFromTau", &ngenLepFromTau, &b_ngenLepFromTau);
     fChain->SetBranchAddress("genLepFromTau_pt", genLepFromTau_pt, &b_genLepFromTau_pt);
     fChain->SetBranchAddress("genLepFromTau_eta", genLepFromTau_eta, &b_genLepFromTau_eta);
     fChain->SetBranchAddress("genLepFromTau_phi", genLepFromTau_phi, &b_genLepFromTau_phi);
     fChain->SetBranchAddress("genLepFromTau_mass", genLepFromTau_mass, &b_genLepFromTau_mass);
     fChain->SetBranchAddress("genLepFromTau_pdgId", genLepFromTau_pdgId, &b_genLepFromTau_pdgId);
     fChain->SetBranchAddress("genLepFromTau_charge", genLepFromTau_charge, &b_genLepFromTau_charge);
     fChain->SetBranchAddress("genLepFromTau_status", genLepFromTau_status, &b_genLepFromTau_status);
     fChain->SetBranchAddress("genLepFromTau_sourceId", genLepFromTau_sourceId, &b_genLepFromTau_sourceId);
     fChain->SetBranchAddress("ngenLep", &ngenLep, &b_ngenLep);
     fChain->SetBranchAddress("genLep_pt", genLep_pt, &b_genLep_pt);
     fChain->SetBranchAddress("genLep_eta", genLep_eta, &b_genLep_eta);
     fChain->SetBranchAddress("genLep_phi", genLep_phi, &b_genLep_phi);
     fChain->SetBranchAddress("genLep_mass", genLep_mass, &b_genLep_mass);
     fChain->SetBranchAddress("genLep_pdgId", genLep_pdgId, &b_genLep_pdgId);
     fChain->SetBranchAddress("genLep_charge", genLep_charge, &b_genLep_charge);
     fChain->SetBranchAddress("genLep_status", genLep_status, &b_genLep_status);
     fChain->SetBranchAddress("genLep_sourceId", genLep_sourceId, &b_genLep_sourceId);
     fChain->SetBranchAddress("nGenP6StatusThree", &nGenP6StatusThree, &b_nGenP6StatusThree);
     fChain->SetBranchAddress("GenP6StatusThree_pt", GenP6StatusThree_pt, &b_GenP6StatusThree_pt);
     fChain->SetBranchAddress("GenP6StatusThree_eta", GenP6StatusThree_eta, &b_GenP6StatusThree_eta);
     fChain->SetBranchAddress("GenP6StatusThree_phi", GenP6StatusThree_phi, &b_GenP6StatusThree_phi);
     fChain->SetBranchAddress("GenP6StatusThree_mass", GenP6StatusThree_mass, &b_GenP6StatusThree_mass);
     fChain->SetBranchAddress("GenP6StatusThree_pdgId", GenP6StatusThree_pdgId, &b_GenP6StatusThree_pdgId);
     fChain->SetBranchAddress("GenP6StatusThree_charge", GenP6StatusThree_charge, &b_GenP6StatusThree_charge);
     fChain->SetBranchAddress("GenP6StatusThree_status", GenP6StatusThree_status, &b_GenP6StatusThree_status);
     fChain->SetBranchAddress("GenP6StatusThree_motherId", GenP6StatusThree_motherId, &b_GenP6StatusThree_motherId);
     fChain->SetBranchAddress("GenP6StatusThree_motherIndex", GenP6StatusThree_motherIndex, &b_GenP6StatusThree_motherIndex);
     fChain->SetBranchAddress("GenP6StatusThree_grandmaId", GenP6StatusThree_grandmaId, &b_GenP6StatusThree_grandmaId);
   }
   fChain->SetBranchAddress("nisoTrack", &nisoTrack, &b_nisoTrack);
   fChain->SetBranchAddress("isoTrack_pt", isoTrack_pt, &b_isoTrack_pt);
   fChain->SetBranchAddress("isoTrack_eta", isoTrack_eta, &b_isoTrack_eta);
   fChain->SetBranchAddress("isoTrack_phi", isoTrack_phi, &b_isoTrack_phi);
   fChain->SetBranchAddress("isoTrack_mass", isoTrack_mass, &b_isoTrack_mass);
   fChain->SetBranchAddress("isoTrack_pdgId", isoTrack_pdgId, &b_isoTrack_pdgId);
   fChain->SetBranchAddress("isoTrack_charge", isoTrack_charge, &b_isoTrack_charge);
   fChain->SetBranchAddress("isoTrack_dz", isoTrack_dz, &b_isoTrack_dz);
   fChain->SetBranchAddress("isoTrack_absIso", isoTrack_absIso, &b_isoTrack_absIso);
   fChain->SetBranchAddress("isoTrack_mcMatchId", isoTrack_mcMatchId, &b_isoTrack_mcMatchId);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_mass", jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_btagCSV", jet_btagCSV, &b_jet_btagCSV);
   fChain->SetBranchAddress("jet_rawPt", jet_rawPt, &b_jet_rawPt);
   fChain->SetBranchAddress("jet_mcPt", jet_mcPt, &b_jet_mcPt);
   fChain->SetBranchAddress("jet_mcFlavour", jet_mcFlavour, &b_jet_mcFlavour);
   fChain->SetBranchAddress("jet_mcMatchId", jet_mcMatchId, &b_jet_mcMatchId);
   fChain->SetBranchAddress("jet_mcMatchFlav", jet_mcMatchFlav, &b_jet_mcMatchFlav);
   fChain->SetBranchAddress("jet_area", jet_area, &b_jet_area);
   fChain->SetBranchAddress("jet_puId", jet_puId, &b_jet_puId);
   fChain->SetBranchAddress("jet_id", jet_id, &b_jet_id);
   fChain->SetBranchAddress("jet_ptd", jet_ptd, &b_jet_ptd);
   fChain->SetBranchAddress("jet_axis2", jet_axis2, &b_jet_axis2);
   fChain->SetBranchAddress("jet_mult", jet_mult, &b_jet_mult);
   fChain->SetBranchAddress("jet_qgl", jet_qgl, &b_jet_qgl);
   fChain->SetBranchAddress("jet_partonId", jet_partonId, &b_jet_partonId);
   fChain->SetBranchAddress("jet_partonMotherId", jet_partonMotherId, &b_jet_partonMotherId);
   fChain->SetBranchAddress("nlep", &nlep, &b_nlep);
   fChain->SetBranchAddress("lep_pt", lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_mass", lep_mass, &b_lep_mass);
   fChain->SetBranchAddress("lep_pdgId", lep_pdgId, &b_lep_pdgId);
   fChain->SetBranchAddress("lep_charge", lep_charge, &b_lep_charge);
   fChain->SetBranchAddress("lep_dxy", lep_dxy, &b_lep_dxy);
   fChain->SetBranchAddress("lep_dz", lep_dz, &b_lep_dz);
   fChain->SetBranchAddress("lep_edxy", lep_edxy, &b_lep_edxy);
   fChain->SetBranchAddress("lep_edz", lep_edz, &b_lep_edz);
   fChain->SetBranchAddress("lep_ip3d", lep_ip3d, &b_lep_ip3d);
   fChain->SetBranchAddress("lep_sip3d", lep_sip3d, &b_lep_sip3d);
   fChain->SetBranchAddress("lep_tightId", lep_tightId, &b_lep_tightId);
   fChain->SetBranchAddress("lep_convVeto", lep_convVeto, &b_lep_convVeto);
   fChain->SetBranchAddress("lep_lostHits", lep_lostHits, &b_lep_lostHits);
   fChain->SetBranchAddress("lep_looseIdSusy", lep_looseIdSusy, &b_lep_looseIdSusy);
   fChain->SetBranchAddress("lep_relIso03", lep_relIso03, &b_lep_relIso03);
   fChain->SetBranchAddress("lep_relIso04", lep_relIso04, &b_lep_relIso04);
   fChain->SetBranchAddress("lep_chargedHadRelIso03", lep_chargedHadRelIso03, &b_lep_chargedHadRelIso03);
   fChain->SetBranchAddress("lep_chargedHadRelIso04", lep_chargedHadRelIso04, &b_lep_chargedHadRelIso04);
   fChain->SetBranchAddress("lep_convVetoFull", lep_convVetoFull, &b_lep_convVetoFull);
   fChain->SetBranchAddress("lep_eleCutIdCSA14_25ns_v1", lep_eleCutIdCSA14_25ns_v1, &b_lep_eleCutIdCSA14_25ns_v1);
   fChain->SetBranchAddress("lep_eleCutIdCSA14_50ns_v1", lep_eleCutIdCSA14_50ns_v1, &b_lep_eleCutIdCSA14_50ns_v1);
   fChain->SetBranchAddress("lep_eleMVAId", lep_eleMVAId, &b_lep_eleMVAId);
   fChain->SetBranchAddress("lep_tightCharge", lep_tightCharge, &b_lep_tightCharge);
   fChain->SetBranchAddress("lep_mvaId", lep_mvaId, &b_lep_mvaId);
   fChain->SetBranchAddress("lep_mvaIdTrig", lep_mvaIdTrig, &b_lep_mvaIdTrig);
   fChain->SetBranchAddress("lep_nStations", lep_nStations, &b_lep_nStations);
   fChain->SetBranchAddress("lep_trkKink", lep_trkKink, &b_lep_trkKink);
   fChain->SetBranchAddress("lep_caloCompatibility", lep_caloCompatibility, &b_lep_caloCompatibility);
   fChain->SetBranchAddress("lep_globalTrackChi2", lep_globalTrackChi2, &b_lep_globalTrackChi2);
   fChain->SetBranchAddress("lep_trackerLayers", lep_trackerLayers, &b_lep_trackerLayers);
   fChain->SetBranchAddress("lep_pixelLayers", lep_pixelLayers, &b_lep_pixelLayers);
   fChain->SetBranchAddress("lep_mvaTTH", lep_mvaTTH, &b_lep_mvaTTH);
   fChain->SetBranchAddress("lep_mvaSusy", lep_mvaSusy, &b_lep_mvaSusy);
   fChain->SetBranchAddress("lep_jetPtRatio", lep_jetPtRatio, &b_lep_jetPtRatio);
   fChain->SetBranchAddress("lep_jetBTagCSV", lep_jetBTagCSV, &b_lep_jetBTagCSV);
   fChain->SetBranchAddress("lep_jetDR", lep_jetDR, &b_lep_jetDR);
   fChain->SetBranchAddress("lep_mcMatchId", lep_mcMatchId, &b_lep_mcMatchId);
   fChain->SetBranchAddress("lep_mcMatchAny", lep_mcMatchAny, &b_lep_mcMatchAny);
   fChain->SetBranchAddress("lep_mcMatchTau", lep_mcMatchTau, &b_lep_mcMatchTau);
   fChain->SetBranchAddress("ngamma", &ngamma, &b_ngamma);
   fChain->SetBranchAddress("gamma_pt", gamma_pt, &b_gamma_pt);
   fChain->SetBranchAddress("gamma_eta", gamma_eta, &b_gamma_eta);
   fChain->SetBranchAddress("gamma_phi", gamma_phi, &b_gamma_phi);
   fChain->SetBranchAddress("gamma_mass", gamma_mass, &b_gamma_mass);
   fChain->SetBranchAddress("gamma_pdgId", gamma_pdgId, &b_gamma_pdgId);
   fChain->SetBranchAddress("gamma_idCutBased", gamma_idCutBased, &b_gamma_idCutBased);
   fChain->SetBranchAddress("gamma_hOverE", gamma_hOverE, &b_gamma_hOverE);
   fChain->SetBranchAddress("gamma_r9", gamma_r9, &b_gamma_r9);
   fChain->SetBranchAddress("gamma_sigmaIetaIeta", gamma_sigmaIetaIeta, &b_gamma_sigmaIetaIeta);
   fChain->SetBranchAddress("gamma_chHadIso", gamma_chHadIso, &b_gamma_chHadIso);
   fChain->SetBranchAddress("gamma_chHadIso04", gamma_chHadIso04, &b_gamma_chHadIso04);
   fChain->SetBranchAddress("gamma_chHadIsoRC"  , gamma_chHadIsoRC,   &b_gamma_chHadIsoRC);
   fChain->SetBranchAddress("gamma_chHadIsoRC04", gamma_chHadIsoRC04, &b_gamma_chHadIsoRC04);
   fChain->SetBranchAddress("gamma_neuHadIso", gamma_neuHadIso, &b_gamma_neuHadIso);
   fChain->SetBranchAddress("gamma_phIso", gamma_phIso, &b_gamma_phIso);
   fChain->SetBranchAddress("gamma_mcMatchId", gamma_mcMatchId, &b_gamma_mcMatchId);
   fChain->SetBranchAddress("gamma_genIso", gamma_genIso, &b_gamma_genIso);
   fChain->SetBranchAddress("gamma_genIso04", gamma_genIso04, &b_gamma_genIso04);
   fChain->SetBranchAddress("ntau", &ntau, &b_ntau);
   fChain->SetBranchAddress("tau_pt", tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_eta", tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_mass", tau_mass, &b_tau_mass);
   fChain->SetBranchAddress("tau_pdgId", tau_pdgId, &b_tau_pdgId);
   fChain->SetBranchAddress("tau_charge", tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("tau_dxy", tau_dxy, &b_tau_dxy);
   fChain->SetBranchAddress("tau_dz", tau_dz, &b_tau_dz);
   fChain->SetBranchAddress("tau_idCI3hit", tau_idCI3hit, &b_tau_idCI3hit);
   fChain->SetBranchAddress("tau_isoCI3hit", tau_isoCI3hit, &b_tau_isoCI3hit);
   fChain->SetBranchAddress("tau_mcMatchId", tau_mcMatchId, &b_tau_mcMatchId);
   fChain->SetBranchAddress("evt_scale1fb", &evt_scale1fb, &b_evt_scale1fb);
   fChain->SetBranchAddress("evt_xsec", &evt_xsec, &b_evt_xsec);
   fChain->SetBranchAddress("evt_kfactor", &evt_kfactor, &b_evt_kfactor);
   fChain->SetBranchAddress("evt_filter", &evt_filter, &b_evt_filter);
   fChain->SetBranchAddress("evt_nEvts", &evt_nEvts, &b_evt_nEvts);
   fChain->SetBranchAddress("evt_id", &evt_id, &b_evt_id);
   Notify();
}

Bool_t MT2Tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MT2Tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MT2Tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mt2_cxx
