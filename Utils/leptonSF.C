/*
How to run:
root -l
.L leptonSF.C+
run()
*/

#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

#include <string>
#include <iostream>
#include "TFile.h"
#include "TFileMerger.h"
#include "TFileCollection.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TBranch.h"
#include "TString.h"
#include "TH1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH2.h"
#include "TH3.h"
#include "TLorentzVector.h"

using namespace std;

int setLeptonSF(std::string inputString="input",
		   std::string inputFolder="./",
		   std::string outputFile="output.root",
		   std::string treeName="mt2"){
  
  TH2D* h_id = 0;
  TH2D* h_iso = 0;
  TH2D* h_trk = 0;
  TH2D* h_elSF = 0;
  TH2D* h_elSF_trk = 0;

  TH2D* h_id_mu = 0;
  TH2D* h_iso_mu = 0;
  TH2D* h_dxyz_mu = 0;
  TH2D* h_muSF = 0;

//  TH1D* h_muTrk_hi = 0;
//  TH1D* h_muTrk_lo = 0;
 
//  TH2D* h_fast_muSF = 0;
//  TH2D* h_fast_elSF = 0;

  TH2D* h_eff_full_mu = 0;
  TH2D* h_eff_full_el = 0;

//  TH2D* h_eff_fast_mu =  0;
//  TH2D* h_eff_fast_el =  0;

  //Getting the lepton scale factor histograms/////////////////
  //Electrons//
  std::string filename = "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/scaleFactors.root";
  std::string filenametrkele = "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/egammaEffi.txt_EGM2D.root";
  TFile * f_ele = new TFile(filename.c_str() );
  TFile * f_eletrk = new TFile(filenametrkele.c_str() );
  if (!f_ele->IsOpen() || !f_eletrk->IsOpen()) std::cout << " ERROR: Could not find scale factor file " << filename << std::endl; 
  //Uncomment for loose Id
  //TH2D* h_id = (TH2D*) f_ele->Get("CutBasedLoose");
  //(TH2D*) f_ele->Get("CutBasedVeto");
  
  h_id = (TH2D*) f_ele->Get("GsfElectronToCutBasedSpring15V");
  h_iso = (TH2D*) f_ele->Get("MVAVLooseElectronToMini");
  h_trk = (TH2D*) f_eletrk->Get("EGamma_SF2D");
//  std::cout << h_id->GetBinContent(1,1) << std::endl;
//  std::cout << h_iso->GetBinContent(1,1) << std::endl;
  if (!h_id || !h_iso) std::cout << "ERROR: Could not find scale factor histogram"<< std::endl;
  h_elSF = (TH2D*) h_id->Clone("h_elSF");
  h_elSF->SetDirectory(0);
  h_elSF->Multiply(h_iso);
  h_elSF_trk = (TH2D*) h_trk->Clone("h_elSF_trk");
  h_elSF_trk->SetDirectory(0);
  
  //Muons//
  std::string filenameID =   "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root";
  std::string filenameISO =  "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/TnP_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root";
  std::string filenamedxyz = "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/TnP_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root";
  TFile * f1 = new TFile(filenameID.c_str() );
  TFile * f2 = new TFile(filenameISO.c_str() );
  TFile * f3 = new TFile(filenamedxyz.c_str() );
  if (!f1->IsOpen()) { std::cout<<" ERROR: Could not find ID scale factor file "<<filenameID<<std::endl; return 0;}
  if (!f2->IsOpen()) { std::cout<<"ERROR: Could not find ISO scale factor file "<<filenameISO<<std::endl; return 0;}
  if (!f3->IsOpen()) { std::cout<<"ERROR: Could not find dxy dz scale factor file "<<filenamedxyz<<std::endl; return 0;}
  h_id_mu = (TH2D*) f1->Get("SF");
  h_iso_mu = (TH2D*) f2->Get("SF");
  h_dxyz_mu = (TH2D*) f3->Get("SF");

  if (!h_id_mu || !h_iso_mu  || !h_dxyz_mu) { std::cout<<"ERROR: Could not find scale factor histogram"<<std::endl; return 0;}
  h_muSF = (TH2D*) h_id_mu->Clone("h_muSF");
  h_muSF->SetDirectory(0);
  h_muSF->Multiply(h_iso_mu);
  h_muSF->Multiply(h_dxyz_mu);
  
  f_ele->Close(); f_eletrk->Close();  f1->Close(); f2->Close(); f3->Close();
  delete f_ele; delete f_eletrk; delete f1; delete f2; delete f3;
  
//  TH1D* h_trk_mu_hi = 0;
//  TH1D* h_trk_mu_lo = 0;
//  
//  std::string filenameTrk = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/general_tracks_and_early_general_tracks_corr_ratio.root";
//  TFile * fTrk = new TFile(filenameTrk.c_str() );
//  if (!fTrk->IsOpen()) { std::cout<<" ERROR: Could not find track ineff scale factor file "<<filenameTrk<<std::endl; return 0;}
//  h_trk_mu_hi = (TH1D*) fTrk->Get("mutrksfptg10");
//  if (!h_trk_mu_hi) { std::cout<<"ERROR: Could not find trk sf histogram"<<std::endl; return 0;}
//  h_muTrk_hi = (TH1D*) h_trk_mu_hi->Clone("h_muTrk_hi");
//  h_muTrk_hi->SetDirectory(0);
//  h_trk_mu_lo = (TH1D*) fTrk->Get("mutrksfptl10");
//  if (!h_trk_mu_lo) { std::cout<<"ERROR: Could not find trk sf histogram"<<std::endl; return 0;}
//  h_muTrk_lo = (TH1D*) h_trk_mu_lo->Clone("h_muTrk_lo");
//  h_muTrk_lo->SetDirectory(0);
//  fTrk->Close(); delete fTrk;
  
  
  std::cout << std::endl;
  std::cout << "Using Loose Muon ID, MiniIso 0.2 lepton scale factors" << std::endl;
  std::cout << "Using Veto Electrons ID, MiniIso 0.1 lepton scale factors" << std::endl;
  std::cout << "Be aware that Veto Electrons are not suited for selecting Electrons." << std::endl;
  std::cout << std::endl;
  
  std::cout << "Also loading the FullSim efficiency map" << std::endl;
  TFile * f_eff_full = new TFile("/mnt/t3nfs01/data01/shome/mmasciov/lepSF/vetoeff_emu_etapt_lostlep.root" );
  if(!f_eff_full->IsOpen()) {std::cout<<" ERROR: Could not find muon Fullsim scale factor file" <<std::endl; return 0;}
  h_eff_full_mu = (TH2D*) f_eff_full->Get("h_mu_comb_eff");
  h_eff_full_el = (TH2D*) f_eff_full->Get("h_ele_comb_eff");
  if(!h_eff_full_mu || !h_eff_full_el ) {std::cout << " ERROR: Could not find the 2D histogram in your files " << std::endl; return 0;}
  h_eff_full_mu->SetDirectory(0);
  h_eff_full_el->SetDirectory(0);
  f_eff_full->Close();
  delete f_eff_full;
 

//  // Add all files in the input folder
//  string dcap = inputFolder.find("pnfs")!=std::string::npos ? "dcap://t3se01.psi.ch:22125/" : "";
//  string fullInputString = dcap + inputFolder + "/" + inputString + ".root";
//  //std::cout << fullInputString << std::endl;                                                                                                                                                                                                                               
  string fullInputString = inputFolder + "/" + inputString + ".root";
  string objectName = "lep";
 
  TFile *inputFile = TFile::Open(fullInputString.c_str());
  TTree* tree_= (TTree*) inputFile->Get("mt2");

  TFile *out = new TFile(outputFile.c_str(), "RECREATE");
  TTree *clone = new TTree("mt2", "post processed baby tree for mt2 analysis");

  clone = tree_->CloneTree(-1, "fast");
  clone->SetName("mt2");

  Float_t weight_lepsf;
  Float_t weight_lepsf_UP;
  Float_t weight_lepsf_DN;

  Float_t weight_lepsf_0l;
  Float_t weight_lepsf_0l_UP;
  Float_t weight_lepsf_0l_DN;


  ////// Lepton Efficiency SF (to be run only for MC)
  Int_t nlep;
  Float_t lep_pt[100];
  Float_t lep_eta[100];
  Int_t lep_pdgId[100];

  int isData;
  clone->SetBranchAddress("isData",&isData);
  clone->SetBranchAddress("nlep", &nlep);
  clone->SetBranchAddress("lep_pt", lep_pt);
  clone->SetBranchAddress("lep_eta", lep_eta);    
  clone->SetBranchAddress("lep_pdgId", lep_pdgId);

  Int_t ngenLep;
  Float_t genLep_pt[100];
  Float_t genLep_eta[100];  
  Int_t genLep_pdgId[100];
  Int_t ngenLepFromTau;
  Float_t genLepFromTau_pt[100];
  Float_t genLepFromTau_eta[100];  
  Int_t genLepFromTau_pdgId[100];
  Int_t nMuons10;
  Int_t nElectrons10;
  Int_t nPFLep5LowMT;
  Int_t nPFHad10LowMT;

  clone->SetBranchAddress("ngenLep", &ngenLep);
  clone->SetBranchAddress("genLep_pt", genLep_pt);
  clone->SetBranchAddress("genLep_eta", genLep_eta);
  clone->SetBranchAddress("genLep_pdgId", genLep_pdgId);
  
  clone->SetBranchAddress("ngenLepFromTau", &ngenLepFromTau);
  clone->SetBranchAddress("genLepFromTau_pt", genLepFromTau_pt);
  clone->SetBranchAddress("genLepFromTau_eta", genLepFromTau_eta);
  clone->SetBranchAddress("genLepFromTau_pdgId", genLepFromTau_pdgId);
  
  clone->SetBranchAddress("nMuons10", &nMuons10);
  clone->SetBranchAddress("nElectrons10", &nElectrons10);
  clone->SetBranchAddress("nPFLep5LowMT", &nPFLep5LowMT);
  clone->SetBranchAddress("nPFHad10LowMT", &nPFHad10LowMT);

 
  TBranch* b1 = 0; 
  TBranch* b2 = 0; 
  TBranch* b3 = 0; 

  TBranch* b4 = 0;
  TBranch* b5 = 0; 
  TBranch* b6 = 0; 

  b1 = clone->Branch("weight_lepsf2017", &weight_lepsf, "weight_lepsf/F");
  b2 = clone->Branch("weight_lepsf2017_UP", &weight_lepsf_UP, "weight_lepsf_UP/F");
  b3 = clone->Branch("weight_lepsf2017_DN", &weight_lepsf_DN, "weight_lepsf_DN/F");

  b4 = clone->Branch("weight_lepsf2017_0l", &weight_lepsf_0l, "weight_lepsf_0l/F");
  b5 = clone->Branch("weight_lepsf2017_0l_UP", &weight_lepsf_0l_UP, "weight_lepsf_0l_UP/F");
  b6 = clone->Branch("weight_lepsf2017_0l_DN", &weight_lepsf_0l_DN, "weight_lepsf_0l_DN/F");


  std::cout << "Entering the final loop over the events" << std::endl;
  
  Long64_t nEntries = clone->GetEntries();

  for( Long64_t i = 0; i < nEntries; i++) {
  
    if( i==0)
      std::cout << "Start of loop over tree entries" << std::endl;


    weight_lepsf = 1.;
    weight_lepsf_UP = 1.;
    weight_lepsf_DN = 1.;

    weight_lepsf_0l = 1;
    weight_lepsf_0l_UP = 1;
    weight_lepsf_0l_DN = 1;
    
    clone->GetEntry(i);

    /////////Add lepton scale factor//////////
    if(nlep>0){
      Float_t uncert = 1; //Place holder for total uncertainty
      Float_t central = 1; 
      Float_t err = 0;
      Float_t uncert_UP = 0; 	Float_t uncert_DN = 0; 
      
      Float_t fast_central = 1; 
      Float_t fast_err = 0;
      
      
      for(int o=0; o<nlep; ++o){
	
	float pt = lep_pt[o];
	float eta = fabs( lep_eta[o] );
	int pdgId = abs( lep_pdgId[o] );
	
	float pt_cutoff = std::max( 10.1, std::min( 100., double(pt) ) );
	float pt_cutoff_eff = std::max( 5.1, std::min( 100., double(pt) ) );
	
	//Electrons
	if( pdgId == 11) {
	  Int_t binx = h_elSF->GetXaxis()->FindBin(pt_cutoff);
	  Int_t biny = h_elSF->GetYaxis()->FindBin(eta);
	  central = h_elSF->GetBinContent(binx,biny);
	  err  = h_elSF->GetBinError(binx,biny);
	  
	  int binx_trk = h_elSF_trk->GetXaxis()->FindBin(eta);
	  int biny_trk = 1; // hardcoding for now - only one bin in pt (hist starts at 20)
	  central *= h_elSF_trk->GetBinContent(binx_trk,biny_trk);
	  float trk_err = h_elSF_trk->GetBinError(binx_trk,biny_trk);
	  if (pt_cutoff < 20. || pt_cutoff > 80.) err = sqrt(err*err + trk_err*trk_err + 0.01*0.01); 
	  else err = sqrt(err*err + trk_err*trk_err);
	  
	  if (central > 1.3 || central < 0.3) 
	    std::cout<<"STRANGE: Electron with pT/eta of "<< pt <<"/"<< eta <<". SF is "<< central <<std::endl;
	  uncert_UP = central + err;
	  uncert_DN = central - err;
	  	  
	} //else Muons
	else if (abs( lep_pdgId[o]) == 13) {
	  Int_t binx = h_muSF->GetXaxis()->FindBin(pt_cutoff);
	  Int_t biny = h_muSF->GetYaxis()->FindBin(fabs(eta));
	  
//	  float central_trk = 1;
//	  Int_t binx_trk = h_muTrk_hi->GetXaxis()->FindBin(  lep_eta[o] );
//	  if( binx_trk>10 ) binx_trk = 10;
//	  else if( binx_trk<1 ) binx_trk = 1;
//	  central_trk = h_muTrk_hi->GetBinContent( binx_trk );
//	  
//	  
//	  if ( binx >7 ) binx = 7; //overflow bin empty for the muons...
//	  central = h_muSF->GetBinContent(binx,biny);
//	  
//	  central *= central_trk;
	  

	  central = h_muSF->GetBinContent(binx,biny);
	  err  = 0.03; //current recomendation is 3% //   err  = 0.014; // adding in quadrature 1% unc. on ID and 1% unc. on ISO
	  if (central > 1.2 || central < 0.8) 
	    std::cout<<"STRANGE: Muon with pT/eta of "<<pt<<"/"<< fabs(eta) <<". SF is "<< central <<std::endl;
	  uncert_UP = central + err;
	  uncert_DN = central - err;
	  
	}//done with one  electron/muon 
	weight_lepsf    *= central;
	weight_lepsf_UP *= uncert_UP;
	weight_lepsf_DN *= uncert_DN;	  	
      }//end of loop over objects
      
    }//end of lepton sf




    //0lepton efficiency correction
    if( nMuons10+nElectrons10+nPFLep5LowMT+nPFHad10LowMT == 0 ){ //otherwise we have caught the lep


      for(int o=0; o < ngenLep+ngenLepFromTau; ++o){

	float central=1.;
	float err=0;
	float uncert=0;
	float fast_central=1.;
	float fast_err=0.;
	float uncert_fast=0;
	
	
	float pt, eta;
	int  pdgId;
	if(o < ngenLep){
	  pt = genLep_pt[o];
	  eta = genLep_eta[o];
	  pdgId = genLep_pdgId[o];
	}else{
	  pt = genLepFromTau_pt[o-ngenLep];
	  eta = genLepFromTau_eta[o-ngenLep];
	  pdgId = genLepFromTau_pdgId[o-ngenLep] ;
	}
	// acceptance of lepton veto: pt>5, |eta| < 2.4
	if( pt < 5.) continue;
	if( fabs(eta) > 2.4) continue; 
	
	
	float pt_cutoff = std::max( 10.1, std::min( 100., double(pt) ) );
	float pt_cutoff_eff = std::max( 5.1, std::min( 100., double(pt) ) );
	float veto_eff = 1.; //just need one veto eff, not one for full and fast sim
	
	if ( abs( pdgId ) == 11) { 
	  Int_t binx = h_elSF->GetXaxis()->FindBin(pt_cutoff);
	  Int_t biny = h_elSF->GetYaxis()->FindBin( fabs(eta) );
	  central = h_elSF->GetBinContent(binx,biny);
	  uncert  = h_elSF->GetBinError(binx,biny);
	  
	  Int_t binx_eff = h_eff_full_el->GetXaxis()->FindBin(pt_cutoff_eff);
	  Int_t biny_eff = h_eff_full_el->GetYaxis()->FindBin( fabs(eta) );
	  veto_eff = h_eff_full_el->GetBinContent(binx_eff,biny_eff);
	  
	} else if( abs( pdgId ) == 13) {
	  Int_t binx = h_muSF->GetXaxis()->FindBin(pt_cutoff);
	  Int_t biny = h_muSF->GetYaxis()->FindBin( fabs(eta) );
	  central = h_muSF->GetBinContent(binx,biny);
	  uncert = 0.03;// uncert  = 0.014;   // adding in quadrature 1% unc. on ID and 1% unc. on ISO
	  
//	  float central_trk = 1;
//	  if( pt >= 10.){
//	    Int_t binx_trk = h_muTrk_hi->GetXaxis()->FindBin(  lep_eta[o] );
//	    if( binx_trk>10 ) binx_trk = 10;
//	    else if( binx_trk<1 ) binx_trk = 1;
//	    central_trk = h_muTrk_hi->GetBinContent( binx_trk );
//	  }else if( pt < 10 ){
//	    Int_t binx_trk = h_muTrk_lo->GetXaxis()->FindBin(  lep_eta[o] );
//	    if( binx_trk>10 ) binx_trk = 10;
//	    else if( binx_trk<1 ) binx_trk = 1;
//	    central_trk = h_muTrk_lo->GetBinContent( binx_trk );
//	  }
//	  central *= central_trk;
	  
	  Int_t binx_eff = h_eff_full_mu->GetXaxis()->FindBin( pt_cutoff_eff );
	  Int_t biny_eff = h_eff_full_mu->GetYaxis()->FindBin( fabs(eta) );
	  veto_eff = h_eff_full_mu->GetBinContent(binx_eff,biny_eff);
	}
	
	//Only full sim correction, only uncertainty, not weight
	    
	float sf = central;
	// float veto_eff_corr = veto_eff* sf; //corrected veto eff with the 
	
	float unc = uncert;
	float veto_eff_unc_UP = veto_eff* (1. + unc);
	
	float unc_UP_0l = (( 1. - veto_eff_unc_UP) / (1. - veto_eff))  -1.;
	weight_lepsf_0l_UP *= ( 1. + unc_UP_0l);
	weight_lepsf_0l_DN *= ( 1. - unc_UP_0l);

      }//end loop over genLep
      
    }//end 0Lep region
    

    if( i==nEntries-1)
      std::cout << "End of loop over tree entries" << std::endl;
    
    b1->Fill();
    b2->Fill();
    b3->Fill();
    b4->Fill();
    b5->Fill();
    b6->Fill();
    

  }//end loop over events
  //-------------------------------------------------------------


  clone->Write();
  delete clone;
  out->Close();
  delete out;
  return 0;
  
}


