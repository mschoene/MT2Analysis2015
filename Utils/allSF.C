/*
How to run:
root -l
.L allSF.C+
allSF(inputString, inputFolder, outputFile, treeName, objectName, scaleFileName)
*/

#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TFile.h"
#include "TFileMerger.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TBranch.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLorentzVector.h"

using namespace std;

void get_SF_btag(float obj_pt, float obj_eta, int obj_mcFlavour, float obj_btagCSV, float &SF, float &SFerr){
  float x = obj_pt; ///the pt of the jet
  float eta = fabs(obj_eta); ///abs(eta)

  if(eta > 2.5){
    std::cout << "warning SF_btag_eta>2.5? " << eta << std::endl;
    exit(1);
  }

  if(abs(obj_mcFlavour)==5 || abs(obj_mcFlavour) == 4){ //for b or c. jet_mcFlavour[indj] refers to the MC-true flavor of the jet

    SF = 0.95;
    if(x < 30. || x > 250.)
      SFerr = 0.05;
    else
      SFerr = 0.02;

    if( abs(obj_mcFlavour) == 4 ) SFerr *= 2;

  }
  else{ ///SFlight

    SF = 1.05;
    SFerr = 0.10;

  }

}


void get_weight_btag(int nobj, float* obj_pt, float* obj_eta, int* obj_mcFlavour, float* obj_btagCSV, float &wtbtag, float &wtbtagErr){

  float mcTag = 1.;
  float mcNoTag = 1.;
  float dataTag = 1.;
  float dataNoTag = 1.;
  float errTag = 0.;
  float errNoTag = 0.;

  float err1 = 0;
  float err2 = 0;
  float err3 = 0;
  float err4 = 0;

  for(int indj=0; indj < nobj; ++indj){ //Here we loop over all selected jets ( for our case, pt>30, PF loose ID, etc )

    float csv = obj_btagCSV[indj]; ////here we get the CSV btag value
    int partonFlavor = abs(obj_mcFlavour[indj]);
    float eta = fabs(obj_eta[indj]);

    if(eta > 2.5) continue;
    if(partonFlavor==0) continue; //for jets with flavor 0, we ignore.

    float eff;
    if( partonFlavor==5 ) {
      ///here one need to provide the pt/eta dependent efficiency for b-tag for "b jet"
      eff = 0.70;
    }
    else if( partonFlavor==4){
      ///here one need to provide the pt/eta dependent efficiency for b-tag for "c jet"
      eff = 0.17;
    }
    else{
      ///here one need to provide the pt/eta dependent efficiency for b-tag for "light jet"
      eff = 0.015;
    }

    bool istag = csv > 0.814 && eta < 2.5 ;
    float SF = 0.;
    float SFerr = 0.;
    get_SF_btag(obj_pt[indj], obj_eta[indj], obj_mcFlavour[indj], obj_btagCSV[indj], SF, SFerr);

    if(istag){
      mcTag *= eff;
      dataTag *= eff*SF;

      if(partonFlavor==5 || partonFlavor ==4)  err1 += SFerr/SF; ///correlated for b/c
      else err3 += SFerr/SF; //correlated for light

    }
    else{
      mcNoTag *= (1 - eff);
      dataNoTag *= (1 - eff*SF);

      if( partonFlavor==5 || partonFlavor==4 ) err2 += (-eff*SFerr)/(1-eff*SF); /// /correlated for b/c
      else err4 +=  (-eff*SFerr)/(1-eff*SF);  ////correlated for light

    }

  }

  wtbtag = (dataNoTag * dataTag ) / ( mcNoTag * mcTag );
  wtbtagErr = sqrt( pow(err1+err2,2) + pow( err3 + err4,2)) * wtbtag;  ///un-correlated for b/c and light

}



int allSF(string inputString,
	  string inputFolder,
	  string outputFile,
	  string treeName,
	  string scaleFileName)
{
  // Add all files in the input folder 
  string dcap = inputFolder.find("pnfs")!=std::string::npos ? "dcap://t3se01.psi.ch:22125/" : "";
  string fullInputString = dcap + inputFolder + "/" + inputString + ".root";
  //std::cout << fullInputString << std::endl;
  TFile *inputFile = TFile::Open(fullInputString.c_str());
  TTree* tree_= (TTree*) inputFile->Get("mt2");
  
  TFile* scaleFile= new TFile(scaleFileName.c_str(), "READ");
  
  TFile *out = new TFile(outputFile.c_str(), "RECREATE");
  TTree *clone = new TTree("mt2", "post processed baby tree for mt2 analysis");

  clone = tree_->CloneTree(-1, "fast"); 
  clone->SetName("mt2");

  Float_t weight_lepsf;
  Float_t weight_lepsf_UP;
  Float_t weight_lepsf_DN;
  Float_t weight_btagsf;
  Float_t weight_btagsf_UP;
  Float_t weight_btagsf_DN;
  Float_t weight_sigtrigsf;
  Float_t weight_dileptrigsf;
  Float_t weight_phottrigsf;
  Float_t weight_pu;
  Float_t weight_isr;
  Float_t weight_scales_UP;
  Float_t weight_scales_DN;
  Float_t weight_pdfs_UP;
  Float_t weight_pdfs_DN;

  Float_t weight_btagsf_err;
  Float_t weight_scales_err;

  int isData; 
  clone->SetBranchAddress("isData",&isData);
  //  clone->GetEntry(0); 
  
  ////// Lepton Efficiency SF
  Int_t nlep;
  clone->SetBranchAddress("nlep", &nlep);
  Float_t lep_pt[100];
  clone->SetBranchAddress("lep_pt", lep_pt);
  Float_t lep_eta[100];
  clone->SetBranchAddress("lep_eta", lep_eta);
  
  ////// b-tag SF
  Int_t njet;
  clone->SetBranchAddress("njet", &njet);
  Float_t jet_pt[100];
  clone->SetBranchAddress("jet_pt", jet_pt);
  Float_t jet_eta[100];
  clone->SetBranchAddress("jet_eta", jet_eta);
  Int_t jet_mcFlavour[100];
  clone->SetBranchAddress("jet_mcFlavour", jet_mcFlavour);
  Float_t jet_btagCSV[100];
  clone->SetBranchAddress("jet_btagCSV", jet_btagCSV);

  ////// isr re-weight
  Int_t nGenPart;
  clone->SetBranchAddress("nGenPart", &nGenPart);
  Float_t GenPart_pt[100];
  clone->SetBranchAddress("GenPart_pt", GenPart_pt);
  Float_t GenPart_eta[100];
  clone->SetBranchAddress("GenPart_eta", GenPart_eta);
  Float_t GenPart_phi[100];
  clone->SetBranchAddress("GenPart_phi", GenPart_phi);
  Float_t GenPart_mass[100];
  clone->SetBranchAddress("GenPart_mass", GenPart_mass);
  Int_t GenPart_status[100];
  clone->SetBranchAddress("GenPart_status", GenPart_status);

  ////// normalization & factorization scales
  Float_t mt2;
  clone->SetBranchAddress("mt2", &mt2);

  ////// Initializing all branches to be filled
  TBranch* b1 = clone->Branch("weight_lepsf", &weight_lepsf, "weight_lepsf/F");
  TBranch* b2 = clone->Branch("weight_lepsf_UP", &weight_lepsf_UP, "weight_lepsf_UP/F");
  TBranch* b3 = clone->Branch("weight_lepsf_DN", &weight_lepsf_DN, "weight_lepsf_DN/F");
  TBranch* b4 = clone->Branch("weight_btagsf", &weight_btagsf, "weight_btagsf/F");
  TBranch* b5 = clone->Branch("weight_btagsf_UP", &weight_btagsf_UP, "weight_btagsf_UP/F");
  TBranch* b6 = clone->Branch("weight_btagsf_DN", &weight_btagsf_DN, "weight_btagsf_DN/F");
  TBranch* b7 = clone->Branch("weight_sigtrigsf", &weight_sigtrigsf, "weight_sigtrigsf/F");
  TBranch* b8 = clone->Branch("weight_dileptrigsf", &weight_dileptrigsf, "weight_dileptrigsf/F");
  TBranch* b9 = clone->Branch("weight_phottrigsf", &weight_phottrigsf, "weight_phottrigsf/F");
  TBranch* b10 = clone->Branch("weight_pu", &weight_pu, "weight_pu/F");
  TBranch* b11 = clone->Branch("weight_isr", &weight_isr, "weight_isr/F");
  TBranch* b12 = clone->Branch("weight_scales_UP", &weight_scales_UP, "weight_scales_UP/F");
  TBranch* b13 = clone->Branch("weight_scales_DN", &weight_scales_DN, "weight_scales_DN/F");
  TBranch* b14 = clone->Branch("weight_pdfs_UP", &weight_pdfs_UP, "weight_pdfs_UP/F");
  TBranch* b15 = clone->Branch("weight_pdfs_DN", &weight_pdfs_DN, "weight_pdfs_DN/F");
  
  ////// Reading histogram for lepton efficiency, from Tag & Probe
  TH2F* hS = (TH2F*) scaleFile->Get("h");
  

  ////// Starting loop over all entries
  int nEntries = clone->GetEntries();
  for(int i = 0; i < nEntries; ++i) {
    
    clone->GetEntry(i);
    
    weight_lepsf=1.;
    weight_lepsf_UP=1.;
    weight_lepsf_DN=1.;
    weight_btagsf=1.;
    weight_btagsf_UP=1.;
    weight_btagsf_DN=1.;
    weight_sigtrigsf=1.;
    weight_dileptrigsf=1.;
    weight_phottrigsf=1.;
    weight_pu=1.;
    weight_isr=1.;
    weight_scales_UP=1.;
    weight_scales_DN=1.;
    weight_pdfs_UP=1.;
    weight_pdfs_DN=1.;
    
    if( isData == 1 ) ;
    else{

      ////// Lepton Efficiency SF
      for(int l=0; l < nlep; ++l){

	Int_t binx= hS->GetXaxis()->FindBin(lep_pt[l]);
	Int_t biny= hS->GetYaxis()->FindBin(fabs(lep_eta[l]));
	
	weight_lepsf *= hS->GetBinContent(binx, biny);
	weight_lepsf_UP *= ( hS->GetBinContent(binx, biny) + hS->GetBinError(binx, biny) );
	weight_lepsf_DN *= ( hS->GetBinContent(binx, biny) - hS->GetBinError(binx, biny) );

      }
      
      ////// b-tag SF
      get_weight_btag(njet, jet_pt, jet_eta, jet_mcFlavour, jet_btagCSV, weight_btagsf, weight_btagsf_err);
      weight_btagsf_UP = weight_btagsf + weight_btagsf_err;
      weight_btagsf_DN = weight_btagsf - weight_btagsf_err;
      
      ////// isr re-weight
      TLorentzVector s;

      for(int o=0; o < nGenPart; ++o){
	
	if(GenPart_status[o] != 62) continue;
	
	TLorentzVector s_;
	s_.SetPtEtaPhiM(GenPart_pt[o], GenPart_eta[o], GenPart_phi[o], GenPart_mass[o]);
	s+=s_;
	
      }
      float pt_hard = s.Pt();
      if( pt_hard < 0. )         continue;
      else if( pt_hard < 120. )  weight_isr = 1.0;
      else if( pt_hard < 150. )  weight_isr = 0.95;
      else if( pt_hard < 250. )  weight_isr = 0.90;
      else if( pt_hard >= 250. ) weight_isr = 0.80;
      
      ////// normalization & factorization scales
      if( mt2 < 0. )         continue;
      else if( mt2 < 600.   )  weight_scales_err = 0.05;
      else if( mt2 < 1000.  )  weight_scales_err = 0.10;
      else if( mt2 >= 1000. )  weight_scales_err = 0.30;
      weight_scales_UP += weight_scales_err;
      weight_scales_DN -= weight_scales_err;

    }
    
    b1->Fill();
    b2->Fill();
    b3->Fill();
    b4->Fill();
    b5->Fill();
    b6->Fill();
    b7->Fill();
    b8->Fill();
    b9->Fill();
    b10->Fill();
    b11->Fill();
    b12->Fill();
    b13->Fill();
    b14->Fill();
    b15->Fill();
   
  }
  //-------------------------------------------------------------


  clone->Write();
  delete clone;
  out->Close();
  return 0;
  
}
