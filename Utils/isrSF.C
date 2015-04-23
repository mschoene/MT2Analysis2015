/*
How to run:
root -l
.L isrSF.C+
isrSF(inputString, inputFolder, outputFile, treeName, objectName)
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


int isrSF(string inputString,
	  string inputFolder,
	  string outputFile,
	  string treeName,
	  string objectName)
{
   // Add all files in the input folder 
  string dcap = inputFolder.find("pnfs")!=std::string::npos ? "dcap://t3se01.psi.ch:22125/" : "";
  string fullInputString = dcap + inputFolder + "/" + inputString + ".root";
  //std::cout << fullInputString << std::endl;
  TFile *inputFile = TFile::Open(fullInputString.c_str());
  TTree* tree_= (TTree*) inputFile->Get("mt2");
    
  TFile *out = new TFile(outputFile.c_str(), "RECREATE");
  TTree *clone = new TTree("mt2", "post processed baby tree for mt2 analysis");

  clone = tree_->CloneTree(-1, "fast"); 
  clone->SetName("mt2");

//  Float_t weight_lepsf;
//  Float_t weight_lepsf_UP;
//  Float_t weight_lepsf_DN;
//  Float_t weight_btagsf;
//  Float_t weight_btagsf_UP;
//  Float_t weight_btagsf_DN;
//  Float_t weight_sigtrigsf;
//  Float_t weight_dileptrigsf;
//  Float_t weight_phottrigsf;
//  Float_t weight_pu;
//  Float_t weight_isr;
//  Float_t weight_scales_UP;
//  Float_t weight_scales_DN;
//  Float_t weight_pdfs_UP;
//  Float_t weight_pdfs_DN;

  Float_t weight_isr;

  int isData; 
  clone->SetBranchAddress("isData",&isData);
  //  clone->GetEntry(0); 
  
  Int_t nobj;
  clone->SetBranchAddress(("n"+objectName).c_str(), &nobj);
  
  Float_t obj_pt[100];
  clone->SetBranchAddress((objectName+"_pt").c_str(), obj_pt);
  Float_t obj_eta[100];
  clone->SetBranchAddress((objectName+"_eta").c_str(), obj_eta);
  Float_t obj_phi[100];
  clone->SetBranchAddress((objectName+"_phi").c_str(), obj_phi);
  Float_t obj_mass[100];
  clone->SetBranchAddress((objectName+"_mass").c_str(), obj_mass);
  Int_t obj_status[100];
  clone->SetBranchAddress((objectName+"_status").c_str(), obj_status);

  TBranch* b1 = clone->Branch("weight_isr", &weight_isr, "weight_isr/F");
    
  int nEntries = clone->GetEntries();
  std::cout << "Starting loop over " << nEntries << " entries..." << std::endl;

  for(int i = 0; i < nEntries; ++i) {
    
    if(i%50000 == 0)
      std::cout << "Entry "<< i << "/" << nEntries << std::endl;
    
    clone->GetEntry(i);
    
    weight_isr=1.;
    
    if( objectName != "GenPart" || isData == 1 ) ;
    else{
      
      TLorentzVector s;
      
      if(nobj>0)
	for(int o=0; o<nobj; ++o){
	  
	  if(obj_status[o] != 62) continue;
	  
	  TLorentzVector s_;
	  s_.SetPtEtaPhiM(obj_pt[o], obj_eta[o], obj_phi[o], obj_mass[o]);
	  s+=s_;
	  
	}
      
      float pt_hard = s.Pt();
      if( pt_hard < 0. )         continue;
      else if( pt_hard < 120. )  weight_isr = 1.0;
      else if( pt_hard < 150. )  weight_isr = 0.95;
      else if( pt_hard < 250. )  weight_isr = 0.90;
      else if( pt_hard >= 250. ) weight_isr = 0.80;
      
    }
    
    b1->Fill();
    
  }
  //-------------------------------------------------------------


  clone->Write();
  delete clone;
  out->Close();
  return 0;
  
}
