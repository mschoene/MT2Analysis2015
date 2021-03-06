/*
How to run:
root -l
.L scaleFactors.C+
scaleFactors(inputString, inputFolder, outputFile, treeName, objectName, scaleFileName)
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

using namespace std;


int scaleFactors(string inputString,
		   string inputFolder,
		   string outputFile,
		   string treeName,
		   string objectName,
		   string scaleFileName)
{
   // Add all files in the input folder 
  string dcap = inputFolder.find("pnfs")!=std::string::npos ? "dcap://t3se01.psi.ch:22125/" : "";
  string fullInputString = dcap + inputFolder + "/" + inputString + ".root";
  //std::cout << fullInputString << std::endl;
  TFile *inputFile = TFile::Open(fullInputString.c_str());
  TTree* tree_= (TTree*) inputFile->Get("mt2");
  
  TFile* scaleFile= new TFile(scaleFileName.c_str(), "READ");
  
  TFile *out = TFile::Open(outputFile.c_str(), "RECREATE");
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
  
  int isData; 
  clone->SetBranchAddress("isData",&isData);
  //  clone->GetEntry(0); 
  
  Int_t nobj;
  clone->SetBranchAddress(("n"+objectName).c_str(), &nobj);
  
  Float_t obj_pt[100];
  clone->SetBranchAddress((objectName+"_pt").c_str(), obj_pt);
  Float_t obj_eta[100];
  clone->SetBranchAddress((objectName+"_eta").c_str(), obj_eta);
  
  TTree *ft = new TTree("mt2_sf", "friend tree containing scale factors");
  
  TBranch* b1 = ft->Branch("weight_lepsf", &weight_lepsf, "weight_lepsf/F");
  TBranch* b2 = ft->Branch("weight_lepsf_UP", &weight_lepsf_UP, "weight_lepsf_UP/F");
  TBranch* b3 = ft->Branch("weight_lepsf_DN", &weight_lepsf_DN, "weight_lepsf_DN/F");
  TBranch* b4 = ft->Branch("weight_btagsf", &weight_btagsf, "weight_btagsf/F");
  TBranch* b5 = ft->Branch("weight_btagsf_UP", &weight_btagsf_UP, "weight_btagsf_UP/F");
  TBranch* b6 = ft->Branch("weight_btagsf_DN", &weight_btagsf_DN, "weight_btagsf_DN/F");
  TBranch* b7 = ft->Branch("weight_sigtrigsf", &weight_sigtrigsf, "weight_sigtrigsf/F");
  TBranch* b8 = ft->Branch("weight_dileptrigsf", &weight_dileptrigsf, "weight_dileptrigsf/F");
  TBranch* b9 = ft->Branch("weight_phottrigsf", &weight_phottrigsf, "weight_phottrigsf/F");
  TBranch* b10 =ft->Branch("weight_pu", &weight_pu, "weight_pu/F");
  TBranch* b11 =ft->Branch("weight_isr", &weight_isr, "weight_isr/F");
  TBranch* b12 =ft->Branch("weight_scales_UP", &weight_scales_UP, "weight_scales_UP/F");
  TBranch* b13 =ft->Branch("weight_scales_DN", &weight_scales_DN, "weight_scales_DN/F");
  TBranch* b14 =ft->Branch("weight_pdfs_UP", &weight_pdfs_UP, "weight_pdfs_UP/F");
  TBranch* b15 =ft->Branch("weight_pdfs_DN", &weight_pdfs_DN, "weight_pdfs_DN/F");
  
  TH2F* hS = (TH2F*) scaleFile->Get("h");
  
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
  
    if( objectName == "lep" )
      if(nobj>0)
	for(int o=0; o<nobj; ++o){
	  Int_t binx= hS->GetXaxis()->FindBin(obj_pt[o]);
	  Int_t biny= hS->GetYaxis()->FindBin(fabs(obj_eta[o]));
	  
	  weight_lepsf *= hS->GetBinContent(binx, biny);
	  weight_lepsf_UP *= ( hS->GetBinContent(binx, biny) + hS->GetBinError(binx, biny) );
	  weight_lepsf_DN *= ( hS->GetBinContent(binx, biny) - hS->GetBinError(binx, biny) );

	  if (weight_lepsf < 0.94)
	    std::cout << weight_lepsf << "\t" << obj_pt[o] << "\t"  << obj_eta[o] << std::endl;
	  
	}
    
    ft->Fill();
    
  }
  //-------------------------------------------------------------


  ft->Write();
  delete ft;
  delete clone;
  out->Close();
  return 0;
  
}
