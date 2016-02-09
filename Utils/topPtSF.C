/*
How to run:
root -l
.L topPtSF.C+
topPtSF(inputString, inputFolder, outputFile, treeName, objectName)
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


int topPtSF(string inputString,
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

  clone = tree_->CloneTree(-1, "mt2"); 
  clone->SetName("mt2");

  //Values from Run1 still valid for now, see here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  float a = 0.156;
  float b = -0.00137;


  Float_t weight_topPt;
  Float_t weight_topPt_uncert;

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
  Int_t obj_pdgId[100];
  clone->SetBranchAddress((objectName+"_pdgId").c_str(), obj_pdgId);


  TBranch* b1 = clone->Branch("weight_topPt", &weight_topPt, "weight_topPt/F");
  TBranch* b2 = clone->Branch("weight_topPt_uncert", &weight_topPt_uncert, "weight_topPt_uncert/F");
    
  int nEntries = clone->GetEntries();
  std::cout << "Starting loop over " << nEntries << " entries..." << std::endl;

  //top pt scale factors should not change the overall normalization
  //first loop get the average of the scale factors
  double average = 0;

  for(int i = 0; i < nEntries; ++i) {
    
    if(i%50000 == 0)
      std::cout << "Entry "<< i << "/" << nEntries << std::endl;
    
    clone->GetEntry(i);
    weight_topPt = 1.;

    if( objectName != "GenPart" || isData == 1 ) ;
    else{
      if(nobj>0)
	for(int o=0; o<nobj; ++o){
   
	  if(obj_status[o] != 62) continue;
    
	  //Only apply weight to top or antitop
	  if( abs(obj_pdgId[o])==6 )
	    weight_topPt *= exp( a + b* obj_pt[o] );
	  else continue;

	}//end loop over objects
    }
    
    average += weight_topPt;   
 
  }
  //-------------------------------------------------------------end of run over events

  average /=  double(nEntries);
  std::cout << "average = " << average << std::endl;

  double average_test = 0;

  for(int i = 0; i < nEntries; ++i) {
    
    if(i%50000 == 0)
      std::cout << "Entry "<< i << "/" << nEntries << std::endl;
    
    clone->GetEntry(i);
    
    weight_topPt = 1.;
    weight_topPt_uncert = 0;

    if( objectName != "GenPart" || isData == 1 ) ;
    else{
      if(nobj>0)
	for(int o=0; o<nobj; ++o){
       
	  if(obj_status[o] != 62) continue;

	  //Only apply weight to top or antitop
	  if( abs(obj_pdgId[o])==6 )
	    weight_topPt *= exp( a + b * obj_pt[o] );
	  else continue;

	}//end loop over objects
    
      weight_topPt /= average;

    }

    b1->Fill();  
    
  }
  //-------------------------------------------------------------end of run over events


  clone->Write();
  delete clone;
  out->Close();
  return 0;
  
}
