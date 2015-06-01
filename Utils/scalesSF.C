/*
How to run:
root -l
.L scalesSF.C+
scalesSF(inputString, inputFolder, outputFile, treeName, objectName)
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


int scalesSF(string inputString,
	  string inputFolder,
	  string outputFile,
	  string treeName)
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


  Float_t weight_scales_UP;
  Float_t weight_scales_DN;

  int isData; 
  clone->SetBranchAddress("isData",&isData);
  //  clone->GetEntry(0); 
  
  Float_t mt2;
  clone->SetBranchAddress("mt2", &mt2);

  TBranch* b1 = clone->Branch("weight_scales_UP", &weight_scales_UP, "weight_scales_UP/F");
  TBranch* b2 = clone->Branch("weight_scales_DN", &weight_scales_DN, "weight_scales_DN/F");
    
  int nEntries = clone->GetEntries();
  std::cout << "Starting loop over " << nEntries << " entries..." << std::endl;

  for(int i = 0; i < nEntries; ++i) {
    
    if(i%50000 == 0)
      std::cout << "Entry "<< i << "/" << nEntries << std::endl;
    
    clone->GetEntry(i);
    
    weight_scales_UP = 1.;
    weight_scales_DN = 1.;
    
    float weight_scales_err = 0.;
    
    if( isData == 1 ) ;
    else{
      
      if( mt2 < 0. )         continue;
      else if( mt2 < 600.   )  weight_scales_err = 0.05;
      else if( mt2 < 1000.  )  weight_scales_err = 0.10;
      else if( mt2 >= 1000. )  weight_scales_err = 0.30;
          
    }
    
    weight_scales_UP += weight_scales_err;
    weight_scales_DN -= weight_scales_err;

    b1->Fill();
    b2->Fill();
    
  }
  //-------------------------------------------------------------

  clone->Write();
  delete clone;
  out->Close();
  return 0;
  
}
