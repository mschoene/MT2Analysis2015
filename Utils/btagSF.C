/*
How to run:
root -l
.L btagSF.C+
btagSF(inputString, inputFolder, outputFile, treeName, objectName)
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

int btagSF(string inputString,
	   string inputFolder,
	   string outputFile,
	   string treeName,
	   string objectName)
{
  //Add all files in the input folder 
  string dcap = inputFolder.find("pnfs")!=std::string::npos ? "dcap://t3se01.psi.ch:22125/" : "";
  string fullInputString = dcap + inputFolder + "/" + inputString + ".root";
  //std::cout << fullInputString << std::endl;
  TFile *inputFile = TFile::Open(fullInputString.c_str());
  TTree* tree_= (TTree*) inputFile->Get("mt2");
      
  TFile *out = new TFile(outputFile.c_str(), "RECREATE");
  TTree *clone = new TTree("mt2", "post processed baby tree for mt2 analysis");

  clone = tree_->CloneTree(-1, "fast"); 
  clone->SetName("mt2");

  Float_t weight_btagsf;
  Float_t weight_btagsf_err;
  Float_t weight_btagsf_UP;
  Float_t weight_btagsf_DN;

  int isData; 
  clone->SetBranchAddress("isData",&isData);
  //clone->GetEntry(0);
  
  Int_t nobj;
  clone->SetBranchAddress(("n"+objectName).c_str(), &nobj);
  
  Float_t obj_pt[100];
  clone->SetBranchAddress((objectName+"_pt").c_str(), obj_pt);
  Float_t obj_eta[100];
  clone->SetBranchAddress((objectName+"_eta").c_str(), obj_eta);
  Int_t obj_mcFlavour[100];
  clone->SetBranchAddress((objectName+"_mcFlavour").c_str(), obj_mcFlavour);
  Float_t obj_btagCSV[100];
  clone->SetBranchAddress((objectName+"_btagCSV").c_str(), obj_btagCSV);
  
  Int_t evt;
  clone->SetBranchAddress("evt", &evt);
  
  TBranch* b1 = clone->Branch("weight_btagsf", &weight_btagsf, "weight_btagsf/F");
  TBranch* b2 = clone->Branch("weight_btagsf_UP", &weight_btagsf_UP, "weight_btagsf_UP/F");
  TBranch* b3 = clone->Branch("weight_btagsf_DN", &weight_btagsf_DN, "weight_btagsf_DN/F");
  
  int nEntries = clone->GetEntries();
  std::cout << "Starting loop over " << nEntries << " entries..." << std::endl;

  for(int i = 0; i < nEntries; ++i) {
    
    if(i%50000 == 0)
      std::cout << "Entry "<< i << "/" << nEntries << std::endl;
    
    clone->GetEntry(i);
    
    weight_btagsf = 1.;
    weight_btagsf_err = 0.;
    weight_btagsf_UP = 1.;
    weight_btagsf_DN = 1.;
    
    
    if( objectName != "jet" || isData == 1 ) ;
    else{
      
      get_weight_btag(nobj, obj_pt, obj_eta, obj_mcFlavour, obj_btagCSV, weight_btagsf, weight_btagsf_err);
      weight_btagsf_UP = weight_btagsf + weight_btagsf_err;
      weight_btagsf_DN = weight_btagsf - weight_btagsf_err;  
	
    }

    b1->Fill();
    b2->Fill();
    b3->Fill();
    
  }
  //-------------------------------------------------------------


  clone->Write();
  delete clone;
  out->Close();
  return 0;
  
}
