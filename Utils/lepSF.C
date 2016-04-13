/*
How to run:
root -l
.L lepSF.C+
lepSF(inputString, inputFolder, outputFile, treeName, objectName)
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


int lepSF(string inputString,
	  string inputFolder,
	  string outputFile,
	  string treeName,
	  string objectName)
{

  //Getting the scale factor histograms/////////////////
  //Electrons//
  std::string filename = "kinematicBinSFele.root";
  TFile * f = new TFile(filename.c_str() );
  if (!f->IsOpen()) std::cout << " ERROR: Could not find scale factor file " << filename << std::endl;
  //TH2D* h_id = (TH2D*) f->Get("CutBasedLoose");
  TH2D* h_id = (TH2D*) f->Get("CutBasedVeto");
  TH2D* h_iso = (TH2D*) f->Get("MiniIso0p1_vs_AbsEta");
  if (!h_id || !h_iso) std::cout << "ERROR: Could not find scale factor histogram"<< std::endl;
  TH2D* h_elSF = (TH2D*) h_id->Clone("h_elSF");
  h_elSF->SetDirectory(0);
  h_elSF->Multiply(h_iso);

  //Muons//
  std::string filenameID = "TnP_MuonID_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root";
  std::string filenameISO = "TnP_MuonID_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root";
  TFile * f1 = new TFile(filenameID.c_str() );
  TFile * f2 = new TFile(filenameISO.c_str() );
  if (!f1->IsOpen()) { std::cout<<" ERROR: Could not find ID scale factor file "<<filenameID<<std::endl; return 0;}
  if (!f2->IsOpen()) { std::cout<<"ERROR: Could not find ISO scale factor file "<<filenameISO<<std::endl; return 0;}
  TH2D* h_id_mu = (TH2D*) f1->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_tag_IsoMu20_pass");
  TH2D* h_iso_mu = (TH2D*) f2->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_PF_pass_&_tag_IsoMu20_pass");
  if (!h_id_mu || !h_iso_mu) { std::cout<<"ERROR: Could not find scale factor histogram"<<std::endl; return 0;}
  TH2D* h_muSF = (TH2D*) h_id_mu->Clone("h_muSF");
  h_muSF->SetDirectory(0);
  h_muSF->Multiply(h_iso_mu);




  std::cout << std::endl;
  std::cout << "Using Loose Muon ID, MiniIso 0.2 lepton scale factors" << std::endl;
  std::cout << "Using Veto Muon ID, MiniIso 0.1 lepton scale factors" << std::endl;
  std::cout << "Be aware that this is not suitable for selecting Electrons." << std::endl;
  std::cout << std::endl;


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

  Float_t weight_lepsf;
  Float_t weight_lepsf_UP;
  Float_t weight_lepsf_DN;
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
  
  int isData; 
  clone->SetBranchAddress("isData",&isData);

  Int_t nobj;
  clone->SetBranchAddress(("n"+objectName).c_str(), &nobj);
  
  Float_t obj_pt[100];
  clone->SetBranchAddress((objectName+"_pt").c_str(), obj_pt);
  Float_t obj_eta[100];
  clone->SetBranchAddress((objectName+"_eta").c_str(), obj_eta);
  Int_t obj_pdgId[100];
  clone->SetBranchAddress((objectName+"_pdgId").c_str(), obj_pdgId);

  TBranch* b1 = clone->Branch("weight_lepsf", &weight_lepsf, "weight_lepsf/F");
  TBranch* b2 = clone->Branch("weight_lepsf_UP", &weight_lepsf_UP, "weight_lepsf_UP/F");
  TBranch* b3 = clone->Branch("weight_lepsf_DN", &weight_lepsf_DN, "weight_lepsf_DN/F");


  int nEntries = clone->GetEntries();
  for(int i = 0; i < nEntries; ++i) {
    
    clone->GetEntry(i);
    
    weight_lepsf=1.;
    weight_lepsf_UP=1.;
    weight_lepsf_DN=1.;

    
    if( objectName == "lep" )
      if(nobj>0){

	Float_t uncert = 0; //Place holder for total uncertainty
	Float_t central = 1; 
	Float_t err = 0;

	for(int o=0; o<nobj; ++o){
	  //Electrons
	  if (abs( obj_pdgId[o]) == 11) {
	    Int_t binx = h_elSF->GetXaxis()->FindBin(obj_pt[o]);
	    Int_t biny = h_elSF->GetYaxis()->FindBin(fabs(obj_eta[o]));
	    central = h_elSF->GetBinContent(binx,biny);
	    err  = h_elSF->GetBinError(binx,biny);
	    if (central > 1.2 || central < 0.8) 
	      std::cout<<"STRANGE: Electron with pT/eta of "<<obj_pt[o]<<"/"<<obj_eta[o]<<". SF is "<< central <<std::endl;
	  } //else Muons
	  else if (abs( obj_pdgId[o]) == 13) {
	    Int_t binx = h_muSF->GetXaxis()->FindBin(obj_pt[o]);
	    Int_t biny = h_muSF->GetYaxis()->FindBin(fabs(obj_eta[o]));
	    if ( binx >7 ) binx = 7; //overflow bin empty for the muons...
	    central = h_muSF->GetBinContent(binx,biny);
	    err  = 0.014; // adding in quadrature 1% unc. on ID and 1% unc. on ISO
	    if (central > 1.2 || central < 0.8) 
	      std::cout<<"STRANGE: Muon with pT/eta of "<<obj_pt[o]<<"/"<< fabs(obj_eta[o]) <<". SF is "<< central <<std::endl;
	  } 
	  weight_lepsf *= central;
	  //uncertainties are supposed to be summed up linearly for the lepton SF
	  uncert += err; 
	
	}//end of loop over objects
    
	weight_lepsf_UP = central + uncert; 
	weight_lepsf_DN = central - uncert;
	  
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
