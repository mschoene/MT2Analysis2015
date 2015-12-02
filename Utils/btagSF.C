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

#include "BTagCalibrationStandalone.h"   // https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration

using namespace std;


// setup calibration readers 
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X  --  official SFs
// https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration  --  calibration reader documentations
BTagCalibration *calib = new BTagCalibration("csvv2", "/shome/casal/btagsf/CSVv2.csv"); // 25 ns official version of SFs
BTagCalibrationReader *reader_heavy    = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "central");  // central
BTagCalibrationReader *reader_heavy_UP = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "up");       // sys up
BTagCalibrationReader *reader_heavy_DN = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "down");     // sys down
BTagCalibrationReader *reader_light    = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "comb"  , "central");  // central
BTagCalibrationReader *reader_light_UP = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "comb"  , "up");       // sys up
BTagCalibrationReader *reader_light_DN = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "comb"  , "down");     // sys down

TFile *f_btag_eff = new TFile("/shome/casal/btagsf/btageff__ttbar_powheg_pythia8_25ns.root"); // Dominick's b-tagging efficiencies
TH2D* h_btag_eff_b    = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_b"   );
TH2D* h_btag_eff_c    = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_c"   );
TH2D* h_btag_eff_udsg = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_udsg");


void get_SF_btag(float pt, float eta, int mcFlavour, float &SF, float &SFup, float &SFdown){

  BTagEntry::JetFlavor flavour = BTagEntry::FLAV_UDSG;
  if      ( abs(mcFlavour)==5 ) flavour = BTagEntry::FLAV_B;
  else if ( abs(mcFlavour)==4 ) flavour = BTagEntry::FLAV_C;
  
  float pt_cutoff  = std::max(30. ,std::min(669., double(pt)));
  float eta_cutoff = std::min(2.39,fabs(double(eta)));

  if ( flavour==BTagEntry::FLAV_UDSG ){
    SF     = reader_light   ->eval(flavour,eta_cutoff, pt_cutoff);
    SFup   = reader_light_UP->eval(flavour,eta_cutoff, pt_cutoff);
    SFdown = reader_light_DN->eval(flavour,eta_cutoff, pt_cutoff);
  }
  else {
    SF     = reader_heavy   ->eval(flavour,eta_cutoff, pt_cutoff);
    SFup   = reader_heavy_UP->eval(flavour,eta_cutoff, pt_cutoff);
    SFdown = reader_heavy_DN->eval(flavour,eta_cutoff, pt_cutoff);
  }

}

float getBtagEffFromFile(float pt, float eta, int mcFlavour){
  if(!h_btag_eff_b || !h_btag_eff_c || !h_btag_eff_udsg) {
    std::cout << "ERROR: missing input hists" << std::endl;
    return 1.;
  }
  
  // only use pt bins up to 400 GeV for charm and udsg
  float pt_cutoff = std::max(20.,std::min(399.,double(pt)));
  TH2D* h(0);
  if (abs(mcFlavour) == 5) {
    h = h_btag_eff_b;
    // use pt bins up to 600 GeV for b
    pt_cutoff = std::max(20.,std::min(599.,double(pt)));
  }
  else if (abs(mcFlavour) == 4) {
    h = h_btag_eff_c;
  }
  else {
    h = h_btag_eff_udsg;
  }
  
  int binx = h->GetXaxis()->FindBin(pt_cutoff);
  int biny = h->GetYaxis()->FindBin(fabs(eta));
  return h->GetBinContent(binx,biny);
}

void get_weight_btag(int nobj, float* obj_pt, float* obj_eta, int* obj_mcFlavour, float* obj_btagCSV, float &wtbtag, float &wtbtagUp_heavy, float &wtbtagDown_heavy, float &wtbtagUp_light, float &wtbtagDown_light){

  float mcTag = 1.;
  float mcNoTag = 1.;
  float dataTag = 1.;
  float dataNoTag = 1.;

  float errHup   = 0;
  float errHdown = 0;
  float errLup   = 0;
  float errLdown = 0;

  for(int indj=0; indj < nobj; ++indj){ //Here we loop over all selected jets ( for our case, pt>30, PF loose ID, etc )
    
    float csv = obj_btagCSV[indj]; ////here we get the CSV btag value
    int   mcFlavour = abs(obj_mcFlavour[indj]);
    float eta = fabs(obj_eta[indj]);
    float pt  = obj_pt[indj];

    if(eta > 2.5) continue;
    if(pt  < 20 ) continue;
    //if(mcFlavour==0) continue; //for jets with flavour 0, we ignore.

    
    // get efficiency
    float eff = getBtagEffFromFile(pt, eta, mcFlavour);


    bool istag = csv > 0.890 && eta < 2.5 && pt>20;
    float SF = 0.;
    float SFup = 0.;
    float SFdown = 0.;

    // get SF
    get_SF_btag(pt, eta, mcFlavour, SF, SFup, SFdown);
 
    if(istag){
      mcTag   *= eff;
      dataTag *= eff*SF;

      if(mcFlavour==5 || mcFlavour ==4) {
	errHup   += (SFup - SF  )/SF;
	errHdown += (SF - SFdown)/SF;
      }
      else {
	errLup   += (SFup - SF  )/SF;
	errLdown += (SF - SFdown)/SF;
      }

    }
    else{
      mcNoTag   *= (1 - eff);
      dataNoTag *= (1 - eff*SF);

      if( mcFlavour==5 || mcFlavour==4 ) {
	errHup   += -eff*(SFup - SF  )/(1-eff*SF);
	errHdown += -eff*(SF - SFdown)/(1-eff*SF);	
      }
      else {
	errLup   += -eff*(SFup - SF  )/(1-eff*SF);
	errLdown += -eff*(SF - SFdown)/(1-eff*SF);	
      }
    }

  }

  wtbtag = (dataNoTag * dataTag ) / ( mcNoTag * mcTag );

  wtbtagUp_heavy   = wtbtag*( 1 + errHup   );
  wtbtagUp_light   = wtbtag*( 1 + errLup   );
  wtbtagDown_heavy = wtbtag*( 1 - errHdown );
  wtbtagDown_light = wtbtag*( 1 - errLdown );

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
  Float_t weight_btagsf_heavy_UP;
  Float_t weight_btagsf_heavy_DN;
  Float_t weight_btagsf_light_UP;
  Float_t weight_btagsf_light_DN;

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
  clone->SetBranchAddress((objectName+"_hadronFlavour").c_str(), obj_mcFlavour);
  //clone->SetBranchAddress((objectName+"_mcFlavour").c_str(), obj_mcFlavour);
  Float_t obj_btagCSV[100];
  clone->SetBranchAddress((objectName+"_btagCSV").c_str(), obj_btagCSV);
  
  
  TBranch* b1 = clone->Branch("weight_btagsf"         , &weight_btagsf         , "weight_btagsf/F"         );
  TBranch* b2 = clone->Branch("weight_btagsf_heavy_UP", &weight_btagsf_heavy_UP, "weight_btagsf_heavy_UP/F");
  TBranch* b3 = clone->Branch("weight_btagsf_heavy_DN", &weight_btagsf_heavy_DN, "weight_btagsf_heavy_DN/F");
  TBranch* b4 = clone->Branch("weight_btagsf_light_UP", &weight_btagsf_light_UP, "weight_btagsf_light_UP/F");
  TBranch* b5 = clone->Branch("weight_btagsf_light_DN", &weight_btagsf_light_DN, "weight_btagsf_light_DN/F");
  
  int nEntries = clone->GetEntries();
  std::cout << "Starting loop over " << nEntries << " entries..." << std::endl;


  for(int i = 0; i < nEntries; ++i) {
    
    if(i%50000 == 0)
      std::cout << "Entry "<< i << "/" << nEntries << std::endl;
    
    clone->GetEntry(i);
    
    weight_btagsf = 1.;
    weight_btagsf_heavy_UP = 1.;
    weight_btagsf_heavy_DN = 1.;
    weight_btagsf_light_UP = 1.;
    weight_btagsf_light_DN = 1.;
    
    
    if( objectName != "jet" || isData == 1 ) ;
    else{
      
      get_weight_btag(nobj, obj_pt, obj_eta, obj_mcFlavour, obj_btagCSV, weight_btagsf, weight_btagsf_heavy_UP, weight_btagsf_heavy_DN, weight_btagsf_light_UP, weight_btagsf_light_DN);
      
    }

    b1->Fill();
    b2->Fill();
    b3->Fill();
    b4->Fill();
    b5->Fill();
    
  }
  //-------------------------------------------------------------


  clone->Write();
  delete clone;
  out->Close();
  return 0;
  
}
