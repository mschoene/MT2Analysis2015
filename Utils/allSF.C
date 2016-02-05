/*
How to run:
root -l
.L allSF.C+
allSF(inputString, inputFolder, outputFile, treeName, objectName)
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

#include "BTagCalibrationStandalone.cc"   // https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration


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






int allSF(string inputString,
	   string inputFolder,
	   string outputFile,
	   string treeName,
	   string objectName)
{
  //Getting the lepton scale factor histograms/////////////////
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
  std::cout << "Using Veto Electrons ID, MiniIso 0.1 lepton scale factors" << std::endl;
  std::cout << "Be aware that Veto Electrons are not suited for selecting Electrons." << std::endl;
  std::cout << std::endl;


  //Top pt reweighting 
  //Values from Run1 still valid for now, see here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  float a = 0.156;
  float b = -0.00137;



  //Add all files in the input folder 
  string dcap = inputFolder.find("pnfs")!=std::string::npos ? "dcap://t3se01.psi.ch:22125/" : "";
  string fullInputString = dcap + inputFolder + "/" + inputString + ".root";
  //std::cout << fullInputString << std::endl;
  TFile *inputFile = TFile::Open(fullInputString.c_str());
  TTree* tree_= (TTree*) inputFile->Get("mt2");
      
  TFile *out = new TFile(outputFile.c_str(), "RECREATE");
  TTree *clone = new TTree("mt2", "post processed baby tree for mt2 analysis");

  clone = tree_->CloneTree( -1, "fast"); 
  clone->SetName("mt2");


  Float_t weight_btagsf;
  Float_t weight_btagsf_heavy_UP;
  Float_t weight_btagsf_heavy_DN;
  Float_t weight_btagsf_light_UP;
  Float_t weight_btagsf_light_DN;
  Float_t weight_lepsf;
  Float_t weight_lepsf_UP;
  Float_t weight_lepsf_DN;
  Float_t weight_toppt;
  Float_t weight_isr;


  int isData; 
  clone->SetBranchAddress("isData",&isData);
 
  int evt_id; 
  clone->SetBranchAddress("evt_id",&evt_id);
  /*
  clone->GetEntry(0); 
  if( evt_id < 1000 ) 
    std::cout << "ISR only calculated for signal samples, skipping filling them here " << std::endl;
  if( isData == 1 ) 
    std::cout << "Scale Factors are only applied to MC, skipping data " << std::endl;
  */
 
  ////// Lepton Efficiency SF
  Int_t nlep;
  clone->SetBranchAddress("nlep", &nlep);
  Float_t lep_pt[100];
  clone->SetBranchAddress("lep_pt", lep_pt);
  Float_t lep_eta[100];
  clone->SetBranchAddress("lep_eta", lep_eta);
  Int_t lep_pdgId[100];
  clone->SetBranchAddress("lep_pdgId", lep_pdgId);
  
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
  Int_t GenPart_pdgId[100];
  clone->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);


  TBranch* b1 = clone->Branch("weight_btagsf"         , &weight_btagsf         , "weight_btagsf/F"         );
  TBranch* b2 = clone->Branch("weight_btagsf_heavy_UP", &weight_btagsf_heavy_UP, "weight_btagsf_heavy_UP/F");
  TBranch* b3 = clone->Branch("weight_btagsf_heavy_DN", &weight_btagsf_heavy_DN, "weight_btagsf_heavy_DN/F");
  TBranch* b4 = clone->Branch("weight_btagsf_light_UP", &weight_btagsf_light_UP, "weight_btagsf_light_UP/F");
  TBranch* b5 = clone->Branch("weight_btagsf_light_DN", &weight_btagsf_light_DN, "weight_btagsf_light_DN/F");
  TBranch* b6 = clone->Branch("weight_lepsf", &weight_lepsf, "weight_lepsf/F");
  TBranch* b7 = clone->Branch("weight_lepsf_UP", &weight_lepsf_UP, "weight_lepsf_UP/F");
  TBranch* b8 = clone->Branch("weight_lepsf_DN", &weight_lepsf_DN, "weight_lepsf_DN/F");
  TBranch* b9 = clone->Branch("weight_toppt", &weight_toppt, "weight_toppt/F");
  TBranch* b10 = clone->Branch("weight_isr", &weight_isr, "weight_isr/F");



  
  int nEntries = clone->GetEntries();
  std::cout << "Starting first loop over " << nEntries << " entries for average weight for top SF..." << std::endl;

  //top pt scale factors should not change the overall normalization
  //first loop to get the average of the scale factors
  double average = 0;
  for(int i = 0; i < nEntries; ++i) { 
    if(i%50000 == 0)
      std::cout << "Entry "<< i << "/" << nEntries << std::endl;
    clone->GetEntry(i);

    if( evt_id >399 || evt_id < 300 ) continue;

    if( isData == 1 ) continue;

    weight_toppt = 1.;
 
    if(nGenPart>0)
      for(int o=0; o<nGenPart; ++o){
	if(GenPart_status[o] != 62) continue;
	//Only apply weight to top or antitop
	if( abs(GenPart_pdgId[o])==6 ){
	  weight_toppt *= exp( a + b * GenPart_pt[o] );
	}
	else continue;


      }//end loop over objects
  
    average += weight_toppt;   
  }
  //-------------------------------------------------------------end of run over events to get the average for the top pt reweighting

  average /=  double(nEntries);
  std::cout << average << std::endl;

  ////////////////Second Loop, filling all the new branches////////////
  std::cout << "Starting second loop over " << nEntries << " entries..." << std::endl;
  for(int i = 0; i < nEntries; ++i) {
    
    if(i%50000 == 0)
      std::cout << "Entry "<< i << "/" << nEntries << std::endl;
    
    clone->GetEntry(i);
    
    weight_btagsf = 1.;
    weight_btagsf_heavy_UP = 1.;
    weight_btagsf_heavy_DN = 1.;
    weight_btagsf_light_UP = 1.;
    weight_btagsf_light_DN = 1.;
    weight_lepsf = 1.;
    weight_lepsf_UP = 1.;
    weight_lepsf_DN = 1.;
    weight_toppt = 1.;
    weight_isr=1.;

    if( isData == 1 ) continue;

    ///////////Add ISR scale factor////////////
    if( evt_id < 1000 ) ;
    else{
      TLorentzVector s;
      
      if(nGenPart>0)
	for(int o=0; o<nGenPart; ++o){
	  if(GenPart_status[o] != 62) continue;
	  TLorentzVector s_;
	  s_.SetPtEtaPhiM(GenPart_pt[o], GenPart_eta[o], GenPart_phi[o], GenPart_mass[o]);
	  s+=s_;	  
	}
      
      float pt_hard = s.Pt();
      if( pt_hard < 0. )         continue;
      else if( pt_hard < 400. )  weight_isr = 1.0;
      else if( pt_hard < 600. )  weight_isr = 1.15;
      else if( pt_hard >= 600. ) weight_isr = 1.30;

    }

    /////////Add Top pt scale factor//////////
    if( evt_id >399 || evt_id < 300 ) ;
    else if(nGenPart>0){
      for(int o=0; o<nGenPart; ++o){
      	if(GenPart_status[o] != 62) continue;
	//Only apply weight to top or antitop
	if( abs(GenPart_pdgId[o])==6 )
	  weight_toppt *= exp( a + b * GenPart_pt[o] );
	else continue;
      }//end loop over objects
      weight_toppt /= average;
    }


    /////////Add b-tagging scale factor//////////     
    get_weight_btag(njet, jet_pt, jet_eta, jet_mcFlavour, jet_btagCSV, weight_btagsf, weight_btagsf_heavy_UP, weight_btagsf_heavy_DN, weight_btagsf_light_UP, weight_btagsf_light_DN);
   

    /////////Add lepton scale factor//////////
    if(nlep>0){
      Float_t uncert = 0; //Place holder for total uncertainty
      Float_t central = 1; 
      Float_t err = 0;

      for(int o=0; o<nlep; ++o){
	//Electrons
	if (abs( lep_pdgId[o]) == 11) {
	  Int_t binx = h_elSF->GetXaxis()->FindBin(lep_pt[o]);
	  Int_t biny = h_elSF->GetYaxis()->FindBin(fabs(lep_eta[o]));
	  central = h_elSF->GetBinContent(binx,biny);
	  err  = h_elSF->GetBinError(binx,biny);
	  if (central > 1.2 || central < 0.8) 
	    std::cout<<"STRANGE: Electron with pT/eta of "<<lep_pt[o]<<"/"<<lep_eta[o]<<". SF is "<< central <<std::endl;
	} //else Muons
	else if (abs( lep_pdgId[o]) == 13) {
	  Int_t binx = h_muSF->GetXaxis()->FindBin(lep_pt[o]);
	  Int_t biny = h_muSF->GetYaxis()->FindBin(fabs(lep_eta[o]));
	  if ( binx >7 ) binx = 7; //overflow bin empty for the muons...
	  central = h_muSF->GetBinContent(binx,biny);
	  err  = 0.014; // adding in quadrature 1% unc. on ID and 1% unc. on ISO
	  if (central > 1.2 || central < 0.8) 
	    std::cout<<"STRANGE: Muon with pT/eta of "<<lep_pt[o]<<"/"<< fabs(lep_eta[o]) <<". SF is "<< central <<std::endl;
	} 
	weight_lepsf *= central;
	//uncertainties are supposed to be summed up linearly for the lepton SF
	uncert += err; 
	
      }//end of loop over objects
    
      weight_lepsf_UP = central + uncert; 
      weight_lepsf_DN = central - uncert;
	 
    }//end of lepton sf


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
    
  }
  //-------------------------------------------------------------


  clone->Write();
  delete clone;
  out->Close();
  return 0;
  
}
