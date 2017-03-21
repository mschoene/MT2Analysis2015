#include "leptonSF.h"

using namespace std;



bool setElHistos( std::string filenameID, std::string filenameTrk, bool useLoose){

  if( filenameID == "")
    filenameID = "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/scaleFactors.root";
  if(filenameTrk == "")
    filenameTrk = "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/egammaEffi.txt_EGM2D.root";

  TFile* f_ele = new TFile(filenameID.c_str() );
  TFile* f_eletrk = new TFile(filenameTrk.c_str() );
 
  if (!f_ele->IsOpen() || !f_eletrk->IsOpen()) 
    std::cout << " ERROR: Could not find scale factor file for ELECTRONS " << filenameID << std::endl; 
 
  if( useLoose)
    std::cout << "You want to use the loose WP, you have to implement it here first" << std::endl;

  TH2D* h_id  = (TH2D*) f_ele->Get("GsfElectronToCutBasedSpring15V");
  TH2D* h_iso = (TH2D*) f_ele->Get("MVAVLooseElectronToMini");
  TH2D* h_trk = (TH2D*) f_eletrk->Get("EGamma_SF2D");
  if (!h_id || !h_iso || !h_trk) std::cout << "ERROR: Could not find scale factor histogram"<< std::endl;

  h_elSF = (TH2D*) h_id->Clone("h_elSF");
  h_elSF->SetDirectory(0);
  h_elSF->Multiply(h_iso);
  h_elSF_trk = (TH2D*) h_trk->Clone("h_elSF_trk");
  h_elSF_trk->SetDirectory(0);

  return true;
}


bool setMuHistos( std::string filenameID, std::string filenameISO, std::string filenamedxyz ){
  
  if( filenameID == "")
    filenameID =   "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root";
  if( filenameISO == "")
    filenameISO =  "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/TnP_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root";
  if( filenamedxyz == "")
    filenamedxyz = "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/TnP_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root";
  
  TFile * f1 = new TFile(filenameID.c_str() );
  TFile * f2 = new TFile(filenameISO.c_str() );
  TFile * f3 = new TFile(filenamedxyz.c_str() );

  if (!f1->IsOpen()) { std::cout<<" ERROR: Could not find ID scale factor file "<<filenameID<<std::endl; return 0;}
  if (!f2->IsOpen()) { std::cout<<"ERROR: Could not find ISO scale factor file "<<filenameISO<<std::endl; return 0;}
  if (!f3->IsOpen()) { std::cout<<"ERROR: Could not find dxy dz scale factor file "<<filenamedxyz<<std::endl; return 0;}

  TH2D* h_id_mu   = (TH2D*) f1->Get("SF");
  TH2D* h_iso_mu  = (TH2D*) f2->Get("SF");
  TH2D* h_dxyz_mu = (TH2D*) f3->Get("SF");

  if (!h_id_mu || !h_iso_mu  || !h_dxyz_mu) { std::cout<<"ERROR: Could not find scale factor histogram"<<std::endl; return 0;}
  h_muSF = (TH2D*) h_id_mu->Clone("h_muSF");
  h_muSF->SetDirectory(0);
  h_muSF->Multiply(h_iso_mu);
  h_muSF->Multiply(h_dxyz_mu);
  

  return true;  
}

// add version for with muon track corr 

// add version for fast sim






bool setVetoEffHistos( std::string filename){

  if( filename == "")
    filename = "/mnt/t3nfs01/data01/shome/mmasciov/lepSF/vetoeff_emu_etapt_lostlep.root";

  TFile* f_eff_full = new TFile(filename.c_str() );
   
  if (!f_eff_full->IsOpen()) std::cout << " ERROR: Could not find scale factor file for veto eff" << filename << std::endl; 
 
  h_eff_full_mu = (TH2D*) f_eff_full->Get("h_mu_comb_eff");
  h_eff_full_el = (TH2D*) f_eff_full->Get("h_ele_comb_eff");
  if(!h_eff_full_mu || !h_eff_full_el ) {std::cout << " ERROR: Could not find the 2D histogram in your files " << std::endl; return 0;}
  h_eff_full_mu->SetDirectory(0);
  h_eff_full_el->SetDirectory(0);

  f_eff_full->Close();

  return true;
}




lepSF getLepSF( float pt, float eta, int pdgId){

  lepSF weights;

  if( !h_muSF || !h_elSF || !h_elSF_trk ){
    std::cout << "There are not scale factor histograms for me to read from " << std::endl;
    return weights;
  }

  Float_t uncert = 1; //Place holder for total uncertainty
  Float_t central = 1; 
  Float_t err = 0;
  Float_t uncert_UP = 0; 	Float_t uncert_DN = 0; 
      
  float pt_cutoff = std::max( 10.1, std::min( 100., double(pt) ) );
 	
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
  else if( pdgId == 13) {
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

  weights.sf = central;
  weights.up = uncert_UP;
  weights.dn = uncert_DN;

  return weights;
}









lepSF getVetoEffLepSF( float pt, float eta, int pdgId){

  lepSF weights;

  if( !h_muSF || !h_elSF || !h_elSF_trk ){
    std::cout << "There are not scale factor histograms for me to read from" << std::endl;
      return weights;
  }

	
  float central=1.;
  float err=0;
  float uncert=0;
  float fast_central=1.;
  float fast_err=0.;
  float uncert_fast=0;
	
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
  // weight_lepsf_0l_UP *= ( 1. + unc_UP_0l);
  // weight_lepsf_0l_DN *= ( 1. - unc_UP_0l);


  weights.sf = central;
  weights.up = ( 1. + unc_UP_0l);
  weights.dn = ( 1. - unc_UP_0l);

 
  return weights;
}












lepSF getLepSF_fast( float pt, float eta, int pdgId){

  lepSF weights;

  if( !h_muSF || !h_elSF || !h_elSF_trk ){
    std::cout << "There are not scale factor histograms for me to read from " << std::endl;
    return weights;
  }

  Float_t uncert = 1; //Place holder for total uncertainty
  Float_t central = 1; 
  Float_t err = 0;
  Float_t uncert_UP = 0; 	Float_t uncert_DN = 0; 
     
  Float_t fast_central = 1; 
  Float_t fast_err = 0;


  float pt_cutoff = std::max( 10.1, std::min( 100., double(pt) ) );
 	
  //Electrons
  if( pdgId == 11) {
 
   // Int_t binx = h_elSF->GetXaxis()->FindBin(pt_cutoff);
   //  Int_t biny = h_elSF->GetYaxis()->FindBin(eta);
   //  central = h_elSF->GetBinContent(binx,biny);
   //  err  = h_elSF->GetBinError(binx,biny);
	  
   //  int binx_trk = h_elSF_trk->GetXaxis()->FindBin(eta);
   //  int biny_trk = 1; // hardcoding for now - only one bin in pt (hist starts at 20)
   //  central *= h_elSF_trk->GetBinContent(binx_trk,biny_trk);
   //  float trk_err = h_elSF_trk->GetBinError(binx_trk,biny_trk);
   //  if (pt_cutoff < 20. || pt_cutoff > 80.) err = sqrt(err*err + trk_err*trk_err + 0.01*0.01); 
   //  else err = sqrt(err*err + trk_err*trk_err);
	  
   //  if (central > 1.3 || central < 0.3) 
   //    std::cout<<"STRANGE: Electron with pT/eta of "<< pt <<"/"<< eta <<". SF is "<< central <<std::endl;
   //  uncert_UP = central + err;
   //  uncert_DN = central - err;
	  
    Int_t fast_binx = h_fast_elSF->GetXaxis()->FindBin(pt_cutoff);
    Int_t fast_biny = h_fast_elSF->GetYaxis()->FindBin(eta);
    fast_central = h_fast_elSF->GetBinContent(fast_binx,fast_biny);
    fast_err = h_fast_elSF->GetBinError(fast_binx,fast_biny);
    fast_err = sqrt(fast_err*fast_err + 0.05*0.05); // 5% systematic uncertainty

    if( fast_central > 1.3 || fast_central < 0.7 )
      std::cout << "Strange FastSim Electron with pT/eta of" << pt <<"/"<< eta <<". SF is "<< fast_central <<std::endl;

    central   = fast_central;
    uncert_UP = ( fast_central + fast_err );
    uncert_DN = ( fast_central - fast_err );
    // uncert_UP += fast_err ;
    //uncert_DN -= fast_err ;



  } //else Muons
  else if( pdgId == 13) {

    // Int_t binx = h_muSF->GetXaxis()->FindBin(pt_cutoff);
    // Int_t biny = h_muSF->GetYaxis()->FindBin(fabs(eta));
	  
    // //	  float central_trk = 1;
    // //	  Int_t binx_trk = h_muTrk_hi->GetXaxis()->FindBin(  lep_eta[o] );
    // //	  if( binx_trk>10 ) binx_trk = 10;
    // //	  else if( binx_trk<1 ) binx_trk = 1;
    // //	  central_trk = h_muTrk_hi->GetBinContent( binx_trk );
    // //	  
    // //	  
    // //	  if ( binx >7 ) binx = 7; //overflow bin empty for the muons...
    // //	  central = h_muSF->GetBinContent(binx,biny);
    // //	  
    // //	  central *= central_trk;
	  

    // central = h_muSF->GetBinContent(binx,biny);
    // err  = 0.03; //current recomendation is 3% //   err  = 0.014; // adding in quadrature 1% unc. on ID and 1% unc. on ISO
    // if (central > 1.2 || central < 0.8) 
    //   std::cout<<"STRANGE: Muon with pT/eta of "<<pt<<"/"<< fabs(eta) <<". SF is "<< central <<std::endl;
    // uncert_UP = central + err;
    // uncert_DN = central - err;
	  
    Int_t fast_binx = h_fast_muSF->GetXaxis()->FindBin(pt);
    Int_t fast_biny = h_fast_muSF->GetYaxis()->FindBin(fabs(eta));
    fast_central = h_fast_muSF->GetBinContent(fast_binx,fast_biny);
    fast_err  = h_fast_muSF->GetBinError(fast_binx,fast_biny);
    if( pt < 20)
      fast_err= sqrt(fast_err*fast_err+ 0.03*0.03); // 3% systematic uncertainty
    else
      fast_err= sqrt(fast_err*fast_err+ 0.01*0.01); // 1% systematic uncertainty

    if( fast_central > 1.3 || fast_central < 0.7 )
      std::cout << "Strange FastSim Muon with pT/eta of" <<pt<<"/"<<eta<<". SF is "<< fast_central <<std::endl;

    central   = fast_central;

    uncert_UP = ( fast_central + fast_err );
    uncert_DN = ( fast_central - fast_err );

    //uncert_UP += fast_err ;
    //uncert_DN -= fast_err ;

  }//done with one  electron/muon 

  weights.sf = central;
  weights.up = uncert_UP;
  weights.dn = uncert_DN;

  return weights;
}
