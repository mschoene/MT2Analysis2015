#ifndef LEPTONSF_H
#define LEPTONSF_H


#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

#include <string>
#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH2.h"




struct lepSF{
  float sf = 1.0;
  float up = 1.0;
  float dn = 1.0;
};



TH2D* h_elSF        = 0;
TH2D* h_elSF_trk    = 0;
TH2D* h_muSF        = 0;

TH2D* h_eff_full_mu = 0;
TH2D* h_eff_full_el = 0;

TH2D* h_fast_elSF   = 0;
TH2D* h_fast_muSF   = 0;

bool setElHistos( std::string filenameID, std::string fileNameTrk, bool useLoose=0);
bool setMuHistos( std::string filenameID, std::string filenameISO, std::string filenamedxyz );
bool setVetoEffHistos( std::string filename);

lepSF getLepSF( float pt, float eta, int pdgId);

lepSF getVetoEffLepSF( float pt, float eta, int pdgId);




#endif
