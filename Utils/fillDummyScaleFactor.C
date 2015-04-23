#include "TFile.h"
#include "TH2F.h"


void fillDummyScaleFactor(){
  
  TFile* f= new TFile("scaleFactors.root", "RECREATE");
  
  int nBinsx=20;
  float minPt=0., maxPt=100.;

  int nBinsy=25;
  float minEta=0., maxEta=2.5.;

  float Pt_t=20.;
  float Eta_t=2.5;
  
  TH2F* h= new TH2F("h", "Lepton scale factor #eta vs p_{T}", nBinsx, minPt, maxPt, nBinsy, minEta, maxEta);
  
  for (int x=1; x<=nBinsx; ++x){
    
    float X=h->GetXaxis()->GetBinCenter(x);
    float Y=0;
    if (X < Pt_t){
      
      for (int y=1; y<=nBinsy; ++y){
	Y=h->GetYaxis()->GetBinCenter(y);
	if (Y < Eta_t){
	  h->Fill(X, Y,  0.95);
	  h->SetBinError(x, y,  0.05);
	}
	else {
	  h->Fill(X, Y, 1.0);
	  h->SetBinError(x, y,  0.0);
	}
	
      }

    }
    
    else 
      for (int y=1; y<=nBinsy; ++y) {
	Y=h->GetYaxis()->GetBinCenter(y);
	if (Y < Eta_t){
          h->Fill(X, Y,  0.99);
          h->SetBinError(x, y,  0.02);
        }
	else{
	  h->Fill(X, Y, 1.0);
	  h->SetBinError(x, y,  0.0);
	}
      }

  }
  
  for (int y=1; y<=nBinsy; ++y){
    X=maxPt+1.;
    Y=h->GetYaxis()->GetBinCenter(y);
    if (Y < Eta_t){
      h->Fill(X, Y,  0.99);
      h->SetBinError(x, y,  0.02);
    }
    else{
      h->Fill(X, Y, 1.0);
      h->SetBinError(x, y,  0.0);
    }
  }
    
  f->cd();
  f->Write();
  f->Close();

}
