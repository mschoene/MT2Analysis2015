#include "TFile.h"
#include "TH2F.h"


void main(){
  
  TFile* f= new TFile("scaleFactors.root", "RECREATE");
  
  TH2F* h= new TH2F("h", "Lepton scale factor #eta vs p_{T}", 20, 0., 100., 50, 0., 5.);
  
  for (int x=1; x<=20; ++x){
    
    float X=h->GetXaxis()->GetBinCenter(x);
    float Y=0;
    if (X < 20.){
      
      for (int y=1; y<=50; ++y){
	Y=h->GetYaxis()->GetBinCenter(y);
	if (Y < 2.5)
	  h->Fill(X, Y,  0.99);
	  
	else h->Fill(X, Y, 0.95);
	
      }

    }
    
    else 
      for (int y=1; y<=50; ++y) {
	Y=h->GetYaxis()->GetBinCenter(y);
	h->Fill(X, Y, 0.95);
      }

  }
  
  for (int y=1; y<=50; ++y){
    X=100.+1.;
    Y=h->GetYaxis()->GetBinCenter(y);
    h->Fill(X, Y, 0.95);
  }
  
  for (int x=1; x<=20; ++x)
    for (int y=1; y<=50; ++y)
      h->SetBinError(x, y, 0.05);
  
  f->cd();
  f->Write();
  f->Close();

}
