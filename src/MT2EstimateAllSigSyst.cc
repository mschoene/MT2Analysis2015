#include "../interface/MT2EstimateAllSigSyst.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TEfficiency.h"








MT2EstimateAllSigSyst::MT2EstimateAllSigSyst( const std::string& aname, const MT2Region& aregion, const std::string& asystName ) : MT2Estimate( aname, aregion ) {

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

//  int nBinsM=81;
//  double binWidthM=25.;
//  double binsM[nBinsM+1];
//  for (int b=0; b<=nBinsM; ++b)
//    binsM[b]=b*binWidthM;

  int nBinsM=93;
  double binWidthM=25.;
  double binsM[nBinsM+1];
  for (int b=0; b<=nBinsM; ++b)
    binsM[b]=b*binWidthM;


////  yield3d_systUp = new TH3D( this->getHistoName("yield3d_"+this->systName+"_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
////  yield3d_systUp->Sumw2();
////  yield3d_systDown = new TH3D( this->getHistoName("yield3d_"+this->systName+"_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
////  yield3d_systDown->Sumw2();

  yield3d_genmet = new TH3D( this->getHistoName("yield3d_genmet").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_genmet->Sumw2();

  yield3d_systUp = new TH3D( this->getHistoName("yield3d_"+asystName+"_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_systUp->Sumw2();
  yield3d_systDown = new TH3D( this->getHistoName("yield3d_"+asystName+"_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_systDown->Sumw2();

  yield3d_isr_UP = new TH3D( this->getHistoName("yield3d_isr_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_isr_UP->Sumw2();
  yield3d_isr_DN = new TH3D( this->getHistoName("yield3d_isr_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_isr_DN->Sumw2();

  yield3d_lepsf_UP = new TH3D( this->getHistoName("yield3d_lepsf_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_lepsf_UP->Sumw2();
  yield3d_lepsf_DN = new TH3D( this->getHistoName("yield3d_lepsf_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_lepsf_DN->Sumw2();

  yield3d_btag_light_UP = new TH3D( this->getHistoName("yield3d_btag_light_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_btag_light_UP->Sumw2();
  yield3d_btag_heavy_DN = new TH3D( this->getHistoName("yield3d_btag_heavy_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_btag_heavy_DN->Sumw2();

  yield3d_btag_light_DN = new TH3D( this->getHistoName("yield3d_btag_light_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_btag_light_DN->Sumw2();
  yield3d_btag_heavy_UP = new TH3D( this->getHistoName("yield3d_btag_heavy_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_btag_heavy_UP->Sumw2();

  yield3d_crsl = new TH3D( this->getHistoName("yield3d_crsl").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_crsl->Sumw2();

  yield3d_alpha = new TH3D( this->getHistoName("yield3d_alpha").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_alpha->Sumw2();



  //  int nBinsMY=161;
  //  double binWidthMY=5.;
  //  double binsMY[nBinsMY+1];
//  for (int b=0; b<=nBinsMY; ++b)
//    binsMY[b]=b*binWidthMY;
//
//  yield3d_systUp = new TH3D( this->getHistoName("yield3d_"+asystName+"_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsMY, binsMY);
//  yield3d_systUp->Sumw2();
//  yield3d_systDown = new TH3D( this->getHistoName("yield3d_"+asystName+"_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsMY, binsMY);
//  yield3d_systDown->Sumw2();

}






MT2EstimateAllSigSyst::MT2EstimateAllSigSyst( const MT2EstimateAllSigSyst& rhs ) : MT2Estimate(rhs) {

  this->yield3d_genmet = new TH3D(*(rhs.yield3d_genmet));

  this->yield3d_systUp = new TH3D(*(rhs.yield3d_systUp));
  this->yield3d_systDown = new TH3D(*(rhs.yield3d_systDown));

  this->yield3d_isr_UP = new TH3D(*(rhs.yield3d_isr_UP));
  this->yield3d_isr_DN = new TH3D(*(rhs.yield3d_isr_DN));

  this->yield3d_lepsf_UP = new TH3D(*(rhs.yield3d_lepsf_UP));
  this->yield3d_lepsf_DN = new TH3D(*(rhs.yield3d_lepsf_DN));

  this->yield3d_btag_light_UP = new TH3D(*(rhs.yield3d_btag_light_UP));
  this->yield3d_btag_light_DN = new TH3D(*(rhs.yield3d_btag_light_DN));
  this->yield3d_btag_heavy_UP = new TH3D(*(rhs.yield3d_btag_heavy_UP));
  this->yield3d_btag_heavy_DN = new TH3D(*(rhs.yield3d_btag_heavy_DN));

  this->yield3d_crsl  = new TH3D(*(rhs.yield3d_crsl));
  this->yield3d_alpha = new TH3D(*(rhs.yield3d_alpha));


}



MT2EstimateAllSigSyst::~MT2EstimateAllSigSyst() {

  delete yield3d_genmet;

  delete yield3d_systUp;
  delete yield3d_systDown;

  delete yield3d_systUp;
  delete yield3d_systDown;

  delete yield3d_isr_UP;
  delete yield3d_isr_DN;

  delete yield3d_lepsf_UP;
  delete yield3d_lepsf_DN;

  delete yield3d_btag_light_UP;
  delete yield3d_btag_light_DN;
  delete yield3d_btag_heavy_UP;
  delete yield3d_btag_heavy_DN;

  delete yield3d_crsl;
  delete yield3d_alpha;

}







void MT2EstimateAllSigSyst::setName( const std::string& newName ) {

  MT2Estimate::setName(newName);

  yield3d_genmet  ->SetName( this->getHistoName("yield3d_genmet").c_str() );

  yield3d_systUp  ->SetName( this->getHistoName("yield3d_"+this->systName+"_UP").c_str() );
  yield3d_systDown->SetName( this->getHistoName("yield3d_"+this->systName+"_DN").c_str() );

  yield3d_isr_UP  ->SetName( this->getHistoName("yield3d_isr_UP").c_str() );
  yield3d_isr_DN  ->SetName( this->getHistoName("yield3d_isr_DN").c_str() );

  yield3d_lepsf_UP  ->SetName( this->getHistoName("yield3d_lepsf_UP").c_str() );
  yield3d_lepsf_DN  ->SetName( this->getHistoName("yield3d_lepsf_DN").c_str() );

  yield3d_btag_heavy_UP  ->SetName( this->getHistoName("yield3d_btag_heavy_UP").c_str() );
  yield3d_btag_heavy_DN  ->SetName( this->getHistoName("yield3d_btag_heavy_DN").c_str() );

  yield3d_btag_light_UP  ->SetName( this->getHistoName("yield3d_btag_light_UP").c_str() );
  yield3d_btag_light_DN  ->SetName( this->getHistoName("yield3d_btag_light_DN").c_str() );


  yield3d_alpha  ->SetName( this->getHistoName("yield3d_alpha").c_str() );
  yield3d_crsl  ->SetName( this->getHistoName("yield3d_crsl").c_str() );


}


void MT2EstimateAllSigSyst::setSystName( const std::string& newSystName ) {

  yield3d_genmet  ->SetName( this->getHistoName("yield3d_genmet").c_str() );

  yield3d_systUp  ->SetName( this->getHistoName("yield3d_"+newSystName+"_UP").c_str() );
  yield3d_systDown->SetName( this->getHistoName("yield3d_"+newSystName+"_DN").c_str() );

  yield3d_isr_UP  ->SetName( this->getHistoName("yield3d_isr_UP").c_str() );
  yield3d_isr_DN  ->SetName( this->getHistoName("yield3d_isr_DN").c_str() );

  yield3d_lepsf_UP  ->SetName( this->getHistoName("yield3d_lepsf_UP").c_str() );
  yield3d_lepsf_DN  ->SetName( this->getHistoName("yield3d_lepsf_DN").c_str() );

  yield3d_btag_heavy_UP  ->SetName( this->getHistoName("yield3d_btag_heavy_UP").c_str() );
  yield3d_btag_heavy_DN  ->SetName( this->getHistoName("yield3d_btag_heavy_DN").c_str() );

  yield3d_btag_light_UP  ->SetName( this->getHistoName("yield3d_btag_light_UP").c_str() );
  yield3d_btag_light_DN  ->SetName( this->getHistoName("yield3d_btag_light_DN").c_str() );


  yield3d_alpha  ->SetName( this->getHistoName("yield3d_alpha").c_str() );
  yield3d_crsl  ->SetName( this->getHistoName("yield3d_crsl").c_str() );

}



void MT2EstimateAllSigSyst::finalize( ) {

  MT2Estimate::finalize();

}



void MT2EstimateAllSigSyst::getShit( TFile* file, const std::string& path ) {

  MT2Estimate::getShit(file, path);
  
  yield3d_genmet   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_genmet->GetName()));

  yield3d_systUp   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_systUp->GetName()));
  yield3d_systDown = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_systDown->GetName()));


  yield3d_isr_UP   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_isr_UP->GetName()));
  yield3d_isr_DN   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_isr_DN->GetName()));

  yield3d_lepsf_UP   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_lepsf_UP->GetName()));
  yield3d_lepsf_DN   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_lepsf_DN->GetName()));

  yield3d_btag_light_UP   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_btag_light_UP->GetName()));
  yield3d_btag_light_DN   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_btag_light_DN->GetName()));
  yield3d_btag_heavy_UP   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_btag_heavy_UP->GetName()));
  yield3d_btag_heavy_DN   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_btag_heavy_DN->GetName()));

  yield3d_crsl    = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_crsl->GetName()));
  yield3d_alpha   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_alpha->GetName()));

}


// void MT2EstimateAllSigSyst::print( std::ofstream& ofs_file, Float_t m1, Float_t m2, Int_t mt2_bin, float k, bool doGenAverage ){

//   int nBinsM=93;
//   double binWidthM=25.;
//   double binsM[nBinsM+1];
//   for (int b=0; b<=nBinsM; ++b)
//     binsM[b]=b*binWidthM;

//   int binY, binZ;

//   TH1D* h_sig;
//   TH3D* h_sig3d;

//   TH3D* h_sig3d_genmet;
//   TH1D* h_sig_genmet;

//   h_sig3d = this->yield3d;

//   if( h_sig3d == 0 ){

//     std::cout << "3d histogram does not exist, initializing empty histogram..." << std::endl;
//     h_sig = new TH1D("emptyHisto", "", nBinsM, binsM);

//   }
//   else{
//     binY = h_sig3d->GetYaxis()->FindBin(m1);
//     binZ = h_sig3d->GetZaxis()->FindBin(m2);

//     h_sig = h_sig3d->ProjectionX("mt2_0", binY, binY, binZ, binZ);

//     if(doGenAverage && this->yield3d_genmet != 0){
//       h_sig3d_genmet = this->yield3d_genmet;
//       h_sig_genmet = h_sig3d_genmet->ProjectionX("mt2_genmet", binY, binY, binZ, binZ);
//       h_sig->Add( h_sig_genmet );
//       h_sig->Scale(0.5);
//     }

//     std::cout << "Printing for m1 = " << m1 << ", m2 = " << m2 << std::endl;

//   }

//   Double_t error = h_sig->GetBinError(mt2_bin);
//   Double_t integral = h_sig->GetBinContent(mt2_bin);

//   if(integral<0){

//     integral=0.;
//     error=0;

//   }

//   integral*=k;
//   error*=k;

//   if(integral >= 10)
//     ofs_file << std::fixed << std::setprecision(1) << " & " << integral << " $\\pm$ " << error;
//   else if(integral < 10)
//     ofs_file << std::fixed << std::setprecision(2) << " & " << integral << " $\\pm$ " << error;

// }


void MT2EstimateAllSigSyst::write() const {

  MT2Estimate::write();
  yield3d_genmet->Write();
  yield3d_systUp->Write();
  yield3d_systDown->Write();

  yield3d_isr_UP->Write();
  yield3d_isr_DN->Write();

  yield3d_lepsf_UP->Write();
  yield3d_lepsf_DN->Write();

  yield3d_btag_light_UP->Write();
  yield3d_btag_light_DN->Write();
  yield3d_btag_heavy_UP->Write();
  yield3d_btag_heavy_DN->Write();

  yield3d_crsl->Write();
  yield3d_alpha->Write();

}





const MT2EstimateAllSigSyst& MT2EstimateAllSigSyst::operator=( const MT2EstimateAllSigSyst& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield3d_genmet   = new TH3D(*(rhs.yield3d_genmet));
  this->yield3d_systUp   = new TH3D(*(rhs.yield3d_systUp));
  this->yield3d_systDown = new TH3D(*(rhs.yield3d_systDown));


  this->yield3d_isr_UP   = new TH3D(*(rhs.yield3d_isr_UP));
  this->yield3d_isr_DN   = new TH3D(*(rhs.yield3d_isr_DN));

  this->yield3d_lepsf_UP   = new TH3D(*(rhs.yield3d_lepsf_UP));
  this->yield3d_lepsf_DN   = new TH3D(*(rhs.yield3d_lepsf_DN));

  this->yield3d_btag_light_UP   = new TH3D(*(rhs.yield3d_btag_light_UP));
  this->yield3d_btag_light_DN   = new TH3D(*(rhs.yield3d_btag_light_DN));
  this->yield3d_btag_heavy_UP   = new TH3D(*(rhs.yield3d_btag_heavy_UP));
  this->yield3d_btag_heavy_DN   = new TH3D(*(rhs.yield3d_btag_heavy_DN));

  this->yield3d_crsl    = new TH3D(*(rhs.yield3d_crsl));
  this->yield3d_alpha   = new TH3D(*(rhs.yield3d_alpha));


  this->systName = rhs.systName;

  this->setName( this->getName() );


  return *this;

}




const MT2EstimateAllSigSyst& MT2EstimateAllSigSyst::operator=( const MT2Estimate& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield3d_genmet   = new TH3D(*(rhs.yield3d));
  this->yield3d_systUp   = new TH3D(*(rhs.yield3d));
  this->yield3d_systDown = new TH3D(*(rhs.yield3d));

  this->yield3d_isr_UP   = new TH3D(*(rhs.yield3d));
  this->yield3d_isr_DN   = new TH3D(*(rhs.yield3d));

  this->yield3d_lepsf_UP   = new TH3D(*(rhs.yield3d));
  this->yield3d_lepsf_DN   = new TH3D(*(rhs.yield3d));

  this->yield3d_btag_light_UP   = new TH3D(*(rhs.yield3d));
  this->yield3d_btag_light_DN   = new TH3D(*(rhs.yield3d));
  this->yield3d_btag_heavy_UP   = new TH3D(*(rhs.yield3d));
  this->yield3d_btag_heavy_DN   = new TH3D(*(rhs.yield3d));

  this->yield3d_crsl    = new TH3D(*(rhs.yield3d));
  this->yield3d_alpha   = new TH3D(*(rhs.yield3d));

  this->systName = "syst";

  this->setName( this->getName() );


  return *this;

}

  



const MT2EstimateAllSigSyst& MT2EstimateAllSigSyst::operator+=( const MT2EstimateAllSigSyst& rhs ) {

  this->yield->Add(rhs.yield);
  this->yield3d->Add(rhs.yield3d);
  this->yield3d_systUp->Add(rhs.yield3d_systUp);
  this->yield3d_systDown->Add(rhs.yield3d_systDown);

  this->yield3d_genmet->Add(rhs.yield3d_genmet );

  this->yield3d_isr_UP->Add(rhs.yield3d_isr_UP);
  this->yield3d_isr_DN->Add(rhs.yield3d_isr_DN );
  this->yield3d_lepsf_UP->Add(rhs.yield3d_lepsf_UP );
  this->yield3d_lepsf_DN->Add(rhs.yield3d_lepsf_DN );
  this->yield3d_btag_light_UP->Add(rhs.yield3d_btag_light_UP );
  this->yield3d_btag_light_DN->Add(rhs.yield3d_btag_light_DN );

  this->yield3d_btag_heavy_UP->Add(rhs.yield3d_btag_heavy_UP );
  this->yield3d_btag_heavy_DN->Add(rhs.yield3d_btag_heavy_DN );

  this->yield3d_crsl->Add(rhs.yield3d_crsl);
  this->yield3d_alpha->Add(rhs.yield3d_alpha);

  return (*this);

}


MT2EstimateAllSigSyst MT2EstimateAllSigSyst::operator*( float k ) const{

  MT2EstimateAllSigSyst result(this->getName(), *(this->region), this->systName );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(k);
  
  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(k);

  result.yield3d_genmet = new TH3D(*(this->yield3d_genmet));
  result.yield3d_genmet->Scale(k);

  result.yield3d_systUp = new TH3D(*(this->yield3d_systUp));
  result.yield3d_systUp->Scale(k);

  result.yield3d_systDown = new TH3D(*(this->yield3d_systDown));
  result.yield3d_systDown->Scale(k);


  result.yield3d_isr_UP   = new TH3D(*(this->yield3d_isr_UP));
  result.yield3d_isr_UP->Scale(k);
  result.yield3d_isr_DN   = new TH3D(*(this->yield3d_isr_DN));
  result.yield3d_isr_DN->Scale(k);

  result.yield3d_lepsf_UP   = new TH3D(*(this->yield3d_lepsf_UP));
  result.yield3d_lepsf_UP->Scale(k);
  result.yield3d_lepsf_DN   = new TH3D(*(this->yield3d_lepsf_DN));
  result.yield3d_lepsf_DN->Scale(k);

  result.yield3d_btag_light_UP   = new TH3D(*(this->yield3d_btag_light_UP));
  result.yield3d_btag_light_UP->Scale(k);
  result.yield3d_btag_light_DN   = new TH3D(*(this->yield3d_btag_light_DN));
  result.yield3d_btag_light_DN->Scale(k);

  result.yield3d_btag_heavy_UP   = new TH3D(*(this->yield3d_btag_heavy_UP));
  result.yield3d_btag_heavy_UP->Scale(k);
  result.yield3d_btag_heavy_DN   = new TH3D(*(this->yield3d_btag_heavy_DN));
  result.yield3d_btag_heavy_DN->Scale(k);

  result.yield3d_crsl    = new TH3D(*(this->yield3d_crsl));  
  result.yield3d_crsl->Scale(k);
  result.yield3d_alpha   = new TH3D(*(this->yield3d_alpha)); 
  result.yield3d_crsl->Scale(k);


  return result;

}



MT2EstimateAllSigSyst MT2EstimateAllSigSyst::operator/( float k ) const{

  MT2EstimateAllSigSyst result(this->getName(), *(this->region), this->systName );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(1./k);

  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(1./k);

  result.yield3d_genmet = new TH3D(*(this->yield3d_genmet));
  result.yield3d_genmet->Scale(1./k);

  result.yield3d_systUp = new TH3D(*(this->yield3d_systUp));
  result.yield3d_systUp->Scale(1./k);

  result.yield3d_systDown = new TH3D(*(this->yield3d_systDown));
  result.yield3d_systDown->Scale(1./k);

  result.yield3d_isr_UP   = new TH3D(*(this->yield3d_isr_UP));
  result.yield3d_isr_UP->Scale(1./k);
  result.yield3d_isr_DN   = new TH3D(*(this->yield3d_isr_DN));
  result.yield3d_isr_DN->Scale(1./k);

  result.yield3d_lepsf_UP   = new TH3D(*(this->yield3d_lepsf_UP));
  result.yield3d_lepsf_UP->Scale(1./k);
  result.yield3d_lepsf_DN   = new TH3D(*(this->yield3d_lepsf_DN));
  result.yield3d_lepsf_DN->Scale(1./k);

  result.yield3d_btag_light_UP   = new TH3D(*(this->yield3d_btag_light_UP));
  result.yield3d_btag_light_UP->Scale(1./k);
  result.yield3d_btag_light_DN   = new TH3D(*(this->yield3d_btag_light_DN));
  result.yield3d_btag_light_DN->Scale(1./k);

  result.yield3d_btag_heavy_UP   = new TH3D(*(this->yield3d_btag_heavy_UP));
  result.yield3d_btag_heavy_UP->Scale(1./k);
  result.yield3d_btag_heavy_DN   = new TH3D(*(this->yield3d_btag_heavy_DN));
  result.yield3d_btag_heavy_DN->Scale(1./k);

  result.yield3d_crsl    = new TH3D(*(this->yield3d_crsl));  
  result.yield3d_crsl->Scale(1./k);
  result.yield3d_alpha   = new TH3D(*(this->yield3d_alpha)); 
  result.yield3d_crsl->Scale(1./k);

  return result;


}




const MT2EstimateAllSigSyst& MT2EstimateAllSigSyst::operator*=( float k ) {

  this->yield->Scale(k);
  this->yield3d->Scale(k);
  this->yield3d_genmet->Scale(k);
  this->yield3d_systUp->Scale(k);
  this->yield3d_systDown->Scale(k);

  this->yield3d_isr_UP->Scale(k);
  this->yield3d_isr_DN->Scale(k);

  this->yield3d_lepsf_UP->Scale(k);
  this->yield3d_lepsf_DN->Scale(k);

  this->yield3d_btag_light_UP->Scale(k);
  this->yield3d_btag_light_DN->Scale(k);
  this->yield3d_btag_heavy_UP->Scale(k);
  this->yield3d_btag_heavy_DN->Scale(k);

  this->yield3d_crsl->Scale(k);
  this->yield3d_alpha->Scale(k);

  return (*this);

}

const MT2EstimateAllSigSyst& MT2EstimateAllSigSyst::operator/=( float k ) {

  this->yield->Scale(1./k);
  this->yield3d->Scale(1./k);
  this->yield3d_genmet->Scale(1./k);
  this->yield3d_systUp->Scale(1./k);
  this->yield3d_systDown->Scale(1./k);

  this->yield3d_isr_UP->Scale(1./k);
  this->yield3d_isr_DN->Scale(1./k);

  this->yield3d_lepsf_UP->Scale(1./k);
  this->yield3d_lepsf_DN->Scale(1./k);

  this->yield3d_btag_light_UP->Scale(1./k);
  this->yield3d_btag_light_DN->Scale(1./k);
  this->yield3d_btag_heavy_UP->Scale(1./k);
  this->yield3d_btag_heavy_DN->Scale(1./k);

  this->yield3d_crsl->Scale(1./k);
  this->yield3d_alpha->Scale(1./k);

  return (*this);

}




// friend functions


MT2EstimateAllSigSyst operator*( float k, const MT2EstimateAllSigSyst& rhs ) {

  return rhs*k;

}


MT2EstimateAllSigSyst operator/( float k, const MT2EstimateAllSigSyst& rhs ) {

  return rhs/k;

}


