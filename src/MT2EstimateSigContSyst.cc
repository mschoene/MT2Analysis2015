#include "../interface/MT2EstimateSigContSyst.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TEfficiency.h"








MT2EstimateSigContSyst::MT2EstimateSigContSyst( const std::string& aname, const MT2Region& aregion, const std::string& asystName ) : MT2Estimate( aname, aregion ) {

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

  yield3d_systUp = new TH3D( this->getHistoName("yield3d_"+asystName+"_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_systUp->Sumw2();
  yield3d_systDown = new TH3D( this->getHistoName("yield3d_"+asystName+"_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_systDown->Sumw2();

  yield_alpha = new TH1D( this->getHistoName("yield_alpha").c_str(), "", nBins, bins);
  yield_alpha->Sumw2();

  yield3d_crsl = new TH3D( this->getHistoName("yield3d_crsl").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_crsl->Sumw2();

  yield3d_genmet = new TH3D( this->getHistoName("yield3d_genmet").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_genmet->Sumw2();

  yield3d_crsl_genmet = new TH3D( this->getHistoName("yield3d_crsl_genmet").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_crsl_genmet->Sumw2();

}






MT2EstimateSigContSyst::MT2EstimateSigContSyst( const MT2EstimateSigContSyst& rhs ) : MT2Estimate(rhs) {

  this->yield3d_systUp = new TH3D(*(rhs.yield3d_systUp));
  this->yield3d_systDown = new TH3D(*(rhs.yield3d_systDown));

  this->yield_alpha = new TH1D(*(rhs.yield_alpha));
  this->yield3d_crsl = new TH3D(*(rhs.yield3d_crsl));

  this->yield3d_genmet = new TH3D(*(rhs.yield3d_genmet));
  this->yield3d_crsl_genmet = new TH3D(*(rhs.yield3d_crsl_genmet));

}



MT2EstimateSigContSyst::~MT2EstimateSigContSyst() {

  delete yield3d_systUp;
  delete yield3d_systDown;

  delete yield_alpha;
  delete yield3d_crsl;

  delete yield3d_genmet;
  delete yield3d_crsl_genmet;

}







void MT2EstimateSigContSyst::setName( const std::string& newName ) {

  MT2Estimate::setName(newName);

  yield3d_systUp  ->SetName( this->getHistoName("yield3d_"+this->systName+"_UP").c_str() );
  yield3d_systDown->SetName( this->getHistoName("yield3d_"+this->systName+"_DN").c_str() );

  yield3d_crsl  ->SetName( this->getHistoName("yield3d_crsl").c_str() );
  yield_alpha->SetName( this->getHistoName("yield_alpha").c_str() );

  yield3d_genmet  ->SetName( this->getHistoName("yield3d_genmet").c_str() );
  yield3d_crsl_genmet  ->SetName( this->getHistoName("yield3d_crsl_genmet").c_str() );

}


void MT2EstimateSigContSyst::setSystName( const std::string& newSystName ) {

  yield3d_systUp  ->SetName( this->getHistoName("yield3d_"+newSystName+"_UP").c_str() );
  yield3d_systDown->SetName( this->getHistoName("yield3d_"+newSystName+"_DN").c_str() );

}



void MT2EstimateSigContSyst::finalize( ) {

  MT2Estimate::finalize();

}



void MT2EstimateSigContSyst::getShit( TFile* file, const std::string& path ) {

  MT2Estimate::getShit(file, path);
  
  yield3d_systUp   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_systUp->GetName()));
  yield3d_systDown = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_systDown->GetName()));

  yield3d_crsl   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_crsl->GetName()));
  yield_alpha = (TH1D*)file->Get(Form("%s/%s", path.c_str(), yield_alpha->GetName()));

  yield3d_genmet   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_genmet->GetName()));
  yield3d_crsl_genmet   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_crsl_genmet->GetName()));

}



void MT2EstimateSigContSyst::write() const {

  MT2Estimate::write();
  yield3d_systUp->Write();
  yield3d_systDown->Write();

  yield3d_crsl->Write();
  yield_alpha->Write();

  yield3d_genmet->Write();
  yield3d_crsl_genmet->Write();

}

void MT2EstimateSigContSyst::print( std::ofstream& ofs_file, Float_t m1, Float_t m2, Int_t mt2_bin, float k, bool doGenAverage ){

  TH1D* h_sig;
  TH3D* h_sig3d;
  h_sig3d = this->yield3d;
  
  if(doGenAverage){
    h_sig3d->Add( (this->yield3d_genmet) );
    h_sig3d->Scale(0.5);
  }

  TH1D* h_sigCont;
  TH3D* h_sig3d_crsl;
  h_sig3d_crsl = this->yield3d_crsl;
  
  if(doGenAverage){
    h_sig3d_crsl->Add( (this->yield3d_crsl_genmet) );
    h_sig3d_crsl->Scale(0.5);
  }

  TH1D* h_alpha;
  h_alpha = this->yield_alpha;

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

  std::cout << "Printing for m1 = " << m1 << ", m2 = " << m2 << std::endl;

  if( h_sig3d == 0 ){

    std::cout << "3d histogram does not exist, initializing empty histogram..." << std::endl;
    h_sig3d = new TH3D("emptyHisto", "", nBinsM, binsM, nBinsM, binsM, nBinsM, binsM);

    h_sig3d_crsl = new TH3D("emptyHisto_crsl", "", nBinsM, binsM, nBinsM, binsM, nBinsM, binsM);
    
  }

  int binY, binZ;
  binY = h_sig3d->GetYaxis()->FindBin(m1);
  binZ = h_sig3d->GetYaxis()->FindBin(m2);
  h_sig = h_sig3d->ProjectionX("mt2_0", binY, binY, binZ, binZ);

  h_sigCont = h_sig3d_crsl->ProjectionX("mt2_crsl", binY, binY, binZ, binZ);
  h_sigCont->Multiply(h_alpha);  
  h_sig->Add(h_sigCont, -1);
  
  Double_t error;
  Double_t integral = h_sig->IntegralAndError(mt2_bin, mt2_bin, error);
  
  if(integral<0){
    
    integral=0.;
    error=0;
    
  }
  
  integral*=k;
  error*=k;
  
  if(integral >= 10)
    ofs_file << std::fixed << std::setprecision(1) << " & " << integral << " $\\pm$ " << error;
  else if(integral < 10)
    ofs_file << std::fixed << std::setprecision(2) << " & " << integral << " $\\pm$ " << error;

}



const MT2EstimateSigContSyst& MT2EstimateSigContSyst::operator=( const MT2EstimateSigContSyst& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield3d_systUp   = new TH3D(*(rhs.yield3d_systUp));
  this->yield3d_systDown = new TH3D(*(rhs.yield3d_systDown));
  this->yield3d_crsl   = new TH3D(*(rhs.yield3d_crsl));
  this->yield_alpha = new TH1D(*(rhs.yield_alpha));

  this->yield3d_genmet   = new TH3D(*(rhs.yield3d_genmet));
  this->yield3d_crsl_genmet   = new TH3D(*(rhs.yield3d_crsl_genmet));

  this->systName = rhs.systName;

  this->setName( this->getName() );


  return *this;

}




const MT2EstimateSigContSyst& MT2EstimateSigContSyst::operator=( const MT2Estimate& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield3d_systUp   = new TH3D(*(rhs.yield3d));
  this->yield3d_systDown = new TH3D(*(rhs.yield3d));
  this->yield_alpha = new TH1D(*(rhs.yield));
  this->yield3d_crsl = new TH3D(*(rhs.yield3d));

  this->yield3d_genmet = new TH3D(*(rhs.yield3d));
  this->yield3d_crsl_genmet = new TH3D(*(rhs.yield3d));

  this->systName = "syst";

  this->setName( this->getName() );


  return *this;

}








MT2EstimateSigContSyst MT2EstimateSigContSyst::operator*( float k ) const{

  MT2EstimateSigContSyst result(this->getName(), *(this->region), this->systName );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(k);
  
  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(k);

  result.yield3d_systUp = new TH3D(*(this->yield3d_systUp));
  result.yield3d_systUp->Scale(k);

  result.yield3d_systDown = new TH3D(*(this->yield3d_systDown));
  result.yield3d_systDown->Scale(k);

  result.yield_alpha = new TH1D(*(this->yield_alpha));
  result.yield_alpha->Scale(k);
  
  result.yield3d_crsl = new TH3D(*(this->yield3d_crsl));
  result.yield3d_crsl->Scale(k);

  result.yield3d_genmet = new TH3D(*(this->yield3d_genmet));
  result.yield3d_genmet->Scale(k);

  result.yield3d_crsl_genmet = new TH3D(*(this->yield3d_crsl_genmet));
  result.yield3d_crsl_genmet->Scale(k);

  return result;

}



MT2EstimateSigContSyst MT2EstimateSigContSyst::operator/( float k ) const{

  MT2EstimateSigContSyst result(this->getName(), *(this->region), this->systName );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(1./k);

  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(1./k);

  result.yield3d_systUp = new TH3D(*(this->yield3d_systUp));
  result.yield3d_systUp->Scale(1./k);

  result.yield3d_systDown = new TH3D(*(this->yield3d_systDown));
  result.yield3d_systDown->Scale(1./k);

  result.yield_alpha = new TH1D(*(this->yield_alpha));
  result.yield_alpha->Scale(1./k);
  
  result.yield3d_crsl = new TH3D(*(this->yield3d_crsl));
  result.yield3d_crsl->Scale(1./k);

  result.yield3d_genmet = new TH3D(*(this->yield3d_genmet));
  result.yield3d_genmet->Scale(1./k);

  result.yield3d_crsl_genmet = new TH3D(*(this->yield3d_crsl_genmet));
  result.yield3d_crsl_genmet->Scale(1./k);


  return result;


}




const MT2EstimateSigContSyst& MT2EstimateSigContSyst::operator*=( float k ) {

  this->yield->Scale(k);
  this->yield3d->Scale(k);
  this->yield3d_systUp->Scale(k);
  this->yield3d_systDown->Scale(k);
  this->yield_alpha->Scale(k);
  this->yield3d_crsl->Scale(k);

  this->yield3d_genmet->Scale(k);
  this->yield3d_crsl_genmet->Scale(k);

  return (*this);

}

const MT2EstimateSigContSyst& MT2EstimateSigContSyst::operator/=( float k ) {

  this->yield->Scale(1./k);
  this->yield3d->Scale(1./k);
  this->yield3d_systUp->Scale(1./k);
  this->yield3d_systDown->Scale(1./k);
  this->yield_alpha->Scale(1./k);
  this->yield3d_crsl->Scale(1./k);

  this->yield3d_genmet->Scale(1./k);
  this->yield3d_crsl_genmet->Scale(1./k);

  return (*this);

}




// friend functions


MT2EstimateSigContSyst operator*( float k, const MT2EstimateSigContSyst& rhs ) {

  return rhs*k;

}


MT2EstimateSigContSyst operator/( float k, const MT2EstimateSigContSyst& rhs ) {

  return rhs/k;

}


