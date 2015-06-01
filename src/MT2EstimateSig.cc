#include "../interface/MT2EstimateSig.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TRandom3.h"



MT2EstimateSig::MT2EstimateSig( const std::string& aname, const MT2Region& aregion ) {

  name = aname;

  region = new MT2Region(aregion);

  int nBins;
  double* bins;
  region->getBins(nBins, bins);
  
  int nBinsM=60;
  double binsM[nBinsM+1];
  for (int b=0; b<=nBinsM; ++b)
    binsM[b]=b*25.;
 
  yield3d = new TH3D(this->getHistoName("yield3d").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d->Sumw2();

  yield = new TH1D(this->getHistoName("yield").c_str(), "", nBins, bins);
  yield->Sumw2();

}


MT2EstimateSig::MT2EstimateSig( const MT2EstimateSig& rhs ) {

  name = rhs.getName();

  region = new MT2Region(*(rhs.region));

  yield3d = new TH3D(*(rhs.yield3d));
  yield = new TH1D(*(rhs.yield));

}



MT2EstimateSig::~MT2EstimateSig() {

  delete region;
  delete yield;
  delete yield3d;

  
}



std::string MT2EstimateSig::getHistoName( const std::string& prefix ) const {

  std::string returnName = prefix + "_" + name + "_" + region->getName();
  
  return returnName;

}



void MT2EstimateSig::setName( const std::string& newName ) {

  name = newName;
  yield3d->SetName( this->getHistoName("yield3d").c_str() );
  yield->SetName( this->getHistoName("yield").c_str() );

}




void MT2EstimateSig::addOverflow() {

  MT2EstimateSig::addOverflowSingleHisto( yield3d );
  MT2EstimateSig::addOverflowSingleHisto( yield );

}


void MT2EstimateSig::addOverflowSingleHisto( TH3D* yield3d ) {
  
  for (int y=1; y<=yield3d->GetNbinsY()+1; ++y)
    for (int z=1; z<=yield3d->GetNbinsZ()+1; ++z){

      yield3d->SetBinContent(yield3d->GetNbinsX(), y, z,
			   yield3d->GetBinContent(yield3d->GetNbinsX(), y, z  )+
			   yield3d->GetBinContent(yield3d->GetNbinsX()+1, y, z)  );
      yield3d->SetBinError(  yield3d->GetNbinsX(), y, z,
			   sqrt(yield3d->GetBinError(yield3d->GetNbinsX(), y, z  )*
				yield3d->GetBinError(yield3d->GetNbinsX(), y, z  )+
				yield3d->GetBinError(yield3d->GetNbinsX()+1, y, z)*
				yield3d->GetBinError(yield3d->GetNbinsX()+1, y, z)  ));
      
      yield3d->SetBinContent(yield3d->GetNbinsX()+1, y, z, 0.);
      yield3d->SetBinError  (yield3d->GetNbinsX()+1, y, z, 0.);
    }
  
}

void MT2EstimateSig::addOverflowSingleHisto( TH1D* yield ) {
  
  yield->SetBinContent(yield->GetNbinsX(),
			 yield->GetBinContent(yield->GetNbinsX()  )+
			 yield->GetBinContent(yield->GetNbinsX()+1)  );
  yield->SetBinError(  yield->GetNbinsX(),
			 sqrt(yield->GetBinError(yield->GetNbinsX() )*
			      yield->GetBinError(yield->GetNbinsX() )+
			      yield->GetBinError(yield->GetNbinsX()+1)*
			      yield->GetBinError(yield->GetNbinsX()+1)  ));
  
  yield->SetBinContent(yield->GetNbinsX()+1, 0.);
  yield->SetBinError  (yield->GetNbinsX()+1, 0.);
      
}


const MT2EstimateSig& MT2EstimateSig::operator=( const MT2EstimateSig& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield = new TH1D(*(rhs.yield));
  
  this->setName(this->getName());

  return *this;

}


MT2EstimateSig MT2EstimateSig::operator+( const MT2EstimateSig& rhs ) const {


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSig::operator+] ERROR! Can't add MT2EstimateSig with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSig result(*this);
  result.yield3d->Add(rhs.yield3d);
  result.yield->Add(rhs.yield);

  return result;

}



MT2EstimateSig MT2EstimateSig::operator-( const MT2EstimateSig& rhs ) const {


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSig::operator-] ERROR! Can't add MT2EstimateSig with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSig result(*this);
  result.yield3d->Add(rhs.yield3d, -1.);
  result.yield->Add(rhs.yield, -1.);

  return result;

}




MT2EstimateSig MT2EstimateSig::operator/( const MT2EstimateSig& rhs ) const {


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSig::operator/] ERROR! Can't divide MT2EstimateSig with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateSig result(*this);
  result.yield3d->Divide(rhs.yield3d);
  result.yield->Divide(rhs.yield);
  //MT2EstimateSig result(name, *(this->region) );
  //result.yield = new TH1D(*(this->yield));
  //result.yield->Divide(rhs.yield);

  return result;

}


MT2EstimateSig MT2EstimateSig::operator*( const MT2EstimateSig& rhs ) const {

  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSig::operator*] ERROR! Can't multiply MT2EstimateSig with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSig result(*this);
  result.yield3d->Multiply(rhs.yield3d);
  result.yield->Multiply(rhs.yield);
  //MT2EstimateSig result(name, *(this->region) );
  //result.yield = new TH1D(*(this->yield));
  //result.yield->Multiply(rhs.yield);

  return result;

}



MT2EstimateSig MT2EstimateSig::operator/( float k ) const {

  MT2EstimateSig result(*this);
  result.yield3d->Scale(1./k);
  result.yield->Scale(1./k);
  //MT2EstimateSig result(name, *(this->region) );
  //result.yield = new TH1D(*(this->yield));
  //result.yield->Scale(1./k);

  return result;

}


MT2EstimateSig MT2EstimateSig::operator*( float k ) const {

  MT2EstimateSig result(*this);
  result.yield3d->Scale(k);
  result.yield->Scale(k);
  //MT2EstimateSig result(name, *(this->region) );
  //result.yield = new TH1D(*(this->yield));
  //result.yield->Scale(k);

  return result;

}




const MT2EstimateSig& MT2EstimateSig::operator/=( const MT2EstimateSig& rhs ) {

  this->yield3d->Divide(rhs.yield3d);
  this->yield->Divide(rhs.yield);
  return (*this);

}



const MT2EstimateSig& MT2EstimateSig::operator+=( const MT2EstimateSig& rhs ) {

  this->yield3d->Add(rhs.yield3d);
  this->yield->Add(rhs.yield);
  return (*this);

}


const MT2EstimateSig& MT2EstimateSig::operator-=( const MT2EstimateSig& rhs ) {

  this->yield3d->Add(rhs.yield3d, -1.);
  this->yield->Add(rhs.yield, -1.);
  return (*this);

}


const MT2EstimateSig& MT2EstimateSig::operator*=( const MT2EstimateSig& rhs ) {

  this->yield3d->Multiply(rhs.yield3d);
  this->yield->Multiply(rhs.yield);
  return (*this);

}


const MT2EstimateSig& MT2EstimateSig::operator*=( float k ) {

  this->yield3d->Scale(k);
  this->yield->Scale(k);
  return (*this);

}


const MT2EstimateSig& MT2EstimateSig::operator/=( float k ) {

  this->yield3d->Scale(1./k);
  this->yield->Scale(1./k);
  return (*this);

}






void MT2EstimateSig::getShit( TFile* file, const std::string& path ) {

  yield3d = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d->GetName()));
  yield = (TH1D*)file->Get(Form("%s/%s", path.c_str(), yield->GetName()));

}



void MT2EstimateSig::print(const std::string& ofs){

  Int_t binXmin=1;
  Int_t binXmax=-1;

  Double_t error;
  Double_t integral = yield->IntegralAndError(binXmin, binXmax, error);

  ofstream ofs_file;
  ofs_file.open( ofs, std::ofstream::app );
  if(integral >= 10)
    ofs_file << std::fixed << std::setprecision(1) << " & " << integral << " $\\pm$ " << error;
  else if(integral < 10)
    ofs_file << std::fixed << std::setprecision(2) << " & " << integral << " $\\pm$ " << error;

}


void MT2EstimateSig::randomizePoisson( float scale ){

  TRandom3 rand(13);
  
  for( int ibin=1; ibin<yield->GetXaxis()->GetNbins()+1; ++ibin ) {
    
    int poisson_data = rand.Poisson(scale * yield->GetBinContent(ibin));
    yield->SetBinContent(ibin, poisson_data);
    yield->SetBinError(ibin, 0.); // it's data 
    
  } 
  
}

void MT2EstimateSig::fillYield( float mt2, int m1, int m2, float weight){

  yield->Fill(mt2, weight);
  yield3d->Fill(mt2, m1, m2, weight);

}


// friend functions

MT2EstimateSig operator*( float k, const MT2EstimateSig& rhs ) {

  return rhs*k;

}


MT2EstimateSig operator/( float k, const MT2EstimateSig& rhs ) {

  return rhs/k;

}






