#include "../interface/MT2Estimate.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TRandom3.h"



MT2Estimate::MT2Estimate( const std::string& aname, const MT2Region& aregion ) {

  name = aname;

  region = new MT2Region(aregion);

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

  yield = new TH1D(this->getHistoName("yield").c_str(), "", nBins, bins);
  yield->Sumw2();

}


MT2Estimate::MT2Estimate( const MT2Estimate& rhs ) {

  name = rhs.getName();

  region = new MT2Region(*(rhs.region));

  yield = new TH1D(*(rhs.yield));

}



MT2Estimate::~MT2Estimate() {

  delete region;
  delete yield;

}



std::string MT2Estimate::getHistoName( const std::string& prefix ) const {

  std::string returnName = prefix + "_" + name + "_" + region->getName();
  
  return returnName;

}



void MT2Estimate::setName( const std::string& newName ) {

  name = newName;
  yield->SetName( this->getHistoName("yield").c_str() );

}




void MT2Estimate::addOverflow() {

  MT2Estimate::addOverflowSingleHisto( yield );

}


void MT2Estimate::addOverflowSingleHisto( TH1D* yield ) {

  yield->SetBinContent(yield->GetNbinsX(),
      yield->GetBinContent(yield->GetNbinsX()  )+
      yield->GetBinContent(yield->GetNbinsX()+1)  );
  yield->SetBinError(  yield->GetNbinsX(),
      sqrt(yield->GetBinError(yield->GetNbinsX()  )*
           yield->GetBinError(yield->GetNbinsX()  )+
           yield->GetBinError(yield->GetNbinsX()+1)*
           yield->GetBinError(yield->GetNbinsX()+1)  ));

  yield->SetBinContent(yield->GetNbinsX()+1, 0.);
  yield->SetBinError  (yield->GetNbinsX()+1, 0.);

}



const MT2Estimate& MT2Estimate::operator=( const MT2Estimate& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));

  this->setName(this->getName());

  return *this;

}


MT2Estimate MT2Estimate::operator+( const MT2Estimate& rhs ) const {


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2Estimate::operator+] ERROR! Can't add MT2Estimate with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2Estimate result(*this);
  result.yield->Add(rhs.yield);

  return result;

}



MT2Estimate MT2Estimate::operator-( const MT2Estimate& rhs ) const {


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2Estimate::operator-] ERROR! Can't add MT2Estimate with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2Estimate result(*this);
  result.yield->Add(rhs.yield, -1.);

  return result;

}




MT2Estimate MT2Estimate::operator/( const MT2Estimate& rhs ) const {


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2Estimate::operator/] ERROR! Can't divide MT2Estimate with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2Estimate result(*this);
  result.yield->Divide(rhs.yield);
  //MT2Estimate result(name, *(this->region) );
  //result.yield = new TH1D(*(this->yield));
  //result.yield->Divide(rhs.yield);

  return result;

}


MT2Estimate MT2Estimate::operator*( const MT2Estimate& rhs ) const {

  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2Estimate::operator*] ERROR! Can't multiply MT2Estimate with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2Estimate result(*this);
  result.yield->Multiply(rhs.yield);
  //MT2Estimate result(name, *(this->region) );
  //result.yield = new TH1D(*(this->yield));
  //result.yield->Multiply(rhs.yield);

  return result;

}



MT2Estimate MT2Estimate::operator/( float k ) const {

  MT2Estimate result(*this);
  result.yield->Scale(1./k);
  //MT2Estimate result(name, *(this->region) );
  //result.yield = new TH1D(*(this->yield));
  //result.yield->Scale(1./k);

  return result;

}


MT2Estimate MT2Estimate::operator*( float k ) const {

  MT2Estimate result(*this);
  result.yield->Scale(k);
  //MT2Estimate result(name, *(this->region) );
  //result.yield = new TH1D(*(this->yield));
  //result.yield->Scale(k);

  return result;

}




const MT2Estimate& MT2Estimate::operator/=( const MT2Estimate& rhs ) {

  this->yield->Divide(rhs.yield);
  return (*this);

}



const MT2Estimate& MT2Estimate::operator+=( const MT2Estimate& rhs ) {

  this->yield->Add(rhs.yield);
  return (*this);

}


const MT2Estimate& MT2Estimate::operator-=( const MT2Estimate& rhs ) {

  this->yield->Add(rhs.yield, -1.);
  return (*this);

}


const MT2Estimate& MT2Estimate::operator*=( const MT2Estimate& rhs ) {

  this->yield->Multiply(rhs.yield);
  return (*this);

}


const MT2Estimate& MT2Estimate::operator*=( float k ) {

  this->yield->Scale(k);
  return (*this);

}


const MT2Estimate& MT2Estimate::operator/=( float k ) {

  this->yield->Scale(1./k);
  return (*this);

}






void MT2Estimate::getShit( TFile* file, const std::string& path ) {

  yield = (TH1D*)file->Get(Form("%s/%s", path.c_str(), yield->GetName()));

}



void MT2Estimate::print(const std::string& ofs){

  Int_t binXmin=1;
  Int_t binXmax=5;

  Double_t error;
  Double_t integral = yield->IntegralAndError(binXmin, binXmax, error);

  ofstream ofs_file;
  ofs_file.open( ofs, std::ofstream::app );
  if(integral >= 10)
    ofs_file << std::fixed << std::setprecision(1) << " & " << integral << " $\\pm$ " << error;
  else if(integral < 10)
    ofs_file << std::fixed << std::setprecision(2) << " & " << integral << " $\\pm$ " << error;

}


void MT2Estimate::randomizePoisson( float scale ){

  TRandom3 rand(13);
  
  for( int ibin=1; ibin<yield->GetXaxis()->GetNbins()+1; ++ibin ) {
    
    int poisson_data = rand.Poisson(scale * yield->GetBinContent(ibin));
    yield->SetBinContent(ibin, poisson_data);
    yield->SetBinError(ibin, 0.); // it's data 
    
  } 
  
}


// friend functions

MT2Estimate operator*( float k, const MT2Estimate& rhs ) {

  return rhs*k;

}


MT2Estimate operator/( float k, const MT2Estimate& rhs ) {

  return rhs/k;

}






