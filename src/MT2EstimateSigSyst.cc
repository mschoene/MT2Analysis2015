#include "../interface/MT2EstimateSigSyst.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TEfficiency.h"








MT2EstimateSigSyst::MT2EstimateSigSyst( const std::string& aname, const MT2Region& aregion, const std::string& asystName ) : MT2Estimate( aname, aregion ) {

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

  int nBinsM=81;
  double binWidthM=25.;
  double binsM[nBinsM+1];
  for (int b=0; b<=nBinsM; ++b)
    binsM[b]=b*binWidthM;

  int nBinsMY=161;
  double binWidthMY=5.;
  double binsMY[nBinsMY+1];
  for (int b=0; b<=nBinsMY; ++b)
    binsMY[b]=b*binWidthMY;

////  yield3d_systUp = new TH3D( this->getHistoName("yield3d_"+this->systName+"_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
////  yield3d_systUp->Sumw2();
////  yield3d_systDown = new TH3D( this->getHistoName("yield3d_"+this->systName+"_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
////  yield3d_systDown->Sumw2();

  yield3d_systUp = new TH3D( this->getHistoName("yield3d_"+asystName+"_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_systUp->Sumw2();
  yield3d_systDown = new TH3D( this->getHistoName("yield3d_"+asystName+"_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsM, binsM);
  yield3d_systDown->Sumw2();

//  yield3d_systUp = new TH3D( this->getHistoName("yield3d_"+asystName+"_UP").c_str(), "", nBins, bins, nBinsM, binsM, nBinsMY, binsMY);
//  yield3d_systUp->Sumw2();
//  yield3d_systDown = new TH3D( this->getHistoName("yield3d_"+asystName+"_DN").c_str(), "", nBins, bins, nBinsM, binsM, nBinsMY, binsMY);
//  yield3d_systDown->Sumw2();

}






MT2EstimateSigSyst::MT2EstimateSigSyst( const MT2EstimateSigSyst& rhs ) : MT2Estimate(rhs) {

  this->yield3d_systUp = new TH3D(*(rhs.yield3d_systUp));
  this->yield3d_systDown = new TH3D(*(rhs.yield3d_systDown));

}



MT2EstimateSigSyst::~MT2EstimateSigSyst() {

  delete yield3d_systUp;
  delete yield3d_systDown;

}







void MT2EstimateSigSyst::setName( const std::string& newName ) {

  MT2Estimate::setName(newName);

  yield3d_systUp  ->SetName( this->getHistoName("yield3d_"+this->systName+"_UP").c_str() );
  yield3d_systDown->SetName( this->getHistoName("yield3d_"+this->systName+"_DN").c_str() );

}


void MT2EstimateSigSyst::setSystName( const std::string& newSystName ) {

  yield3d_systUp  ->SetName( this->getHistoName("yield3d_"+newSystName+"_UP").c_str() );
  yield3d_systDown->SetName( this->getHistoName("yield3d_"+newSystName+"_DN").c_str() );

}



void MT2EstimateSigSyst::finalize( ) {

  MT2Estimate::finalize();

}



void MT2EstimateSigSyst::getShit( TFile* file, const std::string& path ) {

  MT2Estimate::getShit(file, path);
  
  yield3d_systUp   = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_systUp->GetName()));
  yield3d_systDown = (TH3D*)file->Get(Form("%s/%s", path.c_str(), yield3d_systDown->GetName()));


}



void MT2EstimateSigSyst::write() const {

  MT2Estimate::write();
  yield3d_systUp->Write();
  yield3d_systDown->Write();

}





const MT2EstimateSigSyst& MT2EstimateSigSyst::operator=( const MT2EstimateSigSyst& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield3d_systUp   = new TH3D(*(rhs.yield3d_systUp));
  this->yield3d_systDown = new TH3D(*(rhs.yield3d_systDown));

  this->systName = rhs.systName;

  this->setName( this->getName() );


  return *this;

}




const MT2EstimateSigSyst& MT2EstimateSigSyst::operator=( const MT2Estimate& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield3d_systUp   = new TH3D(*(rhs.yield3d));
  this->yield3d_systDown = new TH3D(*(rhs.yield3d));

  this->systName = "syst";

  this->setName( this->getName() );


  return *this;

}








MT2EstimateSigSyst MT2EstimateSigSyst::operator*( float k ) const{

  MT2EstimateSigSyst result(this->getName(), *(this->region), this->systName );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(k);
  
  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(k);

  result.yield3d_systUp = new TH3D(*(this->yield3d_systUp));
  result.yield3d_systUp->Scale(k);

  result.yield3d_systDown = new TH3D(*(this->yield3d_systDown));
  result.yield3d_systDown->Scale(k);

  return result;

}



MT2EstimateSigSyst MT2EstimateSigSyst::operator/( float k ) const{

  MT2EstimateSigSyst result(this->getName(), *(this->region), this->systName );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(1./k);

  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(1./k);

  result.yield3d_systUp = new TH3D(*(this->yield3d_systUp));
  result.yield3d_systUp->Scale(1./k);

  result.yield3d_systDown = new TH3D(*(this->yield3d_systDown));
  result.yield3d_systDown->Scale(1./k);

  return result;


}




const MT2EstimateSigSyst& MT2EstimateSigSyst::operator*=( float k ) {

  this->yield->Scale(k);
  this->yield3d->Scale(k);
  this->yield3d_systUp->Scale(k);
  this->yield3d_systDown->Scale(k);
  return (*this);

}

const MT2EstimateSigSyst& MT2EstimateSigSyst::operator/=( float k ) {

  this->yield->Scale(1./k);
  this->yield3d->Scale(1./k);
  this->yield3d_systUp->Scale(1./k);
  this->yield3d_systDown->Scale(1./k);
  return (*this);

}




// friend functions


MT2EstimateSigSyst operator*( float k, const MT2EstimateSigSyst& rhs ) {

  return rhs*k;

}


MT2EstimateSigSyst operator/( float k, const MT2EstimateSigSyst& rhs ) {

  return rhs/k;

}


