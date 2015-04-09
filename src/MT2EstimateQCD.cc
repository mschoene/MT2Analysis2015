#include "../interface/MT2EstimateQCD.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TRandom3.h"
#include "TMath.h"


MT2EstimateQCD::MT2EstimateQCD( const std::string& aname, const MT2Region& aregion ) : MT2Estimate( aname, aregion ) {

  int nBins;
  double* bins;
  region->getBins_qcdCR(nBins, bins);

  lDphi = new TH1D( this->getHistoName("lDphi").c_str(), "", nBins, bins); lDphi->Sumw2();
  hDphi = new TH1D( this->getHistoName("hDphi").c_str(), "", nBins, bins); hDphi->Sumw2();
  ratio = new TH1D( this->getHistoName("ratio").c_str(), "", nBins, bins); ratio->Sumw2();

  exp      = new TF1(this->getHistoName("exp"     ).c_str(), "expo(0)"                                   , bins[0], bins[nBins]);
  expPlusC = new TF1(this->getHistoName("expPlusC").c_str(), "expo(0)+[2]"                               , bins[0], bins[nBins]);
  expOrC   = new TF1(this->getHistoName("expOrC"  ).c_str(), "(x>=200)*exp([0]+200.*[1])+(x<200)*expo(0)", bins[0], bins[nBins]);

  fitXmin  = 50.;
  fitXmax  = 80.;
  dphi_low = 0.3;

}

MT2EstimateQCD::MT2EstimateQCD( const MT2EstimateQCD& rhs ) : MT2Estimate(rhs) {

  this->lDphi = new TH1D(*(rhs.lDphi));
  this->hDphi = new TH1D(*(rhs.hDphi));
  this->ratio = new TH1D(*(rhs.ratio));

  this->exp      = new TF1(*(rhs.exp));
  this->expPlusC = new TF1(*(rhs.expPlusC));
  this->expOrC   = new TF1(*(rhs.expOrC));

  this->fitXmin  = rhs.fitXmin;
  this->fitXmax  = rhs.fitXmax;
  this->dphi_low = rhs.dphi_low;

}

// destructor
MT2EstimateQCD::~MT2EstimateQCD() {

  delete lDphi;
  delete hDphi;
  delete ratio;

  delete exp;
  delete expPlusC;
  delete expOrC;

}

void MT2EstimateQCD::setName( const std::string& newName ) {

  MT2Estimate::setName(newName);

  lDphi->SetName( this->getHistoName("lDphi").c_str() );
  hDphi->SetName( this->getHistoName("hDphi").c_str() );
  ratio->SetName( this->getHistoName("ratio").c_str() );

  exp     ->SetName( this->getHistoName("exp"     ).c_str() );
  expPlusC->SetName( this->getHistoName("expPlusC").c_str() );
  expOrC  ->SetName( this->getHistoName("expOrC"  ).c_str() );

}

void MT2EstimateQCD::fillDphi(float dphi, float weight, float mt2) {

  if (dphi > 0.3)
    this->hDphi->Fill( mt2, weight );
  else if (dphi < dphi_low)
    this->lDphi->Fill( mt2, weight );
      
}

void MT2EstimateQCD::doFit() {

  ratio->Fit(exp     , "Q0", "", fitXmin, fitXmax             ); // fit to exponential only
  ratio->Fit(expPlusC, "Q0", "", fitXmin, expPlusC->GetXmax() ); // fit to exponential + constant
  expOrC->SetParameter(0, exp->GetParameter(0));
  expOrC->SetParameter(1, exp->GetParameter(1));

}

void MT2EstimateQCD::finalize() {

  MT2Estimate::finalize();

  MT2Estimate::addOverflowSingleHisto( lDphi );
  MT2Estimate::addOverflowSingleHisto( hDphi );

  getRatio();
  doFit   ();

}

void MT2EstimateQCD::getShit( TFile* file, const std::string& path ) {

  MT2Estimate::getShit(file, path);

  lDphi = (TH1D*)file->Get(Form("%s/%s", path.c_str(), lDphi->GetName()));
  hDphi = (TH1D*)file->Get(Form("%s/%s", path.c_str(), hDphi->GetName()));
  ratio = (TH1D*)file->Get(Form("%s/%s", path.c_str(), ratio->GetName()));

  exp      = (TF1*)file->Get(Form("%s/%s", path.c_str(), exp     ->GetName()));
  expPlusC = (TF1*)file->Get(Form("%s/%s", path.c_str(), expPlusC->GetName()));
  expOrC   = (TF1*)file->Get(Form("%s/%s", path.c_str(), expOrC  ->GetName()));

}

void MT2EstimateQCD::print(const std::string& ofs){

  MT2Estimate::print( ofs );

}

void MT2EstimateQCD::randomizePoisson( float scale ){

  MT2Estimate::randomizePoisson( scale );

  TRandom3 rand(13);
  
  for( int ibin=1; ibin<lDphi->GetXaxis()->GetNbins()+1; ++ibin ) {
    
    int poisson_data = rand.Poisson(scale * lDphi->GetBinContent(ibin));
    lDphi->SetBinContent(ibin, poisson_data);
    lDphi->SetBinError( ibin, TMath::Sqrt(poisson_data) ); // here i want an approximation of the Poisson error
    
  } 

  for( int ibin=1; ibin<hDphi->GetXaxis()->GetNbins()+1; ++ibin ) {
    
    int poisson_data = rand.Poisson(scale * hDphi->GetBinContent(ibin));
    hDphi->SetBinContent(ibin, poisson_data);
    hDphi->SetBinError( ibin, TMath::Sqrt(poisson_data) ); // here i want an approximation of the Poisson error
    
  } 
  
}


void MT2EstimateQCD::write() const {

  MT2Estimate::write();

  lDphi->Write();
  hDphi->Write();
  ratio->Write();

  exp     ->Write();
  expPlusC->Write();
  expOrC  ->Write();

}

const MT2EstimateQCD& MT2EstimateQCD::operator=( const MT2EstimateQCD& rhs ) {

  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));

  this->lDphi = new TH1D(*(rhs.lDphi));
  this->hDphi = new TH1D(*(rhs.hDphi));
  this->ratio = new TH1D(*(rhs.ratio));

  this->exp      = new TF1(*(rhs.exp     ));
  this->expPlusC = new TF1(*(rhs.expPlusC));
  this->expOrC   = new TF1(*(rhs.expOrC  ));

  this->fitXmin  = rhs.fitXmin ;
  this->fitXmax  = rhs.fitXmax ;
  this->dphi_low = rhs.dphi_low;

  this->setName(this->getName());

  return *this;

}

MT2EstimateQCD MT2EstimateQCD::operator+( const MT2EstimateQCD& rhs ) const{

  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateQCD::operator+] ERROR! Can't add MT2EstimateQCD with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateQCD result(*this);
  result.yield->Add(rhs.yield);

  result.lDphi->Add(rhs.lDphi);
  result.hDphi->Add(rhs.hDphi);
  
  return result;

}

MT2EstimateQCD MT2EstimateQCD::operator-( const MT2EstimateQCD& rhs ) const{

  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateQCD::operator-] ERROR! Can't add MT2EstimateQCD with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateQCD result(*this);

  result.yield->Add(rhs.yield, -1.);
  result.lDphi->Add(rhs.lDphi, -1.);
  result.hDphi->Add(rhs.hDphi, -1.);
 
  return result;

}

const MT2EstimateQCD& MT2EstimateQCD::operator+=( const MT2EstimateQCD& rhs ) {

  this->yield->Add(rhs.yield);

  this->lDphi->Add(rhs.lDphi);
  this->hDphi->Add(rhs.hDphi);

  return (*this);

}

const MT2EstimateQCD& MT2EstimateQCD::operator-=( const MT2EstimateQCD& rhs ) {

  this->yield->Add(rhs.yield, -1.);

  this->lDphi->Add(rhs.lDphi, -1.);
  this->hDphi->Add(rhs.hDphi, -1.);

  return (*this);

}

MT2EstimateQCD MT2EstimateQCD::operator*( float k ) const{

  MT2EstimateQCD result(*this);
  result.yield->Scale(k);
  result.lDphi->Scale(k);
  result.hDphi->Scale(k);

  return result;

}

const MT2EstimateQCD& MT2EstimateQCD::operator*=( float k ) {

  this->yield->Scale(k);
  this->lDphi->Scale(k);
  this->hDphi->Scale(k);

  return (*this);

}

MT2EstimateQCD MT2EstimateQCD::operator/( float k ) const{


  MT2EstimateQCD result(*this);
  result.yield->Scale(1./k);
  result.lDphi->Scale(1./k);
  result.hDphi->Scale(1./k);

  return result;

}

const MT2EstimateQCD& MT2EstimateQCD::operator/=( float k ) {

  this->yield->Scale(1./k);
  this->lDphi->Scale(1./k);
  this->hDphi->Scale(1./k);

  return (*this);

}

  
// friend functions:

MT2EstimateQCD operator*( float k, const MT2EstimateQCD& rhs ) {

  return rhs*k;

}

MT2EstimateQCD operator/( float k, const MT2EstimateQCD& rhs ) {

  return rhs/k;

}


