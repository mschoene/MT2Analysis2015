#ifndef MT2EstimateSig_h
#define MT2EstimateSig_h


#include "MT2Region.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TFile.h"




// this is the basic Estimate class: 
// it refers to a region, and it has a yield histogram
// other more complex classes (like LostLepton estimates)
// should inherit from this one, 
// and add further specialized data members



class MT2EstimateSig {

 public:

  MT2EstimateSig( const MT2EstimateSig& rhs );  
  MT2EstimateSig( const std::string& aname, const MT2Region& aregion );
  virtual ~MT2EstimateSig();
 
  // the region it refers to
  MT2Region* region;

  // the main data member: the yield histogram
  // classes that inherit from this one will add other data members
  TH3D* yield3d;
  TH1D* yield;

  std::string getName() const { return name; };
  virtual void setName( const std::string& newName );

  std::string getHistoName( const std::string& prefix ) const;

  void getBins( int& nBins, double* bins ) const {
    return region->getBins(nBins, bins);
  }

  MT2HTRegion* htRegion() const {
    return region->htRegion();
  }
  MT2SignalRegion* sigRegion() const {
    return region->sigRegion();
  }

  // this is a univocal identifier of the region
  // regions with the same definition (jet numbers and HT cuts)
  // have the same name
  std::string regionName() const {
    return region->getName();
  }

  const MT2EstimateSig& operator=( const MT2EstimateSig& rhs );
  MT2EstimateSig operator+( const MT2EstimateSig& rhs ) const;
  MT2EstimateSig operator-( const MT2EstimateSig& rhs ) const;
  MT2EstimateSig operator/( const MT2EstimateSig& rhs ) const;
  MT2EstimateSig operator*( const MT2EstimateSig& rhs ) const;
  const MT2EstimateSig& operator+=( const MT2EstimateSig& rhs );
  const MT2EstimateSig& operator-=( const MT2EstimateSig& rhs );
  const MT2EstimateSig& operator/=( const MT2EstimateSig& rhs );
  const MT2EstimateSig& operator*=( const MT2EstimateSig& rhs );

  MT2EstimateSig operator/ ( float k ) const;
  MT2EstimateSig operator* ( float k ) const;
  const MT2EstimateSig& operator/=( float k );
  const MT2EstimateSig& operator*=( float k );

  friend MT2EstimateSig operator*( float k, const MT2EstimateSig& rhs );
  friend MT2EstimateSig operator/( float k, const MT2EstimateSig& rhs );

  virtual void finalize() {
    return this->addOverflow();
  }

  virtual void addOverflow();
  void addOverflowSingleHisto( TH3D* yield3d );
  void addOverflowSingleHisto( TH1D* yield );

  virtual void write() const {
    yield3d->Write();
    yield->Write();
  }

  virtual void getShit( TFile* file, const std::string& path );
  
  virtual void print(const std::string& ofs);
  
  virtual void randomizePoisson( float scale=1. );

  virtual void fillYield( float mt2, int m1, int m2, float weight=1. );
  
 private:
  
  std::string name;

};





#endif
