#ifndef MT2Estimate_h
#define MT2Estimate_h

#include "MT2Analysis.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TFile.h"


class MT2Region;


// this is the basic Estimate class: 
// it refers to a region, and it has a yield histogram
// other more complex classes (like LostLepton estimates)
// should inherit from this one, 
// and add further specialized data members



class MT2Estimate {

 public:

  MT2Estimate( const MT2Estimate& rhs );  
  MT2Estimate( const std::string& aname, const MT2Region& aregion );
  virtual ~MT2Estimate();
 
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

  void getYieldBins( int& nBins, double*& bins ) const;

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

  static MT2Analysis<MT2Estimate>* makeIntegralAnalysisFromEstimate( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2Estimate>* analysis );

  static void rebinYields( MT2Analysis<MT2Estimate>* analysis, int nBins, float xMin, float xMax );

  const MT2Estimate& operator=( const MT2Estimate& rhs );

  MT2Estimate operator+( const MT2Estimate& rhs ) const;
  MT2Estimate operator-( const MT2Estimate& rhs ) const;
  MT2Estimate operator/( const MT2Estimate& rhs ) const;
  MT2Estimate operator*( const MT2Estimate& rhs ) const;
  const MT2Estimate& operator+=( const MT2Estimate& rhs );
  const MT2Estimate& operator-=( const MT2Estimate& rhs );
  const MT2Estimate& operator/=( const MT2Estimate& rhs );
  const MT2Estimate& operator*=( const MT2Estimate& rhs );

  MT2Estimate operator/ ( float k ) const;
  MT2Estimate operator* ( float k ) const;
  const MT2Estimate& operator/=( float k );
  const MT2Estimate& operator*=( float k );

  const MT2Estimate& getMassPoint( const MT2Estimate& rhs, int mParent, int mLSP ); 
  
  friend MT2Estimate operator*( float k, const MT2Estimate& rhs );
  friend MT2Estimate operator/( float k, const MT2Estimate& rhs );

  virtual void finalize() {};

  virtual void write() const {
    yield3d->Write();
    yield->Write();
  }

  virtual void getShit( TFile* file, const std::string& path );
  
  virtual void print( const std::string& ofs );
  virtual void print( std::ofstream& ofs_file );
  virtual void print( std::ofstream& ofs_file, Int_t mt2_bin );
  
  virtual void randomizePoisson( float scale=1., int seed=13 );
  
 private:
  
  std::string name;

};



#endif
