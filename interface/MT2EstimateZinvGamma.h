#ifndef MT2EstimateZinvGamma_h
#define MT2EstimateZinvGamma_h

#include "MT2Estimate.h"

#include <iostream>
#include <vector>
#include "RooRealVar.h"
#include "RooDataSet.h"


class MT2EstimateZinvGamma : public MT2Estimate {

 public:

  MT2EstimateZinvGamma( const MT2EstimateZinvGamma& rhs );
  MT2EstimateZinvGamma( const std::string& aname, const MT2Region& aregion );
  virtual ~MT2EstimateZinvGamma();

  virtual void setName( const std::string& newName );
 
  
  TH1D* sietaieta;

  // integrated over mt2:
  TH1D* iso;
  RooDataSet* isoData() const;

  // this will be used to fill RooDataSets:
  RooRealVar* x_; // iso var
  RooRealVar* w_; // weight var

  // for each bin of mt2:
  std::vector<RooDataSet*> iso_bins;
  std::vector<TH1D*> iso_bins_hist;


  void fakeDatasetsFromHistos( int seed=0 );

  const MT2EstimateZinvGamma& operator=( const MT2EstimateZinvGamma& rhs );
  MT2EstimateZinvGamma operator+( const MT2EstimateZinvGamma& rhs ) const;
  MT2EstimateZinvGamma operator-( const MT2EstimateZinvGamma& rhs ) const;
  //MT2EstimateZinvGamma operator/( const MT2EstimateZinvGamma& rhs ) const;
  const MT2EstimateZinvGamma& operator+=( const MT2EstimateZinvGamma& rhs );
  const MT2EstimateZinvGamma& operator-=( const MT2EstimateZinvGamma& rhs );
  //MT2EstimateZinvGamma operator/=( const MT2EstimateZinvGamma& rhs ) const;

  MT2EstimateZinvGamma operator* ( float k ) const;
  const MT2EstimateZinvGamma& operator*=( float k );


  void fillIso( float iso, float weight=1., float mt2=-1 );

  virtual void finalize();

  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

  virtual void print(const std::string& ofs);

 private:

};





#endif
