#ifndef MT2EstimateQCD_h
#define MT2EstimateQCD_h

#include "MT2Estimate.h"
#include "TF1.h"

#include <iostream>
#include <vector>


class MT2EstimateQCD : public MT2Estimate {

 public:

  MT2EstimateQCD( const MT2EstimateQCD& rhs );
  MT2EstimateQCD( const std::string& aname, const MT2Region& aregion );
  virtual ~MT2EstimateQCD();

  virtual void setName( const std::string& newName );
 
  
  TH1D* lDphi;
  TH1D* hDphi;
  TH1D* ratio;

  TF1* exp;
  TF1* expPlusC;
  TF1* expOrC;

  const MT2EstimateQCD& operator=( const MT2EstimateQCD& rhs );
  MT2EstimateQCD operator+( const MT2EstimateQCD& rhs ) const;
  MT2EstimateQCD operator-( const MT2EstimateQCD& rhs ) const;
  //MT2EstimateQCD operator/( const MT2EstimateQCD& rhs ) const;
  const MT2EstimateQCD& operator+=( const MT2EstimateQCD& rhs );
  const MT2EstimateQCD& operator-=( const MT2EstimateQCD& rhs );
  //MT2EstimateQCD operator/=( const MT2EstimateQCD& rhs ) const;

  MT2EstimateQCD operator* ( float k ) const;
  MT2EstimateQCD operator/ ( float k ) const;
  const MT2EstimateQCD& operator*=( float k );
  const MT2EstimateQCD& operator/=( float k );

  friend MT2EstimateQCD operator*( float k, const MT2EstimateQCD& rhs );
  friend MT2EstimateQCD operator/( float k, const MT2EstimateQCD& rhs );

  void fillDphi( float dphi, float mt2=0., float weight=1. );

  void doFit();

  void getRatio() {
    ratio->Divide(hDphi, lDphi);
  }

  void setFitXmin(float xmin){
    fitXmin = xmin;
  }
  void setFitXmax(float xmax){
    fitXmax = xmax;
  }
  void setDphiLow(float dphi){
    dphi_low = dphi;
  }

  virtual void finalize();

  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

  virtual void print(const std::string& ofs);

 private:

  float fitXmin, fitXmax;  // fit window
  float dphi_low;          // threshold to define low dphi region

};





#endif
