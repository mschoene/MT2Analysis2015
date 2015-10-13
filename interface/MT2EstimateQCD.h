#ifndef MT2EstimateQCD_h
#define MT2EstimateQCD_h

#include "MT2Estimate.h"
#include "TF1.h"

#include <iostream>
#include <vector>


class MT2EstimateTree;



class MT2EstimateQCD : public MT2Estimate {

 public:

  MT2EstimateQCD( const MT2EstimateQCD& rhs );
  MT2EstimateQCD( const MT2EstimateTree& rhs, const std::string& selection="" );
  MT2EstimateQCD( const std::string& aname, const MT2Region& aregion );
  virtual ~MT2EstimateQCD();

  void projectFromTree( const MT2EstimateTree* treeEst, const std::string& selection );

  static MT2Analysis<MT2EstimateQCD>* makeAnalysisFromEstimateTree( const std::string& aname, MT2Analysis<MT2EstimateTree>* analysis, const std::string& selection="" );
  static MT2Analysis<MT2EstimateQCD>* makeAnalysisFromEstimateTreeInclusive( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* analysis, const std::string& selection="" );


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

  float getFitXmin () { return fitXmin; }
  float getFitXmax () { return fitXmax; }
  float getDphiLow() { return dphi_low; }
  virtual void setFitXmin(float xmin){ fitXmin = xmin;  }
  virtual void setFitXmax(float xmax){ fitXmax = xmax;  }
  virtual void setDphiLow(float dphi){ dphi_low = dphi; }

  virtual void finalize();

  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

  virtual void print(const std::string& ofs);

  virtual void randomizePoisson( float scale=1., int seed=13 );
  virtual void sqrtErrors( float scale=1. );

 private:

  float fitXmin, fitXmax;  // fit window
  float dphi_low;          // threshold to define low dphi region

};





#endif
