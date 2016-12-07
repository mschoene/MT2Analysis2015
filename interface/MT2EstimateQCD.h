#ifndef MT2EstimateQCD_h
#define MT2EstimateQCD_h

#include "MT2EstimateTree.h"
#include "TF1.h"

#include <iostream>
#include <vector>





class MT2EstimateQCD : public MT2EstimateTree {

 public:

  MT2EstimateQCD( const MT2EstimateQCD& rhs );
  MT2EstimateQCD( const std::string& aname, const MT2Region& aregion );
  virtual ~MT2EstimateQCD();

  void projectFromTree( const MT2EstimateTree* treeEst, const std::string& selectionTree, const std::string& selectionDphi="" );

  static MT2Analysis<MT2EstimateQCD>* makeAnalysisFromTree( const std::string& aname, MT2Analysis<MT2EstimateTree>* analysis, const std::string& selectionTree="", const std::string& selectionDphi=""  );
  static MT2Analysis<MT2EstimateQCD>* makeAnalysisFromInclusiveTree( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* analysis, const std::string& selectionTree="", const std::string& selectionDphi=""  );

  TH1D* getRatio() const;
  TF1* getFit( const std::string& functionName, float xMin_fit, float xMax_fit, float par0=0., float par1=0. );

  void fillDphi( float dphi, float mt2=0., float weight=1. );

  virtual void setName( const std::string& newName );
 
  
  TH1D* lDphi;
  TH1D* hDphi;

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


  float getDphiLow() { return dphi_low; }
  virtual void setDphiLow(float dphi){ dphi_low = dphi; }

  virtual void finalize();

  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

  virtual void print(const std::string& ofs);

  virtual void randomizePoisson( float scale=1., int seed=13 );
  virtual void sqrtErrors( float scale=1. );

 private:

  float dphi_low;          // threshold to define low dphi region

};





#endif
