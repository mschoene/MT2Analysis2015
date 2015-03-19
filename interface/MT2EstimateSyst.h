#ifndef MT2EstimateSyst_h
#define MT2EstimateSyst_h


#include <iostream>

#include "MT2Estimate.h"
#include "MT2Analysis.h"

#include "TGraphAsymmErrors.h"



class MT2EstimateSyst : public MT2Estimate {

 public:

  MT2EstimateSyst( const MT2Estimate& rhs );
  MT2EstimateSyst( const MT2EstimateSyst& rhs );
  MT2EstimateSyst( const std::string& aname, const MT2Region& aregion );
  MT2EstimateSyst( const std::string& aname, const MT2Region& aregion, const MT2Estimate& pass, const MT2Estimate& tot );
  virtual ~MT2EstimateSyst();

  virtual void setName( const std::string& newName );

  static MT2Analysis<MT2EstimateSyst>* makeEfficiencyAnalysis( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2Estimate>* pass, MT2Analysis<MT2Estimate>* all );
  static MT2Analysis<MT2EstimateSyst>* makeAnalysisFromEstimate( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2Estimate>* analysis );

  TGraphAsymmErrors* getGraph() const;
 
  TH1D* yield_systUp;
  TH1D* yield_systDown;

  const MT2EstimateSyst& operator=( const MT2EstimateSyst& rhs );
  MT2EstimateSyst operator+( const MT2EstimateSyst& rhs ) const;
  MT2EstimateSyst operator/( const MT2EstimateSyst& rhs ) const;
  MT2EstimateSyst operator*( const MT2EstimateSyst& rhs ) const;
  MT2EstimateSyst operator*( const MT2Estimate& rhs ) const;
  const MT2EstimateSyst& operator+=( const MT2EstimateSyst& rhs );
  const MT2EstimateSyst& operator/=( const MT2EstimateSyst& rhs );
  const MT2EstimateSyst& operator*=( const MT2EstimateSyst& rhs );
  const MT2EstimateSyst& operator*=( const MT2Estimate& rhs );

  MT2EstimateSyst operator/ ( float k ) const;
  MT2EstimateSyst operator* ( float k ) const;
  const MT2EstimateSyst& operator/=( float k );
  const MT2EstimateSyst& operator*=( float k );

  friend MT2EstimateSyst operator/( float k, const MT2EstimateSyst& rhs );
  friend MT2EstimateSyst operator*( float k, const MT2EstimateSyst& rhs );


  virtual void finalize() {
    return this->addOverflow();
  }

  virtual void addOverflow();

  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

  virtual void print(const std::string& ofs);

 private:

};





#endif
