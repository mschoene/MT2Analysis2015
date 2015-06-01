#ifndef MT2EstimateSigSyst_h
#define MT2EstimateSigSyst_h


#include <iostream>

#include "MT2EstimateSig.h"
#include "MT2Analysis.h"

#include "TGraphAsymmErrors.h"



class MT2EstimateSigSyst : public MT2EstimateSig {

 public:

  MT2EstimateSigSyst( const MT2EstimateSig& rhs );
  MT2EstimateSigSyst( const MT2EstimateSigSyst& rhs );
  MT2EstimateSigSyst( const std::string& aname, const MT2Region& aregion );
  MT2EstimateSigSyst( const std::string& aname, const MT2Region& aregion, const MT2EstimateSig& pass, const MT2EstimateSig& tot );
  virtual ~MT2EstimateSigSyst();

  virtual void setName( const std::string& newName );

  static MT2Analysis<MT2EstimateSigSyst>* makeEfficiencyAnalysis( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateSig>* pass, MT2Analysis<MT2EstimateSig>* all );
  static MT2Analysis<MT2EstimateSigSyst>* makeAnalysisFromEstimate( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateSig>* analysis );

  TGraphAsymmErrors* getGraph() const;
 
  TH1D* yield_systUp;
  TH1D* yield_systDown;

  const MT2EstimateSigSyst& operator=( const MT2EstimateSigSyst& rhs );
  const MT2EstimateSigSyst& operator=( const MT2EstimateSig& rhs );

  MT2EstimateSigSyst operator+( const MT2EstimateSigSyst& rhs ) const;
  MT2EstimateSigSyst operator-( const MT2EstimateSigSyst& rhs ) const;
  MT2EstimateSigSyst operator*( const MT2EstimateSigSyst& rhs ) const;
  MT2EstimateSigSyst operator/( const MT2EstimateSigSyst& rhs ) const;

  MT2EstimateSigSyst operator+( const MT2EstimateSig& rhs ) const;
  MT2EstimateSigSyst operator-( const MT2EstimateSig& rhs ) const;
  MT2EstimateSigSyst operator*( const MT2EstimateSig& rhs ) const;
  MT2EstimateSigSyst operator/( const MT2EstimateSig& rhs ) const;

  const MT2EstimateSigSyst& operator+=( const MT2EstimateSigSyst& rhs );
  const MT2EstimateSigSyst& operator-=( const MT2EstimateSigSyst& rhs );
  const MT2EstimateSigSyst& operator*=( const MT2EstimateSigSyst& rhs );
  const MT2EstimateSigSyst& operator/=( const MT2EstimateSigSyst& rhs );

  const MT2EstimateSigSyst& operator+=( const MT2EstimateSig& rhs );
  const MT2EstimateSigSyst& operator-=( const MT2EstimateSig& rhs );
  const MT2EstimateSigSyst& operator*=( const MT2EstimateSig& rhs );
  const MT2EstimateSigSyst& operator/=( const MT2EstimateSig& rhs );

  MT2EstimateSigSyst operator/ ( float k ) const;
  MT2EstimateSigSyst operator* ( float k ) const;
  const MT2EstimateSigSyst& operator/=( float k );
  const MT2EstimateSigSyst& operator*=( float k );

  friend MT2EstimateSigSyst operator/( float k, const MT2EstimateSigSyst& rhs );
  friend MT2EstimateSigSyst operator*( float k, const MT2EstimateSigSyst& rhs );


  virtual void finalize();

  virtual void addOverflow();

  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

  virtual void print(const std::string& ofs);

 private:

};





#endif
