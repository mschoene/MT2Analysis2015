#ifndef MT2EstimateSigSyst_h
#define MT2EstimateSigSyst_h



#include "MT2Estimate.h"
#include "TH3D.h"




class MT2EstimateSigSyst : public MT2Estimate {

 public:

  MT2EstimateSigSyst( const MT2EstimateSigSyst& rhs );
  MT2EstimateSigSyst( const std::string& aname, const MT2Region& aregion, const std::string& asystName="" );

  virtual ~MT2EstimateSigSyst();

  virtual void setName( const std::string& newName );
  virtual void setSystName( const std::string& newSystName );

  std::string systName;

  TH3D* yield3d_systUp;
  TH3D* yield3d_systDown;

  const MT2EstimateSigSyst& operator=( const MT2EstimateSigSyst& rhs );
  const MT2EstimateSigSyst& operator=( const MT2Estimate& rhs );

  MT2EstimateSigSyst operator/ ( float k ) const;
  MT2EstimateSigSyst operator* ( float k ) const;
  const MT2EstimateSigSyst& operator/=( float k );
  const MT2EstimateSigSyst& operator*=( float k );

  friend MT2EstimateSigSyst operator/( float k, const MT2EstimateSigSyst& rhs );
  friend MT2EstimateSigSyst operator*( float k, const MT2EstimateSigSyst& rhs );


  virtual void finalize();

  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

 
 private:

};





#endif
