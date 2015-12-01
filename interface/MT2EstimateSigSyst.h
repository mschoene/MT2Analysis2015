#ifndef MT2EstimateSigSyst_h
#define MT2EstimateSigSyst_h



#include "MT2Estimate.h"




class MT2EstimateSigSyst : public MT2Estimate {

 public:

  MT2EstimateSigSyst( const MT2EstimateSigSyst& rhs );
  MT2EstimateSigSyst( const std::string& aname, const std::string& asystName, const MT2Region& aregion );

  virtual ~MT2EstimateSigSyst();

  virtual void setName( const std::string& newName );


  std::string systName;

  TH3D* yield3d_systUp;
  TH3D* yield3d_systDown;


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
