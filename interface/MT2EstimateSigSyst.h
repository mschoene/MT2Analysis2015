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

  TH3D* yield3d_genmet;

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

  virtual void print( std::ofstream& ofs_file, Float_t m1, Float_t m2, Int_t mt2_bin, float k );
  
  virtual void write() const;

 
 private:

};





#endif
