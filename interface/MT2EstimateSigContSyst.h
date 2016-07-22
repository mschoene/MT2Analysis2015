#ifndef MT2EstimateSigContSyst_h
#define MT2EstimateSigContSyst_h



#include "MT2Estimate.h"
#include "TH3D.h"




class MT2EstimateSigContSyst : public MT2Estimate {

 public:

  MT2EstimateSigContSyst( const MT2EstimateSigContSyst& rhs );
  MT2EstimateSigContSyst( const std::string& aname, const MT2Region& aregion, const std::string& asystName="" );

  virtual ~MT2EstimateSigContSyst();

  virtual void setName( const std::string& newName );
  virtual void setSystName( const std::string& newSystName );

  std::string systName;

  TH3D* yield3d_systUp;
  TH3D* yield3d_systDown;

  TH3D* yield3d_crsl;
  TH1D* yield_alpha;

  TH3D* yield3d_genmet;
  TH3D* yield3d_crsl_genmet;


  const MT2EstimateSigContSyst& operator=( const MT2EstimateSigContSyst& rhs );
  const MT2EstimateSigContSyst& operator=( const MT2Estimate& rhs );

  MT2EstimateSigContSyst operator/ ( float k ) const;
  MT2EstimateSigContSyst operator* ( float k ) const;
  
  const MT2EstimateSigContSyst& operator/=( float k );
  const MT2EstimateSigContSyst& operator*=( float k );

  friend MT2EstimateSigContSyst operator/( float k, const MT2EstimateSigContSyst& rhs );
  friend MT2EstimateSigContSyst operator*( float k, const MT2EstimateSigContSyst& rhs );


  virtual void finalize();

  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

  virtual void print( std::ofstream& ofs_file, Float_t m1, Float_t m2, Int_t mt2_bin, float k, bool doGenAverage=true );
 
 private:

};





#endif
