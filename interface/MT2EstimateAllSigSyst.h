#ifndef MT2EstimateAllSigSyst_h
#define MT2EstimateAllSigSyst_h



#include "MT2Estimate.h"
#include "TH3D.h"




class MT2EstimateAllSigSyst : public MT2Estimate {

 public:

  MT2EstimateAllSigSyst( const MT2EstimateAllSigSyst& rhs );
  MT2EstimateAllSigSyst( const std::string& aname, const MT2Region& aregion, const std::string& asystName="" );

  virtual ~MT2EstimateAllSigSyst();

  virtual void setName( const std::string& newName );
  virtual void setSystName( const std::string& newSystName );

  std::string systName;

  TH3D* yield3d_systUp;
  TH3D* yield3d_systDown;


  //  TH3D* yield3d;  
  TH3D* yield3d_genmet;

  TH3D* yield3d_btag_light_UP;
  TH3D* yield3d_btag_light_DN;
  TH3D* yield3d_btag_heavy_UP;
  TH3D* yield3d_btag_heavy_DN;

  TH3D* yield3d_isr_UP;
  TH3D* yield3d_isr_DN;
  
  TH3D* yield3d_lepsf_UP;
  TH3D* yield3d_lepsf_DN;
  
  //For signal contamination
  TH3D* yield3d_crsl;
  TH3D* yield3d_alpha;


  const MT2EstimateAllSigSyst& operator=( const MT2EstimateAllSigSyst& rhs );
  const MT2EstimateAllSigSyst& operator=( const MT2Estimate& rhs );

  MT2EstimateAllSigSyst operator/ ( float k ) const;
  MT2EstimateAllSigSyst operator* ( float k ) const;
  const MT2EstimateAllSigSyst& operator/=( float k );
  const MT2EstimateAllSigSyst& operator*=( float k );

  const MT2EstimateAllSigSyst& operator+=( const MT2EstimateAllSigSyst& rhs );

  friend MT2EstimateAllSigSyst operator/( float k, const MT2EstimateAllSigSyst& rhs );
  friend MT2EstimateAllSigSyst operator*( float k, const MT2EstimateAllSigSyst& rhs );


  virtual void finalize();

  virtual void getShit( TFile* file, const std::string& path );

  //  virtual void print( std::ofstream& ofs_file, Float_t m1, Float_t m2, Int_t mt2_bin, float k, bool doGenAverage=true );
  
  virtual void write() const;

 
 private:

};





#endif
