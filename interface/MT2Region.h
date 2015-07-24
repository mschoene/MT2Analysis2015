#ifndef MT2Region_h
#define MT2Region_h

#include <string>
#include <vector>




class MT2HTRegion {

 public:

  MT2HTRegion( const std::string& name );
  MT2HTRegion( const MT2HTRegion& rhs );
  MT2HTRegion( float ahtMin, float ahtMax );
  ~MT2HTRegion() {};


  // univocal identifier:
  std::string getName() const;

  std::string getNiceName() const;
  std::string getNiceNameLatex() const;


  float htMin;
  float htMax;
  float metMin() const;

  bool isInclusiveHT() const;
  bool metMinInclusiveHT( float ht ) const;

  bool operator==( const MT2HTRegion& rhs ) const;
  bool operator!=( const MT2HTRegion& rhs ) const;
  bool operator<( const MT2HTRegion& rhs ) const;

  bool isIncluded( const MT2HTRegion* htRegion ) const;

 private:

};





class MT2SignalRegion {

 public:

  MT2SignalRegion( const std::string& name );
  MT2SignalRegion( int njmin, int njmax, int nbmin, int nbmax, const std::string& mtCut="" );
  MT2SignalRegion( const MT2SignalRegion& rhs );

  ~MT2SignalRegion() {};

  // univocal identifier:
  std::string getName() const;

  std::string getNameMt() const;
  std::string getNiceName() const;
  std::string getNiceNameLatex() const;
  
  int nJetsMin; 
  int nJetsMax;
  int nBJetsMin;
  int nBJetsMax;

  std::string mtCut; // can be "", "loMT", or "hiMT"

  bool operator==( const MT2SignalRegion& rhs ) const;
  bool operator!=( const MT2SignalRegion& rhs ) const;
  bool operator<( const MT2SignalRegion& rhs ) const;
  bool operator>( const MT2SignalRegion& rhs ) const;
  bool operator<=( const MT2SignalRegion& rhs ) const;
  bool operator>=( const MT2SignalRegion& rhs ) const;

  bool isIncluded( const MT2SignalRegion* sigRegion ) const;
  

 private:

  std::string getNiceJetName( const std::string& pedix, int nmin, int nmax ) const;
  std::string getNiceJetNameLatex( const std::string& pedix, int nmin, int nmax ) const;
  std::string getSingleJetString( const std::string& prefix, int n_min , int n_max=-1 ) const;

};





class MT2Region {

 public:

  MT2Region( const MT2Region& region ) {
    htRegion_ = new MT2HTRegion(*(region.htRegion()));
    sigRegion_ = new MT2SignalRegion(*(region.sigRegion()));
  }

  MT2Region( const MT2HTRegion& htRegion, const MT2SignalRegion& sigRegion ) {
    htRegion_ = new MT2HTRegion(htRegion);
    sigRegion_ = new MT2SignalRegion(sigRegion);
  }


  MT2Region( const std::string& regionName );

  MT2Region( float htMin, float htMax=-1, int njmin=2, int njmax=-1, int nbmin=-1, int nbmax=-1, const std::string& mtCut="" ) {
    htRegion_ = new MT2HTRegion( htMin, htMax );
    sigRegion_ = new MT2SignalRegion( njmin, njmax, nbmin, nbmax, mtCut );
  }

  ~MT2Region() {};


  // this is a univocal identifier (same jet/HT thresholds always correspond to same name)
  std::string getName() const {
    return htRegion_->getName() + "_" + sigRegion_->getName();
  }


  std::vector< std::string > getNiceNames() const;
  std::vector< std::string > getNiceNamesLatex() const;

  void getBins      ( int& nBins, double*& bins ) const;
  void getBins_qcdCR( int& nBins, double*& bins ) const;

  std::string getBinName( double& min, double& max ) const;
  std::vector< std::string > getBinNames() const;

  std::string getBinNameLatex( double& min, double& max ) const;
  std::vector< std::string > getBinNamesLatex() const;

  MT2HTRegion* htRegion() const {
    return htRegion_;
  }

  MT2SignalRegion* sigRegion() const {
    return sigRegion_;
  }

  float htMin()   const { return htRegion_->htMin; };
  float htMax()   const { return htRegion_->htMax; };
  float metMin()  const { return htRegion_->metMin(); };
  
  bool isInclusiveHT()     const { return htRegion_->isInclusiveHT(); };
  bool metMinInclusiveHT( float ht ) const { return htRegion_->metMinInclusiveHT(ht); };
 
  int nJetsMin()  const { return sigRegion_->nJetsMin; };
  int nJetsMax()  const { return sigRegion_->nJetsMax; };
  int nBJetsMin() const { return sigRegion_->nBJetsMin; };
  int nBJetsMax() const { return sigRegion_->nBJetsMax; };
  std::string mtCut() const {return sigRegion_->mtCut; };


  bool operator==( const MT2Region& rhs ) const;
  bool operator!=( const MT2Region& rhs ) const;
  bool operator<( const MT2Region& rhs ) const;
  bool operator>( const MT2Region& rhs ) const;
  bool operator>=( const MT2Region& rhs ) const;
  bool operator<=( const MT2Region& rhs ) const;

  bool isIncluded( const MT2Region* region ) const;
  


 private:
  
  MT2HTRegion* htRegion_;
  MT2SignalRegion* sigRegion_;

};




#endif
