#ifndef MT2Analysis_h
#define MT2Analysis_h


#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "MT2Region.h"
#include "TFile.h"

#include "../interface/mt2.h"


  
template<class T> 
class MT2Analysis {

 public:

  MT2Analysis( const std::string& aname, const std::string& regionsSet="zurich", int id=-1, const std::string& afullName="" );
  MT2Analysis( const std::string& aname, std::set<MT2HTRegion> htRegions, std::set<MT2SignalRegion> signalRegions, int id=-1, const std::string& afullName="" );
  MT2Analysis( const std::string& aname, std::set<MT2Region> regions, int id=-1, const std::string& afullName="" );
  MT2Analysis( const std::string& aname, std::set<T*> data, int id=-1, const std::string& afullName="" );
  MT2Analysis( const MT2Analysis& rhs );
  ~MT2Analysis();


  std::set<MT2Region> multiplyHTandSignal( std::set<MT2HTRegion> htRegions, std::set<MT2SignalRegion> signalRegions );
  std::set<MT2Region> getRegions() const { return regions_; };
  std::set<MT2HTRegion> getHTRegions() const; //{ return htRegions_; };
  std::set<MT2SignalRegion> getSignalRegions() const; // { return signalRegions_; };

  MT2Region* getRegion( float ht, int njets, int nbjets, float mt=-1., float mt2=-1. ) const;
  MT2Region* matchRegion( MT2Region region ) const;

  T* get( const MT2Region& r ) const;
  T* get( const MT2Tree& mt2tree ) const;
  T* get( float ht, int njets, int nbjets, float mt=-1., float mt2=-1. ) const;
  T* getWithMatch( const MT2Region& r ) const;

  std::string getName() const { return name_; };
  std::string getFullName() const { return fullName_; };
  int getColor() const { return color_; };
  int getId() const { return id_; };

  void setName( const std::string& newName );
  void setFullName( const std::string& newName ) { fullName_ = newName; };
  void setColor( const int& newColor ) { color_ = newColor; };
  void setId( const int& newId ) { id_ = newId; };
  

  const MT2Analysis& operator=( const MT2Analysis& rhs);
  //template<class T2>
  //const MT2Analysis<T>& operator=( const MT2Analysis<T2>& rhs);

  template<class T2>
  MT2Analysis<T> operator+( const MT2Analysis<T2>& rhs) const;
  template<class T2>
  MT2Analysis<T> operator-( const MT2Analysis<T2>& rhs) const;
  template<class T2>
  MT2Analysis<T> operator/( const MT2Analysis<T2>& rhs) const;
  template<class T2>
  MT2Analysis<T> operator*( const MT2Analysis<T2>& rhs) const;

  template<class T2>
  const MT2Analysis<T>& operator+=( const MT2Analysis<T2>& rhs);
  template<class T2>
  const MT2Analysis<T>& operator-=( const MT2Analysis<T2>& rhs);
  template<class T2>
  const MT2Analysis<T>& operator/=( const MT2Analysis<T2>& rhs);
  template<class T2>
  const MT2Analysis<T>& operator*=( const MT2Analysis<T2>& rhs);

  void add( const MT2Analysis& rhs );
  void divide( const MT2Analysis& rhs );

  MT2Analysis operator* ( float k ) const;
  MT2Analysis operator/ ( float k ) const;
  const MT2Analysis& operator*=( float k );
  const MT2Analysis& operator/=( float k );

  static MT2Analysis* readFromFile( const std::string& fileName, const std::string& matchName="" );
  static std::vector<MT2Analysis*> readAllFromFile( const std::string& fileName, const std::string& matchName="", bool verbose=true );
  void writeToFile( const std::string& fileName, const std::string& option="UPDATE", bool overwrite=true );
  void addToFile( const std::string& fileName, bool overwrite=true ) {
    return this->writeToFile(fileName,"UPDATE",overwrite);
  }

  static void printFromFile( const std::string& fileName, const std::string& ofs, const std::string& matchName="" );
  static void print( const std::vector<MT2Analysis*> analyses, const std::string& ofs, const std::string& matchName="" );
  void print( const std::string& ofs, MT2Region* matchRegion=0 ) const;
  void print( std::ofstream& ofs_file, MT2Region* matchRegion=0 ) const;
  void print( std::ofstream& ofs_file, MT2HTRegion* thisHTRegion=0 ) const;

  void printRegions() const;

  void finalize();

  void randomizePoisson( float scale=1. );
  void sqrtErrors      ( float scale=1. );

  std::set<T*> data;



 private:

  void createAnalysisStructure() {

    for( std::set<MT2Region>::iterator iR=regions_.begin(); iR!=regions_.end(); ++iR ) {
      MT2Region thisRegion(*iR);
      T* t = new T(name_, thisRegion);
      data.insert(t);
    }

  }

  std::set<MT2Region> regions_;

  std::string name_;
  std::string fullName_;

  int id_;
  int color_;
  
  void setDefaultColor() {
    
    std::string aname = this->name_;

    if( aname == "QCD" )
      this->color_ = 401;
    else if( aname == "WJets" )
      this->color_ = 417;
    else if( aname == "ZJets" )
      this->color_ = 419;
    else if( aname == "Top" )
      this->color_ = 855;
    else if( aname == "Other" )
      this->color_ = 9;
    else if( aname == "qcdEstimate")
      this->color_ = 402;
    else if( aname == "llepEstimate")
      this->color_ = 430;
    else if( aname == "ZinvEstimate")
      this->color_ = 418;
    else
      this->color_ = 1;

  }

};





// constructors

template<class T> 
MT2Analysis<T>::MT2Analysis( const std::string& aname, const std::string& regionsSet, int aid, const std::string& afullname ) {


  name_ = aname;
  fullName_ = (afullname!="") ? afullname : name_;
  id_ = aid;

  this->setDefaultColor();

  if( regionsSet=="8TeV" ) {

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,  750. ));
    htRegions.insert(MT2HTRegion( 750., 1200. ));
    htRegions.insert(MT2HTRegion(1200.,   -1. ));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2, 2, 0, 0));  // 2j0b
    signalRegions.insert(MT2SignalRegion(2, 2, 1, 2));  // 2j1to2b
    signalRegions.insert(MT2SignalRegion(3, 5, 0, 0));  // 3to5j0b
    signalRegions.insert(MT2SignalRegion(3, 5, 1, 1));  // 3to5j1b
    signalRegions.insert(MT2SignalRegion(3, 5, 2, 2));  // 3to5j2b
    signalRegions.insert(MT2SignalRegion(6, -1, 0, 0));  // 6j0b
    signalRegions.insert(MT2SignalRegion(6, -1, 1, 1));  // 6j1b
    signalRegions.insert(MT2SignalRegion(6, -1, 2, 2));  // 6j2b
    signalRegions.insert(MT2SignalRegion(-1, -1, 3, -1));  // 3b

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_noCut" ) {

    regions_.insert(MT2Region( 0. ));

  } else if( regionsSet=="13TeV_inclusive" ) {

    regions_.insert(MT2Region( 200. )); // inclusive 200-inf at least one jet requirement

  } else if( regionsSet=="13TeV_inclusive450" ) {

    regions_.insert(MT2Region( 450., -1., 2., -1. )); // inclusive 450-inf at least two jet requirement

  } else if( regionsSet=="13TeV_inclusive_bjets" ) {

    regions_.insert(MT2Region( 450., -1., 2, -1, 0, 0 )); 
    regions_.insert(MT2Region( 450., -1., 2, -1, 1, 1 )); 
    regions_.insert(MT2Region( 450., -1., 2, -1, 2, 2 )); 
    regions_.insert(MT2Region( 450., -1., 2, -1, 3,-1 )); 

  } else if( regionsSet=="13TeV" ) {

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  2)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 1,  2)); 
    signalRegions.insert(MT2SignalRegion(2, -1, 3, -1)); 

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_CSA14" ) {

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 1,  1)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "loMT" ));
    signalRegions.insert(MT2SignalRegion(4, -1, 2,  2, "loMT" ));
    signalRegions.insert(MT2SignalRegion(2, -1, 3, -1, "loMT" ));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "hiMT" ));
    signalRegions.insert(MT2SignalRegion(4, -1, 2,  2, "hiMT" ));
    signalRegions.insert(MT2SignalRegion(2, -1, 3, -1, "hiMT" ));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_CSA14_noMT" ) {

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  0)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1)); 
    signalRegions.insert(MT2SignalRegion(4, -1, 1,  1)); 
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, -1, 2,  2));
    signalRegions.insert(MT2SignalRegion(3, -1, 3, -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_PHYS14" ){
    
    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, 8, 0,  0));
    signalRegions.insert(MT2SignalRegion(9, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, 8, 1,  1));
    signalRegions.insert(MT2SignalRegion(9, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "hiMT"));
    
    regions_ = multiplyHTandSignal( htRegions, signalRegions );
    
  } else if( regionsSet=="13TeV_PHYS14_noMT" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,    -1.));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0 ));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0 ));
    signalRegions.insert(MT2SignalRegion(7, 8, 0,  0 ));
    signalRegions.insert(MT2SignalRegion(9, -1, 0,  0 ));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1 ));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1 ));
    signalRegions.insert(MT2SignalRegion(7, 8, 1,  1 ));
    signalRegions.insert(MT2SignalRegion(9, -1, 1,  1 ));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2 ));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2 ));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2 ));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2 ));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1 ));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1 ));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

  } else if( regionsSet=="13TeV_PHYS14_hiHT" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, 8, 0,  0));
    signalRegions.insert(MT2SignalRegion(9, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, 8, 1,  1));
    signalRegions.insert(MT2SignalRegion(9, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, 8, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(9, -1, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "hiMT"));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="zurich" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2));
    signalRegions.insert(MT2SignalRegion(7, -1, 2,  2));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

  } else if( regionsSet=="13TeV_inclusive2j" ){

    regions_.insert(MT2Region( 200., -1., 2., -1. )); // inclusive 200-inf at least two jet requirement

   } else if( regionsSet=="zurichPlus_noMonojet" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 200.,  450.));
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2));
    signalRegions.insert(MT2SignalRegion(7, -1, 2,  2));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

  } else if( regionsSet=="zurichPlus" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 200.,  450.));
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2));
    signalRegions.insert(MT2SignalRegion(7, -1, 2,  2));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

    regions_.insert(MT2Region(200., -1., 1, 1, 0, 0));
    regions_.insert(MT2Region(200., -1., 1, 1, 1, -1));

  } else if( regionsSet=="zurichNew" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 200.,  450.));
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2));
    signalRegions.insert(MT2SignalRegion(7, -1, 2,  2));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

  } else if( regionsSet=="zurich_monojet" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2));
    signalRegions.insert(MT2SignalRegion(7, -1, 2,  2));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

    //    regions_.insert(MT2Region(450., -1., 1, 1, 0, -1)); // monojet region
    regions_.insert(MT2Region(200., -1., 1, 1, 0, -1)); // monojet region


  } else if( regionsSet=="monojet" ){ 

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 200.,  450.));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2));
    signalRegions.insert(MT2SignalRegion(7, -1, 2,  2));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1));
    
    regions_ = multiplyHTandSignal( htRegions, signalRegions );
    
    regions_.insert(MT2Region(200., -1., 1, 1, 0, 0));
    regions_.insert(MT2Region(200., -1., 1, 1, 1, -1));

  } else if( regionsSet=="zurich_llep" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  2));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

  } else if( regionsSet=="zurichPlus_llep" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 200.,  450.));
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  2));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

    regions_.insert(MT2Region(200., -1., 1, 1, 0, 0));
    regions_.insert(MT2Region(200., -1., 1, 1, 1, -1));

  } else if( regionsSet=="zurich_monojet_llep" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  2));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

    regions_.insert(MT2Region(450., -1., 1, 1, 0, -1)); // monojet region
    
  } else if( regionsSet=="zurich_onlyHT" ){

    //regions_.insert(MT2Region( 200.,    -1., 1.,  1.));
    regions_.insert(MT2Region( 200.,   450., 2., -1.));
    regions_.insert(MT2Region( 450.,   575., 2., -1.));
    regions_.insert(MT2Region( 575.,  1000., 2., -1.));
    regions_.insert(MT2Region(1000.,  1500., 2., -1.));
    regions_.insert(MT2Region(1500.,    -1., 2., -1.));


  } else if( regionsSet=="zurich_HTtriggers" ){

    regions_.insert(MT2Region( 450.,   575.)); // no cut on jets
    regions_.insert(MT2Region( 575.,  1000.));
    regions_.insert(MT2Region(1000.,    -1.));


  } else if( regionsSet=="zurich_HTtriggers2" ){

    regions_.insert(MT2Region( 450.,   575.)); // no cut on jets
    regions_.insert(MT2Region( 575.,   900.));
    regions_.insert(MT2Region( 900.,    -1.));


  } else if( regionsSet=="zurich_bMerged" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  -1));
    signalRegions.insert(MT2SignalRegion(4,  6, 0,  -1));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="zurich_onlyJets" ){

    regions_.insert(MT2Region(450., -1., 2,  3, 0,  0));
    regions_.insert(MT2Region(450., -1., 4,  6, 0,  0));
    regions_.insert(MT2Region(450., -1., 7, -1, 0,  0));
    regions_.insert(MT2Region(450., -1., 2,  3, 1,  1));
    regions_.insert(MT2Region(450., -1., 4,  6, 1,  1));
    regions_.insert(MT2Region(450., -1., 7, -1, 1,  1));
    regions_.insert(MT2Region(450., -1., 2,  3, 2,  2));
    regions_.insert(MT2Region(450., -1., 4,  6, 2,  2));
    regions_.insert(MT2Region(450., -1., 7, -1, 2,  2));
    regions_.insert(MT2Region(450., -1., 2,  6, 3, -1));
    regions_.insert(MT2Region(450., -1., 7, -1, 3, -1));


  } else if( regionsSet=="zurich_onlyJets_noB" ){

    //regions_.insert(MT2Region(200., -1., 1,  1, 0,  -1));
    regions_.insert(MT2Region(200., -1., 2,  3, 0,  -1));
    regions_.insert(MT2Region(200., -1., 4,  6, 0,  -1));
    regions_.insert(MT2Region(200., -1., 7, -1, 0,  -1));


  } else if( regionsSet=="darkMatter_max1b" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4,  6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4,  6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="darkMatter_max1b_2j" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2, 2, 0,  0));
    signalRegions.insert(MT2SignalRegion(3, 3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2, 2, 1,  1));
    signalRegions.insert(MT2SignalRegion(3, 3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="darkMatter_max1b_2j_4j" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2, 2, 0,  0));
    signalRegions.insert(MT2SignalRegion(3, 3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2, 2, 1,  1));
    signalRegions.insert(MT2SignalRegion(3, 3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, -1, 1,  1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="darkMatter_max1b_all_2j_4j" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2, 2, 0,  1));
    signalRegions.insert(MT2SignalRegion(3, 3, 0,  1));
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="darkMatter_allb" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2, 3, 0,  -1));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  -1));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="darkMatter_allb_2j" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  2, 0,  -1));
    signalRegions.insert(MT2SignalRegion(3,  3, 0,  -1));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  -1));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="darkMatter_allb_2j_4j" || regionsSet=="DMv0" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  2, 0,  -1));
    signalRegions.insert(MT2SignalRegion(3,  3, 0,  -1));
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_PHYS14_hiJet_mergeHT" ){

    regions_.insert(MT2Region(450., 575., 2, 3, 0,  0));
    regions_.insert(MT2Region(575., 1000., 2, 3, 0,  0));
    regions_.insert(MT2Region(1000., 1500., 2, 3, 0,  0));
    regions_.insert(MT2Region(1500., -1, 2, 3, 0,  0));


    regions_.insert(MT2Region(450., 575., 4, 6, 0,  0));
    regions_.insert(MT2Region(575., 1000., 4, 6, 0,  0));
    regions_.insert(MT2Region(1000., 1500., 4, 6, 0,  0));
    regions_.insert(MT2Region(1500., -1, 4, 6, 0,  0));

    regions_.insert(MT2Region(450., 575., 7, 8, 0,  0));
    regions_.insert(MT2Region(575., 1000., 7, 8, 0,  0));
    regions_.insert(MT2Region(1000., 1500., 7, 8, 0,  0));
    regions_.insert(MT2Region(1500., -1, 7, 8, 0,  0));

    regions_.insert(MT2Region(450., -1, 9, -1, 0,  0));

    regions_.insert(MT2Region(450., 575., 2, 3, 1,  1));
    regions_.insert(MT2Region(575., 1000., 2, 3, 1,  1));
    regions_.insert(MT2Region(1000., 1500., 2, 3, 1,  1));
    regions_.insert(MT2Region(1500., -1, 2, 3, 1,  1));

    regions_.insert(MT2Region(450., 575., 4, 6, 1,  1));
    regions_.insert(MT2Region(575., 1000., 4, 6, 1,  1));
    regions_.insert(MT2Region(1000., 1500., 4, 6, 1,  1));
    regions_.insert(MT2Region(1500., -1, 4, 6, 1,  1));

    regions_.insert(MT2Region(450., 575., 7, 8, 1,  1));
    regions_.insert(MT2Region(575., 1000., 7, 8, 1,  1));
    regions_.insert(MT2Region(1000., 1500., 7, 8, 1,  1));
    regions_.insert(MT2Region(1500., -1, 7, 8, 1,  1));

    regions_.insert(MT2Region(450., -1, 9, -1, 1,  1));
    
    regions_.insert(MT2Region(450., 575., 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(575., 1000., 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(1000., 1500., 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(1500., -1, 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(450., 575., 2, 3, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(575., 1000., 2, 3, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(1000., 1500., 2, 3, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(1500., -1, 2, 3, 2,  2, "hiMT" ));

    regions_.insert(MT2Region(450., 575., 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(575., 1000., 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(1000., 1500., 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(1500., -1, 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(450., 575., 4, 6, 2, 2, "hiMT" ));
    regions_.insert(MT2Region(575., 1000., 4, 6, 2, 2, "hiMT" ));
    regions_.insert(MT2Region(1000., 1500., 4, 6, 2, 2, "hiMT" ));
    regions_.insert(MT2Region(1500., -1, 4, 6, 2, 2, "hiMT" ));

    regions_.insert(MT2Region(450., 575., 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(575., 1000., 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(1000., 1500., 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(1500., -1, 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(450., 575., 7, 8, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(575., 1000., 7, 8, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(1000., 1500., 7, 8, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(1500., -1, 7, 8, 2,  2, "hiMT" ));

    regions_.insert(MT2Region(450., -1, 9, -1, 2,  2, "loMT" ));
    regions_.insert(MT2Region(450., -1, 9, -1, 2,  2, "hiMT" ));

    regions_.insert(MT2Region(450., 575., 2,  6, 3,  -1, "loMT"));
    regions_.insert(MT2Region(450., 575., 2,  6, 3,  -1, "hiMT"));
    regions_.insert(MT2Region(450., 575., 7, -1, 3,  -1, "loMT"));
    regions_.insert(MT2Region(450., 575., 7, -1, 3,  -1, "hiMT"));

    regions_.insert(MT2Region(575., 1000., 2,  6, 3,  -1, "loMT"));
    regions_.insert(MT2Region(575., 1000., 2,  6, 3,  -1, "hiMT"));
    regions_.insert(MT2Region(575., 1000., 7, -1, 3,  -1, "loMT"));
    regions_.insert(MT2Region(575., 1000., 7, -1, 3,  -1, "hiMT"));

    regions_.insert(MT2Region(1000., -1, 2,  6, 3,  -1, "loMT"));
    regions_.insert(MT2Region(1000., -1, 2,  6, 3,  -1, "hiMT"));
    regions_.insert(MT2Region(1000., -1, 7, -1, 3,  -1, "loMT"));
    regions_.insert(MT2Region(1000., -1, 7, -1, 3,  -1, "hiMT"));


  } else if( regionsSet=="13TeV_PHYS14_hiJet_mergeHT_noMT" ){

    regions_.insert(MT2Region(450., 575., 2, 3, 0,  0));
    regions_.insert(MT2Region(575., 1000., 2, 3, 0,  0));
    regions_.insert(MT2Region(1000., 1500., 2, 3, 0,  0));
    regions_.insert(MT2Region(1500., -1, 2, 3, 0,  0));


    regions_.insert(MT2Region(450., 575., 4, 6, 0,  0));
    regions_.insert(MT2Region(575., 1000., 4, 6, 0,  0));
    regions_.insert(MT2Region(1000., 1500., 4, 6, 0,  0));
    regions_.insert(MT2Region(1500., -1, 4, 6, 0,  0));

    regions_.insert(MT2Region(450., 575., 7, 8, 0,  0));
    regions_.insert(MT2Region(575., 1000., 7, 8, 0,  0));
    regions_.insert(MT2Region(1000., 1500., 7, 8, 0,  0));
    regions_.insert(MT2Region(1500., -1, 7, 8, 0,  0));

    regions_.insert(MT2Region(450., -1, 9, -1, 0,  0));

    regions_.insert(MT2Region(450., 575., 2, 3, 1,  1));
    regions_.insert(MT2Region(575., 1000., 2, 3, 1,  1));
    regions_.insert(MT2Region(1000., 1500., 2, 3, 1,  1));
    regions_.insert(MT2Region(1500., -1, 2, 3, 1,  1));

    regions_.insert(MT2Region(450., 575., 4, 6, 1,  1));
    regions_.insert(MT2Region(575., 1000., 4, 6, 1,  1));
    regions_.insert(MT2Region(1000., 1500., 4, 6, 1,  1));
    regions_.insert(MT2Region(1500., -1, 4, 6, 1,  1));

    regions_.insert(MT2Region(450., 575., 7, 8, 1,  1));
    regions_.insert(MT2Region(575., 1000., 7, 8, 1,  1));
    regions_.insert(MT2Region(1000., 1500., 7, 8, 1,  1));
    regions_.insert(MT2Region(1500., -1, 7, 8, 1,  1));

    regions_.insert(MT2Region(450., -1, 9, -1, 1,  1));
    
    regions_.insert(MT2Region(450., 575., 2, 3, 2,  2));
    regions_.insert(MT2Region(575., 1000., 2, 3, 2,  2));
    regions_.insert(MT2Region(1000., 1500., 2, 3, 2,  2));
    regions_.insert(MT2Region(1500., -1, 2, 3, 2,  2));
   
    regions_.insert(MT2Region(450., 575., 4, 6, 2, 2));
    regions_.insert(MT2Region(575., 1000., 4, 6, 2, 2));
    regions_.insert(MT2Region(1000., 1500., 4, 6, 2, 2));
    regions_.insert(MT2Region(1500., -1, 4, 6, 2, 2));
   
    regions_.insert(MT2Region(450., 575., 7, 8, 2,  2));
    regions_.insert(MT2Region(575., 1000., 7, 8, 2,  2));
    regions_.insert(MT2Region(1000., 1500., 7, 8, 2,  2));
    regions_.insert(MT2Region(1500., -1, 7, 8, 2,  2));
   
    regions_.insert(MT2Region(450., -1, 9, -1, 2,  2));
   
    regions_.insert(MT2Region(450., 575., 2,  6, 3,  -1));
    regions_.insert(MT2Region(450., 575., 7, -1, 3,  -1));

    regions_.insert(MT2Region(575., 1000., 2,  6, 3,  -1));
    regions_.insert(MT2Region(575., 1000., 7, -1, 3,  -1));
    
    regions_.insert(MT2Region(1000., -1, 2,  6, 3,  -1));
    regions_.insert(MT2Region(1000., -1, 7, -1, 3,  -1));
    
  } else if( regionsSet=="13TeV_PHYS14_hiJet_extremeHT" ){

    regions_.insert(MT2Region(450., 575., 2, 3, 0,  0));
    regions_.insert(MT2Region(575., 1000., 2, 3, 0,  0));
    regions_.insert(MT2Region(1000., 2000., 2, 3, 0,  0));

    regions_.insert(MT2Region(450., 575., 4, 6, 0,  0));
    regions_.insert(MT2Region(575., 1000., 4, 6, 0,  0));
    regions_.insert(MT2Region(1000., 2000., 4, 6, 0,  0));

    regions_.insert(MT2Region(450., 575., 7, 8, 0,  0));
    regions_.insert(MT2Region(575., 1000., 7, 8, 0,  0));
    regions_.insert(MT2Region(1000., 2000., 7, 8, 0,  0));

    regions_.insert(MT2Region(450., 2000, 9, -1, 0,  0));

    regions_.insert(MT2Region(450., 575., 2, 3, 1,  1));
    regions_.insert(MT2Region(575., 1000., 2, 3, 1,  1));
    regions_.insert(MT2Region(1000., 2000., 2, 3, 1,  1));

    regions_.insert(MT2Region(450., 575., 4, 6, 1,  1));
    regions_.insert(MT2Region(575., 1000., 4, 6, 1,  1));
    regions_.insert(MT2Region(1000., 2000., 4, 6, 1,  1));

    regions_.insert(MT2Region(450., 575., 7, 8, 1,  1));
    regions_.insert(MT2Region(575., 1000., 7, 8, 1,  1));
    regions_.insert(MT2Region(1000., 2000., 7, 8, 1,  1));

    regions_.insert(MT2Region(450., 2000, 9, -1, 1,  1));
    
    regions_.insert(MT2Region(450., 575., 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(575., 1000., 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(1000., 2000., 2, 3, 2,  2, "loMT" ));
    regions_.insert(MT2Region(450., 575., 2, 3, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(575., 1000., 2, 3, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(1000., 2000., 2, 3, 2,  2, "hiMT" ));

    regions_.insert(MT2Region(450., 575., 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(575., 1000., 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(1000., 2000., 4, 6, 2, 2, "loMT" ));
    regions_.insert(MT2Region(450., 575., 4, 6, 2, 2, "hiMT" ));
    regions_.insert(MT2Region(575., 1000., 4, 6, 2, 2, "hiMT" ));
    regions_.insert(MT2Region(1000., 2000., 4, 6, 2, 2, "hiMT" ));

    regions_.insert(MT2Region(450., 575., 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(575., 1000., 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(1000., 2000., 7, 8, 2,  2, "loMT" ));
    regions_.insert(MT2Region(450., 575., 7, 8, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(575., 1000., 7, 8, 2,  2, "hiMT" ));
    regions_.insert(MT2Region(1000., 2000., 7, 8, 2,  2, "hiMT" ));

    regions_.insert(MT2Region(450., 2000, 9, -1, 2,  2, "loMT" ));
    regions_.insert(MT2Region(450., 2000, 9, -1, 2,  2, "hiMT" ));

    regions_.insert(MT2Region(450., 575., 2,  6, 3,  -1, "loMT"));
    regions_.insert(MT2Region(450., 575., 2,  6, 3,  -1, "hiMT"));
    regions_.insert(MT2Region(450., 575., 7, -1, 3,  -1, "loMT"));
    regions_.insert(MT2Region(450., 575., 7, -1, 3,  -1, "hiMT"));

    regions_.insert(MT2Region(575., 1000., 2,  6, 3,  -1,
    fitPurity( cfg, thisLoosePurity_ht, thisTightPurity_ht, thisEstimate_ht->x_, thisEstimate_ht->iso_bins, templatePrompt->iso, templateFake->iso); "loMT"));
    regions_.insert(MT2Region(575., 1000., 2,  6, 3,  -1, "hiMT"));
    regions_.insert(MT2Region(575., 1000., 7, -1, 3,  -1, "loMT"));
    regions_.insert(MT2Region(575., 1000., 7, -1, 3,  -1, "hiMT"));

    regions_.insert(MT2Region(1000., 2000., 2,  6, 3,  -1, "loMT"));
    regions_.insert(MT2Region(1000., 2000., 2,  6, 3,  -1, "hiMT"));
    regions_.insert(MT2Region(1000., 2000., 7, -1, 3,  -1, "loMT"));
    regions_.insert(MT2Region(1000., 2000., 7, -1, 3,  -1, "hiMT"));

    regions_.insert(MT2Region(2000., -1, 2,  -1, 0,  -1));

  } else if( regionsSet=="13TeV_PHYS14_hiJet_extremeHT_noMT" ){

    regions_.insert(MT2Region(450., 575., 2, 3, 0,  0));
    regions_.insert(MT2Region(575., 1000., 2, 3, 0,  0));
    regions_.insert(MT2Region(1000., 2000., 2, 3, 0,  0));

    regions_.insert(MT2Region(450., 575., 4, 6, 0,  0));
    regions_.insert(MT2Region(575., 1000., 4, 6, 0,  0));
    regions_.insert(MT2Region(1000., 2000., 4, 6, 0,  0));

    regions_.insert(MT2Region(450., 575., 7, 8, 0,  0));
    regions_.insert(MT2Region(575., 1000., 7, 8, 0,  0));
    regions_.insert(MT2Region(1000., 2000., 7, 8, 0,  0));

    regions_.insert(MT2Region(450., 2000, 9, -1, 0,  0));

    regions_.insert(MT2Region(450., 575., 2, 3, 1,  1));
    regions_.insert(MT2Region(575., 1000., 2, 3, 1,  1));
    regions_.insert(MT2Region(1000., 2000., 2, 3, 1,  1));

    regions_.insert(MT2Region(450., 575., 4, 6, 1,  1));
    regions_.insert(MT2Region(575., 1000., 4, 6, 1,  1));
    regions_.insert(MT2Region(1000., 2000., 4, 6, 1,  1));

    regions_.insert(MT2Region(450., 575., 7, 8, 1,  1));
    regions_.insert(MT2Region(575., 1000., 7, 8, 1,  1));
    regions_.insert(MT2Region(1000., 2000., 7, 8, 1,  1));

    regions_.insert(MT2Region(450., 2000, 9, -1, 1,  1));
    
    regions_.insert(MT2Region(450., 575., 2, 3, 2,  2));
    regions_.insert(MT2Region(575., 1000., 2, 3, 2,  2));
    regions_.insert(MT2Region(1000., 2000., 2, 3, 2,  2));
  
    regions_.insert(MT2Region(450., 575., 4, 6, 2, 2));
    regions_.insert(MT2Region(575., 1000., 4, 6, 2, 2));
    regions_.insert(MT2Region(1000., 2000., 4, 6, 2, 2));
  
    regions_.insert(MT2Region(450., 575., 7, 8, 2,  2));
    regions_.insert(MT2Region(575., 1000., 7, 8, 2,  2));
    regions_.insert(MT2Region(1000., 2000., 7, 8, 2,  2));
  
    regions_.insert(MT2Region(450., 2000, 9, -1, 2,  2));
  
    regions_.insert(MT2Region(450., 575., 2,  6, 3,  -1));
    regions_.insert(MT2Region(450., 575., 7, -1, 3,  -1));
  
    regions_.insert(MT2Region(575., 1000., 2,  6, 3,  -1));
    regions_.insert(MT2Region(575., 1000., 7, -1, 3,  -1));
  
    regions_.insert(MT2Region(1000., 2000., 2,  6, 3,  -1));
    regions_.insert(MT2Region(1000., 2000., 7, -1, 3,  -1));
  
    regions_.insert(MT2Region(2000., -1, 2,  -1, 0,  -1));

  } else if( regionsSet=="13TeV_PHYS14_loJet_hiHT" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));

    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  3, 0,  0));
    signalRegions.insert(MT2SignalRegion(4, 6, 0,  0));
    signalRegions.insert(MT2SignalRegion(7, -1, 0,  0));
    signalRegions.insert(MT2SignalRegion(2,  3, 1,  1));
    signalRegions.insert(MT2SignalRegion(4, 6, 1,  1));
    signalRegions.insert(MT2SignalRegion(7, -1, 1,  1));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  3, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(4, 6, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 2,  2, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 2,  2, "hiMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(2,  6, 3,  -1, "hiMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "loMT"));
    signalRegions.insert(MT2SignalRegion(7, -1, 3,  -1, "hiMT"));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );


  } else if( regionsSet=="13TeV_onlyHT" ) {

    regions_.insert(MT2Region( 450.,   575.)); // no cut on jets
    regions_.insert(MT2Region( 575.,  1000.));
    regions_.insert(MT2Region(1000.,    -1.));


  } else if( regionsSet=="13TeV_onlyJets" ) {

    regions_.insert(MT2Region(450., -1., 2,  3, 0,  0)); 
    regions_.insert(MT2Region(450., -1., 4, -1, 0,  0)); 
    regions_.insert(MT2Region(450., -1., 2,  3, 1,  1)); 
    regions_.insert(MT2Region(450., -1., 4, -1, 1,  1)); 
    regions_.insert(MT2Region(450., -1., 2,  3, 2,  2, "loMT" )); 
    regions_.insert(MT2Region(450., -1., 4, -1, 2,  2, "loMT" )); 
    regions_.insert(MT2Region(450., -1., 3, -1, 3, -1, "loMT" )); 
    regions_.insert(MT2Region(450., -1., 2,  3, 2,  2, "hiMT" )); 
    regions_.insert(MT2Region(450., -1., 4, -1, 2,  2, "hiMT" )); 
    regions_.insert(MT2Region(450., -1., 3, -1, 3, -1, "hiMT" )); 

  } else if( regionsSet=="DMv1" ){

    std::set<MT2HTRegion> htRegions;
    htRegions.insert(MT2HTRegion( 450.,   575.));
    htRegions.insert(MT2HTRegion( 575.,  1000.));
    htRegions.insert(MT2HTRegion(1000.,  1500.));
    htRegions.insert(MT2HTRegion(1500.,    -1 ));
    
    std::set<MT2SignalRegion> signalRegions;
    signalRegions.insert(MT2SignalRegion(2,  2, 0,  -1));
    signalRegions.insert(MT2SignalRegion(3,  3, 0,  -1));
    signalRegions.insert(MT2SignalRegion(4, -1, 0,  -1));

    regions_ = multiplyHTandSignal( htRegions, signalRegions );

    regions_.insert(MT2Region(450., -1., 1, 1, 0, -1)); // monojet region

  } else {

    std::cout << "[MT2Analysis::MT2Analysis] Analysis region set '" << regionsSet << "' not implemented yet. Exiting." << std::endl;
    exit(917);

  }

  this->createAnalysisStructure();

}




template<class T> 
MT2Analysis<T>::MT2Analysis( const std::string& aname, std::set<MT2Region> regions, int aid, const std::string& afullname ) {

  name_ = aname;
  fullName_ = (afullname!="") ? afullname : name_;
  id_ = aid;
  
  this->setDefaultColor();

  regions_ = regions;

  this->createAnalysisStructure();

}


template<class T> 
MT2Analysis<T>::MT2Analysis( const std::string& aname, std::set<MT2HTRegion> htRegions, std::set<MT2SignalRegion> signalRegions, int aid, const std::string& afullname ) {

  name_ = aname;
  fullName_ = (afullname!="") ? afullname : name_;
  id_ = aid;
  
  this->setDefaultColor();

  for( std::set<MT2HTRegion>::iterator iHT=htRegions.begin(); iHT!=htRegions.end(); ++iHT )  {
    for( std::set<MT2SignalRegion>::iterator iSR=signalRegions.begin(); iSR!=signalRegions.end(); ++iSR ) {
      MT2Region thisRegion( *iHT, *iSR );
      regions_.insert(thisRegion);
    }
  }


  this->createAnalysisStructure();

}


template<class T> 
MT2Analysis<T>::MT2Analysis( const std::string& aname, std::set<T*> newdata, int aid, const std::string& afullname ) {

  name_ = aname;
  fullName_ = (afullname!="") ? afullname : name_;
  id_ = aid;
  
  this->setDefaultColor();

  for( typename std::set<T*>::iterator idata=newdata.begin(); idata!=newdata.end(); ++idata ) {

    MT2Region* thisRegion = new MT2Region(*((*idata)->region));
    regions_.insert( *thisRegion );

    T* newdata = new T( *(*idata) );
    this->data.insert( newdata );

  }

  this->setName(aname);

}



template<class T> 
MT2Analysis<T>::MT2Analysis( const MT2Analysis& rhs ) {

  //regions_ = rhs.getRegions();

  name_ = rhs.name_;
  fullName_ = rhs.fullName_;
  id_ = rhs.id_;
  
  color_ = rhs.getColor();

  for( typename std::set<T*>::iterator idata=rhs.data.begin(); idata!=rhs.data.end(); ++idata ) {

    MT2Region* thisRegion = new MT2Region(*((*idata)->region));
    regions_.insert( *thisRegion );

    T* newdata = new T( *(*idata) );
    this->data.insert( newdata );

  }


  //this->createAnalysisStructure();

  //*this = rhs;

}



// destructor

template<class T> 
MT2Analysis<T>::~MT2Analysis() {

  for( typename std::set<T*>::iterator i=data.begin(); i!=data.end(); ++i ) {
    delete *i;
  }

}




// other methods


template<class T>
std::set<MT2Region> MT2Analysis<T>::multiplyHTandSignal( std::set<MT2HTRegion> htRegions, std::set<MT2SignalRegion> signalRegions ) {


  std::set<MT2Region> regions;

  for( std::set<MT2HTRegion>::iterator iHT=htRegions.begin(); iHT!=htRegions.end(); ++iHT ) {

    for( std::set<MT2SignalRegion>::iterator iSR=signalRegions.begin(); iSR!=signalRegions.end(); ++iSR ) {

      MT2Region newregion( *iHT, *iSR );
      regions.insert( newregion );

    }

  }

  return regions;

}




template<class T>
std::set<MT2HTRegion> MT2Analysis<T>::getHTRegions() const {

  std::set<MT2HTRegion> regions;

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {

    regions.insert( *((*it)->region->htRegion()) );

  }

  return regions;

}



template<class T>
std::set<MT2SignalRegion> MT2Analysis<T>::getSignalRegions() const {

  std::set<MT2SignalRegion> regions;

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {

    regions.insert( *((*it)->region->sigRegion()) );

  }

  return regions;

}








template<class T>
void MT2Analysis<T>::printRegions() const {

  std::cout << std::endl;
  std::cout << "-> MT2Analysis '" << name_ << "' has the following regions: " << std::endl;

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) 
    std::cout << "  " << ((*it)->region)->getName() << std::endl;

}



template<class T>
MT2Region* MT2Analysis<T>::getRegion( float ht, int njets, int nbjets, float mt, float mt2 ) const {

  MT2Region* foundRegion = 0;
  
  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {


    float htMin  = (*it)->region->htRegion()->htMin;
    float htMax  = (*it)->region->htRegion()->htMax;
    //float metMin = (*it)->region->htRegion()->metMin( ht );
    
    //float metMin = (*it)->region->htRegion()->metMin();
    //bool isInclusiveHT = (*it)->region->htRegion()->isInclusiveHT();
    //float metMinInclusiveHT = (*it)->region->htRegion()->metMinInclusiveHT( ht );

    
    if( ht<htMin ) continue;
    if( htMax>0. && ht>htMax ) continue;
    //if( metMin>0. && met>0. && met<metMin ) continue;
    
    //if( !(isInclusiveHT) && metMin>0. && met>0. && met<metMin ) continue;
    //if( metMinInclusiveHT>0. && met>0. && met<metMinInclusiveHT ) continue;
    //if( isInclusiveHT && metMinInclusiveHT>0. && met>0. && met<metMinInclusiveHT ) continue;

    int njetsmin  = (*it)->region->sigRegion()->nJetsMin;
    int njetsmax  = (*it)->region->sigRegion()->nJetsMax;
    int nbjetsmin = (*it)->region->sigRegion()->nBJetsMin;
    int nbjetsmax = (*it)->region->sigRegion()->nBJetsMax;

    if( njetsmin >=0 && njets <njetsmin  ) continue;
    if( njetsmax >=0 && njets >njetsmax  ) continue;
    if( nbjetsmin>=0 && nbjets<nbjetsmin ) continue;
    if( nbjetsmax>=0 && nbjets>nbjetsmax ) continue;

    std::string mtCut = (*it)->region->sigRegion()->mtCut;

    if( mtCut == "loMT" ) {

      float mtMin = 0.;
      float mtMax = 200.;
      float mt2Min = 200.;
      float mt2Max = 400.;
      bool insideMT  = ( mt  >= mtMin && mt <= mtMax );
      bool insideMT2 = ( mt2 >= mt2Min && mt2 <= mt2Max );
      bool insideBox = insideMT && insideMT2;
      if( !insideBox ) continue;
 
    } else if( mtCut == "hiMT" ) {

      float mtMin = 0.;
      float mtMax = 200.;
      float mt2Min = 200.;
      float mt2Max = 400.;
      bool insideMT  = ( mt  >= mtMin && mt <= mtMax );
      bool insideMT2 = ( mt2 >= mt2Min && mt2 <= mt2Max );
      bool insideBox = insideMT && insideMT2;
      if( insideBox ) continue; // outside

    } else if( mtCut != "" ) {

      std::cout << "WARNING! Unknown mtCut '" << mtCut << "'! Will ignore." << std::endl;

    }

    foundRegion = (*it)->region;
    break;

  }  // for

  return foundRegion;

}



template<class T>
MT2Region* MT2Analysis<T>::matchRegion( MT2Region region ) const {

  MT2Region* foundRegion = 0;
  
  for( std::set<MT2Region>::iterator iR=regions_.begin(); iR!=regions_.end(); ++iR ) {
    
    MT2Region* thisRegion= new MT2Region( (*iR) );
 
    if(!( region.isIncluded( thisRegion ) ) ) continue;
    foundRegion = ( thisRegion );
    break;

    delete thisRegion;
    
  }
    
  return foundRegion;


}





template<class T>
T* MT2Analysis<T>::get( const MT2Tree& mt2tree ) const {

  return this->get( mt2tree.ht, mt2tree.nJet30, mt2tree.nBJet20, mt2tree.minMTBMet, mt2tree.mt2 );
  //return this->get( mt2tree.ht, mt2tree.nJet30, mt2tree.nBJet20, mt2tree.met_pt, mt2tree.minMTBMet, mt2tree.mt2 );

}


template<class T>
T* MT2Analysis<T>::get( float ht, int njets, int nbjets, float mt, float mt2 ) const {
//T* MT2Analysis<T>::get( float ht, int njets, int nbjets, float met, float mt, float mt2 ) const {

  MT2Region* foundRegion = this->getRegion(ht, njets, nbjets, mt, mt2);
  //MT2Region* foundRegion = this->getRegion(ht, njets, nbjets, met, mt, mt2);

  if( foundRegion==0 ) return 0;

  return this->get( *foundRegion );

}


template<class T>
T* MT2Analysis<T>::get( const MT2Region& r ) const {

  bool found = false;

  T* t = 0;

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {
    if( *((*it)->region) == r ) {
      t = (*it);
      found = true;
      break;
    }
  }

  if( !found ) return 0;

  return t;

}



template<class T>
T* MT2Analysis<T>::getWithMatch( const MT2Region& r ) const {

  MT2Region* matchedRegion = this->matchRegion(r);
  if( matchedRegion==0 ) return 0;
  return this->get(*matchedRegion);

}


template<class T> 
void MT2Analysis<T>::setName( const std::string& newName ) {

  this->name_ = newName;

  for( std::set<MT2Region>::iterator iR=regions_.begin(); iR!=regions_.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t = this->get(thisRegion); 

    t->setName( newName );

  }

}



    

// operator overloading:

//template<class T> 
//template<class T2> 
//const MT2Analysis<T>& MT2Analysis<T>::operator=( const MT2Analysis<T2>& rhs ) {
//
//  regions_ = rhs.getRegions();
//
//  for( std::set<MT2Region>::iterator iR=regions_.begin(); iR!=regions_.end(); ++iR ) {
//
//    MT2Region thisRegion(*iR);
//
//    T* t1 = this->get(thisRegion); 
//    T2* t2 = rhs.get(thisRegion); 
//    if( t2==0 ) {
//      std::cout << "[MT2Analysis::operator=] ERROR! Can't equate MT2Analysis with different regional structures!" << std::endl;
//      exit(111);
//    }
//
//    *t1 = *t2;
//
//    //if( t1==0 ) {
//    //  t1 = new T(*t2);
//    //} else {
//    //  *t1 = *t2;
//    //}
//
//  }
//
//  return *this;
//
//}


// wonder why the above doesnt work (it apparently breaks MT2Estimate *= 0.92 in computeZinvFromGamma)
template<class T> 
const MT2Analysis<T>& MT2Analysis<T>::operator=( const MT2Analysis<T>& rhs ) {

  regions_ = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions_.begin(); iR!=regions_.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator=] ERROR! Can't equate MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    if( t1==0 ) {
      t1 = new T(*t2);
    } else {
      *t1 = *t2;
    }

  }

  return *this;

}



template<class T> 
template<class T2> 
MT2Analysis<T> MT2Analysis<T>::operator+( const MT2Analysis<T2>& rhs ) const {

  std::set<MT2Region> regions = rhs.getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

      MT2Region thisRegion(*iR);

      T* t1 = this->get(thisRegion); 
      T2* t2 = rhs.get(thisRegion); 
      if( t2==0 ) {
        std::cout << "[MT2Analysis::operator+] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
        exit(111);
      }

      T* tnew = new T(*t1 + *t2);
      newdata.insert( tnew ); 

  }

  MT2Analysis<T> result(name_, newdata);

  return result;

}



template<class T> 
template<class T2> 
MT2Analysis<T> MT2Analysis<T>::operator-( const MT2Analysis<T2>& rhs ) const {

  std::set<MT2Region> regions = rhs.getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

      MT2Region thisRegion(*iR);

      T* t1 = this->get(thisRegion); 
      T2* t2 = rhs.get(thisRegion); 
      if( t2==0 ) {
        std::cout << "[MT2Analysis::operator-] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
        exit(111);
      }

      T* tnew = new T(*t1 - *t2);
      newdata.insert( tnew ); 

  }

  MT2Analysis<T> result(name_, newdata);

  return result;

}



template<class T> 
template<class T2> 
const MT2Analysis<T>& MT2Analysis<T>::operator+=( const MT2Analysis<T2>& rhs ) {

  std::set<MT2Region> regions = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator+=] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    *t1 += *t2;

  }


  return *this;

}


template<class T> 
template<class T2> 
const MT2Analysis<T>& MT2Analysis<T>::operator-=( const MT2Analysis<T2>& rhs ) {

  std::set<MT2Region> regions = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T2* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator-=] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    *t1 -= *t2;

  }


  return *this;

}


template<class T> 
template<class T2> 
const MT2Analysis<T>& MT2Analysis<T>::operator/=( const MT2Analysis<T2>& rhs ) {

  std::set<MT2Region> regions = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T2* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator/= ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    *t1 /= *t2;

  }


  return *this;

}



template<class T> 
template<class T2> 
const MT2Analysis<T>& MT2Analysis<T>::operator*=( const MT2Analysis<T2>& rhs ) {


  std::set<MT2Region> regions = rhs.getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T2* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator*= ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    *t1 *= *t2;

  }


  return *this;

}




template<class T> 
const MT2Analysis<T>& MT2Analysis<T>::operator/=( float k ) {

  std::set<MT2Region> regions = this->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    *t1 /= k;

  }


  return *this;

}



template<class T> 
const MT2Analysis<T>& MT2Analysis<T>::operator*=( float k ) {


  std::set<MT2Region> regions = this->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    *t1 *= k;

  }


  return *this;

}





template<class T>
void MT2Analysis<T>::add( const MT2Analysis& rhs ) {

  (*this) = (*this) + rhs;

}




template<class T>
template<class T2>
MT2Analysis<T> MT2Analysis<T>::operator/( const MT2Analysis<T2>& rhs ) const {

  std::set<MT2Region> regions = rhs.getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T2* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator/] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    T* tnew = new T(*t1 / *t2);
    newdata.insert( tnew ); 

  }

  MT2Analysis<T> result(name_, newdata);

  return result;

}





template<class T>
void MT2Analysis<T>::divide( const MT2Analysis& rhs ) {

  (*this) = (*this) / rhs;

}




template<class T>
template<class T2>
MT2Analysis<T> MT2Analysis<T>::operator*( const MT2Analysis<T2>& rhs ) const {

  std::set<MT2Region> regions = rhs.getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = this->get(thisRegion); 
    T2* t2 = rhs.get(thisRegion); 
    if( t2==0 ) {
      std::cout << "[MT2Analysis::operator*] ERROR! Can't add MT2Analysis with different regional structures!" << std::endl;
      exit(111);
    }

    T* tnew = new T(*t1 * (*t2) );
    newdata.insert( tnew ); 

  }

  MT2Analysis<T> result(name_, newdata);

  return result;

}





template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator*( float k ) const {


  std::set<MT2Region> regions = this->getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = new T(*(this->get(thisRegion))); 
    *t1 *= k;

    newdata.insert( t1 );

  }

  MT2Analysis<T> result(name_, newdata);

  return result;

}



template<class T> 
MT2Analysis<T> MT2Analysis<T>::operator/( float k ) const {


  std::set<MT2Region> regions = this->getRegions();

  std::set<T*> newdata;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion(*iR);

    T* t1 = new T(*(this->get(thisRegion))); 
    *t1 /= k;

    newdata.insert( t1 );

  }

  MT2Analysis<T> result(name_, newdata);

  return result;

}






template<class T> 
void MT2Analysis<T>::writeToFile( const std::string& fileName, const std::string& option, bool overwrite ) {

  TFile* file = TFile::Open(fileName.c_str(), option.c_str() );
  file->cd();

  if( file->GetDirectory(this->name_.c_str()) ) {
    file->cd();
    if( overwrite ) {
      file->rmdir(this->name_.c_str());
    } else {
      std::cout << "[MT2Analysis::writeToFile] Directory '" << this->name_ << "' already exists in file '" << fileName << "'. Will not overwrite." << std::endl;
      return;
    }
  }

  file->mkdir(this->name_.c_str());
  file->cd(this->name_.c_str());

  std::set<MT2Region> regions = this->getRegions();
  for( std::set<MT2Region>::iterator it=regions.begin(); it!=regions.end(); ++it ) {
    file->cd();
    file->mkdir(Form("%s/%s", this->name_.c_str(), it->getName().c_str()) );
  }

  
  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) {
    file->cd();
    file->cd(Form("%s/%s", this->name_.c_str(), (*it)->region->getName().c_str()) );
    (*it)->write();
  }

  file->Close();

  std::cout << "-> Wrote '" << this->name_ << "' to file: " << fileName << std::endl;

}




template<class T>
void MT2Analysis<T>::printFromFile( const std::string& fileName, const std::string& ofs, const std::string& matchName ) {

  std::vector<MT2Analysis<T>*> analyses = readAllFromFile(fileName, matchName, false);

  if( analyses.size()==0 ) {
    std::cout << "[MT2Analysis::printFromFile] WARNING!!! Didn't find any MT2Analysis in file " << fileName << std::endl;
    return 0;
  }
  else
    print( analyses, ofs, matchName );

}



template<class T>
void MT2Analysis<T>::print( std::vector<MT2Analysis<T>*> analyses, const std::string& ofs, const std::string& matchName ) {

  MT2Analysis<T>* analysis = 0;
  if( matchName=="" ) {
    analysis = *(analyses.begin());
  } else {
    for( typename std::vector<MT2Analysis<T>*>::iterator iAn=analyses.begin(); iAn!=analyses.end(); ++iAn ) {
      if( (*iAn)->name_ == matchName ) {
        analysis = new MT2Analysis<T>(*(*iAn));
        break;
      }
    }
    if( analysis==0 ) {
      std::cout << "[MT2Analysis::print] WARNING!!! Didn't find any MT2Analysis named '" << matchName << std::endl;
      return 0;
    }
  }

  if( analyses.size()>1 && matchName=="" ) {
    std::cout << "[MT2Analysis::print] WARNING!!! Multiple analyses found, but reading only one ('" << analysis->name_ << "')" << std::endl;
  } else {
    std::cout << "[MT2Analysis::print] Grabbed MT2Analysis '" << analysis->name_ << std::endl;
  }

  analysis->print( ofs );

}


template<class T>
void MT2Analysis<T>::print( const std::string& ofs, MT2Region* matchRegion ) const {

  std::ifstream isExist(ofs);
  if(isExist) {
    
    system( Form("rm %s", ofs.c_str()) );
    //std::cout << "File " << ofs.c_str() << " already exist!" << std::endl;
    //exit(1);

  }

  std::ofstream ofs_file;
  if (ofs_file)
    ofs_file.open( ofs, std::ofstream::app );
  
  std::set<MT2HTRegion> htRegions = this->getHTRegions();
  std::set<MT2SignalRegion> sigRegions = this->getSignalRegions();
  std::set<MT2Region> mt2Regions = this->getRegions();

  std::string oldName="";

  for( std::set<MT2HTRegion>::iterator iHT=htRegions.begin(); iHT!=htRegions.end(); ++iHT ) {
    
    if( iHT->getNiceName() == oldName ) continue;
    if( matchRegion && !(iHT->isIncluded( matchRegion->htRegion() )) ) continue;
    
    std::string htRegionName = iHT->getNiceName();
    ofs_file << htRegionName << std::endl;

    for ( std::set<MT2SignalRegion>::iterator iSR=sigRegions.begin(); iSR!=sigRegions.end(); ++iSR ){
      
      if( matchRegion && !(iSR->isIncluded( matchRegion->sigRegion() )) ) continue;

      std::string sigRegionName = iSR->getNiceName();
      ofs_file << " & "  << sigRegionName;

    }
    
    ofs_file << std::endl;

    for( std::set<MT2Region>::iterator imt2=mt2Regions.begin(); imt2!=mt2Regions.end(); ++imt2 ) {
             
      if( *(imt2->htRegion()) != (*iHT) ) continue;
      if( matchRegion && !(imt2->isIncluded(matchRegion)) ) continue;

      T* thisT = this->get(*imt2);
      thisT->print( ofs );

    } // for mt2 regions                                                                                                                                                          

    ofs_file << std::endl << std::endl;
    
    oldName = iHT->getNiceName();
 
  } // for ht regions

  std::cout << "-> Printed analysis '" << name_ << "' to: " << ofs << std::endl;
  
}


template<class T>
void MT2Analysis<T>::print( std::ofstream& ofs_file, MT2Region* thisRegion ) const {

  int nBins;
  double* bins;
  thisRegion->getBins(nBins, bins);
  
  T* thisT = this->get(*thisRegion);

  for(int i=1; i < nBins+1; ++i){ 
    thisT->print( ofs_file, i );
  }
  
  ofs_file << "\\\\" << std::endl;
  
  //std::cout << "-> Printed analysis '" << name << "', HT region '"<< HTname << "' to: " << ofs << std::endl;
  
}


template<class T>
void MT2Analysis<T>::print( std::ofstream& ofs_file, MT2HTRegion* thisHTRegion ) const {

  std::set<MT2SignalRegion> sigRegions = this->getSignalRegions();
  for( std::set<MT2SignalRegion>::iterator iSig=sigRegions.begin(); iSig!=sigRegions.end(); ++iSig ) {
    
    MT2Region* thisRegion = new MT2Region(*thisHTRegion, *iSig);
    
    T* thisT = this->get(*thisRegion);
    thisT->print( ofs_file );
    
  } // for signal regions                                                                                                                                                          
  
  ofs_file << "\\\\" << std::endl;
  
  //std::cout << "-> Printed analysis '" << name << "', HT region '"<< HTname << "' to: " << ofs << std::endl;
  
}



template<class T> 
std::vector<MT2Analysis<T>*> MT2Analysis<T>::readAllFromFile( const std::string& fileName, const std::string& matchName, bool verbose ) {

  TFile* file = TFile::Open(fileName.c_str());
 
  if( file==0 ) {
    std::cout << "[MT2Analysis::readAllFromFile] ERROR! Can't open file: " << fileName << std::endl;
    exit(1357);
  }

  if( verbose ) std::cout << "[MT2Analysis] Reading analyses from file: " << file->GetName() << std::endl;

  std::vector<MT2Analysis<T>*> analyses;


  TIter next(file->GetListOfKeys());

  // these are the uppermost dirs in the file, 
  // so one dir per analysis name (tyipically only one anyways)
  while(TObject *obj = next()) { 

    TString thisdirname( obj->GetName() );
    if( thisdirname.BeginsWith("ProcessID") ) continue;

    std::string analysisName(obj->GetName());
    file->cd(analysisName.c_str()); 

    std::set<MT2Region> regions;

    TIter next2(gDirectory->GetListOfKeys());

    // these are the directiories inside the analysis dir
    // there will be one directory per region
    // and the dir name is the region identifying name
    while(TObject *obj2 = next2()) { 

      std::string regionName(obj2->GetName());
      MT2Region region(regionName);
      regions.insert(region);

      
      file->cd(analysisName.c_str());

    } // while regions


    //TString analysisName_tstr(analysisName);
    //if( matchExpression!="" && !(analysisName_tstr.Contains(matchExpression)) ) continue;
    TString analysisName_tstr(analysisName);
    if( (matchName=="SMS" || matchName=="DarkMatter" || matchName=="Zprime" || matchName=="Wprime") && !(analysisName_tstr.Contains(matchName)) ) continue;
    else if( matchName!="" && matchName!="SMS" && matchName!="DarkMatter" && matchName!="Zprime" && matchName!="Wprime" && matchName!=analysisName ) continue;

    // now that we know name and region structure we can istantiate an MT2Analysis:
    MT2Analysis<T>* analysis = new MT2Analysis<T>( analysisName, regions );

    // second loop to set everything
    file->cd(analysisName.c_str()); 
    TIter nextAgain(gDirectory->GetListOfKeys());

    while(TObject *obj2 = nextAgain()) { // loop on regions

      std::string regionName(obj2->GetName());
      MT2Region region(regionName);
      
      std::string path = analysisName + "/" + regionName;
      file->cd(path.c_str());

      T* thisT = analysis->get(region);
      thisT->getShit( file, path );

    } // while regions

    analyses.push_back( analysis );

    if( verbose ) std::cout << "  -> added: " << analysis->name_ << std::endl;

  } // while analysis names


  return analyses;
  
}



template<class T> 
MT2Analysis<T>* MT2Analysis<T>::readFromFile( const std::string& fileName, const std::string& matchName ) {


  std::vector<MT2Analysis<T>*> analyses = readAllFromFile(fileName, matchName, false);

  if( analyses.size()==0 ) {
    if( matchName!="" )
      std::cout << "[MT2Analysis::readFromFile] WARNING!!! Didn't find any MT2Analysis matching '" << matchName << "' in file " << fileName << std::endl;
    else
      std::cout << "[MT2Analysis::readFromFile] WARNING!!! Didn't find any MT2Analysis in file " << fileName << std::endl;
    return 0;
  }


  MT2Analysis<T>* analysis = *(analyses.begin());

  if( analyses.size()>1 && matchName=="" ) {
    std::cout << "[MT2Analysis::readFromFile] WARNING!!! Multiple analyses found in file: " << fileName << std::endl;
    std::cout << "[MT2Analysis::readFromFile] but reading only one ('" << analysis->name_ << "')" << std::endl;
    std::cout << "[MT2Analysis::readFromFile] (if you want to read all of them you should use readAllFromFile)" << std::endl;
  } else {
    std::cout << "[MT2Analysis::readFromFile] Grabbed MT2Analysis '" << analysis->name_ << "' from file " << fileName << std::endl;
  }

  return analysis;

}




template<class T> 
void MT2Analysis<T>::finalize() {

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) 
    (*it)->finalize();


}

template<class T> 
void MT2Analysis<T>::randomizePoisson( float scale ) {

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) 
    (*it)->randomizePoisson( scale );

}

template<class T> 
void MT2Analysis<T>::sqrtErrors( float scale ) {

  for( typename std::set<T*>::iterator it=data.begin(); it!=data.end(); ++it ) 
    (*it)->sqrtErrors( scale );

}


// global functions:

template<class T>
MT2Analysis<T> operator*( float k, const MT2Analysis<T>& rhs ) {

  return rhs*k;

}


template<class T>
MT2Analysis<T> operator/( float k, const MT2Analysis<T>& rhs ) {

  return rhs/k;

}



#endif
