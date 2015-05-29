#ifndef MT2EstimateSigTree_h
#define MT2EstimateSigTree_h

#include "MT2EstimateSig.h"
#include "MT2Analysis.h"
#include "TTree.h"
#include "mt2.h"


class MT2EstimateSigTree : public MT2EstimateSig {

 public:

  MT2EstimateSigTree( const std::string& aname, const MT2Region& aregion );
  MT2EstimateSigTree( const MT2EstimateSigTree& rhs );
  virtual ~MT2EstimateSigTree();

  void initTree();
  virtual void setName( const std::string& newName );

  void assignTree( const MT2Tree& mt2tree, float w );
  void fillTree( const MT2Tree& mt2tree, float w );

  void assignTree_gamma( const MT2Tree& mt2tree, float w );
  void fillTree_gamma( const MT2Tree& mt2tree, float w );

  static void addVar( MT2Analysis<MT2EstimateSigTree>* analysis, const std::string& name );
  //static void addVarFloat( MT2Analysis<MT2EstimateSigTree>* analysis, const std::string& name );
  //static void addVarInt( MT2Analysis<MT2EstimateSigTree>* analysis, const std::string& name );
  void assignVars( float aht, int anJets, int anBJets, float amet, float amt2 );
  void assignVar( const std::string& name, float value );
  //void assignVar( const std::string& name, int value );

  int run;
  int lumi;
  int evt;
  float weight;
  int id;

  float mt2;
  float ht;
  float met;
  int nJets;
  int nBJets;

  //std::map< std::string, size_t > extraVars;
  std::map< std::string, float* > extraVars;
 
  TTree* tree;

  const MT2EstimateSigTree& operator=( const MT2EstimateSigTree& rhs );
  MT2EstimateSigTree operator+( const MT2EstimateSigTree& rhs ) const;
  const MT2EstimateSigTree& operator+=( const MT2EstimateSigTree& rhs );

  MT2EstimateSigTree operator/ ( float k ) const;
  MT2EstimateSigTree operator* ( float k ) const;
  const MT2EstimateSigTree& operator/=( float k );
  const MT2EstimateSigTree& operator*=( float k );

  friend MT2EstimateSigTree operator*( float k, const MT2EstimateSigTree& rhs );
  friend MT2EstimateSigTree operator/( float k, const MT2EstimateSigTree& rhs );



  virtual void finalize() {
  }


  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

  virtual void print(const std::string& ofs);

 private:

};





#endif
