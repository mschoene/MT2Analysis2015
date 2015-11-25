#ifndef MT2EstimateTree_h
#define MT2EstimateTree_h

#include "MT2Estimate.h"
#include "MT2Analysis.h"
#include "TTree.h"
#include "mt2.h"


class MT2EstimateTree : public MT2Estimate {

 public:

  MT2EstimateTree( const std::string& aname, const MT2Region& aregion );
  MT2EstimateTree( const MT2EstimateTree& rhs );
  virtual ~MT2EstimateTree();

  void initTree();
  virtual void setName( const std::string& newName );

  void assignTree( const MT2Tree& mt2tree, float w );
  void fillTree( const MT2Tree& mt2tree, float w );

  void assignTree_gamma( const MT2Tree& mt2tree, float w );
  void fillTree_gamma( const MT2Tree& mt2tree, float w );

  void assignTree_zll( const MT2Tree& mt2tree, float w );
  void fillTree_zll( const MT2Tree& mt2tree, float w );

  static void addVar( MT2Analysis<MT2EstimateTree>* analysis, const std::string& name );
  //static void addVarFloat( MT2Analysis<MT2EstimateTree>* analysis, const std::string& name );
  //static void addVarInt( MT2Analysis<MT2EstimateTree>* analysis, const std::string& name );
  void assignVars( float aht, int anJets, int anBJets, float amet, float amt2 );
  void assignVar( const std::string& name, float value );
  //void assignVar( const std::string& name, int value );

  int run;
  int lumi;
  int evt;
  float weight; // = crossSecWeight * puWeight
  float puWeight; // the puWeight part of the above
  int id;

  float mt2;
  float ht;
  float met;

  float deltaPhiMin;
  float diffMetMht;
  int nVert;

  int nJets;
  int nBJets;
  int nElectrons;
  int nMuons;
  int nPFLep;
  int nPFHad;
  int nJetHF;

  int GenSusyMScan1;
  int GenSusyMScan2;

//  float GenSusyMScan1;
//  float GenSusyMScan2;

//  float LHEweight_original;
//  int LHEweight_id[446];
//  float LHEweight_wgt[446];
  

  //std::map< std::string, size_t > extraVars;
  std::map< std::string, float* > extraVars;
 
  TTree* tree;

  const MT2EstimateTree& operator=( const MT2EstimateTree& rhs );
  MT2EstimateTree operator+( const MT2EstimateTree& rhs ) const;
  const MT2EstimateTree& operator+=( const MT2EstimateTree& rhs );

  MT2EstimateTree operator/ ( float k ) const;
  MT2EstimateTree operator* ( float k ) const;
  const MT2EstimateTree& operator/=( float k );
  const MT2EstimateTree& operator*=( float k );

  friend MT2EstimateTree operator*( float k, const MT2EstimateTree& rhs );
  friend MT2EstimateTree operator/( float k, const MT2EstimateTree& rhs );

  virtual void finalize() {
  }

  void projectFromTree( const MT2EstimateTree* treeEst, const std::string& selection, const std::string& variable="mt2" );
  //  void projectFromTree( const MT2EstimateTree* treeEst, const std::string& selection );
  
  //  static MT2Analysis<MT2EstimateTree>* makeRebinnedAnalysisFromInclusiveTree( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* analysis, const std::string& selection="", int nBins=0, float xMin=0., float xMax=0. );
 
  static MT2Analysis<MT2EstimateTree>* makeAnalysisFromInclusiveTree( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* analysis, const std::string& selection="" ) { return makeRebinnedAnalysisFromInclusiveTree( aname, regionsSet, analysis, selection ); };

  static MT2Analysis<MT2EstimateTree>* makeRebinnedAnalysisFromInclusiveTree( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* analysis, const std::string& selection="", int nBins=0, double* bins=0,  const std::string& variable="mt2" );

  virtual void getShit( TFile* file, const std::string& path );

  virtual void write() const;

  virtual void print(const std::string& ofs);

 private:

};





#endif
