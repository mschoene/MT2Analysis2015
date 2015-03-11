#include "../interface/MT2EstimateTree.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>




MT2EstimateTree::MT2EstimateTree( const std::string& aname, const MT2Region& aregion ) : MT2Estimate( aname, aregion ) {

  TH1::AddDirectory(kFALSE);

  tree = new TTree(this->getHistoName("tree").c_str(), "" );
  
  this->initTree();

}




MT2EstimateTree::MT2EstimateTree( const MT2EstimateTree& rhs ) : MT2Estimate( rhs ) {

  TH1::AddDirectory(kFALSE);

  tree = rhs.tree->CloneTree(-1);

  this->initTree();

}



MT2EstimateTree::~MT2EstimateTree() {

  delete tree;

}




void MT2EstimateTree::initTree() {

  tree->Branch( "run", &run, "run/I");
  tree->Branch( "lumi", &lumi, "lumi/I");
  tree->Branch( "evt", &evt, "evt/i");
  tree->Branch( "weight", &weight, "weight/F");
  tree->Branch( "id", &id, "id/I");

  tree->Branch( "mt2", &mt2, "mt2/F");
  tree->Branch( "ht", &ht, "ht/F");
  tree->Branch( "met", &met, "met/F");
  tree->Branch( "nJets", &nJets, "nJets/I");
  tree->Branch( "nBJets", &nBJets, "nBJets/I");

  tree->SetDirectory(0);

}



void MT2EstimateTree::setName( const std::string& newName ) {

  MT2Estimate::setName(newName);

  tree->SetName( this->getHistoName("tree").c_str() );

}



void MT2EstimateTree::addVar( MT2Analysis<MT2EstimateTree>* analysis, const std::string& name ) {


  for( std::set<MT2EstimateTree*>::iterator iD=analysis->data.begin(); iD!=analysis->data.end(); ++iD ) {

    float* x = new float();
    (*iD)->extraVars[name] = x;

    (*iD)->tree->Branch( name.c_str(), x, Form("%s/F", name.c_str()) );

  }

}


//void MT2EstimateTree::addVarFloat( MT2Analysis<MT2EstimateTree>* analysis, const std::string& name ) {
//
//
//  for( std::set<MT2EstimateTree*>::iterator iD=analysis->data.begin(); iD!=analysis->data.end(); ++iD ) {
//
//    float* x = new float();
//    (*iD)->extraVars[name] = (size_t)x;
//
//    (*iD)->tree->Branch( name.c_str(), x, Form("%s/F", name.c_str()) );
//
//  }
//
//}



//void MT2EstimateTree::addVarInt( MT2Analysis<MT2EstimateTree>* analysis, const std::string& name ) {
//
//
//  for( std::set<MT2EstimateTree*>::iterator iD=analysis->data.begin(); iD!=analysis->data.end(); ++iD ) {
//
//    int* x = new int();
//    (*iD)->extraVars[name] = (size_t)x;
//
//    (*iD)->tree->Branch( name.c_str(), x, Form("%s/I", name.c_str()) );
//
//  }
//
//}



void MT2EstimateTree::assignVar( const std::string& name, float value ) {

  //&(extraVars[name]) = value;
  float* x = (float*)extraVars[name];
  *x = value;

}



//void MT2EstimateTree::assignVar( const std::string& name, int value ) {
//
//  int* x = (int*)extraVars[name];
//  *x = value;
//
//}




void MT2EstimateTree::fillTree( const MT2Tree& mt2tree, float w ) {

  this->assignTree( mt2tree, w );

  tree->Fill();

}



void MT2EstimateTree::fillTree_gamma( const MT2Tree& mt2tree, float w ) {

  this->assignTree_gamma( mt2tree, w );

  tree->Fill();

}



void MT2EstimateTree::assignTree( const MT2Tree& mt2tree, float w ) {

  run    = mt2tree.run;
  lumi   = mt2tree.lumi;
  evt    = mt2tree.evt;
  weight = w;
  id     = mt2tree.evt_id;

  mt2    = mt2tree.mt2;
  ht     = mt2tree.ht;
  met    = mt2tree.met_pt;
  nJets  = mt2tree.nJet40;
  nBJets = mt2tree.nBJet40;

}
  

void MT2EstimateTree::assignTree_gamma( const MT2Tree& mt2tree, float w ) {

  run    = mt2tree.run;
  lumi   = mt2tree.lumi;
  evt    = mt2tree.evt;
  weight = w;
  id     = mt2tree.evt_id;

  mt2    = mt2tree.gamma_mt2;
  ht     = mt2tree.gamma_ht;
  met    = mt2tree.gamma_met_pt;
  nJets  = mt2tree.gamma_nJet40;
  nBJets = mt2tree.gamma_nBJet40;

}

  

void MT2EstimateTree::assignVars( float aht, int anJets, int anBJets, float amet, float amt2 ) {


  mt2    = amt2;
  ht     = aht;
  met    = amet;
  nJets  = anJets;
  nBJets = anBJets;

}
  




void MT2EstimateTree::getShit( TFile* file, const std::string& path ) {

  MT2Estimate::getShit(file, path);

  tree = (TTree*)file->Get(Form("%s/%s", path.c_str(), tree->GetName()));

  this->initTree();

}



void MT2EstimateTree::write() const {

  MT2Estimate::write();
  tree->Write();

}


void MT2EstimateTree::print(const std::string& ofs){

  MT2Estimate::print( ofs );

}


const MT2EstimateTree& MT2EstimateTree::operator=( const MT2EstimateTree& rhs ) {

  if( this->yield == 0 ) { // first time

    this->setName(rhs.getName());

    this->region = new MT2Region(*(rhs.region));

    this->yield = new TH1D(*(rhs.yield));

    this->tree = rhs.tree->CloneTree(-1);

  } else { // keep name and histo name, just make histogram identical

    if( this->region!=0 ) delete this->region;
    this->region = new MT2Region(*(rhs.region));

    std::string oldName = this->yield->GetName();
    delete this->yield;
    this->yield = new TH1D(*(rhs.yield));
    this->yield->SetName(oldName.c_str());

    std::string oldTreeName = this->tree->GetName();
    this->tree = rhs.tree->CloneTree(-1);
    this->tree->SetName(oldTreeName.c_str());

  }

  this->initTree();

  return *this;

}




MT2EstimateTree MT2EstimateTree::operator+( const MT2EstimateTree& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateTree::operator+] ERROR! Can't add MT2EstimateTree with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateTree result(*this);
  result.yield->Add(rhs.yield);

  TList* list = new TList;
  list->Add(result.tree);
  list->Add(rhs.tree);
  result.tree = TTree::MergeTrees( list );

  result.initTree();


  return result;
  
}






//MT2EstimateTree MT2EstimateTree::operator*( float k ) const{
//
//  MT2EstimateTree result( this->getName(), *(this->region) );
//  result.yield = new TH1D(*(this->yield));
//  result.yield->Scale(k);
//
//  // loop on entries and multiply weight
//
//  return result;
//
//}



//MT2EstimateTree MT2EstimateTree::operator/( float k ) const{
//
//  MT2EstimateTree result( this->getName(), *(this->region) );
//  result.yield = new TH1D(*(this->yield));
//  result.yield->Scale(1./k);
//  // loop on entries and multiply weight
//
//  return result;
//
//
//}



const MT2EstimateTree& MT2EstimateTree::operator+=( const MT2EstimateTree& rhs ) {

  if( rhs.tree->GetEntries()>0 ) {

    this->yield->Add(rhs.yield);

    std::string oldName(this->tree->GetName());

    TList* list = new TList;
    list->Add(this->tree);
    list->Add(rhs.tree);
    this->tree = TTree::MergeTrees( list );
    this->tree->SetName( oldName.c_str() );

    this->initTree();

  }

  return (*this);

}





//const MT2EstimateTree& MT2EstimateTree::operator*=( float k ) {
//
//  this->yield->Scale(k);
//  // loop on entries and multiply weight
//  return (*this);
//
//}
//
//const MT2EstimateTree& MT2EstimateTree::operator/=( float k ) {
//
//  this->yield->Scale(1./k);
//  // loop on entries and multiply weight
//  return (*this);
//
//}

