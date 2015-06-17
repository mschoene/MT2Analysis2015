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
  tree->Branch( "deltaPhiMin", &deltaPhiMin, "deltaPhiMin/F");
  tree->Branch("diffMetMht", &diffMetMht, "diffMetMht/F");
  tree->Branch( "nVert", &nVert, "nVert/I");
  tree->Branch( "nJets", &nJets, "nJets/I");
  tree->Branch( "nBJets", &nBJets, "nBJets/I");
  tree->Branch( "nElectrons", &nElectrons, "nElectrons/I");
  tree->Branch( "nMuons", &nMuons, "nMuons/I");
  tree->Branch( "nPFLep", &nPFLep, "nPFLep/I");
  tree->Branch( "nPFHad", &nPFHad, "nPFHad/I");

  tree->Branch( "GenSusyMScan1", &GenSusyMScan1, "GenSusyMScan1/I");
  tree->Branch( "GenSusyMScan2", &GenSusyMScan2, "GenSusyMScan2/I");

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

  mt2           = mt2tree.mt2;
  ht            = mt2tree.ht;
  met           = mt2tree.met_pt;
  deltaPhiMin   = mt2tree.deltaPhiMin;
  diffMetMht    = mt2tree.diffMetMht;
  nJets         = mt2tree.nJet40;
  nBJets        = mt2tree.nBJet20;
  nElectrons    = mt2tree.nElectrons10;
  nMuons        = mt2tree.nMuons10;
  nPFLep        = mt2tree.nPFLep5LowMT;
  nPFHad        = mt2tree.nPFHad10LowMT;
  
  GenSusyMScan1 = mt2tree.GenSusyMScan1;
  GenSusyMScan2 = mt2tree.GenSusyMScan2;
  
}
  

void MT2EstimateTree::assignTree_gamma( const MT2Tree& mt2tree, float w ) {

  run    = mt2tree.run;
  lumi   = mt2tree.lumi;
  evt    = mt2tree.evt;
  weight = w;
  id     = mt2tree.evt_id;

  mt2           = mt2tree.gamma_mt2;
  ht            = mt2tree.gamma_ht;
  met           = mt2tree.gamma_met_pt;
  dPhiMin       = mt2tree.gamma_dPhiMin;
  diffMetMht    = mt2tree.gamma_diffMetMht;
  nJets         = mt2tree.gamma_nJet40;
  nBJets        = mt2tree.gamma_nBJet20;
  nElectrons    = mt2tree.nElectrons10;
  nMuons        = mt2tree.nMuons10;
  nPFLep        = mt2tree.nPFLep5LowMT;
  nPFHad        = mt2tree.nPFHad10LowMT;

  GenSusyMScan1 = mt2tree.GenSusyMScan1;
  GenSusyMScan2 = mt2tree.GenSusyMScan2;

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

  this->region = new MT2Region(*(rhs.region));

  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield = new TH1D(*(rhs.yield));

  this->tree = rhs.tree->CloneTree(-1);

  this->setName( this->getName() );

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
  result.yield3d->Add(rhs.yield3d);

  TList* list = new TList;
  list->Add(result.tree);
  list->Add(rhs.tree);
  result.tree = TTree::MergeTrees( list );

  result.initTree();


  return result;
  
}






//MT2EstimateTree MT2EstimateTree::operator*( float k ) const{
//
//std::cout << "aaaaaaaaaaaa here! " << std::endl;
//  MT2EstimateTree result( this->getName(), *(this->region) );
//  result.yield = new TH1D(*(this->yield));
//  result.yield->Scale(k);
//
//  // loop on entries and multiply weight
//
//
//  std::string oldName(tree->GetName());
//  result.tree->SetName("willBeKilled");
//
//  result.tree->SetBranchStatus( "weight", 0 );
//  Float_t oldWeight;
//  result.tree->SetBranchAddress("weight", &oldWeight);
//
//  TTree* newTree = new TTree(oldName.c_str(), "");
//  newTree = result.tree->CloneTree(0);
//  Float_t newWeight;
//  newTree->Branch( "weight", &newWeight, "newWeight/F" );
// 
//  int nentries = result.tree->GetEntries();
//
//  // loop on entries and multiply weight
//  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
//
//    result.tree->GetEntry(iEntry);
//    newWeight = oldWeight*k;
//    newTree->Fill();
//
//  }
//
//
//  delete result.tree;
//  result.tree = newTree;
//  result.tree->SetName( oldName.c_str() );
//
//  return result;
//
//}




MT2EstimateTree MT2EstimateTree::operator*( float k ) const{

  MT2EstimateTree result( this->getName(), *(this->region) );
  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(k);
  
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(k);
  
  // loop on entries and multiply weight

  return result;


}


MT2EstimateTree MT2EstimateTree::operator/( float k ) const{

  MT2EstimateTree result( this->getName(), *(this->region) );
  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(1./k);
  
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(1./k);

  
  // loop on entries and multiply weight

  return result;


}



const MT2EstimateTree& MT2EstimateTree::operator+=( const MT2EstimateTree& rhs ) {


  if( rhs.tree->GetEntries()>0 ) {

    this->yield->Add(rhs.yield);
    this->yield3d->Add(rhs.yield3d);

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
//
//  this->yield->Scale(k);
//
//
//  std::string oldName(tree->GetName());
//  tree->SetName("willBeKilled");
//
//  TTree* newTree = new TTree(oldName.c_str(), "");
//
//  newTree->Branch( "run", &run, "run/I");
//  newTree->Branch( "lumi", &lumi, "lumi/I");
//  newTree->Branch( "evt", &evt, "evt/i");
//  newTree->Branch( "weight", &weight, "weight/F");
//  newTree->Branch( "id", &id, "id/I");
//
//  newTree->Branch( "mt2", &mt2, "mt2/F");
//  newTree->Branch( "ht", &ht, "ht/F");
//  newTree->Branch( "met", &met, "met/F");
//  newTree->Branch( "nJets", &nJets, "nJets/I");
//  newTree->Branch( "nBJets", &nBJets, "nBJets/I");
//
//
//  int nentries = tree->GetEntries();
//
//  // loop on entries and multiply weight
//  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
//
//    tree->GetEntry(iEntry);
//
//    weight *= k;
//
//    newTree->Fill();
//
//  }
//
//
//  delete tree;
//  tree = newTree;
//  tree->SetName( oldName.c_str() );
//
//  return (*this);
//
//}


const MT2EstimateTree& MT2EstimateTree::operator*=( float k ) {

  this->yield->Scale(k);
  this->yield3d->Scale(k);
  // loop on entries and multiply weight
  return (*this);

}


const MT2EstimateTree& MT2EstimateTree::operator/=( float k ) {

  this->yield->Scale(1./k);
  this->yield3d->Scale(1./k);
  // loop on entries and multiply weight
  return (*this);

}



// friend functions:

MT2EstimateTree operator*( float k, const MT2EstimateTree& rhs ) {

  return rhs*k;

}



MT2EstimateTree operator/( float k, const MT2EstimateTree& rhs ) {

  return rhs/k;

}

