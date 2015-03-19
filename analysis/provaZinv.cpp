#include <iostream>

#include "../interface/MT2Analysis.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Sample.h"

#define mt2_cxx
#include "../interface/mt2.h"


float lumi = 5.;


//MT2Analysis<MT2EstimateZinvGamma> computeYield( const MT2Sample& sample, const std::string& regionsSet );
void computeYield( const MT2Sample& sample, MT2Analysis<MT2EstimateZinvGamma>* analysis );


int main() {

  std::string samplesFileName = "../samples/samples_PHYS14_v2_Zinv.dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 200, 210); // only GJet
  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }


  std::string regionsSet = "13TeV_inclusive";

  MT2Analysis<MT2EstimateZinvGamma>*  analysis = new MT2Analysis<MT2EstimateZinvGamma>( "prova", regionsSet );
  for( unsigned i=0; i<fSamples.size(); ++i ) {
    //MT2Analysis<MT2EstimateZinvGamma> newana = ( computeYield( fSamples[i], regionsSet ) );
    //MT2Analysis<MT2EstimateZinvGamma>* newana = new MT2Analysis<MT2EstimateZinvGamma>( computeYield( fSamples[i], regionsSet ) );
    //(*analysis) += newana;
    computeYield( fSamples[i],  analysis );
  }

  analysis->writeToFile("provaZinv.root");


  MT2Analysis<MT2EstimateZinvGamma>* analysis2 = new MT2Analysis<MT2EstimateZinvGamma>( (*analysis)*2. );
  analysis2->setName("prova2");
  analysis2->writeToFile("provaZinv2.root");

  MT2Analysis<MT2EstimateZinvGamma>* analysis3 = new MT2Analysis<MT2EstimateZinvGamma>( "prova3", regionsSet );
  (*analysis3) = (*analysis) + 1000.*(*analysis2);
  analysis3->writeToFile("provaZinv3.root");

  return 0;

}



void computeYield( const MT2Sample& sample, MT2Analysis<MT2EstimateZinvGamma>* analysis ) {
//MT2Analysis<MT2EstimateZinvGamma> computeYield( const MT2Sample& sample, const std::string& regionsSet, MT2Analysis<MT2EstimateZinvGamma>* analysis ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
  

  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);



  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  //MT2Analysis<MT2EstimateZinvGamma> analysis( sample.sname, regionsSet, sample.id );
  //MT2Analysis<MT2EstimateTree>* analysis = new MT2Analysis<MT2EstimateTree>( sample.sname, regionsSet, sample.id );




  int nentries = tree->GetEntries();
  nentries = 10000;

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;


    myTree.GetEntry(iEntry);

    if( myTree.nMuons10 > 0) continue;
    if( myTree.nElectrons10 > 0 ) continue;
    if( myTree.nPFLep5LowMT > 0) continue;
    if( myTree.nPFHad10LowMT > 0) continue;
  
    if( myTree.deltaPhiMin<0.3 ) continue;
    if( myTree.diffMetMht>0.5*myTree.met_pt ) continue;

    if( myTree.nVert==0 ) continue;
    if( myTree.nJet40<2 ) continue;
    if( myTree.njet<2 ) continue;
    if( myTree.jet_pt[1]<100. ) continue;


    MT2EstimateZinvGamma* thisEstimate = analysis->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt, myTree.gamma_minMTBMet, myTree.gamma_mt2 );

    if( thisEstimate==0 ) continue;

    float weight = myTree.evt_scale1fb*lumi;
    thisEstimate->fillIso( myTree.gamma_chHadIso[0], weight, myTree.gamma_mt2 );
    thisEstimate->yield->Fill(myTree.gamma_mt2, weight);

    
  } // for entries

  analysis->finalize();


  delete tree;

  file->Close();
  delete file;
  

}


