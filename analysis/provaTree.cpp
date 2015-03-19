#include <iostream>

#include "../interface/MT2Analysis.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Sample.h"

#define mt2_cxx
#include "../interface/mt2.h"


float lumi = 5.;


MT2Analysis<MT2EstimateTree> computeYield( const MT2Sample& sample, const std::string& regionsSet );


int main() {

  std::string samplesFileName = "../samples/samples_mctest.dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 1000, 2000); // only signal
  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }


  std::string regionsSet = "13TeV_inclusive";

  MT2Analysis<MT2EstimateTree>*  analysis = new MT2Analysis<MT2EstimateTree>( "provaTree", regionsSet );
  for( unsigned i=0; i<1; ++i ) {
    MT2Analysis<MT2EstimateTree> newana = ( computeYield( fSamples[i], regionsSet ) );
    (*analysis) += newana;
  }

  MT2Analysis<MT2EstimateTree>*  analysis2 = new MT2Analysis<MT2EstimateTree>( "provaTree2", regionsSet );
  for( unsigned i=1; i<fSamples.size(); ++i ) {
    MT2Analysis<MT2EstimateTree> newana = ( computeYield( fSamples[i], regionsSet ) );
    (*analysis2) += newana;
  }

  analysis->writeToFile("provaTree.root");
  analysis2->writeToFile("provaTree2.root");


  MT2Analysis<MT2EstimateTree>*  analysis2copy = new MT2Analysis<MT2EstimateTree>( *analysis2 );
  analysis2copy->writeToFile("provaTree2copy.root");
  analysis2copy->setName("provaTree2copy");
  analysis2copy->addToFile("provaTree2copy.root");

  analysis2copy->setName("provaTree2_times2");
  (*analysis2copy) *= 2.;
  analysis2copy->addToFile("provaTree2copy.root");


  MT2Analysis<MT2EstimateTree>*  analysisSum = new MT2Analysis<MT2EstimateTree>( "provaSum", regionsSet );
  (*analysisSum) = 100.*(*analysis) + (*analysis2);
  analysisSum->writeToFile("provaTreeSum.root");


  return 0;

}






MT2Analysis<MT2EstimateTree> computeYield( const MT2Sample& sample, const std::string& regionsSet ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
  

  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);



  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  MT2Analysis<MT2EstimateTree> analysis( sample.sname, regionsSet, sample.id );
  //MT2Analysis<MT2EstimateTree>* analysis = new MT2Analysis<MT2EstimateTree>( sample.sname, regionsSet, sample.id );
  MT2EstimateTree::addVar( &analysis, "puppa" );



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


    MT2EstimateTree* thisEstimate = analysis.get( myTree.ht, myTree.nJet40, myTree.nBJet40, myTree.met_pt, myTree.minMTBMet, myTree.mt2 );

    if( thisEstimate==0 ) continue;

    float weight = lumi*myTree.evt_scale1fb ;

    thisEstimate->yield->Fill( myTree.mt2, weight );

    thisEstimate->assignTree( myTree, weight );
    thisEstimate->assignVar( "puppa", -13. );
    thisEstimate->tree->Fill();
    //thisEstimate->fillTree( myTree, lumi*myTree.evt_scale1fb );

    
  } // for entries

  //ofs.close();

  analysis.finalize();
  
  delete tree;

  file->Close();
  delete file;
  
  return analysis;

}


