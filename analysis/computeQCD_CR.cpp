#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateQCD.h"

#define mt2_cxx
#include "../interface/mt2.h"


float lumi = 4.; //fb-1

bool doDiffMetMht = true;


MT2Analysis<MT2EstimateQCD>* fillCR( const MT2Sample& sample, const std::string& regionsSet);
MT2Analysis<MT2EstimateQCD>* merge ( std::vector<MT2Analysis<MT2EstimateQCD> *> anas, const std::string& regionsSet, const std::string& name, int id_min, int id_max );



int main( int argc, char* argv[] ) {


  if( argc>2 ) {
    std::cout << "USAGE: ./qcdControlRegion [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  std::string samplesFileName = "PHYS14_qcdCR";
  if( argc>1 ) {
    std::string samplesFileName_tmp(argv[1]); 
    samplesFileName = samplesFileName_tmp;
  }


  std::string regionsSet = "zurich";

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string outputdir( Form("QCDcontrolRegion_%s_%s_%.0ffb", samplesFileName.c_str(), regionsSet.c_str(), lumi) );
  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string samplesFile = "../samples/samples_" + samplesFileName + ".dat";

  std::vector<MT2Sample> samples = MT2Sample::loadSamples(samplesFile, 1, 999); // not interested in signal here (see later)
  if( samples.size()==0 ) {
    std::cout << "There must be an error: didn't find any qcd files in " << samplesFile << "!" << std::endl;
    exit(1209);
  }


  std::vector< MT2Analysis<MT2EstimateQCD>* > controlRegions;
  for( unsigned i=0; i<samples.size(); ++i ) 
    controlRegions.push_back( fillCR( samples[i], regionsSet ) );

  MT2Analysis<MT2EstimateQCD>* qcdCR  = merge( controlRegions, regionsSet, "qcdCR", 100, 199 );

  qcdCR->finalize();

  //for( unsigned i=0; i<samples.size(); ++i ) {
  //  controlRegions[i]->finalize();
  //  controlRegions[i]->writeToFile( outputdir + "/mc.root" );
  //}
  //qcdCR->addToFile( outputdir + "/mc.root" );
  
  qcdCR->writeToFile( outputdir + "/mc.root" );

}


MT2Analysis<MT2EstimateQCD>* merge( std::vector<MT2Analysis<MT2EstimateQCD> *> anas, const std::string& regionsSet, const std::string& name, int id_min, int id_max ) {

  if( id_max<0 ) id_max=id_min;

  MT2Analysis<MT2EstimateQCD>* ana = new MT2Analysis<MT2EstimateQCD>(name, regionsSet);

  for( unsigned i=0; i<anas.size(); ++i ) {
    if( anas[i]->id >= id_min && anas[i]->id <= id_max )
      *(ana) += *(anas[i]);
  } 

  return ana;

}

MT2Analysis<MT2EstimateQCD>* fillCR( const MT2Sample& sample, const std::string& regionsSet){

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting to fill QCD CR for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  MT2Analysis<MT2EstimateQCD>* analysis = new MT2Analysis<MT2EstimateQCD>( sample.sname, regionsSet, sample.id );

  int nentries = tree->GetEntries();
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    myTree.GetEntry(iEntry);

    if ( !(myTree.passLeptonVeto() && myTree.passIsoTrackVeto()) )                            continue;
    if ( !(myTree.nVert>0 && myTree.nJet40>=2 && myTree.jet1_pt>40.  && myTree.jet2_pt>40.) ) continue;
    if ( doDiffMetMht && !(myTree.diffMetMht < 0.5*myTree.met_pt) )                           continue;

    float ht   = myTree.ht;
    float mt2  = myTree.mt2;
    float dphi = myTree.deltaPhiMin;
    int njets  = myTree.nJet40;
    int nbjets = myTree.nBJet40;

    Double_t weight = myTree.evt_scale1fb*lumi;

    MT2EstimateQCD* thisEstimate = analysis->get( ht, njets, nbjets );
    if( thisEstimate==0 ) continue;

    thisEstimate->fillDphi(dphi, weight, mt2);

  } // entries loop

  //analysis->finalize();  //will do after merge
  
  delete tree;

  file->Close();
  delete file;
  
  return analysis;

}
