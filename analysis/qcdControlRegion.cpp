#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"




int round(float d) {
  return (int)(floor(d + 0.5));
}


void computeYield( const MT2Sample& sample, const MT2Config& cfg,  MT2Analysis<MT2EstimateTree>* anaTree );




int main( int argc, char* argv[] ) {



  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|          Running qcdControlRegion         |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./qcdControlRegion [configFileName] [data/MC]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  
  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool onlyData = false;
  bool onlyMC   = false;
  if( argc > 2 ) {
    std::string dataMC(argv[2]);
    if( dataMC=="data" ) onlyData = true;
    else if( dataMC=="MC" || dataMC=="mc" ) onlyMC = true;
    else {
      std::cout << "-> You passed a second argument that isn't 'data' or 'MC', so I don't know what to do about it." << std::endl;
    }
  }



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = cfg.getEventYieldDir() + "/qcdControlRegion"; 
  system(Form("mkdir -p %s", outputdir.c_str()));


  if( cfg.useMC() && !onlyData ) { // run on MC

    std::string samplesFile = "../samples/samples_" + cfg.mcSamples() + ".dat";
    
    std::vector<MT2Sample> samples_wjet = MT2Sample::loadSamples(samplesFile, 502, 505);
    std::vector<MT2Sample> samples_top  = MT2Sample::loadSamples(samplesFile, 300, 499);
    std::vector<MT2Sample> samples_qcd  = MT2Sample::loadSamples(samplesFile, 100, 199);


    MT2Analysis<MT2EstimateTree>* qcdCRtree = new MT2Analysis<MT2EstimateTree>( "qcdCRtree", cfg.regionsSet() );
    
    
    for( unsigned i=0; i<samples_wjet.size(); ++i ) 
      computeYield( samples_wjet[i], cfg, qcdCRtree );
    for( unsigned i=0; i<samples_top.size(); ++i ) 
      computeYield( samples_top[i], cfg, qcdCRtree );
    for( unsigned i=0; i<samples_qcd.size(); ++i ) 
      computeYield( samples_qcd[i], cfg, qcdCRtree );
    

   
    std::string mcFile = outputdir + "/mc.root";
    qcdCRtree->writeToFile( mcFile, "RECREATE" );

  }


  if( !(cfg.dummyAnalysis()) && cfg.dataSamples()!="" && !onlyMC  ) {

    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, 1, 2 );
    if( samples_data.size()==0 ) {
      std::cout << "There must be an error: samples_data is empty!" << std::endl;
      exit(1209);
    }

    MT2Analysis<MT2EstimateTree>* data = new MT2Analysis<MT2EstimateTree>( "qcdCRtree", cfg.regionsSet() );
    for( unsigned i=0; i < samples_data.size(); ++i )
      computeYield( samples_data[i], cfg, data );

    std::string dataFile = outputdir + "/data.root";
    data->writeToFile( dataFile, "RECREATE" );


  }



  return 0;

}








void computeYield( const MT2Sample& sample, const MT2Config& cfg, MT2Analysis<MT2EstimateTree>* anaTree ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;


  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  int nentries = tree->GetEntries();


  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( myTree.isData ) {
      if ( myTree.isGolden == 0 ) continue;
      if( !(myTree.Flag_HBHENoiseFilter && myTree.Flag_CSCTightHaloFilter && myTree.Flag_eeBadScFilter) ) continue;
    }
    

    if( !myTree.passSelection("qcd") ) continue;


    float minMTBmet = myTree.minMTBMet;
    //float met       = myTree.met_pt;
    int njets       = myTree.nJet30;
    int nbjets      = myTree.nBJet20;    
    float mt2       = (njets>1) ? myTree.mt2 : myTree.jet1_pt;
    float ht        = myTree.ht;

    if (myTree.isData && ( (ht<575  && myTree.HLT_ht350prescale==0) ||
                 (ht>575 && ht<1000 && myTree.HLT_ht475prescale==0) ) )  continue;



    if( mt2<30. ) continue;

    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi(); 

    MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisTree==0 ) continue;


    thisTree->yield->Fill( mt2, weight );
    thisTree->fillTree( myTree, weight );

    
  } // for entries


  anaTree->finalize();


  delete tree;


  file->Close();
  delete file;
  

}



