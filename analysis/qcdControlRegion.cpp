#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"





bool monojet=false;


int round(float d) {
  return (int)(floor(d + 0.5));
}


void computeYield( const MT2Sample& sample, const MT2Config& cfg,  MT2Analysis<MT2EstimateTree>* anaTree );
float getAverageISRWeight(const int id, const int var);



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
    std::cout << "USAGE: ./qcdControlRegion [configFileName] [data/MC/all] [monojet=false]" << std::endl;
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
  }

  if( onlyData ) {
    std::cout << "-> Will run only on data." << std::endl;
  } else if( onlyMC ) {
    std::cout << "-> Will run only on MC." << std::endl;
  } else {
    std::cout << "-> Will run on both data and MC." << std::endl;
  }
 


  if( argc > 3 ) {
    std::string monojetstr(argv[3]);
    if( monojetstr=="true" || monojetstr=="monojet" ) {
      monojet=true;
      std::cout << "-> Monojet analysis (so releasing the MT2>50. cut)" << std::endl;
    }
  }



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this 

  std::string outputdir = cfg.getEventYieldDir() + "/qcdControlRegion";
  system(Form("mkdir -p %s", outputdir.c_str()));


  if( cfg.useMC() && !onlyData ) { // run on MC

    std::string samplesFile = "../samples/samples_" + cfg.qcdMCSamples() + ".dat";
    if(monojet) samplesFile = "../samples/samples_" + cfg.qcdMonoJetMCSamples() + ".dat";

    std::vector<MT2Sample> samples_zinv = MT2Sample::loadSamples(samplesFile, 602, 699);
    std::vector<MT2Sample> samples_wjet = MT2Sample::loadSamples(samplesFile, 502, 599);
    //std::vector<MT2Sample> samples_top  = MT2Sample::loadSamples(samplesFile, 300, 399); // ignore single top and rares: faster
    std::vector<MT2Sample> samples_top  = MT2Sample::loadSamples(samplesFile, 300, 499);
    std::vector<MT2Sample> samples_qcd  = MT2Sample::loadSamples(samplesFile, 100, 199);


    MT2Analysis<MT2EstimateTree>* qcdCRtree = new MT2Analysis<MT2EstimateTree>( "qcdCRtree", "13TeV_2016_inclusive" );
    MT2EstimateTree::addVar( qcdCRtree, "jet1_pt" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_pt" );
    
    
    for( unsigned i=0; i<samples_zinv.size(); ++i ) 
      computeYield( samples_zinv[i], cfg, qcdCRtree );
    for( unsigned i=0; i<samples_wjet.size(); ++i ) 
      computeYield( samples_wjet[i], cfg, qcdCRtree );
    for( unsigned i=0; i<samples_top.size(); ++i ) 
      computeYield( samples_top[i], cfg, qcdCRtree );
    for( unsigned i=0; i<samples_qcd.size(); ++i ) 
      computeYield( samples_qcd[i], cfg, qcdCRtree );
    

   
    std::string mcFile = outputdir + "/mc";
    if( monojet ) mcFile = mcFile + "_forMonojet";
    mcFile = mcFile + ".root";
    qcdCRtree->writeToFile( mcFile, "RECREATE" );

  }

  if( !(cfg.dummyAnalysis()) && cfg.qcdDataSamples()!="" && !onlyMC  ) {

    std::string samplesFile_data = "../samples/samples_" + cfg.qcdDataSamples() + ".dat";
    if(monojet) samplesFile_data = "../samples/samples_" + cfg.qcdMonoJetDataSamples() + ".dat";

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, 1, 3 );
    if( samples_data.size()==0 ) {
      std::cout << "There must be an error: samples_data is empty!" << std::endl;
      exit(1209);
    }

    MT2Analysis<MT2EstimateTree>* data = new MT2Analysis<MT2EstimateTree>( "qcdCRtree", "13TeV_2016_inclusive" );
    MT2EstimateTree::addVar( data, "jet1_pt" );
    MT2EstimateTree::addVar( data, "jet2_pt" );
    //MT2Analysis<MT2EstimateTree>* data = new MT2Analysis<MT2EstimateTree>( "qcdCRtree", cfg.regionsSet() );
    for( unsigned i=0; i < samples_data.size(); ++i )
      computeYield( samples_data[i], cfg, data );

    std::string dataFile = outputdir + "/data";
    if( monojet ) dataFile = dataFile + "_forMonojet";
    dataFile = dataFile + ".root";
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
      if (  myTree.isGolden == 0 ) continue;
      if ( !myTree.passFilters() ) continue;
    }
    else {
      //if ( !(myTree.Flag_badMuonFilter>0 && myTree.Flag_badChargedHadronFilter>0 && myTree.Flag_EcalDeadCellTriggerPrimitiveFilter>0) ) continue;
      if ( !(myTree.Flag_badChargedHadronFilter>0 && myTree.Flag_EcalDeadCellTriggerPrimitiveFilter>0) ) continue; // americans don't have muon filter
    }
    if (myTree.met_miniaodPt/myTree.met_caloPt > 5.0) continue; // RA2 filter for QCD MC
    

    if( !myTree.passSelection("qcd") ) continue;


    float minMTBmet = myTree.minMTBMet;
    //float met       = myTree.met_pt;
    int njets       = myTree.nJet30;
    int nbjets      = myTree.nBJet20;    
    float mt2       = (njets>1) ? myTree.mt2 : myTree.jet1_pt;
    float ht        = myTree.ht;


   //crazy events! To be piped into a separate txt file
    if(myTree.jet_pt[0] > 13000){
      std::cout << "Rejecting weird event at run:lumi:evt = " << myTree.run << ":" << myTree.lumi << ":" << myTree.evt << std::endl;
      continue;
    }
    //NEW check if is there is a nan
    if( isnan(myTree.ht) || isnan(myTree.met_pt) ||  isinf(myTree.ht) || isinf(myTree.met_pt)  ){
      std::cout << "Rejecting nan/inf event at run:lumi:evt = " << myTree.run << ":" << myTree.lumi << ":" << myTree.evt << std::endl;
      continue;
    }


    if (myTree.isData) {

      if( njets<2 ) continue; // qcd CRs never use 1jet events, not even for monojet CR
      // don't bother about dataset of origin (OR of datasets)
      // fill id with trigger bit logic (0b01 HT-only triggered, 0b10 signal triggered)
      // signal triggered events include the MET>250 selection for HT<1000
      myTree.evt_id = 0;
      if( (ht>1000.&&myTree.HLT_PFHT900) || (ht<1000.&&ht>575.&&myTree.HLT_PFHT475_Prescale) || (ht<575.&&ht>450.&&myTree.HLT_PFHT350_Prescale) || (ht<450.&&ht>250.&&myTree.HLT_PFHT125_Prescale) )
	myTree.evt_id += 0b01;      // HT-only triggered -> turn bit 1 on
      if(ht>250 && (myTree.HLT_PFHT900||myTree.HLT_PFJet450||myTree.HLT_PFHT300_PFMET110||myTree.HLT_PFMET120_PFMHT120) && (ht>1000 || myTree.met_pt>250) )
	myTree.evt_id += 0b10;      // passes signal triggers -> turn bit 2 on
      if( myTree.evt_id==0 )
	continue;                   // doesn't pass any trigger

    } // if is data


    if( monojet ) {
      // only store signal triggered events for monojet CR
      if( !( (njets==2 && myTree.deltaPhiMin<0.3 && myTree.jet1_pt>200. && myTree.met_pt>200. && (myTree.evt_id&0b10)==0b10) ) ) continue;
      if( !myTree.passMonoJetId(0) ) continue;

      //Fix for monojetCR for shape comparsion Data/MC
      if( !myTree.Flag_EcalDeadCellTriggerPrimitiveFilter ) continue;

      if( !myTree.isData && (myTree.met_pt/ myTree.met_caloPt > 5.0) ) continue;


    } else {
      if( mt2<50. || njets<2 ) continue; // remove unnecesary for lighter trees
    }

    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*myTree.weight_btagsf;
    if (myTree.evt_id>300 && myTree.evt_id<310)
      weight *= myTree.weight_isr/getAverageISRWeight(myTree.evt_id,0);

    //float myht = njets==1 ? 201. : ht; // let everything (mt2=ht>40) pass for monojet
    MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisTree==0 ) continue;


    thisTree->yield->Fill( mt2, weight );
    thisTree->assignVar( "jet1_pt", myTree.jet1_pt );
    thisTree->assignVar( "jet2_pt", myTree.jet2_pt );
    thisTree->fillTree( myTree, weight );

    
  } // for entries


  anaTree->finalize();


  delete tree;


  file->Close();
  delete file;
  

}


float getAverageISRWeight(const int evt_id, const int var){

  // madgraph ttsl, from RunIISpring16MiniAODv2
  if (evt_id == 301 || evt_id == 302) {
    if (var == 0) return 0.910; // nominal
    else if (var == 1) return 0.955; // UP
    else if (var == -1) return 0.865; // DN
  }
  // madgraph ttdl, from RunIISpring16MiniAODv2
  else if (evt_id == 303) {
    if (var == 0) return 0.897; // nominal
    else if (var == 1) return 0.948; // UP
    else if (var == -1) return 0.845; // DN
  }

  std::cout << "WARNING: MT2Looper::getAverageISRWeight: didn't recognize either evt_id: " << evt_id
	    << " or variation: " << var << std::endl;
  return 1.;

}

