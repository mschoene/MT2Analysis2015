#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateQCD.h"

#define mt2_cxx
#include "../interface/mt2.h"

#include "TEventList.h"
#include "TVector2.h"

// this analysis:
// fills all necesary histograms for QCD estimation

float lumi = 0.042; //fb-1

// to do: make these as optional arguments
bool doDiffMetMht = true;  // apply the |MET-MHT|/MET cut?
bool doHFjetVeto  = true;  // apply HF jet veto

MT2Analysis<MT2EstimateQCD>* fillCR( const MT2Sample& sample, const std::string& regionsSet);
MT2Analysis<MT2EstimateQCD>* merge ( std::vector<MT2Analysis<MT2EstimateQCD> *> anas, const std::string& regionsSet, const std::string& name, int id_min, int id_max );

int main( int argc, char* argv[] ) {


  if( argc>2 ) {
    std::cout << "USAGE: ./qcdControlRegion [samplesFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  std::string samplesFileName = "74X_jecV4_MET30_QCD";
  if( argc>1 ) {
    std::string samplesFileName_tmp(argv[1]); 
    samplesFileName = samplesFileName_tmp;
  }


  std::string regionsSet = "zurich_HTtriggers";
  //std::string regionsSet = "zurich_HTtriggers2";

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string outputdir( Form("QCDcontrolRegion_%s_%s_%.3ffb%s%s", samplesFileName.c_str(), regionsSet.c_str(), lumi, (doDiffMetMht ? "" : "_noDiffMetMht"), (!doHFjetVeto ? "" : "_HFjetVeto") ) );
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

  controlRegions[0]->finalize();
  controlRegions[0]->writeToFile( outputdir + "/qcdCR.root" );
  for( unsigned i=1; i<samples.size(); ++i ) {
    controlRegions[i]->finalize();
  controlRegions[i]->addToFile( outputdir + "/qcdCR.root" );
  }
  qcdCR->addToFile( outputdir + "/qcdCR.root" );
  
  //qcdCR->writeToFile( outputdir + "/qcdCR.root" );

  MT2Analysis<MT2EstimateQCD>* topCR    = merge( controlRegions, regionsSet, "topCR"  , 300, 499 );
  MT2Analysis<MT2EstimateQCD>* wjetsCR  = merge( controlRegions, regionsSet, "wjetsCR", 500, 599 );
  MT2Analysis<MT2EstimateQCD>* zjetsCR  = merge( controlRegions, regionsSet, "zjetsCR", 600, 699 );
  topCR  ->finalize();
  wjetsCR->finalize();
  zjetsCR->finalize();

  MT2Analysis<MT2EstimateQCD>* dataCR  = merge( controlRegions, regionsSet, "dataCR", 1, 99 );
  dataCR->finalize();

  topCR  ->addToFile( outputdir + "/qcdCR.root" );
  wjetsCR->addToFile( outputdir + "/qcdCR.root" );
  zjetsCR->addToFile( outputdir + "/qcdCR.root" );
  dataCR ->addToFile( outputdir + "/qcdCR.root" );
}


MT2Analysis<MT2EstimateQCD>* merge( std::vector<MT2Analysis<MT2EstimateQCD> *> anas, const std::string& regionsSet, const std::string& name, int id_min, int id_max ) {

  if( id_max<0 ) id_max=id_min;

  MT2Analysis<MT2EstimateQCD>* ana = new MT2Analysis<MT2EstimateQCD>(name, regionsSet);

  for( unsigned i=0; i<anas.size(); ++i ) {
    if( anas[i]->getId() >= id_min && anas[i]->getId() <= id_max )
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
  
  // In absence of skimmed ntuples let's filter to gain some speed
  TString filter = "mt2>30&&ht>1000"; // ht>1000 for 50ns only! dangerous temporary hardcoded cut to be revisited
  tree->Draw(">>selList", filter);
  TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
  tree->SetEventList(myEvtList);

  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  MT2Analysis<MT2EstimateQCD>* analysis = new MT2Analysis<MT2EstimateQCD>( sample.sname, regionsSet, sample.id );

  // int nentries = tree->GetEntries();
  // for( int iEntry=0; iEntry<nentries; ++iEntry ) {

  //   if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
  int nentries = myEvtList->GetN();
  for( int jEntry=0; jEntry<nentries; ++jEntry ) {
    int iEntry = myEvtList->GetEntry(jEntry);

    if( jEntry % 50000 == 0 ) std::cout << "    Entry: " << jEntry << " / " << nentries << std::endl;
    myTree.GetEntry(iEntry);


    if ( !(myTree.met_pt>30 && myTree.nVert>0 && myTree.nJet30>=2)   ) continue;
    if ( !(myTree.passLeptonVeto() && myTree.passIsoTrackVeto())     ) continue;
    if ( doDiffMetMht && !(myTree.diffMetMht < 0.5*myTree.met_pt)    ) continue;
    if ( myTree.isData==1 && !(myTree.Flag_CSCTightHaloFilter==1 && 
			       myTree.Flag_HBHENoiseFilter   ==1 &&
			       myTree.Flag_eeBadScFilter     ==1   ) )  continue;


    if (doHFjetVeto) {
      int nJetHF30 = 0;
      for(int j=0; j<myTree.njet; ++j){
	if( myTree.jet_pt[j] < 30. || fabs(myTree.jet_eta[j]) < 3.0 ) continue;
	else ++nJetHF30;
      }
      if ( !(nJetHF30==0))  continue;
    }

    float ht   = myTree.ht;
    float mt2  = myTree.mt2;
    float dphi = myTree.deltaPhiMin;
    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20;
    int isData = myTree.isData;

    Double_t weight = isData ? 1.0 : myTree.evt_scale1fb*lumi;

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
