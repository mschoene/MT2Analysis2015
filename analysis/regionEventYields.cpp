#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
//#include <cmath>


#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"


#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2EstimateSigSyst.h"
#include "interface/MT2EstimateAllSigSyst.h"
#include "interface/MT2DrawTools.h"
#include "interface/MT2Config.h"

#include "TRandom3.h"


#define mt2_cxx
#include "interface/mt2.h"





void randomizePoisson( MT2Analysis<MT2EstimateTree>* data );
template <class T>
MT2Analysis<T>* computeYield( const MT2Sample& sample, const MT2Config& cfg, std::string otherRegion="" );
template <class T>
MT2Analysis<T>* computeSigYield( const MT2Sample& sample, const MT2Config& cfg );
template <class T>
MT2Analysis<T>* mergeYields( std::vector< MT2Analysis<T> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max=-1, const std::string& legendName="" );
int matchPartonToJet( int index, MT2Tree* myTree );





int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|             Running regionEventYields              |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./regionEventYields [configFileName] [data/MC]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  bool onlyData = false;
  bool onlyMC   = false;
  bool onlySignal = false;
  if( argc > 2 ) {
    std::string dataMC(argv[2]);
    if( dataMC=="data" ) onlyData = true;
    else if( dataMC=="MC" || dataMC=="mc" ) onlyMC = true;
    else if( dataMC=="signal" ) onlySignal = true;
    else {
      std::cout << "-> You passed a second argument that isn't 'data', nor 'MC', nor 'signal', so I don't know what to do about it." << std::endl;
    }
  }

  std::string outputdir = cfg.getEventYieldDir();
  system(Form("mkdir -p %s", outputdir.c_str()));


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::vector<MT2Analysis<MT2EstimateTree>* > yields;
  //MT2Analysis<MT2EstimateTree>* dataYield;  

  if( cfg.useMC() && !onlyData && !onlySignal ) { // use MC BG estimates

    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 101, 999); // not interested in signal here (see later)
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(120);
    }


    
    std::vector< MT2Analysis<MT2EstimateTree>* > EventYield_incl;
    for( unsigned i=0; i<fSamples.size(); ++i ) {
      int this_id = fSamples[i].id;
      if( this_id<600 ) continue; // skip everything that is not ZJets
      if( this_id>=700 ) continue; //
      EventYield_incl.push_back( computeYield<MT2EstimateTree>( fSamples[i], cfg, "13TeV_2016_inclusive"  ));
    }
    MT2Analysis<MT2EstimateTree>* EventYield_zjets_inclusive = mergeYields<MT2EstimateTree>( EventYield_incl, "13TeV_2016_inclusive", "ZJets_inclusive", 600, 699, "Z+jets" );
    EventYield_zjets_inclusive->writeToFile(outputdir + "/ZJetsIncl.root");


    std::vector< MT2Analysis<MT2EstimateTree>* > EventYield;
    for( unsigned i=0; i<fSamples.size(); ++i ) {
      int this_id = fSamples[i].id;
      if( this_id>=200 && this_id<300 ) continue; // skip GJets
      if( this_id>=700 && this_id<800 ) continue; // skip DY
      EventYield.push_back( computeYield<MT2EstimateTree>( fSamples[i], cfg ));
    }



    std::cout << "-> Done looping on samples. Start merging." << std::endl;

    std::cout << "     merging Top..." << std::endl;
    MT2Analysis<MT2EstimateTree>* EventYield_top   = mergeYields<MT2EstimateTree>( EventYield, cfg.regionsSet(), "Top", 300, 499 ); // ttbar, single top, ttW, ttZ...
    std::cout << "     merging QCD..." << std::endl;
    MT2Analysis<MT2EstimateTree>* EventYield_qcd   = mergeYields<MT2EstimateTree>( EventYield, cfg.regionsSet(), "QCD", 100, 199 );
    std::cout << "     merging WJets..." << std::endl;
    MT2Analysis<MT2EstimateTree>* EventYield_wjets = mergeYields<MT2EstimateTree>( EventYield, cfg.regionsSet(), "WJets", 500, 599, "W+jets" );
    std::cout << "     merging ZJets..." << std::endl;
    MT2Analysis<MT2EstimateTree>* EventYield_zjets = mergeYields<MT2EstimateTree>( EventYield, cfg.regionsSet(), "ZJets", 600, 699, "Z+jets" );
    

    //    MT2Analysis<MT2EstimateTree>* EventYield_other = mergeYields<MT2EstimateTree>( EventYield, cfg.regionsSet(), "Other", 700, 999, "Other" );
    std::cout << "-> Done merging." << std::endl;

    yields.push_back( EventYield_qcd );
    yields.push_back( EventYield_wjets );
    yields.push_back( EventYield_zjets );
    yields.push_back( EventYield_top );
    //    yields.push_back( EventYield_other );



    if( cfg.dummyAnalysis() ) {
      MT2Analysis<MT2EstimateTree>* dataYield   = mergeYields<MT2EstimateTree>( EventYield, cfg.regionsSet(), "data", 100, 699 );
      yields.push_back( dataYield );
    } //else {
      //dataYield = new MT2Analysis<MT2EstimateTree>( "data", cfg.regionsSet(), 1 );
    //}

  } // if MC samples



  // load signal samples, if any
  std::vector< MT2Analysis< MT2EstimateAllSigSyst>* > signals;
  if( cfg.mcSamples()!="" && cfg.additionalStuff()!="noSignals" && !onlyData ) {

    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading signal samples from file: " << samplesFileName << std::endl;

    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 1000); // only signal (id>=1000)


    if( fSamples.size()==0 ) {

      std::cout << "No signal samples found, skipping." << std::endl;

    } else {
    
      for( unsigned i=0; i<fSamples.size(); ++i ) 
        signals.push_back( computeSigYield<MT2EstimateAllSigSyst>( fSamples[i], cfg ) );
    
//      std::cout << "     merging T1bbbb full scan..." << std::endl;
//      MT2Analysis<MT2EstimateAllSigSyst>* EventYield_T1bbbb   = mergeYields<MT2EstimateAllSigSyst>( signals, cfg.regionsSet(), "SMS_T1bbbb_fullScan", 1020, 1020 );
//      std::cout << "-> Done merging." << std::endl;
//      signals.push_back( EventYield_T1bbbb );

    } // if samples != 0

  } // if mc samples

  else if ( cfg.sigSamples()!="" && !onlyData ) {

    std::string samplesFileName = "../samples/samples_" + cfg.sigSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading signal samples from file: " << samplesFileName << std::endl;

    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 1000); // only signal (id>=1000)

    if( fSamples.size()==0 ) {

      std::cout << "No signal samples found, skipping." << std::endl;

    } else {

      for( unsigned i=0; i<fSamples.size(); ++i )
        signals.push_back( computeSigYield<MT2EstimateAllSigSyst>( fSamples[i], cfg ) );

//      std::cout << "     merging T1bbbb full scan..." << std::endl;
//      MT2Analysis<MT2EstimateAllSigSyst>* EventYield_T1bbbb   = mergeYields<MT2EstimateAllSigSyst>( signals, cfg.regionsSet(), "SMS_T1bbbb_fullScan", 1020, 1020 );
//      std::cout << "-> Done merging." << std::endl;
//      signals.push_back( EventYield_T1bbbb );

    } // if samples != 0
    
  } // if sig samples
  

  if( !(cfg.dummyAnalysis()) && cfg.dataSamples()!="" && !onlyMC  && !onlySignal ) {

    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    //    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "JetHTMHT"); //, 1, 99 );
    //    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, 1, 3 );
    // std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, -1, 0 );
    // std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "noDuplicates" );
    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "merged" );
    if( samples_data.size()==0 ) {
      std::cout << "There must be an error: samples_data is empty!" << std::endl;
      exit(1209);
    }

    std::vector< MT2Analysis<MT2EstimateTree>* > EventYield_data;
    for( unsigned i=0; i < samples_data.size(); ++i )
      EventYield_data.push_back( computeYield<MT2EstimateTree>( samples_data[i], cfg ) );

    MT2Analysis<MT2EstimateTree>* dataYield;
    //dataYield = EventYield_data[0];
    //dataYield->setName("data");
    //dataYield   = mergeYields<MT2EstimateTree>( EventYield_data, cfg.regionsSet(), "data", 1, 3 );
    dataYield   = mergeYields<MT2EstimateTree>( EventYield_data, cfg.regionsSet(), "data", -1, 10 );

    yields.push_back( dataYield );

  }


  if( yields.size()==0 && signals.size()==0 ) {
    std::cout << "-> Didn't end up with a single yield... something's wrong." << std::endl;
    exit(87);
  }


  // save MT2Analyses:
  if( yields.size()>0 ){
    yields[0]->writeToFile(outputdir + "/analyses.root");
  for( unsigned i=1; i<yields.size(); ++i )
    yields[i]->writeToFile(outputdir + "/analyses.root");
  for( unsigned i=0; i<signals.size(); ++i )
    signals[i]->writeToFile(outputdir + "/analyses.root");
  }
  else if( signals.size()>0 ){
    signals[0]->writeToFile(outputdir + "/analyses.root");
    for( unsigned i=1; i<signals.size(); ++i )
      signals[i]->writeToFile(outputdir + "/analyses.root");
  }
  cfg.saveAs(outputdir + "/config.txt");

  return 0;

}




template <class T>
MT2Analysis<T>* computeYield( const MT2Sample& sample, const MT2Config& cfg, std::string otherRegion ) {

  std::string regionsSet = cfg.regionsSet();
  if(otherRegion!="")
    regionsSet = otherRegion;

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
  

  MT2Tree myTree;
  if( cfg.additionalStuff()=="qgVars" || cfg.additionalStuff()=="hfContent" ) {
     myTree.loadGenStuff = true;
  } else {
    myTree.loadGenStuff = false;
  }
  myTree.Init(tree);



  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  MT2Analysis<T>* analysis = new MT2Analysis<T>( sample.sname, regionsSet, sample.id );
  
 
  if( sample.id < 1000){
    
    if( cfg.additionalStuff()=="qgVars" ) {
      
      T::addVar( analysis, "partId0" );
      T::addVar( analysis, "partId1" );
      T::addVar( analysis, "partId2" );
      T::addVar( analysis, "partId3" );
      T::addVar( analysis, "qgl0" );
      T::addVar( analysis, "qgl1" );
      T::addVar( analysis, "qgl2" );
      T::addVar( analysis, "qgl3" );
      T::addVar( analysis, "qglProd" );
      T::addVar( analysis, "qglAve" );
    }
    
    if( cfg.additionalStuff()=="hfContent" ) {
      T::addVar( analysis, "nTrueB" );
      T::addVar( analysis, "nTrueBJ" );
      T::addVar( analysis, "nTrueC" );
    }
    
    T::addVar( analysis, "jet1_pt" );
    T::addVar( analysis, "jet2_pt" );
    T::addVar( analysis, "mht" );
    
  }
  
  int nentries = tree->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);
    
    if( myTree.isData && !myTree.isGolden ) continue;

    if( regionsSet!="13TeV_noCut" )
      if( !myTree.passSelection(cfg.additionalStuff()) ) continue;

    if ( myTree.nJet30==1 && !myTree.passMonoJetId(0) ) continue;
      

    if( myTree.isData && !(myTree.run<=276811 || ( 278820<=myTree.run && myTree.run<=279931)) ) continue;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20;    
    float mt2  = (njets>1) ? myTree.mt2 : ht;
    //float mt2  = myTree.mt2;

    int GenSusyMScan1=0;
    int GenSusyMScan2=0;
    if(  myTree.evt_id > 999){
      GenSusyMScan1 = myTree.GenSusyMGluino;
      GenSusyMScan2 = myTree.GenSusyMNeutralino;
    }
   
    //Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi()*myTree.puWeight;
    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi();
    //Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi();
    Double_t weight_syst = 1.;

    if( !myTree.isData ){
      weight *= myTree.weight_btagsf;
      weight *= myTree.weight_lepsf;
    }

    if( myTree.evt_id > 1000 )
      weight_syst = myTree.weight_isr;
    
    //The filters to be applied to MC only
    if( !(myTree.nVert>0 && myTree.Flag_HBHENoiseFilter==1 && myTree.Flag_HBHENoiseIsoFilter==1 && myTree.Flag_EcalDeadCellTriggerPrimitiveFilter==1 && myTree.Flag_goodVertices==1 && myTree.Flag_eeBadScFilter==1 && myTree.Flag_badChargedHadronFilter==1)) continue;


    if( myTree.isData ) {
      
      if( !myTree.passFilters() ) continue;

    }

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

    if( myTree.met_miniaodPt/myTree.met_caloPt > 5.0 ) continue;

    if (myTree.isData) {

      if ( !(myTree.HLT_PFMET120_PFMHT120 || myTree.HLT_PFHT900 || myTree.HLT_PFHT300_PFMET110 || myTree.HLT_PFJet450) ) continue;
      //      if ( !(myTree.HLT_PFMET100_PFMHT100 || myTree.HLT_PFHT800 || myTree.HLT_PFHT300_PFMET100) ) continue;

    } // if is data

   
    T* thisEstimate = analysis->get( ht, njets, nbjets, minMTBmet, mt2 );
    //T* thisEstimate = analysis->get( ht, njets, nbjets, met, minMTBmet, mt2 );
    if( thisEstimate==0 ) continue;

    if( sample.id < 1000 ){

      thisEstimate->assignVar( "jet1_pt",  myTree.jet1_pt );
      thisEstimate->assignVar( "jet2_pt",  myTree.jet2_pt );
      thisEstimate->assignVar( "mht",  myTree.mht_pt );
      
      if( cfg.additionalStuff()=="qgVars" ) {
	
	// initialize
	thisEstimate->assignVar( "qgl0", -1. );
	thisEstimate->assignVar( "qgl1", -1. );
	thisEstimate->assignVar( "qgl2", -1. );
	thisEstimate->assignVar( "qgl3", -1. );
	thisEstimate->assignVar( "partId0", 0 );
	thisEstimate->assignVar( "partId1", 0 );
	thisEstimate->assignVar( "partId2", 0 );
	thisEstimate->assignVar( "partId3", 0 );
	
	float qglProd = 1.;
	float qglAve = 0.;
	int denom = 0;
	
	
	if( njets>0 && fabs(myTree.jet_eta[0])<2.5 ) {
	  
	  float qgl0 = myTree.jet_qgl[0];
	  thisEstimate->assignVar( "qgl0", qgl0 );
	  qglProd *= qgl0;
	  qglAve += qgl0;
	  denom++;
	  thisEstimate->assignVar( "partId0", matchPartonToJet( 0, &myTree ) );
	  //thisEstimate->assignVar( "partId0", myTree.jet_mcFlavour[0] );
	  
	}
	
	
	if( njets>1 && fabs(myTree.jet_eta[1])<2.5 ) {
	  
	  float qgl1 = myTree.jet_qgl[1];
	  thisEstimate->assignVar( "qgl1", qgl1 );
	  qglProd *= qgl1;
	  qglAve += qgl1;
	  denom++;
	  
	  thisEstimate->assignVar( "partId1", matchPartonToJet( 1, &myTree ) );
	  //thisEstimate->assignVar( "partId1", myTree.jet_mcFlavour[1] );
	  
	}
        
	if( njets>2 && fabs(myTree.jet_eta[2])<2.5 ) {
	  
	  float qgl2 = myTree.jet_qgl[2];
	  thisEstimate->assignVar( "qgl2", qgl2 );
	  qglProd *= qgl2;
	  qglAve += qgl2;
	  denom++;
	  
	  thisEstimate->assignVar( "partId2", matchPartonToJet( 2, &myTree ) );
	  //thisEstimate->assignVar( "partId2", myTree.jet_mcFlavour[2] );
	  
	}
        
	
	if( njets>3 && fabs(myTree.jet_eta[3])<2.5 ) {
	  
	  float qgl3 = myTree.jet_qgl[3];
	  thisEstimate->assignVar( "qgl3", qgl3 );
	  qglProd *= qgl3;
	  qglAve += qgl3;
	  denom++;
	  
	  thisEstimate->assignVar( "partId3", matchPartonToJet( 3, &myTree ) );
	  //thisEstimate->assignVar( "partId3", myTree.jet_mcFlavour[3] );
	  
	}
	
	qglAve /= (float)denom;
        
	thisEstimate->assignVar( "qglProd", qglProd );
	thisEstimate->assignVar( "qglAve", qglAve );
	
	
      } 
      
      if( cfg.additionalStuff()=="hfContent" ) {
	
	float nTrueB=0.;
	float nTrueC=0.;
	
	for( int ipart=0; ipart<myTree.ngenPart; ++ipart ) {
	  
	  if( myTree.genPart_pt[ipart] < 20. ) continue;
	  if( abs(myTree.genPart_eta[ipart])>2.5 ) continue;
	  if( myTree.genPart_status[ipart] != 23 ) continue;
	  
	  if( abs(myTree.genPart_pdgId[ipart])==5 )
	    nTrueB+=1.;
	  if( abs(myTree.genPart_pdgId[ipart])==4 )
	    nTrueC+=1.;
	  
	}
	
	
	float nTrueBJ=0.;
	
	for( int ijet=0; ijet<myTree.njet; ++ijet ) {
	  
	  if( myTree.jet_pt[ijet] <20. ) continue;
	  if( abs(myTree.jet_eta[ijet])>2.5 ) continue;
	  
	  if( abs(myTree.jet_mcFlavour[ijet])==5 )
	    nTrueBJ+=1.;
	  
	}
	
	thisEstimate->assignVar( "nTrueB", nTrueB );
	thisEstimate->assignVar( "nTrueC", nTrueC );
	thisEstimate->assignVar( "nTrueBJ", nTrueBJ );
	
      }
    }
    
    if( sample.id < 1000 ){
     
      thisEstimate->assignTree( myTree, weight );
      thisEstimate->tree->Fill();
    
    }
    
    thisEstimate->yield->Fill( mt2, weight );
    thisEstimate->yield3d->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight );
    
  } // for entries
    
  //ofs.close();

  analysis->finalize();
  
  delete tree;

  file->Close();
  delete file;
  
  return analysis;

}

template <class T>
MT2Analysis<T>* computeSigYield( const MT2Sample& sample, const MT2Config& cfg ) {

  bool dogenmet = false;

  TString sigSampleName(sample.name);
  TFile* sigXSFile;
  if(sigSampleName.Contains("T1qqqq") || sigSampleName.Contains("T1bbbb") || sigSampleName.Contains("T1tttt"))
    sigXSFile = TFile::Open("/shome/casal/SUSxsecs/SUSYCrossSections13TeVgluglu.root");
  else if(sigSampleName.Contains("T2bb") || sigSampleName.Contains("T2tt"))
    sigXSFile = TFile::Open("/shome/casal/SUSxsecs/SUSYCrossSections13TeVstopstop.root");
  else
    sigXSFile = TFile::Open("/shome/casal/SUSxsecs/SUSYCrossSections13TeVsquarkantisquark.root");

  TH1F* sigXS = (TH1F*) sigXSFile->Get("xs");

  if(sigSampleName.Contains("T2qq"))
    sigXS->Scale(8./10);

  std::string regionsSet = cfg.regionsSet();
  
  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample (computeSigYield): " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
  

  MT2Tree myTree;
  myTree.Init(tree);



  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  MT2Analysis<T>* analysis = new MT2Analysis<T>( sample.sname, regionsSet, sample.id );
  
 
  
  int nentries = tree->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);
    
    bool passGenMET=false;
    if(dogenmet) passGenMET =true;
    
    bool passRecoMET=true;

//    if( regionsSet!="13TeV_noCut" )
//      if( !myTree.passSelection(cfg.additionalStuff()) ) continue;

    if( regionsSet!="13TeV_noCut" ){
     
      if( !myTree.passSelection(cfg.additionalStuff()) ) passRecoMET=false;
      
      if(dogenmet)
	if( !myTree.passSelection("genmet") ) passGenMET=false;      
      
      if (!passGenMET && !passRecoMET) continue;

    }
      

    if ( myTree.nJet30==1 && !myTree.passMonoJetId(0) ) continue;
    if ( myTree.nJet20BadFastsim > 0 ) continue;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float met_gen  = myTree.met_genPt;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20;    
    float mt2  = (njets>1) ? myTree.mt2 : ht;
    float mt2_genmet;

    if(dogenmet)
      mt2_genmet = (njets>1) ? myTree.mt2_genmet : ht;

    float weight_isr = myTree.weight_isr;
    float isr_UP = myTree.weight_isr_UP;
    float isr_DN = myTree.weight_isr_DN;

    float weight_lepsf = myTree.weight_lepsf;
    float lepsf_UP = myTree.weight_lepsf_UP;
    float lepsf_DN = myTree.weight_lepsf_DN;

    float weight_btagsf = myTree.weight_btagsf;
    float btag_heavy_UP = myTree.weight_btagsf_heavy_UP;
    float btag_heavy_DN = myTree.weight_btagsf_heavy_DN;
    float btag_light_UP = myTree.weight_btagsf_light_UP;
    float btag_light_DN = myTree.weight_btagsf_light_DN;
    


    int GenSusyMScan1=0;
    int GenSusyMScan2=0;
    if(  myTree.evt_id > 999){
     
      if(sigSampleName.Contains("T2qq")){
	
	GenSusyMScan1 = myTree.GenSusyMSquark;
	GenSusyMScan2 = myTree.GenSusyMNeutralino;

      }
      else if(sigSampleName.Contains("T2bb")){
	
	GenSusyMScan1 = myTree.GenSusyMSbottom;
	GenSusyMScan2 = myTree.GenSusyMNeutralino;
	
      }
      else if(sigSampleName.Contains("T2tt")){
	
	GenSusyMScan1 = myTree.GenSusyMStop;
	GenSusyMScan2 = myTree.GenSusyMNeutralino;

      }
      else{
      
      GenSusyMScan1 = myTree.GenSusyMGluino;
      GenSusyMScan2 = myTree.GenSusyMNeutralino;
      
//      GenSusyMScan1 = myTree.GenSusyMScan1;
//      GenSusyMScan2 = myTree.GenSusyMScan2;
      
      }

    }
 
    //Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi()*myTree.puWeight;
    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi();
    //Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi(); //Keeping normalization to luminosity for signal?
    Double_t weight_syst = 1.;

    //weight = 1000./nentries*cfg.lumi(); //Exceptionally for signal from muricans 

    if( !myTree.isData ){
      weight *= weight_btagsf;
      weight *= weight_lepsf;
      weight *= weight_isr;
    }


//    if( myTree.evt_id > 1000 )
//      weight_syst = myTree.weight_isr;
    

    float sig_xs=0.;
    if( myTree.evt_id >= 1000  && myTree.evt_id < 2000){
   
      int thisBinX = sigXS->FindBin( GenSusyMScan1 );
     
      sig_xs = sigXS->GetBinContent(thisBinX);
     
      weight *= sig_xs;
   
    }
    
    
    T* thisEstimate = analysis->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisEstimate==0 ) continue;
    
    if(passRecoMET){
      
      thisEstimate->yield->Fill( mt2, weight );
      thisEstimate->yield3d->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight );
      
      //    thisEstimate->yield3d_systUp->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight*(1.+(weight_syst-1.)));
      //    thisEstimate->yield3d_systDown->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight*(1.-(weight_syst-1.)));
      
//      thisEstimate->yield3d_isr_UP->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight*(1.+(isr-1.)));
//      thisEstimate->yield3d_isr_DN->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight*(1.-(isr-1.)));
//      
//      thisEstimate->yield3d_btag_heavy_UP->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight*(1.+(btag_heavy_UP-1.)));
//      thisEstimate->yield3d_btag_heavy_DN->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight*(1.-(btag_heavy_DN-1.)));
//      
//      thisEstimate->yield3d_btag_light_UP->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight*(1.+(btag_light_UP-1.)));
//      thisEstimate->yield3d_btag_light_DN->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight*(1.-(btag_light_DN-1.)));


      thisEstimate->yield3d_isr_UP->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight/weight_isr*(isr_UP) );
      thisEstimate->yield3d_isr_DN->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight/weight_isr*(isr_DN) );

      thisEstimate->yield3d_lepsf_UP->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight/weight_lepsf*(lepsf_UP) );
      thisEstimate->yield3d_lepsf_DN->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight/weight_lepsf*(lepsf_DN) );

      thisEstimate->yield3d_btag_heavy_UP->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight/weight_btagsf*(btag_heavy_UP) );
      thisEstimate->yield3d_btag_heavy_DN->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight/weight_btagsf*(btag_heavy_DN) );
            
      thisEstimate->yield3d_btag_light_UP->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight/weight_btagsf*(btag_light_UP) );
      thisEstimate->yield3d_btag_light_DN->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight/weight_btagsf*(btag_light_DN) );
                  
    }
    
    if(dogenmet && passGenMET){
      
      thisEstimate->yield3d_genmet->Fill( mt2_genmet, GenSusyMScan1, GenSusyMScan2, weight );

    }
    
  } // for entries
    
  //ofs.close();

  analysis->finalize();
  
  delete tree;

  file->Close();
  delete file;
  
  return analysis;

}


template <class T>
MT2Analysis<T>* mergeYields( std::vector<MT2Analysis<T> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max, const std::string& legendName ) {

  if( id_max<0 ) id_max=id_min;

  MT2Analysis<T>* return_EventYield = new MT2Analysis<T>(name, regionsSet, id_min, legendName);

  for( unsigned i=0; i<EventYield.size(); ++i ) {

    if( EventYield[i]->getId() >= id_min && EventYield[i]->getId() <= id_max ) {

       *(return_EventYield) += *(EventYield[i]);

    }

  } // for EventYield


  return return_EventYield;

}





void randomizePoisson( MT2Analysis<MT2EstimateTree>* data ) {

  TRandom3 rand(13);


//  std::set<MT2HTRegion> HTRegions = data->getHTRegions();
//  std::set<MT2SignalRegion> signalRegions = data->getSignalRegions();

  std::set<MT2Region> MT2Regions = data->getRegions();

//  for( std::set<MT2HTRegion>::iterator iHT = HTRegions.begin(); iHT!=HTRegions.end(); ++iHT ) {
//    for( std::set<MT2SignalRegion>::iterator iSR = signalRegions.begin(); iSR!=signalRegions.end(); ++iSR ) {

      for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

	MT2Region thisRegion( (*iMT2) );
	
	TH1D* h1_data = data->get(thisRegion)->yield;

	for( int ibin=1; ibin<h1_data->GetXaxis()->GetNbins()+1; ++ibin ) {

	  int poisson_data = rand.Poisson(h1_data->GetBinContent(ibin, 0, 0));
	  h1_data->SetBinContent(ibin, 0, 0, poisson_data);
	  h1_data->SetBinError(ibin, 0, 0,  0.);
	  
	}  // for bins

      }// for MT2 regions

//    }// for signal regions
//  }// for HT regions

}


int matchPartonToJet( int index, MT2Tree* myTree ) {

  float deltaRMin = 0.3; // at least 0.3
  int foundId = 0;

  TLorentzVector jet;
  jet.SetPtEtaPhiM( myTree->jet_pt[index], myTree->jet_eta[index], myTree->jet_phi[index], myTree->jet_mass[index] );


  for( int i=0; i<myTree->ngenPart; ++i ) {

    if( myTree->genPart_status[i]!=23 ) continue;
    if( !(myTree->genPart_pdgId[i]==21 || abs(myTree->genPart_pdgId[i]<6)) ) continue;

    TLorentzVector thisPart;
    thisPart.SetPtEtaPhiM( myTree->genPart_pt[i], myTree->genPart_eta[i], myTree->genPart_phi[i], myTree->genPart_mass[i] );

    float thisDeltaR = jet.DeltaR(thisPart);
    if( thisDeltaR<deltaRMin ) {
      deltaRMin = thisDeltaR;
      foundId = myTree->genPart_pdgId[i];
    }
  }

  return foundId;

}
