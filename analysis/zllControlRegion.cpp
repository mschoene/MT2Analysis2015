#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "../interface/MT2Analysis.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Sample.h"
#include "../interface/MT2DrawTools.h"
#include "../interface/MT2Config.h"

#define mt2_cxx
#include "../interface/mt2.h"

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"


int round(float d) {
  return (int)(floor(d + 0.5));
}



bool do_bg = true;



void computeYieldSnO( const MT2Sample& sample, const MT2Config& cfg,   
		      MT2Analysis<MT2EstimateTree>* anaTree,  
		      MT2Analysis<MT2EstimateTree>* anaTree_of );
void addVariables(MT2Analysis<MT2EstimateTree>* anaTree);
void roundLikeData( MT2Analysis<MT2EstimateTree>* data );


int main(int argc, char* argv[]) {


  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|            Running zllControlRegion                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./zllControlRegion [configFileName] [data/MC]" << std::endl;
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
    else if( dataMC=="MC" ) onlyMC = true;
    else {
      std::cout << "-> You passed a second argument that isn't 'data' nor 'MC', so I don't know what to do about it." << std::endl;
    }
  }


  TH1::AddDirectory(kFALSE); //stupid ROOT memory allocation needs this

  std::string outputdir = cfg.getEventYieldDir() + "/zllControlRegion";
  system(Form("mkdir -p %s", outputdir.c_str()));

   std::string regionsSet;// = "13TeV_inclusive";
  regionsSet=cfg.regionsSet();
  std::cout << "-> Using regions: " << regionsSet << std::endl;


  if( cfg.useMC() && !onlyData ) { // run on MC  
  
    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;


    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, "DYJetsToLL", 699, 800); // not interested in signal here
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(1209);
    }

    MT2Analysis<MT2EstimateTree>* mcTree = new MT2Analysis<MT2EstimateTree>( "zllCR", cfg.regionsSet() );
    addVariables(mcTree); //Adds some additional variables Zpt,Zmass, raw MT2...


    MT2Analysis<MT2EstimateTree>* mcTree_of = new MT2Analysis<MT2EstimateTree>( "zllCR_of", cfg.regionsSet() );
    addVariables(mcTree_of);

  

    // std::vector< MT2Analysis<MT2EstimateTree>* > EventYield;
    for( unsigned i=0; i<fSamples.size(); ++i ) 
      computeYieldSnO( fSamples[i], cfg, mcTree, mcTree_of);
    //      EventYield.push_back( computeYield( fSamples[i], cfg, cfg.lumi() ) );
    
    // MT2Analysis<MT2EstimateTree>* zllCR = mergeYields( EventYield, cfg.regionsSet(), "DYJets", 700, 799, "DYJets" );
    // zllCR->setName("zllCR");
    // zllCR->writeToFile(outputdir+"/mc.root");
    mcTree->writeToFile(outputdir+"/mc.root");
    mcTree_of->writeToFile(outputdir+"/mc_of.root");
    
    if( cfg.dummyAnalysis() ) {
    
      roundLikeData(mcTree); 
      mcTree->addToFile(outputdir+"/data.root");
      roundLikeData(mcTree_of); 
      mcTree_of->addToFile(outputdir+"/data_of.root");
     
    }
    
    if(do_bg==true){
      //MC
      MT2Analysis<MT2EstimateTree>* mc_top = new MT2Analysis<MT2EstimateTree>( "Top", cfg.regionsSet(),300, "Top" );
      MT2Analysis<MT2EstimateTree>* mc_top_of = new MT2Analysis<MT2EstimateTree>( "Top", cfg.regionsSet(),300, "Top" );
      addVariables(mc_top);      addVariables(mc_top_of);
      std::vector<MT2Sample> fSamples_top = MT2Sample::loadSamples(samplesFileName, 300, 499);   
      for( unsigned i=0; i<fSamples_top.size(); ++i )
	computeYieldSnO( fSamples_top[i], cfg, mc_top, mc_top_of);
   
      MT2Analysis<MT2EstimateTree>* mc_qcd = new MT2Analysis<MT2EstimateTree>( "QCD", cfg.regionsSet(),100, "QCD" );
      MT2Analysis<MT2EstimateTree>* mc_qcd_of = new MT2Analysis<MT2EstimateTree>( "QCD", cfg.regionsSet(),100, "QCD");
      addVariables(mc_qcd);      addVariables(mc_qcd_of);
      std::vector<MT2Sample> fSamples_qcd = MT2Sample::loadSamples(samplesFileName, 100, 199);   
      for( unsigned i=0; i<fSamples_qcd.size(); ++i )
	computeYieldSnO( fSamples_qcd[i], cfg, mc_qcd, mc_qcd_of);

      MT2Analysis<MT2EstimateTree>* mc_wjets = new MT2Analysis<MT2EstimateTree>( "WJets", cfg.regionsSet(),500, "W+jets"  );
      MT2Analysis<MT2EstimateTree>* mc_wjets_of = new MT2Analysis<MT2EstimateTree>( "WJets", cfg.regionsSet(),500, "W+jets");
      addVariables(mc_wjets);      addVariables(mc_wjets_of);
      std::vector<MT2Sample> fSamples_wjets = MT2Sample::loadSamples(samplesFileName, 500, 599);   
      for( unsigned i=0; i<fSamples_wjets.size(); ++i )
	computeYieldSnO( fSamples_wjets[i], cfg, mc_wjets, mc_wjets_of);

      MT2Analysis<MT2EstimateTree>* mc_zjets = new MT2Analysis<MT2EstimateTree>( "ZJets", cfg.regionsSet(),600, "Z+jets" );
      MT2Analysis<MT2EstimateTree>* mc_zjets_of = new MT2Analysis<MT2EstimateTree>( "ZJets", cfg.regionsSet(),600, "Z+jets");
      addVariables(mc_zjets);      addVariables(mc_zjets_of);
      std::vector<MT2Sample> fSamples_zjets = MT2Sample::loadSamples(samplesFileName, 600, 699);   
      for( unsigned i=0; i<fSamples_zjets.size(); ++i )
	computeYieldSnO( fSamples_zjets[i], cfg, mc_zjets, mc_zjets_of); 

 
  
      std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
      MT2Analysis<MT2EstimateTree>* mc_zll   = mcTree;
      mc_zll->setName("DYJets");
      mc_zll->setFullName("DY+jets");
      
      bgYields.push_back( mc_qcd );
      bgYields.push_back( mc_wjets );
      bgYields.push_back( mc_zjets );
      bgYields.push_back( mc_top );
      
      std::string outFile = outputdir + "/ZllPurityTrees.root";
      mc_zll->writeToFile( outFile );
      mc_top->addToFile( outFile );
      mc_qcd->addToFile( outFile );
      mc_wjets->addToFile( outFile );
      mc_zjets->addToFile( outFile );


      //For the OPPOSITE FLAVOR EVENTS:
      std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_of;   
      MT2Analysis<MT2EstimateTree>* mc_zll_of   = mcTree_of;
      mc_zll_of->setName("DYJets");
      mc_zll_of->setFullName("DY+jets");
      
      bgYields_of.push_back( mc_qcd_of );
      bgYields_of.push_back( mc_wjets_of );
      bgYields_of.push_back( mc_zjets_of );
      bgYields_of.push_back( mc_top_of );
      
      std::string outFile_of = outputdir + "/ZllPurityTrees_of.root";
      mc_zll_of->writeToFile( outFile_of );
      mc_top_of->addToFile( outFile_of );
      mc_qcd_of->addToFile( outFile_of );
      mc_wjets_of->addToFile( outFile_of );
      mc_zjets_of->addToFile( outFile_of );
     
    } //End do background trees
    
  } //if only MC
  
  if( !onlyMC ) {
    
    //DATA
    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;
    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "Double");  
    std::vector<MT2Sample> samples_data_of = MT2Sample::loadSamples(samplesFile_data, "MuonEG");
  
    MT2Analysis<MT2EstimateTree>* dataTree = new MT2Analysis<MT2EstimateTree>( "data", cfg.regionsSet() );
    MT2Analysis<MT2EstimateTree>* dataTree_of = new MT2Analysis<MT2EstimateTree>( "data_of", cfg.regionsSet() );
  
    //Filler Tree so that I don't have to rewrite the function
    MT2Analysis<MT2EstimateTree>* dataTree_filler = new MT2Analysis<MT2EstimateTree>( "data_filler", cfg.regionsSet() );

    addVariables(dataTree);      addVariables(dataTree_of); addVariables(dataTree_filler);


    if( samples_data.size()==0 ) {
      std::cout << std::endl;
      std::cout << "-> WARNING!! Didn't find any data in file: " << samplesFile_data << "!" << std::endl;
      std::cout << "-> Exiting." << std::endl;
      std::cout << std::endl;
    } else {
      for( unsigned i=0; i<samples_data.size(); ++i ) {
	computeYieldSnO( samples_data[i], cfg, dataTree, dataTree_filler);
      }
      for( unsigned i=0; i<samples_data_of.size(); ++i ) {
	computeYieldSnO( samples_data_of[i], cfg, dataTree_filler, dataTree_of);
      }
    }

    dataTree->addToFile(outputdir+"/data.root");
    dataTree_of->writeToFile(outputdir+"/data_of.root");

  } // if DATA
  
  return 0;
  
}















void roundLikeData( MT2Analysis<MT2EstimateTree>* data ) {

  std::set<MT2Region> regions = data->getRegions();
  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {
    TH1D* thisYield = data->get(*iR)->yield;
    for( int iBin=1; iBin<thisYield->GetNbinsX()+1; ++iBin ) {
      float yield = thisYield->GetBinContent(iBin);
      int yield_rounded = round(yield);
      thisYield->SetBinContent(iBin, yield_rounded  );
      thisYield->SetBinError(iBin, 0. );
    } // for bins
  } // for regions

}




void addVariables(MT2Analysis<MT2EstimateTree>* anaTree){

  MT2EstimateTree::addVar( anaTree, "Z_pt" );
  MT2EstimateTree::addVar( anaTree, "Z_phi" );
  MT2EstimateTree::addVar( anaTree, "Z_mass" );
  MT2EstimateTree::addVar( anaTree, "Z_lepId" );
  MT2EstimateTree::addVar( anaTree, "nLep" );

  MT2EstimateTree::addVar( anaTree, "lep_pt0");
  MT2EstimateTree::addVar( anaTree, "lep_pt1");
  MT2EstimateTree::addVar( anaTree, "lep_eta0");
  MT2EstimateTree::addVar( anaTree, "lep_eta1");
  MT2EstimateTree::addVar( anaTree, "raw_mt2");

  MT2EstimateTree::addVar( anaTree, "nJetHF30" );
  
}







//Loop over same and oppsite flavor just once
void computeYieldSnO( const MT2Sample& sample, const MT2Config& cfg, 
		      MT2Analysis<MT2EstimateTree>* anaTree,
		      MT2Analysis<MT2EstimateTree>* anaTree_of) {

  std::string regionsSet = cfg.regionsSet();
  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
 
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  int nentries = tree->GetEntries();
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "   Entry: " << iEntry << " / " << nentries << std::endl;
    myTree.GetEntry(iEntry);

    if(!( myTree.nlep==2 )) continue; 

    //baseline single cuts for flexibility, change back
    //if( !(myTree.passSelection("zll")) ) continue; 
    if(myTree.nVert < 0) continue;
  
  if(myTree.nJet30 < 1  ) continue;
  /*
   if(myTree.nJet30 < 2  ) continue;
   if(myTree.zll_deltaPhiMin < 0.3) continue;
    if(myTree.zll_diffMetMht > 0.5*myTree.zll_met_pt) continue;
    if( myTree.mt2 >200 ) continue;  
  */
    if(( myTree.lep_pdgId[0]*myTree.lep_pdgId[1])>0 )   continue;

    //FILTERS
    if( myTree.isData &&( myTree.Flag_HBHENoiseFilter==0 || myTree.Flag_CSCTightHaloFilter==0 || myTree.Flag_goodVertices==0 ||  myTree.Flag_eeBadScFilter==0 ) ) continue;
    if(myTree.isData && myTree.isGolden == 0) continue;

    if(myTree.lep_pt[0]<25) continue;
    if(myTree.lep_pt[1]<20) continue;

    //Need the lorentz vectors of the leptons first
    TLorentzVector *LVec = new TLorentzVector[3];
    for(int i=0; i< 2; i++){
      LVec[i].SetPtEtaPhiM(myTree.lep_pt[i], myTree.lep_eta[i],myTree.lep_phi[i], myTree.lep_mass[i]);
    }

    TLorentzVector z = LVec[0] + LVec[1]; //leptons invariant mass

    float ht   = myTree.zll_ht;
    float met  = myTree.zll_met_pt;
    float mt2  = myTree.zll_mt2;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20;

    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi(); 

    bool isSF = false;
    bool isOF = false;

    //SAME OR OPPOSITE FLAVOR selections
    if(  (myTree.lep_pdgId[0] == -myTree.lep_pdgId[1]) ) isSF = true;
    if( !(myTree.lep_pdgId[0] == -myTree.lep_pdgId[1]) ) isOF = true;
    
    if(isSF){ //////////SAME FLAVOR//////////////////////////////////////////
      if(  !(myTree.HLT_DoubleMu || myTree.HLT_DoubleEl) ) continue;
      MT2EstimateTree* thisTree = anaTree->get( myTree.zll_ht, njets, nbjets, myTree.zll_met_pt, minMTBmet, myTree.zll_mt2 );
      if (thisTree==0) continue;

      int nJetHF30_ = 0;
      for(int j=0; j<myTree.njet; ++j){
	if( myTree.jet_pt[j] < 30. || fabs(myTree.jet_eta[j]) < 3.0 ) continue;
	else ++nJetHF30_;
      }

      thisTree->assignVar("Z_pt", z.Perp() );
      thisTree->assignVar("Z_phi", z.Phi() );
      thisTree->assignVar("Z_mass", z.M() );
      thisTree->assignVar("Z_lepId", abs(myTree.lep_pdgId[0])  );

      thisTree->assignVar("nLep", myTree.nlep );
      thisTree->assignVar("lep_pt0", myTree.lep_pt[0] );
      thisTree->assignVar("lep_pt1", myTree.lep_pt[1] );
      thisTree->assignVar("lep_eta0", myTree.lep_eta[0] );
      thisTree->assignVar("lep_eta1", myTree.lep_eta[1] );
      thisTree->assignVar("raw_mt2", myTree.mt2 );

      thisTree->assignVar( "nJetHF30",  nJetHF30_ );

      thisTree->fillTree_zll(myTree, weight );
      thisTree->yield->Fill(myTree.zll_mt2, weight );

    }else if(isOF){ //////////Opposite FLAVOR//////////////////////////////////////////
      if(  myTree.isData && !(myTree.HLT_MuX_Ele12 || myTree.HLT_Mu8_EleX) ) continue;
      MT2EstimateTree* thisTree_of = anaTree_of->get( myTree.zll_ht, njets, nbjets, myTree.zll_met_pt, minMTBmet, myTree.zll_mt2 );
      if(thisTree_of==0) continue;

      int nJetHF30_ = 0;
      for(int j=0; j<myTree.njet; ++j){
	if( myTree.jet_pt[j] < 30. || fabs(myTree.jet_eta[j]) < 3.0 ) continue;
	else ++nJetHF30_;
      }

      thisTree_of->assignVar("Z_pt", z.Perp() );
      thisTree_of->assignVar("Z_phi", z.Phi() );
      thisTree_of->assignVar("Z_mass", z.M() );
      thisTree_of->assignVar("Z_lepId", abs(myTree.lep_pdgId[0])  );

      thisTree_of->assignVar("nLep", myTree.nlep );
      thisTree_of->assignVar("lep_pt0", myTree.lep_pt[0] );
      thisTree_of->assignVar("lep_pt1", myTree.lep_pt[1] );
      thisTree_of->assignVar("lep_eta0", myTree.lep_eta[0] );
      thisTree_of->assignVar("lep_eta1", myTree.lep_eta[1] );
      thisTree_of->assignVar("raw_mt2", myTree.mt2 );

      thisTree_of->assignVar( "nJetHF30",  nJetHF30_ );

      thisTree_of->fillTree_zll(myTree, weight );
      thisTree_of->yield->Fill(myTree.zll_mt2, weight );

    }else
      continue;
 
  } // for entries

  anaTree->finalize();
  anaTree_of->finalize();
  
  delete tree;

  file->Close();
  delete file;
   
}
