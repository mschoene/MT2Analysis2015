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

bool do_of = true;
bool do_bg = true;


MT2Analysis<MT2EstimateTree>* mergeYields( std::vector< MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max=-1, const std::string& legendName="" );
MT2Analysis<MT2EstimateTree>* computeYield( const MT2Sample& sample, const MT2Config& cfg, float lumi=1., bool doSameFlavor=1 );

void roundLikeData( MT2Analysis<MT2Estimate>* data );


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
  //  bool onlyMC   = false;
  bool onlyMC   = true;
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


  std::string outputdir_purity = cfg.getEventYieldDir() + "/zllControlRegionPurity";
  system(Form("mkdir -p %s", outputdir_purity.c_str()));


  std::string regionsSet = "13TeV_inclusive";
  std::cout << "-> Using regions: " << regionsSet << std::endl;

  
  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";

  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;


  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, "DYJetsToLL", 700, 799  ); // not interested in signal here

  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }


  //DATA
  std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;
  std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "Double");

  std::vector< MT2Analysis<MT2EstimateTree>* > dataTree;

  if( samples_data.size()==0 ) {
    std::cout << std::endl;
    std::cout << "-> WARNING!! Didn't find any data in file: " << samplesFile_data << "!" << std::endl;
    std::cout << "-> Exiting." << std::endl;
    std::cout << std::endl;
  } else {
    for( unsigned i=0; i<samples_data.size(); ++i ) {
      dataTree.push_back( computeYield( samples_data[i], cfg, cfg.lumi() ));
    }
  }

    MT2Analysis<MT2EstimateTree>* EventYield_data = mergeYields( dataTree, cfg.regionsSet(), "data", 0, 2000, "" );



  std::vector< MT2Analysis<MT2EstimateTree>* > EventYield;
  for( unsigned i=0; i<fSamples.size(); ++i ) 
    EventYield.push_back( computeYield( fSamples[i], cfg, cfg.lumi() ) );
   

  MT2Analysis<MT2EstimateTree>* EventYield_zll = mergeYields( EventYield, cfg.regionsSet(), "DYJets", 700, 799, "DYJets" );

  /*
  MT2Analysis<MT2EstimateTree>* Zinv = MT2Analysis<MT2EstimateTree>::readFromFile(cfg.getEventYieldDir() + "/analyses.root", "ZJets");
  if( Zinv==0 ) {
    std::cout << "-> Please run regionEventYields on MC first. I need to get the Z->vv MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }
  */

  MT2Analysis<MT2Estimate>* alpha = new MT2Analysis<MT2Estimate>( "alpha", regionsSet );
  MT2Analysis<MT2Estimate>* yield_zll = new MT2Analysis<MT2Estimate>( "Zll", regionsSet );
  *yield_zll = (* (MT2Analysis<MT2Estimate>*) EventYield_zll);

  /*  MT2Analysis<MT2Estimate>* yield_zinv = new MT2Analysis<MT2Estimate>( "ZJets", regionsSet );
  *yield_zinv = (* (MT2Analysis<MT2Estimate>*) Zinv);
  */

  EventYield_zll->writeToFile(outputdir+"/Zll_analyses.root");
  EventYield_data->addToFile(outputdir+"/Zll_analyses.root");
  /*
  yield_zll->writeToFile(outputdir+"/mc.root");
  yield_zinv->addToFile(outputdir+"/mc.root");

  (*alpha) = (*yield_zinv) / (*yield_zll);
  alpha->writeToFile(outputdir+"/data.root");

  roundLikeData(yield_zll); 
  yield_zll->addToFile(outputdir+"/data.root");
  */


  if(do_bg==true){
  //MC
  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 1, 999); // not interested in signal here
  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }
 std::vector< MT2Analysis<MT2EstimateTree>* > EventYield;
  for( unsigned i=0; i<fSamples.size(); ++i ) 
    EventYield.push_back( computeYield( fSamples[i], cfg, cfg.lumi(), 1 ) );
    
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 
  MT2Analysis<MT2EstimateTree>* EventYield_top   = mergeYields( EventYield, cfg.regionsSet(), "Top", 300, 499 ); // ttbar, single top, ttW, ttZ...
  MT2Analysis<MT2EstimateTree>* EventYield_qcd   = mergeYields( EventYield, cfg.regionsSet(), "QCD", 100, 199 );
  MT2Analysis<MT2EstimateTree>* EventYield_wjets = mergeYields( EventYield, cfg.regionsSet(), "WJets", 500, 599, "W+jets" );
  MT2Analysis<MT2EstimateTree>* EventYield_zjets = mergeYields( EventYield, cfg.regionsSet(), "ZJets", 600, 699, "Z+jets" );

  bgYields.push_back( EventYield_qcd );
  bgYields.push_back( EventYield_wjets );
  bgYields.push_back( EventYield_zjets );
  bgYields.push_back( EventYield_top );
 
  std::string outFile = outputdir + "/ZllPurityTrees.root";
  EventYield_zll->writeToFile( outFile );
  EventYield_top->addToFile( outFile );
  EventYield_qcd->addToFile( outFile );
  EventYield_wjets->addToFile( outFile );
  EventYield_zjets->addToFile( outFile );



  if(do_of==true){
    std::vector<MT2Sample> samples_data_of = MT2Sample::loadSamples(samplesFile_data, "MuonEG");
    std::vector< MT2Analysis<MT2EstimateTree>* > dataTree_of;
    for( unsigned i=0; i<samples_data.size(); ++i ) {
      dataTree_of.push_back( computeYield( samples_data_of[i], cfg, cfg.lumi(),0 ));
    }
    MT2Analysis<MT2EstimateTree>* EventYield_data_of = mergeYields( dataTree_of, cfg.regionsSet(), "data_of", 0, 2000, "" );

  std::vector< MT2Analysis<MT2EstimateTree>* > EventYield_of;
  for( unsigned i=0; i<fSamples.size(); ++i ) 
    EventYield_of.push_back( computeYield( fSamples[i], cfg, cfg.lumi(), 0 ) );
    
  MT2Analysis<MT2EstimateTree>* EventYield_zll_of = mergeYields( EventYield_of, cfg.regionsSet(), "DYJets", 700, 799, "DYJets" );

  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_of; 

  MT2Analysis<MT2EstimateTree>* EventYield_top_of   = mergeYields( EventYield_of, cfg.regionsSet(), "Top", 300, 499 ); // ttbar, single top, ttW, ttZ...
  MT2Analysis<MT2EstimateTree>* EventYield_qcd_of   = mergeYields( EventYield_of, cfg.regionsSet(), "QCD", 100, 199 );
  MT2Analysis<MT2EstimateTree>* EventYield_wjets_of = mergeYields( EventYield_of, cfg.regionsSet(), "WJets", 500, 599, "W+jets" );
  MT2Analysis<MT2EstimateTree>* EventYield_zjets_of = mergeYields( EventYield_of, cfg.regionsSet(), "ZJets", 600, 699, "Z+jets" );
 
  bgYields_of.push_back( EventYield_qcd_of );
  bgYields_of.push_back( EventYield_wjets_of );
  bgYields_of.push_back( EventYield_zjets_of );
  bgYields_of.push_back( EventYield_top_of );

  std::string outFile_of = outputdir + "/ZllPurityTrees_of.root";
  EventYield_zll_of->writeToFile( outFile_of );
  EventYield_top_of->addToFile( outFile_of );
  EventYield_qcd_of->addToFile( outFile_of );
  EventYield_wjets_of->addToFile( outFile_of );
  EventYield_zjets_of->addToFile( outFile_of );

  std::string outFile_data_of = outputdir + "/ZllPurityTrees_data_of.root";
  EventYield_data_of->writeToFile(outFile_data_of);
  




  }//etz langets aber mit of

  }//gnueg bg abgschetzt
  return 0;

}





void roundLikeData( MT2Analysis<MT2Estimate>* data ) {

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

MT2Analysis<MT2EstimateTree>* computeYield( const MT2Sample& sample, const MT2Config& cfg, float lumi, bool doSameFlavor ) {


  std::string regionsSet = cfg.regionsSet();

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
  

  MT2Tree myTree;
  if( cfg.additionalStuff()=="qgVars" ) {
     myTree.loadGenStuff = true;
  } else {
    myTree.loadGenStuff = false;
  }
  myTree.Init(tree);



  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  MT2Analysis<MT2EstimateTree>* analysis = new MT2Analysis<MT2EstimateTree>( sample.sname, regionsSet, sample.id );


  MT2EstimateTree::addVar( analysis, "Z_pt" );
  MT2EstimateTree::addVar( analysis, "Z_phi" );
  MT2EstimateTree::addVar( analysis, "Z_mass" );
  MT2EstimateTree::addVar( analysis, "Z_lepId" );
  MT2EstimateTree::addVar( analysis, "nLep" );

  MT2EstimateTree::addVar( analysis, "lep_pt0");
  MT2EstimateTree::addVar( analysis, "lep_pt1");
  MT2EstimateTree::addVar( analysis, "lep_eta0");
  MT2EstimateTree::addVar( analysis, "lep_eta1");
  MT2EstimateTree::addVar( analysis, "raw_mt2");
  
  MT2EstimateTree::addVar( analysis, "HLT_DoubleMu");
  MT2EstimateTree::addVar( analysis, "HLT_DoubleEl");
  


  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "   Entry: " << iEntry << " / " << nentries << std::endl;
    myTree.GetEntry(iEntry);

    //if( !(myTree.passSelection("zll")) ) continue; 

    if(!( myTree.nlep==2 )) continue; 
    //baseline

    if(myTree.nVert < 0) continue;
    if(myTree.nJet30 < 2  ) continue;
    if(myTree.zll_deltaPhiMin < 0.3) continue;
    if(myTree.zll_diffMetMht > 0.5*myTree.zll_met_pt) continue;


    if( myTree.mt2 >200 ) continue; //change back when more data

    if(doSameFlavor==1 && !(myTree.lep_pdgId[0] == -myTree.lep_pdgId[1]) ) continue;
    if(doSameFlavor==0 &&  (myTree.lep_pdgId[0] == -myTree.lep_pdgId[1]) ) continue;
 
    if(( myTree.lep_pdgId[0]*myTree.lep_pdgId[1])>0 )   continue;
    
    if(  doSameFlavor==1 && !(myTree.HLT_DoubleMu || myTree.HLT_DoubleEl) ) continue;
    if(  doSameFlavor==0 && !(myTree.HLT_DoubleMu || myTree.HLT_DoubleEl ) ) continue;

    //Implemented in computeZllGammaRatio for now to allow flexibility
    //  if(myTree.lep_pt[0]<25) continue;
    //  if(myTree.lep_pt[1]<20) continue;

    //Need the lorentz vectors of the leptons first
    TLorentzVector *LVec = new TLorentzVector[5];
    for(int i=0; i< 2; i++){
      LVec[i].SetPtEtaPhiM(myTree.lep_pt[i], myTree.lep_eta[i],myTree.lep_phi[i], myTree.lep_mass[i]);
    }

    double Z_invM_true = 91.19;
    TLorentzVector z = LVec[0] + LVec[1]; //leptons invariant mass
    double M_ll = z.M(); //Z mass

    //  if( abs(M_ll - Z_invM_true)>20.) continue;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = myTree.mt2;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20;

    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi(); 

    MT2EstimateTree* thisEstimate = analysis->get( myTree.zll_ht, njets, nbjets, myTree.zll_met_pt, minMTBmet, myTree.zll_mt2 );
    if( thisEstimate==0 ) continue; 


    //initialize
    thisEstimate->assignVar("Z_pt", z.Perp() );
    thisEstimate->assignVar("Z_phi", z.Phi() );
    thisEstimate->assignVar("Z_mass", z.M() );
    thisEstimate->assignVar("Z_lepId", abs(myTree.lep_pdgId[0])  );

    thisEstimate->assignVar("nLep", myTree.nlep );
    thisEstimate->assignVar("lep_pt0", myTree.lep_pt[0] );
    thisEstimate->assignVar("lep_pt1", myTree.lep_pt[1] );
    thisEstimate->assignVar("lep_eta0", myTree.lep_eta[0] );
    thisEstimate->assignVar("lep_eta1", myTree.lep_eta[1] );
    thisEstimate->assignVar("raw_mt2", myTree.mt2 );

    thisEstimate->assignVar("HLT_DoubleMu", myTree.HLT_DoubleMu );
    thisEstimate->assignVar("HLT_DoubleEl", myTree.HLT_DoubleEl );

    //JUGJUGJUG
    thisEstimate->fillTree_zll(myTree, weight );

    thisEstimate->yield->Fill(myTree.zll_mt2, weight );
  
  } // for entries

  analysis->finalize();
  
  delete tree;

  file->Close();
  delete file;
  
  return analysis;

}

MT2Analysis<MT2EstimateTree>* mergeYields( std::vector<MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max, const std::string& legendName ) {

  if( id_max<0 ) id_max=id_min;

  MT2Analysis<MT2EstimateTree>* return_EventYield = new MT2Analysis<MT2EstimateTree>(name, regionsSet, id_min, legendName);

  for( unsigned i=0; i<EventYield.size(); ++i ) {

    if( EventYield[i]->id >= id_min && EventYield[i]->id <= id_max ) {

       *(return_EventYield) += *(EventYield[i]);

     }

  } // for EventYield


  return return_EventYield;

}




























