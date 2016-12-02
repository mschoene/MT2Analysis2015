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

bool doZinvEst = true;
bool do_bg = true;


TH1D*  h_muTrk_hi = 0;
TH2D*  h_elTrk = 0;


void computeYieldSnO( const MT2Sample& sample, const MT2Config& cfg,   
		      MT2Analysis<MT2EstimateTree>* anaTree,  
		      MT2Analysis<MT2EstimateTree>* anaTree_of,
		      TH2D* h_elSF, TH2D* h_muSF, bool doZinvEst );
void addVariables(MT2Analysis<MT2EstimateTree>* anaTree);
void roundLikeData( MT2Analysis<MT2EstimateTree>* data );


int main(int argc, char* argv[]) {

  std::cout << std::endl << std::endl;
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

  bool onlyData   = false;
  bool onlyMC     = false;
  bool onlySignal = false;
  if( argc > 2 ) {
    std::string dataMC(argv[2]);
    if( dataMC=="data" ) onlyData = true;
    else if( dataMC=="MC" ) onlyMC = true;
    else if( dataMC=="signal" ) onlySignal = true;
    else {
      std::cout << "-> You passed a second argument that isn't 'data' nor 'MC', so I don't know what to do about it." << std::endl;
    }
  }

  TH1::AddDirectory(kFALSE); //stupid ROOT memory allocation needs this

  std::string outputdir = cfg.getEventYieldDir() + "/zllControlRegion";
  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string regionsSet;// = "13TeV_inclusive";
  regionsSet=cfg.crRegionsSet();
  // regionsSet=cfg.zllRegions();
  // std::string regionsSet = cfg.zllRegions();

  std::cout << "-> Using regions: " << regionsSet << std::endl;

  //Getting the scale factor histogram/////////////////
  //Electrons//
  std::string filename = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/scaleFactors.root";
  TFile * f_ele = new TFile(filename.c_str() );
  if (!f_ele->IsOpen()) std::cout << " ERROR: Could not find scale factor file " << filename << std::endl; 
  //Uncomment for loose Id
  //TH2D* h_id = (TH2D*) f_ele->Get("CutBasedLoose");
  //(TH2D*) f_ele->Get("CutBasedVeto");
  TH2D* h_id = (TH2D*) f_ele->Get("GsfElectronToVeto");
  TH2D* h_iso = (TH2D*) f_ele->Get("MVAVLooseElectronToMini");
  if (!h_id || !h_iso) std::cout << "ERROR: Could not find scale factor histogram"<< std::endl;
  TH2D* h_elSF = (TH2D*) h_id->Clone("h_elSF");
  h_elSF->SetDirectory(0);
  h_elSF->Multiply(h_iso);

  std::string filenameElTrk = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/egammaEffi_SF2D.root";
  TFile * f_eleTrk = new TFile(filenameElTrk.c_str() );
  if (!f_eleTrk->IsOpen()) std::cout << " ERROR: Could not find scale factor file " << filenameElTrk << std::endl; 
  h_elTrk = (TH2D*) f_eleTrk->Get("EGamma_SF2D");
  h_elTrk->SetDirectory(0);
  f_eleTrk->Close(); delete f_eleTrk; 




  //Muons//
  std::string filenameID = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/TnP_MuonID_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root";
  std::string filenameISO = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/TnP_MuonID_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root";
  std::string filenamedxyz = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/TnP_MuonID_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root";
  TFile * f1 = new TFile(filenameID.c_str() );
  TFile * f2 = new TFile(filenameISO.c_str() );
  TFile * f3 = new TFile(filenamedxyz.c_str() );
  if (!f1->IsOpen()) { std::cout<<" ERROR: Could not find ID scale factor file "<<filenameID<<std::endl; return 0;}
  if (!f2->IsOpen()) { std::cout<<"ERROR: Could not find ISO scale factor file "<<filenameISO<<std::endl; return 0;}
  if (!f3->IsOpen()) { std::cout<<"ERROR: Could not find dxy dz scale factor file "<<filenamedxyz<<std::endl; return 0;}
  TH2D* h_id_mu = (TH2D*) f1->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0");
  TH2D* h_iso_mu = (TH2D*) f2->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_PF_pass");
  TH2D* h_dxyz_mu = (TH2D*) f3->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_PF_pass");
  if (!h_id_mu || !h_iso_mu  || !h_dxyz_mu) { std::cout<<"ERROR: Could not find scale factor histogram"<<std::endl; return 0;}
  TH2D* h_muSF = (TH2D*) h_id_mu->Clone("h_muSF");
  h_muSF->SetDirectory(0);
  h_muSF->Multiply(h_iso_mu);
  h_muSF->Multiply(h_dxyz_mu);
 


  TH1D* h_trk_mu_hi = 0;

  std::string filenameTrk = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/general_tracks_and_early_general_tracks_corr_ratio.root";
  TFile * fTrk = new TFile(filenameTrk.c_str() );
  if (!fTrk->IsOpen()) { std::cout<<" ERROR: Could not find track ineff scale factor file "<<filenameTrk<<std::endl; return 0;}
  h_trk_mu_hi = (TH1D*) fTrk->Get("mutrksfptg10");
  if (!h_trk_mu_hi) { std::cout<<"ERROR: Could not find trk sf histogram"<<std::endl; return 0;}
  h_muTrk_hi = (TH1D*) h_trk_mu_hi->Clone("h_muTrk_hi");
  h_muTrk_hi->SetDirectory(0);
  // fTrk->Close(); delete fTrk;


  if( onlySignal ){
    onlyData=true;
    onlyMC=true;
  }


   if( cfg.useMC() && !onlyData ) { // run on MC  
  
    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;


    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 700, 799); // DY signal only
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(1209);
    }

    MT2Analysis<MT2EstimateTree>* mcTree = new MT2Analysis<MT2EstimateTree>( "zllCR", cfg.crRegionsSet() );
    addVariables(mcTree); //Adds some additional variables Zpt,Zmass, raw MT2...

    MT2Analysis<MT2EstimateTree>* mcTree_of = new MT2Analysis<MT2EstimateTree>( "zllCR_of", cfg.crRegionsSet() );
    addVariables(mcTree_of);
  
    for( unsigned i=0; i<fSamples.size(); ++i ) 
      computeYieldSnO( fSamples[i], cfg, mcTree, mcTree_of, h_elSF, h_muSF, false);
 
    mcTree->writeToFile(outputdir+"/mc.root");
    mcTree_of->writeToFile(outputdir+"/mc_of.root");

    if(doZinvEst){
      MT2Analysis<MT2EstimateTree>* mcTree_forZinvEst = new MT2Analysis<MT2EstimateTree>( "zllCR", cfg.regionsSet() );
      addVariables(mcTree_forZinvEst); //Adds some additional variables Zpt,Zmass, raw MT2...

      MT2Analysis<MT2EstimateTree>* mcTree_of_forZinvEst = new MT2Analysis<MT2EstimateTree>( "zllCR_of", cfg.regionsSet() );
      addVariables(mcTree_of_forZinvEst);
  
      for( unsigned i=0; i<fSamples.size(); ++i ) 
	computeYieldSnO( fSamples[i], cfg, mcTree_forZinvEst, mcTree_of_forZinvEst, h_elSF, h_muSF, true);
 
      mcTree_forZinvEst->writeToFile(outputdir+"/mc_forZinvEst.root");
      mcTree_of_forZinvEst->writeToFile(outputdir+"/mc_of_forZinvEst.root");
    }
    
    if( cfg.dummyAnalysis() ) {
      roundLikeData(mcTree); 
      mcTree->addToFile(outputdir+"/data.root");
      roundLikeData(mcTree_of); 
      mcTree_of->addToFile(outputdir+"/data_of.root");    
    }
    
    if(do_bg==true){
      //MC
      MT2Analysis<MT2EstimateTree>* mc_top = new MT2Analysis<MT2EstimateTree>( "Top", cfg.crRegionsSet(),300, "Top" );
      MT2Analysis<MT2EstimateTree>* mc_top_of = new MT2Analysis<MT2EstimateTree>( "Top", cfg.crRegionsSet(),300, "Top" );
      addVariables(mc_top);      addVariables(mc_top_of);
      std::vector<MT2Sample> fSamples_top = MT2Sample::loadSamples(samplesFileName, 300, 499);   
      for( unsigned i=0; i<fSamples_top.size(); ++i )
	computeYieldSnO( fSamples_top[i], cfg, mc_top, mc_top_of, h_elSF, h_muSF, false);

      /* 
	 MT2Analysis<MT2EstimateTree>* mc_qcd = new MT2Analysis<MT2EstimateTree>( "QCD", cfg.crRegionsSet(),100, "QCD" );
	 MT2Analysis<MT2EstimateTree>* mc_qcd_of = new MT2Analysis<MT2EstimateTree>( "QCD", cfg.crRegionsSet(),100, "QCD");
	 addVariables(mc_qcd);      addVariables(mc_qcd_of);
	 std::vector<MT2Sample> fSamples_qcd = MT2Sample::loadSamples(samplesFileName, 100, 199);   
	 for( unsigned i=0; i<fSamples_qcd.size(); ++i )
	 computeYieldSnO( fSamples_qcd[i], cfg, mc_qcd, mc_qcd_of, h_elSF, h_muSF, false);

	 MT2Analysis<MT2EstimateTree>* mc_wjets = new MT2Analysis<MT2EstimateTree>( "WJets", cfg.crRegionsSet(),500, "W+jets"  );
	 MT2Analysis<MT2EstimateTree>* mc_wjets_of = new MT2Analysis<MT2EstimateTree>( "WJets", cfg.crRegionsSet(),500, "W+jets");
	 addVariables(mc_wjets);      addVariables(mc_wjets_of);
	 std::vector<MT2Sample> fSamples_wjets = MT2Sample::loadSamples(samplesFileName, 500, 599);   
	 for( unsigned i=0; i<fSamples_wjets.size(); ++i )
	 computeYieldSnO( fSamples_wjets[i], cfg, mc_wjets, mc_wjets_of, h_elSF, h_muSF, false);
      */

      MT2Analysis<MT2EstimateTree>* mc_zll   = mcTree;
      mc_zll->setName("DYJets");
      mc_zll->setFullName("DY+jets");
      
      std::string outFile = outputdir + "/ZllPurityTrees.root";
      mc_zll->writeToFile( outFile );
      mc_top->addToFile( outFile );
      //      mc_qcd->addToFile( outFile );
      //      mc_wjets->addToFile( outFile );

      //For the OPPOSITE FLAVOR EVENTS:
      MT2Analysis<MT2EstimateTree>* mc_zll_of   = mcTree_of;
      mc_zll_of->setName("DYJets");
      mc_zll_of->setFullName("DY+jets");
 
      std::string outFile_of = outputdir + "/ZllPurityTrees_of.root";
      mc_zll_of->writeToFile( outFile_of );
      mc_top_of->addToFile( outFile_of );
      //     mc_qcd_of->addToFile( outFile_of );
      //     mc_wjets_of->addToFile( outFile_of );

      if(doZinvEst){
	MT2Analysis<MT2EstimateTree>* mc_top_forZinvEst = new MT2Analysis<MT2EstimateTree>( "zllCR", cfg.regionsSet() );
	addVariables(mc_top_forZinvEst); //Adds some additional variables Zpt,Zmass, raw MT2...

	MT2Analysis<MT2EstimateTree>* mc_top_of_forZinvEst = new MT2Analysis<MT2EstimateTree>( "zllCR_of", cfg.regionsSet() );
	addVariables(mc_top_of_forZinvEst);
  
	for( unsigned i=0; i<fSamples_top.size(); ++i ) 
	  computeYieldSnO( fSamples_top[i], cfg, mc_top_forZinvEst, mc_top_of_forZinvEst, h_elSF, h_muSF, true);
 
	mc_top_forZinvEst->writeToFile(outputdir+"/mc_Top_forZinvEst.root");
	mc_top_of_forZinvEst->writeToFile(outputdir+"/mc_Top_of_forZinvEst.root");
      }
    
    } //End do background trees
    
  } //if only MC
  
  if( !onlyMC ) {
    
    //DATA
    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;
    //    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "merged");  
    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "Double");  

    // std::vector<MT2Sample> samples_data_of = MT2Sample::loadSamples(samplesFile_data, "merged");
  
    MT2Analysis<MT2EstimateTree>* dataTree = new MT2Analysis<MT2EstimateTree>( "data", cfg.crRegionsSet() );
    MT2Analysis<MT2EstimateTree>* dataTree_of = new MT2Analysis<MT2EstimateTree>( "data_of", cfg.crRegionsSet() );
  
    //Filler Tree so that I don't have to rewrite the function
    // MT2Analysis<MT2EstimateTree>* dataTree_filler = new MT2Analysis<MT2EstimateTree>( "data_filler", cfg.crRegionsSet() );

    addVariables(dataTree);      addVariables(dataTree_of); // addVariables(dataTree_filler);

    if( samples_data.size()==0 ) {
      std::cout << std::endl;
      std::cout << "-> WARNING!! Didn't find any data in file: " << samplesFile_data << "!" << std::endl;
      std::cout << "-> Exiting." << std::endl;
      std::cout << std::endl;
    } else {
      for( unsigned i=0; i<samples_data.size(); ++i ) 
	computeYieldSnO( samples_data[i], cfg, dataTree, dataTree_of, h_elSF, h_muSF, false);
      //      computeYieldSnO( samples_data[i], cfg, dataTree, dataTree_filler, h_elSF, h_muSF, false);

      //for( unsigned i=0; i<samples_data_of.size(); ++i )
      //	computeYieldSnO( samples_data_of[i], cfg, dataTree_filler, dataTree_of, h_elSF, h_muSF, false);
    }

    dataTree->addToFile(outputdir+"/data.root");
    dataTree_of->writeToFile(outputdir+"/data_of.root");




    if(doZinvEst){
      MT2Analysis<MT2EstimateTree>* dataTree_forZinvEst = new MT2Analysis<MT2EstimateTree>( "data", cfg.regionsSet() );
      MT2Analysis<MT2EstimateTree>* dataTree_of_forZinvEst = new MT2Analysis<MT2EstimateTree>( "data_of", cfg.regionsSet() );
      addVariables(dataTree_forZinvEst); addVariables(dataTree_of_forZinvEst);

      for( unsigned i=0; i<samples_data.size(); ++i ) 
	computeYieldSnO( samples_data[i], cfg, dataTree_forZinvEst, dataTree_of_forZinvEst, h_elSF, h_muSF, true);

      dataTree_forZinvEst->addToFile(outputdir+"/data_forZinvEst.root");
      dataTree_of_forZinvEst->addToFile(outputdir+"/data_of_forZinvEst.root");
    }



  } // if DATA


  if( onlySignal){
    std::cout << "DOING THE SIGNAL" << std::endl;

    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;
    
    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 999, 2000); // signal only
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(1209);
    }

    MT2Analysis<MT2EstimateTree>* signalTree = new MT2Analysis<MT2EstimateTree>( "zllSigCR", cfg.regionsSet() );
    addVariables(signalTree); //Adds some additional variables Zpt,Zmass, raw MT2...

    MT2Analysis<MT2EstimateTree>* signalTree_of = new MT2Analysis<MT2EstimateTree>( "zllSigCR_of", cfg.regionsSet() );
    addVariables(signalTree_of);
  
    for( unsigned i=0; i<fSamples.size(); ++i ) 
      computeYieldSnO( fSamples[i], cfg, signalTree, signalTree_of, h_elSF, h_muSF, true);
 
    signalTree->addToFile(outputdir+"/signal_forZinvEst.root");

    signalTree_of->addToFile(outputdir+"/signal_of_forZinvEst.root");


  }//end of only signal

  
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

  MT2EstimateTree::addVar( anaTree, "ID" );

  MT2EstimateTree::addVar( anaTree, "Z_pt" );
  MT2EstimateTree::addVar( anaTree, "Z_phi" );
  MT2EstimateTree::addVar( anaTree, "Z_eta" );
  MT2EstimateTree::addVar( anaTree, "Z_mass" );
  MT2EstimateTree::addVar( anaTree, "Z_lepId" );
  MT2EstimateTree::addVar( anaTree, "nLep" );

  MT2EstimateTree::addVar( anaTree, "lep_tightId0" );
  MT2EstimateTree::addVar( anaTree, "lep_tightId1" );

  MT2EstimateTree::addVar( anaTree, "lep_pdgId0");
  MT2EstimateTree::addVar( anaTree, "lep_pdgId1");
  MT2EstimateTree::addVar( anaTree, "lep_pt0");
  MT2EstimateTree::addVar( anaTree, "lep_pt1");
  MT2EstimateTree::addVar( anaTree, "lep_eta0");
  MT2EstimateTree::addVar( anaTree, "lep_eta1");
  MT2EstimateTree::addVar( anaTree, "lep_phi0");
  MT2EstimateTree::addVar( anaTree, "lep_phi1");
  MT2EstimateTree::addVar( anaTree, "raw_mt2"); // = mt2 with the two leptons
  MT2EstimateTree::addVar( anaTree, "raw_met"); // = met with the two leptons

  MT2EstimateTree::addVar( anaTree, "weight_lep0"); 
  MT2EstimateTree::addVar( anaTree, "weight_lep1");
  MT2EstimateTree::addVar( anaTree, "weight_lep_err");

  MT2EstimateTree::addVar( anaTree, "HLT_weight");

  MT2EstimateTree::addVar( anaTree, "nJetHF30" );
  MT2EstimateTree::addVar( anaTree, "jet1_pt" );
  
}







//Loop over same and oppsite flavor just once
void computeYieldSnO( const MT2Sample& sample, const MT2Config& cfg, 
		      MT2Analysis<MT2EstimateTree>* anaTree,
		      MT2Analysis<MT2EstimateTree>* anaTree_of,
		      TH2D* h_elSF, TH2D* h_muSF, bool doZinvEst) {

  std::string regionsSet = cfg.crRegionsSet();
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

    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20;
    float ht   = myTree.zll_ht;
    //float met  = myTree.zll_met_pt;
    float mt2  = (njets>1) ? myTree.zll_mt2 : myTree.zll_ht;
    float minMTBmet = myTree.minMTBMet;

    int ID = myTree.evt_id;

    //Minimal selection for the standard model Z/Gamma ratio
    if(myTree.nVert < 1) continue;

    if( cfg.analysisType() == "mt2"){
      if( regionsSet!="13TeV_noCut" )
        if( !myTree.passSelection("zll") ) continue;
    }

    if(( myTree.lep_pdgId[0]*myTree.lep_pdgId[1])>0 )   continue;
    
    if(  myTree.nJet30==1 && !(myTree.jet_id[0]>=4)) continue;    



    //FILTERS
    if( myTree.isData && !myTree.passFilters() ) continue;
    if( myTree.isData &&  myTree.isGolden == 0 ) continue;

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

    // if(myTree.lep_pt[0]<35) continue;
    // if(myTree.lep_pt[1]<35) continue; 
    if(myTree.lep_pt[0]<100) continue;
    if(myTree.lep_pt[1]<30) continue; 

    //    if(myTree.lep_pt[0]<25) continue;
    //    if(myTree.lep_pt[1]<20) continue; 

    
    //Need the lorentz vectors of the leptons first
    TLorentzVector *LVec = new TLorentzVector[3];
    for(int i=0; i< 2; i++){
      LVec[i].SetPtEtaPhiM(myTree.lep_pt[i], myTree.lep_eta[i],myTree.lep_phi[i], myTree.lep_mass[i]);
    }

    TLorentzVector z = LVec[0] + LVec[1]; //leptons invariant mass
    


    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi(); 
    //Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi(); 

    if(ID >999)
      weight = 1000.*myTree.evt_xsec / nentries;
    
    //get the lepton scale factors
    Double_t weight_lep0 = 1.;
    Double_t weight_lep1 = 1.;
    Double_t weight_lep_err = 1.;





    if( !myTree.isData ){ 
      weight *= myTree.weight_btagsf;
      
      /* //ONLY for 74MC which has xsec corrections //temporarily scaling by hand the cross sections
      if( myTree.evt_id == 702) weight = weight * 1.0573;
      if( myTree.evt_id == 703) weight = weight * 0.9588;
      if( myTree.evt_id == 704) weight = weight * 1.0329;
      if( myTree.evt_id == 705) weight = weight * 0.9945;
      */

      //weight *= myTree.weight_toppt;
    
      /////////Add lepton scale factor//////////
      Float_t central = 1; 
      Float_t err = 0;
      Float_t uncert_UP = 0; //	Float_t uncert_DN = 0; 

      Float_t weight_lepsf = 1.;
      Float_t weight_lepsf_UP = 1.;

      for(int o=0; o < 2; ++o){

	float pt = myTree.lep_pt[o];
	float eta = fabs( myTree.lep_eta[o] );
	int pdgId = abs( myTree.lep_pdgId[o] );

	float pt_cutoff = std::max( 10.1, std::min( 100., double(pt) ) );
 
	//Electrons
	if( pdgId == 11) {
	  Int_t binx = h_elSF->GetXaxis()->FindBin(pt_cutoff);
	  Int_t biny = h_elSF->GetYaxis()->FindBin(eta);
	  central = h_elSF->GetBinContent(binx,biny);
	  err  = h_elSF->GetBinError(binx,biny);
	  if (central > 1.2 || central < 0.8) 
	    std::cout<<"STRANGE: Electron with pT/eta of "<< pt <<"/"<< eta <<". SF is "<< central <<std::endl;

	  float central_trk = 1;
	  Int_t binx_trk = h_elTrk->GetXaxis()->FindBin(  myTree.lep_eta[o] );
	  if( binx_trk>28 ) binx_trk = 28; //though we shouldn't really get many electrons with eta = 2.5
	  else if( binx_trk<1 ) binx_trk = 1;
	  central_trk = h_elTrk->GetBinContent( binx_trk,1  ); // y=pT range of 20-200

	  central *= central_trk;

	  uncert_UP = central + err;
	  //	uncert_DN = central - err;

	} //else Muons
	else if ( pdgId == 13) {
	  Int_t binx = h_muSF->GetXaxis()->FindBin(pt);
	  Int_t biny = h_muSF->GetYaxis()->FindBin(fabs(eta));
	  if ( binx >7 ) binx = 7; //overflow bin empty for the muons...
	  central = h_muSF->GetBinContent(binx,biny);
	  err  = 0.03; // adding in quadrature 1% unc. on ID and 1% unc. on ISO
	  if (central > 1.3 || central < 0.7) 
	    std::cout<<"STRANGE: Muon with pT/eta of "<<pt<<"/"<< fabs(eta) <<". SF is "<< central <<std::endl;

	  float central_trk = 1;
	  Int_t binx_trk = h_muTrk_hi->GetXaxis()->FindBin(  myTree.lep_eta[o] );
	  if( binx_trk>10 ) binx_trk = 10;
	  else if( binx_trk<1 ) binx_trk = 1;
	  central_trk = h_muTrk_hi->GetBinContent( binx_trk );

	  central *= central_trk;

	  uncert_UP = central + err;
	  //	uncert_DN = central - err;

	}//done with one  electron/muon 
	weight_lepsf    *= central;
	weight_lepsf_UP *= uncert_UP;
	//  weight_lepsf_DN *= uncert_DN;	

	//Backwards compatible, don't make me think now, it's too warm
	weight_lep0 = central;
	weight_lep_err = uncert_UP;
	weight_lep1 = 1.;  	
      }//end of loop over objects

    }//end of applying SF




    bool isSF = false;
    bool isOF = false;

    //SAME OR OPPOSITE FLAVOR selections
    if(  (myTree.lep_pdgId[0] == -myTree.lep_pdgId[1]) ) isSF = true;
    if( !(myTree.lep_pdgId[0] == -myTree.lep_pdgId[1]) ) isOF = true;
    
    if(isSF){ //////////SAME FLAVOR//////////////////////////////////////////
      if(  myTree.isData && !( myTree.HLT_DoubleMu || myTree.HLT_DoubleMu_NonIso || myTree.HLT_SingleMu_NonIso || myTree.HLT_DoubleEl || myTree.HLT_DoubleEl33 || myTree.HLT_Photon165_HE10 ) )continue;

      if(doZinvEst){
	//SF part
	if( fabs(myTree.zll_mass-91.19)>=20 ) continue;
	if( myTree.zll_pt <= 200. ) continue;
	//if( fabs(z.M()-91.19)>=20 ) continue;
	//if( z.Perp() <= 180. ) continue;
	//if( fabs(z.M()-91.19)>10 ) continue;
      }

      if( abs(myTree.lep_pdgId[0])==11 && myTree.lep_tightId[0]< 0.5 ) continue;
      if( abs(myTree.lep_pdgId[1])==11 && myTree.lep_tightId[1]< 0.5 ) continue;
      

      //      if(myTree.nVert > 0 && myTree.nJet30 >= 1 && myTree.zll_ht>250. && ((myTree.nJet30>1 && myTree.zll_mt2>200.)||myTree.nJet30==1) && myTree.nJet30FailId == 0 && myTree.zll_deltaPhiMin > 0.3 && ( (myTree.nJet30>1 && myTree.zll_ht<1000. && myTree.zll_met_pt>250.) || (myTree.nJet30>1 && myTree.zll_ht>=1000. && myTree.zll_met_pt>30.) || (myTree.nJet30==1 &&myTree.zll_met_pt>250.)) &&  myTree.zll_diffMetMht < 0.5*myTree.zll_met_pt &&  myTree.nlep==2 && (myTree.lep_pdgId[0]*myTree.lep_pdgId[1])<0 && ((myTree.nJet30==1 && myTree.jet_id[0]>=4) || myTree.nJet30>1) && myTree.lep_pt[0]>=25. && myTree.lep_pt[1]>=20. && myTree.Flag_HBHENoiseFilter>0 && myTree.Flag_HBHENoiseIsoFilter>0 && myTree.Flag_globalTightHalo2016Filter>0 && myTree.Flag_EcalDeadCellTriggerPrimitiveFilter>0 && myTree.Flag_goodVertices>0 && myTree.Flag_eeBadScFilter>0 && myTree.Flag_badMuonFilter>0 && myTree.Flag_badChargedHadronFilter>0 && (myTree.HLT_DoubleMu || myTree.HLT_DoubleMu_NonIso || myTree.HLT_SingleMu_NonIso || myTree.HLT_DoubleEl || myTree.HLT_DoubleEl33 || myTree.HLT_Photon165_HE10 ) && (abs(myTree.lep_pdgId[0]) == abs(myTree.lep_pdgId[1])) && ((abs(myTree.lep_pdgId[0])==11 && myTree.lep_tightId[0]> 0.5) || (abs(myTree.lep_pdgId[0])==13)) && ((abs(myTree.lep_pdgId[1])==11 && myTree.lep_tightId[1]> 0.5) || (abs(myTree.lep_pdgId[1])==13)) && myTree.jet_pt[0]<13000. && myTree.met_pt/myTree.met_caloPt < 5. );
      //      else continue;


      MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
      if (thisTree==0) continue;

      int nJetHF30_ = 0;
      for(int j=0; j<myTree.njet; ++j){
	if( myTree.jet_pt[j] < 30. || fabs(myTree.jet_eta[j]) < 3.0 ) continue;
	else ++nJetHF30_;
      }

      float HLT_weight = 1;
      if( !myTree.isData){
	if( abs(myTree.lep_pdgId[0])==11 )
	  HLT_weight = 0.993;
	else if(abs(myTree.lep_pdgId[0])==13 )
	  HLT_weight = 0.969;
	
	weight = weight* (HLT_weight * weight_lep0);
      }


      thisTree->assignVar("ID", ID );

      //      thisTree->assignVar("Z_pt", z.Perp() );
      thisTree->assignVar("Z_pt", myTree.zll_pt );
      thisTree->assignVar("Z_phi", z.Phi() );
      thisTree->assignVar("Z_eta", z.Eta() );
      //      thisTree->assignVar("Z_mass", z.M() );
      thisTree->assignVar("Z_mass", myTree.zll_mass );
      thisTree->assignVar("Z_lepId", abs(myTree.lep_pdgId[0]) );

      thisTree->assignVar("nLep", myTree.nlep );
      thisTree->assignVar("lep_pdgId0", myTree.lep_pdgId[0] );
      thisTree->assignVar("lep_pdgId1", myTree.lep_pdgId[1] );
      thisTree->assignVar("lep_pt0", myTree.lep_pt[0] );
      thisTree->assignVar("lep_pt1", myTree.lep_pt[1] );
      thisTree->assignVar("lep_eta0", myTree.lep_eta[0] );
      thisTree->assignVar("lep_eta1", myTree.lep_eta[1] );
      thisTree->assignVar("lep_phi0", myTree.lep_phi[0] );
      thisTree->assignVar("lep_phi1", myTree.lep_phi[1] );
      thisTree->assignVar("raw_mt2", myTree.mt2 );
      thisTree->assignVar("raw_met", myTree.met_pt );

      thisTree->assignVar("weight_lep0", weight_lep0);
      thisTree->assignVar("weight_lep1", weight_lep1);
      thisTree->assignVar("weight_lep_err", weight_lep_err);

      thisTree->assignVar("HLT_weight", HLT_weight );
  
      thisTree->assignVar( "nJetHF30", nJetHF30_ );
      thisTree->assignVar( "jet1_pt", myTree.jet1_pt );

      thisTree->assignVar("lep_tightId0", myTree.lep_tightId[0] );
      thisTree->assignVar("lep_tightId1", myTree.lep_tightId[1] );

 
      thisTree->fillTree_zll(myTree, weight );
      thisTree->yield->Fill( mt2, weight );

    } else if(isOF){ //////////Opposite FLAVOR//////////////////////////////////////////
      if(  myTree.isData && !( (myTree.HLT_MuX_Ele12 || myTree.HLT_Mu8_EleX || myTree.HLT_Mu33_Ele33_NonIso || myTree.HLT_Mu30_Ele30_NonIso) ) ) continue;
      //if(  myTree.isData && !( (myTree.HLT_MuX_Ele12 || myTree.HLT_Mu8_EleX || myTree.HLT_Mu30_Ele30_NonIso )) ) continue;

      if(doZinvEst){
	//SF part
	if( fabs(myTree.zll_mass-91.19)>=20. ) continue;
	if( myTree.zll_pt <= 200. ) continue;
//	if( fabs(z.M()-91.19)>=20 ) continue;
//	if( z.Perp() <= 180. ) continue;
	//if( fabs(z.M()-91.19)>10 ) continue;
      }


      if( abs(myTree.lep_pdgId[0])==11 && myTree.lep_tightId[0]< 0.5 ) continue;
      if( abs(myTree.lep_pdgId[1])==11 && myTree.lep_tightId[1]< 0.5 ) continue;


      //      if(myTree.nVert > 0 && myTree.nJet30 >= 1 && myTree.zll_ht>250. && ((myTree.nJet30>1 && myTree.zll_mt2>200.)||myTree.nJet30==1) && myTree.nJet30FailId == 0 && myTree.zll_deltaPhiMin > 0.3 && ( (myTree.nJet30>1 && myTree.zll_ht<1000. && myTree.zll_met_pt>250.) || (myTree.nJet30>1 && myTree.zll_ht>=1000. && myTree.zll_met_pt>30.) || (myTree.nJet30==1 &&myTree.zll_met_pt>250.)) &&  myTree.zll_diffMetMht < 0.5*myTree.zll_met_pt &&  myTree.nlep==2 && (myTree.lep_pdgId[0]*myTree.lep_pdgId[1])<0 && ((myTree.nJet30==1 && myTree.jet_id[0]>=4) || myTree.nJet30>1) && myTree.lep_pt[0]>=25. && myTree.lep_pt[1]>=20. && myTree.Flag_HBHENoiseFilter>0 && myTree.Flag_HBHENoiseIsoFilter>0 && myTree.Flag_globalTightHalo2016Filter>0 && myTree.Flag_EcalDeadCellTriggerPrimitiveFilter>0 && myTree.Flag_goodVertices>0 && myTree.Flag_eeBadScFilter>0 && myTree.Flag_badMuonFilter>0 && myTree.Flag_badChargedHadronFilter>0 && (myTree.HLT_MuX_Ele12 || myTree.HLT_Mu8_EleX || myTree.HLT_Mu33_Ele33_NonIso || myTree.HLT_Mu30_Ele30_NonIso) && (abs(myTree.lep_pdgId[0]) != abs(myTree.lep_pdgId[1])) && ((abs(myTree.lep_pdgId[0])==11 && myTree.lep_tightId[0]> 0.5) || (abs(myTree.lep_pdgId[0])==13)) && ((abs(myTree.lep_pdgId[1])==11 && myTree.lep_tightId[1]> 0.5) || (abs(myTree.lep_pdgId[1])==13)) && myTree.jet_pt[0]<13000. && myTree.met_pt/myTree.met_caloPt < 5. );
      //      else continue;
      //if( myTree.zll_ht>=250. && ((myTree.nJet30>1 && myTree.zll_mt2>=200.) || (myTree.zll_met_pt>250.)) );
      //else continue;


      MT2EstimateTree* thisTree_of = anaTree_of->get( myTree.zll_ht, njets, nbjets, minMTBmet, mt2 );
      if(thisTree_of==0) continue;

      int nJetHF30_ = 0;
      for(int j=0; j<myTree.njet; ++j){
	if( myTree.jet_pt[j] < 30. || fabs(myTree.jet_eta[j]) < 3.0 ) continue;
	else ++nJetHF30_;
      }


      if( !myTree.isData){
	  weight = weight*(weight_lep0);
      }

      thisTree_of->assignVar("ID", sample.id );

      //      thisTree_of->assignVar("Z_pt", z.Perp() );
      thisTree_of->assignVar("Z_pt", myTree.zll_pt );
      thisTree_of->assignVar("Z_phi", z.Phi() );
      thisTree_of->assignVar("Z_eta", z.Eta() );
      //      thisTree_of->assignVar("Z_mass", z.M() );
      thisTree_of->assignVar("Z_mass", myTree.zll_mass );
      thisTree_of->assignVar("Z_lepId", abs(myTree.lep_pdgId[0]) );

      thisTree_of->assignVar("nLep", myTree.nlep );
      thisTree_of->assignVar("lep_pdgId0", myTree.lep_pdgId[0] );
      thisTree_of->assignVar("lep_pdgId1", myTree.lep_pdgId[1] );
      thisTree_of->assignVar("lep_pt0", myTree.lep_pt[0] );
      thisTree_of->assignVar("lep_pt1", myTree.lep_pt[1] );
      thisTree_of->assignVar("lep_eta0", myTree.lep_eta[0] );
      thisTree_of->assignVar("lep_eta1", myTree.lep_eta[1] );
      thisTree_of->assignVar("lep_phi0", myTree.lep_phi[0] );
      thisTree_of->assignVar("lep_phi1", myTree.lep_phi[1] );
      thisTree_of->assignVar("raw_mt2", myTree.mt2 );
      thisTree_of->assignVar("raw_met", myTree.met_pt );

      thisTree_of->assignVar("weight_lep0", weight_lep0);
      thisTree_of->assignVar("weight_lep1", weight_lep1);
      thisTree_of->assignVar("weight_lep_err", weight_lep_err);

      thisTree_of->assignVar( "nJetHF30", nJetHF30_ );
      thisTree_of->assignVar( "jet1_pt", myTree.jet1_pt );

      thisTree_of->assignVar("lep_tightId0", myTree.lep_tightId[0] );
      thisTree_of->assignVar("lep_tightId1", myTree.lep_tightId[1] );



//      thisTree_of->assignVar("ID", sample.id );
//      //      thisTree_of->assignVar("Z_pt", z.Perp() );
//      thisTree_of->assignVar("Z_pt", myTree.zll_pt );
//      thisTree_of->assignVar("Z_phi", z.Phi() );
//      //      thisTree_of->assignVar("Z_mass", z.M() );
//      thisTree_of->assignVar("Z_mass", myTree.zll_mass );
//      thisTree_of->assignVar("Z_lepId", abs(myTree.lep_pdgId[0])  );
//
//      thisTree_of->assignVar("nLep", myTree.nlep );
//      thisTree_of->assignVar("lep_pt0", myTree.lep_pt[0] );
//      thisTree_of->assignVar("lep_pt1", myTree.lep_pt[1] );
//      thisTree_of->assignVar("lep_eta0", myTree.lep_eta[0] );
//      thisTree_of->assignVar("lep_eta1", myTree.lep_eta[1] );
//      thisTree_of->assignVar("lep_phi0", myTree.lep_phi[0] );
//      thisTree_of->assignVar("lep_phi1", myTree.lep_phi[1] );
//      thisTree_of->assignVar("raw_mt2", myTree.mt2 );
//      thisTree_of->assignVar("raw_met", myTree.met_pt );
//
//      thisTree_of->assignVar("weight_lep0", weight_lep0);
//      thisTree_of->assignVar("weight_lep1", weight_lep1);
//      thisTree_of->assignVar("weight_lep_err", weight_lep_err);
//
//      thisTree_of->assignVar( "nJetHF30",  nJetHF30_ );
//      thisTree_of->assignVar( "jet1_pt",  myTree.jet1_pt );
//  
//      thisTree_of->assignVar("lep_tightId0", myTree.lep_tightId[0] );
//      thisTree_of->assignVar("lep_tightId1", myTree.lep_tightId[1] );

 
      thisTree_of->fillTree_zll(myTree, weight );
      thisTree_of->yield->Fill(mt2, weight );

    } else
      continue;
 
  } // for entries

  anaTree->finalize();
  anaTree_of->finalize();
  
  delete tree;

  file->Close();
  delete file;
   
}
