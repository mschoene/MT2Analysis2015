#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

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



int round(float d) {
  return (int)(floor(d + 0.5));
}

void computeYield( const MT2Sample& sample, const MT2Config& cfg, MT2Analysis<MT2EstimateTree>* anaTree );
void roundLikeData( MT2Analysis<MT2EstimateTree>* data );

float DeltaR(float eta1, float eta2, float phi1, float phi2);
float DeltaPhi(float phi1, float phi2);

template <class T>
MT2Analysis<T>* computeSigYield( const MT2Sample& sample, const MT2Config& cfg );



int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running llepControlRegion                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc < 2 ) {
    std::cout << "USAGE: ./llepControlRegion [configFileName] [data/mc]" << std::endl;
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
  
  std::string outputdir = cfg.getEventYieldDir() + "/llepControlRegion";
  system(Form("mkdir -p %s", outputdir.c_str()));

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this  

  std::string regionsSet = cfg.regionsSet();
  std::cout << "Using region set: " << regionsSet << std::endl;

  if( cfg.useMC() && !onlyData && !onlySignal ) { // use MC BG estimates 
    
    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;
    
    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 300, 599); // only top (tt, t, ttW, ttZ) and W+jets
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(1209);
    }
    
    MT2Analysis<MT2EstimateTree>* llepCR = new MT2Analysis<MT2EstimateTree> ( "llepCR", regionsSet );  
        
    for( unsigned i=0; i < fSamples.size(); ++i )
      computeYield( fSamples[i], cfg, llepCR );
    
    llepCR->writeToFile( outputdir + "/mc.root" );
    
    if( cfg.dummyAnalysis() ) {
      
      // emulate data:
      roundLikeData(llepCR);
      
      llepCR->writeToFile( outputdir + "/data.root" );
      
    }
  }


  // load signal samples, if any
  std::vector< MT2Analysis< MT2EstimateAllSigSyst >* > signals;

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

      if( signals.size()>0 ){

	signals[0]->writeToFile( outputdir + "/mc.root" );
	for( unsigned i=1; i<signals.size(); ++i )
	  signals[i]->writeToFile( outputdir + "/mc.root" );
      }
      
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

      if( signals.size()>0 ){
	signals[0]->writeToFile(outputdir + "/mc.root");
	for( unsigned i=1; i<signals.size(); ++i )
	  signals[i]->writeToFile(outputdir + "/mc.root");
      }

    } // if samples != 0
    
  } // if sig samples
  

  if( !(cfg.dummyAnalysis()) && cfg.dataSamples()!="" && !onlyMC && !onlySignal) {

    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "merged" );
    //    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "noDuplicates" );
    if( samples_data.size()==0 ) {
      std::cout << "There must be an error: samples_data is empty!" << std::endl;
      exit(1209);
    }

    
    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this
  
    MT2Analysis<MT2EstimateTree>* dataCR = new MT2Analysis<MT2EstimateTree> ( "llepCR", regionsSet );

    for( unsigned i=0; i < samples_data.size(); ++i )
      computeYield( samples_data[i], cfg, dataCR );

    dataCR->writeToFile( outputdir + "/data.root" );

  }


  return 0;
  
}



void computeYield( const MT2Sample& sample, const MT2Config& cfg, MT2Analysis<MT2EstimateTree>* anaTree ){

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  MT2Tree myTree;
  myTree.Init(tree);

  std::string regionsSet = cfg.regionsSet();
  std::cout << "Using region set: " << regionsSet << std::endl;
  
  int nentries = tree->GetEntries();
    
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( myTree.isData && !myTree.isGolden ) continue;
    if( myTree.isData ) {

      if( !myTree.passFilters() ) continue;

    }

    if( !myTree.passBaseline() ) continue;
    if( myTree.nLepLowMT==1 ) ; // For lost lepton CR
    else continue;
    
    if ( myTree.nJet30==1 && !myTree.passMonoJetId(0) ) continue;

    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20; 
    //    float ht   = myTree.ht;
    float ht   = (njets>1) ? myTree.ht : myTree.jet1_pt;
    //    float met  = myTree.met_pt;
    float mt2  = (njets>1) ? myTree.mt2 : ht;
    float minMTBmet = myTree.minMTBMet;
    
    //    int nMuons10 = myTree.nMuons10;
    //    int nElectrons10 = myTree.nElectrons10;
    //    int nPFLep5LowMT = myTree.nPFLep5LowMT;
    //    int nPFHad10LowMT = myTree.nPFHad10LowMT;
    
    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi();
    // if(!myTree.isData) {
    //   weight *= myTree.weight_btagsf;
    //   weight *= myTree.weight_lepsf;
    // }

    if( !myTree.isData ){
      weight *= myTree.weight_btagsf;
      weight *= myTree.weight_lepsf2017;
      //      weight *= myTree.weight_lepsf;

      // // ETH has a branch witht he average weight stored:
      // // Also we have a different numbering scheme...
      // if (myTree.evt_id == 302 || myTree.evt_id == 303 || myTree.evt_id == 304) //singleLep T/Tbar, Dilep
      // 	weight *= myTree.weight_isr / myTree.weight_isr_norm;

      ///AMERICAN WAY
      if (myTree.evt_id == 301 || myTree.evt_id == 302)
      	weight *= myTree.weight_isr/0.909; // nominal
      else if (myTree.evt_id == 303) 
      	weight *= myTree.weight_isr/0.895;

    }

    if (myTree.isData) {

      if ( !(myTree.HLT_PFMET120_PFMHT120 || myTree.HLT_PFHT900 || myTree.HLT_PFHT300_PFMET110 || myTree.HLT_PFJet450  || myTree.HLT_PFMETNoMu120_PFMHTNoMu120 ) ) continue;
      //      if ( !(myTree.HLT_PFMET100_PFMHT100 || myTree.HLT_PFHT800 || myTree.HLT_PFHT300_PFMET100) ) continue; //ICHEP 2016
      //OLD if( !(myTree.HLT_PFMETNoMu90_PFMHTNoMu90 || myTree.HLT_PFHT350_PFMET100 || myTree.HLT_PFHT800) ) continue;

    } // if is data


    if( myTree.nJet200MuFrac50DphiMet > 0 ) continue; // new RA2 filter

    if( myTree.met_miniaodPt/myTree.met_caloPt > 5.0 ) continue;

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

    MT2EstimateTree* thisEstimate;

    if( regionsSet=="zurich" || regionsSet=="zurichPlus" || regionsSet=="zurich2016" ){ // To avoid signal contamination in 7j 2b and 7j 3b
    //    if( false ){ // To avoid signal contamination in 7j 2b and 7j 3b
      
      if( ht>450. && njets>=7 && nbjets>2 ) continue;
      
      else if( ht<450 || njets<7 || nbjets<1 ) {
	
	thisEstimate = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
	if( thisEstimate==0 ) continue;
	
	thisEstimate->assignTree( myTree, weight );
	thisEstimate->tree->Fill();
	thisEstimate->yield->Fill(mt2, weight );
	
      }
      else {
	
	thisEstimate = anaTree->get( ht, njets, 1, minMTBmet, mt2 );
	if( thisEstimate==0 ) continue;
	thisEstimate->assignTree( myTree, weight );
	thisEstimate->tree->Fill();
	thisEstimate->yield->Fill(mt2, weight );
	
	thisEstimate = anaTree->get( ht, njets, 2, minMTBmet, mt2 );
	if( thisEstimate==0 ) continue;
	thisEstimate->assignTree( myTree, weight );
	thisEstimate->tree->Fill();
	thisEstimate->yield->Fill(mt2, weight );
	
	thisEstimate = anaTree->get( ht, njets, 3, minMTBmet, mt2 );
	if( thisEstimate==0 ) continue;
	thisEstimate->assignTree( myTree, weight );
	thisEstimate->tree->Fill();
	thisEstimate->yield->Fill(mt2, weight );
	
      }
    
    }
    else {

	thisEstimate = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
	if( thisEstimate==0 ) continue;
	
	thisEstimate->assignTree( myTree, weight );
	thisEstimate->tree->Fill();
	thisEstimate->yield->Fill(mt2, weight );

    }
    

  } // for entries

  anaTree->finalize();

  delete tree;

  file->Close();
  delete file;

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
    
    bool passGenMET = false;
    if(dogenmet)
      passGenMET = true;

    bool passRecoMET=true;

    if( !myTree.passBaseline() ) passRecoMET=false;

    if(dogenmet)
      if( !myTree.passBaseline("genmet") ) passGenMET=false;

    if (!passGenMET && !passRecoMET) continue;
 
    if( myTree.nLepLowMT==1 ); 
    else continue;
    
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

      }

    }
 
    //Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi()*myTree.puWeight;
    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi();
    //Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi();
    Double_t weight_syst = 1.;

    /////    weight = 1000.* myTree.evt_xsec/nentries; //Exceptionally for signal from muricans 

    if( !myTree.isData ){
      weight *= myTree.weight_btagsf;
      weight *= myTree.weight_lepsf2017;
   
      if (myTree.evt_id == 301 || myTree.evt_id == 302)
	weight *= myTree.weight_isr/0.909; // central/average
      else if (myTree.evt_id == 303) 
	weight *= myTree.weight_isr/0.895;  
    }

    float sig_xs=0.;
    if( myTree.evt_id >= 1000  && myTree.evt_id < 2000){

      int thisBinX = sigXS->FindBin( GenSusyMScan1 );

      sig_xs = sigXS->GetBinContent(thisBinX);

      weight *= sig_xs;

    }

    T* thisEstimate = analysis->get( ht, njets, nbjets, minMTBmet, 201. );
    if( thisEstimate==0 ) continue;

    if( regionsSet=="zurich" || regionsSet=="zurichPlus" || regionsSet=="zurich2016" ){ // To avoid signal contamination in 7j 2b and 7j 3b                                                                                

      if( ht>450 && njets>=7 && nbjets>2 ) continue;

      else if( ht<450 || njets<7 || nbjets<1) {
	
	if(passRecoMET){

	  thisEstimate  = analysis->get( ht, njets, nbjets, minMTBmet, mt2 );
	  if( thisEstimate==0 ) continue;
		  
	  thisEstimate->yield->Fill( mt2, weight );
	  thisEstimate->yield3d->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight );
	  
	}
	
	if(dogenmet && passGenMET){

	  thisEstimate  = analysis->get( ht, njets, nbjets, minMTBmet, mt2_genmet );
	  if( thisEstimate==0 ) continue;
		  	  
	  thisEstimate->yield3d_genmet->Fill( mt2_genmet, GenSusyMScan1, GenSusyMScan2, weight );
	  
	}
    
      }
      else {

        if(passRecoMET){

	  thisEstimate  = analysis->get( ht, njets, 1, minMTBmet, mt2 );
	  if( thisEstimate==0 ) continue;
	
          thisEstimate->yield->Fill( mt2, weight );
          thisEstimate->yield3d->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight );

        }
	
        if(dogenmet && passGenMET){

	  thisEstimate  = analysis->get( ht, njets, 1, minMTBmet, mt2_genmet );
	  if( thisEstimate==0 ) continue;
	
          thisEstimate->yield3d_genmet->Fill( mt2_genmet, GenSusyMScan1, GenSusyMScan2, weight );

        }

        if(passRecoMET){

	  thisEstimate  = analysis->get( ht, njets, 2, minMTBmet, mt2 );
	  if( thisEstimate==0 ) continue;

          thisEstimate->yield->Fill( mt2, weight );
          thisEstimate->yield3d->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight );

        }

        if(dogenmet && passGenMET){

	  thisEstimate  = analysis->get( ht, njets, 2, minMTBmet, mt2_genmet );
	  if( thisEstimate==0 ) continue;

          thisEstimate->yield3d_genmet->Fill( mt2_genmet, GenSusyMScan1, GenSusyMScan2, weight );

        }
	
        if(passRecoMET){

	  thisEstimate  = analysis->get( ht, njets, 3, minMTBmet, mt2 );
	  if( thisEstimate==0 ) continue;

          thisEstimate->yield->Fill( mt2, weight );
          thisEstimate->yield3d->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight );

        }

        if(dogenmet && passGenMET){

	  thisEstimate  = analysis->get( ht, njets, 3, minMTBmet, mt2_genmet );
	  if( thisEstimate==0 ) continue;

          thisEstimate->yield3d_genmet->Fill( mt2_genmet, GenSusyMScan1, GenSusyMScan2, weight );

        }

	
      }

    }
    else {

      
      if(passRecoMET){

	thisEstimate  = analysis->get( ht, njets, nbjets, minMTBmet, mt2 );
	if( thisEstimate==0 ) continue;
      
	thisEstimate->yield->Fill( mt2, weight );
	thisEstimate->yield3d->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight );

      }

      if(dogenmet && passGenMET){

	thisEstimate  = analysis->get( ht, njets, nbjets, minMTBmet, mt2_genmet );
	if( thisEstimate==0 ) continue;
      
	thisEstimate->yield3d_genmet->Fill( mt2_genmet, GenSusyMScan1, GenSusyMScan2, weight );

      }
      
    }

  } // for entries

  analysis->finalize();

  delete tree;

  file->Close();
  delete file;

  return analysis;

}




void roundLikeData( MT2Analysis<MT2EstimateTree>* data ) {
  
   std::set<MT2Region> regions = data->getRegions();

   for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

     TH1D* thisYield = data->get(*iR)->yield;

     for( int iBin=1; iBin<thisYield->GetNbinsX()+1; ++iBin ) {

       float yield = thisYield->GetBinContent(iBin);
       thisYield->SetBinContent(iBin, round( yield ));
       thisYield->SetBinError(iBin, 0. );
       
     } // for bins
     
   } // for regions

}


float DeltaR(float eta1, float eta2, float phi1, float phi2){
  float dEta = eta1 - eta2;
  float dPhi = DeltaPhi(phi1, phi2);
  return TMath::Sqrt(dEta*dEta + dPhi*dPhi);
}

float DeltaPhi(float phi1, float phi2){
  float dPhi = phi1 - phi2;
  while (dPhi  >  TMath::Pi()) dPhi -= 2*TMath::Pi();
  while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
  return fabs(dPhi);
}
