#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TVector2.h"

#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"

#define mt2_cxx
#include "interface/mt2.h"

int round(float d) {
  return (int)(floor(d + 0.5));
}

void computeYield( const MT2Sample& sample, const MT2Config& cfg, MT2Analysis<MT2EstimateTree>* anaTree );
void roundLikeData( MT2Analysis<MT2EstimateTree>* data );

float DeltaR(float eta1, float eta2, float phi1, float phi2);
float DeltaPhi(float phi1, float phi2);



int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeLostLepton_CR             |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc < 2 ) {
    std::cout << "USAGE: ./computeLostLepton_CR [configFileName] [data/mc]" << std::endl;
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
      std::cout << "-> You passed a second argument that isn't 'data', nor 'MC', nor 'signal', so I don't know what to do about it." << std::endl;
    }
  }

  std::string outputdir = cfg.getEventYieldDir() + "/llepControlRegion";
  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string regionsSet = cfg.regionsSet();
  std::cout << "Using region set: " << regionsSet << std::endl;

  if( cfg.useMC() && !onlyData ) { // use MC BG estimates 
    
    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;
    
    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 300, 599); // only top (tt, t, ttW, ttZ) and W+jets
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(1209);
    }
    
    
    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this
    
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

  if( !(cfg.dummyAnalysis()) && cfg.dataSamples()!="" && !onlyMC ) {

    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    //    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, 1, 3 );
    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "noDuplicates" );
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

    //    if( myTree.isData && !myTree.isGolden ) continue;
    if( myTree.isData ) {

      if( !myTree.passFilters() ) continue;

    }

    if( !myTree.passBaseline() ) continue;
    if( myTree.nLepLowMT==1 ) ; // For lost lepton CR
    else continue;
    
    if ( myTree.nJet30==1 && !myTree.passMonoJetId(0) ) continue;

    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20csv; 
    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = (njets>1) ? myTree.mt2 : ht;
    float minMTBmet = myTree.minMTBMet;
    
    myTree.nBJet20csv=nbjets;
    //    if( myTree.isData && myTree.run>275125.) continue;
    
    int nMuons10 = myTree.nMuons10;
    int nElectrons10 = myTree.nElectrons10;
    int nPFLep5LowMT = myTree.nPFLep5LowMT;
    int nPFHad10LowMT = myTree.nPFHad10LowMT;
    
    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi();
    //Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi();
    if(!myTree.isData) {
      weight *= myTree.weight_btagsf;
      weight *= myTree.weight_lepsf;
    }

    if (myTree.isData) {

      int id = sample.id;
      // sample IDs for data:
      // JetHT = 1
      // HTMHT = 2
      // MET   = 3

//      if( njets==1 ) {
//
//        if( !( id==3 && myTree.HLT_PFMETNoMu90_PFMHTNoMu90) ) continue;
//
//      } else { // njets>=2
//
//	if( ht>1000. ) {
//          if( !( id==1 && myTree.HLT_PFHT800) ) continue;
//	} else if( ht>575. ) {
//          if( !( id==2 && myTree.HLT_PFHT350_PFMET100 )  ) continue;
//        } else if( ht>450. ) {
//          if( !( id==2 && myTree.HLT_PFHT350_PFMET100 )  ) continue;
//	} else if( ht>200. ) {
//          if( !( id==3 && myTree.HLT_PFMETNoMu90_PFMHTNoMu90  )  ) continue;
//	}
//
//      }

      if ( !(myTree.HLT_PFMET100_PFMHT100 || myTree.HLT_PFHT800 || myTree.HLT_PFHT300_PFMET100) ) continue;

      //OOOOLD if( !(myTree.HLT_PFMETNoMu90_PFMHTNoMu90 || myTree.HLT_PFHT350_PFMET100 || myTree.HLT_PFHT800) ) continue;

    } // if is data


    MT2EstimateTree* thisEstimate;

    if( regionsSet=="zurich" || regionsSet=="zurichPlus" ){ // To avoid signal contamination in 7j 2b and 7j 3b
      
      if( njets>=7 && nbjets>2 ) continue;
      
      else if( njets<7 || nbjets<1) {
	
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
