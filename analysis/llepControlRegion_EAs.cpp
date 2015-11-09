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

//  std::string regionsSet;
//  if( cfg.regionsSet()=="zurich" )
//    regionsSet="zurich_llep";
//  else if( cfg.regionsSet()=="zurichPlus" )
//    regionsSet="zurichPlus_llep";
//  else
//    regionsSet=cfg.regionsSet();
  
  std::string regionsSet = cfg.regionsSet();
  std::cout << "Using region set: " << regionsSet << std::endl;

  if( cfg.useMC() && !onlyData ) { // use MC BG estimates 
    
    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;
    
    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 100, 999); // only top (tt, t, ttW, ttZ) and W+jets
    //    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 300, 499); // only top (tt, t, ttW, ttZ)
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(1209);
    }
    
    
    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this
    
    MT2Analysis<MT2EstimateTree>* llepCR = new MT2Analysis<MT2EstimateTree> ( "llepCR", regionsSet );  
    MT2EstimateTree::addVar( llepCR, "lep_pt" );
        
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

    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, 1, 3 );
    if( samples_data.size()==0 ) {
      std::cout << "There must be an error: samples_data is empty!" << std::endl;
      exit(1209);
    }

    
    TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this
  
    MT2Analysis<MT2EstimateTree>* dataCR = new MT2Analysis<MT2EstimateTree> ( "llepCR", regionsSet );
    MT2EstimateTree::addVar( dataCR, "lep_pt" );

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

//  std::string regionsSet;
//  if( cfg.regionsSet()=="zurich" )
//    regionsSet="zurich_llep";
//  else if( cfg.regionsSet()=="zurichPlus" )
//    regionsSet="zurichPlus_llep";
//  else
//    regionsSet=cfg.regionsSet();

  std::string regionsSet = cfg.regionsSet();
  std::cout << "Using region set: " << regionsSet << std::endl;
  
  int nentries = tree->GetEntries();
    
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 500 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( myTree.isData && !myTree.isGolden ) continue;
    if( myTree.isData ) {

      if( !( myTree.Flag_HBHENoiseFilter && myTree.Flag_CSCTightHaloFilter &&  myTree.Flag_eeBadScFilter ) ) continue;

    }

    if( myTree.nVert<1 ) continue;
    if( myTree.nJet40<2 ) continue;
    if( myTree.met_pt<250 ) continue;
    if( myTree.ht<500 ) continue;

    if( myTree.nElectrons10>0 ) continue;
    if( myTree.nMuons10<=0 ) continue;
    if( !(myTree.nlep==1) ) continue;
    
  
    int nMuons25LowMT=0;

    if( fabs(myTree.lep_pdgId[0])!=13 || myTree.lep_pt[0]<25. || fabs(myTree.lep_eta[0])>2.5 ) continue;
    
    float mt=0;

    float met  = myTree.met_pt;
    float met_phi = myTree.met_phi;
    
    float lepPt=myTree.lep_pt[0];
    float lepPhi=myTree.lep_phi[0];
    
    float dPhi =DeltaPhi(lepPhi, met_phi);
    mt = TMath::Sqrt( 2*lepPt*met*( 1 - TMath::Cos(dPhi) ) );
    if( mt < 120. ) 
      nMuons25LowMT=1;

//    for(int l=0; l<myTree.nlep; ++l){
//      
//      if( abs(myTree.lep_pdgId[l])!=13 || myTree.lep_pt[l]<25. || fabs(myTree.lep_eta[l])>2.5 ) continue;
//      
//      float mt=0;
//      float lepPt=myTree.lep_pt[l];
//      float lepPhi=myTree.lep_phi[l];
//      
//      float dPhi =DeltaPhi(lepPhi, met_phi);
//      
//      mt = TMath::Sqrt( 2*lepPt*met*( 1-TMath::Cos(dPhi) ) );
//      if( mt > 120. ) continue;
//      else {
//	++nMuons25LowMT;
//	lepIndex=l;
//      }
//
//    }
    
    if( nMuons25LowMT < 1 ) continue;
    float thisLep_pt = lepPt;

//    float ht=0.;
//    for(int j=0; j<myTree.njet; ++j){
//
//      if(myTree.jet_pt[j]<40. || fabs(myTree.jet_eta[j])>2.5) continue;
//      ht+=myTree.jet_pt[j];
//
//    }

    int njets  = myTree.nJet40;
    int nbjets = myTree.nBJet40; 
    float ht   = myTree.ht;
    float mt2  = 201.;
    float minMTBmet = myTree.minMTBMet;
        
    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi();

    if (myTree.isData) {

      int id = sample.id;
      // sample IDs for data:
      // JetHT = 1
      // HTMHT = 2
      // MET   = 3

      if( njets==1 ) {

        if( !( id==3 && myTree.HLT_PFMET90_PFMHT90) ) continue;

      } else { // njets>=2

	if( ht>1000. ) {
          if( !( id==1 && myTree.HLT_PFHT800) ) continue;
	} else if( ht>575. ) {
          if( !( id==2 && myTree.HLT_PFHT350_PFMET100 )  ) continue;
        } else if( ht>450. ) {
          if( !( id==2 && myTree.HLT_PFHT350_PFMET100 )  ) continue;
	} else if( ht>200. ) {
          if( !( id==3 && myTree.HLT_PFMET90_PFMHT90  )  ) continue;
	}

      }

    } // if is data


    MT2EstimateTree* thisEstimate;
    thisEstimate = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisEstimate==0 ) continue;
    
    thisEstimate->yield->Fill( mt2, weight );

    thisEstimate->assignVar( "lep_pt", thisLep_pt );
    thisEstimate->assignTree( myTree, weight );
    
    thisEstimate->tree->Fill();
    

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
