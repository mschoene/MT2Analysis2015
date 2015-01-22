#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>



#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateSyst.h"

#define mt2_cxx
#include "../interface/mt2_float.h"


#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLorentzVector.h"




float lumi = 5.; // fb-1



MT2Analysis<MT2EstimateSyst> computeYield( const MT2Sample& sample, const std::string& regionsSet, const std::string& prefix="" );
void addPoissonError( MT2Analysis<MT2EstimateSyst>* analysis );
MT2Analysis<MT2EstimateSyst>* combineDataAndMC( MT2Analysis<MT2EstimateSyst>* data, MT2Analysis<MT2EstimateSyst>* mc );


int main( int argc, char* argv[] ) {


  //if( argc!=2 ) {
  //  std::cout << "USAGE: ./computeZinvFromGamma [samplesFileName]" << std::endl;
  //  std::cout << "Exiting." << std::endl;
  //  exit(11);
  //}


  std::string samplesFileName = "CSA14_Zinv";
  //std::string samplesFileName = "CSA14_skimprune_Zinv";
  if( argc>1 ) {
    std::string samplesFileName_tmp(argv[1]); 
    samplesFileName = samplesFileName_tmp;
  }

  std::string samplesFile = "../samples/samples_" + samplesFileName + ".dat";

  std::cout << std::endl << std::endl;
  std::cout << "-> Loading gamma+jet samples" << std::endl;

  std::vector<MT2Sample> samples_gammaJet = MT2Sample::loadSamples(samplesFile, "GJets");
  if( samples_gammaJet.size()==0 ) {
    std::cout << "There must be an error: didn't find any gamma+jet files in " << samplesFile << "!" << std::endl;
    exit(1209);
  }

  std::cout << std::endl << std::endl;
  std::cout << "-> Loading QCD samples" << std::endl;

  std::vector<MT2Sample> samples_qcd = MT2Sample::loadSamples(samplesFile, "QCD");
  if( samples_qcd.size()==0 ) {
    std::cout << "There must be an error: didn't find any QCD files in " << samplesFile << "!" << std::endl;
    exit(1205);
  }
  

  std::cout << std::endl << std::endl;
  std::cout << "-> Loading Zinv samples" << std::endl;

  std::vector<MT2Sample> samples_Zinv = MT2Sample::loadSamples(samplesFile, "ZJetsToNuNu");
  if( samples_Zinv.size()==0 ) {
    std::cout << "There must be an error: didn't find any Zinv files in " << samplesFile << "!" << std::endl;
    exit(1207);
  }


  //std::string regionsSet = "13TeV_inclusive";
  std::string regionsSet = "13TeV_CSA14";

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string outputdir = "ZinvEstimateFromGamma_" + samplesFileName + "_" + regionsSet;
  system(Form("mkdir -p %s", outputdir.c_str()));

  
  MT2Analysis<MT2EstimateSyst>* gammaJet = new MT2Analysis<MT2EstimateSyst>( "gammaJet", regionsSet );
  for( unsigned i=0; i<samples_gammaJet.size(); ++i ) {
    (*gammaJet) += (computeYield( samples_gammaJet[i], regionsSet, "gamma_" ));
  }

  
  MT2Analysis<MT2EstimateSyst>* qcd = new MT2Analysis<MT2EstimateSyst>( "qcd", regionsSet );
  for( unsigned i=0; i<samples_qcd.size(); ++i ) {
    (*qcd) += (computeYield( samples_qcd[i], regionsSet, "gamma_" ));
  }

  MT2Analysis<MT2EstimateSyst>* gamma_plus_qcd = new MT2Analysis<MT2EstimateSyst>( "gamma_plus_qcd", regionsSet );
  *gamma_plus_qcd = *gammaJet;
  *gamma_plus_qcd += *qcd;
  
  MT2Analysis<MT2EstimateSyst>* purity = new MT2Analysis<MT2EstimateSyst>( "purity", regionsSet );
  (*purity) = (*gammaJet) / (*gamma_plus_qcd);
  purity->setName( "purity" );

  MT2Analysis<MT2EstimateSyst>* Zinv = new MT2Analysis<MT2EstimateSyst>( "Zinv", regionsSet );
  for( unsigned i=0; i<samples_Zinv.size(); ++i ) {
    (*Zinv) += (computeYield( samples_Zinv[i], regionsSet ));
  }
  


  MT2Analysis<MT2EstimateSyst>* ZgammaRatio = new MT2Analysis<MT2EstimateSyst>( "ZgammaRatio", regionsSet );
  (*ZgammaRatio) = (*Zinv) / (*gammaJet);

  // now that the MC ratio is done, add poisson error to gammajet sample:
  addPoissonError(gammaJet);


  MT2Analysis<MT2EstimateSyst>* ZinvEstimateFromGamma = new MT2Analysis<MT2EstimateSyst>( "ZinvEstimateFromGamma", regionsSet );
  (*ZinvEstimateFromGamma) = (*ZgammaRatio) * (*gammaJet);


  MT2Analysis<MT2EstimateSyst>* ZJets = MT2Analysis<MT2EstimateSyst>::readFromFile( "EventYields_mc_PHYS14_dummy_5fb/analyses.root", "ZJets" );

  MT2Analysis<MT2EstimateSyst>* ZinvEstimate = combineDataAndMC( ZinvEstimateFromGamma, ZJets );
  //MT2Analysis<MT2EstimateSyst>* ZinvEstimate = combineDataAndMC( ZinvEstimateFromGamma, Zinv );
  ZinvEstimate->writeToFile( outputdir + "/MT2ZinvEstimate.root" );


  std::string outputdirPlots = outputdir + "/plots";
  system(Form("mkdir -p %s", outputdirPlots.c_str()));


  std::string mcFile = outputdir + "/mc.root";
  gammaJet->writeToFile( mcFile );
  qcd->addToFile( mcFile );
  Zinv->addToFile( mcFile );
  ZgammaRatio->addToFile( mcFile );

  purity->writeToFile( outputdir + "/MT2GammaPurity.root" );


  return 0;

}




MT2Analysis<MT2EstimateSyst> computeYield( const MT2Sample& sample, const std::string& regionsSet, const std::string& prefix ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;



  MT2Analysis<MT2EstimateSyst> analysis( sample.sname, regionsSet, sample.id );

  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  int nentries = tree->GetEntries();

  float mt2;
  tree->SetBranchAddress( Form("%smt2"   , prefix.c_str()), &mt2 );
  float ht;
  tree->SetBranchAddress( Form("%sht"    , prefix.c_str()), &ht );
  float met;
  tree->SetBranchAddress( Form("%smet_pt", prefix.c_str()), &met );
  int njets;
  tree->SetBranchAddress( Form("%snJet40", prefix.c_str()), &njets );
  int nbjets;
  tree->SetBranchAddress( Form("%snBJet40", prefix.c_str()), &nbjets );
  float deltaPhiMin;
  tree->SetBranchAddress( Form("%sdeltaPhiMin", prefix.c_str()), &deltaPhiMin );
  float diffMetMht;
  tree->SetBranchAddress( Form("%sdiffMetMht", prefix.c_str()), &diffMetMht );
  float minMTBMet;
  tree->SetBranchAddress( Form("%sminMTBMet", prefix.c_str()), &minMTBMet );




  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( myTree.nMuons10 > 0) continue;
    if( myTree.nElectrons10 > 0 ) continue;
    if( myTree.nPFLep5LowMT > 0) continue;
    if( myTree.nPFHad10LowMT > 0) continue;

    if( deltaPhiMin<0.3 ) continue;
    if( diffMetMht>0.5*met ) continue;
  
    if( myTree.nVert==0 ) continue;

    if( njets<2 ) continue;

    if( myTree.ngamma>0 && prefix=="gamma_" ) {

      if( myTree.gamma_idCutBased==0 ) continue;
      if( myTree.gamma_chHadIso[0]+myTree.gamma_neuHadIso[0] > 10. ) continue;

      TLorentzVector gamma;
      gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );
      float found_pt = 0.;
      int foundjet = 0;
      for( unsigned i=0; i<myTree.njet; ++i ) {
        TLorentzVector thisjet;
        thisjet.SetPtEtaPhiM( myTree.jet_pt[i], myTree.jet_eta[i], myTree.jet_phi[i], myTree.jet_mass[i] );
        if( gamma.DeltaR(thisjet)>0.4 ) foundjet++;
        if( foundjet==2 ) {
          found_pt = thisjet.Pt();
          break;
        }
      }
      if( found_pt<100. ) continue;

    } else {

      if( myTree.jet_pt[1]<100. ) continue;

    }


    Double_t weight = myTree.evt_scale1fb*lumi; 

    MT2EstimateSyst* thisEstimate = analysis.get( ht, njets, nbjets, met, minMTBMet, mt2 );
    if( thisEstimate==0 ) continue;

    thisEstimate->yield->Fill(mt2, weight );

    
  } // for entries


  analysis.finalize();
  

  delete tree;

  file->Close();
  delete file;
  
  return analysis;

}






void addPoissonError( MT2Analysis<MT2EstimateSyst>* analysis ) {


  std::set<MT2Region> regions = analysis->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

      TH1D* h1 = analysis->get(*iR)->yield;

      for( unsigned ibin=1; ibin<h1->GetXaxis()->GetNbins()+1; ++ibin ) {

        int nData = (int) h1->GetBinContent(ibin);
        h1->SetBinContent(ibin, nData);
        if( nData==0 )
          h1->SetBinError(ibin, 0.);
        else
          h1->SetBinError(ibin, sqrt((float)nData));

      }  // for bins

  }// for regions

}



MT2Analysis<MT2EstimateSyst>* combineDataAndMC( MT2Analysis<MT2EstimateSyst>* data, MT2Analysis<MT2EstimateSyst>* mc ) {

  std::string dataname = data->getName();
  std::string mcname = mc->getName();

  // temporarily set all names to the output name so that returned MT2Analysis has consistent naming in all regions:
  std::string estimateName = "ZinvEstimate";
  data->setName( estimateName );
  mc->setName( estimateName );

  std::set<MT2Region> regions = data->getRegions();

  std::set<MT2EstimateSyst*> newData;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateSyst* dataEst = data->get(*iR);
    MT2EstimateSyst* mcEst = mc->get(*iR);

    MT2EstimateSyst* thisNewEstimate;
    if( iR->nBJetsMin()>1 ) {
      thisNewEstimate =  new MT2EstimateSyst(*mcEst);
      for( unsigned ibin=1; ibin<thisNewEstimate->yield->GetNbinsX()+1; ++ibin )
        thisNewEstimate->yield->SetBinError( ibin, thisNewEstimate->yield->GetBinContent(ibin) );
    } else {
      thisNewEstimate =  new MT2EstimateSyst(*dataEst);
    }
    newData.insert( thisNewEstimate );

  }

  MT2Analysis<MT2EstimateSyst>* analysis = new MT2Analysis<MT2EstimateSyst>( estimateName, newData );

  // set names back to original:
  data->setName( dataname );
  mc->setName( mcname );


  return analysis;

}
