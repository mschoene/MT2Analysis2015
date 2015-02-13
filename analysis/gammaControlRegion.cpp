#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"


#define mt2_cxx
#include "../interface/mt2_float.h"


#include "TLorentzVector.h"
#include "TH1F.h"
#include "TRandom3.h"


float lumi = 5.; //fb-1



MT2Analysis<MT2EstimateZinvGamma> computeYield( const MT2Sample& sample, const std::string& regionsSet, int selectID=-1 );
void randomizePoisson( MT2Analysis<MT2EstimateZinvGamma>* data );
void randomizeSingleHisto( TRandom3 rand, TH1D* histo );




int main( int argc, char* argv[] ) {


  std::string samplesFileName = "CSA14_Zinv";
  if( argc>1 ) {
    std::string samplesFileName_tmp(argv[1]); 
    samplesFileName = samplesFileName_tmp;
  }

  std::string samplesFile = "../samples/samples_" + samplesFileName + ".dat";
  
  std::vector<MT2Sample> samples_gammaJet = MT2Sample::loadSamples(samplesFile, "GJets");
  if( samples_gammaJet.size()==0 ) {
    std::cout << "There must be an error: didn't find any gamma+jet files in " << samplesFile << "!" << std::endl;
    exit(1209);
  }

  std::cout << std::endl << std::endl;
  std::cout << "-> Loading QCD samples" << std::endl;

  std::vector<MT2Sample> samples_qcd = MT2Sample::loadSamples(samplesFile, "QCD");
  


  std::string regionsSet = "13TeV_CSA14";
  //std::string regionsSet = "13TeV_onlyHT";
  //std::string regionsSet = "13TeV_inclusive";
  //std::string regionsSet = "13TeV_ZinvGammaPurity";

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string outputdir = "GammaControlRegion_" + samplesFileName + "_" + regionsSet;
  system(Form("mkdir -p %s", outputdir.c_str()));

  

  MT2Analysis<MT2EstimateZinvGamma>* prompt = new MT2Analysis<MT2EstimateZinvGamma>( "prompt", regionsSet );
  //MT2Analysis<MT2EstimateZinvGamma>* fake = new MT2Analysis<MT2EstimateZinvGamma>( "fake", regionsSet );

  MT2Analysis<MT2EstimateZinvGamma>* gammaJet = new MT2Analysis<MT2EstimateZinvGamma>( "gammaJet", regionsSet );
  for( unsigned i=0; i<samples_gammaJet.size(); ++i ) {
    (*gammaJet) += (computeYield( samples_gammaJet[i], regionsSet ));
    (*prompt) += (computeYield( samples_gammaJet[i], regionsSet, 22 ));
    //(*fake) += (computeYield( samples_gammaJet[i], regionsSet, 0 ));
  }

  MT2Analysis<MT2EstimateZinvGamma>* qcd = new MT2Analysis<MT2EstimateZinvGamma>( "qcd", regionsSet );
  for( unsigned i=0; i<samples_qcd.size(); ++i ) {
    (*qcd) += (computeYield( samples_qcd[i], regionsSet ));
    (*prompt) += (computeYield( samples_qcd[i], regionsSet, 22 ));
    //(*fake) += (computeYield( samples_qcd[i], regionsSet, 0 ));
  }

  MT2Analysis<MT2EstimateZinvGamma>* gammaCR = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR", regionsSet );
  (*gammaCR) = (*gammaJet) + (*qcd);


  MT2Analysis<MT2Estimate>* purity = new MT2Analysis<MT2Estimate>( "purityMC", regionsSet );
  //(*purity) = *( (MT2Analysis<MT2Estimate>*) (gammaJet) );
  (*purity) = *( (MT2Analysis<MT2Estimate>*) (prompt) );
  (*purity) /= *( (MT2Analysis<MT2Estimate>*) (gammaCR) );


  gammaCR->writeToFile( outputdir + "/mc.root" );
  qcd->addToFile( outputdir + "/mc.root" );
  gammaJet->addToFile( outputdir + "/mc.root" );
  purity->writeToFile( outputdir + "/purityMC.root" );

  // emulate data:
  randomizePoisson(gammaCR);
  gammaCR->writeToFile( outputdir + "/data.root" );


  return 0;

}









MT2Analysis<MT2EstimateZinvGamma> computeYield( const MT2Sample& sample, const std::string& regionsSet, int selectID ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;



  MT2Analysis<MT2EstimateZinvGamma> analysis( sample.sname, regionsSet, sample.id );

  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  int nentries = tree->GetEntries();




  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( myTree.gamma_ht>1000. && sample.id==204 ) continue; // remove high-weight spikes (remove GJet_400to600 leaking into HT>1000)

    if( myTree.gamma_mt2 < 200.) continue;

    if( myTree.nMuons10 > 0) continue;
    if( myTree.nElectrons10 > 0 ) continue;
    if( myTree.nPFLep5LowMT > 0) continue;
    if( myTree.nPFHad10LowMT > 0) continue;

    if( myTree.gamma_deltaPhiMin<0.3 ) continue;
    if( myTree.gamma_diffMetMht>0.5*myTree.gamma_met_pt ) continue;
  
    if( myTree.nVert==0 ) continue;

    if( myTree.gamma_nJet40<2 ) continue;

    if( myTree.ngamma==0 ) continue;

    if( myTree.gamma_idCutBased[0]==0 ) continue;
    //if( myTree.gamma_chHadIso[0]+myTree.gamma_neuHadIso[0] > 10. ) continue;

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

    int mcMatchId = myTree.gamma_mcMatchId[0];
    if( selectID>0 && mcMatchId==0 ) continue;
    //if( selectID>=0 && selectID!=mcMatchId ) continue;

    Double_t weight = myTree.evt_scale1fb*lumi; 

    MT2EstimateZinvGamma* thisEstimate = analysis.get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
    if( thisEstimate==0 ) continue;

    thisEstimate->yield->Fill(myTree.gamma_mt2, weight );

    float iso = myTree.gamma_chHadIso[0]/myTree.gamma_pt[0];

    thisEstimate->fillIso( iso, weight, myTree.gamma_mt2 );

    
  } // for entries


  analysis.finalize();
  

  delete tree;


  file->Close();
  delete file;
  
  return analysis;

}



void randomizePoisson( MT2Analysis<MT2EstimateZinvGamma>* data ) {

  TRandom3 rand(13);


  std::set<MT2Region> regions = data->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Region thisRegion( (*iR) );

    randomizeSingleHisto(rand, data->get(thisRegion)->yield);
    randomizeSingleHisto(rand, data->get(thisRegion)->iso);

    for( unsigned i=0; i < data->get(thisRegion)->iso_bins_hist.size(); ++i ) {
      randomizeSingleHisto(rand, data->get(thisRegion)->iso_bins_hist[i]);
    }

    data->get( thisRegion)->fakeDatasetsFromHistos();

  }// for regions



}



void randomizeSingleHisto( TRandom3 rand, TH1D* histo ) {

  for( unsigned ibin=1; ibin<histo->GetXaxis()->GetNbins()+1; ++ibin ) {

    int poisson_data = rand.Poisson(histo->GetBinContent(ibin));
    histo->SetBinContent(ibin, poisson_data);
    histo->SetBinError(ibin, 0.);

  }  // for bins

}
