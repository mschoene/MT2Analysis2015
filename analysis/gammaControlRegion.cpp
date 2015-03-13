#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"
#include "interface/MT2EstimateSyst.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"
#include "TRandom3.h"


float lumi = 4.; //fb-1



void computeYield( const MT2Sample& sample, const std::string& regionsSet, MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* prompt_pass, MT2Analysis<MT2EstimateZinvGamma>* nip, MT2Analysis<MT2EstimateZinvGamma>* nip_pass, MT2Analysis<MT2EstimateZinvGamma>* fake, MT2Analysis<MT2EstimateZinvGamma>* fake_pass, float isoCut );
void randomizePoisson( MT2Analysis<MT2EstimateZinvGamma>* data );
void randomizeSingleHisto( TRandom3 rand, TH1D* histo );




int main( int argc, char* argv[] ) {


  std::string samplesFileName = "PHYS14_v2_Zinv";
  //std::string samplesFileName = "CSA14_Zinv";
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
  MT2Analysis<MT2EstimateZinvGamma>* prompt_pass = new MT2Analysis<MT2EstimateZinvGamma>( "prompt_pass", regionsSet );

  MT2Analysis<MT2EstimateZinvGamma>* fake = new MT2Analysis<MT2EstimateZinvGamma>( "fake", regionsSet );
  MT2Analysis<MT2EstimateZinvGamma>* fake_pass = new MT2Analysis<MT2EstimateZinvGamma>( "fake_pass", regionsSet );

  MT2Analysis<MT2EstimateZinvGamma>* nip = new MT2Analysis<MT2EstimateZinvGamma>( "nip", regionsSet );
  MT2Analysis<MT2EstimateZinvGamma>* nip_pass = new MT2Analysis<MT2EstimateZinvGamma>( "nip_pass", regionsSet );


  float isoCut = 2.5; // GeV (absolute iso)

  for( unsigned i=0; i<samples_gammaJet.size(); ++i ) {
    computeYield( samples_gammaJet[i], regionsSet, prompt, prompt_pass, nip, nip_pass, fake, fake_pass, isoCut );
  }

  for( unsigned i=0; i<samples_qcd.size(); ++i ) {
    computeYield( samples_qcd[i], regionsSet, prompt, prompt_pass, nip, nip_pass, fake, fake_pass, isoCut );
  }


  MT2Analysis<MT2EstimateZinvGamma>* gammaCR = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR_loose", regionsSet );
  (*gammaCR) = (*prompt) + (*nip) + (*fake);

 
  MT2Analysis<MT2EstimateZinvGamma>* fake_plus_prompt_pass = new MT2Analysis<MT2EstimateZinvGamma>( (*prompt_pass) + (*fake_pass) ); 
  fake_plus_prompt_pass->setName("fake_plus_prompt_pass");

  MT2Analysis<MT2EstimateZinvGamma>* nip_pass_times2 = new MT2Analysis<MT2EstimateZinvGamma>( *nip );
  nip_pass_times2->setName("nip_pass_times2");
  (*nip_pass_times2) *= 2.;

  MT2Analysis<MT2EstimateZinvGamma>* nip_pass_div2 = new MT2Analysis<MT2EstimateZinvGamma>( *nip );
  nip_pass_div2->setName("nip_pass_div2");
  (*nip_pass_div2) /= 2.;

 
  MT2Analysis<MT2EstimateZinvGamma>* gammaCR_pass         = new MT2Analysis<MT2EstimateZinvGamma>( *fake_plus_prompt_pass + *nip_pass );
  gammaCR_pass->setName("gammaCR");
  MT2Analysis<MT2EstimateZinvGamma>* gammaCR_nipUp_pass   = new MT2Analysis<MT2EstimateZinvGamma>( *fake_plus_prompt_pass + *nip_pass_times2 );
  gammaCR_nipUp_pass->setName("gammaCR_nipUp");
  MT2Analysis<MT2EstimateZinvGamma>* gammaCR_nipDown_pass = new MT2Analysis<MT2EstimateZinvGamma>( *fake_plus_prompt_pass + *nip_pass_div2 );
  gammaCR_nipDown_pass->setName("gammaCR_nipDown");
  
 

  MT2Analysis<MT2EstimateSyst>* eff_isoCut = MT2EstimateSyst::makeEfficiencyAnalysis( "eff_isoCut", regionsSet, (MT2Analysis<MT2Estimate>*)prompt_pass, (MT2Analysis<MT2Estimate>*)prompt );

  MT2Analysis<MT2EstimateSyst>* purityTight = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", regionsSet, (MT2Analysis<MT2Estimate>*)prompt_pass, (MT2Analysis<MT2Estimate>*)gammaCR_pass );

  MT2Analysis<MT2EstimateSyst>* purityLoose = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", regionsSet, (MT2Analysis<MT2Estimate>*)prompt, (MT2Analysis<MT2Estimate>*)gammaCR );


  gammaCR->writeToFile( outputdir + "/mc.root" );
  prompt->addToFile( outputdir + "/mc.root" );
  fake->addToFile( outputdir + "/mc.root" );

  purityTight->writeToFile( outputdir + "/purityMC.root" );
  purityLoose->addToFile( outputdir + "/purityMC.root" );
  eff_isoCut->addToFile( outputdir + "/purityMC.root" );

  // emulate data:
  //randomizePoisson(gammaCR);
  gammaCR_pass->writeToFile( outputdir + "/data.root" );
  gammaCR_nipUp_pass->addToFile( outputdir + "/data.root" );
  gammaCR_nipDown_pass->addToFile( outputdir + "/data.root" );
 

  return 0;

}









void computeYield( const MT2Sample& sample, const std::string& regionsSet, MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* prompt_pass, MT2Analysis<MT2EstimateZinvGamma>* nip, MT2Analysis<MT2EstimateZinvGamma>* nip_pass, MT2Analysis<MT2EstimateZinvGamma>* fake, MT2Analysis<MT2EstimateZinvGamma>* fake_pass, float isoCut ) {

  

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;



  bool isQCD  = sample.id>=100 && sample.id<200;
  bool isGJet = sample.id>=200 && sample.id<300;


  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);


  int nentries = tree->GetEntries();


  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( myTree.met_pt>100. ) continue; // orthogonal to signal regions

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


    int mcMatchId = myTree.gamma_mcMatchId[0];
    bool isMatched = (mcMatchId==22 || mcMatchId==7);
    bool isGenIso = (myTree.gamma_genIso[0]<5.);

    bool isPrompt = ( isMatched &&  isGenIso);
    bool isNIP    = ( isMatched && !isGenIso);
    bool isFake   = (!isMatched);

    if( isPrompt && isQCD  ) continue; //isolated prompts taken from GJet only
    if( isNIP    && isGJet ) continue; //non-isolated prompts taken from QCD only
    if( isFake   && isGJet ) continue; //fakes from QCD only


    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );

    float hOverE = myTree.gamma_hOverE[0];
    float sietaieta = myTree.gamma_sigmaIetaIeta[0];
    if( fabs( gamma.Eta() )<1.479 ) {
      if( hOverE > 0.058 ) continue;
      if( sietaieta > 0.01 ) continue;
    } else {  
      if( hOverE > 0.020 ) continue;
      if( sietaieta > 0.03 ) continue;
    }


    //// relative iso:
    //float iso = myTree.gamma_chHadIso[0]/myTree.gamma_pt[0];
    //if( iso>isoCut ) continue;
    //if( iso>0.1 ) continue; // preselection anyways in there

    // absolute iso:
    float iso = myTree.gamma_chHadIso[0];
    if( iso>20. ) continue; // preselection anyways in there



    int closestJet = -1;
    float deltaRmin = 0.4;
    for( unsigned i=0; i<myTree.njet; ++i ) {
      if( fabs(myTree.jet_eta[i])>2.5 ) continue;
      if( myTree.jet_pt[i]<40. ) continue;
      TLorentzVector thisjet;
      thisjet.SetPtEtaPhiM( myTree.jet_pt[i], myTree.jet_eta[i], myTree.jet_phi[i], myTree.jet_mass[i] );
      float thisDeltaR = gamma.DeltaR(thisjet);
      if( thisDeltaR<deltaRmin ) {
        deltaRmin = thisDeltaR;
        closestJet = i;
      }
    }
    float found_pt = 0.;
    int jet_counter = 0;
    for( unsigned i=0; i<myTree.njet; ++i ) {
      if( i==closestJet ) continue;
      if( fabs(myTree.jet_eta[i])>2.5 ) continue;
      if( myTree.jet_pt[i]<40. ) continue;
      jet_counter++;
      if( jet_counter==2 ) {
        found_pt = myTree.jet_pt[i];
        break;
      }
    }

    if( found_pt<100. ) continue;


    Double_t weight = myTree.evt_scale1fb*lumi; 

    if( isPrompt ) {

      MT2EstimateZinvGamma* thisPrompt = prompt->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
      if( thisPrompt==0 ) continue;

      thisPrompt->yield->Fill(myTree.gamma_mt2, weight );
      thisPrompt->fillIso( iso, weight, myTree.gamma_mt2 );

      if( iso<isoCut ) {

        MT2EstimateZinvGamma* thisPrompt_pass = prompt_pass->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
        if( thisPrompt_pass==0 ) continue;

        thisPrompt_pass->yield->Fill(myTree.gamma_mt2, weight );
        thisPrompt_pass->fillIso( iso, weight, myTree.gamma_mt2 );

      }

    } else if( isNIP ) { 
      
      MT2EstimateZinvGamma* thisnip = nip->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
      if( thisnip==0 ) continue;

      thisnip->yield->Fill(myTree.gamma_mt2, weight );
      thisnip->fillIso( iso, weight, myTree.gamma_mt2 );

      if( iso<isoCut ) {

        MT2EstimateZinvGamma* thisnip_pass = nip_pass->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
        if( thisnip_pass==0 ) continue;

        thisnip_pass->yield->Fill(myTree.gamma_mt2, weight );
        thisnip_pass->fillIso( iso, weight, myTree.gamma_mt2 );

      }

    } else if( isFake ) {

      MT2EstimateZinvGamma* thisFake = fake->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
      if( thisFake==0 ) continue;

      thisFake->yield->Fill(myTree.gamma_mt2, weight );
      thisFake->fillIso( iso, weight, myTree.gamma_mt2 );

      if( iso<isoCut ) {

        MT2EstimateZinvGamma* thisFake_pass = fake_pass->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
        if( thisFake_pass==0 ) continue;

        thisFake_pass->yield->Fill(myTree.gamma_mt2, weight );
        thisFake_pass->fillIso( iso, weight, myTree.gamma_mt2 );

      }

    } // is prompt/nip/fake

    
  } // for entries


  prompt->finalize();
  prompt_pass->finalize();
  nip->finalize();
  nip_pass->finalize();
  fake->finalize();
  fake_pass->finalize();
  

  delete tree;


  file->Close();
  delete file;
  

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

