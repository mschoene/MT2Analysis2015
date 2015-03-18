#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"
#include "TRandom3.h"


float lumi = 5.; //fb-1



void computeYield( const MT2Sample& sample, const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* prompt, MT2Analysis<MT2EstimateTree>* nip, MT2Analysis<MT2EstimateTree>* fake );




int main( int argc, char* argv[] ) {


  std::string samplesFileName = "PHYS14_v2_Zinv";
  //std::string samplesFileName = "CSA14_Zinv";
  if( argc>1 ) {
    std::string samplesFileName_tmp(argv[1]); 
    samplesFileName = samplesFileName_tmp;
  }

  std::string samplesFile = "../samples/samples_" + samplesFileName + ".dat";
  
  std::vector<MT2Sample> samples = MT2Sample::loadSamples(samplesFile, 100, 299); // GJet + QCD
  if( samples.size()==0 ) {
    std::cout << "There must be an error: didn't find any samples in " << samplesFile << "!" << std::endl;
    exit(1209);
  }

  


  //std::string regionsSet = "13TeV_CSA14";
  //std::string regionsSet = "13TeV_onlyHT";
  std::string regionsSet = "13TeV_inclusive";
  //std::string regionsSet = "13TeV_ZinvGammaPurity";

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string outputdir = "GenIsoCheck_" + samplesFileName + "_" + regionsSet;
  system(Form("mkdir -p %s", outputdir.c_str()));

  

  MT2Analysis<MT2EstimateTree>* prompt = new MT2Analysis<MT2EstimateTree>( "prompt", regionsSet );
  MT2Analysis<MT2EstimateTree>* nip = new MT2Analysis<MT2EstimateTree>( "nip", regionsSet );
  MT2Analysis<MT2EstimateTree>* fake = new MT2Analysis<MT2EstimateTree>( "fake", regionsSet );

  MT2EstimateTree::addVar( prompt, "genIso" );
  MT2EstimateTree::addVar( prompt, "iso" );
  MT2EstimateTree::addVar( prompt, "isoRC" );
  MT2EstimateTree::addVar( prompt, "isoCN" );
  MT2EstimateTree::addVar( prompt, "isoCP" );
  MT2EstimateTree::addVar( prompt, "isoCPN" );
  MT2EstimateTree::addVar( prompt, "sietaieta" );
  MT2EstimateTree::addVar( prompt, "mcMatchId" );
  MT2EstimateTree::addVar( prompt, "ptGamma" );
  MT2EstimateTree::addVar( prompt, "etaGamma" );

  MT2EstimateTree::addVar( nip, "genIso" );
  MT2EstimateTree::addVar( nip, "iso" );
  MT2EstimateTree::addVar( nip, "isoRC" );
  MT2EstimateTree::addVar( nip, "isoCN" );
  MT2EstimateTree::addVar( nip, "isoCP" );
  MT2EstimateTree::addVar( nip, "isoCPN" );
  MT2EstimateTree::addVar( nip, "sietaieta" );
  MT2EstimateTree::addVar( nip, "mcMatchId" );
  MT2EstimateTree::addVar( nip, "ptGamma" );
  MT2EstimateTree::addVar( nip, "etaGamma" );

  MT2EstimateTree::addVar( fake, "genIso" );
  MT2EstimateTree::addVar( fake, "iso" );
  MT2EstimateTree::addVar( fake, "isoRC" );
  MT2EstimateTree::addVar( fake, "isoCN" );
  MT2EstimateTree::addVar( fake, "isoCP" );
  MT2EstimateTree::addVar( fake, "isoCPN" );
  MT2EstimateTree::addVar( fake, "sietaieta" );
  MT2EstimateTree::addVar( fake, "mcMatchId" );
  MT2EstimateTree::addVar( fake, "ptGamma" );
  MT2EstimateTree::addVar( fake, "etaGamma" );


  for( unsigned i=0; i<samples.size(); ++i ) {
    computeYield( samples[i], regionsSet, prompt, nip, fake );
  }


  prompt->writeToFile( outputdir + "/genIso.root" );
  nip->addToFile( outputdir + "/genIso.root" );
  fake->addToFile( outputdir + "/genIso.root" );


  return 0;

}








void computeYield( const MT2Sample& sample, const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* prompt, MT2Analysis<MT2EstimateTree>* nip, MT2Analysis<MT2EstimateTree>* fake ) {

  
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

    if( myTree.mt2>200. ) continue; // orthogonal to signal regions

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
    if( myTree.gamma_pt[0]<160. ) continue;


    //// aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa REMOVE SOON
    if( myTree.evt_scale1fb > 1. ) continue;

    //if( myTree.gamma_idCutBased[0]==0 ) continue;


    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );

    float hOverE = myTree.gamma_hOverE[0];
    float sietaieta = myTree.gamma_sigmaIetaIeta[0];
    if( fabs( gamma.Eta() )<1.479 ) {
      if( hOverE > 0.058 ) continue;
      //if( sietaieta > 0.014 ) continue;
      //if( sietaieta > 0.0099 ) continue;
    } else {  
      if( hOverE > 0.020 ) continue;
      if( sietaieta > 0.035 ) continue;
      //if( sietaieta > 0.0268 ) continue;
    }

    //if( myTree.gamma_chHadIso[0]+myTree.gamma_neuHadIso[0] > 10. ) continue;
    float iso = myTree.gamma_chHadIso[0]/myTree.gamma_pt[0];

    float effAreaCH = 0.; // not used now
    float effAreaNH = 0.;
    float effAreaPH = 0.;
    if( fabs(gamma.Eta())< 1. ) {
      effAreaCH = 0.0130;
      effAreaNH = 0.0056;
      effAreaPH = 0.0896;
    } else if( fabs(gamma.Eta())< 1.479 ) {
      effAreaCH = 0.0096;
      effAreaNH = 0.0107;
      effAreaPH = 0.0762;
    } else if( fabs(gamma.Eta())< 2.0   ) {
      effAreaCH = 0.0107;
      effAreaNH = 0.0019;
      effAreaPH = 0.0383;
    } else if( fabs(gamma.Eta())< 2.2   ) {
      effAreaCH = 0.0077;
      effAreaNH = 0.0011;
      effAreaPH = 0.0534;
    } else if( fabs(gamma.Eta())< 2.3   ) {
      effAreaCH = 0.0088;
      effAreaNH = 0.0077;
      effAreaPH = 0.0846;
    } else if( fabs(gamma.Eta())< 2.4   ) {
      effAreaCH = 0.0065;
      effAreaNH = 0.0178;
      effAreaPH = 0.1032;
    } else                                {
      effAreaCH = 0.0030;
      effAreaNH = 0.1675;
      effAreaPH = 0.1598;
    }

    float isoNH = myTree.gamma_neuHadIso[0] - effAreaNH*myTree.rho;
    float isoPH = myTree.gamma_phIso[0]     - effAreaPH*myTree.rho;
    if( isoNH<0. ) isoNH = 0.;
    if( isoPH<0. ) isoPH = 0.;
    float isoCN = (myTree.gamma_chHadIso[0]+isoNH)/myTree.gamma_pt[0];
    float isoCP = (myTree.gamma_chHadIso[0]+isoPH)/myTree.gamma_pt[0];
    float isoCPN = (myTree.gamma_chHadIso[0]+isoNH+isoPH)/myTree.gamma_pt[0];
    //if( iso>0.1 ) continue; // preselection anyways in there

    float isoRC = myTree.gamma_chHadIsoRC[0];

    int mcMatchId = myTree.gamma_mcMatchId[0];
    bool isMatched = (mcMatchId==22 || mcMatchId==7);
    bool isGenIso = (myTree.gamma_genIso[0]<5.);

    bool isPrompt = ( isMatched &&  isGenIso);
    bool isNIP    = ( isMatched && !isGenIso);
    bool isFake   = (!isMatched);

    if( isPrompt && isQCD  ) continue; //isolated prompts taken from GJet only
    if( isNIP    && isGJet ) continue; //non-isolated prompts taken from QCD only
    if( isFake   && isGJet ) continue; //fakes from QCD only


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

      MT2EstimateTree* thisPrompt = prompt->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
      if( thisPrompt==0 ) continue;

      thisPrompt->yield->Fill(myTree.gamma_mt2, weight );

      thisPrompt->assignTree(myTree, lumi*myTree.evt_scale1fb);
      thisPrompt->assignVars( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt, myTree.gamma_mt2 );
      thisPrompt->assignVar( "genIso", myTree.gamma_genIso[0] );
      thisPrompt->assignVar( "sietaieta", sietaieta );
      thisPrompt->assignVar( "iso", iso );
      thisPrompt->assignVar( "isoRC", isoRC );
      thisPrompt->assignVar( "isoCN", isoCN );
      thisPrompt->assignVar( "isoCP", isoCP );
      thisPrompt->assignVar( "isoCPN", isoCPN );
      thisPrompt->assignVar( "mcMatchId", mcMatchId );
      thisPrompt->assignVar( "ptGamma", gamma.Pt() );
      thisPrompt->assignVar( "etaGamma", gamma.Eta() );
      thisPrompt->tree->Fill();

    } else if( isNIP ) {

      MT2EstimateTree* thisNIP = nip->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
      if( thisNIP==0 ) continue;

      thisNIP->yield->Fill(myTree.gamma_mt2, weight );

      thisNIP->assignTree(myTree, lumi*myTree.evt_scale1fb);
      thisNIP->assignVars( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt, myTree.gamma_mt2 );
      thisNIP->assignVar( "genIso", myTree.gamma_genIso[0] );
      thisNIP->assignVar( "sietaieta", sietaieta );
      thisNIP->assignVar( "iso", iso );
      thisNIP->assignVar( "isoRC", isoRC );
      thisNIP->assignVar( "isoCN", isoCN );
      thisNIP->assignVar( "isoCP", isoCP );
      thisNIP->assignVar( "isoCPN", isoCPN );
      thisNIP->assignVar( "mcMatchId", mcMatchId );
      thisNIP->assignVar( "ptGamma", gamma.Pt() );
      thisNIP->assignVar( "etaGamma", gamma.Eta() );
      thisNIP->tree->Fill();

    } else if( isFake ) {

      MT2EstimateTree* thisFake = fake->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
      if( thisFake==0 ) continue;

      thisFake->yield->Fill(myTree.gamma_mt2, weight );

      thisFake->assignTree(myTree, lumi*myTree.evt_scale1fb);
      thisFake->assignVars( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt, myTree.gamma_mt2 );
      thisFake->assignVar( "genIso", myTree.gamma_genIso[0] );
      thisFake->assignVar( "sietaieta", sietaieta );
      thisFake->assignVar( "iso", iso );
      thisFake->assignVar( "isoRC", isoRC );
      thisFake->assignVar( "isoCN", isoCN );
      thisFake->assignVar( "isoCP", isoCP );
      thisFake->assignVar( "isoCPN", isoCPN );
      thisFake->assignVar( "mcMatchId", mcMatchId );
      thisFake->assignVar( "ptGamma", gamma.Pt() );
      thisFake->assignVar( "etaGamma", gamma.Eta() );
      thisFake->tree->Fill();

    }


  } // for entries


  prompt->finalize();
  nip->finalize();
  fake->finalize();
  

  delete tree;


  file->Close();
  delete file;
  
}



