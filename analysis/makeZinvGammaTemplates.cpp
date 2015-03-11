#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"





float lumi = 5.; //fb-1





void computeYield( const MT2Sample& sample, const std::string& regionsSet, MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* fake, bool useMC );
void setPoissonError( MT2Analysis<MT2EstimateZinvGamma>* data );
//void randomizePoisson( MT2Analysis<MT2EstimateZinvGamma>* data );




int main( int argc, char* argv[] ) {


  if( argc==1 ) {
    std::cout << "-> You need to pass me the regions set name. Here are some suggestions: " << std::endl;
    std::cout << "  13TeV_CSA14" << std::endl;
    std::cout << "  13TeV_onlyHT" << std::endl;
    std::cout << "  13TeV_onlyJets" << std::endl;
    std::cout << "  13TeV_inclusive" << std::endl;
    exit(101);
  }

  //std::string regionsSet = "13TeV_onlyHT";
  std::string regionsSet = "13TeV_CSA14";
  //std::string regionsSet = "13TeV_inclusive";
  //std::string regionsSet = "13TeV_ZinvGammaPurity";
  if( argc>1 ) {
    std::string regionsSet_tmp(argv[1]); 
    regionsSet = regionsSet_tmp;
  }


  bool useMC = true;
  if( argc>2 ) {
    std::string useMC_str(argv[2]); 
    if( useMC_str!="data" && useMC_str!="MC" ) {
      std::cout << "ERROR! Second argument may only be 'MC' or 'data'" << std::endl;
      exit(1111);
    }
    if( useMC_str=="data" ) useMC = false;
  }



  std::string samplesFileName = "PHYS14_v2_Zinv";
  std::string samplesFile = "../samples/samples_" + samplesFileName + ".dat";
  
  std::vector<MT2Sample> samples = MT2Sample::loadSamples(samplesFile, 100, 299); // GJet and QCD
  if( samples.size()==0 ) {
    std::cout << "There must be an error: didn't find any good files in " << samplesFile << "!" << std::endl;
    exit(1209);
  }



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string outputdir = "ZinvGammaPurity_" + samplesFileName + "_" + regionsSet;
  system(Form("mkdir -p %s", outputdir.c_str()));


  
  MT2Analysis<MT2EstimateZinvGamma>* templatesPrompt = new MT2Analysis<MT2EstimateZinvGamma>( "templatesPrompt", regionsSet );
  MT2Analysis<MT2EstimateZinvGamma>* templatesFake   = new MT2Analysis<MT2EstimateZinvGamma>( "templatesFake", regionsSet );

  for( unsigned i=0; i<samples.size(); ++i ) {
    computeYield( samples[i], regionsSet, templatesPrompt, templatesFake, useMC );
  }



  std::string templateFileName = "gammaTemplates";
  if( useMC ) templateFileName += "MC";
  else        templateFileName += "Data";
  templateFileName = templateFileName + "_" + samplesFileName + "_" + regionsSet + ".root";



  if( !useMC ) {
    setPoissonError( templatesFake );
    setPoissonError( templatesPrompt );
  }
  templatesFake->writeToFile(templateFileName);

  if( !useMC ) 
    templatesPrompt->setName( "templatesPromptRaw" ); // still contaminated from fakes with good sietaieta
  templatesPrompt->addToFile(templateFileName, true);


  return 0;

}









void computeYield( const MT2Sample& sample, const std::string& regionsSet, MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* fake, bool useMC ) {


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

    //if( myTree.gamma_ht>1000. && sample.id==204 ) continue; // remove high-weight spikes (remove GJet_400to600 leaking into HT>1000)

    if( myTree.met_pt > 100.) continue;
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


    // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH REMOVE THIS SOOOOON
    if( myTree.evt_scale1fb>1. ) continue;
    

    int mcMatchId = myTree.gamma_mcMatchId[0];
    bool isMatched = (mcMatchId==22 || mcMatchId==7);
    bool isGenIso = (myTree.gamma_genIso[0]<5.);

    if( isMatched  &&  isGenIso && isQCD  ) continue; //isolated prompts taken from GJet only
    if( isMatched  && !isGenIso && isGJet ) continue; //non-isolated prompts taken from QCD only
    if( !isMatched &&              isGJet ) continue; //fakes from QCD only


    float iso = myTree.gamma_chHadIso[0];
    if( iso > 30. ) continue;
    //float iso = myTree.gamma_chHadIso[0]/myTree.gamma_pt[0];
    //if( iso > 0.1 ) continue;


    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );


    float hOverE = myTree.gamma_hOverE[0];
    float sietaieta = myTree.gamma_sigmaIetaIeta[0];
    bool sietaietaOK = false;
    if( fabs( gamma.Eta() )<1.479 ) {
      if( hOverE > 0.058 ) continue;
      if( sietaieta>0.010 && sietaieta<0.011 ) continue; // no man's land
      if( sietaieta>0.015 ) continue; // end of sidebands
      sietaietaOK = (sietaieta < 0.01);
    } else {  
      if( hOverE > 0.020 ) continue;
      if( sietaieta>0.035 ) continue; // end of sidebands
      sietaietaOK = (sietaieta < 0.03);
    }

    bool isPrompt = false;

    if( useMC ) {

      if( !sietaietaOK ) continue;
      isPrompt = isMatched;

    } else { // is data

      isPrompt = sietaietaOK;

      if( isPrompt ) {
        // for prompts use only low-sensitivity regions:
        if( myTree.gamma_mt2>300. ) continue;
        if( myTree.gamma_ht>1000. ) continue;
        if( myTree.gamma_nBJet40>0 ) continue;
      }

    }



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
      thisPrompt->sietaieta->Fill(myTree.gamma_sigmaIetaIeta[0], weight );
      thisPrompt->fillIso( iso, weight, myTree.gamma_mt2 );

    } else {

      MT2EstimateZinvGamma* thisFake = fake->get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
      if( thisFake==0 ) continue;

      thisFake->yield->Fill(myTree.gamma_mt2, weight );
      thisFake->sietaieta->Fill(myTree.gamma_sigmaIetaIeta[0], weight );
      thisFake->fillIso( iso, weight, myTree.gamma_mt2 );

    }

    
  } // for entries


  //prompt->finalize();
  //fake->finalize();
  

  delete tree;


  file->Close();
  delete file;
  

}



//void randomizePoisson( MT2Analysis<MT2EstimateZinvGamma>* data ) {
//
//  TRandom3 rand(13);
//
//
//  std::set<MT2Region> MT2Regions = data->getRegions();
//
//  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
//
//    TH1D* h1_iso = data->get(*iMT2)->iso;
//    
//    for( unsigned ibin=1; ibin<h1_iso->GetXaxis()->GetNbins()+1; ++ibin ) {
//    
//      int poisson_data = rand.Poisson(h1_iso->GetBinContent(ibin));
//      //h1_iso->SetBinContent(ibin, poisson_data);
//      h1_iso->SetBinError(ibin, sqrt(poisson_data));
//      
//    }  // for bins
//
//  }// for MT2 regions
//
//
//}


void setPoissonError( MT2Analysis<MT2EstimateZinvGamma>* data ) {


  std::set<MT2Region> MT2Regions = data->getRegions();

  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    TH1D* h1_iso = data->get(*iMT2)->iso;
    
    for( unsigned ibin=1; ibin<h1_iso->GetXaxis()->GetNbins()+1; ++ibin ) {
    
      h1_iso->SetBinError(ibin, sqrt(h1_iso->GetBinContent(ibin)));
      
    }  // for bins

  }// for MT2 regions


}
