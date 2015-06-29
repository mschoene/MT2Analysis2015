#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"





float lumi = 20.; //fb-1





void computeYield( const MT2Sample& sample, const std::string& regionsSet, MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* fake, const std::string& useMC );
void setPoissonError( MT2Analysis<MT2EstimateZinvGamma>* data );
//void randomizePoisson( MT2Analysis<MT2EstimateZinvGamma>* data );




int main( int argc, char* argv[] ) {



  std::string useMC = "dataRC";

  if( argc>1 ) {

    useMC = std::string(argv[1]); 

    if( useMC=="dataFR" ) useMC="DataFR"; // data Fake Removal
    if( useMC=="dataRC" ) useMC="DataRC"; // data Random Cone
    if( useMC=="data"   ) {
      std::cout << std::endl;
      std::cout << "-> Asking for 'data': will use data Random Cone. (default)" << std::endl;
      std::cout << std::endl;
      useMC="DataRC"; // (default for data)
    }

    if( useMC!="data" && useMC!="DataFR" && useMC!="MC" && useMC!="DataRC" ) {
      std::cout << "ERROR! Second argument may only be 'MC' or 'dataFR' or 'dataRC'" << std::endl;
      exit(1111);
    }

  }


  //std::string regionsSet = "13TeV_onlyHT";
  std::string regionsSet = "13TeV_inclusive";
  //std::string regionsSet = "13TeV_inclusive";
  //std::string regionsSet = "13TeV_ZinvGammaPurity";
  if( argc>2 ) {
    std::string regionsSet_tmp(argv[2]); 
    regionsSet = regionsSet_tmp;
  }



  std::string samplesFileName = "PHYS14_v5_skimprune";
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


  if( useMC=="DataFR" || useMC=="DataRC" ) {
    setPoissonError( templatesFake );
    setPoissonError( templatesPrompt );
    if( useMC=="DataFR" ) templatesPrompt->setName("templatesPromptRaw");
  }



  std::string templateFileName = "gammaTemplates" + useMC;
  templateFileName = templateFileName + "_" + samplesFileName + "_" + regionsSet + ".root";


  templatesFake->writeToFile(templateFileName);
  templatesPrompt->addToFile(templateFileName, true);


  return 0;

}









void computeYield( const MT2Sample& sample, const std::string& regionsSet, MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* fake, const std::string& useMC ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;


  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  int nentries = tree->GetEntries();




  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( !myTree.passSelection("gamma") ) continue;
    if( !myTree.passGammaAdditionalSelection(sample.id) ) continue;


    // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH REMOVE THIS SOOOOON
    if( myTree.evt_scale1fb>1. ) continue;



    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );


    float ht        = myTree.gamma_ht;
    float met       = myTree.gamma_met_pt;
    float mt2       = myTree.gamma_mt2;
    float minMTBmet = myTree.gamma_minMTBMet;
    int njets       = myTree.gamma_nJet40;
    int nbjets      = myTree.gamma_nBJet20;    



    float iso = myTree.gamma_chHadIso[0];

    float sietaieta = myTree.gamma_sigmaIetaIeta[0];
    float id_thresh;
    bool sietaietaOK = false;
    if( fabs( gamma.Eta() )<1.479 ) {
      id_thresh = 0.0106;
      if( sietaieta>id_thresh && sietaieta<0.011 ) continue; // no man's land
      if( sietaieta>0.015 ) continue; // end of sidebands
    } else {  
      id_thresh = 0.0266;
      if( sietaieta>id_thresh && sietaieta<0.030 ) continue; // no man's land
      if( sietaieta>0.035 ) continue; // end of sidebands
    }
    sietaietaOK = (sietaieta < id_thresh);

    bool isWorkingPrompt = false;

    if( useMC=="MC" ) {

      if( !sietaietaOK ) continue;
      isWorkingPrompt = myTree.gamma_mcMatchId[0]==22; // prompt = matched

    } else if( useMC=="DataFR" ) { 

      isWorkingPrompt = sietaietaOK;

      if( isWorkingPrompt ) {
        // for prompts use only low-sensitivity regions:
        if( mt2>300. ) continue;
        if( ht>1000. ) continue;
        if( njets>6 ) continue;
        if( nbjets>0 ) continue;
      }


    } else if( useMC=="DataRC" ) { 

      isWorkingPrompt = sietaietaOK;

      if( isWorkingPrompt ) iso = myTree.gamma_chHadIsoRC[0]; // random cone


    } else {

      // this shouldnt be possible
      std::cout << "-> NOPE. Don't know anything about useMC='" << useMC << "'." << std::endl;
      exit(191);

    }


    //// preselection
    //if( iso > 20. ) continue;
    ////if( iso > 10. ) continue;

    Double_t weight = myTree.evt_scale1fb*lumi; 


    if( isWorkingPrompt ) {

      MT2EstimateZinvGamma* thisPrompt = prompt->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      if( thisPrompt==0 ) continue;

      thisPrompt->yield->Fill(mt2, weight );
      thisPrompt->sietaieta->Fill(myTree.gamma_sigmaIetaIeta[0], weight );
      thisPrompt->fillIso( iso, weight, mt2 );

    } else {

      MT2EstimateZinvGamma* thisFake = fake->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      if( thisFake==0 ) continue;

      thisFake->yield->Fill(mt2, weight );
      thisFake->sietaieta->Fill(myTree.gamma_sigmaIetaIeta[0], weight );
      thisFake->fillIso( iso, weight, mt2 );

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
    
    for( int ibin=1; ibin<h1_iso->GetXaxis()->GetNbins()+1; ++ibin ) {
    
      h1_iso->SetBinError(ibin, sqrt(h1_iso->GetBinContent(ibin)));
      
    }  // for bins

  }// for MT2 regions


}
