#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
//#include <cmath>


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
#include "interface/MT2DrawTools.h"
#include "interface/MT2Config.h"

#include "TRandom3.h"


#define mt2_cxx
#include "interface/mt2.h"





void randomizePoisson( MT2Analysis<MT2EstimateTree>* data );
template <class T>
MT2Analysis<T>* computeYield( const MT2Sample& sample, const MT2Config& cfg );
MT2Analysis<MT2EstimateTree>* mergeYields( std::vector< MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max=-1, const std::string& legendName="" );
int matchPartonToJet( int index, MT2Tree* myTree );





int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|             Running regionEventYields              |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./regionEventYields [configFileName] [data/MC]" << std::endl;
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
    else if( dataMC=="MC" ) onlyMC = true;
    else {
      std::cout << "-> You passed a second argument that isn't 'data' nor 'MC', so I don't know what to do about it." << std::endl;
    }
  }

  std::string outputdir = cfg.getEventYieldDir();
  system(Form("mkdir -p %s", outputdir.c_str()));


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::vector<MT2Analysis<MT2EstimateTree>* > yields;
  //MT2Analysis<MT2EstimateTree>* dataYield;  

  if( cfg.useMC() && !onlyData ) { // use MC BG estimates

    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 100, 999); // not interested in signal here (see later)
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(120);
    }


    
    std::vector< MT2Analysis<MT2EstimateTree>* > EventYield;
    for( unsigned i=0; i<fSamples.size(); ++i ) {
      int this_id = fSamples[i].id;
      if( this_id>=200 && this_id<300 ) continue; // skip GJets
      if( this_id>=700 && this_id<800 ) continue; // skip DY
      //if( !(this_id>=600 && this_id<700) ) continue; // only Zinv
      EventYield.push_back( computeYield<MT2EstimateTree>( fSamples[i], cfg ));
    }


    std::cout << "-> Done looping on samples. Start merging." << std::endl;

    std::cout << "     merging Top..." << std::endl;
    MT2Analysis<MT2EstimateTree>* EventYield_top   = mergeYields( EventYield, cfg.regionsSet(), "Top", 300, 499 ); // ttbar, single top, ttW, ttZ...
    std::cout << "     merging QCD..." << std::endl;
    MT2Analysis<MT2EstimateTree>* EventYield_qcd   = mergeYields( EventYield, cfg.regionsSet(), "QCD", 100, 199 );
    std::cout << "     merging WJets..." << std::endl;
    MT2Analysis<MT2EstimateTree>* EventYield_wjets = mergeYields( EventYield, cfg.regionsSet(), "WJets", 500, 599, "W+jets" );
    std::cout << "     merging ZJets..." << std::endl;
    MT2Analysis<MT2EstimateTree>* EventYield_zjets = mergeYields( EventYield, cfg.regionsSet(), "ZJets", 600, 699, "Z+jets" );
    //MT2Analysis<MT2EstimateTree>* EventYield_other = mergeYields( EventYield, cfg.regionsSet(), "Diboson", 700, 899, "Other" );
    std::cout << "-> Done merging." << std::endl;

    yields.push_back( EventYield_qcd );
    yields.push_back( EventYield_wjets );
    yields.push_back( EventYield_zjets );
    yields.push_back( EventYield_top );
    //yields.push_back( EventYield_other );

    if( cfg.dummyAnalysis() ) {
      MT2Analysis<MT2EstimateTree>* dataYield   = mergeYields( EventYield, cfg.regionsSet(), "data", 100, 699 );
      yields.push_back( dataYield );
    } //else {
      //dataYield = new MT2Analysis<MT2EstimateTree>( "data", cfg.regionsSet(), 1 );
    //}

  } // if MC samples



  // load signal samples, if any
  std::vector< MT2Analysis< MT2EstimateTree>* > signals;
  if( cfg.mcSamples()!="" && cfg.additionalStuff()!="noSignals" && !onlyData ) {

    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading signal samples from file: " << samplesFileName << std::endl;

    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 1000); // only signal (id>=1000)


    if( fSamples.size()==0 ) {

      std::cout << "No signal samples found, skipping." << std::endl;

    } else {
    
      for( unsigned i=0; i<fSamples.size(); ++i ) 
        signals.push_back( computeYield<MT2EstimateTree>( fSamples[i], cfg ) );
    
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
        signals.push_back( computeYield<MT2EstimateTree>( fSamples[i], cfg ) );

    } // if samples != 0
    
  } // if sig samples
  

  if( !(cfg.dummyAnalysis()) && cfg.dataSamples()!="" && !onlyMC ) {

    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "JetHTMHT"); //, 1, 99 );
    if( samples_data.size()==0 ) {
      std::cout << "There must be an error: samples_data is empty!" << std::endl;
      exit(1209);
    }

    std::vector< MT2Analysis<MT2EstimateTree>* > EventYield_data;
    for( unsigned i=0; i < samples_data.size(); ++i )
      EventYield_data.push_back( computeYield<MT2EstimateTree>( samples_data[i], cfg ) );

    //dataYield   = mergeYields( EventYield_data, cfg.regionsSet(), "data"); //, 1, 99 );
    MT2Analysis<MT2EstimateTree>* dataYield   = EventYield_data[0];
    dataYield->setName("data");

    yields.push_back( dataYield );

  }


  if( yields.size()==0 ) {
    std::cout << "-> Didn't end up with a single yield... something's wrong." << std::endl;
    exit(87);
  }


  // save MT2Analyses:
  yields[0]->writeToFile(outputdir + "/analyses.root");
  for( unsigned i=1; i<yields.size(); ++i )
    yields[i]->writeToFile(outputdir + "/analyses.root");
  for( unsigned i=0; i<signals.size(); ++i )
    signals[i]->writeToFile(outputdir + "/analyses.root");

  cfg.saveAs(outputdir + "/config.txt");

  return 0;

}




template <class T>
MT2Analysis<T>* computeYield( const MT2Sample& sample, const MT2Config& cfg ) {


  std::string regionsSet = cfg.regionsSet();

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
  

  MT2Tree myTree;
  if( cfg.additionalStuff()=="qgVars" ) {
     myTree.loadGenStuff = true;
  } else {
    myTree.loadGenStuff = false;
  }
  myTree.Init(tree);



  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  MT2Analysis<T>* analysis = new MT2Analysis<T>( sample.sname, regionsSet, sample.id );
  
 
  if( cfg.additionalStuff()=="qgVars" ) {
    T::addVar( analysis, "partId0" );
    T::addVar( analysis, "partId1" );
    T::addVar( analysis, "partId2" );
    T::addVar( analysis, "partId3" );
    T::addVar( analysis, "qgl0" );
    T::addVar( analysis, "qgl1" );
    T::addVar( analysis, "qgl2" );
    T::addVar( analysis, "qgl3" );
    T::addVar( analysis, "qglProd" );
    T::addVar( analysis, "qglAve" );
  }
  



  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( regionsSet!="13TeV_noCut" )
      if( !myTree.passSelection(cfg.additionalStuff()) ) continue;

    
    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20;    
    float mt2  = (njets>1) ? myTree.mt2 : myTree.jet_pt[0];
    //float mt2  = myTree.mt2;
    
    float GenSusyMScan1 = myTree.GenSusyMScan1;
    float GenSusyMScan2 = myTree.GenSusyMScan2;
    
    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi()*myTree.puWeight;
    //weight *= myTree.weight_lepsf;

    if( myTree.isData ) {
      if( !(  (myTree.HLT_PFHT800 && ht>=900.) || (myTree.HLT_PFHT350_PFMET100 && ht<900.)  ) ) continue;
    }

   
    T* thisEstimate = analysis->get( ht, njets, nbjets, met, minMTBmet, mt2 );
    if( thisEstimate==0 ) continue;


//    //////QCD
//    if( ht > 575. && ht < 1000. && mt2 < 300. ) continue;
//    else if( ht > 1000 && ht < 1500. && mt2 < 300. ) continue;
//    else if( ht > 1500. && mt2 < 400. ) continue;
//    //////

    if( cfg.additionalStuff()=="qgVars" ) {

      // initialize
      thisEstimate->assignVar( "qgl0", -1. );
      thisEstimate->assignVar( "qgl1", -1. );
      thisEstimate->assignVar( "qgl2", -1. );
      thisEstimate->assignVar( "qgl3", -1. );
      thisEstimate->assignVar( "partId0", 0 );
      thisEstimate->assignVar( "partId1", 0 );
      thisEstimate->assignVar( "partId2", 0 );
      thisEstimate->assignVar( "partId3", 0 );

      float qglProd = 1.;
      float qglAve = 0.;
      int denom = 0;


      if( njets>0 && fabs(myTree.jet_eta[0])<2.5 ) {

        float qgl0 = myTree.jet_qgl[0];
        thisEstimate->assignVar( "qgl0", qgl0 );
        qglProd *= qgl0;
        qglAve += qgl0;
        denom++;
        thisEstimate->assignVar( "partId0", matchPartonToJet( 0, &myTree ) );
        //thisEstimate->assignVar( "partId0", myTree.jet_mcFlavour[0] );

      }


      if( njets>1 && fabs(myTree.jet_eta[1])<2.5 ) {

        float qgl1 = myTree.jet_qgl[1];
        thisEstimate->assignVar( "qgl1", qgl1 );
        qglProd *= qgl1;
        qglAve += qgl1;
        denom++;

        thisEstimate->assignVar( "partId1", matchPartonToJet( 1, &myTree ) );
        //thisEstimate->assignVar( "partId1", myTree.jet_mcFlavour[1] );

      }
        
      if( njets>2 && fabs(myTree.jet_eta[2])<2.5 ) {

        float qgl2 = myTree.jet_qgl[2];
        thisEstimate->assignVar( "qgl2", qgl2 );
        qglProd *= qgl2;
        qglAve += qgl2;
        denom++;

        thisEstimate->assignVar( "partId2", matchPartonToJet( 2, &myTree ) );
        //thisEstimate->assignVar( "partId2", myTree.jet_mcFlavour[2] );

      }
        

      if( njets>3 && fabs(myTree.jet_eta[3])<2.5 ) {

        float qgl3 = myTree.jet_qgl[3];
        thisEstimate->assignVar( "qgl3", qgl3 );
        qglProd *= qgl3;
        qglAve += qgl3;
        denom++;

        thisEstimate->assignVar( "partId3", matchPartonToJet( 3, &myTree ) );
        //thisEstimate->assignVar( "partId3", myTree.jet_mcFlavour[3] );

      }

      qglAve /= (float)denom;
        
      thisEstimate->assignVar( "qglProd", qglProd );
      thisEstimate->assignVar( "qglAve", qglAve );

      thisEstimate->assignTree( myTree, weight );
      thisEstimate->tree->Fill();

    } else {

      thisEstimate->fillTree( myTree, weight );

    }

    thisEstimate->yield->Fill( mt2, weight );
    thisEstimate->yield3d->Fill( mt2, GenSusyMScan1, GenSusyMScan2, weight );

  } // for entries
    
  //ofs.close();

  analysis->finalize();
  
  delete tree;

  file->Close();
  delete file;
  
  return analysis;

}



MT2Analysis<MT2EstimateTree>* mergeYields( std::vector<MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max, const std::string& legendName ) {

  if( id_max<0 ) id_max=id_min;

  MT2Analysis<MT2EstimateTree>* return_EventYield = new MT2Analysis<MT2EstimateTree>(name, regionsSet, id_min, legendName);

  for( unsigned i=0; i<EventYield.size(); ++i ) {

    if( EventYield[i]->getId() >= id_min && EventYield[i]->getId() <= id_max ) {

       *(return_EventYield) += *(EventYield[i]);

    }

  } // for EventYield


  return return_EventYield;

}





void randomizePoisson( MT2Analysis<MT2EstimateTree>* data ) {

  TRandom3 rand(13);


//  std::set<MT2HTRegion> HTRegions = data->getHTRegions();
//  std::set<MT2SignalRegion> signalRegions = data->getSignalRegions();

  std::set<MT2Region> MT2Regions = data->getRegions();

//  for( std::set<MT2HTRegion>::iterator iHT = HTRegions.begin(); iHT!=HTRegions.end(); ++iHT ) {
//    for( std::set<MT2SignalRegion>::iterator iSR = signalRegions.begin(); iSR!=signalRegions.end(); ++iSR ) {

      for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

	MT2Region thisRegion( (*iMT2) );
	
	TH1D* h1_data = data->get(thisRegion)->yield;

	for( int ibin=1; ibin<h1_data->GetXaxis()->GetNbins()+1; ++ibin ) {

	  int poisson_data = rand.Poisson(h1_data->GetBinContent(ibin, 0, 0));
	  h1_data->SetBinContent(ibin, 0, 0, poisson_data);
	  h1_data->SetBinError(ibin, 0, 0,  0.);
	  
	}  // for bins

      }// for MT2 regions

//    }// for signal regions
//  }// for HT regions

}


int matchPartonToJet( int index, MT2Tree* myTree ) {

  float deltaRMin = 0.3; // at least 0.3
  int foundId = 0;

  TLorentzVector jet;
  jet.SetPtEtaPhiM( myTree->jet_pt[index], myTree->jet_eta[index], myTree->jet_phi[index], myTree->jet_mass[index] );


  for( int i=0; i<myTree->ngenPart; ++i ) {

    if( myTree->genPart_status[i]!=23 ) continue;
    if( !(myTree->genPart_pdgId[i]==21 || abs(myTree->genPart_pdgId[i]<6)) ) continue;

    TLorentzVector thisPart;
    thisPart.SetPtEtaPhiM( myTree->genPart_pt[i], myTree->genPart_eta[i], myTree->genPart_phi[i], myTree->genPart_mass[i] );

    float thisDeltaR = jet.DeltaR(thisPart);
    if( thisDeltaR<deltaRMin ) {
      deltaRMin = thisDeltaR;
      foundId = myTree->genPart_pdgId[i];
    }
  }

  return foundId;

}
