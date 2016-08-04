#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"





bool monojet=false;


int round(float d) {
  return (int)(floor(d + 0.5));
}


void computeYield( const MT2Sample& sample, const MT2Config& cfg,  MT2Analysis<MT2EstimateTree>* anaTree );




int main( int argc, char* argv[] ) {



  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|          Running qcdControlRegion         |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./qcdControlRegion [configFileName] [data/MC/all] [monojet=false]" << std::endl;
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
  }

  if( onlyData ) {
    std::cout << "-> Will run only on data." << std::endl;
  } else if( onlyMC ) {
    std::cout << "-> Will run only on MC." << std::endl;
  } else {
    std::cout << "-> Will run on both data and MC." << std::endl;
  }
 


  if( argc > 3 ) {
    std::string monojetstr(argv[3]);
    if( monojetstr=="true" || monojetstr=="monojet" ) {
      monojet=true;
      std::cout << "-> Monojet analysis (so releasing the MT2>50. cut)" << std::endl;
    }
  }



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = cfg.getEventYieldDir() + "/qcdControlRegion"; 
  system(Form("mkdir -p %s", outputdir.c_str()));


  if( cfg.useMC() && !onlyData ) { // run on MC

    std::string samplesFile = "../samples/samples_" + cfg.qcdMCSamples() + ".dat";
    if(monojet) samplesFile = "../samples/samples_" + cfg.qcdMonoJetMCSamples() + ".dat";
    
    std::vector<MT2Sample> samples_zinv = MT2Sample::loadSamples(samplesFile, 602, 699);
    std::vector<MT2Sample> samples_wjet = MT2Sample::loadSamples(samplesFile, 502, 599);
    std::vector<MT2Sample> samples_top  = MT2Sample::loadSamples(samplesFile, 300, 399); // ignore single top and rares: faster
    //std::vector<MT2Sample> samples_top  = MT2Sample::loadSamples(samplesFile, 300, 499);
    std::vector<MT2Sample> samples_qcd  = MT2Sample::loadSamples(samplesFile, 100, 199);


    MT2Analysis<MT2EstimateTree>* qcdCRtree = new MT2Analysis<MT2EstimateTree>( "qcdCRtree", "13TeV_inclusive" );
    MT2EstimateTree::addVar( qcdCRtree, "jet1_pt" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_pt" );
    MT2EstimateTree::addVar( qcdCRtree, "caloMet" );
    MT2EstimateTree::addVar( qcdCRtree, "met_phi" );
    MT2EstimateTree::addVar( qcdCRtree, "jet1_eta" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_eta" );
    MT2EstimateTree::addVar( qcdCRtree, "jet1_phi" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_phi" );


    MT2EstimateTree::addVar( qcdCRtree, "jet1_chHEF" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_chHEF" );

    MT2EstimateTree::addVar( qcdCRtree, "jet1_neHEF" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_neHEF" );

    MT2EstimateTree::addVar( qcdCRtree, "jet1_phEF" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_phEF" );

    MT2EstimateTree::addVar( qcdCRtree, "jet1_muEF" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_muEF" );
    
    MT2EstimateTree::addVar( qcdCRtree, "jet1_eEF" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_eEF" );
    
    MT2EstimateTree::addVar( qcdCRtree, "jet1_energy" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_energy" );

    MT2EstimateTree::addVar( qcdCRtree, "badFilter" );
    MT2EstimateTree::addVar( qcdCRtree, "caloMetFilter" );
    
    MT2EstimateTree::addVar( qcdCRtree, "jet2_id" );
 

    for( unsigned i=0; i<samples_zinv.size(); ++i ) 
      computeYield( samples_zinv[i], cfg, qcdCRtree );
    for( unsigned i=0; i<samples_wjet.size(); ++i ) 
      computeYield( samples_wjet[i], cfg, qcdCRtree );
    for( unsigned i=0; i<samples_top.size(); ++i ) 
      computeYield( samples_top[i], cfg, qcdCRtree );
    for( unsigned i=0; i<samples_qcd.size(); ++i ) 
      computeYield( samples_qcd[i], cfg, qcdCRtree );
    

   
    std::string mcFile = outputdir + "/mc";
    if( monojet ) mcFile = mcFile + "_forMonojet";
    mcFile = mcFile + ".root";
    qcdCRtree->writeToFile( mcFile, "RECREATE" );

  }


  if( !(cfg.dummyAnalysis()) && cfg.qcdDataSamples()!="" && !onlyMC  ) {

    std::string samplesFile_data = "../samples/samples_" + cfg.qcdDataSamples() + ".dat";
    if(monojet) samplesFile_data = "../samples/samples_" + cfg.qcdMonoJetDataSamples() + ".dat";

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, 1, 3 );
    if( samples_data.size()==0 ) {
      std::cout << "There must be an error: samples_data is empty!" << std::endl;
      exit(1209);
    }

    MT2Analysis<MT2EstimateTree>* data = new MT2Analysis<MT2EstimateTree>( "qcdCRtree", "13TeV_inclusive" );
    MT2EstimateTree::addVar( data, "jet1_pt" );
    MT2EstimateTree::addVar( data, "jet2_pt" );
    MT2EstimateTree::addVar( data, "jet1_eta" );
    MT2EstimateTree::addVar( data, "jet2_eta" );
    MT2EstimateTree::addVar( data, "jet1_phi" );
    MT2EstimateTree::addVar( data, "jet2_phi" );
    MT2EstimateTree::addVar( data, "caloMet" );
    MT2EstimateTree::addVar( data, "met_phi" );

    MT2EstimateTree::addVar( data, "jet1_chHEF" );
    MT2EstimateTree::addVar( data, "jet2_chHEF" );

    MT2EstimateTree::addVar( data, "jet1_neHEF" );
    MT2EstimateTree::addVar( data, "jet2_neHEF" );

    MT2EstimateTree::addVar( data, "jet1_phEF" );
    MT2EstimateTree::addVar( data, "jet2_phEF" );

    MT2EstimateTree::addVar( data, "jet1_muEF" );
    MT2EstimateTree::addVar( data, "jet2_muEF" );

    MT2EstimateTree::addVar( data, "jet1_eEF" );
    MT2EstimateTree::addVar( data, "jet2_eEF" );

    MT2EstimateTree::addVar( data, "jet1_energy" );
    MT2EstimateTree::addVar( data, "jet2_energy" );

    MT2EstimateTree::addVar( data, "badFilter" );
    MT2EstimateTree::addVar( data, "caloMetFilter" );
 
    MT2EstimateTree::addVar( data, "jet2_id" );
    
 

    //MT2Analysis<MT2EstimateTree>* data = new MT2Analysis<MT2EstimateTree>( "qcdCRtree", cfg.regionsSet() );
    for( unsigned i=0; i < samples_data.size(); ++i )
      computeYield( samples_data[i], cfg, data );

    std::string dataFile = outputdir + "/data";
    if( monojet ) dataFile = dataFile + "_forMonojet";
    dataFile = dataFile + ".root";
    data->writeToFile( dataFile, "RECREATE" );


  }



  return 0;

}








void computeYield( const MT2Sample& sample, const MT2Config& cfg, MT2Analysis<MT2EstimateTree>* anaTree ) {


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

    if( myTree.isData ) {
      if (  myTree.isGolden == 0 ) continue;
      if ( !myTree.passFilters() ) continue;
      //  if( !(myTree.Flag_badMuonFilter && myTree.Flag_badChargedHadronFilter) ) continue;
    }
    else {
      if ( !(myTree.Flag_badMuonFilter>0 && myTree.Flag_badChargedHadronFilter>0 && myTree.Flag_EcalDeadCellTriggerPrimitiveFilter>0) ) continue;
      if (myTree.met_pt/myTree.met_caloPt > 5.0) continue; // RA2 filter for QCD MC
    }
    

    if( !myTree.passSelection("qcd") ) continue;


    float minMTBmet = myTree.minMTBMet;
    //float met       = myTree.met_pt;
    int njets       = myTree.nJet30;
    int nbjets      = myTree.nBJet20;    
    //int nbjets      = (sample.id>=600&&sample.id<=610) ? myTree.nBJet20 : myTree.nBJet20csv;
    //myTree.nBJet20 = nbjets; // this is becuase we are still using 2015 Zinv 
    float mt2       = (njets>1) ? myTree.mt2 : myTree.jet1_pt;
    float ht        = myTree.ht;

    if (myTree.isData) {

      int id = sample.id;
      // sample IDs for data:
      // JetHT = 1
      // HTMHT = 2
      // MET   = 3
      //myTree.evt_id = id; // useful when using american trees

      if( njets==1 ) {

        if( !( id==3 && myTree.HLT_PFMET100_PFMHT100) ) continue;

      } else { // njets>=2


	if( !monojet ){
	  if( ht>1000. ) {
	    if( !( id==1 && myTree.HLT_PFHT800) ) continue;
	  } else if( ht>575. ) {
	    if( !( (id==2 && myTree.HLT_PFHT300_PFMET100 ) || (id==1 && myTree.HLT_PFHT475_Prescale))  ) continue;
	  } else if( ht>450. ) {
	    if( !( (id==2 && myTree.HLT_PFHT300_PFMET100 ) || (id==1 && myTree.HLT_PFHT350_Prescale))  ) continue;
	  } else if( ht>200. ) {
	    if( !( (id==3 && myTree.HLT_PFMET100_PFMHT100) || (id==1 && (myTree.HLT_PFHT125_Prescale)))  ) continue;
	  }
	}else{
	  if( !( (id==1 && myTree.HLT_PFHT800) || (id==2 &&  myTree.HLT_PFHT300_PFMET100 && !myTree.HLT_PFHT800 ) || ( id==3 && myTree.HLT_PFMET100_PFMHT100 && !myTree.HLT_PFHT800 && !myTree.HLT_PFHT300_PFMET100) ) ) continue;
	}
      }

    } // if is data


    if( monojet ) {
      if( !( (njets==2 && myTree.deltaPhiMin<0.3 && myTree.jet1_pt>200. && myTree.met_pt>200.) ) ) continue;
      if( !myTree.passMonoJetId(0) ) continue;
      //if( !myTree.passMonoJetId(1) ) continue;

      //Fix for monojetCR
      if( !myTree.Flag_EcalDeadCellTriggerPrimitiveFilter ) continue;

      if( !myTree.isData && (myTree.met_pt/ myTree.met_caloPt > 5.0) ) continue;
      // if( (myTree.met_pt/ myTree.met_caloPt > 5.0) ) continue;

    } else {
      if( mt2<50. ) continue;
    }

    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi(); 

    if(!myTree.isData){
      weight *= myTree.weight_btagsf;
      weight *= myTree.weight_lepsf;
    }
    //float myht = njets==1 ? 201. : ht; // let everything (mt2=ht>40) pass for monojet
    MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisTree==0 ) continue;


    thisTree->yield->Fill( mt2, weight );
    thisTree->assignVar( "jet1_pt", myTree.jet1_pt );
    thisTree->assignVar( "jet2_pt", myTree.jet2_pt );

    double jet1_energy= -99; double  jet2_energy = -99;

    double jet1_eta= -99; double  jet2_eta = -99;
    double jet1_phi= -99; double  jet2_phi= -99;
    double jet1_chHEF= -99; double  jet2_chHEF = -99;
    double jet1_neHEF= -99; double  jet2_neHEF = -99;
    double jet1_phEF= -99; double  jet2_phEF = -99;
    double jet1_muEF= -99; double  jet2_muEF = -99;
    double jet1_eEF= -99; double  jet2_eEF = -99;


 
    bool isLeading=false;
    bool isSubLeading=false;
    for(int i=0; i<myTree.njet; i++){
      
      if( fabs(myTree.jet_eta[i]) < 2.5 ){

      TLorentzVector *LVec = new TLorentzVector;

        if( isLeading==false ){
          jet1_eta = myTree.jet_eta[i];
          jet1_phi = myTree.jet_phi[i];

          jet1_chHEF = myTree.jet_chHEF[i];
          jet1_neHEF = myTree.jet_neHEF[i];
          jet1_phEF = myTree.jet_phEF[i];
          jet1_muEF = myTree.jet_muEF[i];
	  jet1_eEF = myTree.jet_eEF[i];

          LVec->SetPtEtaPhiM(myTree.jet_pt[i], myTree.jet_eta[i],myTree.jet_phi[i], myTree.jet_mass[i]);

	  jet1_energy = LVec->E();

	  isLeading=true;

	}else if( isSubLeading==false ){
	  jet2_eta = myTree.jet_eta[i];
	  jet2_phi = myTree.jet_phi[i];

	  jet2_chHEF = myTree.jet_chHEF[i];
	  jet2_neHEF = myTree.jet_neHEF[i];
	  jet2_phEF = myTree.jet_phEF[i];
	  jet2_muEF = myTree.jet_muEF[i];
	  jet2_eEF = myTree.jet_eEF[i];

          LVec->SetPtEtaPhiM(myTree.jet_pt[i], myTree.jet_eta[i],myTree.jet_phi[i], myTree.jet_mass[i]);

	  jet2_energy = LVec->E();

	  isSubLeading=true;

	}else continue;
	
      }else continue;
    }


    thisTree->assignVar( "jet1_eta", jet1_eta );
    thisTree->assignVar( "jet2_eta", jet2_eta );
    thisTree->assignVar( "jet1_phi", jet1_phi );
    thisTree->assignVar( "jet2_phi", jet2_phi );
    thisTree->assignVar( "caloMet", myTree.met_caloPt );
    thisTree->assignVar( "met_phi", myTree.met_phi );


    thisTree->assignVar( "jet1_chHEF", jet1_chHEF );
    thisTree->assignVar( "jet2_chHEF", jet2_chHEF );

    thisTree->assignVar( "jet1_neHEF", jet1_neHEF );
    thisTree->assignVar( "jet2_neHEF", jet2_neHEF );

    thisTree->assignVar( "jet1_phEF", jet1_phEF );
    thisTree->assignVar( "jet2_phEF", jet2_phEF );

    thisTree->assignVar( "jet1_muEF", jet1_muEF );
    thisTree->assignVar( "jet2_muEF", jet2_muEF );

    thisTree->assignVar( "jet1_eEF", jet1_eEF );
    thisTree->assignVar( "jet2_eEF", jet2_eEF );

    thisTree->assignVar( "jet1_energy", jet1_energy );
    thisTree->assignVar( "jet2_energy", jet2_energy );

    thisTree->assignVar( "badFilter", myTree.Flag_badMuonFilter && myTree.Flag_badChargedHadronFilter );

    thisTree->assignVar( "caloMetFilter", (myTree.met_pt/ myTree.met_caloPt) < 5.0 );
    


    bool jet2_id =  (sample.id >=600 &&sample.id<700) ?  !(myTree.jet_id[2]>=3 && myTree.jet_chHEF[2]>0.05 && myTree.jet_neHEF[2]<0.8 && myTree.jet_phEF[2]<0.7) :  myTree.passMonoJetId(1)  ;   
    thisTree->assignVar( "jet2_id",  jet2_id );


    // jet_chHEF[j]>0.05 && jet_neHEF[j]<0.8 && jet_phEF[j]<0.7;

    thisTree->fillTree( myTree, weight );

    
  } // for entries


  anaTree->finalize();


  delete tree;


  file->Close();
  delete file;
  

}



