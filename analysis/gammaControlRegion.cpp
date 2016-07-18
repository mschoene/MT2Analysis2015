#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2EstimateTree.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"
#include "TF1.h"


bool alsoSignals = false;


int round(float d) {
  return (int)(floor(d + 0.5));
}


void computeYield( const MT2Sample& sample, const MT2Config& cfg, 
                   MT2Analysis<MT2EstimateTree>* anaTree,
                   MT2Analysis<MT2EstimateTree>* anaTree_pass,
                   MT2Analysis<MT2EstimateZinvGamma>* prompt=0, MT2Analysis<MT2EstimateZinvGamma>* prompt_pass=0, 
                   MT2Analysis<MT2EstimateZinvGamma>* nip=0, MT2Analysis<MT2EstimateZinvGamma>* nip_pass=0, 
                   MT2Analysis<MT2EstimateZinvGamma>* fake=0, MT2Analysis<MT2EstimateZinvGamma>* fake_pass=0 );
void roundLikeData( MT2Analysis<MT2EstimateZinvGamma>* data );
void fillYields( MT2Analysis<MT2EstimateZinvGamma>* est, float weight, float ht, int njets, int nbjets, float met, float minMTBmet, float mt2, float iso );
void fillOneTree( MT2EstimateTree* thisTree, const MT2Tree& myTree, float weight, float ht, int njets, int nbjets, float met, float minMTBmet, float mt2, float iso, int nTrueB, int nTrueC );

float DeltaR(float eta1, float eta2, float phi1, float phi2);

float DeltaPhi(float phi1, float phi2);


int main( int argc, char* argv[] ) {



  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|            Running gammaControlRegion              |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./gammaControlRegion [configFileName] [data/MC]" << std::endl;
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
      std::cout << "-> You passed a second argument that isn't 'data' or 'MC', so I don't know what to do about it." << std::endl;
    }
  }



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = cfg.getGammaCRdir();
  system(Form("mkdir -p %s", outputdir.c_str()));


  if( cfg.useMC() && !onlyData ) { // run on MC

    std::string samplesFile = "../samples/samples_" + cfg.mcSamples() + ".dat";
    
    //std::vector<MT2Sample> samples_gammaJet = MT2Sample::loadSamples(samplesFile, "gjet");
    //std::vector<MT2Sample> samples_gammaJet = MT2Sample::loadSamples(samplesFile, 200, 299);
    std::vector<MT2Sample> samples_gammaJet = MT2Sample::loadSamples(samplesFile, "GJets");
    if( samples_gammaJet.size()==0 ) {
      std::cout << "There must be an error: didn't find any gamma+jet files in " << samplesFile << "!" << std::endl;
      exit(1209);
    }

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading QCD samples" << std::endl;

    //std::vector<MT2Sample> samples_qcd = MT2Sample::loadSamples(samplesFile, "qcd");
    std::vector<MT2Sample> samples_qcd = MT2Sample::loadSamples(samplesFile, "QCD");


    MT2Analysis<MT2EstimateTree>* tree = new MT2Analysis<MT2EstimateTree>( "gammaCRtree_loose", cfg.crRegionsSet() );
    MT2EstimateTree::addVar( tree, "prompt" );
    MT2EstimateTree::addVar( tree, "iso" );
    //MT2EstimateTree::addVar( tree, "isoRC" );
    MT2EstimateTree::addVar( tree, "sietaieta" );
    MT2EstimateTree::addVar( tree, "ptGamma" );
    MT2EstimateTree::addVar( tree, "etaGamma" );
    MT2EstimateTree::addVar( tree, "jet1_pt" );
    MT2EstimateTree::addVar( tree, "jet2_pt" );
    MT2EstimateTree::addVar( tree, "drMinParton" );

    
    MT2Analysis<MT2EstimateTree>* tree_pass = new MT2Analysis<MT2EstimateTree>( "gammaCRtree", cfg.crRegionsSet() );
    MT2EstimateTree::addVar( tree_pass, "prompt" );
    MT2EstimateTree::addVar( tree_pass, "iso" );
    //MT2EstimateTree::addVar( tree_pass, "isoRC" );
    MT2EstimateTree::addVar( tree_pass, "sietaieta" );
    MT2EstimateTree::addVar( tree_pass, "ptGamma" );
    MT2EstimateTree::addVar( tree_pass, "etaGamma" );
    MT2EstimateTree::addVar( tree_pass, "jet1_pt" );
    MT2EstimateTree::addVar( tree_pass, "jet2_pt" );
    MT2EstimateTree::addVar( tree_pass, "drMinParton" );
  

    if(cfg.analysisType() == "ZG"){
      MT2EstimateTree::addVar( tree, "raw_mt2" );
      MT2EstimateTree::addVar( tree_pass, "raw_mt2" );
    }
      
    
    MT2EstimateTree::addVar( tree, "gamma_chHadIsoRC" );
    MT2EstimateTree::addVar( tree_pass, "gamma_chHadIsoRC" );
    
    MT2Analysis<MT2EstimateZinvGamma>* prompt = new MT2Analysis<MT2EstimateZinvGamma>( "prompt", cfg.regionsSet() );
    MT2Analysis<MT2EstimateZinvGamma>* prompt_pass = new MT2Analysis<MT2EstimateZinvGamma>( "prompt_pass", cfg.regionsSet() );

    MT2Analysis<MT2EstimateZinvGamma>* fake = new MT2Analysis<MT2EstimateZinvGamma>( "fake", cfg.regionsSet() );
    MT2Analysis<MT2EstimateZinvGamma>* fake_pass = new MT2Analysis<MT2EstimateZinvGamma>( "fake_pass", cfg.regionsSet() );

    MT2Analysis<MT2EstimateZinvGamma>* nip = new MT2Analysis<MT2EstimateZinvGamma>( "nip", cfg.regionsSet() );
    MT2Analysis<MT2EstimateZinvGamma>* nip_pass = new MT2Analysis<MT2EstimateZinvGamma>( "nip_pass", cfg.regionsSet() );

    for( unsigned i=0; i<samples_gammaJet.size(); ++i ) {
      computeYield( samples_gammaJet[i], cfg, tree, tree_pass, prompt, prompt_pass, nip, nip_pass, fake, fake_pass );
    }

    for( unsigned i=0; i<samples_qcd.size(); ++i ) {
      computeYield( samples_qcd[i], cfg, tree, tree_pass, prompt, prompt_pass, nip, nip_pass, fake, fake_pass );
    }


    MT2Analysis<MT2EstimateZinvGamma>* gammaCR_loose = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR_loose",  cfg.regionsSet() );
    (*gammaCR_loose) = (*prompt) + (*nip) + (*fake);
   
   
    MT2Analysis<MT2EstimateZinvGamma>* gammaCR         = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR", cfg.regionsSet() );
    MT2Analysis<MT2EstimateZinvGamma>* gammaCR_nipUp   = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR_nipUp", cfg.regionsSet() );
    MT2Analysis<MT2EstimateZinvGamma>* gammaCR_nipDown = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR_nipDown", cfg.regionsSet() );
    (*gammaCR        ) = (*prompt_pass) + (*fake_pass) +     (*nip_pass);
    (*gammaCR_nipUp  ) = (*prompt_pass) + (*fake_pass) + 2. *(*nip_pass);
    (*gammaCR_nipDown) = (*prompt_pass) + (*fake_pass) + 0.5*(*nip_pass);
    
   
    MT2Analysis<MT2EstimateZinvGamma>* matched = new MT2Analysis<MT2EstimateZinvGamma>( "matched", cfg.regionsSet() ); 
    (*matched) = (*prompt) + (*nip);

    MT2Analysis<MT2EstimateZinvGamma>* matched_pass = new MT2Analysis<MT2EstimateZinvGamma>( "matched_pass", cfg.regionsSet() ); 
    (*matched_pass) = (*prompt_pass) + (*nip_pass);
  

    MT2Analysis<MT2EstimateSyst>* f = MT2EstimateSyst::makeEfficiencyAnalysis( "f", (MT2Analysis<MT2Estimate>*)prompt, (MT2Analysis<MT2Estimate>*)matched );
    MT2Analysis<MT2EstimateSyst>* f_pass = MT2EstimateSyst::makeEfficiencyAnalysis( "f_pass", (MT2Analysis<MT2Estimate>*)prompt_pass, (MT2Analysis<MT2Estimate>*)matched_pass );



    MT2Analysis<MT2EstimateSyst>* eff = MT2EstimateSyst::makeEfficiencyAnalysis( "eff", (MT2Analysis<MT2Estimate>*)prompt_pass, (MT2Analysis<MT2Estimate>*)prompt );

    MT2Analysis<MT2EstimateSyst>* purityTight = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", (MT2Analysis<MT2Estimate>*)matched_pass, (MT2Analysis<MT2Estimate>*)gammaCR);

    MT2Analysis<MT2EstimateSyst>* purityLoose = MT2EstimateSyst::makeEfficiencyAnalysis( "purityLoose", (MT2Analysis<MT2Estimate>*)matched, (MT2Analysis<MT2Estimate>*)gammaCR_loose );


    std::string mcFile = outputdir + "/mc.root";

    gammaCR_loose->writeToFile( mcFile, "RECREATE" );
    gammaCR      ->writeToFile( mcFile );
    prompt       ->writeToFile( mcFile );
    fake         ->writeToFile( mcFile );
    nip          ->writeToFile( mcFile );
    prompt_pass  ->writeToFile( mcFile );
    fake_pass    ->writeToFile( mcFile );
    nip_pass     ->writeToFile( mcFile );
    f            ->writeToFile( mcFile );
    f_pass       ->writeToFile( mcFile );
    tree         ->writeToFile( mcFile );
    tree_pass    ->writeToFile( mcFile );



    std::vector< MT2Analysis< MT2EstimateTree>* > signals;
    if( alsoSignals ) {

      std::cout << std::endl << std::endl;
      std::cout << "-> Loading signal samples" << std::endl;
      std::string samplesFile_signals = "samples_signalsSNT.txt";
      std::vector<MT2Sample> samples_sig = MT2Sample::loadSamples(samplesFile_signals,1000 );

      for( unsigned i=0; i<samples_sig.size(); ++i ) {

        MT2Analysis<MT2EstimateTree>* thisSig = new MT2Analysis<MT2EstimateTree>( samples_sig[i].sname + "_loose", cfg.regionsSet(), samples_sig[i].id );
        MT2EstimateTree::addVar( thisSig, "prompt" );
        MT2EstimateTree::addVar( thisSig, "iso" );

        MT2Analysis<MT2EstimateTree>* thisSig_pass = new MT2Analysis<MT2EstimateTree>( samples_sig[i].sname, cfg.regionsSet(), samples_sig[i].id );
        MT2EstimateTree::addVar( thisSig_pass, "prompt" );
        MT2EstimateTree::addVar( thisSig_pass, "iso" );

        computeYield( samples_sig[i], cfg, thisSig, thisSig_pass );

        signals.push_back( thisSig_pass );

      }  // for signals

    }

    for( unsigned i=0; i<signals.size(); ++i )
      signals[i]->writeToFile( mcFile);

    //comment here again
    if( cfg.dummyAnalysis() ) {

      // emulate data:
      roundLikeData(gammaCR);
      roundLikeData(gammaCR_loose);
      gammaCR->writeToFile( outputdir + "/data.root" );
      gammaCR_loose->writeToFile( outputdir + "/data.root" );
      gammaCR_nipUp->writeToFile( outputdir + "/data.root" );
      gammaCR_nipDown->writeToFile( outputdir + "/data.root" );
      tree->writeToFile( outputdir + "/data.root" );
      tree_pass->writeToFile( outputdir + "/data.root" );

    } // if dummy

    purityTight->writeToFile( outputdir + "/purityMC.root", "RECREATE" );
    purityLoose->writeToFile( outputdir + "/purityMC.root" );
    eff->writeToFile( outputdir + "/purityMC.root" );

  } // if run on MC



  if( !onlyMC && !cfg.dummyAnalysis() ) {

    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";
 
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "SinglePhoton");

    if( samples_data.size()==0 ) {

      std::cout << std::endl;
      std::cout << "-> WARNING!! Didn't find any data in file: " << samplesFile_data << "!" << std::endl;
      std::cout << "-> Exiting." << std::endl;
      std::cout << std::endl;

    } else {


      MT2Analysis<MT2EstimateZinvGamma>* dataCR_loose = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR_loose",  cfg.regionsSet() );
      MT2Analysis<MT2EstimateZinvGamma>* dataCR       = new MT2Analysis<MT2EstimateZinvGamma>( "gammaCR",  cfg.regionsSet() );

      MT2Analysis<MT2EstimateTree>* tree = new MT2Analysis<MT2EstimateTree>( "gammaCRtree_loose", cfg.crRegionsSet() );
      MT2EstimateTree::addVar( tree, "iso" );
      MT2EstimateTree::addVar( tree, "sietaieta" );
      MT2EstimateTree::addVar( tree, "ptGamma" );
      MT2EstimateTree::addVar( tree, "etaGamma" );
      MT2EstimateTree::addVar( tree, "jet1_pt" );
      MT2EstimateTree::addVar( tree, "jet2_pt" );
      MT2EstimateTree::addVar( tree, "drMinParton" );
      
      MT2Analysis<MT2EstimateTree>* tree_pass = new MT2Analysis<MT2EstimateTree>( "gammaCRtree", cfg.crRegionsSet() );
      MT2EstimateTree::addVar( tree_pass, "iso" );
      MT2EstimateTree::addVar( tree_pass, "sietaieta" );
      MT2EstimateTree::addVar( tree_pass, "ptGamma" );
      MT2EstimateTree::addVar( tree_pass, "etaGamma" );
      MT2EstimateTree::addVar( tree_pass, "jet1_pt" );
      MT2EstimateTree::addVar( tree_pass, "jet2_pt" );
      MT2EstimateTree::addVar( tree_pass, "drMinParton" );

      MT2EstimateTree::addVar( tree, "gamma_chHadIsoRC" );
      MT2EstimateTree::addVar( tree_pass, "gamma_chHadIsoRC" );

      if( cfg.additionalStuff()=="hfContent" ) {
        MT2EstimateTree::addVar( tree, "nTrueB" );
        MT2EstimateTree::addVar( tree, "nTrueC" );

        MT2EstimateTree::addVar( tree_pass, "nTrueB" );
        MT2EstimateTree::addVar( tree_pass, "nTrueC" );
      }

 
      if(cfg.analysisType() == "ZG"){
	MT2EstimateTree::addVar( tree, "raw_mt2" );
	MT2EstimateTree::addVar( tree_pass, "raw_mt2" );
      }
      

      for( unsigned i=0; i<samples_data.size(); ++i ) {
        computeYield( samples_data[i], cfg, tree, tree_pass, dataCR_loose, dataCR );
      }

      tree        ->writeToFile( outputdir + "/data.root", "RECREATE" );
      tree_pass   ->writeToFile( outputdir + "/data.root" );
      dataCR_loose->writeToFile( outputdir + "/data.root" );
      dataCR      ->writeToFile( outputdir + "/data.root" );

    } // if data samples > 0

  } // if run on data


  return 0;

}








void computeYield( const MT2Sample& sample, const MT2Config& cfg,
                   MT2Analysis<MT2EstimateTree>* anaTree,
                   MT2Analysis<MT2EstimateTree>* anaTree_pass,
                   MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* prompt_pass, 
                   MT2Analysis<MT2EstimateZinvGamma>* nip, MT2Analysis<MT2EstimateZinvGamma>* nip_pass, 
                   MT2Analysis<MT2EstimateZinvGamma>* fake, MT2Analysis<MT2EstimateZinvGamma>* fake_pass ) {



  float isoCut = cfg.gammaIsoCut();  

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;


  TF1* f1_ratio_2b1b = 0;

  if( cfg.gamma2bMethod()=="2b1bRatio" ) {
    TFile* ratioFile = TFile::Open( Form( "%s/2bRatio/mc.root", cfg.getEventYieldDir().c_str() ) );
    if( ratioFile==0 ) {
      std::cout << "-> Didn't find 2b/1b ratio file. Please produce it first with computeZinv2b" << std::endl;
      exit(1919191);
    }
    TH1D* h1_ratio = (TH1D*)ratioFile->Get("r_vs_nJets_zllMC");
    f1_ratio_2b1b = h1_ratio->GetFunction("line");
  }


  bool isQCD  = sample.id>=100 && sample.id<200;

  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  if( cfg.additionalStuff()=="qgVars" || cfg.additionalStuff()=="hfContent" ) 
     myTree.loadGenStuff = true;
  myTree.Init(tree);

  int nentries = tree->GetEntries();


  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    //    if( myTree.isData )
    //      if ( myTree.isGolden == 0 ) continue;

    if(cfg.analysisType() == "mt2"){
   
      if( !myTree.passSelection("gamma") ) continue;

      if( myTree.mt2>200. ) continue; // orthogonal to signal region
      if( myTree.gamma_pt[0]<180. ) continue;
      if( (myTree.gamma_nJet30>1 && myTree.gamma_mt2<200.) || (myTree.gamma_nJet30==1 && myTree.gamma_ht<200.) ) continue;
    
    }

    // REMOVE THIS SOOOOON
    // or maybe not due to spiky behavior in qcd samples
    
// if( myTree.evt_scale1fb>10. ) continue;


//    if( myTree.gamma_nJet30==1 ){
//      
//      float maxDR=0;
//      int J=0;
//      for( int j=0; j<myTree.njet; ++j ){
//
//	if( fabs( myTree.jet_eta[j] ) > 2.5 || myTree.jet_pt[j] < 30. || j>1) continue;
//	
//	float thisDR = DeltaR( myTree.gamma_eta[0], myTree.jet_eta[j], myTree.gamma_phi[0], myTree.jet_phi[j] );
//	maxDR = ( thisDR > maxDR ) ? thisDR : maxDR; 
//	J = ( thisDR > maxDR ) ? j : J;
//
//      }
//      
//      if(!(myTree.passMonoJetId(J)));
//    
//    }
  
    if( myTree.isData ) {
      if( !myTree.passGammaAdditionalSelection(1) ) continue;
    } else {
      if( !myTree.passGammaAdditionalSelection(sample.id) ) continue;
    }
    
  
   
    if( !( myTree.HLT_Photon165_HE10 ) ) continue; 

    if( cfg.additionalStuff()=="gammaNoSietaieta" ) {
      //if( myTree.gamma_hOverE[0]>0.1 ) continue;
    } else {
      if( myTree.gamma_idCutBased[0]==0 ) continue;
    }

     
    if( myTree.isData ) {
      
      if( !( myTree.isGolden ) ) continue;


      if( !(myTree.passFilters() ) ) continue;
    
      if( !(myTree.HLT_Photon165_HE10) ) continue;
    
      //if( !(nVert>0. && Flag_HBHENoiseFilter==1 && Flag_eeBadScFilter==1 ) ) continue;
      //if( !(myTree.nVert>0. && myTree.Flag_HBHENoiseFilter==1 ) ) continue;

    }


    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );


    // absolute iso:
    float iso = myTree.gamma_chHadIso[0];
    if( cfg.additionalStuff()=="gammaNoSietaieta" ) {
      if( iso>20. ) continue; // loose preselection for plots
    } else {
      if( iso>10. ) continue; // preselection
    }

    float minMTBmet = myTree.gamma_minMTBMet;
    float met       = myTree.gamma_met_pt;
    int njets       = myTree.gamma_nJet30;
    int nbjets      = (sample.id>10) ?  myTree.gamma_nBJet20 :  myTree.gamma_nBJet20csv;    
    //int nbjets      = myTree.gamma_nBJet20csv;    
    float ht        = myTree.gamma_ht;
    float mt2       = (njets>1) ? myTree.gamma_mt2 : ht;
    //float mt2       = myTree.gamma_mt2;

    //TEEEEEMPORARYYYYY
    //    if( myTree.isData && myTree.run>275125.) continue;

    //temp fix
    myTree.nBJet20csv=nbjets;

    if( cfg.gamma2bMethod()=="2b1bRatio" && nbjets==2 )
      continue; // will take 2b from reweighted 1b so skip

    //    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi()*myTree.puWeight; 
    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi(); 
    //  Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi(); 
    
    if( !myTree.isData ){
      weight *= myTree.weight_btagsf;

      float SF = 1.0;
      float pt = myTree.gamma_pt[0];
      
      if ( fabs(myTree.gamma_eta[0])<1.479 ) {
	if (pt > 150 && pt < 160) SF = 0.03576134;
	else if (pt < 170.) SF = 0.2844761;
	else if (pt < 180.) SF = 0.9541034;
	else if (pt < 200.) SF = 0.9432587;
	else if (pt < 250.) SF = 0.9998274;
	else if (pt < 300.) SF = 0.9997936;
	else if (pt < 350.) SF = 0.9997995;
	else if (pt < 400.) SF = 0.9727327;
	else if (pt < 450.) SF = 0.9686112;
	else if (pt < 500.) SF = 0.9689014;
	else if (pt < 600.) SF = 0.965793;
	else if (pt < 700.) SF = 0.9734969;
	else if (pt < 800.) SF = 0.9776452;
	else SF = 0.9367322;
      }else {
	if (pt > 150 && pt < 160) SF = 0.02972659;
	else if (pt < 170.) SF = 0.1191435;
	else if (pt < 180.) SF = 0.5990605;
	else if (pt < 200.) SF = 0.9105484;
	else if (pt < 250.) SF = 0.9995712;
	else if (pt < 300.) SF = 0.9915838;
	else if (pt < 350.) SF = 0.9545882;
	else if (pt < 400.) SF = 0.9673124;
	else if (pt < 450.) SF = 0.960792;
	else if (pt < 500.) SF = 0.9519475;
	else if (pt < 600.) SF = 0.950766;
	else if (pt < 700.) SF = 0.9538072;
	else if (pt < 800.) SF = 0.9777396;
	else SF = 0.9862047;
      }
    
      weight *= SF;
      
    }

    bool passIso = iso<isoCut;

    MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
    MT2EstimateTree* thisTree_pass = anaTree_pass->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisTree==0 ) continue;


    MT2Analysis<MT2EstimateZinvGamma>* theEstimate = 0;
    MT2Analysis<MT2EstimateZinvGamma>* theEstimate_pass = 0;

    if( !myTree.isData ) {

      int mcMatchId = myTree.gamma_mcMatchId[0];
      bool isMatched = (mcMatchId==22 || mcMatchId==7);

      // bool isPrompt = isMatched && !isQCD;
      //bool isNIP    = isMatched && isQCD;
      //bool isFake   = !isMatched;
      bool isPrompt = isMatched && !isQCD && myTree.gamma_drMinParton[0]>0.4;
      bool isNIP    = isMatched && isQCD && myTree.gamma_drMinParton[0]<0.4;
      bool isFake   = !isMatched && isQCD;

  
      int promptLevel = -1; 
      if( isPrompt ) promptLevel = 2;
      else if( isNIP ) promptLevel = 1;
      else if( isFake ) promptLevel = 0;
      else std::cout << "WARNING!!! This photon is neither prompt, nor fragmentation nor fake!" << std::endl;

      thisTree->assignVar( "prompt", promptLevel );
      if( passIso ) thisTree_pass->assignVar( "prompt", promptLevel );


      if( promptLevel==2 ) {
        theEstimate = prompt;
        theEstimate_pass = prompt_pass;
      } else if( promptLevel==1 ) {
        theEstimate = nip;
        theEstimate_pass = nip_pass;
      } else if( promptLevel==0 ) {
        theEstimate = fake;
        theEstimate_pass = fake_pass;
      }



    } else { // this is data:

      // so nevermind that it's called prompt, it's actually the full CR
      theEstimate = prompt;
      theEstimate_pass = prompt_pass;
      
      //MT2EstimateZinvGamma* thisEstimate = prompt->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      //if( thisEstimate==0 ) continue;

      //thisEstimate->yield->Fill(mt2, weight );
      //thisEstimate->fillIso( iso, weight, mt2 );

      //if( passIso ) {

      //  MT2EstimateZinvGamma* thisEstimate_pass = prompt_pass->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      //  if( thisEstimate_pass==0 ) continue;

      //  thisEstimate_pass->yield->Fill(mt2, weight );
      //  thisEstimate_pass->fillIso( iso, weight, mt2 );

      //} // if pass iso

    }  // if is data


    fillYields( theEstimate, weight, ht, njets, nbjets, met, minMTBmet, mt2, iso );
    if( passIso ) {
      fillYields( theEstimate_pass, weight, ht, njets, nbjets, met, minMTBmet, mt2, iso );
    }


  


    int nTrueB=-1;
    int nTrueC=-1;

    if( cfg.additionalStuff()=="hfContent" ) {

      nTrueB = 0;
      nTrueC = 0;

      for( int ipart=0; ipart<myTree.ngenPart; ++ipart ) {

        if( myTree.genPart_pt[ipart] < 20. ) continue;
        if( abs(myTree.genPart_eta[ipart])>2.5 ) continue;
        if( myTree.genPart_status[ipart] != 23 ) continue;

        if( abs(myTree.genPart_pdgId[ipart])==5 )
          nTrueB++;
        if( abs(myTree.genPart_pdgId[ipart])==4 )
          nTrueC++;

      }
    }

    if(cfg.analysisType() == "ZG")
      thisTree->assignVar( "raw_mt2", myTree.mt2 );
      
    //    thisTree_pass->fillTree_gamma(myTree, weight );


    fillOneTree( thisTree, myTree, weight, ht, njets, nbjets, met, minMTBmet, mt2, iso, nTrueB, nTrueC );

    if( passIso ) 
      fillOneTree( thisTree_pass, myTree, weight, ht, njets, nbjets, met, minMTBmet, mt2, iso, nTrueB, nTrueC );



    if( cfg.gamma2bMethod()=="2b1bRatio" && nbjets==1 ) { // fill again to reweight 1b events to get 2b

      float corr = TMath::Max( f1_ratio_2b1b->Eval(njets), 0. );

      fillYields( theEstimate, corr*weight, ht, njets, 2, met, minMTBmet, mt2, iso );
      if( passIso ) {
        fillYields( theEstimate_pass, corr*weight, ht, njets, 2, met, minMTBmet, mt2, iso );
      }

 

      MT2EstimateTree* thisTree_2b = anaTree->get( ht, njets, 2, minMTBmet, mt2 );
      MT2EstimateTree* thisTree_2b_pass = anaTree_pass->get( ht, njets, 2, minMTBmet, mt2 );
      if( thisTree_2b==0 ) continue;

      fillOneTree( thisTree_2b, myTree, corr*weight, ht, njets, 2, met, minMTBmet, mt2, iso, nTrueB, nTrueC );
      if( passIso ) 
        fillOneTree( thisTree_2b_pass, myTree, corr*weight, ht, njets, 2, met, minMTBmet, mt2, iso, nTrueB, nTrueC );
      
    } // if ratio method

    //    std::cout << "ht " << ht << std::endl;
    
    
  } // for entries


  if( prompt!=0 ) prompt->finalize();
  if( prompt_pass!=0 ) prompt_pass->finalize();
  if( nip!=0 ) nip->finalize();
  if( nip_pass!=0 ) nip_pass->finalize();
  if( fake!=0 ) fake->finalize();
  if( fake_pass!=0 ) fake_pass->finalize();
  anaTree->finalize();
  anaTree_pass->finalize();
  

  delete tree;


  file->Close();
  delete file;
  

}



void roundLikeData( MT2Analysis<MT2EstimateZinvGamma>* data ) {



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



void fillOneTree( MT2EstimateTree* thisTree, const MT2Tree& myTree, float weight, float ht, int njets, int nbjets, float met, float minMTBmet, float mt2, float iso, int nTrueB, int nTrueC ) {

  thisTree->yield->Fill(mt2, weight );
  thisTree->assignVar( "iso", iso );
  thisTree->assignVar( "sietaieta", myTree.gamma_sigmaIetaIeta[0] );
  thisTree->assignVar( "ptGamma", myTree.gamma_pt[0] );
  thisTree->assignVar( "etaGamma", myTree.gamma_eta[0] );
  thisTree->assignVar( "jet1_pt", myTree.gamma_jet1_pt );
  thisTree->assignVar( "jet2_pt", myTree.gamma_jet2_pt );
  thisTree->assignVar( "gamma_chHadIsoRC",  myTree.gamma_chHadIsoRC[0] );
  thisTree->assignVar( "drMinParton", myTree.gamma_drMinParton[0] );


  if( nTrueB>=0 ) {

    thisTree->assignVar( "nTrueB", nTrueB );
    thisTree->assignVar( "nTrueC", nTrueC );

  }

  thisTree->assignTree_gamma(myTree, weight ); // assign all to defaults
  thisTree->assignVars( ht, njets, nbjets, met, mt2 ); //manually override these (so change nbjets)
  thisTree->tree->Fill();

}


void fillYields( MT2Analysis<MT2EstimateZinvGamma>* est, float weight, float ht, int njets, int nbjets, float met, float minMTBmet, float mt2, float iso ) {

  if( est!=0 ) {

    MT2EstimateZinvGamma* thisEst = est->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisEst==0 ) return;

    thisEst->yield->Fill(mt2, weight );
    thisEst->fillIso( iso, weight, mt2 );

  }

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
  return dPhi;
}
