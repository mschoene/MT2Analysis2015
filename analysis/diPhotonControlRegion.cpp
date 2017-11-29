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


bool alsoSignals = true;


int round(float d) {
  return (int)(floor(d + 0.5));
}

void addVariables(MT2Analysis<MT2EstimateTree>* anaTree);

void computeYield( const MT2Sample& sample, const MT2Config& cfg, 
                   MT2Analysis<MT2EstimateTree>* anaTree,
		   int mass_neut2=-1,
		   int mass_neut1=-1 );
void fillYields( MT2Analysis<MT2EstimateZinvGamma>* est, float weight, float ht, int njets, int nbjets, float met, float minMTBmet, float mt2, float iso );
void fillOneTree( MT2EstimateTree* thisTree, const MT2Tree& myTree, float weight, float ht, int njets, int nbjets, float met, float minMTBmet, float mt2);

float DeltaR(float eta1, float eta2, float phi1, float phi2);

float DeltaPhi(float phi1, float phi2);


int main( int argc, char* argv[] ) {



  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|            Running diPhotonControlRegion           |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./diPhotonControlRegion [configFileName] [data/MC/qcd/diPhoton/higgs]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  
  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool onlyData = false;
  bool onlyMC   = false;
  std::string  process;
  if( argc > 2 ) {
    std::string dataMC(argv[2]);
    if( dataMC=="data" ) onlyData = true;
    else if( dataMC=="MC" || dataMC=="mc" || dataMC=="gjets" || dataMC=="qcd" || dataMC=="diPhoton" || dataMC=="higgs" ) {
      onlyMC = true;
      process = dataMC;
    }  
    else {
      std::cout << "-> You passed a second argument that isn't 'data' or 'MC', so I don't know what to do about it." << std::endl;
    }
  }

  if( onlyData ) {
    std::cout << "-> Will run only on data." << std::endl;
  } else if( onlyMC ) {
    std::cout << "-> Will run only on " << process << std::endl;
  } else {
    std::cout << "-> Will run on both data and MC." << std::endl;
  }

  int mass_neutralino2 = -1;
  if( argc > 3 ){
    mass_neutralino2 = std::stof( argv[3] );
    std::cout << "working on the mass point " << mass_neutralino2 << std::endl;
  }


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = cfg.getDiPhotonCRdir();
  system(Form("mkdir -p %s", outputdir.c_str()));


  if( cfg.useMC() && !onlyData ) { // run on MC

    std::string samplesFile = "../samples/samples_" + cfg.mcSamples() + ".dat";
    
    //std::vector<MT2Sample> samples_gammaJet = MT2Sample::loadSamples(samplesFile, "gjet");
    
    std::vector<MT2Sample> samples_gammaJet = MT2Sample::loadSamples(samplesFile, 200, 294);
    
    // std::vector<MT2Sample> samples_gammaJet = MT2Sample::loadSamples(samplesFile, "GJets");
    // if( samples_gammaJet.size()==0 ) {
    //   std::cout << "There must be an error: didn't find any gamma+jet files in " << samplesFile << "!" << std::endl;
    //   exit(1209);
    // }

    std::vector<MT2Sample> samples_diPhoton = MT2Sample::loadSamples(samplesFile, 295, 299);

    std::vector<MT2Sample> samples_higgs = MT2Sample::loadSamples(samplesFile, 900, 999 );

    // std::cout << std::endl << std::endl;
    // std::cout << "-> Loading QCD samples" << std::endl;

    // //std::vector<MT2Sample> samples_qcd = MT2Sample::loadSamples(samplesFile, "qcd");
    std::vector<MT2Sample> samples_qcd = MT2Sample::loadSamples(samplesFile, "QCD");


    MT2Analysis<MT2EstimateTree>* tree = new MT2Analysis<MT2EstimateTree>( "gjets", cfg.regionsSet() );
    //   MT2Analysis<MT2EstimateTree>* tree_diPhoton = new MT2Analysis<MT2EstimateTree>( "diPhotonCRtree", cfg.regionsSet() );
    MT2Analysis<MT2EstimateTree>* tree_higgs = new MT2Analysis<MT2EstimateTree>( "higgs", cfg.regionsSet() );
    MT2Analysis<MT2EstimateTree>* tree_qcd = new MT2Analysis<MT2EstimateTree>( "qcd", cfg.regionsSet() );

    addVariables( tree );
    //   addVariables( tree_diPhoton );
    addVariables( tree_higgs );
    addVariables( tree_qcd );


    if ( process=="mc" || process=="MC" || process=="gjets")
      for( unsigned i=0; i<samples_gammaJet.size(); ++i ) {
    	computeYield( samples_gammaJet[i], cfg, tree );
      }

    // if ( process=="mc" || process=="MC" || process=="diPhoton")
    //   for( unsigned i=0; i<samples_diPhoton.size(); ++i ) {
    // 	computeYield( samples_diPhoton[i], cfg, tree_diPhoton );
    //   }

    if ( process=="mc" || process=="MC" || process=="higgs")
      for( unsigned i=0; i<samples_higgs.size(); ++i ) {
    	computeYield( samples_higgs[i], cfg, tree_higgs );
      }
    if ( process=="mc" || process=="MC" || process=="qcd")
      for( unsigned i=0; i<samples_qcd.size(); ++i ) {
    	computeYield( samples_qcd[i], cfg, tree_qcd );
      }



    std::string mcFile = outputdir + "/mc.root";

    if ( process=="mc" || process=="MC" || process=="gjets"){
      std::string mcFile_gjets = outputdir + "/mc_gjets.root";
      tree         ->writeToFile( mcFile_gjets, "RECREATE"  );
    }
    if ( process=="mc" || process=="MC" || process=="qcd"){
      std::string mcFile_qcd = outputdir + "/mc_qcd.root";
      tree_qcd     ->writeToFile( mcFile_qcd, "RECREATE"  );
    }

    if ( process=="mc" || process=="MC" || process=="higgs"){
      std::string mcFile_higgs = outputdir + "/mc_higgs.root";
      tree_higgs   ->writeToFile( mcFile_higgs, "RECREATE"  );
    }

    // if ( process=="mc" || process=="MC" || process=="diPhoton"){
    //   std::string mcFile_diPhoton = outputdir + "/mc_diPhoton.root";
    //   tree_diPhoton->writeToFile( mcFile_diPhoton, "RECREATE"  );
    // }

    if ( process=="mc" || process=="MC" || process=="diPhoton"){
      if ( process=="mc" || process=="MC" || process=="diPhoton")
    	for( unsigned i=0; i<samples_diPhoton.size(); ++i ) {
    	  MT2Analysis<MT2EstimateTree>* tree_diPhoton = new MT2Analysis<MT2EstimateTree>( "diPhoton", cfg.regionsSet() );
    	  addVariables( tree_diPhoton );
    	  computeYield( samples_diPhoton[i], cfg, tree_diPhoton );
    	  std::string i_str = Form("%d",i);
    	  std::string mcFile_diPhoton = outputdir + "/mc_diPhoton_" + i_str + ".root";
    	  tree_diPhoton->writeToFile( mcFile_diPhoton, "RECREATE"  );
    	}
    }


    std::vector< MT2Analysis< MT2EstimateTree>* > signals;
    if( alsoSignals ) {

      std::cout << std::endl << std::endl;
      std::cout << "-> Loading signal samples" << std::endl;
      //      std::string samplesFile_signals = "samples_signalsSNT.txt";
      std::vector<MT2Sample> samples_sig = MT2Sample::loadSamples(samplesFile, 1000,5000 );

      std::cout << "size " << samples_sig.size() << std::endl;

      if( mass_neutralino2 > 0 ){   
	for( int mass2 =125; mass2<410; mass2+=25){
	  for( int mass1 =0 ; mass1< (mass2 - 120); mass1+=25){

	    int massNeut2 = mass2;
	    if( mass2 == 125 ) massNeut2 = 127;
	    int massNeut1 = mass1;
	    if( mass1 == (mass2-125)) massNeut1 = (mass1-1);
	    if( mass1 == 0 ) massNeut1 = 1;

	    MT2Analysis<MT2EstimateTree>* treeWH = new MT2Analysis<MT2EstimateTree>( Form("SMS_TChiWH_HToGG_%d_%d", massNeut2, massNeut1 ), cfg.regionsSet() );
	    addVariables( treeWH );

	    for( unsigned i=0; i<samples_sig.size(); ++i ) {
	      TString sigSampleName(samples_sig[i].name);
	      if( sigSampleName.Contains("TChiWH")  ){

		computeYield( samples_sig[i], cfg, treeWH, massNeut2, massNeut1 );
    
	      }
	    }
	    signals.push_back( treeWH );
	  }	  

	}//done loop over masses for tchwh
      }

      for( unsigned i=0; i<samples_sig.size(); ++i ) {

	TString sigSampleName(samples_sig[i].name);

	std::string properName = samples_sig[i].sname;
	std::string properNameSave = samples_sig[i].sname;

	if( mass_neutralino2 > 0 || !( sigSampleName.Contains("TChi") ) ){
	  if( sigSampleName.Contains("TChi") && !sigSampleName.Contains("TChiWH")  ){

	    //make a loop from 1 to 1000 over the mass of the neutralino to get the full scan
	    for( int mass =125; mass<220; mass+=25){
	    //	    for( int mass =125; mass<1001; mass+=25){

	      mass_neutralino2 = mass;
	      if(mass == 0  ) mass_neutralino2 = 1;
	      if(mass == 125) mass_neutralino2 = 127;

	      properName = Form("%s_m%d", properNameSave.c_str(), mass_neutralino2 );

	      MT2Analysis<MT2EstimateTree>* thisSig = new MT2Analysis<MT2EstimateTree>( properName, cfg.regionsSet(), samples_sig[i].id );
	      addVariables( thisSig );

	      computeYield( samples_sig[i], cfg, thisSig, mass_neutralino2 );

	      signals.push_back( thisSig );

	      if( sigSampleName.Contains("TChiHH")){
		  properName = Form("%s_m%d_HH0p25", properNameSave.c_str(), mass_neutralino2 );
		  MT2Analysis<MT2EstimateTree>* thisSig_HH0p25 = new MT2Analysis<MT2EstimateTree>( properName, cfg.regionsSet(), samples_sig[i].id );
		  addVariables( thisSig_HH0p25 );
		  computeYield( samples_sig[i], cfg, thisSig_HH0p25, mass_neutralino2 );
		  signals.push_back( thisSig_HH0p25 ); 
		}


	    }



	  }else if( !sigSampleName.Contains("TChiWH")  ){
	    MT2Analysis<MT2EstimateTree>* thisSig = new MT2Analysis<MT2EstimateTree>( samples_sig[i].sname, cfg.regionsSet(), samples_sig[i].id );
	    addVariables( thisSig );
	    computeYield( samples_sig[i], cfg, thisSig);
	    signals.push_back( thisSig );
	  }

	}
	
      }  // for signals

    }

    std::cout << "signal size is " << signals.size() << std::endl;

    std::string mcFile_sig = outputdir + "/mc_sig.root";
    for( unsigned i=0; i<signals.size(); ++i )
      //   if( i==0)
      //	signals[i]->writeToFile( mcFile_sig, "RECREATE");
      //      else
      signals[i]->addToFile( mcFile_sig);

    // //comment here again
    // if( cfg.dummyAnalysis() ) {

    //   // emulate data:
    //   tree->writeToFile( outputdir + "/data.root" );

    // } // if dummy

  } // if run on MC



  if( !onlyMC && !cfg.dummyAnalysis() ) {

    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";
 
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "DoubleEG");
    //  std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "merged");
    //    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "SinglePhoton");

    if( samples_data.size()==0 ) {

      std::cout << std::endl;
      std::cout << "-> WARNING!! Didn't find any data in file: " << samplesFile_data << "!" << std::endl;
      std::cout << "-> Exiting." << std::endl;
      std::cout << std::endl;

    } else {


      MT2Analysis<MT2EstimateTree>* tree = new MT2Analysis<MT2EstimateTree>( "diPhoton_data", cfg.regionsSet() );

      addVariables( tree );
      
      for( unsigned i=0; i<samples_data.size(); ++i ) {
        computeYield( samples_data[i], cfg, tree );
      }

      tree        ->writeToFile( outputdir + "/data.root", "RECREATE" );
            
    } // if data samples > 0

  } // if run on data


  return 0;

}








void computeYield( const MT2Sample& sample, const MT2Config& cfg,
                   MT2Analysis<MT2EstimateTree>* anaTree, 
		   int mass_neut2, int mass_neut1 ) {


  float isoCut = cfg.gammaIsoCut();  

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;

  bool isQCD  =   (sample.id>=100 && sample.id<200) ||  (sample.id>=3100 && sample.id<3200 );
  bool isGJets = sample.id>=200 && sample.id<295;
  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);


  TString sigSampleName(sample.name);
  bool isEW = 0;
  if(sigSampleName.Contains("TChi"))
    isEW = 1;

  bool isEW_WH = 0;
  if(sigSampleName.Contains("TChiWH"))
    isEW_WH = 1;


  int nentries = tree->GetEntries();

  //  for( int iEntry=0; iEntry< 5000 ; ++iEntry ) {
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( isEW && mass_neut2>0 && ( mass_neut2 != myTree.GenSusyMNeutralino2) ) {
      if ( !isEW )
	std::cout << "Why are you doing this" << std::endl;

      // if( mass_neut2 != myTree.GenSusyMNeutralino2) 
      //	std::cout << "not the mass i want" << std::endl;

      continue;
    }


    if( isEW_WH && mass_neut2>0 && ( mass_neut1 != myTree.GenSusyMNeutralino) ) {
      if ( !isEW_WH )
	std::cout << "Why are you doing this" << std::endl;
      // if( mass_neut2 != myTree.GenSusyMNeutralino2) 
      //	std::cout << "not the mass i want" << std::endl;
      continue;
    }



    //crazy events! To be piped into a separate txt file
    if(myTree.jet_pt[0] > 13000){
      std::cout << "Rejecting weird event at run:lumi:evt = " << myTree.run << ":" << myTree.lumi << ":" << myTree.evt << std::endl;
      continue;
    }

    //NEW check if is there is a nan
    if( isnan(myTree.ht) || isnan(myTree.met_pt) ||  isinf(myTree.ht) || isinf(myTree.met_pt)  ){
      std::cout << "Rejecting nan/inf event at run:lumi:evt = " << myTree.run << ":" << myTree.lumi << ":" << myTree.evt << std::endl;
      continue;
    }

    //    if( myTree.met_miniaodPt/myTree.met_caloPt > 5. ) continue;
     
    if( myTree.isData ) {
      if( !( myTree.isGolden ) ) continue;
      if( !(myTree.passFilters() ) ) continue;
      if( !(myTree.HLT_DiPhoton30) ) continue;
    }

    if( myTree.ngamma<2 ) continue;



    // if( myTree.gamma_pt[0]<30 ) continue;
    // if( myTree.gamma_pt[1]<20 ) continue;


    // TLorentzVector gamma;
    // gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );

    // //Need the lorentz vectors of the photons first
    TLorentzVector *LVec = new TLorentzVector[myTree.ngamma];

    std::vector<int> indexG;     

    for(int i=0; i< myTree.ngamma; i++){
      if(myTree.gamma_chHadIso[i]<=2.5){

	bool isNotALep = true;
    	if( myTree.nlep>0){
    	  for(int lep=0; lep< myTree.nlep; lep++){
    	    float dR = sqrt( ((myTree.lep_eta[lep]-myTree.gamma_eta[i])*(myTree.lep_eta[lep]-myTree.gamma_eta[i])) + ((myTree.lep_phi[lep]-myTree.gamma_phi[i])*(myTree.lep_phi[lep]-myTree.gamma_phi[i])) );
    	    //Clean out photons from electrons with dR<1, from muons wiht dR<0.5
    	    if( !( (abs(myTree.lep_pdgId[lep])==11 && dR>1.0 ) || (abs(myTree.lep_pdgId[lep])==13 && dR>0.5) ) ){ 
	      isNotALep = false;
	    }//only take the nice photons
    	  }//done loop over leptons
    	}

	if( myTree.nlep==0 || isNotALep==true ) {
    	  LVec[i].SetPtEtaPhiM(myTree.gamma_pt[i], myTree.gamma_eta[i],myTree.gamma_phi[i], myTree.gamma_mass[i]);
	  indexG.push_back( i );
 
    	}
      }
    }

    if( indexG.size()<2 ){
      //    if( ngam<2 || gotSecond==false ){
      
    //   std::cout << "Kicking out event with " << myTree.ngamma << " of with only " << ngam << " are clean and got a second = " << << std::endl;     
      continue;
    }

    if( indexG.size() > myTree.ngamma ){
      std::cout << "FUCK how do you find " << indexG.size() << " when there are only " << myTree.ngamma << " in total " << std::endl;
      continue;
    }

    if( (isQCD || isGJets) && myTree.gamma_mcMatchId[indexG[0]]==22 && myTree.gamma_mcMatchId[indexG[1]]==22 ) continue; //remove promptprom from QCD & GJets as it is already in the BOX

    TLorentzVector Hvec;
    if( indexG.size()>=2 )
      Hvec = LVec[indexG[0]] + LVec[indexG[1]]; //gamma invariant mass
    else 
      continue; //should already be taken care of in the skim


    float minMTBmet = myTree.minMTBMet;
    float met       = myTree.gg_met_pt;
    int njets       = myTree.gg_nJet30;
    int nbjets      = myTree.gg_nBJet20;    
    float ht        = myTree.gg_ht;
    float mt2       = (njets>1) ? myTree.mt2 : ht;




    TLorentzVector gamma0;
    gamma0.SetPtEtaPhiM( myTree.gamma_pt[indexG[0]], myTree.gamma_eta[indexG[0]], myTree.gamma_phi[indexG[0]], myTree.gamma_mass[indexG[0]] );

    TLorentzVector gamma1;
    gamma1.SetPtEtaPhiM( myTree.gamma_pt[indexG[1]], myTree.gamma_eta[indexG[1]], myTree.gamma_phi[indexG[1]], myTree.gamma_mass[indexG[1]] );

    TLorentzVector pseudoJet1;
    pseudoJet1.SetPtEtaPhiM( myTree.pseudoJet1_pt, myTree.pseudoJet1_eta, myTree.pseudoJet1_phi, myTree.pseudoJet1_mass );

    // TLorentzVector pseudoJet2;
    // pseudoJet2.SetPtEtaPhiM( myTree.pseudoJet2_pt, myTree.pseudoJet2_eta, myTree.pseudoJet2_phi, myTree.pseudoJet2_mass );


    float scalarProd11 = gamma0.Vect()*pseudoJet1.Vect();
    float scalarProd12 = gamma1.Vect()*pseudoJet1.Vect();



    // float minMTBmet = myTree.minMTBMet;
    // float met       = myTree.met_pt;
    // int njets       = myTree.nJet30;
    // int nbjets      = myTree.nBJet20;    
    // float ht        = myTree.ht;
    // float mt2       = (njets>1) ? myTree.mt2 : ht;


    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi(); 

    //std::string anaTreeName = anaTree->getName();
    TString anaTreeName( anaTree->getName() );

    if(sigSampleName.Contains("HZ"))
      weight *= 0.5;
    else if( anaTreeName.Contains("HH0p25"))
      weight *= 0.25;

  

    // if(sigSampleName.Contains("T2bH")){ // needed for first production of T2bH, fixed for oct15 version
    //   weight = weight / (0.00227*0.5824 ) * ( 0.00227*2. - 0.00227*0.00227) ;
    // }

    // std::string whatPt;
    // if( myTree.h_pt/myTree.h_mass > 0.8 )
    //   whatPt = "hiPt";
    // else
    //   whatPt = "loPT";

    // MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, minMTBmet, myTree.h_mass, myTree.h_pt/myTree.h_mass );
    
    MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, Hvec.Perp()/Hvec.M() , Hvec.M());
    //MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, myTree.h_pt/myTree.h_mass , myTree.h_mass );
    //    MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisTree==0 ) continue;



    // if(  myTree.gg_nJet30 == 0 )
    //   std::cout << "you do actually have one of those in here ffs" << std::endl;


    float h_mass = -999;
    float h_pt = -999;

    // if( myTree.ngamma>=2 ){
    //     h_mass = Hve.Mc();
    //     h_pt   = Hvec.Perp();
    // }

    // if( fabs(myTree.h_mass - Hvec.M() ) > 0.1 )
    //   std::cout << myTree.h_mass << " vs " << Hvec.M() << std::endl;


    if( (fabs(myTree.gamma_eta[indexG[0]])<=1.4) && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.010 ) continue;
    if( (fabs(myTree.gamma_eta[indexG[1]])<=1.4) && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.010 ) continue;
    if( (fabs(myTree.gamma_eta[indexG[0]])>1.4)  && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.0267 ) continue;
    if( (fabs(myTree.gamma_eta[indexG[1]])>1.4)  && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.0267 ) continue;


    h_mass = Hvec.M();
    h_pt = Hvec.Perp();

    // h_mass = myTree.h_mass;
    // h_pt = myTree.h_pt;


    float massEG00 = -99.;
    float massEG01 = -99.;
    float massEG10 = -99.;
    float massEG11 = -99.;

    float diLepMass = -99.;
    bool isDiLepZ = 0;
    //di lep
    if( ( myTree.nlep>1 )){

      if( ((myTree.lep_pdgId[0]*myTree.lep_pdgId[1])< 0 )  && myTree.lep_pt[0]> 20 && myTree.lep_pt[1]>20  && (abs(myTree.lep_pdgId[0])==13 || ( abs(myTree.lep_pdgId[0])==11 && myTree.lep_tightId[0]> 0.5 && abs(myTree.lep_pdgId[1])==11 && myTree.lep_tightId[1]> 0.5  )) ) {
	//why use 30 if you can use 20 in pt      if( ((myTree.lep_pdgId[0]*myTree.lep_pdgId[1])< 0 )  && myTree.lep_pt[0]> 30 && myTree.lep_pt[1]>30  && (abs(myTree.lep_pdgId[0])==13 || ( abs(myTree.lep_pdgId[0])==11 && myTree.lep_tightId[0]> 0.5 && abs(myTree.lep_pdgId[1])==11 && myTree.lep_tightId[1]> 0.5  )) ) {

	//Need the lorentz vectors of the leptons first
	TLorentzVector *LVec = new TLorentzVector[2];
	for(int i=0; i< 2; i++){
	  LVec[i].SetPtEtaPhiM(myTree.lep_pt[i], myTree.lep_eta[i],myTree.lep_phi[i], myTree.lep_mass[i]);
	}

	TLorentzVector Zvec = LVec[0] + LVec[1]; //leptons invariant mass
	diLepMass = Zvec.M(); // compare with myTree.zll_mass

    	if( fabs( diLepMass -91.19) <= 20 ) 
	  isDiLepZ  = 1;

	//build M egamma for purity, ie get those dirty photon electrons out
	TLorentzVector vecEG00 = LVec[0] + gamma0 ; //leptons invariant mass
	TLorentzVector vecEG01 = LVec[0] + gamma1 ; //leptons invariant mass
	TLorentzVector vecEG10 = LVec[1] + gamma0 ; //leptons invariant mass
	TLorentzVector vecEG11 = LVec[1] + gamma1 ; //leptons invariant mass

	massEG00 = vecEG00.M();
	massEG01 = vecEG01.M();
	massEG10 = vecEG10.M();
	massEG11 = vecEG11.M();

      }
    }


    bool is1Mu = 0;
    bool is1El = 0;

    //SINGLE LEPTON region
    if( myTree.nlep==1 ){

      if(  myTree.lep_pt[0]> 20  && (abs(myTree.lep_pdgId[0])==13 || ( abs(myTree.lep_pdgId[0])==11 && myTree.lep_tightId[0]> 0.5  )) ) {
	if( (abs(myTree.lep_pdgId[0])==11 ))
	  is1El = 1;
	else if((abs(myTree.lep_pdgId[0])==13 ))
	  is1Mu = 1;
     
      }
    }
  

    float diBMass = -99.;     bool isDiBZ = 0;bool isDiBH = 0;    
    std::vector<float> bMassVector;
    std::vector<float> bPTVector;
    float diBpT = -99;

    if( myTree.gg_nBJet20 > 1 ){
      std::vector<int> indexB;       //   int indexB0 = -1;      int indexB1 = -1;
      for( int i=0; i< myTree.njet; i++){
	if ( fabs(myTree.jet_eta[i]) > 2.4 ) 
	  continue;
	if( myTree.jet_btagCSV[i]>0.8484 ){
	  indexB.push_back( i );
	}
      }

      if( indexB.size()>1 ){
	for(std::vector<int>::iterator bIndex = indexB.begin(); bIndex != indexB.end()-1; ++bIndex){
	  for(std::vector<int>::iterator bIndex2 = indexB.begin()+1; bIndex2 != indexB.end(); ++bIndex2){
	    //  std::cout << "Working on b jets of jets at indeces " << *bIndex << ", " << *bIndex2 << std::endl;
	    TLorentzVector *LVec = new TLorentzVector[2];
	    LVec[0].SetPtEtaPhiM(myTree.jet_pt[*bIndex], myTree.jet_eta[*bIndex], myTree.jet_phi[*bIndex], myTree.jet_mass[*bIndex]);
	    LVec[1].SetPtEtaPhiM(myTree.jet_pt[*bIndex2], myTree.jet_eta[*bIndex2], myTree.jet_phi[*bIndex2], myTree.jet_mass[*bIndex2]);
	    TLorentzVector Bvec = LVec[0] + LVec[1]; //bb invariant mass
	    diBMass = Bvec.M();
	    bMassVector.push_back( diBMass );
	    bPTVector.push_back( Bvec.Perp() );
	    //	    std::cout << "FOUNG A BB pair with mass "  << diBMass << std::endl;
	    diBpT  = Bvec.Perp();
	    //	if( (fabs( diBMass-91.19) <= 20.) ||  (fabs(diBMass-125) <= 20.)  ) 
	    //isDiB = 1;
	  }
	}

	diBMass  = bMassVector[0]; //set to first one then see if there is a better one
	//	diBMass0 = bMassVector[0]; // first in list
	diBpT  = bPTVector[0]; //set to first one then see if there is a better one
	//	diBpT0 = bPTVector[0]; // first in list
	//	std::cout << bMassVector.size()  << std::endl;
	//	std::cout << "Original "  << diBMass << std::endl;

	for( unsigned int bbCand=0; bbCand<bMassVector.size(); bbCand++){
	  if( (bMassVector[ bbCand ] <140) && (bMassVector[ bbCand ] >= 95) )
	    isDiBH = 1;
	}
	for( unsigned int bbCand=0; bbCand<bMassVector.size(); bbCand++){
	  if( (bMassVector[ bbCand ] <95) && (bMassVector[ bbCand ] >= 60)  && (isDiBH==0) )
	    isDiBZ = 1;
	}
	//   //	  std::cout << "current Original "  << diBMass << std::endl;
	//   for( unsigned int bbCand2=bbCand+1; bbCand2<bMassVector.size(); bbCand2++){
	//     if( ((fabs( bMassVector[bbCand2]-125.)) < (fabs( diBMass-125.) )) || ((fabs( bMassVector[bbCand2]- 90.)) < (fabs( diBMass - 90.)) ) )
	//       //  std::cout << "Diff up = " <<  fabs( bMassVector[bbCand2]- 125.) << " vs " << fabs( diBMass- 125.) << std::endl;
	//       //std::cout << "Diff up = " <<  fabs( bMassVector[bbCand2]- 90.) << " vs " << fabs( diBMass  - 90.) << std::endl;
	//       diBMass = bMassVector[bbCand2];
	//       diBpT   = bPTVector[bbCand2];
	//       //  std::cout << "changed to "  << diBMass << std::endl;
	//   }  
	// }

      } // only do all this if there are 2 bjets

    }

  

    thisTree->assignVar("massEG00", massEG00 );
    thisTree->assignVar("massEG01", massEG01 );
    thisTree->assignVar("massEG10", massEG10 );
    thisTree->assignVar("massEG11", massEG11 );


    float lep_id = -99;
    float lep_tightId0 = -99;
    float lep_tightId1 = -99;
    float lep_eta0 = -99;
    float lep_eta1 = -99;
    float lep_phi0 = -99;
    float lep_phi1 = -99;

    if( myTree.nlep >1 ){
      lep_id = myTree.lep_pdgId[0];
      lep_tightId0 = myTree.lep_tightId[0];
      lep_tightId1 = myTree.lep_tightId[1];
      lep_eta0 =  myTree.lep_eta[0];
      lep_eta1 =  myTree.lep_eta[1];
      lep_phi0 = myTree.lep_phi[0];
      lep_phi1 = myTree.lep_phi[1];
    }

    thisTree->assignVar("lep_pdgId0", lep_id );
    thisTree->assignVar( "lep_tightId0", lep_tightId0 );
    thisTree->assignVar( "lep_tightId1", lep_tightId1 );
    thisTree->assignVar( "lep_eta0", lep_eta0 );
    thisTree->assignVar( "lep_eta1", lep_eta1 );
    thisTree->assignVar( "lep_phi0", lep_phi0 );
    thisTree->assignVar( "lep_phi1", lep_phi1 );


    thisTree->assignVar( "diLepMass", diLepMass );
    thisTree->assignVar( "isDiLepZ", isDiLepZ );

    thisTree->assignVar( "is1El", is1El);
    thisTree->assignVar( "is1Mu", is1Mu);

    thisTree->assignVar( "diBMass", diBMass );
    thisTree->assignVar( "isDiBH", isDiBH);
    thisTree->assignVar( "isDiBZ", isDiBZ );
    thisTree->assignVar( "diBpT", diBpT );

    thisTree->assignVar( "scProd11", scalarProd11 );
    thisTree->assignVar( "scProd12", scalarProd12 );


    thisTree->assignVar( "angleGammas", LVec[0].Angle(LVec[1].Vect()) );

    thisTree->assignVar( "met_phi", myTree.met_phi );

    thisTree->assignVar( "h_phi", Hvec.Phi() );
    thisTree->assignVar( "h_eta", Hvec.Eta() );
    thisTree->assignVar( "h_mass", h_mass );
    thisTree->assignVar( "h_pt",   h_pt   );

    thisTree->assignVar( "ptGamma0", myTree.gamma_pt[indexG[0]] );
    thisTree->assignVar( "etaGamma0", myTree.gamma_eta[indexG[0]] );
    thisTree->assignVar( "phiGamma0", myTree.gamma_phi[indexG[0]] );
    thisTree->assignVar( "r9Gamma0", myTree.gamma_r9[indexG[0]] );
    thisTree->assignVar( "sigmaIetaIetaGamma0", myTree.gamma_sigmaIetaIeta[indexG[0]] );

    // thisTree->assignVar( "chHadIsoGamma0", myTree.gamma_chHadIso[indexG[0]] );
    // thisTree->assignVar( "chHadIsoRCGamma0", myTree.gamma_chHadIsoRC[indexG[0]] );
    // thisTree->assignVar( "chHadIsoRC04Gamma0", myTree.gamma_chHadIsoRC04[indexG[0]] );

    thisTree->assignVar( "ptGamma1", myTree.gamma_pt[indexG[1]] );
    thisTree->assignVar( "etaGamma1", myTree.gamma_eta[indexG[1]] );
    thisTree->assignVar( "phiGamma1", myTree.gamma_phi[indexG[1]] );
    thisTree->assignVar( "r9Gamma1", myTree.gamma_r9[indexG[1]] );
    thisTree->assignVar( "sigmaIetaIetaGamma1", myTree.gamma_sigmaIetaIeta[indexG[1]] );

    // thisTree->assignVar( "chHadIsoGamma1", myTree.gamma_chHadIso[indexG[1]] );
    // thisTree->assignVar( "chHadIsoRCGamma1", myTree.gamma_chHadIsoRC[indexG[1]] );
    // thisTree->assignVar( "chHadIsoRC04Gamma1", myTree.gamma_chHadIsoRC04[indexG[1]] );

    thisTree->assignVar( "gamma_jet1_pt", myTree.gamma_jet1_pt );
    thisTree->assignVar( "gamma_jet2_pt", myTree.gamma_jet2_pt );
    // thisTree->assignVar( "gamma_chHadIsoRC",  myTree.gamma_chHadIsoRC[0] );
    thisTree->assignVar( "drMinParton", myTree.gamma_drMinParton[0] );
    //thisTree->assignVars( ht, njets, nbjets, met, mt2 ); //manually override these (so change nbjets)
    thisTree->assignVar( "raw_mt2", myTree.mt2 );


    thisTree->assignVar( "gg_nJets", myTree.gg_nJet30 );
    thisTree->assignVar( "gg_nBJets", myTree.gg_nBJet20 );
    thisTree->assignVar( "gg_mt2", myTree.gg_mt2 );
    thisTree->assignVar( "gg_ht", myTree.gg_ht );
    thisTree->assignVar( "gg_deltaPhiMin", myTree.gg_deltaPhiMin );
    thisTree->assignVar( "gg_diffMetMht", myTree.gg_diffMetMht );
    thisTree->assignVar( "gg_mht", myTree.gg_mht_pt );
    thisTree->assignVar( "gg_met", myTree.gg_met_pt );
    thisTree->assignVar( "gg_met_phi", myTree.gg_met_phi );
    thisTree->assignVar( "gg_jet1_pt", myTree.gg_jet1_pt );
    thisTree->assignVar( "gg_jet2_pt", myTree.gg_jet2_pt );
    thisTree->assignVar( "jet1_pt", myTree.jet1_pt );
    thisTree->assignVar( "jet2_pt", myTree.jet2_pt );

    if(sigSampleName.Contains("T2bH")){ 
      thisTree->assignVar( "weight_isr", myTree.weight_isr/myTree.weight_isr_av );
      thisTree->assignVar( "weight_isr_UP", myTree.weight_isr_UP/myTree.weight_isr_UP_av );
      thisTree->assignVar( "weight_isr_DN", myTree.weight_isr_DN/myTree.weight_isr_DN_av );

      weight *=  myTree.weight_isr/myTree.weight_isr_av;
    }

    thisTree->assignTree(myTree, weight ); // assign all to defaults
 
    thisTree->yield->Fill( h_mass , weight );
    //    thisTree->yield->Fill(mt2, weight );
    thisTree->tree->Fill();
    

    // fillOneTree( thisTree, myTree, weight, ht, njets, nbjets, met, minMTBmet, mt2 );
    
  } // for entries


  anaTree->finalize();
    

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



void fillOneTree( MT2EstimateTree* thisTree, const MT2Tree& myTree, float weight, float ht, int njets, int nbjets, float met, float minMTBmet, float mt2 ) {

  thisTree->yield->Fill(mt2, weight );
  //  thisTree->assignVar( "iso", iso );
  thisTree->assignVar( "sietaieta", myTree.gamma_sigmaIetaIeta[0] );
  thisTree->assignVar( "ptGamma", myTree.gamma_pt[0] );
  thisTree->assignVar( "etaGamma", myTree.gamma_eta[0] );
  thisTree->assignVar( "jet1_pt", myTree.gamma_jet1_pt );
  thisTree->assignVar( "jet2_pt", myTree.gamma_jet2_pt );
  thisTree->assignVar( "gamma_chHadIsoRC",  myTree.gamma_chHadIsoRC[0] );
  thisTree->assignVar( "drMinParton", myTree.gamma_drMinParton[0] );


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




void addVariables(MT2Analysis<MT2EstimateTree>* tree){

  MT2EstimateTree::addVar( tree, "gamma_jet1_pt" );
  MT2EstimateTree::addVar( tree, "gamma_jet2_pt" );

  MT2EstimateTree::addVar( tree, "drMinParton" );
  MT2EstimateTree::addVar( tree, "raw_mt2" );

  MT2EstimateTree::addVar( tree, "met_phi" );

  MT2EstimateTree::addVar( tree, "scProd12" );
  MT2EstimateTree::addVar( tree, "scProd11" );

  MT2EstimateTree::addVar( tree, "angleGammas" );
  MT2EstimateTree::addVar( tree, "h_phi" );
  MT2EstimateTree::addVar( tree, "h_eta" );
  MT2EstimateTree::addVar( tree, "h_mass");
  MT2EstimateTree::addVar( tree, "h_pt" );

  MT2EstimateTree::addVar( tree, "ptGamma0");
  MT2EstimateTree::addVar( tree, "etaGamma0");
  MT2EstimateTree::addVar( tree, "phiGamma0");
  MT2EstimateTree::addVar( tree, "r9Gamma0");
  MT2EstimateTree::addVar( tree, "sigmaIetaIetaGamma0");

  // MT2EstimateTree::addVar( tree, "chHadIsoGamma0");
  // MT2EstimateTree::addVar( tree, "chHadIsoRCGamma0");
  // MT2EstimateTree::addVar( tree, "chHadIsoRC04Gamma0");

  MT2EstimateTree::addVar( tree, "ptGamma1");
  MT2EstimateTree::addVar( tree, "etaGamma1");
  MT2EstimateTree::addVar( tree, "phiGamma1");
  MT2EstimateTree::addVar( tree, "r9Gamma1");
  MT2EstimateTree::addVar( tree, "sigmaIetaIetaGamma1");

  // MT2EstimateTree::addVar( tree, "chHadIsoGamma1");
  // MT2EstimateTree::addVar( tree, "chHadIsoRCGamma1");
  // MT2EstimateTree::addVar( tree, "chHadIsoRC04Gamma1");

  MT2EstimateTree::addVar( tree, "gg_nJets" );
  MT2EstimateTree::addVar( tree, "gg_nBJets" );
  MT2EstimateTree::addVar( tree, "gg_ht" );
  MT2EstimateTree::addVar( tree, "gg_mt2" );
  MT2EstimateTree::addVar( tree, "gg_deltaPhiMin" );
  MT2EstimateTree::addVar( tree, "gg_diffMetMht" );
  MT2EstimateTree::addVar( tree, "gg_mht" );
  MT2EstimateTree::addVar( tree, "gg_met" );
  MT2EstimateTree::addVar( tree, "gg_met_phi" );
  MT2EstimateTree::addVar( tree, "gg_jet1_pt" );
  MT2EstimateTree::addVar( tree, "gg_jet2_pt" );

  MT2EstimateTree::addVar( tree, "jet1_pt" );
  MT2EstimateTree::addVar( tree, "jet2_pt" );

  MT2EstimateTree::addVar( tree, "weight_isr" );
  MT2EstimateTree::addVar( tree, "weight_isr_UP" );
  MT2EstimateTree::addVar( tree, "weight_isr_DN" );

  MT2EstimateTree::addVar( tree, "isDiBH" );
  MT2EstimateTree::addVar( tree, "isDiBZ" );
  MT2EstimateTree::addVar( tree, "diBMass" );
  MT2EstimateTree::addVar( tree, "isDiLepZ" );
  MT2EstimateTree::addVar( tree, "diLepMass" ); 

  MT2EstimateTree::addVar( tree, "is1El" );
  MT2EstimateTree::addVar( tree, "is1Mu" );

  MT2EstimateTree::addVar( tree, "diBpT" );

  MT2EstimateTree::addVar( tree, "lep_pdgId0");
  // MT2EstimateTree::addVar( tree, "lep_pdgId1");

  MT2EstimateTree::addVar( tree, "lep_tightId0" );
  MT2EstimateTree::addVar( tree, "lep_tightId1" );

  MT2EstimateTree::addVar( tree, "lep_eta0" );
  MT2EstimateTree::addVar( tree, "lep_eta1" );

  MT2EstimateTree::addVar( tree, "lep_phi0" );
  MT2EstimateTree::addVar( tree, "lep_phi1" );

  MT2EstimateTree::addVar( tree, "massEG00" );
  MT2EstimateTree::addVar( tree, "massEG01" );
  MT2EstimateTree::addVar( tree, "massEG10" );
  MT2EstimateTree::addVar( tree, "massEG11" );

}
