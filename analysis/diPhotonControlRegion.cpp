#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2EstimateTree.h"

#include "interface/Davismt2.h"
#include "interface/Hemisphere.h"
#define mt2_cxx
#include "../interface/mt2.h"

#include "TLorentzVector.h"
#include "TH1F.h"
#include "TF1.h"

#include <memory>
using std::unique_ptr;
using std::string;

bool alsoSignals = true;
bool is2017=false; //true for 2017 data&mc, otherwise it's 2016 data&mc (for trigger and filters, evlt btag id)


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

// void getHemispheres( std::vector<float> input_pt, std::vector<float> input_eta, std::vector<float> input_phi, std::vector<float> input_mass, TLorentzVector *v1, TLorentzVector *v2);
// Float_t calcMT2(std::vector<float> input_pt, std::vector<float> input_eta, std::vector<float> input_phi, std::vector<float> input_mass, float met, float metphi);

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
    std::cout << "USAGE: ./diPhotonControlRegion [configFileName] [data/MC/qcd/diPhoton/higgs/diH]" << std::endl;
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
    else if( dataMC=="MC" || dataMC=="mc" || dataMC=="gjets" || dataMC=="qcd" || dataMC=="diPhoton" || dataMC=="diH" || dataMC=="higgs" || dataMC=="HZ" || dataMC=="HH" || dataMC=="WH" || dataMC=="T2bH" ) {
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


  TString name2017( argv[1] );
  if( name2017.Contains("2017")  ){
    std::cout << "changing to 2017 settings " << std::endl;
    is2017 = true;
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

    std::vector<MT2Sample> samples_higgs = MT2Sample::loadSamples(samplesFile, 900, 969 );
    std::vector<MT2Sample> samples_diH = MT2Sample::loadSamples(samplesFile, 970, 980 );

    // std::cout << std::endl << std::endl;
    // std::cout << "-> Loading QCD samples" << std::endl;

    // //std::vector<MT2Sample> samples_qcd = MT2Sample::loadSamples(samplesFile, "qcd");
    std::vector<MT2Sample> samples_qcd = MT2Sample::loadSamples(samplesFile, "QCD");


    MT2Analysis<MT2EstimateTree>* tree = new MT2Analysis<MT2EstimateTree>( "gjets", cfg.regionsSet() );
    //   MT2Analysis<MT2EstimateTree>* tree_diPhoton = new MT2Analysis<MT2EstimateTree>( "diPhotonCRtree", cfg.regionsSet() );
    MT2Analysis<MT2EstimateTree>* tree_higgs = new MT2Analysis<MT2EstimateTree>( "higgs", cfg.regionsSet() );
    MT2Analysis<MT2EstimateTree>* tree_diH = new MT2Analysis<MT2EstimateTree>( "diH", cfg.regionsSet() );
    MT2Analysis<MT2EstimateTree>* tree_qcd = new MT2Analysis<MT2EstimateTree>( "qcd", cfg.regionsSet() );

    addVariables( tree );
    //   addVariables( tree_diPhoton );
    addVariables( tree_higgs );
    addVariables( tree_diH );
    addVariables( tree_qcd );


    if ( process=="mc" || process=="MC" || process=="gjets")
      for( unsigned i=0; i<samples_gammaJet.size(); ++i ) {
    	computeYield( samples_gammaJet[i], cfg, tree );
      }

    // if ( process=="mc" || process=="MC" || process=="diPhoton")
    //   for( unsigned i=0; i<samples_diPhoton.size(); ++i ) {
    // 	computeYield( samples_diPhoton[i], cfg, tree_diPhoton );
    //   }

    if ( process=="mc" || process=="MC" || process=="diH")
      for( unsigned i=0; i<samples_diH.size(); ++i ) {
    	computeYield( samples_diH[i], cfg, tree_diH );
      }
    // if ( process=="mc" || process=="MC" || process=="higgs")
    //   for( unsigned i=0; i<samples_higgs.size(); ++i ) {
    // 	computeYield( samples_higgs[i], cfg, tree_higgs );
    //   }
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
    if ( process=="mc" || process=="MC" || process=="diH"){
      std::string mcFile_diH = outputdir + "/mc_diH.root";
      tree_diH     ->writeToFile( mcFile_diH, "RECREATE"  );
    }

    // if ( process=="mc" || process=="MC" || process=="higgs"){
    //   std::string mcFile_higgs = outputdir + "/mc_higgs.root";
    //   tree_higgs   ->writeToFile( mcFile_higgs, "RECREATE"  );
    // }

    if ( process=="higgs"){
      for( unsigned i=0; i<samples_higgs.size(); ++i ) {
	MT2Analysis<MT2EstimateTree>* tree_higgs = new MT2Analysis<MT2EstimateTree>( "higgs", cfg.regionsSet() );
	addVariables( tree_higgs );
	computeYield( samples_higgs[i], cfg, tree_higgs );
	std::string i_str = Form("%d",i);
	std::string mcFile_higgs = outputdir + "/mc_higgs_" + i_str + ".root";
	tree_higgs->writeToFile( mcFile_higgs, "RECREATE"  );
      }
    }

   
    if ( process=="mc" || process=="MC" || process=="diPhoton"){
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

      if( process== "WH"  )
	if( mass_neutralino2 > 0 ){   
	  for( int mass2 =125; mass2<701; mass2+=25){
	  //	  for( int mass2 =125; mass2<210; mass2+=25){
	    for( int mass1 =0 ; mass1< (mass2 - 120); mass1+=25){

	      if( mass1 > 302 ) continue; //scan stops there, cut of triangular shape (stop was at 300

	      int massNeut2 = mass2;
	      if( mass2 == 125 ) massNeut2 = 127;
	      int massNeut1 = mass1;
	      if( mass1 == (mass2-125)) massNeut1 = (mass1-1);
	      if( mass1 == 0 ) massNeut1 = 1;

	      std::cout << "Working on point " << mass1 << " : " << mass2 << std::endl;
	      std::cout << "Working on mass point " << massNeut1 << " : " << massNeut2 << std::endl;

	      std::string  bound = "";
	      if(massNeut2 < 201)
		bound = "200";
	      else if(massNeut2 < 301)
		bound = "300";
	      else if(massNeut2 < 401)
		bound = "400";
	      else if(massNeut2 < 501)
		bound = "500";
	      else if(massNeut2 < 601)
		bound = "600";
	      else
		bound = "INF";

	      std::cout << bound << std::endl;

	      MT2Analysis<MT2EstimateTree>* treeWH = new MT2Analysis<MT2EstimateTree>( Form("SMS_TChiWH_HToGG_m%d_m%d", massNeut2, massNeut1 ), cfg.regionsSet() );
	      addVariables( treeWH );

	      for( unsigned i=0; i<samples_sig.size(); ++i ) {

		TString sigSampleFile(samples_sig[i].file);
		//		if( ! sigSampleFile.Contains(Form("WH%s", bound.c_str() ))  )
		if( (! sigSampleFile.Contains(Form("WH%s", bound.c_str() )) ) &&  !(bound=="200" && sigSampleFile.Contains("TChiWH_HToGG_175")  )     )
		  continue;		

		TString sigSampleName(samples_sig[i].name);
		std::cout << "sample name = " << sigSampleName  << std::endl;
		if( sigSampleName.Contains("TChiWH")  ){

		  computeYield( samples_sig[i], cfg, treeWH, massNeut2, massNeut1 );
    
		}
	      }
	      signals.push_back( treeWH );
	    }	  

	  }//done loop over masses for tchwh
	}


      if( process== "HH"  ){ 
	//make a loop from 1 to 1000 over the mass of the neutralino to get the full scan
	//	for( int mass =125; mass<1001; mass+=25){
	for( int mass =125; mass<531; mass+=25){
	  //	    for( int mass =125; mass<1001; mass+=25){
	  mass_neutralino2 = mass;
	  if(mass == 0  ) mass_neutralino2 = 1;
	  if(mass == 125) mass_neutralino2 = 127;

	  std::string properName = Form("SMS_TChiHH_HToGG_m%d", mass_neutralino2 );

	  MT2Analysis<MT2EstimateTree>* treeHH = new MT2Analysis<MT2EstimateTree>( properName, cfg.regionsSet() );
	  addVariables( treeHH );

	  for( unsigned i=0; i<samples_sig.size(); ++i ) {
	    TString sigSampleName(samples_sig[i].name);
	    if( sigSampleName.Contains("TChiHH") && !sigSampleName.Contains("TChiWH")  ){
	      computeYield( samples_sig[i], cfg, treeHH, mass_neutralino2 );
	    }//created masses
	  }
	  signals.push_back( treeHH );
	}
      }
      


      if( process=="HZ"  ){ 
	//make a loop from 1 to 1000 over the mass of the neutralino to get the full scan
	for( int mass =125; mass<531; mass+=25){
	  //	for( int mass =125; mass<1001; mass+=25){
	  //	    for( int mass =125; mass<1001; mass+=25){
	  mass_neutralino2 = mass;
	  if(mass == 0  ) mass_neutralino2 = 1;
	  if(mass == 125) mass_neutralino2 = 127;

	  std::string properName = Form("SMS_TChiHZ_HToGG_m%d", mass_neutralino2 );
	  MT2Analysis<MT2EstimateTree>* treeHZ = new MT2Analysis<MT2EstimateTree>( properName, cfg.regionsSet() );
	  addVariables( treeHZ );

	  std::string properName_HH0p25 = Form("SMS_TChiHZ_HToGG_m%d_HH0p25", mass_neutralino2 );
	  MT2Analysis<MT2EstimateTree>* treeHZ_HH0p25 = new MT2Analysis<MT2EstimateTree>( properName_HH0p25, cfg.regionsSet() );
	  addVariables( treeHZ_HH0p25 );

	  for( unsigned i=0; i<samples_sig.size(); ++i ) {

	    TString sigSampleName(samples_sig[i].name);

	    if( sigSampleName.Contains("TChi") && !sigSampleName.Contains("TChiWH")  ){
	      
	      if( sigSampleName.Contains("TChiHZ"))
		computeYield( samples_sig[i], cfg, treeHZ, mass_neutralino2 );
		
	      if( sigSampleName.Contains("TChiHH"))
		computeYield( samples_sig[i], cfg, treeHZ_HH0p25, mass_neutralino2 );

	    }
	  
	  
	  }//done loop over sig
	  signals.push_back( treeHZ );
	  signals.push_back( treeHZ_HH0p25 ); 
	}//done masses
      }
	

      if( process== "T2bH"  ){ // normal T2bH
	  
	if( mass_neutralino2 > 0 ){   
	  for( int mass2 =250; mass2<901; mass2+=50){
	  //	  for( int mass2 =125; mass2<210; mass2+=25){
	    for( int mass1 =0 ; mass1< (mass2 - 140); mass1+=50){

	      //  if( mass1 > 320 ) continue; //scan stops there, cut of triangular shape (stop was at 300

	      int massNeut2 = mass2;
	      //	      if( mass2 == 125 ) massNeut2 = 127;
	      int massNeut1 = mass1;
	      //if( mass1 == (mass2-125)) massNeut1 = (mass1-1);
	      if( mass1 == 0 ) massNeut1 = 1;

	      std::cout << "Working on point " << mass1 << " : " << mass2 << std::endl;
	      std::cout << "Working on mass point " << massNeut1 << " : " << massNeut2 << std::endl;

	      std::string  bound = "";
	      if(massNeut2 < 401) 		bound = "400";
	      else if(massNeut2 < 551)		bound = "550";
	      else if(massNeut2 < 601)		bound = "600";
	      else if(massNeut2 < 651)		bound = "650";
	      else if(massNeut2 < 701)		bound = "700";
	      else if(massNeut2 < 751)		bound = "750";
	      else if(massNeut2 < 801)		bound = "800";
	      else if(massNeut2 < 851)		bound = "850";
	      else              		bound = "900";

	      std::cout << bound << std::endl;

	      MT2Analysis<MT2EstimateTree>* treeT2bH = new MT2Analysis<MT2EstimateTree>( Form("SMS_T2bH_mSbottom%d_mLSP%d", massNeut2, massNeut1 ), cfg.regionsSet() );
	      addVariables( treeT2bH );

	      for( unsigned i=0; i<samples_sig.size(); ++i ) {
		TString sigSampleFile(samples_sig[i].file);
		if( (! sigSampleFile.Contains(Form("T2bHto%s", bound.c_str() )) ) &&  !(bound=="400" && sigSampleFile.Contains("T2bH_HToGG_250")  )     )
		  continue;		
		TString sigSampleName(samples_sig[i].name);
		std::cout << "sample name = " << sigSampleName  << std::endl;
		if( sigSampleName.Contains("T2bH")  ){

		  computeYield( samples_sig[i], cfg, treeT2bH, massNeut2, massNeut1 );
    
		}
	      }
	      signals.push_back( treeT2bH );
	    }	  

	  }//done loop over masses for tchwh
	}

	// for( unsigned i=0; i<samples_sig.size(); ++i ) {

	//   TString sigSampleName(samples_sig[i].name);

	//   if( sigSampleName.Contains("T2bH")){

	//     std::string properName = samples_sig[i].sname;
	//     std::string properNameSave = samples_sig[i].sname;


	//     MT2Analysis<MT2EstimateTree>* thisSig = new MT2Analysis<MT2EstimateTree>( samples_sig[i].sname, cfg.regionsSet(), samples_sig[i].id );
	//     addVariables( thisSig );
	//     computeYield( samples_sig[i], cfg, thisSig);
	//     signals.push_back( thisSig );

	// }

      }
   

    }  // for signals



    std::cout << "signal size is " << signals.size() << std::endl;

    //    std::string mcFile_sig = outputdir + "/mc_sig.root";
    std::string mcFile_sig = outputdir + "/mc_"+process+".root";

    for( unsigned i=0; i<signals.size(); ++i )
      if( i==0)
      	signals[i]->writeToFile( mcFile_sig, "RECREATE");
      else
	signals[i]->addToFile( mcFile_sig);

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


  //  float isoCut = cfg.gammaIsoCut();  

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;

  bool isQCD  =   (sample.id>=100 && sample.id<200) ||  (sample.id>=3100 && sample.id<3200 );
  bool isGJets = sample.id>=200 && sample.id<294;
  
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

  bool isT2bH = 0;
  if(sigSampleName.Contains("T2bH"))
    isT2bH = 1;


  TFile* r9corrFile = TFile::Open("transformation_Moriond17_AfterPreApr_v1.root");
  TGraph* r9corr = (TGraph*)r9corrFile->Get("transffull5x5R9EB");
  r9corrFile->Close();

  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( myTree.ngamma<2 ) continue;

    if( isEW && mass_neut2>0 && ( mass_neut2 != myTree.GenSusyMNeutralino2) ) {
      if ( !isEW )
	std::cout << "Why are you doing this" << std::endl;
      continue;
    }

    if( isEW_WH && mass_neut2>0 && ( mass_neut1 != myTree.GenSusyMNeutralino) ) {
      if ( !isEW_WH )
	std::cout << "Why are you doing this" << std::endl;
      continue;
    }

    if( isT2bH && mass_neut2>0 && ( mass_neut2 != myTree.GenSusyMSbottom) ) {
      continue;
    }
    if( isT2bH && mass_neut2>0 && ( mass_neut1 != myTree.GenSusyMNeutralino) ) {
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
      
      if( !is2017 && !(myTree.passFilters()  ) ) continue;
      if( !is2017 && !(myTree.HLT_DiPhoton30 ) ) continue;

      if( is2017 && !(myTree.passFilters2017()   ) ) continue;
      if( is2017 && !(myTree.HLT_DiPhoton30_2017 ) ) continue;
    }else{
      //is MC but different filters for full and fast sim, not trigger in fastsim, ie sig

      if( !(mass_neut2>0) ){ //trigger decision apparently not in signal babies
	if( !is2017 && !(myTree.HLT_DiPhoton30 ) ) continue;
	if( is2017 && !(myTree.HLT_DiPhoton30_2017 ) ) continue;
      }

      if( !(mass_neut2>0) ){
	if( !is2017 && !(myTree.passFiltersMC()  ) ) continue;
	if( is2017 && !(myTree.passFilters2017MC()   ) ) continue;
      }else{
	if( !is2017 && !(myTree.passFiltersMCFastSim()  ) ) continue;
	if( is2017 && !(myTree.passFilters2017MCFastSim()   ) ) continue;
      }

    }


    // if( myTree.gamma_pt[0]<30 ) continue;
    // if( myTree.gamma_pt[1]<20 ) continue;

    // TLorentzVector gamma;
    // gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );

    // //Need the lorentz vectors of the photons first
    TLorentzVector *LVec = new TLorentzVector[myTree.ngamma];

    std::vector<int> indexG;     

    for(int i=0; i< myTree.ngamma; i++){
      //      if(myTree.gamma_chHadIso[i]<=2.5){

      bool isNotALep = true; //already taken care of at heppy level now
      // if( myTree.nlep>0){
      //   for(int lep=0; lep< myTree.nlep; lep++){
      //     //    	    float dR = sqrt( ((myTree.lep_eta[lep]-myTree.gamma_eta[i])*(myTree.lep_eta[lep]-myTree.gamma_eta[i])) + ((myTree.lep_phi[lep]-myTree.gamma_phi[i])*(myTree.lep_phi[lep]-myTree.gamma_phi[i])) );
      //     float dR = DeltaR ( myTree.lep_eta[lep], myTree.gamma_eta[i], myTree.lep_phi[lep], myTree.gamma_phi[i] );
      //     //Clean out photons from electrons with dR<1, from muons wiht dR<0.5
      //     if( !( (abs(myTree.lep_pdgId[lep])==11 && dR>1.0 ) || (abs(myTree.lep_pdgId[lep])==13 && dR>0.5) ) ){ 
      //       isNotALep = false;
      //     }//only take the nice photons
      //   }//done loop over leptons
      // }

      if( myTree.nlep==0 || isNotALep==true ) {
	LVec[i].SetPtEtaPhiM(myTree.gamma_pt[i], myTree.gamma_eta[i],myTree.gamma_phi[i], myTree.gamma_mass[i]);
	indexG.push_back( i );
 
	//    	}
      }
    }

    if( indexG.size()<2 ){
      //    if( ngam<2 || gotSecond==false ){
      //   std::cout << "Kicking out event with " << myTree.ngamma << " of with only " << ngam << " are clean and got a second = " << << std::endl;     
      continue;
    }

    //    if( (isQCD || isGJets) && myTree.gamma_mcMatchId[indexG[0]]==22 && myTree.gamma_mcMatchId[1]==22 )
    if( (isQCD || isGJets) && myTree.gamma_mcMatchId[0]==22 && myTree.gamma_mcMatchId[1]==22 )
      {
	// 	std::cout << "removing pp event from QCD and gJets " << std::endl;
	continue; //remove promptprompt from QCD & GJets as it is already in the BOX
      }
    TLorentzVector Hvec;
    if( indexG.size()>=2 )
      Hvec = LVec[0] + LVec[1]; //gamma invariant mass
    else 
      continue; //should already be taken care of in the skim

    // std::vector<float> input_pt;    // std::vector<float> input_eta;    // std::vector<float> input_phi;   
    // std::vector<float> input_mass;
    // for(int jj=0; jj < myTree.njet; jj++){
    //   if( fabs(myTree.jet_eta[jj]) < 2.4 ){
    // 	input_pt.push_back( myTree.jet_pt[jj] );    // 	input_eta.push_back( myTree.jet_eta[jj] );    
    // 	input_phi.push_back( myTree.jet_phi[jj] );    // 	input_mass.push_back( myTree.jet_mass[jj] );
    //   }
    // }
    // if( ( myTree.nlep>0 )){
    //   for(int jj=0; jj < myTree.nlep; jj++){
    // 	input_pt.push_back( myTree.lep_pt[jj] );    // 	input_eta.push_back( myTree.lep_eta[jj] );
    // 	input_phi.push_back( myTree.lep_phi[jj] );    // 	input_mass.push_back( myTree.lep_mass[jj] );
    //   }
    // }
    // if( ( myTree.nisoTrack>0 )){
    //   for(int jj=0; jj < myTree.nisoTrack; jj++){
    // 	input_pt.push_back( myTree.isoTrack_pt[jj] );    // 	input_eta.push_back( myTree.isoTrack_eta[jj] );
    // 	input_phi.push_back( myTree.isoTrack_phi[jj] );    // 	input_mass.push_back( myTree.isoTrack_mass[jj] );
    //   }
    // }

    //    float minMTBmet = myTree.minMTBMet;
    //    float met       = myTree.gg_met_pt;
    float met       = myTree.met_pt;
    int njets       = myTree.gg_nJet30;
    int nbjets      = myTree.gg_nBJet20;    
    //    int nbjets      = (is2017)  ?  myTree.gg_nBJet20deepcsv : myTree.gg_nBJet20;    
    float ht        = myTree.gg_ht;
    float mt2       = (njets>1) ?  myTree.gg_mt2 : ht;



    TLorentzVector gamma0;
    gamma0.SetPtEtaPhiM( myTree.gamma_pt[indexG[0]], myTree.gamma_etaSc[indexG[0]], myTree.gamma_phi[indexG[0]], myTree.gamma_mass[indexG[0]] );

    TLorentzVector gamma1;
    gamma1.SetPtEtaPhiM( myTree.gamma_pt[indexG[1]], myTree.gamma_etaSc[indexG[1]], myTree.gamma_phi[indexG[1]], myTree.gamma_mass[indexG[1]] );

    // TLorentzVector pseudoJet1;
    // pseudoJet1.SetPtEtaPhiM( myTree.pseudoJet1_pt, myTree.pseudoJet1_eta, myTree.pseudoJet1_phi, myTree.pseudoJet1_mass );
    // TLorentzVector pseudoJet2;
    // pseudoJet2.SetPtEtaPhiM( myTree.pseudoJet2_pt, myTree.pseudoJet2_eta, myTree.pseudoJet2_phi, myTree.pseudoJet2_mass );
    // float scalarProd11 = gamma0.Vect()*pseudoJet1.Vect();
    // float scalarProd12 = gamma1.Vect()*pseudoJet1.Vect();

    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi(); 

    TString anaTreeName( anaTree->getName() );

    if(sigSampleName.Contains("HZ"))
      weight *= 0.5;
    else if( anaTreeName.Contains("HH0p25"))
      weight *= 0.25;

    MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, Hvec.Perp()/Hvec.M() , Hvec.M());
    if( thisTree==0 ) continue;

    float h_mass = -999;
    float h_pt = -999;

    // //This is now already taken care of by heppy
    // if( is2017 && (fabs(myTree.gamma_etaSc[indexG[0]])<=1.4442) && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.0103 ) continue;
    // if( is2017 && (fabs(myTree.gamma_etaSc[indexG[1]])<=1.4442) && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.0103 ) continue;
    // if( is2017 && (fabs(myTree.gamma_etaSc[indexG[0]])>1.4442)  && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.0276 ) continue;
    // if( is2017 && (fabs(myTree.gamma_etaSc[indexG[1]])>1.4442)  && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.0276 ) continue;

    // if( !is2017 && (fabs(myTree.gamma_etaSc[indexG[0]])<=1.4442) && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.01031 ) continue;
    // if( !is2017 && (fabs(myTree.gamma_etaSc[indexG[1]])<=1.4442) && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.01031 ) continue;
    // if( !is2017 && (fabs(myTree.gamma_etaSc[indexG[0]])>1.4442)  && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.03013 ) continue;
    // if( !is2017 && (fabs(myTree.gamma_etaSc[indexG[1]])>1.4442)  && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.03013 ) continue;

    // if( is2017 && (fabs(myTree.gamma_etaSc[indexG[0]])<=1.4442) && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.0105 ) continue;
    // if( is2017 && (fabs(myTree.gamma_etaSc[indexG[1]])<=1.4442) && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.0105 ) continue;
    // if( is2017 && (fabs(myTree.gamma_etaSc[indexG[0]])>1.4442)  && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.0356 ) continue;
    // if( is2017 && (fabs(myTree.gamma_etaSc[indexG[1]])>1.4442)  && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.0356 ) continue;

    // if( !is2017 && (fabs(myTree.gamma_etaSc[indexG[0]])<=1.4442) && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.011 ) continue;
    // if( !is2017 && (fabs(myTree.gamma_etaSc[indexG[1]])<=1.4442) && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.011 ) continue;
    // if( !is2017 && (fabs(myTree.gamma_etaSc[indexG[0]])>1.4442)  && myTree.gamma_sigmaIetaIeta[indexG[0]]> 0.0314 ) continue;
    // if( !is2017 && (fabs(myTree.gamma_etaSc[indexG[1]])>1.4442)  && myTree.gamma_sigmaIetaIeta[indexG[1]]> 0.0314 ) continue;

    h_mass = Hvec.M();
    h_pt = Hvec.Perp();

    if( h_mass > 200 ) continue;
    if( h_mass < 99 ) continue;

    // h_mass = myTree.h_mass;
    // h_pt = myTree.h_pt;

    float diLepMass = -99.;
    bool isDiLepZ = 0;
    int diLepPdgId = -99.;
    //di lep
    if( ( myTree.nlep>1 )){
      if( ((myTree.lep_pdgId[0]*myTree.lep_pdgId[1])< 0 )  && myTree.lep_pt[0]>= 20 && myTree.lep_pt[1]>=20  && ((abs(myTree.lep_pdgId[0])==13 && abs(myTree.lep_pdgId[1])==13) || ( abs(myTree.lep_pdgId[0])==11 && abs(myTree.lep_pdgId[1])==11 )  ) ) {
	//Need the lorentz vectors of the leptons first
	TLorentzVector *LVec = new TLorentzVector[2];
	for(int i=0; i< 2; i++){
	  LVec[i].SetPtEtaPhiM(myTree.lep_pt[i], myTree.lep_eta[i],myTree.lep_phi[i], myTree.lep_mass[i]);
	}
	TLorentzVector Zvec = LVec[0] + LVec[1]; //leptons invariant mass
	diLepMass = Zvec.M(); // compare with myTree.zll_mass
    	if( fabs( diLepMass -91.19) <= 20 ) 
	  isDiLepZ  = 1;
	diLepPdgId = abs ( myTree.lep_pdgId[0] );
      }
    }

    bool is1Mu = 0;
    bool is1El = 0;

    //SINGLE LEPTON region
    if( myTree.nlep==1 ){
      if(  myTree.lep_pt[0]>= 20 ) {
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

    int numberOfCombos = -1;

    float btag_wp = 0.5426; //2016 value 
    if( is2017 )
      btag_wp = 0.5803; //2017 value 

    float bjet0_pt = -99;
    float bjet0_eta = -99;
    float bjet0_phi= -99;
    float bjet0_mass = -99;
    float bjet0_btagCSV = -99;

    float bjet1_pt = -99;
    float bjet1_eta = -99;
    float bjet1_phi= -99;
    float bjet1_mass = -99;
    float bjet1_btagCSV = -99;

    float bjet2_pt = -99;
    float bjet2_eta = -99;
    float bjet2_phi= -99;
    float bjet2_mass = -99;
    float bjet2_btagCSV = -99;

    float bjet3_pt = -99;
    float bjet3_eta = -99;
    float bjet3_phi= -99;
    float bjet3_mass = -99;
    float bjet3_btagCSV = -99;


    if( nbjets > 1 ){
      std::vector<int>* indexB = new std::vector<int>;
      //      int counter = 0;      //      std::cout << "working on event = " << myTree.evt << std::endl;
      //     std::vector<std::unique_ptr<int>> indexB= new std::vector<unique_ptr<s>>;;   
      for( int i=0; i< myTree.njet; i++){
      // for( int i=0; i< (njets+nbjets) ; i++){		
	float dR0 = DeltaR ( myTree.jet_eta[i], myTree.gamma_eta[0], myTree.jet_phi[i], myTree.gamma_phi[0] );
	float dR1 = DeltaR ( myTree.jet_eta[i], myTree.gamma_eta[1], myTree.jet_phi[i], myTree.gamma_phi[1] );
	//	std::cout<< "njets " << myTree.gg_nJet30 << " nbjets " << myTree.gg_nBJet20  <<  " pt = " << myTree.jet_pt[i]  <<  " eta = " << myTree.jet_eta[i]  <<  " btag = " << myTree.jet_btagCSV[i] << " dr0, dr1 =" << dR0 << ", " << dR1 << std::endl;
	if( dR1<0.4  || dR0<0.4)
	  continue;	
	if ( fabs(myTree.jet_eta[i]) > 2.4 ) 
	  continue;

	if( myTree.jet_btagCSV[i] > btag_wp ){
	  indexB->push_back(  i  );
	  // std::cout << " pushing back index " << i << std::endl;
	  //	  counter++;
	}
	//	  indexB->push_back( std::make_unique<int>( i ) );
      }

      // if (indexB->size()> nbjets )
      // 	std::cout << "found in evt " << myTree.evt << " more bjets "<< indexB->size() << " than bjets in the tree " << nbjets << std::endl;
      // if (indexB->size()< nbjets )
      // 	std::cout << "found in evt " << myTree.evt << " less bjets "<< indexB->size() << " than bjets in the tree " << nbjets << std::endl;

      int count = 0;
      for(auto bIndex = indexB->begin(); bIndex != indexB->end(); ++bIndex){
	if(count ==0){
	  bjet0_pt = myTree.jet_pt[*bIndex];
	  bjet0_eta = myTree.jet_eta[*bIndex];
	  bjet0_mass = myTree.jet_mass[*bIndex];
	  bjet0_phi = myTree.jet_phi[*bIndex];
	  bjet0_btagCSV = myTree.jet_btagCSV[*bIndex];
	}else if( count ==1 ){
	  bjet1_pt = myTree.jet_pt[*bIndex];
	  bjet1_eta = myTree.jet_eta[*bIndex];
	  bjet1_mass = myTree.jet_mass[*bIndex];
	  bjet1_phi = myTree.jet_phi[*bIndex];
	  bjet1_btagCSV = myTree.jet_btagCSV[*bIndex];
	}else if( count ==2 ){
	  bjet2_pt = myTree.jet_pt[*bIndex];
	  bjet2_eta = myTree.jet_eta[*bIndex];
	  bjet2_mass = myTree.jet_mass[*bIndex];
	  bjet2_phi = myTree.jet_phi[*bIndex];
	  bjet2_btagCSV = myTree.jet_btagCSV[*bIndex];
	}else if( count ==3 ){
	  bjet3_pt = myTree.jet_pt[*bIndex];
	  bjet3_eta = myTree.jet_eta[*bIndex];
	  bjet3_mass = myTree.jet_mass[*bIndex];
	  bjet3_phi = myTree.jet_phi[*bIndex];
	  bjet3_btagCSV = myTree.jet_btagCSV[*bIndex];
	}

	++count;
      }

      //      if( counter>1 ){
      if( indexB->size()>1  ){
	for(auto bIndex = indexB->begin(); bIndex != indexB->end()-1; ++bIndex){
	  for(auto bIndex2 = indexB->begin()+1; bIndex2 != indexB->end(); ++bIndex2){

	    if( *bIndex != *bIndex2 && ( myTree.jet_pt[*bIndex] > myTree.jet_pt[*bIndex2]  ) ){
	    // std::cout << "Working on b jets of jets at indeces " << *bIndex << ", " << *bIndex2 << std::endl;
	    // std::cout << "Working on b jets of jets at indeces " << myTree.jet_pt[*bIndex]<< ", " << myTree.jet_pt[*bIndex2] << std::endl;
	    // std::cout << "Working on b jets of jets at indeces " << myTree.jet_eta[*bIndex]<< ", " << myTree.jet_eta[*bIndex2] << std::endl;
	    // // std::cout << "Working on b jets of jets at indeces " << myTree.jet_btagCSV[*bIndex]<< ", " << myTree.jet_btagCSV[*bIndex2] << std::endl;

	    TLorentzVector *LVec = new TLorentzVector[2];
	    LVec[0].SetPtEtaPhiM(myTree.jet_pt[*bIndex], myTree.jet_eta[*bIndex], myTree.jet_phi[*bIndex], myTree.jet_mass[*bIndex]);
	    LVec[1].SetPtEtaPhiM(myTree.jet_pt[*bIndex2], myTree.jet_eta[*bIndex2], myTree.jet_phi[*bIndex2], myTree.jet_mass[*bIndex2]);
	    TLorentzVector Bvec = LVec[0] + LVec[1]; //bb invariant mass
	    diBMass = Bvec.M();
	    bMassVector.push_back( diBMass );
	    bPTVector.push_back( Bvec.Perp() );
	    //	    std::cout << "FOUND A BB pair with mass "  << diBMass << std::endl;
	    diBpT  = Bvec.Perp();
	    //	if( (fabs( diBMass-91.19) <= 20.) ||  (fabs(diBMass-125) <= 20.)  ) 
	    //isDiB = 1;

	    //	    std::cout << "yielding mass =  " << diBMass << std::endl;
	    }
	  }
	}

	// diBMass  = bMassVector[0]; //set to first one then see if there is a better one	// //	diBMass0 = bMassVector[0]; // first in list 	// diBpT  = bPTVector[0]; //set to first one then see if there is a better one 	// //	diBpT0 = bPTVector[0]; // first in list
	// //alternative method: only check the leading pair	// if( (diBMass <140) && (diBMass >= 95) //   isDiBH = 1;	// else  if( (diBMass <95) && (diBMass >= 60) )	//   isDiBZ = 1;

	numberOfCombos = bMassVector.size();
  
	for( unsigned int bbCand=0; bbCand<bMassVector.size(); bbCand++){
	  if( (bMassVector[ bbCand ] <140) && (bMassVector[ bbCand ] >= 95) ){
	    isDiBH = 1;
	    diBMass  = bMassVector[bbCand];
	    diBpT  = bPTVector[bbCand];
	  }
	}
	for( unsigned int bbCand=0; bbCand<bMassVector.size(); bbCand++){
	  if( (bMassVector[ bbCand ] <95) && (bMassVector[ bbCand ] >= 60)  && (isDiBH==0) ){
	    isDiBZ = 1;
	    diBMass  = bMassVector[bbCand];
	    diBpT  = bPTVector[bbCand];
	  }
	}

      } 

      //      indexB->erase( std::next( indexB->begin() ));
      indexB->clear();

    }// only do all this if there are 2 bjets

    //    int nwideEtaGamma = ngamma;

    // ECAL prefiring: filter out events that have a jet in 2.25<eta<3 && pt>100 or photon pt>50
    bool passPrefire = 1;

    if( !myTree.isData ){

      if( is2017 )
	for(int i=0; i< 20; i++){
	  if ( fabs(myTree.gammaWideEta_eta[i]) > 2.24 && fabs(myTree.gammaWideEta_eta[i]) < 3.0 && myTree.gammaWideEta_pt[i]>50 ) 
	    passPrefire = 0;
	}
      else
	for(int i=0; i< myTree.ngamma; i++){
	  if ( fabs(myTree.gamma_eta[i]) > 2.24 && fabs(myTree.gamma_eta[i]) < 3.0 && myTree.gamma_pt[i]>50 ) 
	    passPrefire = 0;
	}

      if ( passPrefire ) // only need to loop over jets if it wasn't already flagged as bad
	for( int i=0; i< myTree.njet; i++){
	  if ( fabs(myTree.jet_eta[i]) > 2.24 && fabs(myTree.jet_eta[i]) < 3.0 && myTree.jet_pt[i]>100 ) 
	    passPrefire = 0;
	}
    }


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


    thisTree->assignVar( "passPrefire", passPrefire );
    thisTree->assignVar( "ngamma", myTree.ngamma );

    thisTree->assignVar( "hgg_mt2", myTree.hgg_mt2 );
    thisTree->assignVar( "hgg_mt2_genMET", myTree.hgg_mt2_genMET );



    thisTree->assignVar( "bjet0_pt", bjet0_pt );
    thisTree->assignVar( "bjet0_eta", bjet0_eta );
    thisTree->assignVar( "bjet0_mass", bjet0_mass );
    thisTree->assignVar( "bjet0_phi", bjet0_phi );
    thisTree->assignVar( "bjet0_btagCSV", bjet0_btagCSV );

    thisTree->assignVar( "bjet1_pt", bjet1_pt );
    thisTree->assignVar( "bjet1_eta", bjet1_eta );
    thisTree->assignVar( "bjet1_mass", bjet1_mass );
    thisTree->assignVar( "bjet1_phi", bjet1_phi );
    thisTree->assignVar( "bjet1_btagCSV", bjet1_btagCSV );

    thisTree->assignVar( "bjet2_pt", bjet2_pt );
    thisTree->assignVar( "bjet2_eta", bjet2_eta );
    thisTree->assignVar( "bjet2_mass", bjet2_mass );
    thisTree->assignVar( "bjet2_phi", bjet2_phi );
    thisTree->assignVar( "bjet2_btagCSV", bjet2_btagCSV );

    thisTree->assignVar( "bjet3_pt", bjet3_pt );
    thisTree->assignVar( "bjet3_eta", bjet3_eta );
    thisTree->assignVar( "bjet3_mass", bjet3_mass );
    thisTree->assignVar( "bjet3_phi", bjet3_phi );
    thisTree->assignVar( "bjet3_btagCSV", bjet3_btagCSV );


    thisTree->assignVar("nVert", myTree.nVert );

    thisTree->assignVar( "lep_pdgId0", lep_id );
    thisTree->assignVar( "lep_tightId0", lep_tightId0 );
    thisTree->assignVar( "lep_tightId1", lep_tightId1 );
    thisTree->assignVar( "lep_eta0", lep_eta0 );
    thisTree->assignVar( "lep_eta1", lep_eta1 );
    thisTree->assignVar( "lep_phi0", lep_phi0 );
    thisTree->assignVar( "lep_phi1", lep_phi1 );

    thisTree->assignVar( "diLepMass", diLepMass );
    thisTree->assignVar( "isDiLepZ", isDiLepZ );
    thisTree->assignVar( "diLepId", diLepPdgId );

    thisTree->assignVar( "is1El", is1El);
    thisTree->assignVar( "is1Mu", is1Mu);

    thisTree->assignVar( "diBMass", diBMass );
    thisTree->assignVar( "isDiBH", isDiBH);
    thisTree->assignVar( "isDiBZ", isDiBZ );
    thisTree->assignVar( "diBpT", diBpT );
    thisTree->assignVar( "diBCombos", numberOfCombos );


    thisTree->assignVar( "met_phi", myTree.met_phi );

    thisTree->assignVar( "h_phi", Hvec.Phi() );
    thisTree->assignVar( "h_eta", Hvec.Eta() );
    thisTree->assignVar( "h_mass", h_mass );
    thisTree->assignVar( "h_pt",   h_pt   );

    thisTree->assignVar( "ptGamma0", myTree.gamma_pt[0] );
    thisTree->assignVar( "etaGamma0", myTree.gamma_etaSc[0] );
    thisTree->assignVar( "phiGamma0", myTree.gamma_phi[0] );

    float gamma0_r9_corr = (myTree.isData) ? myTree.gamma_r9[0] : r9corr->Eval(myTree.gamma_r9[0]);
    float gamma1_r9_corr = (myTree.isData) ? myTree.gamma_r9[1] : r9corr->Eval(myTree.gamma_r9[1]);

    thisTree->assignVar( "r9Gamma0", gamma0_r9_corr );
    thisTree->assignVar( "r9Gamma1", gamma1_r9_corr );
    thisTree->assignVar( "sigmaIetaIetaGamma0", myTree.gamma_sigmaIetaIeta[0] );

    // thisTree->assignVar( "chHadIsoGamma0", myTree.gamma_chHadIso[0] );
    // thisTree->assignVar( "chHadIsoRCGamma0", myTree.gamma_chHadIsoRC[0] );
    // thisTree->assignVar( "chHadIsoRC04Gamma0", myTree.gamma_chHadIsoRC04[0] );

    thisTree->assignVar( "ptGamma1", myTree.gamma_pt[1] );
    thisTree->assignVar( "etaGamma1", myTree.gamma_etaSc[1] );
    thisTree->assignVar( "phiGamma1", myTree.gamma_phi[1] );
    //    thisTree->assignVar( "r9Gamma1", myTree.gamma_r9[1] );
    thisTree->assignVar( "sigmaIetaIetaGamma1", myTree.gamma_sigmaIetaIeta[1] );




    // thisTree->assignVar( "ptGamma0", myTree.gamma_pt[indexG[0]] );
    // thisTree->assignVar( "etaGamma0", myTree.gamma_etaSc[indexG[0]] );
    // thisTree->assignVar( "phiGamma0", myTree.gamma_phi[indexG[0]] );

    // float gamma0_r9_corr = (myTree.isData) ? myTree.gamma_r9[indexG[0]] : r9corr->Eval(myTree.gamma_r9[indexG[0]]);
    // float gamma1_r9_corr = (myTree.isData) ? myTree.gamma_r9[indexG[1]] : r9corr->Eval(myTree.gamma_r9[indexG[1]]);

    // thisTree->assignVar( "r9Gamma0", gamma0_r9_corr );
    // thisTree->assignVar( "r9Gamma1", gamma1_r9_corr );
    // thisTree->assignVar( "sigmaIetaIetaGamma0", myTree.gamma_sigmaIetaIeta[indexG[0]] );

    // // thisTree->assignVar( "chHadIsoGamma0", myTree.gamma_chHadIso[indexG[0]] );
    // // thisTree->assignVar( "chHadIsoRCGamma0", myTree.gamma_chHadIsoRC[indexG[0]] );
    // // thisTree->assignVar( "chHadIsoRC04Gamma0", myTree.gamma_chHadIsoRC04[indexG[0]] );

    // thisTree->assignVar( "ptGamma1", myTree.gamma_pt[indexG[1]] );
    // thisTree->assignVar( "etaGamma1", myTree.gamma_etaSc[indexG[1]] );
    // thisTree->assignVar( "phiGamma1", myTree.gamma_phi[indexG[1]] );
    // //    thisTree->assignVar( "r9Gamma1", myTree.gamma_r9[indexG[1]] );
    // thisTree->assignVar( "sigmaIetaIetaGamma1", myTree.gamma_sigmaIetaIeta[indexG[1]] );

    // thisTree->assignVar( "chHadIsoGamma1", myTree.gamma_chHadIso[indexG[1]] );
    // thisTree->assignVar( "chHadIsoRCGamma1", myTree.gamma_chHadIsoRC[indexG[1]] );
    // thisTree->assignVar( "chHadIsoRC04Gamma1", myTree.gamma_chHadIsoRC04[indexG[1]] );

    //    thisTree->assignVar( "gamma_jet1_pt", myTree.gamma_jet1_pt );
    //    thisTree->assignVar( "gamma_jet2_pt", myTree.gamma_jet2_pt );
    // thisTree->assignVar( "gamma_chHadIsoRC",  myTree.gamma_chHadIsoRC[0] );
    thisTree->assignVar( "drMinParton", myTree.gamma_drMinParton[0] );
    //thisTree->assignVars( ht, njets, nbjets, met, mt2 ); //manually override these (so change nbjets)
    thisTree->assignVar( "raw_mt2", myTree.mt2 );

    thisTree->assignVar( "gg_nJets", myTree.gg_nJet30 );
    thisTree->assignVar( "gg_nBJets", nbjets );
    thisTree->assignVar( "gg_nBJetsCSV", myTree.gg_nBJet20 );

    thisTree->assignVar( "gg_mt2", myTree.gg_mt2 );
    thisTree->assignVar( "gg_ht", myTree.gg_ht );
    thisTree->assignVar( "gg_deltaPhiMin", myTree.gg_deltaPhiMin );
    thisTree->assignVar( "gg_met", myTree.gg_met_pt );
    thisTree->assignVar( "gg_met_phi", myTree.gg_met_phi );

    float normalMT2 = myTree.hgg_mt2;
    float normalMET = myTree.met_pt;
    if( is2017 ){
      normalMET = myTree.met_stdMET_pt;
      normalMT2 = myTree.hgg_stdMET_mt2;
    }
    thisTree->assignVar( "std_met", normalMET );
    thisTree->assignVar( "std_mt2", normalMT2 );



    thisTree->assignVar( "met_pt", myTree.met_pt );


    // thisTree->assignVar( "gg_jet1_pt", myTree.gg_jet1_pt );
    // thisTree->assignVar( "gg_jet2_pt", myTree.gg_jet2_pt );

    // thisTree->assignVar( "jet1_pt", myTree.jet1_pt );
    // thisTree->assignVar( "jet2_pt", myTree.jet2_pt );

    if(sigSampleName.Contains("T2bH") || sigSampleName.Contains("TChi") ){ 

      //only apply the 2016 weight correction to 2016, for 2017, weight =1 and uncert =10
      if( !is2017 ){
	thisTree->assignVar( "weight_isr", myTree.weight_isr/myTree.weight_isr_av );
	thisTree->assignVar( "weight_isr_UP", myTree.weight_isr_UP/myTree.weight_isr_UP_av );
	thisTree->assignVar( "weight_isr_DN", myTree.weight_isr_DN/myTree.weight_isr_DN_av );
	weight *=  myTree.weight_isr/myTree.weight_isr_av;
      }else{
	thisTree->assignVar( "weight_isr", 1. );
	thisTree->assignVar( "weight_isr_UP", 1.1);
	thisTree->assignVar( "weight_isr_DN", 1.1);
      }

    }

    if( !(myTree.isData) ){

      // scale diphoton xsec to nlo, only affects data to mc plots
      if( sample.id==295)
	weight *= 135.1/84.4;

      float weight_btag_av = myTree.weight_btagsf_av;
      if( myTree.weight_btagsf_av == 0 )
	weight_btag_av = 1.;

      float weight_btag = myTree.weight_btagsf;
      if( myTree.weight_btagsf == 0 )
	weight_btag = 1.;

      float weight_lepsf = myTree.weight_lepsf2017;
      float weight_lepsf_UP = myTree.weight_lepsf2017_UP;
      float weight_lepsf_DN = myTree.weight_lepsf2017_DN;
      if( myTree.weight_lepsf2017 == 0 ){
	weight_lepsf = 1.;
	weight_lepsf_UP = 1. + weight_lepsf_UP;
	weight_lepsf_DN = 1. + weight_lepsf_DN;
      }

      float weight_gammasf = myTree.weight_gammasf2017;
      float weight_gammasf_UP = myTree.weight_gammasf2017_UP;
      float weight_gammasf_DN = myTree.weight_gammasf2017_DN;
      if( myTree.weight_gammasf2017 == 0 ){	
	weight_gammasf = 1.;
	weight_gammasf_UP = 1. + weight_gammasf_UP;
	weight_gammasf_DN = 1. + weight_gammasf_DN;
      }


      //only need it for SMH and signal
      if( sample.id>950 )
	for( int i=0; i<110; i++)
	  thisTree->assignVar( Form("weight_scales_%d",i), myTree.weight_scales[i]/myTree.weight_scales_av[i] );

      //since we do yield = 1/2(yield_met + yield_genMET) we will just halve the weight before so that I can just add the datasets
      if( sample.id>=1000 )
	weight *= 0.5;

      //      weight *= myTree.weight_btagsf / myTree.weight_btagsf_av * myTree.weight_lepsf2017 ;
      //      weight *= myTree.weight_btagsf / myTree.weight_btagsf_av * myTree.weight_lepsf2017  * myTree.weight_gammasf2017;
      weight *= weight_btag / weight_btag_av * weight_lepsf  * weight_gammasf;

      thisTree->assignVar( "weight_lepsf", weight_lepsf );
      thisTree->assignVar( "weight_lepsf_UP", weight_lepsf_UP );
      thisTree->assignVar( "weight_lepsf_DN", weight_lepsf_DN );

      thisTree->assignVar( "weight_gammasf", weight_gammasf );
      thisTree->assignVar( "weight_gammasf_UP", weight_gammasf_UP );
      thisTree->assignVar( "weight_gammasf_DN", weight_gammasf_DN );

      thisTree->assignVar( "weight_btagsf", weight_btag/myTree.weight_btagsf_av );
      thisTree->assignVar( "weight_btagsf_light_UP", myTree.weight_btagsf_light_UP/myTree.weight_btagsf_light_UP_av );
      thisTree->assignVar( "weight_btagsf_light_DN", myTree.weight_btagsf_light_DN/myTree.weight_btagsf_light_DN_av );
      thisTree->assignVar( "weight_btagsf_heavy_UP", myTree.weight_btagsf_heavy_UP/myTree.weight_btagsf_heavy_UP_av );
      thisTree->assignVar( "weight_btagsf_heavy_DN", myTree.weight_btagsf_heavy_DN/myTree.weight_btagsf_heavy_DN_av );

    }


    thisTree->assignTree(myTree, weight ); // assign all to defaults
 
    thisTree->yield->Fill( h_mass , weight );
    thisTree->tree->Fill();
    
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
  thisTree->assignVar( "etaGamma", myTree.gamma_etaSc[0] );
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


  MT2EstimateTree::addVar( tree, "passPrefire" );

  MT2EstimateTree::addVar( tree, "ngamma" );

  MT2EstimateTree::addVar( tree, "hgg_mt2" );
  MT2EstimateTree::addVar( tree, "hgg_mt2_genMET" );

  MT2EstimateTree::addVar( tree, "bjet0_pt" );
  MT2EstimateTree::addVar( tree, "bjet1_pt" );
  MT2EstimateTree::addVar( tree, "bjet2_pt" );
  MT2EstimateTree::addVar( tree, "bjet3_pt" );

  MT2EstimateTree::addVar( tree, "bjet0_eta" );
  MT2EstimateTree::addVar( tree, "bjet1_eta" );
  MT2EstimateTree::addVar( tree, "bjet2_eta" );
  MT2EstimateTree::addVar( tree, "bjet3_eta" );

  MT2EstimateTree::addVar( tree, "bjet0_phi" );
  MT2EstimateTree::addVar( tree, "bjet1_phi" );
  MT2EstimateTree::addVar( tree, "bjet2_phi" );
  MT2EstimateTree::addVar( tree, "bjet3_phi" );
 
  MT2EstimateTree::addVar( tree, "bjet0_mass" );
  MT2EstimateTree::addVar( tree, "bjet1_mass" );
  MT2EstimateTree::addVar( tree, "bjet2_mass" );
  MT2EstimateTree::addVar( tree, "bjet3_mass" );
 
  MT2EstimateTree::addVar( tree, "bjet0_btagCSV" );
  MT2EstimateTree::addVar( tree, "bjet1_btagCSV" );
  MT2EstimateTree::addVar( tree, "bjet2_btagCSV" );
  MT2EstimateTree::addVar( tree, "bjet3_btagCSV" );

  // MT2EstimateTree::addVar( tree, "gamma_jet1_pt" );
  // MT2EstimateTree::addVar( tree, "gamma_jet2_pt" );

  MT2EstimateTree::addVar( tree, "drMinParton" );
  MT2EstimateTree::addVar( tree, "raw_mt2" );

  MT2EstimateTree::addVar( tree, "met_phi" );
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
  MT2EstimateTree::addVar( tree, "gg_nBJetsCSV" );

  MT2EstimateTree::addVar( tree, "gg_ht" );
  MT2EstimateTree::addVar( tree, "gg_mt2" );
  MT2EstimateTree::addVar( tree, "gg_deltaPhiMin" );
  //  MT2EstimateTree::addVar( tree, "gg_diffMetMht" );
  // MT2EstimateTree::addVar( tree, "gg_mht" );
  MT2EstimateTree::addVar( tree, "gg_met" );
  MT2EstimateTree::addVar( tree, "gg_met_phi" );
  //  MT2EstimateTree::addVar( tree, "gg_jet1_pt" );
  //  MT2EstimateTree::addVar( tree, "gg_jet2_pt" );

  // MT2EstimateTree::addVar( tree, "jet1_pt" );
  // MT2EstimateTree::addVar( tree, "jet2_pt" );

  //   Float_t weight_scales[110];

  for( int i=0; i<110; i++)
    MT2EstimateTree::addVar( tree, Form("weight_scales_%d",i) );


  MT2EstimateTree::addVar( tree, "weight_isr" );
  MT2EstimateTree::addVar( tree, "weight_isr_UP" );
  MT2EstimateTree::addVar( tree, "weight_isr_DN" );

  MT2EstimateTree::addVar( tree, "weight_lepsf" );
  MT2EstimateTree::addVar( tree, "weight_lepsf_UP" );
  MT2EstimateTree::addVar( tree, "weight_lepsf_DN" );

  MT2EstimateTree::addVar( tree, "weight_gammasf" );
  MT2EstimateTree::addVar( tree, "weight_gammasf_UP" );
  MT2EstimateTree::addVar( tree, "weight_gammasf_DN" );

  MT2EstimateTree::addVar( tree, "weight_btagsf" );
  MT2EstimateTree::addVar( tree, "weight_btagsf_light_UP" );
  MT2EstimateTree::addVar( tree, "weight_btagsf_light_DN" );
  MT2EstimateTree::addVar( tree, "weight_btagsf_heavy_UP" );
  MT2EstimateTree::addVar( tree, "weight_btagsf_heavy_DN" );
    

  MT2EstimateTree::addVar( tree, "isDiBH" );
  MT2EstimateTree::addVar( tree, "isDiBZ" );
  MT2EstimateTree::addVar( tree, "diBMass" );
  MT2EstimateTree::addVar( tree, "isDiLepZ" );
  MT2EstimateTree::addVar( tree, "diLepMass" );  
  MT2EstimateTree::addVar( tree, "diLepId" ); 

  MT2EstimateTree::addVar( tree, "is1El" );
  MT2EstimateTree::addVar( tree, "is1Mu" );

  MT2EstimateTree::addVar( tree, "diBpT" );
  MT2EstimateTree::addVar( tree, "diBCombos");

  MT2EstimateTree::addVar( tree, "lep_pdgId0");
  // MT2EstimateTree::addVar( tree, "lep_pdgId1");

  MT2EstimateTree::addVar( tree, "lep_tightId0" );
  MT2EstimateTree::addVar( tree, "lep_tightId1" );

  MT2EstimateTree::addVar( tree, "lep_eta0" );
  MT2EstimateTree::addVar( tree, "lep_eta1" );

  MT2EstimateTree::addVar( tree, "lep_phi0" );
  MT2EstimateTree::addVar( tree, "lep_phi1" );

  MT2EstimateTree::addVar( tree, "nVert" );

  MT2EstimateTree::addVar( tree, "std_met" );
  MT2EstimateTree::addVar( tree, "std_mt2" );
  MT2EstimateTree::addVar( tree, "met_pt" );

}






// void getHemispheres(std::vector<float> input_pt, std::vector<float> input_eta, std::vector<float> input_phi, std::vector<float> input_mass, TLorentzVector *v1, TLorentzVector *v2){
  
//   std::vector<float> px, py, pz, E;
//   int NJets = input_pt.size();
//   for(int j=0; j<NJets; ++j){
//     TLorentzVector jet;
//     //    if (input_pt.at(j) > 30 && fabs(jEtaReb[j])<2.5)
    
//     // selected relevant jets, b-jets 
//     jet.SetPtEtaPhiM(input_pt.at(j), input_eta.at(j), input_phi.at(j), input_mass.at(j) );
//     //else
//     //  continue;
//     px.push_back(jet.Px());
//     py.push_back(jet.Py());
//     pz.push_back(jet.Pz());
//     E .push_back(jet.E ());
//   }

//   Hemisphere* hemisp = new Hemisphere(px, py, pz, E, 2, 3);
//   std::vector<int> grouping = hemisp->getGrouping();

//   v1->SetPxPyPzE(0.,0.,0.,0.);
//   v2->SetPxPyPzE(0.,0.,0.,0.);
//   for(unsigned int i=0; i<px.size(); ++i){
// 	if(grouping[i]==1){
// 		v1->SetPx(v1->Px() + px[i]);
// 		v1->SetPy(v1->Py() + py[i]);
// 		v1->SetPz(v1->Pz() + pz[i]);
// 		v1->SetE (v1->E () + E [i]);	
// 	}else if(grouping[i] == 2){
// 		v2->SetPx(v2->Px() + px[i]);
// 		v2->SetPy(v2->Py() + py[i]);
// 		v2->SetPz(v2->Pz() + pz[i]);
// 		v2->SetE (v2->E () + E [i]);
// 	}
//   }
//   delete hemisp;

// }






// Float_t calcMT2(std::vector<float> input_pt, std::vector<float> input_eta, std::vector<float> input_phi, std::vector<float> input_mass, float met, float metphi){

//   TLorentzVector *visible1 = new TLorentzVector(0.,0.,0.,0.);
//   TLorentzVector *visible2 = new TLorentzVector(0.,0.,0.,0.);

//   getHemispheres(input_pt, input_eta, input_phi, input_mass, visible1, visible2);

//   double pa[3];
//   double pb[3];
//   double pmiss[3];

//   pmiss[0] = 0;
//   pmiss[1] = static_cast<double> (met*TMath::Cos(metphi));
//   pmiss[2] = static_cast<double> (met*TMath::Sin(metphi));

//   pa[0] = 0;
//   pa[1] = static_cast<double> (visible1->Px());
//   pa[2] = static_cast<double> (visible1->Py());
  
//   pb[0] = 0;
//   pb[1] = static_cast<double> (visible2->Px());
//   pb[2] = static_cast<double> (visible2->Py());

//   Davismt2 *mt2 = new Davismt2();
//   mt2->set_momenta(pa, pb, pmiss);
//   mt2->set_mn(0);
//   Float_t MT2=mt2->get_mt2();
//   delete mt2;
//   delete visible1;
//   delete visible2;
//   return MT2;

// }
