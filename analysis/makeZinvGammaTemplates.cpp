#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"








void computeYield( const MT2Sample& sample, const MT2Config& cfg, MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* fake );
void setPoissonError( MT2Analysis<MT2EstimateZinvGamma>* data );
//void randomizePoisson( MT2Analysis<MT2EstimateZinvGamma>* data );




int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|          Running makeZinvGammaTemplates            |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


 
  if( argc<2 ) {
    std::cout << "USAGE: ./makeZinvGammaTemplates [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool useMC = false;

  if( argc>2 ) {
    std::string data_or_mc = std::string(argv[2]); 
    if( data_or_mc=="data" || data_or_mc=="Data" || data_or_mc=="DATA" ) {
      useMC=false;
    } else if( data_or_mc=="mc" || data_or_mc=="MC" ) {
      useMC=true;
    } else {
      std::cout << std::endl;
      std::cout << "-> WARNING! Second argument should be 'data' or 'MC'." << std::endl;
      std::cout << "Exiting." << std::endl;
      std::cout << std::endl;
      exit(817);
    }
  } 


  std::string templateType = cfg.gammaTemplateType();
 
  if( argc>3 ) {

    templateType = std::string(argv[3]); 
    std::cout << std::endl;
    std::cout << "-> Will disobey the cfg and use templateType = " << argv[2] << std::endl;
    std::cout << std::endl;

  } 
 cfg.set_gammaTemplateType(templateType);



  if( templateType!="FR" && templateType!="MC" && templateType!="RC" ) {
    std::cout << "ERROR! templateType may only be 'MC' or 'FR' or 'RC'" << std::endl;
    exit(1111);
  }




  std::cout << std::endl;
  std::cout << "-> Starting to build templates with:" << std::endl;
  std::cout << "      type   : " << templateType << std::endl;
  std::cout << "      regions: " << cfg.gammaTemplateRegions() << std::endl;
  if( useMC )
    std::cout << "      using MC " << std::endl;
  else
    std::cout << "      using data " << std::endl;
  std::cout << std::endl << std::endl;



  std::string regionsSet = cfg.gammaTemplateRegions();


  std::string samplesName = (useMC) ? cfg.mcSamples() : cfg.dataSamples();
  std::string samplesFile = "../samples/samples_" + samplesName + ".dat";
  
  std::vector<MT2Sample> samples = (useMC) ? MT2Sample::loadSamples(samplesFile, 100, 299) : MT2Sample::loadSamples(samplesFile, "SinglePhoton");
  if( samples.size()==0 ) {
    std::cout << "There must be an error: didn't find any good files in " << samplesFile << "!" << std::endl;
    exit(1209);
  }



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  //std::string outputdir = cfg.getEventYieldDir() + "/gammaControlRegion/purity"; //"ZinvGammaPurity_" + samplesFileName + "_" + regionsSet;
  //system(Form("mkdir -p %s", outputdir.c_str()));


  
  MT2Analysis<MT2EstimateZinvGamma>* templatesPrompt = new MT2Analysis<MT2EstimateZinvGamma>( "templatesPrompt", regionsSet );
  MT2Analysis<MT2EstimateZinvGamma>* templatesFake   = new MT2Analysis<MT2EstimateZinvGamma>( "templatesFake", regionsSet );

  for( unsigned i=0; i<samples.size(); ++i ) {
    computeYield( samples[i], cfg, templatesPrompt, templatesFake );
  }


  if( templateType=="FR" || templateType=="RC" ) {
    setPoissonError( templatesFake );
    setPoissonError( templatesPrompt );
    if( templateType=="FR" ) templatesPrompt->setName("templatesPromptRaw");
  }


  std::string templateFileName = cfg.getGammaCRdir() + "/gammaTemplates" + templateType;
  if( useMC ) templateFileName = templateFileName + "_MC";
  else        templateFileName = templateFileName + "_data";
  templateFileName = templateFileName + ".root";


  templatesFake->writeToFile(templateFileName);
  templatesPrompt->addToFile(templateFileName, true);


  return 0;

}









void computeYield( const MT2Sample& sample, const MT2Config& cfg, MT2Analysis<MT2EstimateZinvGamma>* prompt, MT2Analysis<MT2EstimateZinvGamma>* fake ) {


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

    if( !(myTree.HLT_Photon165_HE10) ) continue;

    if( myTree.isData ) {
      if( !( myTree.isGolden ) ) continue;
      if( !(myTree.passFilters() ) ) continue;
    }
    //   if( myTree.isData && !(myTree.HLT_Photon165_HE10) ) continue;

    if( myTree.mt2>200. ) continue; // orthogonal to signal region                                                                                                                                               
    if( myTree.gamma_pt[0]<180. ) continue;
    if( (myTree.gamma_nJet30>1 && myTree.gamma_mt2<200.) || (myTree.gamma_nJet30==1 && myTree.gamma_ht<200.) ) continue;



    // REMOVE THIS SOOOOON
    // or maybe not due to spiky behavior
    if( myTree.evt_scale1fb>10. ) continue;



    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );


    int njets       = myTree.gamma_nJet30;
    int nbjets      = (sample.id>10) ?  myTree.gamma_nBJet20 :  myTree.gamma_nBJet20csv; 
    //  int nbjets      = myTree.gamma_nBJet20csv;    
    float ht        = myTree.gamma_ht;
    float mt2       = (njets>1) ? myTree.gamma_mt2 : ht;
    float minMTBmet = myTree.gamma_minMTBMet;

    //TEEEEEMPORARYYYYY
    //    if( myTree.isData && myTree.run>275125.) continue;

    //temp fix
    myTree.nBJet20csv=nbjets;



    int nJetHF30_ = 0;
    for(int j=0; j<myTree.njet; ++j){
      
      if( myTree.jet_pt[j] < 30. || fabs(myTree.jet_eta[j]) < 3.0 ) continue;
      else ++nJetHF30_;

    }
//    //HF Veto
//    if( nJetHF30_ >0 ) continue;

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
    bool isWorkingFake   = false;

    if( cfg.gammaTemplateType()=="MC" ) {

      if( !sietaietaOK ) continue;
      
      isWorkingPrompt = ( myTree.gamma_mcMatchId[0]==22 && myTree.gamma_drMinParton[0]>0.4); // prompt = matched //no fakes
      isWorkingFake   = ( myTree.gamma_mcMatchId[0]!=22 && myTree.gamma_mcMatchId[0]!=7 && sample.id>=100 && sample.id<199 );
      // isWorkingPrompt = myTree.gamma_mcMatchId[0]==22; // prompt = matched

    } else if( cfg.gammaTemplateType()=="FR" ) { 

      isWorkingPrompt = sietaietaOK;
      isWorkingFake = !isWorkingPrompt;

      if( isWorkingPrompt ) {
        // for prompts use only low-sensitivity regions:
        if( mt2>300. ) continue;
        if( ht>1000. ) continue;
        if( njets>6 ) continue;
        if( nbjets>0 ) continue;
      }

    } else if( cfg.gammaTemplateType()=="RC" ) { 

      isWorkingPrompt = sietaietaOK;
      isWorkingFake = !isWorkingPrompt;

      if( isWorkingPrompt ) iso = myTree.gamma_chHadIsoRC[0]; // random cone


    } else {

      // this shouldnt be possible
      std::cout << "-> NOPE. Don't know anything about cfg.gammaTemplateType()='" << cfg.gammaTemplateType() << "'." << std::endl;
      exit(191);

    }


    //// preselection
    //if( iso > 20. ) continue;
    ////if( iso > 10. ) continue;

    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi(); 


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

    if( isWorkingPrompt ) {

      //      MT2EstimateZinvGamma* thisPrompt = prompt->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      MT2EstimateZinvGamma* thisPrompt = prompt->get( ht, njets, nbjets, minMTBmet, mt2 );
      if( thisPrompt==0 ) continue;

      thisPrompt->yield->Fill(mt2, weight );
      thisPrompt->sietaieta->Fill(myTree.gamma_sigmaIetaIeta[0], weight );
      thisPrompt->fillIso( iso, weight, mt2 );

    } else if( isWorkingFake ) {

      //      MT2EstimateZinvGamma* thisFake = fake->get( ht, njets, nbjets, met, minMTBmet, mt2 );
      MT2EstimateZinvGamma* thisFake = fake->get( ht, njets, nbjets, minMTBmet, mt2 );
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
