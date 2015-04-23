#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>


#include "../interface/MT2Analysis.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Sample.h"

#define mt2_cxx
#include "../interface/mt2.h"


#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"


class MT2Config {
 public:
  MT2Config( const std::string& configFileName );
  std::string regionsSet()      const { return regionsSet_; };
  std::string mcSamples()       const { return mcSamples_; };
  std::string sigSamples()      const { return sigSamples_; };
  std::string dataSamples()     const { return dataSamples_; };
  std::string lostLeptonTag()   const { return lostLeptonTag_; };
  std::string qcdTag()          const { return qcdTag_; };
  std::string zinvTag()         const { return zinvTag_; };
  std::string additionalStuff() const { return additionalStuff_; };

  bool useMC() {
    bool useEstimates = lostLeptonTag_!="" && qcdTag_!="" && zinvTag_!="";
    return !useEstimates; }

 private:

  std::string regionsSet_;
  std::string mcSamples_;
  std::string sigSamples_;
  std::string dataSamples_;
  std::string lostLeptonTag_;
  std::string qcdTag_;
  std::string zinvTag_;
  std::string additionalStuff_;
};


float lumi = 4.;


MT2Analysis<MT2EstimateTree>* computeYield( const MT2Sample& sample, const MT2Config& cfg, float lumi=1. );
MT2Analysis<MT2EstimateTree>* mergeYields( std::vector< MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max=-1, const std::string& legendName="" );

int main(int argc, char* argv[]) {


  // std::string regionsSet = "zurich";
  std::string regionsSet = "13TeV_inclusive";
 if( argc>2 ) {
    regionsSet = std::string(argv[2]);
  }

  std::cout << "-> Using regions: " << regionsSet << std::endl;

  
  if( argc<2 ) {
    std::cout << "USAGE: ./regionEventYields [configFileName] regionSet" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  

  std::string configFileName(argv[1]);


  MT2Config cfg("cfgs/" + configFileName + ".txt");


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
  //  std::string samplesFileName = "../samples/samples_samples_PHYS14_skimprune.dat";

  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;



  std::string outputdir = "Zll_" + cfg.mcSamples() + "_" + regionsSet;
  //  std::string outputdir = "Zll_" + configFileName;
    double intpart;
    double fracpart = modf(lumi, &intpart);
    std::string suffix;
    if( fracpart>0. )
      suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
    else
      suffix = std::string( Form("_%.0ffb", intpart ) );
    outputdir += suffix;
  
  system(Form("mkdir -p %s", outputdir.c_str()));
 



  //Zll
  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, "DYJetsToLL", 700, 799  ); // not interested in signal here

  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }


  std::vector< MT2Analysis<MT2EstimateTree>* > EventYield;
  for( unsigned i=0; i<fSamples.size(); ++i ) 
    EventYield.push_back( computeYield( fSamples[i], cfg, lumi ) );
   

  MT2Analysis<MT2EstimateTree>* EventYield_zll = mergeYields( EventYield, cfg.regionsSet(), "DYJets", 700, 799, "DYJets" );



  EventYield_zll->writeToFile(outputdir+"/Zll_analyses.root","UPDATE");


  return 0;

}













MT2Analysis<MT2EstimateTree>* computeYield( const MT2Sample& sample, const MT2Config& cfg, float lumi ) {


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
  MT2Analysis<MT2EstimateTree>* analysis = new MT2Analysis<MT2EstimateTree>( sample.sname, regionsSet, sample.id );


  MT2EstimateTree::addVar( analysis, "Z_pt" );
  MT2EstimateTree::addVar( analysis, "Z_phi" );
  MT2EstimateTree::addVar( analysis, "Z_mass" );
  MT2EstimateTree::addVar( analysis, "Z_lepId" );
  MT2EstimateTree::addVar( analysis, "nLep" );

  


  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "   Entry: " << iEntry << " / " << nentries << std::endl;
    myTree.GetEntry(iEntry);

     if( !(myTree.passSelection("zll"))) continue; 
    /*
    if(!( myTree.nVert > 0 )) continue;
    //  if(!( myTree.nJet40 >= 2 )) continue;
    if(!( myTree.zll_deltaPhiMin > 0.3 )) continue;
    if(!( myTree.zll_diffMetMht < 0.5*myTree.zll_met_pt )) continue;
       //  if(!( myTree.jet1_pt>40. )) continue;
       //   if(!( myTree.jet2_pt>40. )) continue;
    if(!( myTree.nlep>1 )) continue; 

    */

     if(!( myTree.nlep==2 )) continue; 



    //Sample  are the Z leptons
    //and thus that if they don't have the same flavor they are rejected
    if( !(myTree.lep_pdgId[0] == -myTree.lep_pdgId[1]) ) continue;

    //Need the lorentz vectors of the leptons first
    TLorentzVector *LVec = new TLorentzVector[50];
    for(int i=0; i< 2; i++){
      LVec[i].SetPtEtaPhiM(myTree.lep_pt[i], myTree.lep_eta[i],myTree.lep_phi[i], myTree.lep_mass[i]);
    }



    double Z_invM_true = 91.19;
    TLorentzVector z = LVec[0] + LVec[1]; //leptons invariant mass
    double M_ll = z.M(); //Z mass

    //  if( abs(M_ll - Z_invM_true)>20.) continue;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = myTree.mt2;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet40;
    int nbjets = myTree.nBJet40;

    Double_t weight = myTree.evt_scale1fb*lumi;

    MT2EstimateTree* thisEstimate = analysis->get( myTree.zll_ht, njets, nbjets, myTree.zll_met_pt, minMTBmet, myTree.zll_mt2 );
    if( thisEstimate==0 ) continue; 

    //Fills the variables defined in MT2EstimateTree to the tree
    //at leatst partially...
    thisEstimate->fillTree(myTree, weight ,"zll");

    //initialize
    thisEstimate->assignVar("Z_pt", z.Perp() );
    thisEstimate->assignVar("Z_phi", z.Phi() );
    thisEstimate->assignVar("Z_mass", z.M() );
    thisEstimate->assignVar("Z_lepId", abs(myTree.lep_pdgId[0])  );



    //Fills the above variables into the tree
    thisEstimate->tree->Fill(); 


    thisEstimate->yield->Fill(myTree.zll_mt2, weight );
  
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

    if( EventYield[i]->id >= id_min && EventYield[i]->id <= id_max ) {

       *(return_EventYield) += *(EventYield[i]);

     }

  } // for EventYield


  return return_EventYield;

}































//Configuration file stuff
MT2Config::MT2Config( const std::string& configFileName ) {

  std::cout << std::endl;
  std::cout << "-> Reading config file: " << configFileName << std::endl;
  std::cout << std::endl;

  regionsSet_ = "13TeV_inclusive"; 
  //  regionsSet_ = ""; 
  mcSamples_ = "";
  sigSamples_ = "";
  dataSamples_ = "";
  lostLeptonTag_ = "";
  qcdTag_ = "";
  zinvTag_ = "";
  additionalStuff_ = "";

  ifstream IN(configFileName.c_str());
  char buffer[200];
  char StringValue[1000];


  while( IN.getline(buffer, 200, '\n') ) {

    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'                        
    }

    std::cout << buffer << std::endl;

    char name_c[200];
    sscanf(buffer, "%s %s", name_c, StringValue);
    std::string name(name_c);

    if( name=="regionsSet" )
      regionsSet_ = std::string(StringValue);
    else if( name=="mcSamples" )
      mcSamples_ = std::string(StringValue);
    else if( name=="sigSamples" )
      sigSamples_ = std::string(StringValue);
    else if( name=="dataSamples" )
      dataSamples_ = std::string(StringValue);
    else if( name=="lostLeptonTag" )
      lostLeptonTag_ = std::string(StringValue);
    else if( name=="qcdTag" )
      qcdTag_ = std::string(StringValue);
    else if( name=="zinvTag" )
      zinvTag_ = std::string(StringValue);
    else if( name=="additionalStuff" )
      additionalStuff_ = std::string(StringValue);

  } // while getline

  if( mcSamples_=="" && lostLeptonTag_=="" && qcdTag_=="" && zinvTag_=="" ) {
    std::cout << "[MT2Config] ERROR! Config file missing BG estimates!" << std::endl;
    exit(333);
  }

  if( mcSamples_!="" && ( lostLeptonTag_!="" || qcdTag_!="" || zinvTag_!="" ) ) {
    std::cout << "[MT2Config] ERROR! Config file must have either a mcSamples line OR the lostLeptonTag/qcdTag/zinvTag lines. Not both!" << std::endl;
    exit(335);
  }

  if( mcSamples_=="" && !( lostLeptonTag_!="" || qcdTag_!="" || zinvTag_!="" ) ) {
    std::cout << "[MT2Config] ERROR! All three data-driven BG estimate tags need to be specified in the config (lostLeptonTag/qcdTag/zinvTag)!" << std::endl;
    exit(337);
  }

  if( mcSamples_!="" && sigSamples_!="" ) {
    std::cout << "[MT2Config] ERROR! Config file must have either a mcSamples line OR (exclusive OR) a sigSamples line together with BG estimate tags." << std::endl;
    exit(339);
  }

  std::cout << std::endl;
     
}


