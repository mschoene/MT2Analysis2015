#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateZinvGamma.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"
#include "TRandom3.h"


float lumi = 4.; //fb-1
bool dummyData = true;



MT2Analysis<MT2EstimateZinvGamma> computeYield( const MT2Sample& sample, const std::string& regionsSet, bool onlyPrompt );
void randomizePoisson( MT2Analysis<MT2EstimateZinvGamma>* data );
void randomizeSingleHisto( TRandom3* rand, TH1D* h1 );



int main( int argc, char* argv[] ) {


  std::string regionsSet = "13TeV_inclusive";
  //if( argc>1 ) {
  //  std::string regionsSet_tmp(argv[1]); 
  //  regionsSet = regionsSet_tmp;
  //}

  std::string samplesFileName = "PHYS14_v2_Zinv";

  std::string samplesFile = "../samples/samples_" + samplesFileName + ".dat";
  
  std::vector<MT2Sample> samples = MT2Sample::loadSamples(samplesFile, 100, 299); // QCD + GJet
  if( samples.size()==0 ) {
    std::cout << "There must be an error: didn't find any files in " << samplesFile << "!" << std::endl;
    exit(1209);
  }


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string outputdir = "ZinvGammaPurityFromData_" + samplesFileName + "_" + regionsSet;
  system(Form("mkdir -p %s", outputdir.c_str()));


  
  MT2Analysis<MT2EstimateZinvGamma>* templatesPromptRaw  = new MT2Analysis<MT2EstimateZinvGamma>( "templatesPromptRaw", regionsSet );
  MT2Analysis<MT2EstimateZinvGamma>* templatesFake       = new MT2Analysis<MT2EstimateZinvGamma>( "templatesFake", regionsSet );


  for( unsigned i=0; i<samples.size(); ++i ) {
    (*templatesPromptRaw) += (computeYield( samples[i], regionsSet, true  ));
    (*templatesFake) += (computeYield( samples[i], regionsSet, false  ));
  }

  
  randomizePoisson( templatesPromptRaw );
  randomizePoisson( templatesFake );

  MT2Analysis<MT2EstimateZinvGamma>* templatesPrompt  = new MT2Analysis<MT2EstimateZinvGamma>( "templatesPrompt", regionsSet );
  (*templatesPrompt) = (*templatesPromptRaw) - (*templatesFake);


  std::string templateFileName = "gammaTemplatesDummy_" + samplesFileName + "_" + regionsSet + ".root";

  templatesPrompt->writeToFile(templateFileName);
  templatesFake->addToFile(templateFileName);

  return 0;

}






MT2Analysis<MT2EstimateZinvGamma> computeYield( const MT2Sample& sample, const std::string& regionsSet, bool onlyPrompt ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;



  MT2Analysis<MT2EstimateZinvGamma> analysis( sample.sname, regionsSet, sample.id );

  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  int nentries = tree->GetEntries();




  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    // template control region (use low-significance regions):
    if( myTree.gamma_ht > 1000. ) continue;
    if( myTree.gamma_mt2 < 200.) continue;
    if( myTree.gamma_mt2 > 300. ) continue;
    if( myTree.met_pt > 100.) continue;
    //if( myTree.gamma_nJet40 > 4) continue;
    if( myTree.gamma_nBJet40 > 0) continue;
    //if( myTree.gamma_ht>1000. && sample.id==204 ) continue; // remove high-weight spikes (remove GJet_400to600 leaking into HT>1000)



    if( myTree.nMuons10 > 0) continue;
    if( myTree.nElectrons10 > 0 ) continue;
    if( myTree.nPFLep5LowMT > 0) continue;
    if( myTree.nPFHad10LowMT > 0) continue;

    if( myTree.gamma_deltaPhiMin<0.3 ) continue;
    if( myTree.gamma_diffMetMht>0.5*myTree.gamma_met_pt ) continue;
  
    if( myTree.nVert==0 ) continue;

    if( myTree.gamma_nJet40<2 ) continue;

    if( myTree.ngamma==0 ) continue;



    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );
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


    float sietaieta = myTree.gamma_sigmaIetaIeta[0];
    bool isPrompt = false;
    if( fabs( gamma.Eta() )<1.4445 ) {
      isPrompt = sietaieta<0.01;
    } else {
      isPrompt = sietaieta<0.027;
    }
    
    if(  onlyPrompt && !isPrompt ) continue;
    if( !onlyPrompt &&  isPrompt ) continue;


    Double_t weight = myTree.evt_scale1fb*lumi; 

    MT2EstimateZinvGamma* thisEstimate = analysis.get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
    if( thisEstimate==0 ) continue;

    thisEstimate->yield->Fill(myTree.gamma_mt2, weight );
    thisEstimate->sietaieta->Fill(sietaieta, weight );

    float iso = myTree.gamma_chHadIso[0]/myTree.gamma_pt[0];

    thisEstimate->fillIso( iso, weight, myTree.gamma_mt2 );

    
  } // for entries


  analysis.finalize();
  

  delete tree;


  file->Close();
  delete file;
  
  return analysis;

}





void randomizePoisson( MT2Analysis<MT2EstimateZinvGamma>* data ) {

  TRandom3 rand(13);


  std::set<MT2Region> MT2Regions = data->getRegions();


  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
      
    randomizeSingleHisto( &rand, data->get(thisRegion)->yield );
    randomizeSingleHisto( &rand, data->get(thisRegion)->iso );
    randomizeSingleHisto( &rand, data->get(thisRegion)->sietaieta );

  } // for regions

}




void randomizeSingleHisto( TRandom3* rand, TH1D* h1 ) {

  for( unsigned ibin=1; ibin<h1->GetXaxis()->GetNbins()+1; ++ibin ) {

    int poisson_data = rand->Poisson(h1->GetBinContent(ibin));
    h1->SetBinContent(ibin, poisson_data);
    h1->SetBinError(ibin, 0.);
    
  }  // for bins

}
