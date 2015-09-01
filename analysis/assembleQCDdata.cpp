#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateQCD.h"

#define mt2_cxx
#include "../interface/mt2.h"

#include "TEventList.h"

float lumi = 0.042; //fb-1

bool doDiffMetMht = true;
bool doHFjetVeto = true;

bool doRand = false; // randomize pseudoData


// the following is done in this analysis:
// - simulate data sample from MC
// - subtract non-QCD from (pseudo)data
// - apply prescales to (pseudo)data from prescaled triggers
// - redo the fits after steps above

int main( int argc, char* argv[] ) {


  if( argc>4 ) {
    std::cout << "USAGE: ./qcdControlRegion [samplesFileName] [regionsSet] [postfix]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  std::string samplesFileName = "74X_jecV4_MET30_QCD";
  //std::string regionsSet = "zurich_HTtriggers2";
  std::string regionsSet = "zurich_HTtriggers";
  std::string postfix = "";

  if( argc>1 ) {
    std::string samplesFileName_tmp(argv[1]); 
    samplesFileName = samplesFileName_tmp;
    if( argc>2 ) {
      std::string regionsSet_tmp(argv[2]); 
      regionsSet = regionsSet_tmp; 
      if( argc>3 ) {
	std::string postfix_tmp(argv[3]); 
	postfix = postfix_tmp; 
      }
    }
  }


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string CRdir( Form("QCDcontrolRegion_%s_%s_%.3ffb%s%s%s", samplesFileName.c_str(), regionsSet.c_str(), lumi,(doDiffMetMht ? "" : "_noDiffMetMht"), (!doHFjetVeto ? "" : "_HFjetVeto"), postfix.c_str()) );


  MT2Analysis<MT2EstimateQCD>* qcdCR   = MT2Analysis<MT2EstimateQCD>::readFromFile(CRdir + "/qcdCR.root", "qcdCR"  );
  MT2Analysis<MT2EstimateQCD>* topCR   = MT2Analysis<MT2EstimateQCD>::readFromFile(CRdir + "/qcdCR.root", "topCR"  );
  MT2Analysis<MT2EstimateQCD>* wjetsCR = MT2Analysis<MT2EstimateQCD>::readFromFile(CRdir + "/qcdCR.root", "wjetsCR");
  MT2Analysis<MT2EstimateQCD>* zjetsCR = MT2Analysis<MT2EstimateQCD>::readFromFile(CRdir + "/qcdCR.root", "zjetsCR");
  MT2Analysis<MT2EstimateQCD>* dataCR  = MT2Analysis<MT2EstimateQCD>::readFromFile(CRdir + "/qcdCR.root", "dataCR" );

  std::string rand_text = doRand ? "_rand" : "";

  qcdCR->writeToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  MT2Analysis<MT2EstimateQCD>* allCR = new MT2Analysis<MT2EstimateQCD>("MCallCR", regionsSet);
  MT2Analysis<MT2EstimateQCD>* allCRps700 = new MT2Analysis<MT2EstimateQCD>("MCallCRps700", regionsSet);
  MT2Analysis<MT2EstimateQCD>* allCRps175  = new MT2Analysis<MT2EstimateQCD>("MCallCRps175" , regionsSet);

  *allCR = *qcdCR + *topCR + *wjetsCR + *zjetsCR;
  *allCRps700 = *allCR;
  *allCRps175  = *allCR;
  allCR->finalize();
  allCR->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCR->setName("pdataCR");
  if (doRand)
      allCR->randomizePoisson();
  else
      allCR->sqrtErrors();
  allCR->finalize();
  allCR->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  dataCR->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCR->setName("pdataSubCR");
  *allCR = *allCR - (*topCR + *wjetsCR + *zjetsCR);
  allCR->finalize();
  allCR->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  dataCR->setName("dataSubCR");
  *dataCR = *dataCR - (*topCR + *wjetsCR + *zjetsCR);
  dataCR->finalize();
  dataCR->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCRps700->setName("pdataCRps700");
  if (doRand)
      allCRps700->randomizePoisson(1./700);
  else
      allCRps700->sqrtErrors(1./700);
  allCRps700->finalize();
  allCRps700->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCRps700->setName("pdataSubCRps700");
  *allCRps700 = *allCRps700 - (1./700)*(*topCR + *wjetsCR + *zjetsCR);
  allCRps700->finalize();
  allCRps700->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCRps175->setName("pdataCRps175");
  if (doRand)
      allCRps175->randomizePoisson(1./175);
  else
      allCRps175->sqrtErrors(1./175);
  allCRps175->finalize();
  allCRps175->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCRps175->setName("pdataSubCRps175");
  *allCRps175 = *allCRps175 - (1./175)*(*topCR + *wjetsCR + *zjetsCR);
  allCRps175->finalize();
  allCRps175->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

}


