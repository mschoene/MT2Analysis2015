#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateQCD.h"

#define mt2_cxx
#include "../interface/mt2.h"

#include "TEventList.h"

float lumi = 0.454; //fb-1

bool doDiffMetMht = true;
bool doHFjetVeto = false;

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

  //std::string samplesFileName = "74X_jecV4_MET30_QCD";
  std::string samplesFileName = "Spring15_25ns_qcdSkim";
  //std::string regionsSet = "zurich_HTtriggers2";
  std::string regionsSet = "zurich_onlyHT";
  //std::string regionsSet = "zurich_HTtriggers";
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
  MT2Analysis<MT2EstimateQCD>* allCRps180 = new MT2Analysis<MT2EstimateQCD>("MCallCRps180", regionsSet);
  MT2Analysis<MT2EstimateQCD>* allCRps60  = new MT2Analysis<MT2EstimateQCD>("MCallCRps60" , regionsSet);
  MT2Analysis<MT2EstimateQCD>* dataCRps180 = new MT2Analysis<MT2EstimateQCD>("dataSubCRps180", regionsSet);
  MT2Analysis<MT2EstimateQCD>* dataCRps60  = new MT2Analysis<MT2EstimateQCD>("dataSubCRps60" , regionsSet);

  *allCR = *qcdCR + *topCR + *wjetsCR + *zjetsCR;
  *allCRps180 = *allCR;
  *allCRps60  = *allCR;
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
  *dataCRps180 = *dataCR; dataCRps180->setName("dataSubCRps180");
  *dataCRps60  = *dataCR; dataCRps60 ->setName("dataSubCRps60");
  *dataCR = *dataCR - (*topCR + *wjetsCR + *zjetsCR);
  *dataCRps180 = *dataCRps180 - (1./180)*(*topCR + *wjetsCR + *zjetsCR);
  *dataCRps60  = *dataCRps60  - (1./60)*(*topCR + *wjetsCR + *zjetsCR);
  dataCR->finalize();
  dataCRps180->finalize();
  dataCRps60->finalize();
  dataCR->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );
  dataCRps180->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );
  dataCRps60->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCRps180->setName("pdataCRps180");
  if (doRand)
      allCRps180->randomizePoisson(1./180);
  else
      allCRps180->sqrtErrors(1./180);
  allCRps180->finalize();
  allCRps180->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCRps180->setName("pdataSubCRps180");
  *allCRps180 = *allCRps180 - (1./180)*(*topCR + *wjetsCR + *zjetsCR);
  allCRps180->finalize();
  allCRps180->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCRps60->setName("pdataCRps60");
  if (doRand)
      allCRps60->randomizePoisson(1./60);
  else
      allCRps60->sqrtErrors(1./60);
  allCRps60->finalize();
  allCRps60->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

  allCRps60->setName("pdataSubCRps60");
  *allCRps60 = *allCRps60 - (1./60)*(*topCR + *wjetsCR + *zjetsCR);
  allCRps60->finalize();
  allCRps60->addToFile( CRdir + "/qcdCRmod"+rand_text+".root" );

}


