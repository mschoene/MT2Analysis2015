#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateQCD.h"

#define mt2_cxx
#include "../interface/mt2.h"

#include "TEventList.h"

float lumi = 4.; //fb-1


int main( int argc, char* argv[] ) {

  std::string regionsSet = "zurich";

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string samplesFileName = "PHYS14_qcdCR";
  std::string CRdir( Form("QCDcontrolRegion_%s_%s_%.0ffb", samplesFileName.c_str(), regionsSet.c_str(), lumi) );


  MT2Analysis<MT2EstimateQCD>* qcdCR   = MT2Analysis<MT2EstimateQCD>::readFromFile(CRdir + "/mc.root", "qcdCR");
  MT2Analysis<MT2EstimateQCD>* topCR   = MT2Analysis<MT2EstimateQCD>::readFromFile(CRdir + "/mc.root", "topCR");
  MT2Analysis<MT2EstimateQCD>* wjetsCR = MT2Analysis<MT2EstimateQCD>::readFromFile(CRdir + "/mc.root", "wjetsCR");
  MT2Analysis<MT2EstimateQCD>* zjetsCR = MT2Analysis<MT2EstimateQCD>::readFromFile(CRdir + "/mc.root", "zjetsCR");

  qcdCR->writeToFile( CRdir + "/pseudodata.root" );

  MT2Analysis<MT2EstimateQCD>* allCR = new MT2Analysis<MT2EstimateQCD>("MCallCR", regionsSet);
  MT2Analysis<MT2EstimateQCD>* allCRps1200 = new MT2Analysis<MT2EstimateQCD>("MCallCRps1200", regionsSet);
  MT2Analysis<MT2EstimateQCD>* allCRps125  = new MT2Analysis<MT2EstimateQCD>("MCallCRps125" , regionsSet);

  *allCR = *qcdCR + *topCR + *wjetsCR + *zjetsCR;
  *allCRps1200 = *allCR;
  *allCRps125  = *allCR;
  allCR->finalize();
  allCR->addToFile( CRdir + "/pseudodata.root" );

  allCR->setName("dataCR");
  allCR->randomizePoisson();
  allCR->finalize();
  allCR->addToFile( CRdir + "/pseudodata.root" );

  allCR->setName("dataSubCR");
  *allCR = *allCR - (*topCR + *wjetsCR + *zjetsCR);
  allCR->finalize();
  allCR->addToFile( CRdir + "/pseudodata.root" );

  allCRps1200->setName("dataCRps1200");
  allCRps1200->randomizePoisson(1./1200);
  allCRps1200->finalize();
  allCRps1200->addToFile( CRdir + "/pseudodata.root" );

  allCRps1200->setName("dataSubCRps1200");
  *allCRps1200 = *allCRps1200 - (1./1200)*(*topCR + *wjetsCR + *zjetsCR);
  allCRps1200->finalize();
  allCRps1200->addToFile( CRdir + "/pseudodata.root" );

  allCRps125->setName("dataCRps125");
  allCRps125->randomizePoisson(1./125);
  allCRps125->finalize();
  allCRps125->addToFile( CRdir + "/pseudodata.root" );

  allCRps125->setName("dataSubCRps125");
  *allCRps125 = *allCRps125 - (1./125)*(*topCR + *wjetsCR + *zjetsCR);
  allCRps125->finalize();
  allCRps125->addToFile( CRdir + "/pseudodata.root" );

}


