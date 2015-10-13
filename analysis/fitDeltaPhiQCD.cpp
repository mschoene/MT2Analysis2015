#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"


#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateQCD.h"
#include "../interface/MT2EstimateTree.h"










int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|              Running fitDeltaPhiQCD                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./fitDeltaPhiQCD [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool useMC = true;

  if( argc>2 ) {

    std::string mc_or_data_templates = std::string(argv[2]); 
    if( mc_or_data_templates=="mc" ) mc_or_data_templates="MC";
    std::cout << std::endl;
    std::cout << "-> Will disobey the cfg and use mc_or_data_templates = " << mc_or_data_templates << std::endl;
    cfg.set_gammaTemplateType(mc_or_data_templates);
    if( mc_or_data_templates=="MC" ) useMC = true;
    else useMC=false;
    std::cout << std::endl;

  } 



  //MT2DrawTools::setStyle();


  //TH1::AddDirectory(kTRUE);


  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion"; 
  MT2Analysis<MT2EstimateTree>* qcdTree_mc = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/mc.root", "qcdCRtree" );


  std::string outputdir = qcdCRdir;
  //std::string outputdir = qcdCRdir + "/QCDFits";
  //system( Form( "mkdir -p %s", outputdir.c_str()) );

  std::set<MT2Region> regions = qcdTree_mc->getRegions();

  //MT2Analysis<MT2EstimateQCD>* theFits     = MT2EstimateQCD::makeAnalysisFromEstimateTree( "mc"     , qcdTree_mc, "" );
  MT2Analysis<MT2EstimateQCD>* theFits     = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "mc"     , cfg.qcdRegionsSet(), qcdTree_mc, "" );
  theFits->finalize();

  //MT2Analysis<MT2EstimateQCD>* theFits_qcd = MT2EstimateQCD::makeAnalysisFromEstimateTree( "qcdOnly", qcdTree_mc, "id>=100 && id<200" );
  MT2Analysis<MT2EstimateQCD>* theFits_qcd = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "qcdOnly", cfg.qcdRegionsSet(), qcdTree_mc, "id>=100 && id<200" );
  theFits_qcd->finalize();

  theFits_qcd->writeToFile( outputdir + "/mcFits.root", "recreate" );
  theFits->writeToFile( outputdir + "/mcFits.root" );

  return 0;

}





