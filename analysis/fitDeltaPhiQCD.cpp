#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TProfile.h"


#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateQCD.h"
#include "../interface/MT2EstimateTree.h"







void fillFromTreeAndRatio( MT2Estimate* estimate, MT2Estimate* r_effective, TTree* tree, TF1* f1_ratio );



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




  // always start from inclusive qcd tree:
  MT2Analysis<MT2EstimateTree>* qcdTree_mc = MT2Analysis<MT2EstimateTree>::readFromFile( "EventYields_data_Run2015D_25nsGolden_v4/qcdControlRegion/mc.root", "qcdCRtree" );
  //MT2Analysis<MT2EstimateTree>* qcdTree_mc = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/mc.root", "qcdCRtree" );


  MT2Analysis<MT2Estimate>* estimate     = new MT2Analysis<MT2Estimate>("qcdEstimate", cfg.regionsSet());
  MT2Analysis<MT2Estimate>* r_effective  = new MT2Analysis<MT2Estimate>("r_effective", cfg.regionsSet());
  //MT2Analysis<MT2Estimate>* r_hat        = new MT2Analysis<MT2Estimate>("r_hat"      , cfg.regionsSet());
  //MT2Analysis<MT2Estimate>* f_jets       = new MT2Analysis<MT2Estimate>("f_jets"     , cfg.regionsSet());


  TH1D::AddDirectory(kTRUE);

  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion/";
  std::string outputdir = qcdCRdir;

  MT2Analysis<MT2EstimateQCD>* mc_all  = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "mc"     , cfg.qcdRegionsSet(), qcdTree_mc, "" );
  MT2Analysis<MT2EstimateQCD>* qcdOnly = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "qcdOnly", cfg.qcdRegionsSet(), qcdTree_mc, "id>=100 && id<200" );



  std::set<MT2Region> regions = estimate->getRegions();


  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* this_estimate     = estimate    ->get( *iR );
    MT2Estimate* this_r_effective  = r_effective ->get( *iR );
    //MT2Estimate* this_r_hat        = r_hat       ->get( *iR );
    //MT2Estimate* this_f_jets       = f_jets      ->get( *iR );

    MT2Region* matchedRegion = qcdOnly->matchRegion(*iR);
    MT2EstimateQCD* matchedQCDEstimate = qcdOnly->get( *matchedRegion );
    TF1* f1_ratio = matchedQCDEstimate->getFit( "pow", 70., 100. );

    fillFromTreeAndRatio( this_estimate, this_r_effective, matchedQCDEstimate->tree, f1_ratio );

  }  // for regions
      
    

  mc_all ->writeToFile( outputdir + "/mcFits.root", "recreate" );
  qcdOnly->writeToFile( outputdir + "/mcFits.root" );


  estimate ->writeToFile( outputdir + "/qcdEstimate.root", "recreate" );
  r_effective ->writeToFile( outputdir + "/qcdEstimate.root" );
  //r_hat ->writeToFile( outputdir + "/qcdEstimate.root" );
  //f_jets ->writeToFile( outputdir + "/qcdEstimate.root" );

  return 0;

}




void fillFromTreeAndRatio( MT2Estimate* estimate, MT2Estimate* r_effective, TTree* tree, TF1* f1_ratio ) {

  int nBins = estimate->yield->GetNbinsX();
  float xMin = estimate->yield->GetXaxis()->GetXmin();
  float xMax = estimate->yield->GetXaxis()->GetXmax();

  TProfile* rprof = new TProfile( "rprof", "", nBins, xMin, xMax );

  float weight;
  tree->SetBranchAddress( "weight", &weight );
  float mt2;
  tree->SetBranchAddress( "mt2", &mt2 );

  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    float r = f1_ratio->Eval( mt2 );
    estimate->yield->Fill( mt2, weight*r );

    rprof->Fill( mt2, r, weight );

  } // for entries

  for( int iBin=1; iBin<nBins+1; ++iBin ) {
    r_effective->yield->SetBinContent( iBin, rprof->GetBinContent(iBin) );
    r_effective->yield->SetBinError( iBin, rprof->GetBinError(iBin) );
  }

  delete rprof;

}
