#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>



#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2DrawTools.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2EstimateSyst.h"

#define mt2_cxx
#include "../interface/mt2.h"


#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLorentzVector.h"



bool do_dummyMC = false;
bool use_extrapolation = false;
bool do_hybrid = true;


MT2Analysis<MT2Estimate>* getInclusiveRatioMC( const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* Zinv, MT2Analysis<MT2EstimateTree>* gammaCRtree );
MT2Analysis<MT2EstimateSyst>* combineDataAndMC( MT2Analysis<MT2EstimateSyst>* data, MT2Analysis<MT2Estimate>* mc );

MT2Analysis<MT2EstimateSyst>* computePurityOF( MT2Analysis<MT2Estimate>* SF, MT2Analysis<MT2Estimate>* OF, bool do_Runcert=0 );


void extrapolToTopoRegion( MT2Analysis<MT2Estimate>* shape_TR, MT2Analysis<MT2Estimate>* shape, bool isMC=0 );
//void buildHybrid( MT2Analysis<MT2Estimate>* shape_hybrid, MT2Analysis<MT2Estimate>* shape_data, MT2Analysis<MT2Estimate>* shape_MC, MT2Analysis<MT2Estimate>* MC_ratio, MT2Analysis<MT2Estimate>* bin_extrapol );
void buildHybrid( MT2Analysis<MT2Estimate>* shape_hybrid, MT2Analysis<MT2Estimate>* shape_data, MT2Analysis<MT2Estimate>* shape_MCsr, MT2Analysis<MT2Estimate>* shape_MCcr, MT2Analysis<MT2Estimate>* bin_extrapol );

int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeZinvFromZll               |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc!=2 ) {
    std::cout << "USAGE: ./computeZinvFromZll [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  float lumi = cfg.lumi();

  std::string zllControlRegionDir = cfg.getEventYieldDir() + "/zllControlRegion"; //(Form("ZllControlRegion_%s_%s_%.0ffb", samples.c_str(), regionsSet.c_str(), lumi));




  MT2Analysis<MT2Estimate>* TopMC_ = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_Top_forZinvEst.root", "zllCR");
  (*TopMC_) = (*TopMC_) * lumi;

  MT2Analysis<MT2Estimate>* TopMC_of = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_Top_of_forZinvEst.root", "zllCR_of");
  (*TopMC_of) = (*TopMC_of) * lumi;


  MT2Analysis<MT2Estimate>* TopMC2_ = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_Top_forZinvEst.root", "zllCR");
  (*TopMC2_) = (*TopMC2_) * lumi;


  // MT2Analysis<MT2Estimate>* signal = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/signal_forZinvEst.root", "zllSigCR");
  // (*signal) = (*signal) * lumi;

 
 MT2Analysis<MT2Estimate>* zllData;
  if( !do_dummyMC )
    zllData = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/data_forZinvEst.root", "data");
  else{
    zllData = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_forZinvEst.root", "zllCR");
    (*zllData) = (*zllData) * lumi;
    (*zllData) = (*zllData) +  (*TopMC_);
  }

 
  MT2Analysis<MT2Estimate>* zllData_of;
  if( !do_dummyMC )
    zllData_of = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/data_of_forZinvEst.root", "data_of");
  else{
    zllData_of = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_of_forZinvEst.root", "zllCR_of");
    (*zllData_of) = (*zllData_of) * lumi;
    (*zllData_of) = (*zllData_of) +  (*TopMC_of);
  }
 

  MT2Analysis<MT2EstimateTree>* TopTree_mc = MT2Analysis<MT2EstimateTree>::readFromFile(zllControlRegionDir + "/ZllPurityTrees_of.root", "Top");

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this
  MT2Analysis<MT2EstimateTree>* TTZ = MT2EstimateTree::makeAnalysisFromInclusiveTree( "TTZ" , cfg.regionsSet(), TopTree_mc,  "(id==413 && Z_pt>200.)");
  (*TTZ) *= lumi;
  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


   MT2Analysis<MT2Estimate>* TopMC_woTTZ = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_Top_forZinvEst.root", "zllCR");
  (*TopMC_woTTZ) = (*TopMC_woTTZ) * lumi;
  (*TopMC_woTTZ) = (*TopMC_woTTZ) - (*TTZ);


   MT2Analysis<MT2Estimate>* MC_of = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_of_forZinvEst.root", "zllCR_of");
  (*MC_of) = (*MC_of) * lumi;
  (*MC_of) = (*MC_of) +  (*TopMC_of);

  
  // //static MT2Analysis<MT2EstimateTree>* makeAnalysisFromInclusiveTree( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateTree>* analysis, const std::string& selection="" ) { return makeRebinnedAnalysisFromInclusiveTree( aname, regionsSet, analysis, selection ); };


  // MT2Analysis<MT2Estimate>* bg_tot;
  // (*bg_tot) = (*zllData_of) +  (*TTZ);

 
  MT2Analysis<MT2Estimate>* zllMC_ = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_forZinvEst.root", "zllCR");

  MT2Analysis<MT2Estimate>* zllMC2_ = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_forZinvEst.root", "zllCR");

  MT2Analysis<MT2Estimate>* zllMC_temp = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_forZinvEst.root", "zllCR");


 


  MT2Analysis<MT2Estimate>* zllMC;
  if ( !use_extrapolation )  {
    zllMC = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_forZinvEst.root", "");
    (*zllMC) = (*zllMC) * lumi;
    (*zllMC_) = (*zllMC_) * lumi; 
    (*zllMC2_) = (*zllMC2_) * lumi; 
  }else{
    (*zllMC_) = (*zllMC_) * lumi;
    (*zllMC2_) = (*zllMC2_) * lumi;
    (*zllMC_temp) = (*zllMC_temp) * lumi;
  }

  MT2Analysis<MT2Estimate>* zllMC_integral;
  if( use_extrapolation )
    zllMC_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "zllMC_integral", cfg.regionsSet(), zllMC_temp );




  if( zllData==0 || zllMC_==0  ) {
    std::cout << "-> Please run gammaControlRegion first. I need to get the zllData yields from there." << std::endl;
    std::cout << "->Thank you for your cooperation." << std::endl;
    exit(193);
  }

  
  MT2Analysis<MT2EstimateTree>* Zinv = MT2Analysis<MT2EstimateTree>::readFromFile(cfg.getEventYieldDir() + "/analyses.root", "ZJets");
  if( Zinv==0 ) {
    std::cout << "-> Please run regionEventYields on MC first. I need to get the Z->vv MC yields from there." << std::endl;
    std::cout << "-> Thank you for your cooperation." << std::endl;
    exit(197);
  }
  (* (MT2Analysis<MT2Estimate>*) Zinv) = (* (MT2Analysis<MT2Estimate>*)Zinv) * lumi;


  //MT2Analysis<MT2EstimateTree>* gammaCRtree = MT2Analysis<MT2EstimateTree>::readFromFile(gammaControlRegionDir + "/data.root", "gammaCRtree");
  //MT2Analysis<MT2Estimate>* ZgammaRatioMC = getInclusiveRatioMC( regionsSet, Zinv, gammaCRtree );

  MT2Analysis<MT2Estimate>* ZinvZllRatioMC = new MT2Analysis<MT2Estimate>( "ZinvZllRatioMC", cfg.regionsSet() );
  if( !use_extrapolation )
    (*ZinvZllRatioMC) = ( (* (MT2Analysis<MT2Estimate>*)Zinv) / (*zllMC) );
  else
    (*ZinvZllRatioMC) = ( (* (MT2Analysis<MT2Estimate>*)Zinv) / (*zllMC_integral) );
  //  (*ZinvZllRatioMC) = ( (* (MT2Analysis<MT2Estimate>*)Zinv) / (*zllMC) );

  //MT2Analysis<MT2Estimate>* ZinvZllRatio = MT2EstimateSyst::makeAnalysisFromEstimate( "ZinvZllRatio", regionsSet, ZinvZllRatioMC );
  MT2Analysis<MT2Estimate>* ZinvZllRatio = new MT2Analysis<MT2Estimate>( "ZinvZllRatio", cfg.regionsSet() );
  (*ZinvZllRatio) = (*ZinvZllRatioMC);



  MT2Analysis<MT2EstimateSyst>* zll_est_ = MT2EstimateSyst::makeAnalysisFromEstimate( "zll_est_", cfg.regionsSet(), zllData );

  MT2Analysis<MT2EstimateSyst>* zll_est;
  if( !use_extrapolation ) {
    zll_est = MT2EstimateSyst::makeAnalysisFromEstimate( "zll_est", cfg.regionsSet(), zllData );
  }



  //  MT2Analysis<MT2Estimate>* TopMC_integral;
  //  TopMC_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "TopMC_integral", cfg.regionsSet(), TopMC_ );


  //PurityStuff
  MT2Analysis<MT2Estimate>* diLep = new MT2Analysis<MT2Estimate>( "diLep", cfg.regionsSet() ); 
  (*diLep) = (*TopMC_) + (*zllMC_);

  MT2Analysis<MT2EstimateSyst>* purity = MT2EstimateSyst::makeEfficiencyAnalysis( "purity", (MT2Analysis<MT2Estimate>*)zllMC_, (MT2Analysis<MT2Estimate>*)diLep);

  (*zll_est_) = (*zll_est_ ) *  (*purity);

  //Purity Uncert comes from uncert on Top subtraction = 50% Top 
  MT2Analysis<MT2Estimate>* purity_err = new MT2Analysis<MT2Estimate>( "purity_err", cfg.regionsSet() ); 
  (*purity_err) = (*zllMC_ )*0.5*(*TopMC_)/ ((*diLep)*(*diLep)) ;
  

  std::cout << "Creating Shapes..." << std::endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Creating the shapes //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  MT2Analysis<MT2Estimate>* Zinv_forHybrid_notIntegral = MT2Analysis<MT2Estimate>::readFromFile(cfg.getEventYieldDir() + "/analyses.root", "ZJets");
  (* (MT2Analysis<MT2Estimate>*) Zinv_forHybrid_notIntegral) = (* (MT2Analysis<MT2Estimate>*)Zinv_forHybrid_notIntegral) * lumi;

  MT2Analysis<MT2Estimate>* zllMC_forHybrid_notIntegral = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_forZinvEst.root", "zllCR");
  (*zllMC_forHybrid_notIntegral) *= lumi;

  // MT2Analysis<MT2Estimate>* ZinvZllRatioHybrid_notIntegral = new MT2Analysis<MT2Estimate>( "ZinvZllRatioHybrid_notIntegral", cfg.regionsSet() );
  // (*ZinvZllRatioHybrid_notIntegral) = (( *Zinv_forHybrid_notIntegral) / (*zllMC_forHybrid_notIntegral) );


  MT2Analysis<MT2Estimate>* zllMC_forHybrid = MT2Estimate::makeIntegralAnalysisFromEstimate( "zllMC_forHybrid", cfg.regionsSet(), zllMC_forHybrid_notIntegral );

  MT2Analysis<MT2Estimate>* Zinv_forHybrid = MT2Estimate::makeIntegralAnalysisFromEstimate( "Zinv_forHybrid", cfg.regionsSet(), Zinv_forHybrid_notIntegral );
  //Integrated per topological region
  MT2Analysis<MT2Estimate>* ZinvZllRatioHybrid = new MT2Analysis<MT2Estimate>( "ZinvZllRatioHybrid", cfg.regionsSet() );
  (*ZinvZllRatioHybrid) = ( *Zinv_forHybrid) / (*zllMC_forHybrid);


  //Normalize the ratio
  //  (*ZinvZllRatioHybrid_notIntegral) /= (*ZinvZllRatioHybrid) ;



  MT2Analysis<MT2EstimateTree>* zllData_tree;
  if( !do_dummyMC )
    zllData_tree = MT2Analysis<MT2EstimateTree>::readFromFile(zllControlRegionDir + "/data.root", "data");
  else{
    zllData_tree = MT2Analysis<MT2EstimateTree>::readFromFile(zllControlRegionDir + "/mc.root", "zllCR");
  }

  MT2Analysis<MT2EstimateTree>* zllMC_tree = MT2Analysis<MT2EstimateTree>::readFromFile(zllControlRegionDir + "/mc.root", "zllCR");

  MT2Analysis<MT2EstimateTree>* zllData_tree_of = MT2Analysis<MT2EstimateTree>::readFromFile(zllControlRegionDir + "/data_of.root", "data_of");

  MT2Analysis<MT2EstimateTree>* zinvMC_tree = MT2Analysis<MT2EstimateTree>::readFromFile(cfg.getEventYieldDir() + "/ZJetsIncl.root", "ZJets_inclusive"); 
 
  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this
  TFile* temp = new TFile("temp.root","RECREATE");
  ///////////////////////////////////////////////
  ///////////FILLING THE SHAPE///////////////////
  ///////////////////////////////////////////////
  MT2Analysis<MT2EstimateTree>* zinvMC_forShape = MT2EstimateTree::makeAnalysisFromInclusiveTree( "zinv_forShape" , "zurich2016_forExtrapol", zinvMC_tree,  "");
  MT2Analysis<MT2EstimateTree>* zllMC_forShape = MT2EstimateTree::makeAnalysisFromInclusiveTree( "zllMC_forShape" , "zurich2016_forExtrapol", zllMC_tree,  "(fabs(Z_mass-91.19)<=20.) && Z_pt>200.");
   
  MT2Analysis<MT2EstimateTree>* zllData_shape = MT2EstimateTree::makeAnalysisFromInclusiveTree( "shape" , "zurich2016_forExtrapol", zllData_tree,  "(fabs(Z_mass-91.19)<=20.) && Z_pt>200.");
  
  MT2Analysis<MT2EstimateTree>* zllData_shape_of = MT2EstimateTree::makeAnalysisFromInclusiveTree( "shape_of" , "zurich2016_forExtrapol", zllData_tree_of,  "(fabs(Z_mass-91.19)<=20.) && Z_pt>200.");

  // MT2Analysis<MT2EstimateTree>* zllData_shape_noTTZ = MT2EstimateTree::makeAnalysisFromInclusiveTree( "shape_noTTZ" , "zurich2016_forExtrapol", zllData_tree,  "(fabs(Z_mass-91.19)<=20.) && Z_pt>200.");
  MT2Analysis<MT2EstimateTree>* TTZ_shape = MT2EstimateTree::makeAnalysisFromInclusiveTree( "TTZ_shape" , "zurich2016_forExtrapol", TopTree_mc,  "(id==413 && (fabs(Z_mass-91.19)<=20.) && Z_pt>200.)");

  MT2Analysis<MT2EstimateTree>* zllMC_shape = MT2EstimateTree::makeAnalysisFromInclusiveTree( "zllMC_shape" , "zurich2016_forExtrapol", zllMC_tree,  "(fabs(Z_mass-91.19)<=20.) && Z_pt>200.");

  MT2Analysis<MT2EstimateTree>* zllMC_shape_withTTZ = MT2EstimateTree::makeAnalysisFromInclusiveTree( "zllMC_shape_withTTZ" , "zurich2016_forExtrapol", zllMC_tree,  "(fabs(Z_mass-91.19)<=20.) && Z_pt>200.");

  // MT2Analysis<MT2EstimateTree>* zllMC_shape_TR = MT2EstimateTree::makeAnalysisFromInclusiveTree( "zllMC_shape_TR" , cfg.regionsSet(), zllMC_tree,  "(fabs(Z_mass-91.19)<=20.) && Z_pt>200.");

  std::cout << "got the cr in the regions wanted" << std::endl;
  //  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this
  ///////////////////////////////////////////////
  (*zinvMC_forShape) *= lumi; 
  (*zllMC_forShape) *= lumi; 

  (*TTZ_shape) *= lumi; 

  (*zllMC_shape_withTTZ) *= lumi;
  (*zllMC_shape_withTTZ) += (*TTZ_shape);

  (*zllMC_shape) = (*zllMC_shape) * lumi;
 
  if( do_dummyMC )
    (*zllData_shape) = (*zllData_shape) * lumi;

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  MT2Analysis<MT2Estimate>* bin_extrapol  = new MT2Analysis<MT2Estimate>( "bin_extrapol", cfg.regionsSet() );
 
  MT2Analysis<MT2Estimate>* zllData_shape_TR = new MT2Analysis<MT2Estimate>("zllData_shape_TR", cfg.regionsSet() ); 
  MT2Analysis<MT2Estimate>* zllData_of_shape_TR = new MT2Analysis<MT2Estimate>("zllData_of_shape_TR", cfg.regionsSet()); 
 
  MT2Analysis<MT2Estimate>* zllMC_shape_TR = new MT2Analysis<MT2Estimate>( "zllMC_shape_TR", cfg.regionsSet() );
  MT2Analysis<MT2Estimate>* zllHybrid_shape_TR = new MT2Analysis<MT2Estimate>( "zllHybrid_shape_TR", cfg.regionsSet() ); 

  MT2Analysis<MT2Estimate>* zinvMC_forShape_TR = new MT2Analysis<MT2Estimate>( "zinvMC_forShape_TR", cfg.regionsSet() );

  std::cout << "initialized the shape analyses" << std::endl;


  /////////DATA///////////////////
  MT2Analysis<MT2Estimate>* zllData_forHybrid_notIntegral;
  MT2Analysis<MT2Estimate>* zllData_of_forHybrid_notIntegral;
  if( !do_dummyMC ){
    zllData_forHybrid_notIntegral = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/data_forZinvEst.root", "data");
    //OF background:////////////////////////
    zllData_of_forHybrid_notIntegral = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/data_of_forZinvEst.root", "data_of");
  }else{
    zllData_forHybrid_notIntegral = MT2Analysis<MT2Estimate>::readFromFile(zllControlRegionDir + "/mc_forZinvEst.root", "zllCR");
    (*zllData_forHybrid_notIntegral) = (*zllData_forHybrid_notIntegral) * lumi;
  }

  
  MT2Analysis<MT2Estimate>* zllData_of_forHybrid = MT2Estimate::makeIntegralAnalysisFromEstimate( "zllData_of_forHybrid", cfg.regionsSet(), zllData_of_forHybrid_notIntegral );

  MT2Analysis<MT2Estimate>* zllData_forHybrid = MT2Estimate::makeIntegralAnalysisFromEstimate( "zllData_forHybrid", cfg.regionsSet(), zllData_forHybrid_notIntegral );

  MT2Analysis<MT2EstimateSyst>* purity_forHybrid = computePurityOF(zllData_forHybrid, zllData_of_forHybrid, 1);  

  std::cout << "filling the shapes from the extrapolation region to the TR ..." << std::endl;
  extrapolToTopoRegion( zllData_shape_TR, (MT2Analysis<MT2Estimate>*)zllData_shape );
  extrapolToTopoRegion( zllData_of_shape_TR, (MT2Analysis<MT2Estimate>*)zllData_shape_of );
  extrapolToTopoRegion( zllMC_shape_TR, (MT2Analysis<MT2Estimate>*)zllMC_shape, 1 ); //1 means it is mc
  extrapolToTopoRegion( zinvMC_forShape_TR, (MT2Analysis<MT2Estimate>*)zinvMC_forShape, 1 ); //1 means it is mc
 
  std::cout << "Filled the shapes into the TR" << std::endl;

  MT2Analysis<MT2EstimateSyst>* purity_shape_TR = computePurityOF(zllData_shape_TR, zllData_of_shape_TR );
  (*zllData_shape_TR) *= (*purity_shape_TR);
//  MT2Analysis<MT2Estimate>* zinvZllRatio_forShape_TR = new MT2Analysis<MT2Estimate>( "ZinvZllRatio_forShape_TR", cfg.regionsSet() );
//  (*zinvZllRatio_forShape_TR) = (*zinvMC_forShape_TR) / (*zllMC_shape_TR);

  /////////////////Building the hybrid shape from data & MC, and the ratio of the MC/////////////////////////
  buildHybrid( zllHybrid_shape_TR, zllData_shape_TR, zinvMC_forShape_TR, zllMC_shape_TR, bin_extrapol );
  //  buildHybrid( zllHybrid_shape_TR, zllData_shape_TR,zllMC_shape_TR ,zinvZllRatio_forShape_TR, bin_extrapol );

  //The ALPHA
  MT2Analysis<MT2Estimate>* alpha = new MT2Analysis<MT2Estimate>( "alpha", cfg.regionsSet() );
  (*alpha) =  (*zllHybrid_shape_TR) * (*ZinvZllRatioHybrid);
 
  MT2Analysis<MT2Estimate>* ZinvEstimateFromZll_hybrid = new MT2Analysis<MT2Estimate>( "ZinvEstimateFromZll_hybrid", cfg.regionsSet() );
  (*ZinvEstimateFromZll_hybrid) = (*zllData_forHybrid) * (*purity_forHybrid) * (*alpha) ;


  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////

  MT2Analysis<MT2EstimateSyst>* zll_est_integral;
  MT2Analysis<MT2Estimate>* diLep_dummy = new MT2Analysis<MT2Estimate>( "diLep_dummy", cfg.regionsSet() ); 
  (*diLep_dummy) = (*TopMC_) + (*zllMC_);
  MT2Analysis<MT2Estimate>* diLep_integral;
  MT2Analysis<MT2Estimate>* TopMC_integral;

  if( use_extrapolation ){
    zll_est_integral = MT2EstimateSyst::makeIntegralAnalysisFromEstimate("zll_est_integral", cfg.regionsSet(), zll_est_);
    TopMC_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "TopMC_integral", cfg.regionsSet(), TopMC_ );
    diLep_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "diLep_integral", cfg.regionsSet(), diLep_dummy ); }

  MT2Analysis<MT2EstimateSyst>* purity_integral;

  MT2Analysis<MT2Estimate>* purity_err_integral = new MT2Analysis<MT2Estimate>( "purity_err", cfg.regionsSet() );   if( use_extrapolation ){
    purity_integral = MT2EstimateSyst::makeEfficiencyAnalysis( "purity_integral", (MT2Analysis<MT2Estimate>*)zllMC_integral, (MT2Analysis<MT2Estimate>*)diLep_integral);

    (*purity_err_integral) = (*zllMC_integral )*0.5*(*TopMC_integral)/ ((*diLep_integral)*(*diLep_integral)) ;
  }
  //if( use_extrapolation )
  //   (*diLep_integral) = (*TopMC_integral) + (*zllMC_integral);
 



  // //SIGNAL contamination STUFF
  // MT2Analysis<MT2Estimate>* diLep_dummy_signal = new MT2Analysis<MT2Estimate>( "diLep2_dummy", cfg.regionsSet() ); 
  // (*diLep_dummy_signal) = (*TopMC2_) + (*zllMC2_) + (*signal);

  // MT2Analysis<MT2Estimate>* signal_integral;
  // MT2Analysis<MT2Estimate>* diLep_signal_integral;
  
  // signal_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "signal_integral", cfg.regionsSet(), signal );
  // diLep_signal_integral = MT2Estimate::makeIntegralAnalysisFromEstimate( "diLep_signal_integral", cfg.regionsSet(), diLep_dummy_signal );

  // MT2Analysis<MT2EstimateSyst>* purity_signal = MT2EstimateSyst::makeEfficiencyAnalysis( "purity_signal", (MT2Analysis<MT2Estimate>*)signal_integral, (MT2Analysis<MT2Estimate>*)diLep_signal_integral);
  
 


  MT2Analysis<MT2EstimateSyst>* ZinvEstimateFromZll = new MT2Analysis<MT2EstimateSyst>( "ZinvEstimateFromZll", cfg.regionsSet() );
  if( !use_extrapolation )
    (*ZinvEstimateFromZll) = (*zll_est) * (*ZinvZllRatio);
  else
    (*ZinvEstimateFromZll) = (*zll_est_integral) * (*ZinvZllRatio);

  MT2Analysis<MT2EstimateSyst>* ZinvEstimate = combineDataAndMC( ZinvEstimateFromZll, (MT2Analysis<MT2Estimate>*)Zinv );

  std::string outFile = cfg.getEventYieldDir() + "/zinvFromZll";
  outFile += ".root";

  ZinvEstimate->writeToFile( outFile, "recreate" );
  ZinvZllRatio->addToFile( outFile );
  Zinv->setName("Zinv");
  Zinv->addToFile( outFile );
  purity->addToFile( outFile );
 
  //  purity_signal->addToFile( outFile );
  // diLep->addToFile( outFile );
  //  TopMC_->addToFile( outFile );

  zllData_shape	      ->addToFile( outFile );				     
  //  zllData_shape_noTTZ ->addToFile( outFile );				     
  zllMC_shape	      ->addToFile( outFile );				     
  zllMC_shape_withTTZ ->addToFile( outFile );				     

  //  TTZ	->addToFile( outFile );				     
  TopMC_woTTZ->setName("TopMC_woTTZ");
  MC_of->setName("MC_of");
  zllData_of->addToFile( outFile );				     
  TopMC_woTTZ->addToFile( outFile );				     
  MC_of	->addToFile( outFile );				     

  if( do_hybrid ){
  purity_shape_TR->setName("purity_shape_TR");
  purity_shape_TR->addToFile( outFile);
  
  purity_forHybrid->setName("purity_forHybrid");
  purity_forHybrid->addToFile( outFile);	

//  zinvZllRatio_forShape_TR->setName("ZinvZllRatio_mt2Binned");
//  zinvZllRatio_forShape_TR->addToFile( outFile );

  bin_extrapol->addToFile( outFile );  

  ZinvEstimateFromZll_hybrid->addToFile( outFile );
  zllData_forHybrid->addToFile( outFile );

  zllData_shape_TR->addToFile( outFile );
  zinvMC_forShape_TR->addToFile( outFile );
  zllMC_shape_TR->addToFile( outFile );
  zllHybrid_shape_TR->addToFile( outFile );

  alpha->addToFile( outFile);
  }
			     
  if( !use_extrapolation && !do_hybrid ){
    zllMC->addToFile( outFile );
    zll_est->addToFile( outFile );
    purity->addToFile( outFile );
    purity_err->addToFile( outFile );
    diLep->addToFile( outFile );
    TopMC_->setName("TopMC");
    TopMC_->addToFile( outFile );
  }
  else if( use_extrapolation){
    zllMC_integral->addToFile( outFile );
    zll_est_integral->addToFile( outFile );
    purity_integral->addToFile( outFile );
    purity_err_integral->addToFile( outFile );
    diLep_integral->addToFile( outFile );
    //  TopMC_integral->addToFile( outFile );
  }

  return 0;

}






MT2Analysis<MT2Estimate>* getInclusiveRatioMC( const MT2Config& cfg, MT2Analysis<MT2EstimateTree>* Zinv, MT2Analysis<MT2EstimateTree>* gammaCRtree ) {


  std::cout << "THIS FUNCTION IS BUGGED!!! DID YOU FIX THE BUG???" << std::endl;

  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this


  MT2Analysis<MT2Estimate>* inclusiveZinv  = new MT2Analysis<MT2Estimate>( "inclusiveZinv" , cfg.regionsSet() );
  MT2Analysis<MT2Estimate>* inclusiveGamma = new MT2Analysis<MT2Estimate>( "inclusiveGamma", cfg.regionsSet() );

  std::set<MT2Region> regions = Zinv->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()>1 ) continue;

    TTree* thisTree_Zinv  = Zinv->get(*iR)->tree;
    TTree* thisTree_gamma = gammaCRtree->get(*iR)->tree;


    // loop on all regions and fill all histograms:
    for( std::set<MT2Region>::iterator iR2 = regions.begin(); iR2!=regions.end(); ++iR2 ) {

      TH1D* thisZinv  = inclusiveZinv ->get(*iR2)->yield;
      TH1D* thisGamma = inclusiveGamma->get(*iR2)->yield;

      TH1D* tmp_gamma = new TH1D(*thisGamma); // so that it gets the same binning
      tmp_gamma->SetName("tmp_gamma");
      TH1D* tmp_Zinv  = new TH1D(*thisZinv); // so that it gets the same binning
      tmp_Zinv->SetName("tmp_Zinv");

      tmp_gamma->Fill( 500., 10000.);

      thisTree_Zinv ->Project( "tmp_Zinv" , "mt2", "weight" );
      thisTree_gamma->Project( "tmp_gamma", "mt2", "weight*(prompt>1.5)" );

      thisZinv ->Add( tmp_Zinv  );
      thisGamma->Add( tmp_gamma );

      delete tmp_gamma;
      delete tmp_Zinv;

    } // for regions 2

  } // for regions 


  MT2Analysis<MT2Estimate>* inclusiveRatio = new MT2Analysis<MT2Estimate>( "ZgammaRatioMC", cfg.regionsSet() );
  (*inclusiveRatio) = (*inclusiveZinv) / (*inclusiveGamma);


  return inclusiveRatio;

}




MT2Analysis<MT2EstimateSyst>* computePurityOF( MT2Analysis<MT2Estimate>* SF, MT2Analysis<MT2Estimate>* OF, bool do_Runcert ) {

  std::string SFname = SF->getName();
  std::string OFname = OF->getName();

  // temporarily set all names to the output name so that returned MT2Analysis has consistent naming in all regions:
  std::string estimateName = "purity";
  SF->setName( estimateName );
  OF->setName( estimateName );

  std::set<MT2Region> regions = SF->getRegions();

  std::set<MT2EstimateSyst*> newData;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    // MT2EstimateSyst* Est = SF->get(*iR);
    MT2Estimate* SFEst = SF->get(*iR);
    MT2Estimate* OFEst = OF->get(*iR);
   
    MT2EstimateSyst* thisNewEstimate;
    thisNewEstimate =  new MT2EstimateSyst(*OFEst);
    for( int ibin=1; ibin<thisNewEstimate->yield->GetNbinsX()+1; ++ibin ) {
      float contentSF = SFEst->yield->GetBinContent(ibin);
      float contentOF = OFEst->yield->GetBinContent(ibin);

      float R_sfof = 1.11;
      float R_sfof_err = 0.15;

      float purity = 0.0;
      if(contentSF>(R_sfof * contentOF))
	purity = (contentSF - R_sfof * contentOF)/ contentSF;
      float purity_err = 1.0;
      if( purity > 0.){
	if(do_Runcert)
	  purity_err = sqrt( sqrt( contentOF )*sqrt( contentOF )*R_sfof*R_sfof + R_sfof_err*R_sfof_err* contentOF * contentOF  )/ contentSF;
	else
	  purity_err = sqrt( contentOF )/ contentSF;
      }
      if( contentOF == 0 )
	purity_err = 1.8/ contentSF; //uncert = 1.8 in case of 0 events

      thisNewEstimate->yield->SetBinContent( ibin, purity );
     
      thisNewEstimate->yield->SetBinError( ibin, purity_err );
      
      thisNewEstimate->yield_systUp->SetBinContent( ibin, 1+purity_err );
      thisNewEstimate->yield_systDown->SetBinContent( ibin, 1-purity_err );
    }
    newData.insert( thisNewEstimate );
  }
  MT2Analysis<MT2EstimateSyst>* analysis = new MT2Analysis<MT2EstimateSyst>( estimateName, newData );

  // set names back to original:
  SF->setName( SFname );
  OF->setName( OFname );

  return analysis;
}






MT2Analysis<MT2EstimateSyst>* combineDataAndMC( MT2Analysis<MT2EstimateSyst>* data, MT2Analysis<MT2Estimate>* mc ) {

  std::string dataname = data->getName();
  std::string mcname = mc->getName();

  // temporarily set all names to the output name so that returned MT2Analysis has consistent naming in all regions:
  std::string estimateName = "ZinvEstimate";
  data->setName( estimateName );
  mc->setName( estimateName );

  std::set<MT2Region> regions = data->getRegions();

  std::set<MT2EstimateSyst*> newData;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateSyst* dataEst = data->get(*iR);
    MT2Estimate* mcEst = mc->get(*iR);

    MT2EstimateSyst* thisNewEstimate;
    if( (!(use_extrapolation) && iR->nBJetsMin()>1) || iR->nBJetsMin()>3 /*|| iR->nJetsMin()==1*/ ) {
    //    if( (!(use_extrapolation) && iR->nBJetsMin()>1) || iR->nBJetsMin()>2 /*|| iR->nJetsMin()==1*/ ) {
      
      thisNewEstimate =  new MT2EstimateSyst(*mcEst);
      for( int ibin=1; ibin<thisNewEstimate->yield->GetNbinsX()+1; ++ibin ) {
        thisNewEstimate->yield_systUp->SetBinContent( ibin, 2.*thisNewEstimate->yield->GetBinContent(ibin) );
        thisNewEstimate->yield_systDown->SetBinContent( ibin, 0. );
      }
    } else {
      thisNewEstimate =  new MT2EstimateSyst(*dataEst);
    }
    newData.insert( thisNewEstimate );

  }

  MT2Analysis<MT2EstimateSyst>* analysis = new MT2Analysis<MT2EstimateSyst>( estimateName, newData );

  // set names back to original:
  data->setName( dataname );
  mc->setName( mcname );


  return analysis;

}





void extrapolToTopoRegion( MT2Analysis<MT2Estimate>* shape_TR, MT2Analysis<MT2Estimate>* shape, bool isMC ) {

  std::set<MT2Region> regions       = shape_TR->getRegions();
  std::set<MT2Region> regions_shape = shape   ->getRegions();
  //////////////loop over the SR regions, fill if contained//////////////
  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
 
    MT2Region* regionToMatch = new MT2Region( *iR );

    TH1D* this_shape_TR    = shape_TR   ->get( *iR)->yield;
    int nBins = this_shape_TR->GetNbinsX();

    //Loop over the shape regions///////////FILL CORRECTLY THE DATA SHAPE (NON NORMALIZED)
    for( std::set<MT2Region>::iterator iR_shape=regions_shape.begin(); iR_shape!=regions_shape.end(); ++iR_shape ) { 
      MT2Region* regionToMatch_shape = new MT2Region( *iR_shape );
      if( !(regionToMatch->MT2Region::isIncluded(regionToMatch_shape)) && !(regionToMatch->nJetsMin()==2 && (regionToMatch->nJetsMax()==6 || regionToMatch->nJetsMax()==-1)  )  ) continue; //not contained, doesn't matter
      if( regionToMatch->nJetsMin()==2 && regionToMatch->nJetsMax()==6 ){
	if( regionToMatch->htMin() != regionToMatch_shape->htMin() ) continue; //HT has to match both high and low
	if( regionToMatch->htMax() != regionToMatch_shape->htMax() ) continue;
	if( (regionToMatch->nJetsMin() != regionToMatch_shape->nJetsMin()) && (regionToMatch->nJetsMax() != regionToMatch_shape->nJetsMax() ) ) continue; //lower or upper njets has to match, but not both
      }

      TH1D* this_shape     = shape   ->get( *iR_shape)->yield;
      int nBins_shape = this_shape->GetNbinsX();

      for(int iBin=1; iBin<= nBins; iBin++){

	if(iBin == nBins && nBins_shape > nBins){
	  if( regionToMatch->nJetsMin()==2 && (regionToMatch->nJetsMax()==6 || regionToMatch->nJetsMax()==-1) && this_shape_TR->GetBinContent( iBin )!=0 ){
	    	  
	    //If MC
	    if( isMC ){
	      double int_err= 0.;
	      this_shape->IntegralAndError(iBin,-1,int_err);
	      double int_err_previous_Region= 0.;
	      this_shape_TR->IntegralAndError(iBin,-1,int_err_previous_Region);
	      this_shape_TR->SetBinContent( iBin, this_shape->Integral(iBin,-1)+this_shape_TR->Integral( iBin,-1 ) );
	      this_shape_TR->SetBinError( iBin, sqrt( int_err*int_err +int_err_previous_Region *int_err_previous_Region ) );
	    }else{ 
	      this_shape_TR->SetBinError( iBin, sqrt( this_shape->Integral(iBin,-1)*this_shape->Integral(iBin,-1)+this_shape_TR->Integral( iBin,-1 )*this_shape_TR->Integral( iBin,-1 )) );
	      this_shape_TR->SetBinContent( iBin, this_shape->Integral(iBin,-1)+this_shape_TR->Integral( iBin,-1 ) );
	    }
	    
	  }else{
	    if( isMC ){
	      double int_err= 0.;
	      this_shape->IntegralAndError(iBin,-1,int_err);
	      this_shape_TR->SetBinContent( iBin, this_shape->Integral(iBin, -1) );
	      this_shape_TR->SetBinError( iBin, int_err );
	    }else{
	      this_shape_TR->SetBinContent( iBin, this_shape->Integral(iBin, -1) );
	      this_shape_TR->SetBinError( iBin, sqrt(this_shape->Integral(iBin, -1)) );
	    }
	  }
	}else{
	  if( regionToMatch->nJetsMin()==2 && (regionToMatch->nJetsMax()==6|| regionToMatch->nJetsMax()==-1) && this_shape_TR->GetBinContent( iBin )!=0 ){

	    if( !isMC ){
	      this_shape_TR->SetBinError( iBin, sqrt(this_shape->GetBinContent(iBin)+this_shape_TR->GetBinContent(iBin)));
	      this_shape_TR->SetBinContent(iBin,this_shape->GetBinContent(iBin)+this_shape_TR->GetBinContent( iBin ) );
	    }else{ 
	      this_shape_TR->SetBinError( iBin, sqrt(this_shape->GetBinError(iBin)*this_shape->GetBinError(iBin)+this_shape_TR->GetBinError(iBin)*this_shape_TR->GetBinError(iBin)));
	      this_shape_TR->SetBinContent(iBin,this_shape->GetBinContent(iBin)+this_shape_TR->GetBinContent( iBin ) );
	    }  
	  }else{
	    if( !isMC ){
	      this_shape_TR->SetBinContent( iBin, this_shape->GetBinContent(iBin));
	      this_shape_TR->SetBinError( iBin, sqrt(this_shape->GetBinContent(iBin)));
	    }else{
	      this_shape_TR->SetBinContent( iBin, this_shape->GetBinContent(iBin));
	      this_shape_TR->SetBinError( iBin, this_shape->GetBinError(iBin) );
	    }
	  }

	} 
      }//filled and fine
    }//end loop over shape regions
  }//end loop over TR 

  return;

}









void buildHybrid( MT2Analysis<MT2Estimate>* shape_hybrid, MT2Analysis<MT2Estimate>* shape_data, MT2Analysis<MT2Estimate>* shape_MCsr, MT2Analysis<MT2Estimate>* shape_MCcr, MT2Analysis<MT2Estimate>* bin_extrapol ) {

  std::cout << "In build hybrid:" << std::endl;

  std::set<MT2Region> regions       = shape_data->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
    MT2Region* region = new MT2Region( *iR );

    TH1D* this_shape_data    = (TH1D*)shape_data   ->get( *iR)->yield;
    TH1D* this_shape_MCsr    = (TH1D*)shape_MCsr   ->get( *iR)->yield;
    TH1D* this_shape_MCcr    = (TH1D*)shape_MCcr   ->get( *iR)->yield;
    TH1D* this_shape_hybrid  = (TH1D*)shape_hybrid ->get( *iR)->yield;
    //    TH1D* this_shape_hybrid  = (TH1D*)shape_hybrid ->get( *iR)->yield->Clone("hybrid_helpH");
    TH1D* this_binExtrapol   = (TH1D*)bin_extrapol ->get( *iR)->yield;

    int nBins = this_shape_data->GetNbinsX();
    //for each topo region will have a bin number indicating where we extrapolate
    int bin_extrapol = 1;
    double integral = 0.;
    double errData = 0.;
    double integralMC = 0.;
    double errMC = 0.;

    std::vector< std::string > niceNames = region->getNiceNames();
    std::cout << niceNames[0] << " " << niceNames[1] << std::endl;
    for( int iBin=nBins; iBin>= 1; iBin-- ){
      // std::cout << this_shape_data->Integral( iBin, -1)  << std::endl;
      // if( this_shape_data->Integral( iBin, -1) >= 10. ){
      if( this_shape_MCcr->Integral( iBin, -1) >= 50. ){
	if( iBin == nBins ){ //We take the full shape from data!
	  bin_extrapol = iBin+1;
	  integral = 1.;//we don't have to do a special normalization in this case
	  errData = 0.;
	  integralMC = this_shape_MCcr->IntegralAndError( iBin, -1, errMC);
	}else{
	  bin_extrapol = iBin;
	  integral = this_shape_data->IntegralAndError( iBin, -1, errData);
	  integralMC = this_shape_MCcr->IntegralAndError( iBin, -1, errMC);
	}
	break;
      }
      else{
	
	bin_extrapol = 1;
	integralMC = this_shape_MCcr->IntegralAndError( bin_extrapol, -1, errMC);
	integral  = this_shape_data->IntegralAndError( bin_extrapol, -1, errData);

      }
    }

    //Filling the histo that knows where we extrapolate
    this_binExtrapol->SetBinContent( 1, bin_extrapol );
    
    double errZinv;
    double integralZinv = this_shape_MCsr->IntegralAndError( bin_extrapol, -1, errZinv);
    double relativeErrZinv = 1.0;
    if(integralZinv>0)
      relativeErrZinv = errZinv/integralZinv;

    double relativeErrData = 1.0;
    if(integral>0)
      relativeErrData = errData/integral;

    double errShapeExt = sqrt(relativeErrData*relativeErrData+relativeErrZinv*relativeErrZinv);

    double relativeErrMC = 1.0;
    if(integralMC>0)
      relativeErrMC = errMC/integralMC;

    double errShapeMCExt = sqrt(relativeErrMC*relativeErrMC+relativeErrZinv*relativeErrZinv);

    std::cout << "extrapol bin / total bins= " << bin_extrapol << " / " << nBins << " : " << integral << " : " << integralMC << " : " << integralZinv <<std::endl;

    for(int iBin=1; iBin<= nBins; iBin++){
      double MCsr_cont;
      double MCcr_cont;
      MCsr_cont = this_shape_MCsr->GetBinContent(iBin);
      MCcr_cont = this_shape_MCcr->GetBinContent(iBin);

      if(MCcr_cont == 0) 
	std::cout << std::endl << std::endl << "EMPTY CR bin!!!" << std::endl << std::endl;
      
//      if(iBin<bin_extrapol){
//	MCsr_cont = this_shape_MCsr->GetBinContent(iBin);
//	MCcr_cont = this_shape_MCcr->GetBinContent(iBin);
//      }
//      else{
//	MCsr_cont = this_shape_MCsr->GetBinContent(iBin);
//	MCcr_cont = integralMC;
//      }

      double ratioMC_cont = MCsr_cont/MCcr_cont;
      double kMT2_zinv = MCsr_cont/integralZinv;
      
//      if(ratioMC_cont > 3){
//	std::cout << ratioMC_cont << std::endl;
//	ratioMC_cont = 3.;
//      }
      
      if( iBin<bin_extrapol && (bin_extrapol != nBins) ){
	this_shape_data ->SetBinContent(iBin, this_shape_data->GetBinContent(iBin)*ratioMC_cont);
	this_shape_MCcr ->SetBinContent(iBin, this_shape_MCcr->GetBinContent(iBin)*ratioMC_cont);
	this_shape_data ->SetBinError(iBin, this_shape_data->GetBinError(iBin)*ratioMC_cont);
	this_shape_MCcr ->SetBinError(iBin, this_shape_MCcr->GetBinError(iBin)*ratioMC_cont);
      }else{
	this_shape_data ->SetBinContent(iBin, integral*ratioMC_cont*kMT2_zinv);
	this_shape_MCcr ->SetBinContent(iBin, integralMC*ratioMC_cont*kMT2_zinv);
	this_shape_data ->SetBinError(iBin, integral*ratioMC_cont*kMT2_zinv*errShapeExt);
	this_shape_MCcr ->SetBinError(iBin, integralMC*ratioMC_cont*kMT2_zinv*errShapeMCExt);
      }

    }

    std::cout << "extrapol bin / total bins= " << bin_extrapol << " / " << nBins << " : " << this_shape_data->Integral(bin_extrapol,-1) << " : " << this_shape_MCcr->Integral(bin_extrapol,-1) << " : " << this_shape_MCsr->Integral(bin_extrapol,-1) <<std::endl;

    //And now it has to be normalized
    this_shape_MCcr  ->Scale( 1./this_shape_MCcr->Integral());
    this_shape_data->Scale( 1./this_shape_data->Integral());

//     //Normalized
//    this_shape_MCcr->Scale(this_shape_data->Integral(bin_extrapol,-1)/this_shape_MCcr->Integral(bin_extrapol,-1) );

    for(int iBin=1; iBin<= nBins; iBin++){
	this_shape_hybrid->SetBinContent(iBin, this_shape_data->GetBinContent(iBin) );
	this_shape_hybrid->SetBinError  (iBin, this_shape_data->GetBinError(iBin) );

//      if( ( bin_extrapol==nBins+1 ) || ( iBin<bin_extrapol && (bin_extrapol != nBins) ) ){
//	this_shape_hybrid->SetBinContent(iBin, this_shape_data->GetBinContent(iBin) );
//	this_shape_hybrid->SetBinError  (iBin, this_shape_data->GetBinError(iBin) );
//      }else{
//	this_shape_hybrid->SetBinContent(iBin, this_shape_MCcr->GetBinContent(iBin) );
//	this_shape_hybrid->SetBinError  (iBin, (1./sqrt(integral))*this_shape_MCcr->GetBinContent(iBin) );
//      }

    }
    if( nBins == 1) this_shape_hybrid->SetBinError( 1, 0.0 );

//    if( this_shape_hybrid->Integral() != 0 ){
//      this_shape_hybrid->Scale( 1./ this_shape_hybrid->Integral() );
//    }

//    this_shape_hybrid_->Reset();
//    for(int iBin = 1; iBin<=nBins; ++iBin){
//      this_shape_hybrid_->SetBinContent(iBin,this_shape_hybrid->GetBinContent(iBin));
//      this_shape_hybrid_->SetBinError  (iBin,this_shape_hybrid->GetBinError(iBin));
//    }
  
  }//end loop over final estimate loops

  return;

}






//////void buildHybrid( MT2Analysis<MT2Estimate>* shape_hybrid, MT2Analysis<MT2Estimate>* shape_data, MT2Analysis<MT2Estimate>* shape_MC, MT2Analysis<MT2Estimate>* MC_ratio, MT2Analysis<MT2Estimate>* bin_extrapol ) {
//////
//////  std::set<MT2Region> regions       = shape_data->getRegions();
//////
//////  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
//////    MT2Region* region = new MT2Region( *iR );
//////
//////    TH1D* this_shape_data    = shape_data   ->get( *iR)->yield;
//////    TH1D* this_shape_MC      = shape_MC     ->get( *iR)->yield;
//////    TH1D* this_shape_hybrid  = shape_hybrid ->get( *iR)->yield;
//////    TH1D* this_MC_ratio      = MC_ratio     ->get( *iR)->yield;
//////    TH1D* this_binExtrapol   = bin_extrapol ->get( *iR)->yield;
//////
//////    int nBins = this_shape_data->GetNbinsX();
//////    //for each topo region will have a bin number indicating where we extrapolate
//////    int bin_extrapol = 1;
//////    float integral = 0.;
//////
//////    std::vector< std::string > niceNames = region->getNiceNames();
//////    std::cout << niceNames[0] << " " << niceNames[1] << std::endl;
//////    for(int iBin=nBins; iBin>= 1; iBin--){
//////      //      std::cout << this_shape_data->Integral( iBin, -1)  << std::endl;
//////      //      if( this_shape_data->Integral( iBin, -1) >= 10. ){
//////      if( this_shape_MC->Integral( iBin, -1) >= 50. ){
//////	if( iBin == nBins ){ //We take the full shape from data!
//////	  bin_extrapol = iBin+1;
//////	  integral = 1.;    //we don't have to do a special normalization in this case
//////	}else{
//////	  bin_extrapol = iBin;
//////	  integral = this_shape_data->Integral( iBin, -1);
//////	}
//////	break;
//////      }
//////    }
//////
//////    //Filling the histo that knows where we extrapolate
//////    this_binExtrapol->SetBinContent( 1, bin_extrapol );
//////    
//////    std::cout << "extrapol bin / total bins= " << bin_extrapol << " / " << nBins << std::endl;
//////
//////    for(int iBin=1; iBin<= nBins; iBin++){
//////      if( iBin< bin_extrapol && (bin_extrapol != nBins) ){
//////	this_shape_data->SetBinContent(iBin, this_shape_data->GetBinContent(iBin)*this_MC_ratio->GetBinContent(iBin));
//////	this_shape_MC  ->SetBinContent(iBin, this_shape_MC->GetBinContent(iBin)*this_MC_ratio->GetBinContent(iBin));
//////	this_shape_data->SetBinError(iBin, this_shape_data->GetBinError(iBin)*this_MC_ratio->GetBinContent(iBin));
//////	this_shape_MC  ->SetBinError(iBin, this_shape_MC->GetBinError(iBin)*this_MC_ratio->GetBinContent(iBin));
//////      }else{
//////	this_shape_data->SetBinContent(iBin, this_shape_data->GetBinContent(iBin)*this_MC_ratio->GetBinContent(iBin));
//////	this_shape_MC  ->SetBinContent(iBin, this_shape_MC->GetBinContent(iBin)*this_MC_ratio->GetBinContent(iBin));
//////	this_shape_data->SetBinError(iBin, this_shape_data->GetBinError(iBin)*this_MC_ratio->GetBinContent(iBin));
//////	this_shape_MC  ->SetBinError(iBin, this_shape_MC->GetBinError(iBin)*this_MC_ratio->GetBinContent(iBin));
//////      }
//////    }
//////
//////    //And now it has to be normalized
//////    this_shape_MC  ->Scale( 1./this_shape_MC->Integral());
//////    this_shape_data->Scale( 1./this_shape_data->Integral());
//////    //Normalized
//////    this_shape_MC->Scale(this_shape_data->Integral(bin_extrapol,-1)/this_shape_MC->Integral(bin_extrapol,-1) );
//////
//////    for(int iBin=1; iBin<= nBins; iBin++){
//////      if( (bin_extrapol==nBins+1) || ( iBin< bin_extrapol && (bin_extrapol != nBins)) ){
//////	this_shape_hybrid->SetBinContent(iBin,this_shape_data->GetBinContent(iBin) );
//////	this_shape_hybrid->SetBinError(iBin, this_shape_data->GetBinError(iBin) );
//////      }else{
//////	this_shape_hybrid->SetBinContent(iBin, this_shape_MC->GetBinContent(iBin) );
//////	this_shape_hybrid->SetBinError(iBin, (1./sqrt(integral))*this_shape_MC->GetBinContent(iBin) );
//////      }
//////    }
//////    if( nBins == 1) this_shape_hybrid->SetBinError(nBins, 0.0 );
//////
//////    if( this_shape_hybrid->Integral() != 0 )
//////      this_shape_hybrid->Scale( 1./ this_shape_hybrid->Integral() );
//////
//////  }//end loop over final estimate loops
//////
//////  return;
//////
//////}


















//First version of the filling from extrapolation to topological region
// std::set<MT2Region> regions_shape = zllData_shape    ->getRegions();
// //////////////loop over the SR regions, fill if contained//////////////
// for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
//   MT2Region* regionToMatch = new MT2Region( *iR );
//   TH1D* this_zllData_shape_TR    = zllData_shape_TR   ->get( *iR)->yield;
//   TH1D* this_zllData_of_shape_TR = zllData_of_shape_TR->get( *iR)->yield;
//   TH1D* this_zllMC_shape_TR      = zllMC_shape_TR     ->get( *iR)->yield;
//   // TH1D* this_zllHybrid_shape_TR  = zllHybrid_shape_TR ->get( *iR)->yield;
//   // TH1D* this_binExtrapol = bin_extrapol ->get( *iR)->yield;
//   int nBins = this_zllData_shape_TR->GetNbinsX();
//   //Loop over the shape regions///////////FILL CORRECTLY THE DATA SHAPE (NON NORMALIZED)
//   for( std::set<MT2Region>::iterator iR_shape=regions_shape.begin(); iR_shape!=regions_shape.end(); ++iR_shape ) { 
//     MT2Region* regionToMatch_shape = new MT2Region( *iR_shape );
//     if( !(regionToMatch->MT2Region::isIncluded(regionToMatch_shape)) && !(regionToMatch->nJetsMin()==2 && (regionToMatch->nJetsMax()==6 || regionToMatch->nJetsMax()==-1)  )  ) continue; //not contained, doesn't matter
//     if( regionToMatch->nJetsMin()==2 && regionToMatch->nJetsMax()==6 ){
// 	if( regionToMatch->htMin() != regionToMatch_shape->htMin() ) continue; //HT has to match both high and low
// 	if( regionToMatch->htMax() != regionToMatch_shape->htMax() ) continue;
// 	if( (regionToMatch->nJetsMin() != regionToMatch_shape->nJetsMin()) && (regionToMatch->nJetsMax() != regionToMatch_shape->nJetsMax() ) ) continue; //lower or upper njets has to match, but not both
//     }
//     TH1D* this_zllData_shape     = zllData_shape   ->get( *iR_shape)->yield;
//     TH1D* this_zllData_of_shape  = zllData_shape_of->get( *iR_shape)->yield;
//     TH1D* this_zllMC_shape     = zllMC_shape   ->get( *iR_shape)->yield;
//     int nBins_shape = this_zllData_shape->GetNbinsX();
//     for(int iBin=1; iBin<= nBins; iBin++){
// 	if(iBin == nBins && nBins_shape > nBins){
// 	  if( regionToMatch->nJetsMin()==2 && (regionToMatch->nJetsMax()==6 || regionToMatch->nJetsMax()==-1) && this_zllData_shape_TR->GetBinContent( iBin )!=0 ){
// 	    this_zllData_shape_TR->SetBinError( iBin, sqrt( this_zllData_shape->Integral(iBin,-1)*this_zllData_shape->Integral(iBin,-1)+this_zllData_shape_TR->Integral( iBin,-1 )*this_zllData_shape_TR->Integral( iBin,-1 )) );
// 	    this_zllData_of_shape_TR->SetBinError( iBin, sqrt(this_zllData_of_shape->Integral(iBin,-1)*this_zllData_of_shape->Integral(iBin,-1)+this_zllData_of_shape_TR->Integral( iBin,-1 )*this_zllData_of_shape_TR->Integral( iBin,-1 )) );
// 	    this_zllData_shape_TR->SetBinContent( iBin, this_zllData_shape->Integral(iBin,-1)+this_zllData_shape_TR->Integral( iBin,-1 ) );
// 	    this_zllData_of_shape_TR->SetBinContent( iBin, this_zllData_of_shape->Integral(iBin,-1)+this_zllData_of_shape_TR->Integral( iBin,-1 ) );
// 	    double int_err= 0.;
// 	    this_zllMC_shape->IntegralAndError(iBin,-1,int_err);
// 	    this_zllMC_shape_TR->SetBinContent( iBin, this_zllMC_shape->Integral(iBin,-1)+this_zllMC_shape_TR->Integral( iBin,-1 ) );
// 	    this_zllMC_shape_TR->SetBinError( iBin, int_err );

// 	  }else{
	    
// 	    this_zllData_shape_TR->SetBinContent( iBin, this_zllData_shape->Integral(iBin, -1) );
// 	    this_zllData_of_shape_TR->SetBinContent( iBin, this_zllData_of_shape->Integral(iBin, -1) );
// 	    this_zllData_shape_TR->SetBinError( iBin, sqrt(this_zllData_shape->Integral(iBin, -1)) );
// 	    this_zllData_of_shape_TR->SetBinError( iBin, sqrt(this_zllData_of_shape->Integral(iBin, -1) ));   

// 	    this_zllMC_shape_TR->SetBinContent( iBin, this_zllMC_shape->Integral(iBin, -1) );
// 	    this_zllMC_shape_TR->SetBinError( iBin, sqrt(this_zllMC_shape->Integral(iBin, -1)) );
	    
// 	  }
// 	}else{

// 	  if( regionToMatch->nJetsMin()==2 && (regionToMatch->nJetsMax()==6|| regionToMatch->nJetsMax()==-1) && this_zllData_shape_TR->GetBinContent( iBin )!=0 ){
	    
// 	    this_zllData_shape_TR->SetBinError( iBin, sqrt(this_zllData_shape->GetBinContent(iBin)+this_zllData_shape_TR->GetBinContent(iBin)));
// 	    this_zllData_of_shape_TR->SetBinError(iBin,sqrt(this_zllData_of_shape->GetBinContent(iBin)+this_zllData_of_shape_TR->GetBinContent(iBin )));	 
	  
// 	    this_zllData_shape_TR->SetBinContent(iBin,this_zllData_shape->GetBinContent(iBin)+this_zllData_shape_TR->GetBinContent( iBin ) );
// 	    this_zllData_of_shape_TR->SetBinContent(iBin,this_zllData_of_shape->GetBinContent(iBin)+this_zllData_of_shape_TR->GetBinContent(iBin));
	   
// 	    this_zllMC_shape_TR->SetBinError( iBin, sqrt(this_zllMC_shape->GetBinError(iBin)*this_zllMC_shape->GetBinError(iBin)+this_zllMC_shape_TR->GetBinError(iBin)*this_zllMC_shape_TR->GetBinError(iBin)));
// 	    this_zllMC_shape_TR->SetBinContent(iBin,this_zllMC_shape->GetBinContent(iBin)+this_zllMC_shape_TR->GetBinContent( iBin ) );
	    	  
// 	  }else{
// 	    this_zllData_shape_TR->SetBinContent( iBin, this_zllData_shape->GetBinContent(iBin));
// 	    this_zllData_of_shape_TR->SetBinContent( iBin, this_zllData_of_shape->GetBinContent(iBin));
// 	    this_zllData_shape_TR->SetBinError( iBin, sqrt(this_zllData_shape->GetBinContent(iBin)));
// 	    this_zllData_of_shape_TR->SetBinError( iBin, sqrt( this_zllData_of_shape->GetBinContent(iBin)));

// 	    this_zllMC_shape_TR->SetBinContent( iBin, this_zllMC_shape->GetBinContent(iBin));
// 	    this_zllMC_shape_TR->SetBinError( iBin, this_zllMC_shape->GetBinError(iBin) );
// 	  }

// 	} 
//     }//filled and fine
//   }//end loop over shape regions
// }//end loop over TR 
