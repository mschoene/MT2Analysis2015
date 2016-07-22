#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip> 
#include <stdlib.h> 

#include "TSystem.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TList.h"
#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "RooHistError.h"

#include "interface/MT2Config.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2EstimateSigSyst.h"
#include "interface/MT2EstimateSigContSyst.h"


bool use_extrapolation = true;
bool doSignalContamination = true;
bool doSimultaneousFit = false;
bool includeSignalUnc = true; // signal lep eff commented out till available
bool copy2SE = false; // copy datacards to SE
bool doGenAverage = true;

int round(float d) {
  return (int)(floor(d + 0.5));
}

std::string getSimpleSignalName( const std::string& longName );
std::string gammaConvention( float yieldSR, int yieldCR, int position, const std::string& corrName, const std::string& uncorrName="", float testAlpha=1. );


int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|              Running createDatacards               |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc != 3 && argc != 7 && argc != 8) {
    std::cout << "USAGE: ./createDatacards [configFileName] [model] [m1] [m2] [m11] [m22] ([label])" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  float m1=0.;
  float m2=2000.;
  
  float m11=0.;
  float m22=2000;

  std::string model(argv[2]);

  if( argc == 7 || argc == 8 ){
    
    m1  = std::stof( argv[3] );
    m2  = std::stof( argv[4] );

    m11  = std::stof( argv[5] );
    m22  = std::stof( argv[6] );

  }
  
  std::string label;
  if( argc == 8 )
    label = argv[7];
  else
    label = "";

  std::cout << "Will produce datacards for parent mass between " << m1 << " and " << m2 << ", and LSP mass between" << m11 <<  " and " << m22 << std::endl;

  std::string dir = cfg.getEventYieldDir();
  std::string mc_fileName = dir + "/analyses.root";
  std::string data_fileName = dir + "/analyses.root";


  bool useMC_qcd  = false;
  bool useMC_zinv = false;
  bool useMC_llep = false;

  //  float err_qcd_uncorr  = 1.0; // 100% of QCD MC yield, if use MC for QCD

  float err_llep_shape = 0.40;
  float err_llep_lepEff = 0.07; // Uncertainty on llep estimate from lepton efficiency (7%)

  float err_zinv_shape = 0.40;
  float err_zinv_uncorr_2b = 1.0;

  float err_zinv_puritySyst = 0.1; // 10%, including 5% on purity + 8% on fragmentation
  float err_zinv_doubleRatioOffset = 0.10; // 10%, fully correlated, on zinv // TO BE UPDATED IN CASE IT CHANGES WITH FULL LUMI
  float zinv_doubleRatioOffset = 0.90; // 90%, used to correct the Z/G ratio in zinv estimate // TO BE UPDATED IN CASE IT CHANGES WITH FULL LUMI
  
  float err_lumi_corr   = 0.027; // Uncertainty on luminosity (2.7% for 2016 public results)


  // Reading data analysis (in search region)
  MT2Analysis<MT2Estimate>* data  = MT2Analysis<MT2Estimate>::readFromFile( data_fileName, "data" );
  
  // Reading QCD estimate
  MT2Analysis<MT2Estimate>* qcd;
  //  MT2Analysis<MT2Estimate>* qcd_mc;
  MT2Analysis<MT2Estimate>* qcdCR;
  MT2Analysis<MT2Estimate>* qcd_ratio;
  MT2Analysis<MT2Estimate>* qcd_purity;
  MT2Analysis<MT2Estimate>* qcd_fjets;
  // MT2Analysis<MT2Estimate>* qcd_fjets_vlht;
  MT2Analysis<MT2Estimate>* qcd_rb;
  MT2Analysis<MT2Estimate>* qcd_ratioSystFit;

  MT2Analysis<MT2Estimate>* qcd_monojet;
  MT2Analysis<MT2Estimate>* qcdCR_monojet;
  MT2Analysis<MT2Estimate>* qcd_ratio_monojet;


  if( useMC_qcd )
    qcd = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "QCD"  );
  else{

    qcd = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateData.root", "qcdEstimate" );
    qcdCR = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateData.root", "nCR" );
    qcd_ratio = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateData.root", "r_effective" );
    qcd_purity = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateData.root", "qcdPurity" );
    qcd_fjets = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateData.root", "f_jets_data" );
    //qcd_fjets_vlht = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateData.root", "f_jets_data_noPS" );
    qcd_rb = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateData.root", "r_hat_data" );
    qcd_ratioSystFit = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateData.root", "r_systFit" );
    //    qcd_mc = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "QCD"  );

    qcd_monojet = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateMonojet.root", "monojet_qcdEstimate" );
    qcdCR_monojet = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateMonojet.root", "monojet_nCR" );
    qcd_ratio_monojet = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdControlRegion/qcdEstimateMonojet.root", "monojet_r" );
    
  }



  // Reading invisible Z estimate
  MT2Analysis<MT2Estimate>* zinv;
  MT2Analysis<MT2Estimate>* zinvCR;
  MT2Analysis<MT2Estimate>* zinv_ratio;
  MT2Analysis<MT2EstimateSyst>* purity;

  MT2Analysis<MT2EstimateSyst>* zllG_ht;
  MT2Analysis<MT2EstimateSyst>* zllG_nJets;
  MT2Analysis<MT2EstimateSyst>* zllG_nBJets;
  MT2Analysis<MT2EstimateSyst>* zllG_ht_monojet;

  MT2Analysis<MT2EstimateSyst>* zllG_mc_ht;
  MT2Analysis<MT2EstimateSyst>* zllG_mc_nJets;
  MT2Analysis<MT2EstimateSyst>* zllG_mc_nBJets;
  MT2Analysis<MT2EstimateSyst>* zllG_mc_ht_monojet;

  if( useMC_zinv )
    zinv = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "ZJets");
  else {
    
    zinvCR      = MT2Analysis<MT2Estimate>    ::readFromFile( dir + "/gammaControlRegion/data.root", "gammaCR");
    
    zinv        = MT2Analysis<MT2Estimate>    ::readFromFile( dir + "/zinvFromGamma.root", "ZinvEstimate");
    
    zinv_ratio  = MT2Analysis<MT2Estimate>    ::readFromFile( dir + "/zinvFromGamma.root", "ZgammaRatio");
    
    purity      = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zinvFromGamma.root", "purity");
    
    zllG_ht_monojet  = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_data_ratio.root", "zllG_data_mono_ht");
    zllG_ht          = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_data_ratio.root", "zllG_data_ht");
    zllG_nJets       = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_data_ratio.root", "zllG_data_nJets");
    zllG_nBJets      = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_data_ratio.root", "zllG_data_nBJets");
    
    zllG_mc_ht_monojet  = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_mc_ratio.root", "zllG_mc_mono_ht");
    zllG_mc_ht          = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_mc_ratio.root", "zllG_mc_ht");
    zllG_mc_nJets       = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_mc_ratio.root", "zllG_mc_nJets");
    zllG_mc_nBJets      = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_mc_ratio.root", "zllG_mc_nBJets");
    
  }
  zinv->setName("zinv");
  //zinv->addToFile( mc_fileName, true ); // Optionally, to add estimate used for invisible Z estimate to analyses.root


  
  // Reading lost lepton estimate
  MT2Analysis<MT2Estimate>* llep;
  MT2Analysis<MT2Estimate>* llepCR;
  MT2Analysis<MT2Estimate>* llep_ratio;

  if( useMC_llep ) {
    
    MT2Analysis<MT2Estimate>* wjets = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "WJets");
    MT2Analysis<MT2Estimate>* top   = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "Top");
    llep = new MT2Analysis<MT2Estimate>( (*wjets) + (*top) );
    
    llepCR = MT2Analysis<MT2Estimate>::readFromFile( cfg.getEventYieldDir() + "/llepControlRegion/mc.root", "llepCR" );

  } 
  else {
   
    llep = MT2Analysis<MT2Estimate>::readFromFile( dir + "/llepEstimate.root", "llepEstimate" );
    llep_ratio = MT2Analysis<MT2Estimate>::readFromFile( dir + "/llepEstimate.root", "llepRatioMC" );
    llepCR = MT2Analysis<MT2Estimate>::readFromFile( cfg.getEventYieldDir() + "/llepControlRegion/data.root", "llepCR" );
  
  }
  llep->setName( "llep" );
  //llep->addToFile( mc_fileName, true ); // Optionally, to add estimate used for invisible Z estimate to analyses.root



  // Getting region set used (from data)
  std::set<MT2Region> regions = data->getRegions();

  // Getting inclusive region set (as used for Zll/Gamma ratio)
  std::set<MT2Region> inclusiveRegions=  zllG_ht->getRegions();
  MT2Region inclusiveRegion( (*inclusiveRegions.begin() ) );


  // Reading Zll/G ratio as a function of HT/NJ/NB/HT(monojet)
  TH1D* thisUp_zllG_ht     = zllG_ht    ->get(inclusiveRegion)->yield_systUp;
  TH1D* thisUp_zllG_nJets  = zllG_nJets ->get(inclusiveRegion)->yield_systUp;
  TH1D* thisUp_zllG_nBJets = zllG_nBJets->get(inclusiveRegion)->yield_systUp;
  TH1D* thisUp_zllG_ht_monojet = zllG_ht_monojet->get(inclusiveRegion)->yield_systUp;

  TH1D* thisDn_zllG_ht     = zllG_ht    ->get(inclusiveRegion)->yield_systDown;
  TH1D* thisDn_zllG_nJets  = zllG_nJets ->get(inclusiveRegion)->yield_systDown;
  TH1D* thisDn_zllG_nBJets = zllG_nBJets->get(inclusiveRegion)->yield_systDown;
  TH1D* thisDn_zllG_ht_monojet = zllG_ht_monojet->get(inclusiveRegion)->yield_systDown;

  TH1D* this_zllG_ht     = zllG_ht    ->get(inclusiveRegion)->yield;
  TH1D* this_zllG_nJets  = zllG_nJets ->get(inclusiveRegion)->yield;
  TH1D* this_zllG_nBJets = zllG_nBJets->get(inclusiveRegion)->yield;
  TH1D* this_zllG_ht_monojet = zllG_ht_monojet->get(inclusiveRegion)->yield;

  TH1D* this_zllG_mc_ht         = zllG_mc_ht    ->get(inclusiveRegion)->yield;
  TH1D* this_zllG_mc_nJets      = zllG_mc_nJets ->get(inclusiveRegion)->yield;
  TH1D* this_zllG_mc_nBJets     = zllG_mc_nBJets->get(inclusiveRegion)->yield;
  TH1D* this_zllG_mc_ht_monojet = zllG_mc_ht_monojet->get(inclusiveRegion)->yield;

  
  // First create template datacards
  std::string path_templ = dir + "/datacard_templates";
  system(Form("mkdir -p %s", path_templ.c_str()));

  
  // Start loop over topological regions
  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) { 
    // Getting data yield histogram
    TH1D* this_data = data->get(*iR)->yield;
    
    // Getting QCD yield histograms, plus uncertainties from QCD estimate
    TH1D* this_qcd;
    TH1D* this_qcdCR;
    TH1D* this_qcd_ratio;
    TH1D* this_qcd_purity;
    TH1D* this_qcd_fjets;
    TH1D* this_qcd_rb;
    TH1D* this_qcd_ratioSystFit;
    
    std::string qcd_fjetsCR_name;       
    std::string qcd_rbCR_name;       
    
    if( iR->nJetsMax()==1 ){ // First QCD estimate from monojet
      
      this_qcd = qcd_monojet->get(*iR)->yield;
      this_qcdCR = qcdCR_monojet->get(*iR)->yield;
      this_qcd_ratio = qcd_ratio_monojet->get(*iR)->yield;
      
    }
    else{ // Then QCD estimate for multijet
      
      this_qcd = qcd->get(*iR)->yield;
      this_qcdCR = qcdCR->get(*iR)->yield;
      this_qcd_ratio = qcd_ratio->get(*iR)->yield;
      this_qcd_purity = qcd_purity->get(*iR)->yield;
      this_qcd_ratioSystFit = qcd_ratioSystFit->get(*iR)->yield;
      
      // not anymore due to new PFHT125 trigger 
      // if( iR->htMin() < 450 ){ // For very low HT, read different analysis for F(J)
      // 	MT2Region* thisQCDCR = qcd_fjets_vlht->matchRegion(*iR);
      // 	this_qcd_fjets = qcd_fjets_vlht->get(*thisQCDCR)->yield;
      // 	qcd_fjetsCR_name = thisQCDCR->getName();
      // }
      // else{ // From low HT, read from main analysis
	
      MT2Region* thisQCDCR = qcd_fjets->matchRegion(*iR);
      this_qcd_fjets = qcd_fjets->get(*thisQCDCR)->yield;
      
      qcd_fjetsCR_name = thisQCDCR->getName();
	
      // }
      
      // For region 2-6j, 3b, take R(B) from 4-6 jets region
      //MT2Region* thisQCDCR;
      if( iR->nBJetsMin()==3 && iR->nJetsMin()==2 ) 
	thisQCDCR = new MT2Region( 200, -1, 4, 6, 0, -1 );
      else  
	thisQCDCR = qcd_rb->matchRegion(*iR);
      
      this_qcd_rb = qcd_rb->get(*thisQCDCR)->yield; // Get R(B) fraction
      
      qcd_rbCR_name = thisQCDCR->getName();
      
    }
    
    // Get histograms for invisible Z estimate
    TH1D* this_zinv = zinv->get(*iR)->yield;
    TH1D* this_zinv_ratio =  zinv_ratio->get(*iR)->yield;
    TH1D* this_zinvCR;
    this_zinvCR = zinvCR->get(*iR)->yield;
    
    TGraphAsymmErrors* this_zinv_purity;
    this_zinv_purity = purity->get(*iR)->getGraph();
    
    int Ngamma=0; // Initialize variable for number of events in photon CR
    // Get histograms for lost lepton estimate
    TH1D* this_llep = llep->get(*iR)->yield;
    TH1D* this_llep_ratio = llep_ratio->get(*iR)->yield;
    TH1D* this_llepCR = llepCR->get(*iR)->yield;
    
    float N_llep_CR = this_llepCR->Integral();
    // Regions with N(J)>=7 and N(B)>=1 share same control region (1-2 b jets)
    std::string llepCR_name;
    if(iR->nJetsMin()>=7 && iR->nBJetsMin()>=1){
      MT2Region* thisCR = new MT2Region(iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), 1, 2);
      llepCR_name = thisCR->getName();
    }
    else
      llepCR_name = iR->getName();
    
    
    int nBins = this_data->GetNbinsX(); // Getting total number of bins for this topological region
    
    // Calculating shape uncertainty for invisible Z (default: linear extrapolation)
    float shapeErr_zinv=0.;
    for( int iBin=1; iBin<this_data->GetNbinsX()+1; ++iBin ) {
      
      float relativeErr;
      if( fabs(this_zinv->GetBinContent(iBin))>0 )
	relativeErr = err_zinv_shape / (nBins-1) * (iBin-1);
      //	 relativeErr = err_zinv_shape / ((nBins-1) * (nBins-1)) * (iBin-1) * (iBin-1); // Parabolic shape uncertainty
      else
	relativeErr = 0.0;
      
      shapeErr_zinv+=relativeErr*fabs(this_zinv->GetBinContent(iBin));
      
    }

    // Calculating shape uncertainty for lost lepton (default: linear extrapolation)   
    float shapeErr_llep=0.;
    for( int iBin=1; iBin<this_data->GetNbinsX()+1; ++iBin ) {
      
      float relativeErr;
      if( this_llep->GetBinContent(iBin)>0 )
	relativeErr = err_llep_shape / (nBins-1) * (iBin-1);
      //	 relativeErr = err_llep_shape / ((nBins-1) * (nBins-1)) * (iBin-1) * (iBin-1); // Parabolic shape uncertainty
      else
	relativeErr = 0.0;
      
      shapeErr_llep+=relativeErr*this_llep->GetBinContent(iBin);
      
    }
    
    float lastR_zinv; // To keep information on last non-zero ratio. Will use last non-zero ratio if ratio for one bin is zero. 
    float lastR_llep; // To keep information on last non-zero ratio. Will use last non-zero ratio if ratio for one bin is zero. 
    float lastR_qcd_mono; // To keep information on last non-zero ratio. Will use last non-zero ratio if ratio for one bin is zero. 
    
    // Start loop over MT2 bins for this topological region
    for( int iBin=1; iBin<this_data->GetNbinsX()+1; ++iBin ) {
      
      bool includeCR=false; // In case do simultaneous fit of control region (for signal contamination), need to include CR into datacard
      if(iBin==1) includeCR=true;
      if(iR->nJetsMin()>=7 && iR->nBJetsMin()>1) includeCR=false;
      
      
      // If MT2 bin low edge > max HT continue (should never happen by construction)
      if(this_data->GetBinLowEdge( iBin ) > iR->htMax() && iR->htMax()>0 ) continue;
      
      // Getting MT2 bin edges
      float mt2Min = this_data->GetBinLowEdge( iBin );
      float mt2Max = (iBin==this_data->GetNbinsX()) ?  -1. : this_data->GetBinLowEdge( iBin+1 );
      
      // Getting bin name
      std::string binName;
//      if( iR->nJetsMax() ==1 ){ 
//	binName = std::string( Form("%s_m0toInf", iR->getName().c_str() ) ); // Needed to cope with SnT input, to remove for ETH input
//      }
//      else{
	
	if( mt2Max>=0. )
	  binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
	else
	  binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );
	
	//      }
      
      // Getting HT region name
      std::string htName;
      htName = iR->htRegion()->getName();

      // Getting NJ name
      int nJetsMin = iR->nJetsMin();
      int nJetsMax = iR->nJetsMax();
      std::string jName;
      jName = iR->sigRegion()->getSingleJetString( "j", nJetsMin,  nJetsMax  );
      
      // Getting NB name
      int nBJetsMin = iR->nBJetsMin();
      int nBJetsMax = iR->nBJetsMax();
      std::string bName;
      bName = iR->sigRegion()->getSingleJetString( "b", nBJetsMin,  nBJetsMax  );
      
      // Set template datacard name for this bin
      std::string datacardName( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
      
      std::ifstream thisDatacard( datacardName.c_str() );
      if( thisDatacard.good() ) continue; // If template already exists, move on
      
      std::ofstream datacard( datacardName.c_str() );
      
      // Initializing table for estimates + observation
      std::string tableName( Form("%s/table_%s.txt", path_templ.c_str(), binName.c_str()) );
      std::ofstream table( tableName.c_str() );
      table << std::setprecision(3);
      
      // Set name for 1L CR (needed for simultaneous fit)
      std::string binName1L( Form("%s_1L", binName.c_str()) );

      // Start writing template card

      if(doSimultaneousFit && includeCR){ 
	
	datacard << "imax 2" << std::endl; //Number of bins
	datacard << "jmax 3" << std::endl; //Number of backgrounds
	datacard << "kmax *" << std::endl; //Number of nuisances
	datacard << "-------------" << std::endl;
	datacard << std::endl << std::endl;
	
	
	datacard << std::fixed;
	datacard << std::setprecision(3) << std::endl << std::endl;
	datacard << "bin  " << binName << "\t" << llepCR_name << std::endl;
	//       datacard << "observation  " << (fabs(this_qcd->GetBinContent(iBin))+fabs(this_zinv->GetBinContent(iBin))+fabs(this_llep->GetBinContent(iBin))) << std::endl;
	datacard << "observation  " << this_data->GetBinContent(iBin) << "\t" << N_llep_CR << std::endl;
	datacard << "-------------" << std::endl;
	datacard << std::endl << std::endl;
	
      }
      
      else{
	
	datacard << "imax 1" << std::endl; //Number of bins
	datacard << "jmax 3" << std::endl; //Number of backgrounds
	datacard << "kmax *" << std::endl; //Number of nuisances
	datacard << "-------------" << std::endl;
	datacard << std::endl << std::endl;
	
	
	datacard << std::fixed;
	datacard << std::setprecision(3) << std::endl << std::endl;
	datacard << "bin  " << binName << std::endl;
	//       datacard << "observation  " << (fabs(this_qcd->GetBinContent(iBin))+fabs(this_zinv->GetBinContent(iBin))+fabs(this_llep->GetBinContent(iBin))) << std::endl;
	datacard << "observation  " << this_data->GetBinContent(iBin) << std::endl;
	datacard << "-------------" << std::endl;
	datacard << std::endl << std::endl;	 
	
      }
      // Read background estimates values
      float yield_llep = fabs(this_llep->GetBinContent(iBin));
      float yield_qcd = fabs(this_qcd ->GetBinContent(iBin));
      float yield_zinv = fabs(this_zinv->GetBinContent(iBin));
      
      
      if(doSimultaneousFit && includeCR){
	
	// sig qcd zinv llep sig1L llep1L
	datacard << "bin \t" << binName << "\t" << binName << "\t" << binName << "\t" << binName << "\t" << llepCR_name << "\t" << llepCR_name <<std::endl;
	datacard << "process \t sig \t zinv \t llep \t qcd \t sig \t llep" << std::endl;
	datacard << "process \t 0 \t 1 \t 2 \t 3 \t 0 \t 2" << std::endl;
	datacard << "rate \t XXX";
	datacard << " \t " << yield_zinv << " \t " << yield_llep << " \t " << yield_qcd << " \t YYY \t " << ( (N_llep_CR>0.001) ? N_llep_CR : 0.001 ) << std::endl;
	datacard << "-------------" << std::endl;
	
	datacard << "lumi_syst    lnN    " << 1.+err_lumi_corr << " - - - " << 1.+err_lumi_corr << " -" << std::endl;
	datacard << "sig_MCstat_" << binName << " lnN UUU - - - - -" << std::endl;
	if(doSignalContamination && doSimultaneousFit)
	  datacard << "sig_MCstat_1L_" << llepCR_name << " lnN - - - - VVV -" << std::endl;
	if(!includeSignalUnc)
	  datacard << "sig_syst_" << binName << " lnN 1.2 - - - - -" << std::endl;
	else{
	  datacard << "sig_isrSyst lnN III - - - - -" << std::endl;
	  datacard << "sig_bTagHeavySyst lnN HHH - - - - -" << std::endl;
	  datacard << "sig_bTagLightSyst lnN LLL - - - - -" << std::endl;
	  if(model=="T2tt" || model=="T1tttt")
	    datacard << "sig_lepEffSyst lnN EEE - - - - -" << std::endl; // Include lepton eff. uncertainty only for T2tt and T1tttt
	}
      }
      
      else {
	
	// sig qcd zinv llep
	datacard << "bin \t" << binName << "\t" << binName << "\t" << binName << "\t" << binName << std::endl;
	datacard << "process \t sig \t zinv \t llep \t qcd" << std::endl;
	datacard << "process \t 0 \t 1 \t 2 \t 3" << std::endl;
	datacard << "rate \t XXX";
	datacard << " \t " << yield_zinv << " \t " << yield_llep << " \t " << yield_qcd << std::endl;
	datacard << "-------------" << std::endl;
	
	datacard << "lumi_syst    lnN    " << 1.+err_lumi_corr << " - - -" << std::endl;
	datacard << "sig_MCstat_" << binName << " lnN UUU - - -" << std::endl;
	if (doGenAverage)
	  datacard << "sig_gensyst lnU SSS - - -" << std::endl;
	if(!includeSignalUnc)
	  datacard << "sig_syst_" << binName << " lnN 1.2 - - -" << std::endl;
	else{
	  datacard << "sig_isrSyst lnN III - - -" << std::endl;
	  datacard << "sig_bTagHeavySyst lnN HHH - - -" << std::endl;
	  datacard << "sig_bTagLightSyst lnN LLL - - -" << std::endl;
	  if(model=="T2tt" || model=="T1tttt")
	    datacard << "sig_lepEffSyst lnN EEE - - -" << std::endl; // Include lepton eff. uncertainty only for T2tt and T1tttt 
	}
      }
      
      // Initialize variables for tables
      float zinv_statUp = 0.;
      float zinv_statDn = 0.;
      float zinv_systUp = 0.;
      float zinv_systDn = 0.;
      
      float qcd_statUp = 0.;
      float qcd_statDn = 0.;
      float qcd_systUp = 0.;
      float qcd_systDn = 0.;
      
      float llep_statUp = 0.;
      float llep_statDn = 0.;
      float llep_systUp = 0.;
      float llep_systDn = 0.;
      
      int zinv_nCR = 0;
      int qcd_nCR  = 0;
      int llep_nCR = 0;
      
   
      std::string zinvCR_name;
      zinvCR_name = iR->getName();
      
      // Z INVISIBLE SYSTEMATICS:
      if( yield_zinv>=0.) {
	
	if( (!use_extrapolation && iR->nBJetsMin()<2) || iR->nBJetsMin()<=3 ) { // 0 and 1 btag if not using extrapolation, all otherwise
	  
	  // Get uncertainty from Zll/G ratio in data
	  int thisBinNJ;
	  int thisBinNB;
	  int thisBinHT;
	  
	  float thisErrNJUp;
	  float thisErrNBUp;
	  float thisErrHTUp;
	  
	  float thisErrNJDn;
	  float thisErrNBDn;
	  float thisErrHTDn;
	  
	  float thisCentralNJ;
	  float thisCentralNJ_mc;

	  float thisCentralNB;
	  float thisCentralNB_mc;

	  float thisCentralHT;
	  float thisCentralHT_mc;

	  if( iR->nJetsMax()>1 || iR->nJetsMax()<0){
	    
	    thisBinNJ = this_zllG_nJets->FindBin(iR->nJetsMin());
	    if( iR->nBJetsMin() < 3 )
	      thisBinNB = this_zllG_nBJets->FindBin(iR->nBJetsMin());
	    else
	      thisBinNB = this_zllG_nBJets->FindBin(2); //If NB>=3 take uncertainty from 2B (stats)
	    thisBinHT = this_zllG_ht->FindBin(iR->htMin()+1.);
	    
	    thisErrNJUp = ( this_zllG_nJets->GetBinContent(thisBinNJ) > 0 )  ? ( thisUp_zllG_nJets->GetBinContent(thisBinNJ) - this_zllG_nJets->GetBinContent(thisBinNJ) )  / this_zllG_nJets->GetBinContent(thisBinNJ)  : 1.0;
	    thisErrNBUp = ( this_zllG_nBJets->GetBinContent(thisBinNB) > 0 ) ? ( thisUp_zllG_nBJets->GetBinContent(thisBinNB) - this_zllG_nBJets->GetBinContent(thisBinNB) )/ this_zllG_nBJets->GetBinContent(thisBinNB) : 1.0;
	    thisErrHTUp = ( this_zllG_ht->GetBinContent(thisBinHT) > 0 )     ? ( thisUp_zllG_ht->GetBinContent(thisBinHT) - this_zllG_ht->GetBinContent(thisBinHT) )           / this_zllG_ht->GetBinContent(thisBinHT)     : 1.0;
	    
	    thisErrNJDn = ( this_zllG_nJets->GetBinContent(thisBinNJ) > 0 )  ? ( this_zllG_nJets->GetBinContent(thisBinNJ)  - thisDn_zllG_nJets->GetBinContent(thisBinNJ) ) / this_zllG_nJets->GetBinContent(thisBinNJ)  : 1.0;
	    thisErrNBDn = ( this_zllG_nBJets->GetBinContent(thisBinNB) > 0 ) ? ( this_zllG_nBJets->GetBinContent(thisBinNB) - thisDn_zllG_nBJets->GetBinContent(thisBinNB) )/ this_zllG_nBJets->GetBinContent(thisBinNB) : 1.0;
	    thisErrHTDn = ( this_zllG_ht->GetBinContent(thisBinHT) > 0 )     ? ( this_zllG_ht->GetBinContent(thisBinHT)     - thisDn_zllG_ht->GetBinContent(thisBinHT) )    / this_zllG_ht->GetBinContent(thisBinHT)     : 1.0;
	    
	    if( iR->nBJetsMin() >= 3 ){ // If NB>=3 take uncertainty from 2B (stats), THEN DOUBLE IT
	      thisErrNBUp*=2;
	      thisErrNBDn*=2;
	    }
	    
	    thisCentralNJ = this_zllG_nJets->GetBinContent(thisBinNJ);
	    thisCentralNJ_mc = this_zllG_mc_nJets->GetBinContent(thisBinNJ);

	    thisCentralNB = this_zllG_nBJets->GetBinContent(thisBinNB);
	    thisCentralNB_mc = this_zllG_mc_nBJets->GetBinContent(thisBinNB);

	    thisCentralHT = this_zllG_ht->GetBinContent(thisBinHT);
	    thisCentralHT_mc = this_zllG_mc_ht->GetBinContent(thisBinHT);

	  }
	  else{
	    
	    thisBinNJ = this_zllG_nJets->FindBin(iR->nJetsMin());
	    thisBinNB = this_zllG_nBJets->FindBin(iR->nBJetsMin());
	    
	    thisErrNJUp = ( this_zllG_nJets->GetBinContent(thisBinNJ) > 0 )  ? ( thisUp_zllG_nJets->GetBinContent(thisBinNJ) - this_zllG_nJets->GetBinContent(thisBinNJ) )  / this_zllG_nJets->GetBinContent(thisBinNJ)  : 1.0;
	    thisErrNBUp = ( this_zllG_nBJets->GetBinContent(thisBinNB) > 0 ) ? ( thisUp_zllG_nBJets->GetBinContent(thisBinNB) - this_zllG_nBJets->GetBinContent(thisBinNB) )/ this_zllG_nBJets->GetBinContent(thisBinNB) : 1.0;
	    
	    thisErrNJDn = ( this_zllG_nJets->GetBinContent(thisBinNJ) > 0 )  ? ( this_zllG_nJets->GetBinContent(thisBinNJ)  - thisDn_zllG_nJets->GetBinContent(thisBinNJ) ) / this_zllG_nJets->GetBinContent(thisBinNJ)  : 1.0;
	    thisErrNBDn = ( this_zllG_nBJets->GetBinContent(thisBinNB) > 0 ) ? ( this_zllG_nBJets->GetBinContent(thisBinNB) - thisDn_zllG_nBJets->GetBinContent(thisBinNB) )/ this_zllG_nBJets->GetBinContent(thisBinNB) : 1.0;
	    
	    thisBinHT = this_zllG_ht_monojet->FindBin(iR->htMin()+1.);
	    thisErrHTUp = ( this_zllG_ht_monojet->GetBinContent(thisBinHT) > 0 )     ? ( thisUp_zllG_ht_monojet->GetBinContent(thisBinHT) - this_zllG_ht_monojet->GetBinContent(thisBinHT) ) / this_zllG_ht_monojet->GetBinContent(thisBinHT)     : 1.0;
	    thisErrHTDn = ( this_zllG_ht_monojet->GetBinContent(thisBinHT) > 0 )     ? ( this_zllG_ht_monojet->GetBinContent(thisBinHT)   - thisDn_zllG_ht_monojet->GetBinContent(thisBinHT) ) / this_zllG_ht_monojet->GetBinContent(thisBinHT)     : 1.0;
	    
	    
	    thisCentralNJ = this_zllG_nJets->GetBinContent(thisBinNJ);
	    thisCentralNJ_mc = this_zllG_mc_nJets->GetBinContent(thisBinNJ);

	    thisCentralNB = this_zllG_nBJets->GetBinContent(thisBinNB);
	    thisCentralNB_mc = this_zllG_mc_nBJets->GetBinContent(thisBinNB);

	    thisCentralHT = this_zllG_ht_monojet->GetBinContent(thisBinHT);
	    thisCentralHT_mc = this_zllG_mc_ht_monojet->GetBinContent(thisBinHT);
	    
	  }
	  
	  
	  ////// NEW in 2016:
	  //If stat. uncertainty in Z/G ratio (from data) is smaller than data-MC (AFTER correction), then take data-MC difference as uncertainty
	  thisCentralNJ_mc *= zinv_doubleRatioOffset;
	  thisCentralNB_mc *= zinv_doubleRatioOffset;
	  thisCentralHT_mc *= zinv_doubleRatioOffset;
	  
	  if( thisCentralNJ_mc > thisCentralNJ*(1+thisErrNJUp) ) 
	    thisErrNJUp = (thisCentralNJ>0) ? (thisCentralNJ_mc - thisCentralNJ)/thisCentralNJ : 1.0;
	  else if ( thisCentralNJ_mc < thisCentralNJ*(1-thisErrNJDn) ) 
	    thisErrNJDn = (thisCentralNJ>0) ? (thisCentralNJ - thisCentralNJ_mc)/thisCentralNJ : 1.0;

	  if( thisCentralNB_mc > thisCentralNB*(1+thisErrNBUp) ) 
	    thisErrNBUp = (thisCentralNB>0) ? (thisCentralNB_mc - thisCentralNB)/thisCentralNB : 1.0;
	  else if ( thisCentralNB_mc < thisCentralNB*(1-thisErrNBDn) ) 
	    thisErrNBDn = (thisCentralNB>0) ? (thisCentralNB - thisCentralNB_mc)/thisCentralNB : 1.0;

	  if( thisCentralHT_mc > thisCentralHT*(1+thisErrHTUp) ) 
	    thisErrHTUp = (thisCentralHT>0) ? (thisCentralHT_mc - thisCentralHT)/thisCentralHT : 1.0;
	  else if ( thisCentralHT_mc < thisCentralHT*(1-thisErrHTDn) ) 
	    thisErrHTDn = (thisCentralHT>0) ? (thisCentralHT - thisCentralHT_mc)/thisCentralHT : 1.0;
	  //////

	  if(doSimultaneousFit && includeCR)
	    datacard << "zinv_doubleRatioOffset lnN   - " << 1.+err_zinv_doubleRatioOffset << " - - - -" << std::endl;
	  else
	    datacard << "zinv_doubleRatioOffset lnN   - " << 1.+err_zinv_doubleRatioOffset << " - -" << std::endl;
	  
	  zinv_systUp += err_zinv_doubleRatioOffset*err_zinv_doubleRatioOffset;
	  zinv_systDn += err_zinv_doubleRatioOffset*err_zinv_doubleRatioOffset;
	  
	  float err_ZGUp = TMath::Sqrt( thisErrNJUp*thisErrNJUp + thisErrNBUp*thisErrNBUp + thisErrHTUp*thisErrHTUp );
	  float err_ZGDn = TMath::Sqrt( thisErrNJDn*thisErrNJDn + thisErrNBDn*thisErrNBDn + thisErrHTDn*thisErrHTDn );
	  
	  err_ZGDn = (err_ZGDn < 1) ? err_ZGDn : 0.99; //No sense to have DN uncertainty > 100%
	  
	  zinv_systUp += err_ZGUp*err_ZGUp;
	  zinv_systDn += err_ZGDn*err_ZGDn;
	  
	  if(doSimultaneousFit && includeCR){
	    
	    datacard << "zinv_ZGratio_" << htName << " lnN   - " << 1.+thisErrHTUp << "/" << 1.-thisErrHTDn << " - - - -" << std::endl;
	    datacard << "zinv_ZGratio_" << jName << " lnN   - " << 1.+thisErrNJUp << "/" << 1.-thisErrNJDn << " - - - -" << std::endl;
	    datacard << "zinv_ZGratio_" << bName << " lnN   - " << 1.+thisErrNBUp << "/" << 1.-thisErrNBDn << " - - - -" << std::endl;
	    
	  }
	  else{
	    
	    datacard << "zinv_ZGratio_" << htName << " lnN   - " << 1.+thisErrHTUp << "/" << 1.-thisErrHTDn << " - -" << std::endl;
	    datacard << "zinv_ZGratio_" << jName << " lnN   - " << 1.+thisErrNJUp << "/" << 1.-thisErrNJDn << " - -" << std::endl;
	    datacard << "zinv_ZGratio_" << bName << " lnN   - " << 1.+thisErrNBUp << "/" << 1.-thisErrNBDn << " - -" << std::endl;
	    
	   }
	  
	}
	
	// If NOT using extraploation, use MC for NB>=3 and take uncorrelated uncertainty
	if( (!use_extrapolation && iR->nBJetsMin()>=2) || iR->nBJetsMin()>3 ) {
	  
	  if(includeCR)
	    datacard << "zinv_MC_" << binName << " lnN - " << 1.+err_zinv_uncorr_2b << " - - - -" << std::endl;
	  else
	    datacard << "zinv_MC_" << binName << " lnN - " << 1.+err_zinv_uncorr_2b << " - -" << std::endl;
	  
	  zinv_systUp += err_zinv_uncorr_2b*err_zinv_uncorr_2b;
	  zinv_systDn += err_zinv_uncorr_2b*err_zinv_uncorr_2b;
	  
	} else {
	  
	  if ( !use_extrapolation )
	    Ngamma = round(this_zinvCR->GetBinContent(iBin)); //If NOT using extrapolation, photon CR is corresponding bin
	  else
	    Ngamma = round(this_zinvCR->Integral()); // If using extrapolation (DEFAULT), photon CR is integrated over MT2
	  
	  // Read photon purity and corresponding uncertainty
	  Double_t x_tmp, p, p_errUp, p_errDown;
	  this_zinv_purity->GetPoint( 0, x_tmp, p );
	  p_errUp   = this_zinv_purity->GetErrorYhigh( 0 );
	  p_errDown = this_zinv_purity->GetErrorYlow ( 0 ); 

	  // Uncerainty on purity
	  if( Ngamma>0 && p>0 ) {
	    
	    if(doSimultaneousFit && includeCR)
	      datacard << "zinv_puritySyst_" << zinvCR_name << " lnN  - " << 1.+err_zinv_puritySyst << " - - - -" << std::endl;
	    else
	      datacard << "zinv_puritySyst_" << zinvCR_name << " lnN  - " << 1.+err_zinv_puritySyst << " - -" << std::endl;
	    zinv_systUp += err_zinv_puritySyst*err_zinv_puritySyst;
	    zinv_systDn += err_zinv_puritySyst*err_zinv_puritySyst;
	    
	    p = fabs(p);
	    
	    if(doSimultaneousFit && includeCR){
	      if( fabs(p_errDown)/p < 1. )  
		datacard << "zinv_purity_" << zinvCR_name << " lnN  - " << 1.+p_errUp/p << "/" << 1.-p_errDown/p << " - - - -" << std::endl;
	      else 
		datacard << "zinv_purity_" << zinvCR_name << " lnN  - " << 1.+p_errUp/p << "/0.01 - - - -" << std::endl;
	    }
	    else{
	      if( fabs(p_errDown)/p < 1. )  
		datacard << "zinv_purity_" << zinvCR_name << " lnN  - " << 1.+p_errUp/p << "/" << 1.-p_errDown/p << " - -" << std::endl;
	      else 
		datacard << "zinv_purity_" << zinvCR_name << " lnN  - " << 1.+p_errUp/p << "/0.01 - -" << std::endl;
	    }
	    zinv_systUp += (p_errUp/p)*(p_errUp/p);
	    zinv_systDn += (p_errDown/p)*(p_errDown/p);
	    
	  }
	  else if( Ngamma>0 ){ //If purity is 0 (due to missing stats), take 100% stat. uncerainty on purity
	    
	    if(doSimultaneousFit && includeCR){
	      datacard << "zinv_puritySyst_" << zinvCR_name << " lnN  - " << 1.+err_zinv_puritySyst << " - - - -" << std::endl;
	      datacard << "zinv_purity_" << zinvCR_name << " lnN  - 1.0/0.01 - - - -" << std::endl;
	    }
	    else{
	      datacard << "zinv_puritySyst_" << zinvCR_name << " lnN  - " << 1.+err_zinv_puritySyst << " - -" << std::endl;
	      datacard << "zinv_purity_" << zinvCR_name << " lnN  - 1.0/0.01 - -" << std::endl;
	    }
	    
	    zinv_systUp += 1.0;
	    zinv_systDn += 1.0;
	    
	  }
	  
          
	  // Get Z/G ratio
	  float R = fabs(this_zinv_ratio->GetBinContent(iBin));
	  float relativeErr;
	  
	  // Keep information on last non-zero ratio. Will use last non-zero ratio if ratio for one bin is zero.
	  if (R>0) lastR_zinv = R;
	  
	  if( fabs(this_zinv->GetBinContent(iBin))>0 )
	    relativeErr = err_zinv_shape / (nBins-1) * (iBin-1);
	  //	     relativeErr = err_zinv_shape / ((nBins-1) * (nBins-1)) * (iBin-1) * (iBin-1); //Parabolic shape uncertainty
	  else
	    relativeErr = 0.0;
	  
	  if( !use_extrapolation ){
	    if(includeCR)
	      datacard << "zinv_CRstat_" << std::setprecision(5) << gammaConvention( yield_zinv, Ngamma, 1, binName, binName, (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv ) << " - -" << std::setprecision(3) << std::endl;
	    else
	      datacard << "zinv_CRstat_" << std::setprecision(5) << gammaConvention( yield_zinv, Ngamma, 1, binName, binName, (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv ) << std::setprecision(3) << std::endl;
	  }
	  else {
	    if(doSimultaneousFit && includeCR)
	      datacard << "zinv_CRstat_" << std::setprecision(5) << gammaConvention( yield_zinv, Ngamma, 1, zinvCR_name, binName, (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv ) << " - -" << std::setprecision(3) << std::endl;
	    else
	      datacard << "zinv_CRstat_" << std::setprecision(5) << gammaConvention( yield_zinv, Ngamma, 1, zinvCR_name, binName, (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv ) << std::setprecision(3) << std::endl;
	    
	    if( nBins>1 ){
	      if( iBin==1 && fabs(yield_zinv)>0 ){ //For shape uncertainty on 1st MT2 bin, take 1-tot. unc. form other bins (PIVOT at second bin)
		
		if(doSimultaneousFit && includeCR)
		  datacard << "zinv_shape_" << zinvCR_name << " lnN  - " << 1.-shapeErr_zinv/fabs(yield_zinv) << " - - - -" << std::endl;
		else
		  datacard << "zinv_shape_" << zinvCR_name << " lnN  - " << 1.-shapeErr_zinv/fabs(yield_zinv) << " - - " << std::endl;
		zinv_systUp += (shapeErr_zinv/fabs(yield_zinv))*(shapeErr_zinv/fabs(yield_zinv));
		zinv_systDn += (shapeErr_zinv/fabs(yield_zinv))*(shapeErr_zinv/fabs(yield_zinv));
		
	      }
	      else{
		if(doSimultaneousFit && includeCR)
		  datacard << "zinv_shape_" << zinvCR_name << " lnN  - " << 1.+relativeErr << " - - - -" << std::endl;
		else
		  datacard << "zinv_shape_" << zinvCR_name << " lnN  - " << 1.+relativeErr << " - - " << std::endl;
		zinv_systUp += relativeErr*relativeErr;
		zinv_systDn += relativeErr*relativeErr;
	      }
	    }
	  }

	  // Get Poisson uncertainty for table
	  double yield_zinv_up, yield_zinv_dn;
	  RooHistError::instance().getPoissonInterval(Ngamma,yield_zinv_dn,yield_zinv_up,1.);
	  yield_zinv_up *= (Ngamma>0) ? yield_zinv/Ngamma : (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv;
	  yield_zinv_dn *= (Ngamma>0) ? yield_zinv/Ngamma : (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv;
	  
	  zinv_statUp = yield_zinv_up-yield_zinv;
	  zinv_statDn = yield_zinv-yield_zinv_dn;
	  
	  // Uncertainty on transfer factor
	  float alphaErr;
	  if ( R>0 ) alphaErr= this_zinv_ratio->GetBinError(iBin)/R;
	  else alphaErr=1.0;
	  if( doSimultaneousFit && includeCR )
	    datacard << "zinv_alphaErr_" << binName << " lnN  - " << 1.+alphaErr << " - - - -" << std::endl;
	  else
	    datacard << "zinv_alphaErr_" << binName << " lnN  - " << 1.+alphaErr << " - -" << std::endl;

	  zinv_systUp += alphaErr*alphaErr;
	  zinv_systDn += alphaErr*alphaErr;
	  
	  
	} // if extrapolation
	
	zinv_nCR = Ngamma; // Just for table
	
      } // if zinv
      
      
      
      
      // LOST LEPTON SYSTEMATICS:
      if( yield_llep>=0. ) {
	
	// Get TF for lost-lepton
	float Rllep = this_llep_ratio->GetBinContent(iBin);
	if( Rllep > 0 ) lastR_llep = Rllep; // Keep track of last non-zero TF. Will use it for 'next' bin if TF will be zero.
	
//	float alphaErr_llep;
//	if ( Rllep>0 ) alphaErr_llep= this_llep_ratio->GetBinError(iBin)/Rllep;
//	else alphaErr_llep=1.0;
	
	// Shape uncertainty
	float relativeErr_llep;
	 if( this_llep->GetBinContent(iBin)>0 )
	   relativeErr_llep = err_llep_shape / (nBins-1) * (iBin-1);
	 //relativeErr_llep = err_llep_shape / ((nBins-1) * (nBins-1)) * (iBin-1) * (iBin-1); // Parabolic shape extrapolation
	 else
	   relativeErr_llep = 0.0;
	 
	 // Uncertainty from lepton efficiency
	 if(doSimultaneousFit && includeCR)
	   datacard << "llep_lepeff_" << llepCR_name << "  lnN  - - " << 1.+err_llep_lepEff << " - - -" << std::endl;
	 else
	   datacard << "llep_lepeff_" << llepCR_name << "  lnN  - - " << 1.+err_llep_lepEff << " -" << std::endl;	   

	 llep_systUp += err_llep_lepEff*err_llep_lepEff;
	 llep_systDn += err_llep_lepEff*err_llep_lepEff;
	 
	 // Gamma function, or lnU for simultaneous fit
	 if(doSimultaneousFit && includeCR)
	   datacard << "llep_" << llepCR_name << "  lnU  - - 5.0 - - 5.0" << std::endl;
	 else if (doSimultaneousFit)
	   datacard << "llep_" << llepCR_name << "  lnU  - - 5.0 -" << std::endl;
	 else
	   datacard << "llep_CRstat_" << gammaConvention( yield_llep, round(N_llep_CR), 2, llepCR_name, binName, (Rllep>0) ? ( (Rllep>3) ? 2 : Rllep ) : lastR_llep ) << std::endl;  

	 // Get Poisson uncertainty for table
	 double yield_llep_up, yield_llep_dn;
	 RooHistError::instance().getPoissonInterval(round(N_llep_CR),yield_llep_dn,yield_llep_up,1.);
	 yield_llep_up *= (round(N_llep_CR)>0) ? yield_llep/round(N_llep_CR) : (Rllep>0) ? ( (Rllep>3) ? 2 : Rllep ) : lastR_llep;
	 yield_llep_dn *= (round(N_llep_CR)>0) ? yield_llep/round(N_llep_CR) : (Rllep>0) ? ( (Rllep>3) ? 2 : Rllep ) : lastR_llep;
	 llep_statUp = yield_llep_up-yield_llep;
	 llep_statDn = yield_llep-yield_llep_dn;

	 if( yield_llep>=0. ) {

	   // MC stat. uncertainty
	   float err_llep_mcstat;
	   if(yield_llep>0) err_llep_mcstat = this_llep->GetBinError(iBin)/yield_llep;
	   else err_llep_mcstat = 1.0;
	 
	   if(doSimultaneousFit && includeCR)
	     datacard << "llep_MCstat_" << binName << " lnN  - - " << 1.+err_llep_mcstat << " - - -" << std::endl;
	   else
	     datacard << "llep_MCstat_" << binName << " lnN  - - " << 1.+err_llep_mcstat << " -" << std::endl;
	   llep_systUp += err_llep_mcstat*err_llep_mcstat;
	   llep_systDn += err_llep_mcstat*err_llep_mcstat;
	     
	   // Shape uncertainty
	   if( nBins > 1 ){
	     if( iBin==1 && yield_llep>0 ){
	       
	       if(doSimultaneousFit && includeCR)
		 datacard << "llep_shape_" << llepCR_name << " lnN - - " << 1.-shapeErr_llep/yield_llep << " - - -" << std::endl;
	       else
		 datacard << "llep_shape_" << llepCR_name << " lnN - - " << 1.-shapeErr_llep/yield_llep << " - " << std::endl;
	       llep_systUp += (shapeErr_llep/yield_llep)*(shapeErr_llep/yield_llep);
	       llep_systDn += (shapeErr_llep/yield_llep)*(shapeErr_llep/yield_llep);
	     
	     }
	     else{
	       
	       if(doSimultaneousFit && includeCR)
		 datacard << "llep_shape_" << llepCR_name << " lnN - - " << 1+relativeErr_llep << " - - -" << std::endl;
	       else
		 datacard << "llep_shape_" << llepCR_name << " lnN - - " << 1+relativeErr_llep << " - " << std::endl;
	       llep_systUp += relativeErr_llep*relativeErr_llep;
	       llep_systDn += relativeErr_llep*relativeErr_llep;
	     
	     }
	   }

	   
	   // TF uncertainty
	   if( iR->htMin()==200 && iR->nJetsMin()>=7 ){ // 40% for very-low-HT and NJ>=7 (due to JECs)
	     
	     if(doSimultaneousFit && includeCR)
	       datacard << "llep_alpha_" << llepCR_name << " lnN - - " << 1.+0.4 << " - - -" << std::endl; 
	     else
	       datacard << "llep_alpha_" << llepCR_name << " lnN - - " << 1.+0.4 << " - " << std::endl;
	     llep_systUp += 0.4*0.4;
	     llep_systDn += 0.4*0.4;
	   
	   } else if( iR->nJetsMin()>=7 && iR->nBJetsMin()>=3 ){ // 15% for NJ>=7 and NB>=3 (due to b-tag SF)

	     if(doSimultaneousFit && includeCR)
	       datacard << "llep_alpha_" << llepCR_name << " lnN - - " << 1.+0.15 << " - - -" << std::endl;
	     else
	       datacard << "llep_alpha_" << llepCR_name << " lnN - - " << 1.+0.15 << " - " << std::endl;
             llep_systUp += 0.15*0.15;
             llep_systDn += 0.15*0.15;

	   } else{ // 10% elsewhere
	     
	     if(doSimultaneousFit && includeCR)
	       datacard << "llep_alpha_" << llepCR_name << " lnN - - " << 1.+0.1 << " - - -" << std::endl;
	     else
	       datacard << "llep_alpha_" << llepCR_name << " lnN - - " << 1.+0.1 << " - " << std::endl;
             llep_systUp += 0.1*0.1;
             llep_systDn += 0.1*0.1;

	   }
	   
	   
	 }
	 
	 llep_nCR = N_llep_CR; // Just for table (CR counts)
	 
       }




       // QCD SYSTEMATICS:
       int NQCD_cr;

       if( yield_qcd>=0. ) {

	 // QCD estimate for monojet
	 if( iR->nJetsMax()==1 ){
	   
	   NQCD_cr = round(this_qcdCR->GetBinContent(iBin));
	   float thisRqcd = this_qcd_ratio->GetBinContent(iBin);
	   
	   if( thisRqcd > 0 ) lastR_qcd_mono = thisRqcd;
	   
	   if(doSimultaneousFit && includeCR)
	     datacard << "qcd_CRstat_" << std::setprecision(6) << gammaConvention( yield_qcd, NQCD_cr, 3, binName, binName, (thisRqcd>0) ? thisRqcd : lastR_qcd_mono ) << " - -" << std::setprecision(3) << std::endl;
	   else
	     datacard << "qcd_CRstat_" << std::setprecision(6) << gammaConvention( yield_qcd, NQCD_cr, 3, binName, binName, (thisRqcd>0) ? thisRqcd : lastR_qcd_mono ) << std::setprecision(3) << std::endl;

	   // Get Poisson uncertainty (from CR) for table
           double yield_qcd_up, yield_qcd_dn;
	   RooHistError::instance().getPoissonInterval(NQCD_cr,yield_qcd_dn,yield_qcd_up,1.);
           yield_qcd_up *= (NQCD_cr>0.) ? yield_qcd/NQCD_cr : lastR_qcd_mono;
           yield_qcd_dn *= (NQCD_cr>0.) ? yield_qcd/NQCD_cr : lastR_qcd_mono;
           qcd_statUp = yield_qcd_up-yield_qcd;
           qcd_statDn = yield_qcd-yield_qcd_dn;

	   // Get uncertainty on TF
	   float thisRQCDErr_monojet = ( thisRqcd > 0 ) ? this_qcd_ratio->GetBinError(iBin)/thisRqcd : 1.0;
	   if(doSimultaneousFit && includeCR)
	     datacard << "qcd_alphaErr_" << binName << " lnN - - - " << 1.+thisRQCDErr_monojet << " - -" << std::endl;
	   else
	     datacard << "qcd_alphaErr_" << binName << " lnN - - - " << 1.+thisRQCDErr_monojet << std::endl;
           qcd_systUp += thisRQCDErr_monojet*thisRQCDErr_monojet;
           qcd_systDn += thisRQCDErr_monojet*thisRQCDErr_monojet;
	 
	 }
	 
	 else{ // QCD estimate for multi-jet

	   // Get events in CR
	   NQCD_cr = round(this_qcdCR->GetBinContent(iBin));
	   
	   // Get dPhi ratio
	   float thisRqcd = this_qcd_ratio->GetBinContent(iBin);
	   float thisRqcdErr = this_qcd_ratio->GetBinError(iBin);
	   
	   // Get Purity
	   float thisQCDPurity = this_qcd_purity->GetBinContent(iBin);
	   float thisQCDPurityErr = this_qcd_purity->GetBinError(iBin);
	   
	   if( thisRqcd>0 && thisRqcdErr>0) thisRqcdErr=thisRqcdErr/thisRqcd;
	   else thisRqcdErr=0.0;

	   if( thisQCDPurity>0 && thisQCDPurityErr>0) thisQCDPurityErr=thisQCDPurityErr/thisQCDPurity;
	   else thisQCDPurityErr=0.0;
	   
	   // Get uncertainty on dPhi ratio from fit window variation
	   float thisRFitVarErr = this_qcd_ratioSystFit->GetBinContent(iBin);
	   
	   // Get F(J) and R(B)
	   int thisFJetsBin = this_qcd_fjets->FindBin( iR->nJetsMin() );
	   int thisRBBin = this_qcd_rb->FindBin( iR->nBJetsMin() );
	   
	   float thisFJets = this_qcd_fjets->GetBinContent(thisFJetsBin);
	   float thisRB = this_qcd_rb->GetBinContent(thisRBBin);
	   
	   float thisFJetsErr = this_qcd_fjets->GetBinError(thisFJetsBin);
	   float thisRBErr = this_qcd_rb->GetBinError(thisRBBin);
	   
	   // For TR 2-6j, 3b get F(J) as 2-3j + 4-6j
	   if( iR->nBJetsMin()==3 && iR->nJetsMin()==2 ) {
	     
	     int otherFJetsBin  = this_qcd_fjets->FindBin(4);
	     thisFJets += this_qcd_fjets->GetBinContent( otherFJetsBin );
	     thisFJetsErr *= thisFJetsErr;
	     thisFJetsErr += this_qcd_fjets->GetBinError( otherFJetsBin )*this_qcd_fjets->GetBinError( otherFJetsBin );
	     thisFJetsErr  = sqrt(thisFJetsErr);

	   }

	   if( thisFJets>0 && thisFJetsErr>0 ) thisFJetsErr=thisFJetsErr/thisFJets;
	   else thisFJetsErr=0.0;
	   
	   if( thisRB>0  && thisRBErr>0 ) thisRBErr=thisRBErr/thisRB;
	   else thisRBErr=0.0;
	   
	   float thisFractionsErr = TMath::Sqrt(thisRBErr*thisRBErr+thisFJetsErr*thisFJetsErr+thisQCDPurityErr*thisQCDPurityErr);
	   float thisFractions = thisFJets*thisRB;
	   if( thisQCDPurity > 0 ) thisFractions *= thisQCDPurity;

	   if(doSimultaneousFit && includeCR)
	     datacard << "qcd_CRstat_" << std::setprecision(6) << gammaConvention( yield_qcd, NQCD_cr, 3, binName, binName, (thisRqcd*thisFractions>0) ? thisRqcd*thisFractions : 1.0 ) << " - -" << std::setprecision(3) << std::endl;
	   else
	     datacard << "qcd_CRstat_" << std::setprecision(6) << gammaConvention( yield_qcd, NQCD_cr, 3, binName, binName, (thisRqcd*thisFractions>0) ? thisRqcd*thisFractions : 1.0 ) << std::setprecision(3) << std::endl;   

	   // Get Poisson uncertainty from CR for tables
	   double yield_qcd_up, yield_qcd_dn;
	   RooHistError::instance().getPoissonInterval(NQCD_cr,yield_qcd_dn,yield_qcd_up,1.);
	   yield_qcd_up *= (NQCD_cr>0.) ? yield_qcd/NQCD_cr : thisRqcd*thisFJets*thisRB;
	   yield_qcd_dn *= (NQCD_cr>0.) ? yield_qcd/NQCD_cr : thisRqcd*thisFJets*thisRB;
	   qcd_statUp = yield_qcd_up-yield_qcd;
	   qcd_statDn = yield_qcd-yield_qcd_dn;

	   // Uncertainty on F(J) and R(B)
	   if(doSimultaneousFit && includeCR)
	     datacard << "qcd_FJRBsyst_" << binName << " lnN - - - " <<  1.+thisFractionsErr  << " - -" << std::endl;
	   else
	     datacard << "qcd_FJRBsyst_" << binName << " lnN - - - " <<  1.+thisFractionsErr  << std::endl;
	   qcd_systUp += thisFractionsErr*thisFractionsErr;
           qcd_systDn += thisFractionsErr*thisFractionsErr;

	   // Stat. uncertainty on r(dPhi)
	   if(doSimultaneousFit && includeCR)
	     datacard << "qcd_RPHIstat_" << htName << " lnN - - - " <<  1.+thisRqcdErr  << " - -" << std::endl;
	   else
	     datacard << "qcd_RPHIstat_" << htName << " lnN - - - " <<  1.+thisRqcdErr  << std::endl;
	   qcd_systUp += thisRqcdErr*thisRqcdErr;
	   qcd_systDn += thisRqcdErr*thisRqcdErr;
	 
	   // Syst. ucnertainty on r(dPhi)
	   if(doSimultaneousFit && includeCR)
	     datacard << "qcd_RPHIsyst_" << htName << " lnN - - - " <<  1.+thisRFitVarErr  << " - -" << std::endl;
	   else
	     datacard << "qcd_RPHIsyst_" << htName << " lnN - - - " <<  1.+thisRFitVarErr  << std::endl;
	   qcd_systUp += thisRFitVarErr*thisRFitVarErr;
	   qcd_systDn += thisRFitVarErr*thisRFitVarErr;

	 }
	
	 qcd_nCR = NQCD_cr; // Events in CR for table
	 
       }


       datacard.close();

       std::cout << "-> Created template datacard: " << datacardName << std::endl;



       // Make absolute uncertainties for table
       zinv_systUp = yield_zinv*sqrt(zinv_systUp);
       zinv_systDn = yield_zinv*sqrt(zinv_systDn);

       llep_systUp = yield_llep*sqrt(llep_systUp);
       llep_systDn = yield_llep*sqrt(llep_systDn);

       qcd_systUp = yield_qcd*sqrt(qcd_systUp);
       qcd_systDn = yield_qcd*sqrt(qcd_systDn);

       // Print the table:
       table << "### bg_name yield statUp statDown systUp systDown" << std::endl;
       table << "zinv " << yield_zinv << " " << zinv_statUp << " " << zinv_statDn << "  " << zinv_systUp << " " << zinv_systDn << std::endl;
       table << "llep " << yield_llep << " " << llep_statUp << " " << llep_statDn << "  " << llep_systUp << " " << llep_systDn << std::endl;
       table << "qcd  " << yield_qcd << " " << qcd_statUp << " " << qcd_statDn << "  " << qcd_systUp << " " << qcd_systDn << std::endl;
       table << "data " << std::setprecision(6) << this_data->GetBinContent(iBin) << std::setprecision(3) <<std::endl;
       table << "zinv_nCR " << std::setprecision(6) << zinv_nCR << std::setprecision(3) <<std::endl;
       table << "llep_nCR " << std::setprecision(6) << llep_nCR << std::setprecision(3) <<std::endl;
       table << "qcd_nCR " << std::setprecision(6) << qcd_nCR << std::setprecision(3)   <<std::endl;
       table.close();

       std::cout << "-> Created BG table: " << tableName << std::endl;

       
     }// for bins
     
  } // for regions
  
  
  
  // now create datacards for all signals
  //////  std::vector<MT2Analysis<MT2Estimate>*> signals = MT2Analysis<MT2Estimate>::readAllFromFile( mc_fileName, "SMS" );
  //////  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( mc_fileName, "SMS", "" );

  std::vector<MT2Analysis<MT2EstimateSigContSyst>*> signals;
  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_isr;
  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagHeavy;
  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagLight;
  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_lepEff;

  std::string modelName = model;
  if( model == "T2tt" || model == "T1tttt" )
    modelName += "_sigcontam";

  signals       = MT2Analysis<MT2EstimateSigContSyst>::readAllSystFromFile( "./signalScansFromDominick/"+modelName+"_eth.root", modelName, "isr" );

  if( includeSignalUnc ){
    signals_isr       = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "./signalScansFromDominick/"+modelName+"_eth.root", modelName, "isr" );
    signals_bTagHeavy = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "./signalScansFromDominick/"+modelName+"_eth.root", modelName, "btagsf_heavy" );
    signals_bTagLight = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "./signalScansFromDominick/"+modelName+"_eth.root", modelName, "btagsf_light" );
    
    if( model == "T2tt" || model == "T1tttt" )
      signals_lepEff = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "./signalScansFromDominick/"+modelName+"_eth.root", modelName, "lepeff" );
  }
  
  
  ////// To replace signal to existing one (for T2tt corridor studies)
  //std::vector<MT2Analysis<MT2Estimate>*> signalVeto = MT2Analysis<MT2Estimate>::readAllFromFile( dir+"/T2tt_175_0_evtVeto.root", "T2tt_175_0" );
  //std::vector<MT2Analysis<MT2Estimate>*> signalVeto_1lCR = MT2Analysis<MT2Estimate>::readAllFromFile( dir+"/T2tt_175_0_evtVeto_1LCR.root", "llepCR" );
  //std::vector<MT2Analysis<MT2Estimate>*> signalVeto = MT2Analysis<MT2Estimate>::readAllFromFile( dir+"/t2tt.root", "T2tt" );
  //std::vector<MT2Analysis<MT2Estimate>*> signalVeto_1lCR = MT2Analysis<MT2Estimate>::readAllFromFile( dir+"/llepControlRegion/t2tt.root", "llepCR" );
  

  for( unsigned  isig=0; isig<signals.size(); ++isig ) {

    
    // signals[isig]           ->setName(model.c_str());
    // if(includeSignalUnc){
    //   signals_isr[isig]       ->setName(model.c_str());
    //   signals_bTagHeavy[isig] ->setName(model.c_str());
    //   signals_bTagLight[isig] ->setName(model.c_str());
    //   // if( model == "T2tt" || model == "T1tttt" )
    //   // 	signals_lepEff[isig]    ->setName(model.c_str());
    // }

    // Name convention
    std::string sigName;
    sigName = signals[isig]->getName();
    //sigName = getSimpleSignalName( signals[isig]->getName() );

    std::string scont = "_sigcontam";
    std::string::size_type pos = sigName.find(scont);
    if(pos != std::string::npos) sigName.erase(pos,scont.length());

    
    // Local path for datacards
    std::string path = dir + "/datacards_" + sigName;
    system(Form("mkdir -p %s", path.c_str()));

    // SE path for datacards
    std::string pathSE = "";
    if (label=="")
      pathSE = dir + "/datacards_" + sigName;
    else
      pathSE = dir + "/datacards_" + sigName + "_" + label;
    
    std::string path_mass = path;

    float xs_norm=1.;
    ////// If you need to renormalize T2qq xsec, uncomment
//    if( signals[isig]->getName().find("T2qq") != std::string::npos ) 
//      xs_norm=8./10.;

    // Start loop over topological regions
    for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
      
      TH3D* this_signal3d_central;
      TH1D* this_signalParent;
      
      // Read signal analysis, and take 3D histogram if it exists for this region
      MT2EstimateSigContSyst* thisSigSystCentral = signals[isig]->get(*iR);
      //MT2Estimate* thisSigSystCentral = signalVeto[isig]->get(*iR);
      if( thisSigSystCentral->yield3d!=0 ){
	
	//this_signal3d_central = signalVeto[isig]->get(*iR)->yield3d;
	this_signal3d_central        = signals[isig]->get(*iR)->yield3d;
	
      }
      else continue;
      

      // Project SUSY parent mass on 1D histogram (to be used to loop over scan masses)
      this_signalParent = this_signal3d_central->ProjectionY("mParent");
      
      // Initialize histograms for signal systematics
      TH3D* this_signal3d_bTagHeavy_Up;
      TH3D* this_signal3d_bTagLight_Up;
      TH3D* this_signal3d_isr_Up;
      TH3D* this_signal3d_lepEff_Up;
      
      // Read signal systematic analysis, and take 3D histogrms if they exist for this region
      MT2EstimateSigSyst* thisSigSyst_isr;
      MT2EstimateSigSyst* thisSigSyst_bTagHeavy;
      MT2EstimateSigSyst* thisSigSyst_bTagLight;
      MT2EstimateSigSyst* thisSigSyst_lepEff;

      if( includeSignalUnc ){
	
	thisSigSyst_isr = signals_isr[isig]->get(*iR);
	if( thisSigSyst_isr->yield3d_systUp!=0 )
	  this_signal3d_isr_Up       = (TH3D*) signals_isr[isig]->get(*iR)->yield3d_systUp->Clone();
	else
	  this_signal3d_isr_Up       = (TH3D*) signals[isig]->get(*iR)->yield3d->Clone();

	thisSigSyst_bTagHeavy = signals_bTagHeavy[isig]->get(*iR);
	if( thisSigSyst_bTagHeavy->yield3d_systUp!=0 )
	  this_signal3d_bTagHeavy_Up       = (TH3D*) signals_bTagHeavy[isig]->get(*iR)->yield3d_systUp->Clone();
	else
	  this_signal3d_bTagHeavy_Up = (TH3D*) signals[isig]->get(*iR)->yield3d->Clone(); 

	thisSigSyst_bTagLight = signals_bTagLight[isig]->get(*iR);
	if( thisSigSyst_bTagLight->yield3d_systUp!=0 )
	  this_signal3d_bTagLight_Up       = (TH3D*) signals_bTagLight[isig]->get(*iR)->yield3d_systUp->Clone();
	else
	  this_signal3d_bTagLight_Up = (TH3D*) signals[isig]->get(*iR)->yield3d->Clone();
	
	if( model == "T2tt" || model == "T1tttt" ){
	  thisSigSyst_lepEff = signals_lepEff[isig]->get(*iR);
	  if( thisSigSyst_lepEff->yield3d_systUp!=0 )
	    this_signal3d_lepEff_Up       = (TH3D*) signals_lepEff[isig]->get(*iR)->yield3d_systUp->Clone();
	  else
	    this_signal3d_lepEff_Up = (TH3D*) signals[isig]->get(*iR)->yield3d->Clone();
	}
	
      }
      
      // Start loop over SUSY parent mass
      for( int iBinY=1; iBinY<this_signalParent->GetNbinsX()+1; ++iBinY ){
	
	float mParent = this_signalParent->GetBinLowEdge(iBinY);
	if( !(mParent >= m1-1 && mParent < m2-1) ) continue;
	
	// Get LSP masses for this parent mass
	TH1D* this_signalLSP = this_signal3d_central->ProjectionZ("mLSP", 0, -1, iBinY, iBinY);
	
	// Start loop over LSP masses
	for( int iBinZ=1; iBinZ < iBinY; ++iBinZ ) {
	  
	  float mLSP = this_signalLSP->GetBinLowEdge(iBinZ);
	  if( !(mLSP >= m11-1 && mLSP < m22-1) ) continue;
	  
	  // Get MT2 yield histogram for this mass point
	  TH1D* this_signal      = this_signal3d_central->ProjectionX("mt2"     , iBinY, iBinY, iBinZ, iBinZ);
	  TH1D* this_signal_syst = this_signal3d_central->ProjectionX("mt2_syst", iBinY, iBinY, iBinZ, iBinZ);
	  
	  if (doGenAverage) {
	    TH3D* this_signal3d_central_genmet = signals[isig]->get(*iR)->yield3d_genmet;
	    TH1D* this_signal_genmet = this_signal3d_central_genmet->ProjectionX("mt2_genmet", iBinY, iBinY, iBinZ, iBinZ);
	    this_signal->Add(this_signal_genmet);
	    this_signal->Scale(0.5);

	    this_signal_syst->Add(this_signal_genmet, -1.0);
	    this_signal_syst->Scale(0.5); // half difference between gen and reco

	  }
	  

//	  // If want to replace central yield
//	  TH1D* this_signal_veto = (TH1D*) signalVeto[isig]->get(*iR)->yield->Clone();
	  
	  if( this_signal->Integral() <=0 ) continue;
	  
	  ////// Signal contamination
	  TH1D* this_signalContamination;
	  TH1D* this_signalContamination_syst;
	 
	  TH3D* this_signal3d_crsl;
	  TH1D* this_signal_crsl;
	  TH1D* this_signal_crsl_syst;
	  TH1D* this_signal_alpha;
	  
	  if( (model == "T2tt" || model == "T1tttt") && doSignalContamination ){
	    this_signal3d_crsl        = (TH3D*) signals[isig]->get(*iR)->yield3d_crsl       ->Clone();
	    //this_signal3d_crsl = (TH3D*) signalVeto_1lCR[isig]->get(*iR)->yield3d->Clone();
	    this_signal_crsl        = this_signal3d_crsl       ->ProjectionX("mt2_crsl", iBinY, iBinY, iBinZ, iBinZ);
	    this_signal_alpha  = (TH1D*) signals[isig]->get(*iR)->yield_alpha->Clone();

	    if (doGenAverage){
	      TH3D *this_signal3d_crsl_genmet = (TH3D*) signals[isig]->get(*iR)->yield3d_crsl_genmet->Clone();
	      TH1D *this_signal_crsl_genmet = this_signal3d_crsl_genmet->ProjectionX("mt2_crsl", iBinY, iBinY, iBinZ, iBinZ);
	      this_signal_crsl->Add(this_signal_crsl_genmet);
	      this_signal_crsl->Scale(0.5);

	      this_signal_crsl_syst= this_signal3d_crsl->ProjectionX("mt2_crsl_syst", iBinY, iBinY, iBinZ, iBinZ);
	      this_signal_crsl_syst->Add(this_signal_crsl_genmet,-1.0);
	      this_signal_crsl_syst->Scale(0.5);
	    }
	    
	  //	  //If want to replace yield in 1l CR for signal contamination
	  //	  TH1D* this_signal_veto_1lCR = (TH1D*) signalVeto_1lCR[isig]->get(*iR)->yield->Clone();
	  //	  // Then, below replace this_signal_crsl with this histogram;
	    
	    if( doSimultaneousFit )
	      this_signalContamination = (TH1D*) this_signal_crsl->Clone();
	    else{
	      this_signalContamination = (TH1D*) this_signal_alpha->Clone();
	      this_signalContamination->Scale( this_signal_crsl->Integral() );
	      if(doGenAverage){
		this_signalContamination_syst = (TH1D*) this_signal_alpha->Clone("mt2_cont_syst");
		this_signalContamination_syst->Scale( this_signal_crsl_syst->Integral() );
	      }
	    }
	  }
	  else{
	    
	    this_signalContamination = (TH1D*) this_signal->Clone();
	    this_signalContamination->Scale(0.);
	    if(doGenAverage){
	      this_signalContamination_syst = (TH1D*) this_signal_syst->Clone("mt2_cont_syst");
	      this_signalContamination_syst->Scale(0.);
	    }

	  }
	    
	  //////
	  
	  //Start loop over MT2 bins
	  for( int iBin=1; iBin<this_signal->GetNbinsX()+1; ++iBin ) {
	    
	    bool includeCR=false;
	    if(iBin==1) includeCR=true;
	    if(iR->nJetsMin()>=7 && iR->nBJetsMin()>1) includeCR=false;
	    
	    if( this_signal->GetBinLowEdge( iBin ) > iR->htMax() && iR->htMax()>0 ) continue;
	    
	    float mt2Min = this_signal->GetBinLowEdge( iBin );
	    float mt2Max = (iBin==this_signal->GetNbinsX()) ?  -1. : this_signal->GetBinLowEdge( iBin+1 );
	    
	    // If bin is empty, do not create card
	    if( this_signal->GetBinContent(iBin) <=0 );
	    else{
	      
	      std::string binName;
	      if( mt2Max>=0. )
		binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
	      else
		binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );
	      
	      // If datacard exists already on SE, do not create it again
	      Long_t id;
	      Long_t flags; 
	      Long_t modtime;
	      Long_t size;
	      std::string fullPathSE;
	      int checkFileSE;
	      std::string rmOnSE;
	      if( copy2SE ){
		fullPathSE = Form("/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/%s/datacards_%.0f_%.0f/datacard_%s_%s_%.0f_%.0f.txt", pathSE.c_str(), mParent, mLSP, binName.c_str(), sigName.c_str(), mParent, mLSP);
		checkFileSE = (int) gSystem->GetPathInfo(fullPathSE.c_str(), &id, &size, &flags, &modtime);
		std::cout << fullPathSE << "\t" << checkFileSE << "\t" <<size<< std::endl;
		
		//std::string rmOnSE( Form("env --unset=LD_LIBRARY_PATH gfal-rm srm://t3se01.psi.ch/%s", fullPathSE.c_str()) );
		rmOnSE = Form("gfal-rm srm://t3se01.psi.ch/%s", fullPathSE.c_str()) ;
		
		if( checkFileSE==0 && (size)==0 ){
		  
		  std::cout << "Removing. File " << fullPathSE << " exists and has zero size " << (size) << ". Removing." << std::endl;
		  system( rmOnSE.c_str() );
		  
		}
		else if ( checkFileSE==0 && (size)>0 ){
		  
		  std::cout << "Skipping. File " << fullPathSE << " exists and has non-zero size  " << (size) << ". Skipping." << std::endl;
		  
		  continue;
		
		}
	      }
	      
	      // Get template card for this bin
	      std::string templateDatacard( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
	      
	      // Create new card for this bin
	      std::string newDatacard( Form("%s/datacard_%s_%s_%.0f_%.0f.txt", path_mass.c_str(), binName.c_str(), sigName.c_str(), mParent, mLSP) );
	      std::string helpDatacard( Form("%s/datacard_%s_%s_%.0f_%.0f_forSed.txt", path_mass.c_str(), binName.c_str(), sigName.c_str(), mParent, mLSP) );
	      
	      std::ifstream thisNewDatacard( newDatacard.c_str() );
	      if( thisNewDatacard.good() ) continue;

	      float sig = this_signal->GetBinContent(iBin);
	      float sigErr = this_signal->GetBinError(iBin)/sig;
	      float sig_syst = 0;
	      
	      ////// If you want to replace central yield
	      //	      float sig_veto = this_signal_veto->GetBinContent(iBin);
	      //	      float sigErr_veto = this_signal_veto->GetBinError(iBin)/sig_veto;

	      sig*=xs_norm; // To eventually rescale xsec.
	      
	      // Siganl Contamination
	      double sigContErr = 0.0;
	      double sigCont = 0.0;
	      if(doSignalContamination && doSimultaneousFit)
		sigCont = this_signalContamination->IntegralAndError(1, -1, sigContErr);
	      else if(doSignalContamination && !doSimultaneousFit){
		sigCont = this_signalContamination->GetBinContent(iBin);
		sigContErr = this_signalContamination->GetBinError(iBin);
	      }
	      sigContErr = (sigCont > 0) ? fabs(sigContErr)/sigCont : 0.0;
	      //
	      
	      
	      float isrErr;
	      float bTagErr_heavy;
	      float bTagErr_light;
	      float lepEffErr;
	      
	      if( includeSignalUnc ) {
		
		isrErr = this_signal3d_isr_Up->GetBinContent(iBin, iBinY, iBinZ);
		//isrErr = 2 - isrErr/sig;
		isrErr = isrErr/sig;
	      
		bTagErr_heavy = this_signal3d_bTagHeavy_Up->GetBinContent(iBin, iBinY, iBinZ);
		bTagErr_heavy = bTagErr_heavy/sig;
	      
		bTagErr_light = this_signal3d_bTagLight_Up->GetBinContent(iBin, iBinY, iBinZ);
		bTagErr_light = bTagErr_light/sig;
		
		lepEffErr = this_signal3d_lepEff_Up->GetBinContent(iBin, iBinY, iBinZ);
		lepEffErr = lepEffErr/sig;
		
	      }
	      
	      float totUncorrErr = 1.+sqrt(sigErr*sigErr+2*0.05*0.05+0.1*0.1); // MC stat + scales (5%) + JEC (10%)
	      float totUncorrErrCont = 1.+sqrt(sigContErr*sigContErr+2*0.05*0.05+0.1*0.1); // MC stat + scales (5%) + JEC (10%)

	      if(doSignalContamination && !doSimultaneousFit) sig=sig-sigCont;
	      else if(!doSignalContamination) {
		sigCont=0.;
		totUncorrErrCont=0.;
	      }
	      if(sig<0.) sig=0.;

	      if(doGenAverage){
		sig_syst = 1 + fabs((this_signal_syst->GetBinContent(iBin)-this_signalContamination_syst->GetBinContent(iBin))/(sig !=0 ? sig : 1.0)); // cont_syst=0 if no doSignalCont
		if ( (this_signal_syst->GetBinContent(iBin)-this_signalContamination_syst->GetBinContent(iBin))*sig < 0 )
		  sig_syst = 1/sig_syst; // to account for negative variation
	      }

	      
	      std::string mvCommand( Form("mv %s %s", newDatacard.c_str(), helpDatacard.c_str()) );
	      std::string rmCommand( Form("rm -f %s", helpDatacard.c_str()) );
	      

	      std::string sedCommand( Form("sed 's/XXX/%.3f/' %s > %s", sig, templateDatacard.c_str(), newDatacard.c_str()) );
	      system( sedCommand.c_str() );
	      
	      std::string sedCommand_sigCont( Form("sed -i 's/YYY/%.3f/' %s", sigCont, newDatacard.c_str()) );
	      if(doSimultaneousFit && includeCR)
		system( sedCommand_sigCont.c_str() );
	      
	      std::string sedCommand_uncErr( Form("sed -i 's/UUU/%.3f/' %s", totUncorrErr, newDatacard.c_str()) );
	      system( sedCommand_uncErr.c_str() );
	      
	      std::string sedCommand_uncErrCR( Form("sed -i 's/VVV/%.3f/' %s", totUncorrErrCont, newDatacard.c_str()) );
	      if(doSimultaneousFit && includeCR)
		system( sedCommand_uncErrCR.c_str() );
	      
	      std::string sedCommand_isrErr( Form("sed -i 's/III/%.3f/' %s", isrErr, newDatacard.c_str()) );
	      std::string sedCommand_bTagHErr( Form("sed -i 's/HHH/%.3f/' %s", bTagErr_heavy, newDatacard.c_str()) );
	      std::string sedCommand_bTagLErr( Form("sed -i 's/LLL/%.3f/' %s", bTagErr_light, newDatacard.c_str()) );
	      std::string sedCommand_lepEffErr( Form("sed -i 's/EEE/%.3f/' %s", lepEffErr, newDatacard.c_str()) );

	      std::string sedCommand_genErr( Form("sed -i 's/SSS/%.3f/' %s", sig_syst, newDatacard.c_str()) );

	      if (doGenAverage)
		system( sedCommand_genErr.c_str() );


	      if( includeSignalUnc ){
		
		system( sedCommand_isrErr.c_str() );
		system( sedCommand_bTagHErr.c_str() );
		system( sedCommand_bTagLErr.c_str() );
		
		if( model == "T2tt" || model == "T1tttt" )
		  system( sedCommand_lepEffErr.c_str() );
	      
	      }
	      
	      if( copy2SE ){
		// Copying on SE
		std::string mkdirOnSE( Form("env --unset=LD_LIBRARY_PATH gfal-mkdir -p srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/%s/datacards_%.0f_%.0f", pathSE.c_str(), mParent, mLSP) );
		std::string copyOnSE( Form("xrdcp -v %s root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/`whoami`/%s/datacards_%.0f_%.0f/datacard_%s_%s_%.0f_%.0f.txt", newDatacard.c_str(), pathSE.c_str(), mParent, mLSP, binName.c_str(), sigName.c_str(), mParent, mLSP) );
		system( mkdirOnSE.c_str() );
		system( copyOnSE.c_str() );
		
		// Attempt copying 3 times (to maximize efficiency)
		for(int c=0; c<3; ++c){
		  
		  
		  checkFileSE = (int) gSystem->GetPathInfo(fullPathSE.c_str(),&id, &size, &flags, &modtime);
		  
		  if( checkFileSE==0 && (size)==0 ){
		    
		    std::cout << "Copy did not work. Trying again: " << c << std::endl;
		    
		    system( rmOnSE.c_str() );
		    system( copyOnSE.c_str() );
		    
		  }
		  else{
		    
		    std::cout << "Copy succeded. Exiting." << std::endl;
		    
		    system( rmCommand.c_str() );
		    break;
		    
		  }
		}
	      }
	      
	      
	    }
	    
	  } // for bins X (MT2)
	} // for bins Z (mLSP)
      }// for bins Y (mParent)      
    } // for regions
    
//////    For simple MT2Estimate signal (no systematics), and one only mass point

//      for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
//	
//	TH1D* this_signal = signals[isig]->get(*iR)->yield;
//	
//	//	if( this_signal->Integral() < 0.01 ) continue; 
//
//	for( int iBin=1; iBin<this_signal->GetNbinsX()+1; ++iBin ) {
//	  
//	  if( this_signal->GetBinLowEdge( iBin ) > iR->htMax() && iR->htMax()>0 ) continue;
//	  
//	  float mt2Min = this_signal->GetBinLowEdge( iBin );
//	  float mt2Max = (iBin==this_signal->GetNbinsX()) ?  -1. : this_signal->GetBinLowEdge( iBin+1 );
//	  
//	  if( this_signal->GetBinContent(iBin) < 0 );
//	  //if( this_signal->GetBinContent(iBin) < 0.01 );
//	  else{
//	    
//	    std::string binName;
//	    if( mt2Max>=0. )
//	      binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
//	    else
//	      binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );
//	    
//	    std::string templateDatacard( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
//	    
//	    std::string newDatacard( Form("%s/datacard_%s_%s.txt", path.c_str(), binName.c_str(), sigName.c_str()) );
//	    
//	    float sig = this_signal->GetBinContent(iBin);
//	    sig*=xs_norm;
//
//	    std::string sedCommand( Form("sed 's/XXX/%.3f/g' %s > %s", sig, templateDatacard.c_str(), newDatacard.c_str()) );
//	    system( sedCommand.c_str() );
//	    
//	  }
//	  
//	} // for bins X (MT2)
//      } // for regions
    
    std::cout << "-> Created datacards in " << path_mass << std::endl;
       
  } // for signals

  return 0;

} 



std::string getSimpleSignalName( const std::string& longName ) {

  TString longName_tstr(longName);

  longName_tstr.ReplaceAll( "_", " " );
  longName_tstr.ReplaceAll( "mStop", " " );
  longName_tstr.ReplaceAll( "mGluino", " " );
  longName_tstr.ReplaceAll( "mLSP", " " );

  std::istringstream iss(longName_tstr.Data());
  std::vector<std::string> parts;
  do {
    std::string sub;
    iss >> sub;
    parts.push_back(sub);
  } while (iss);

  // parts should be:
  // [0]: SMS
  // [1]: model
  // [2]: parent mass
  // [3]: lsp mass


  std::string simpleName = parts[1] + "_" + parts[2] + "_" + parts[3];

  return simpleName;

}


std::string gammaConvention( float yieldSR, int yieldCR, int position, const std::string& corrName, const std::string& uncorrName, float testAlpha ) {
  
  std::string use_uncorrName(uncorrName);
  if( uncorrName=="" ) 
    use_uncorrName = corrName;

  std::stringstream line;
  line << std::fixed;
  line << std::setprecision(3);

  int precision = 3;
  float syst = -1.;
  if( yieldCR==0 && yieldSR==0. ) {
    line << corrName << "  gmN " << yieldCR << "   ";
    syst = testAlpha;
  } else if( yieldCR==0 && yieldSR>0. ) {
    line << use_uncorrName << "  lnN  ";
    syst = 2.;
  } else if( yieldCR>0 && yieldSR==0. ) {
    line << use_uncorrName << "  gmN 0  ";
    syst = testAlpha;
  } else {
    float alpha = yieldSR/((float)yieldCR);
    line << corrName << "  gmN " << yieldCR << "   ";
    syst = alpha;
    precision = 5;
  }
  line << std::setprecision(precision);

  for( int i=0; i<position; ++i )
    line << " - ";

  line << syst;

  for( int i=position+1; i<4; ++i )
    line << " - ";

  std::string line_str = line.str();
  return line_str;

}
