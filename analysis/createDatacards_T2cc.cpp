#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip> 

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


bool use_extrapolation = true;

int round(float d) {
  return (int)(floor(d + 0.5));
}


void writeToTemplateFile( TFile* file, MT2Analysis<MT2Estimate>* analysis, float err_uncorr );
void writeToTemplateFile_poisson( TFile* file, MT2Analysis<MT2Estimate>* analysis, const std::string& name="stat" );
MT2Analysis<MT2Estimate>* get( const std::string& name, std::vector< MT2Analysis<MT2Estimate>* > analyses, const std::string& name1, const std::string& name2="", const std::string& name3="", const std::string& name4="" );
std::string getSimpleSignalName( const std::string& longName );
std::string gammaConvention( float yieldSR, int yieldCR, int position, const std::string& corrName, const std::string& uncorrName="", float testAlpha=1. );
void getQCDestimate( float htMin, float mt2Min, float mt2Max, float nB, int& NQCD_cr, float& r );


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


  if( argc != 2 && argc != 6 ) {
    std::cout << "USAGE: ./createDatacards [configFileName] [m1] [m2] [m11] [m22]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  float m1=0.;
  float m2=2000.;
  
  float m11=0.;
  float m22=2000;

  if( argc == 6 ){
    
    m1  = std::stof( argv[2] );
    m2  = std::stof( argv[3] );

    m11  = std::stof( argv[4] );
    m22  = std::stof( argv[5] );

  }
  
  std::cout << "Will only produce datacards for parent mass and LSP mass " << m1 <<  " and " << m11 << std::endl;

  std::string dir = cfg.getEventYieldDir();
  std::string mc_fileName = dir + "/analyses.root";
  
  std::string data_fileName = dir + "/analyses.root";


  bool useMC_qcd  = false;
  bool useMC_zinv = false;
  bool useMC_llep = false;

  float err_qcd_corr    = 0.0;
  float err_qcd_uncorr  = 1.0; // 100% of QCD MC yield

  float err_llep_corr   = 0.;
  float err_llep_shape = 0.20;
  float err_llep_alpha = 0.10;
  float err_llep_lepEff = 0.07;

  float err_zinv_corr   = 0.21; // 20% on Z/gamma ratio plus added in quadrature syst on templates (2%) and on f (4%) and MC stat on Rzg (5%) -> sqrt( 20*20 + 2*2 + 4*4 +5*5 ) = 21
  float err_zinv_shape = 0.40;
  float err_zinv_alpha = 0.05;
  float err_zinv_uncorr = -1.; // will take histogram bin error
  float err_zinv_alpha_extra  = 0.2; // 20% extra uncertainty on alpha if using lower MT2 as CR
  float err_zinv_uncorr_2b = 1.0;

  float err_zinv_puritySyst = 0.1; // 10%, including 5% on purity + 8% on fragmentation
  float err_zinv_doubleRatioOffset = 0.11; // 7%, fully correlated, on zinv
  
  float err_lumi_corr   = 0.046;

  float err_sig_corr    = 0.1;
  float err_sig_uncorr  = 0.;

  MT2Analysis<MT2Estimate>* data  = MT2Analysis<MT2Estimate>::readFromFile( data_fileName, "data" );
  
  MT2Analysis<MT2Estimate>* qcd;
  MT2Analysis<MT2Estimate>* qcd_mc;
  MT2Analysis<MT2Estimate>* qcdCR;
  MT2Analysis<MT2Estimate>* qcd_ratio;
  MT2Analysis<MT2Estimate>* qcd_purity;
  MT2Analysis<MT2Estimate>* qcd_fjets;
  MT2Analysis<MT2Estimate>* qcd_fjets_vlht;
  MT2Analysis<MT2Estimate>* qcd_rb;
  MT2Analysis<MT2Estimate>* qcd_ratioSystFit;

  MT2Analysis<MT2Estimate>* qcd_monojet;
  MT2Analysis<MT2Estimate>* qcdCR_monojet;
  MT2Analysis<MT2Estimate>* qcd_ratio_monojet;


  if( useMC_qcd )
    qcd = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "QCD"  );
  else{

    qcd = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateData.root", "qcdEstimate" );
    qcdCR = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateData.root", "nCR" );
    qcd_ratio = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateData.root", "r_effective" );
    qcd_purity = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateData.root", "qcdPurity" );
    qcd_fjets = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateData.root", "f_jets_data" );
    qcd_fjets_vlht = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateData.root", "f_jets_data_noPS" );
    qcd_rb = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateData.root", "r_hat_data" );
    qcd_ratioSystFit = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateData.root", "r_systFit" );
    qcd_mc = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "QCD"  );

    qcd_monojet = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateMonojet.root", "monojet_qcdEstimate" );
    qcdCR_monojet = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateMonojet.root", "monojet_nCR" );
    qcd_ratio_monojet = MT2Analysis<MT2Estimate>::readFromFile( dir + "/qcdEstimateMonojet.root", "monojet_r" );
    
  }


  
  MT2Analysis<MT2Estimate>* zinv;
  MT2Analysis<MT2Estimate>* zinvCR;
  MT2Analysis<MT2Estimate>* zinv_ratio;
  MT2Analysis<MT2EstimateSyst>* purity;


  MT2Analysis<MT2EstimateSyst>* zllG_ht;
  MT2Analysis<MT2EstimateSyst>* zllG_nJets;
  MT2Analysis<MT2EstimateSyst>* zllG_nBJets;
  MT2Analysis<MT2EstimateSyst>* zllG_ht_monojet;

  if( useMC_zinv )
    zinv = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "ZJets");
  else {
    
    zinvCR      = MT2Analysis<MT2Estimate>    ::readFromFile( dir + "/gammaControlRegion/data.root", "gammaCR");
    
    zinv        = MT2Analysis<MT2Estimate>    ::readFromFile( dir + "/zinvFromGamma.root", "ZinvEstimate");
    
    zinv_ratio  = MT2Analysis<MT2Estimate>    ::readFromFile( dir + "/zinvFromGamma.root", "ZgammaRatio");
    
    purity      = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zinvFromGamma.root", "purity");
    //purity      = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/gammaControlRegion/purityMC.root", "purity");
    
    zllG_ht_monojet     = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_data_ratio.root", "zllG_data_mono_ht");
    zllG_ht     = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_data_ratio.root", "zllG_data_ht");
    zllG_nJets  = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_data_ratio.root", "zllG_data_nJets");
    zllG_nBJets = MT2Analysis<MT2EstimateSyst>::readFromFile( dir + "/zllGammaRatio/zllG_data_ratio.root", "zllG_data_nBJets");
    
  }
  zinv->setName("zinv");
  //zinv->addToFile( mc_fileName, true );


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
  //llep->addToFile( mc_fileName, true );


  std::set<MT2Region> regions = data->getRegions();

  std::set<MT2Region> inclusiveRegions=  zllG_ht->getRegions();
  MT2Region inclusiveRegion( (*inclusiveRegions.begin() ) );

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

  
  // first create template datacards
  std::string path_templ = dir + "/datacard_templates";
  system(Form("mkdir -p %s", path_templ.c_str()));

  
  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

     TH1D* this_data = data->get(*iR)->yield;

     TH1D* this_qcd;
     TH1D* this_qcdCR;
     TH1D* this_qcd_ratio;
     TH1D* this_qcd_purity;
     TH1D* this_qcd_fjets;
     TH1D* this_qcd_rb;
     TH1D* this_qcd_ratioSystFit;
     
     std::string qcd_fjetsCR_name;       
     std::string qcd_rbCR_name;       
     
     if( iR->nJetsMax()==1 ){

       this_qcd = qcd_monojet->get(*iR)->yield;
       this_qcdCR = qcdCR_monojet->get(*iR)->yield;
       this_qcd_ratio = qcd_ratio_monojet->get(*iR)->yield;

     }
     else{
       
       this_qcd = qcd->get(*iR)->yield;
       this_qcdCR = qcdCR->get(*iR)->yield;
       this_qcd_ratio = qcd_ratio->get(*iR)->yield;
       this_qcd_purity = qcd_purity->get(*iR)->yield;
       this_qcd_ratioSystFit = qcd_ratioSystFit->get(*iR)->yield;

       if( iR->htMin() < 450 ){
	 
	 MT2Region* thisQCDCR = qcd_fjets_vlht->matchRegion(*iR);
	 
	 this_qcd_fjets = qcd_fjets_vlht->get(*thisQCDCR)->yield;
	 
	 qcd_fjetsCR_name = thisQCDCR->getName();
	 
       }
       else{
	 
	 MT2Region* thisQCDCR = qcd_fjets->matchRegion(*iR);
	 this_qcd_fjets = qcd_fjets->get(*thisQCDCR)->yield;
	 
	 qcd_fjetsCR_name = thisQCDCR->getName();
       
       }
       
       MT2Region* thisQCDCR;
       if( iR->nBJetsMin()==3 && iR->nJetsMin()==2 ) 
	 thisQCDCR = new MT2Region( 200, -1, 4, 6, 0, -1 );
       else  
	 thisQCDCR = qcd_rb->matchRegion(*iR);
       
       this_qcd_rb = qcd_rb->get(*thisQCDCR)->yield;
       
       qcd_rbCR_name = thisQCDCR->getName();

     }


     TH1D* this_zinv = zinv->get(*iR)->yield;
     TH1D* this_zinv_ratio =  zinv_ratio->get(*iR)->yield;
     TH1D* this_zinvCR;
     this_zinvCR = zinvCR->get(*iR)->yield;

     TGraphAsymmErrors* this_zinv_purity;
     this_zinv_purity = purity->get(*iR)->getGraph();

     int Ngamma=0;

     TH1D* this_llep = llep->get(*iR)->yield;
     TH1D* this_llep_ratio = llep_ratio->get(*iR)->yield;
     TH1D* this_llepCR = llepCR->get(*iR)->yield;

     float N_llep_CR = this_llepCR->Integral();
     std::string llepCR_name;
     if(iR->nJetsMin()>=7 && iR->nBJetsMin()>=1){
       MT2Region* thisCR = new MT2Region(iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), 1, 2);
       llepCR_name = thisCR->getName();
     }
     else
       llepCR_name = iR->getName();

     
         
     unsigned iEmptyZinvBin=this_data->GetNbinsX()+1;
     int nEmptyCR=0;

     int nBins = this_data->GetNbinsX();

     float shapeErr_zinv=0.;
     float shapeErr_llep=0.;

     for( int iBin=1; iBin<this_data->GetNbinsX()+1; ++iBin ) {
       
       float relativeErr;
       if( fabs(this_zinv->GetBinContent(iBin))>0 )
	 relativeErr = 0.4 / (nBins-1) * (iBin-1);
////       //relativeErr = 1.0 / (nBins-1) * (iBin-1);
//	 relativeErr = 0.2 / ((nBins-1) * (nBins-1)) * (iBin-1) * (iBin-1);
       else
	 relativeErr = 0.0;
       
       shapeErr_zinv+=relativeErr*fabs(this_zinv->GetBinContent(iBin));
     
     }

     for( int iBin=1; iBin<this_data->GetNbinsX()+1; ++iBin ) {
       
       float relativeErr;
       if( this_llep->GetBinContent(iBin)>0 )
	 relativeErr = 0.4 / (nBins-1) * (iBin-1);
////       //relativeErr = 1.0 / (nBins-1) * (iBin-1);
//	 relativeErr = 0.2 / ((nBins-1) * (nBins-1)) * (iBin-1) * (iBin-1);
       else
	 relativeErr = 0.0;
       
      shapeErr_llep+=relativeErr*this_llep->GetBinContent(iBin);
     
     }
     
     float lastR_zinv;
     float lastR_llep;
     float lastR_qcd_mono;

     for( int iBin=1; iBin<this_data->GetNbinsX()+1; ++iBin ) {
       
       if(this_data->GetBinLowEdge( iBin ) > iR->htMax() && iR->htMax()>0 ) continue;
       
       float mt2Min = this_data->GetBinLowEdge( iBin );
       float mt2Max = (iBin==this_data->GetNbinsX()) ?  -1. : this_data->GetBinLowEdge( iBin+1 );
       
       std::string binName;
       if( iR->nJetsMax() ==1 ){ 
	   binName = std::string( Form("%s_m0toInf", iR->getName().c_str() ) );
       }
       else{
	 
	 if( mt2Max>=0. )
	   binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
	 else
	   binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );

       }

       std::string htName;
       htName = iR->htRegion()->getName();

       int nJetsMin = iR->nJetsMin();
       int nJetsMax = iR->nJetsMax();
       int nBJetsMin = iR->nBJetsMin();
       int nBJetsMax = iR->nBJetsMax();
       
       std::string jName;
       jName = iR->sigRegion()->getSingleJetString( "j", nJetsMin,  nJetsMax  );

       //       std::cout << jName << std::endl;

       std::string bName;
       bName = iR->sigRegion()->getSingleJetString( "b", nBJetsMin,  nBJetsMax  );
       
       //       std::cout << bName << std::endl;
       
       std::string datacardName( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
       
       std::ifstream thisDatacard( datacardName.c_str() );
       if( thisDatacard.good() ) continue;
       
       std::ofstream datacard( datacardName.c_str() );
       
       std::string tableName( Form("%s/table_%s.txt", path_templ.c_str(), binName.c_str()) );
       std::ofstream table( tableName.c_str() );
       table << std::setprecision(3);
       
       
       datacard << "imax 1" << std::endl;
       datacard << "jmax 3" << std::endl;
       datacard << "kmax *" << std::endl;
       datacard << "-------------" << std::endl;
       datacard << std::endl << std::endl;
       
       
       datacard << std::fixed;
       datacard << std::setprecision(3) << std::endl << std::endl;
       datacard << "bin  " << binName<< std::endl;
       //       datacard << "observation  " << (fabs(this_qcd->GetBinContent(iBin))+fabs(this_zinv->GetBinContent(iBin))+fabs(this_llep->GetBinContent(iBin))) << std::endl;
       datacard << "observation  " << this_data->GetBinContent(iBin) << std::endl;
       datacard << "-------------" << std::endl;
       datacard << std::endl << std::endl;
       
       
       float yield_llep = fabs(this_llep->GetBinContent(iBin));
       float yield_qcd = fabs(this_qcd ->GetBinContent(iBin));
       //       if( !(iR->nJetsMax()==1) && this_qcd_purity->GetBinContent(iBin) > 0 ) yield_qcd *= (this_qcd_purity->GetBinContent(iBin));

//       ////// Suicidal attempt
//       if( iR->nJetsMax() >  1 || iR->nJetsMax() < 0 ){
//	 
//	 if( iR->htMin() < 1500 )
//	   yield_qcd = yield_qcd*0.3/(this_qcd_ratio->GetBinContent(iBin));
//	 else
//	   yield_qcd = yield_qcd*0.6/(this_qcd_ratio->GetBinContent(iBin));
//       
//       }
//       /////

       float yield_zinv = fabs(this_zinv->GetBinContent(iBin));
       
       // sig qcd zinv llep
       datacard << "bin \t" << binName << "\t" << binName << "\t" << binName << "\t" << binName << std::endl;
       datacard << "process \t sig \t zinv \t llep \t qcd" << std::endl;
       datacard << "process \t 0 \t 1 \t 2 \t 3" << std::endl;
       datacard << "rate \t XXX";
       datacard << " \t " << yield_zinv << " \t " << yield_llep << " \t " << yield_qcd << std::endl;
       datacard << "-------------" << std::endl;
       
       datacard << "lumi_syst    lnN    " << 1.+err_lumi_corr << " - - -" << std::endl;
       datacard << "sig_MCstat_" << binName << " lnN UUU - - -" << std::endl;
       datacard << "sig_isrSyst lnN III - - -" << std::endl;
       datacard << "sig_bTagHeavySyst lnN HHH - - -" << std::endl;
       datacard << "sig_bTagLightSyst lnN LLL - - -" << std::endl;
       
       
       // these needed for table
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


           /*// correlated:
           datacard << "zinv_ZGratio_" << zinvCR_name << " lnN   - " << 1.+err_zinv_corr << " - -" << std::endl;
           zinv_systUp += err_zinv_corr*err_zinv_corr;
           zinv_systDn += err_zinv_corr*err_zinv_corr;
	   */
	   
	   int thisBinNJ;
	   int thisBinNB;
	   int thisBinHT;

	   float thisErrNJUp;
	   float thisErrNBUp;
	   float thisErrHTUp;

	   float thisErrNJDn;
	   float thisErrNBDn;
	   float thisErrHTDn;

//	   float thisErrNJ;
//	   float thisErrNB;
//	   float thisErrHT;

	   if( iR->nJetsMax()>1 || iR->nJetsMax()<0){
	     
	      thisBinNJ = this_zllG_nJets->FindBin(iR->nJetsMin());
	      if( iR->nBJetsMin() < 3 )
		thisBinNB = this_zllG_nBJets->FindBin(iR->nBJetsMin());
	      else
		thisBinNB = this_zllG_nBJets->FindBin(2);
	      thisBinHT = this_zllG_ht->FindBin(iR->htMin()+1.);
	     
	      thisErrNJUp = ( this_zllG_nJets->GetBinContent(thisBinNJ) > 0 )  ? ( thisUp_zllG_nJets->GetBinContent(thisBinNJ) - this_zllG_nJets->GetBinContent(thisBinNJ) )  / this_zllG_nJets->GetBinContent(thisBinNJ)  : 1.0;
	      thisErrNBUp = ( this_zllG_nBJets->GetBinContent(thisBinNB) > 0 ) ? ( thisUp_zllG_nBJets->GetBinContent(thisBinNB) - this_zllG_nBJets->GetBinContent(thisBinNB) )/ this_zllG_nBJets->GetBinContent(thisBinNB) : 1.0;
	      thisErrHTUp = ( this_zllG_ht->GetBinContent(thisBinHT) > 0 )     ? ( thisUp_zllG_ht->GetBinContent(thisBinHT) - this_zllG_ht->GetBinContent(thisBinHT) )           / this_zllG_ht->GetBinContent(thisBinHT)     : 1.0;

	      thisErrNJDn = ( this_zllG_nJets->GetBinContent(thisBinNJ) > 0 )  ? ( this_zllG_nJets->GetBinContent(thisBinNJ)  - thisDn_zllG_nJets->GetBinContent(thisBinNJ) ) / this_zllG_nJets->GetBinContent(thisBinNJ)  : 1.0;
	      thisErrNBDn = ( this_zllG_nBJets->GetBinContent(thisBinNB) > 0 ) ? ( this_zllG_nBJets->GetBinContent(thisBinNB) - thisDn_zllG_nBJets->GetBinContent(thisBinNB) )/ this_zllG_nBJets->GetBinContent(thisBinNB) : 1.0;
	      thisErrHTDn = ( this_zllG_ht->GetBinContent(thisBinHT) > 0 )     ? ( this_zllG_ht->GetBinContent(thisBinHT)     - thisDn_zllG_ht->GetBinContent(thisBinHT) )    / this_zllG_ht->GetBinContent(thisBinHT)     : 1.0;

	      if( iR->nBJetsMin() >= 3 ){
		
		thisErrNBUp*=2;
		thisErrNBDn*=2;

	      }

//	      thisErrNJ = ( this_zllG_nJets->GetBinContent(thisBinNJ) > 0 )  ? this_zllG_nJets->GetBinError(thisBinNJ) / this_zllG_nJets->GetBinContent(thisBinNJ)  : 1.0;
//	      thisErrNB = ( this_zllG_nBJets->GetBinContent(thisBinNB) > 0 ) ? this_zllG_nBJets->GetBinError(thisBinNB)/ this_zllG_nBJets->GetBinContent(thisBinNB) : 1.0;
//	      thisErrHT = ( this_zllG_ht->GetBinContent(thisBinHT) > 0 )     ? this_zllG_ht->GetBinError(thisBinHT)    / this_zllG_ht->GetBinContent(thisBinHT)     : 1.0;
	      
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
	     thisErrHTDn = ( this_zllG_ht_monojet->GetBinContent(thisBinHT) > 0 )     ? ( this_zllG_ht_monojet->GetBinContent(thisBinHT)     - thisDn_zllG_ht_monojet->GetBinContent(thisBinHT) )    / this_zllG_ht_monojet->GetBinContent(thisBinHT)     : 1.0;
	     
	   }


           datacard << "zinv_doubleRatioOffset lnN   - " << 1.+err_zinv_doubleRatioOffset << " - -" << std::endl;
           zinv_systUp += err_zinv_doubleRatioOffset*err_zinv_doubleRatioOffset;
           zinv_systDn += err_zinv_doubleRatioOffset*err_zinv_doubleRatioOffset;

	   float err_ZGUp = TMath::Sqrt( thisErrNJUp*thisErrNJUp + thisErrNBUp*thisErrNBUp + thisErrHTUp*thisErrHTUp );
	   float err_ZGDn = TMath::Sqrt( thisErrNJDn*thisErrNJDn + thisErrNBDn*thisErrNBDn + thisErrHTDn*thisErrHTDn );
	   
	   err_ZGDn = (err_ZGDn < 1) ? err_ZGDn : 0.99;

//           datacard << "zinv_ZGratio_" << zinvCR_name << " lnN   - " << 1.+err_ZGUp << "/" << 1.-err_ZGDn << " - -" << std::endl;
           zinv_systUp += err_ZGUp*err_ZGUp;
           zinv_systDn += err_ZGDn*err_ZGDn;

           datacard << "zinv_ZGratio_" << htName << " lnN   - " << 1.+thisErrHTUp << "/" << 1.-thisErrHTDn << " - -" << std::endl;
           datacard << "zinv_ZGratio_" << jName << " lnN   - " << 1.+thisErrNJUp << "/" << 1.-thisErrNJDn << " - -" << std::endl;
           datacard << "zinv_ZGratio_" << bName << " lnN   - " << 1.+thisErrNBUp << "/" << 1.-thisErrNBDn << " - -" << std::endl;

         }
	 
         // uncorrelated:
         float thisError_zinv_uncorr_rel = this_zinv->GetBinError(iBin)/yield_zinv;
	 
	 if( (!use_extrapolation && iR->nBJetsMin()>=2) || iR->nBJetsMin()>3 ) {
	     
	   datacard << "zinv_MC_" << binName << " lnN - " << 1.+err_zinv_uncorr_2b << " - -" << std::endl;
	   zinv_systUp += err_zinv_uncorr_2b*err_zinv_uncorr_2b;
	   zinv_systDn += err_zinv_uncorr_2b*err_zinv_uncorr_2b;
	   
	 } else {
	   
	   if ( !use_extrapolation )
	     Ngamma = round(this_zinvCR->GetBinContent(iBin));
	   else
	     Ngamma = round(this_zinvCR->Integral());
	   
	   Double_t x_tmp, p, p_errUp, p_errDown;
	   this_zinv_purity->GetPoint( 0, x_tmp, p );
	   p_errUp   = this_zinv_purity->GetErrorYhigh( 0 );
	   p_errDown = this_zinv_purity->GetErrorYlow ( 0 ); 
//	   this_zinv_purity->GetPoint( 0, x_tmp, p );
//	   p_errUp   = this_zinv_purity->GetErrorYhigh( 0 );
//	   p_errDown = this_zinv_purity->GetErrorYlow ( 0 ); 
	   
	   if( Ngamma>0 && p>0 ) {
	     
	     datacard << "zinv_puritySyst_" << zinvCR_name << " lnN  - " << 1.+err_zinv_puritySyst << " - -" << std::endl;
	     zinv_systUp += err_zinv_puritySyst*err_zinv_puritySyst;
             zinv_systDn += err_zinv_puritySyst*err_zinv_puritySyst;
	     
	     p = fabs(p);
	     
	     if( fabs(p_errDown)/p < 1. )  
	       datacard << "zinv_purity_" << zinvCR_name << " lnN  - " << 1.+p_errUp/p << "/" << 1.-p_errDown/p << " - -" << std::endl;
	     else 
	       datacard << "zinv_purity_" << zinvCR_name << " lnN  - " << 1.+p_errUp/p << "/0.01 - -" << std::endl;
	   
	     zinv_systUp += (p_errUp/p)*(p_errUp/p);
	     zinv_systDn += (p_errDown/p)*(p_errDown/p);
	   
	   }
	   else if( Ngamma>0 ){

	     datacard << "zinv_puritySyst_" << zinvCR_name << " lnN  - " << 1.+err_zinv_puritySyst << " - -" << std::endl;
  
	     datacard << "zinv_purity_" << zinvCR_name << " lnN  - 1.0/0.01 - -" << std::endl;
	     
	     zinv_systUp += 1.0;
	     zinv_systDn += 1.0;
	     
	   }
	   
             
	   float R = fabs(this_zinv_ratio->GetBinContent(iBin));
	   float relativeErr;
	   
	   if (R>0) lastR_zinv = R;

	   if( fabs(this_zinv->GetBinContent(iBin))>0 )
	     relativeErr = 0.4 / (nBins-1) * (iBin-1);
	     //	     relativeErr = 0.2 / ((nBins-1) * (nBins-1)) * (iBin-1) * (iBin-1);
	   else
	     relativeErr = 0.0;
	   
	   if( !use_extrapolation )
	     datacard << "zinv_CRstat_" << std::setprecision(5) << gammaConvention( yield_zinv, Ngamma, 1, binName, binName, (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv ) << std::setprecision(3) << std::endl;
	   else {
	     datacard << "zinv_CRstat_" << std::setprecision(5) << gammaConvention( yield_zinv, Ngamma, 1, zinvCR_name, binName, (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv ) << std::setprecision(3) << std::endl;
	     if( nBins>1 ){
	       if( iBin==1 && fabs(yield_zinv)>0 ){
		 
		 datacard << "zinv_shape_" << zinvCR_name << " lnN  - " << 1.-shapeErr_zinv/fabs(yield_zinv) << " - - " << std::endl;
		 zinv_systUp += (shapeErr_zinv/fabs(yield_zinv))*(shapeErr_zinv/fabs(yield_zinv));
		 zinv_systDn += (shapeErr_zinv/fabs(yield_zinv))*(shapeErr_zinv/fabs(yield_zinv));
	       }
	       else{
		 datacard << "zinv_shape_" << zinvCR_name << " lnN  - " << 1.+relativeErr << " - - " << std::endl;
		 zinv_systUp += relativeErr*relativeErr;
		 zinv_systDn += relativeErr*relativeErr;
	       }
	     }
	   }

	   double yield_zinv_up, yield_zinv_dn;
	   RooHistError::instance().getPoissonInterval(Ngamma,yield_zinv_dn,yield_zinv_up,1.);
	   yield_zinv_up *= (Ngamma>0) ? yield_zinv/Ngamma : (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv;
	   yield_zinv_dn *= (Ngamma>0) ? yield_zinv/Ngamma : (R>0) ? ((R>0.5) ? 0.5 : R) : lastR_zinv;
	   
	   zinv_statUp = yield_zinv_up-yield_zinv;
	   zinv_statDn = yield_zinv-yield_zinv_dn;
	   
	   float alphaErr;
	   if ( R>0 ) alphaErr= this_zinv_ratio->GetBinError(iBin)/R;
	   else alphaErr=1.0;
	   datacard << "zinv_alphaErr_" << binName << " lnN  - " << 1.+alphaErr << " - -" << std::endl;
	   //datacard << "zinv_alphaErr_" << zinvCR_name << " lnN  - " << 1.+1. << " - -" << std::endl;
	   zinv_systUp += alphaErr*alphaErr;
	   zinv_systDn += alphaErr*alphaErr;
	   
	   
	 } // if extrapolation
	 
	 zinv_nCR = Ngamma;

       } // if zinv
       
       
       
       
       // LOST LEPTON SYSTEMATICS:
       if( yield_llep>=0. ) {
	 
         // correlated within the SR (stat-like):
         float llep_stat_err = (N_llep_CR>0) ? 1./sqrt((float)N_llep_CR) : 0.;
         float llep_tot_err = sqrt( llep_stat_err*llep_stat_err + err_llep_lepEff*err_llep_lepEff );
	 
	 float Rllep = this_llep_ratio->GetBinContent(iBin);
	 
	 if( Rllep > 0 ) lastR_llep = Rllep;

	 float alphaErr_llep;
	 if ( Rllep>0 ) alphaErr_llep= this_llep_ratio->GetBinError(iBin)/Rllep;
	 else alphaErr_llep=1.0;
	
	 float relativeErr_llep;
	 if( this_llep->GetBinContent(iBin)>0 )
	   relativeErr_llep = 0.4 / (nBins-1) * (iBin-1);
	 //relativeErr_llep = 0.2 / ((nBins-1) * (nBins-1)) * (iBin-1) * (iBin-1);
	 else
	   relativeErr_llep = 0.0;

	 datacard << "llep_lepeff_" << llepCR_name << "  lnN  - - " << 1.+err_llep_lepEff << " -" << std::endl;
	 llep_systUp += err_llep_lepEff*err_llep_lepEff;
	 llep_systDn += err_llep_lepEff*err_llep_lepEff;
	 
	 datacard << "llep_CRstat_" << gammaConvention( yield_llep, round(N_llep_CR), 2, llepCR_name, binName, (Rllep>0) ? ( (Rllep>3) ? 2 : Rllep ) : lastR_llep ) << std::endl;
        
	 double yield_llep_up, yield_llep_dn;
	 RooHistError::instance().getPoissonInterval(round(N_llep_CR),yield_llep_dn,yield_llep_up,1.);
	 yield_llep_up *= (round(N_llep_CR)>0) ? yield_llep/round(N_llep_CR) : (Rllep>0) ? ( (Rllep>3) ? 2 : Rllep ) : lastR_llep;
	 yield_llep_dn *= (round(N_llep_CR)>0) ? yield_llep/round(N_llep_CR) : (Rllep>0) ? ( (Rllep>3) ? 2 : Rllep ) : lastR_llep;
	 llep_statUp = yield_llep_up-yield_llep;
	 llep_statDn = yield_llep-yield_llep_dn;

//	 if( iR->nJetsMin()==7 && iR->nBJetsMin()>=1 )
//	   datacard << "llep_bTag_" << llepCR_name << " lnN - - 1.2 -" << std::endl;
	     
	 if( yield_llep>=0. ) {

	   float err_llep_mcstat;
	   if(yield_llep>0) err_llep_mcstat = this_llep->GetBinError(iBin)/yield_llep;
	   else err_llep_mcstat = 1.0;
	 
	   datacard << "llep_MCstat_" << binName << " lnN  - - " << 1.+err_llep_mcstat << " -" << std::endl;
	   llep_systUp += err_llep_mcstat*err_llep_mcstat;
	   llep_systDn += err_llep_mcstat*err_llep_mcstat;
	     
	   if( nBins > 1 ){
	     if( iBin==1 && yield_llep>0 ){
	       
	       datacard << "llep_shape_" << llepCR_name << " lnN - - " << 1.-shapeErr_llep/yield_llep << " - " << std::endl;
	       llep_systUp += (shapeErr_llep/yield_llep)*(shapeErr_llep/yield_llep);
	       llep_systDn += (shapeErr_llep/yield_llep)*(shapeErr_llep/yield_llep);
	     
	     }
	     else{
	     
	       datacard << "llep_shape_" << llepCR_name << " lnN - - " << 1+relativeErr_llep << " - " << std::endl;
	       llep_systUp += relativeErr_llep*relativeErr_llep;
	       llep_systDn += relativeErr_llep*relativeErr_llep;
	     
	     }
	   }

	   if( iR->htMin()==200 && iR->nJetsMin()>=7 ){
	 
	     datacard << "llep_alpha_" << llepCR_name << " lnN - - " << 1.+0.4 << " - " << std::endl;
	     llep_systUp += 0.4*0.4;
	     llep_systDn += 0.4*0.4;
	   
	   } else if( iR->nJetsMin()>=7 && iR->nBJetsMin()>=3 ){

	     datacard << "llep_alpha_" << llepCR_name << " lnN - - " << 1.+0.15 << " - " << std::endl;
             llep_systUp += 0.15*0.15;
             llep_systDn += 0.15*0.15;

	   } else{

	     datacard << "llep_alpha_" << llepCR_name << " lnN - - " << 1.+0.1 << " - " << std::endl;
             llep_systUp += 0.1*0.1;
             llep_systDn += 0.1*0.1;

	   }
	   
	   
	 }
	 
	 llep_nCR = N_llep_CR;
	 
       }



       // QCD SYSTEMATICS:
       int NQCD_cr;

       if( yield_qcd>=0. ) {

	 if( iR->nJetsMax()==1 ){
	   
	   NQCD_cr = round(this_qcdCR->GetBinContent(iBin));
	   float thisRqcd = this_qcd_ratio->GetBinContent(iBin);
	   
	   if( thisRqcd > 0 ) lastR_qcd_mono = thisRqcd;
	   
	   datacard << "qcd_CRstat_" << std::setprecision(6) << gammaConvention( yield_qcd, NQCD_cr, 3, binName, binName, (thisRqcd>0) ? thisRqcd : lastR_qcd_mono ) << std::setprecision(3) << std::endl;

           double yield_qcd_up, yield_qcd_dn;
	   RooHistError::instance().getPoissonInterval(NQCD_cr,yield_qcd_dn,yield_qcd_up,1.);
           yield_qcd_up *= (NQCD_cr>0.) ? yield_qcd/NQCD_cr : lastR_qcd_mono;
           yield_qcd_dn *= (NQCD_cr>0.) ? yield_qcd/NQCD_cr : lastR_qcd_mono;
           qcd_statUp = yield_qcd_up-yield_qcd;
           qcd_statDn = yield_qcd-yield_qcd_dn;

	   float thisRQCDErr_monojet = ( thisRqcd > 0 ) ? this_qcd_ratio->GetBinError(iBin)/thisRqcd : 1.0;
	   datacard << "qcd_alphaErr_" << binName << " lnN - - - " << 1.+thisRQCDErr_monojet << std::endl;
           qcd_systUp += thisRQCDErr_monojet*thisRQCDErr_monojet;
           qcd_systDn += thisRQCDErr_monojet*thisRQCDErr_monojet;

//	   datacard << "qcd_syst_" << binName << " lnN - - - " << 1.+err_qcd_uncorr << std::endl;
//	   qcd_systUp += err_qcd_uncorr*err_qcd_uncorr;
//	   qcd_systDn += err_qcd_uncorr*err_qcd_uncorr;
	 
	 }
	 
	 else{

	   NQCD_cr = round(this_qcdCR->GetBinContent(iBin));
	   float thisRqcd = this_qcd_ratio->GetBinContent(iBin);
	   float thisRqcdErr = this_qcd_ratio->GetBinError(iBin);
	   
	   float thisQCDPurity = this_qcd_purity->GetBinContent(iBin);
	   float thisQCDPurityErr = this_qcd_purity->GetBinError(iBin);
	   
	   if( thisRqcd>0 && thisRqcdErr>0) thisRqcdErr=thisRqcdErr/thisRqcd;
	   else thisRqcdErr=0.0;

	   if( thisQCDPurity>0 && thisQCDPurityErr>0) thisQCDPurityErr=thisQCDPurityErr/thisQCDPurity;
	   else thisQCDPurityErr=0.0;
	   
	   float thisRFitVarErr = this_qcd_ratioSystFit->GetBinContent(iBin);
	   
//	   ////// Suicidal attempt
//	   if( iR->htMin()<1500 ){
//	     thisRqcd = 0.3;
//	     thisRqcdErr = 0.1/0.3;
//	   }
//	   else{
//	     thisRqcd = 0.6;
//             thisRqcdErr = 0.2/0.6;
//	   }
//	   /////

	   int thisFJetsBin = this_qcd_fjets->FindBin( iR->nJetsMin() );
	   int thisRBBin = this_qcd_rb->FindBin( iR->nBJetsMin() );
	   
	   float thisFJets = this_qcd_fjets->GetBinContent(thisFJetsBin);
	   float thisRB = this_qcd_rb->GetBinContent(thisRBBin);
	   
	   float thisFJetsErr = this_qcd_fjets->GetBinError(thisFJetsBin);
	   float thisRBErr = this_qcd_rb->GetBinError(thisRBBin);
	   
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

	   datacard << "qcd_CRstat_" << std::setprecision(6) << gammaConvention( yield_qcd, NQCD_cr, 3, binName, binName, (thisRqcd*thisFractions>0) ? thisRqcd*thisFractions : 1.0 ) << std::setprecision(3) << std::endl;   

	   double yield_qcd_up, yield_qcd_dn;
	   RooHistError::instance().getPoissonInterval(NQCD_cr,yield_qcd_dn,yield_qcd_up,1.);
	   yield_qcd_up *= (NQCD_cr>0.) ? yield_qcd/NQCD_cr : thisRqcd*thisFJets*thisRB;
	   yield_qcd_dn *= (NQCD_cr>0.) ? yield_qcd/NQCD_cr : thisRqcd*thisFJets*thisRB;
	   qcd_statUp = yield_qcd_up-yield_qcd;
	   qcd_statDn = yield_qcd-yield_qcd_dn;

	   datacard << "qcd_FJRBsyst_" << binName << " lnN - - - " <<  1.+thisFractionsErr  << std::endl;
	   qcd_systUp += thisFractionsErr*thisFractionsErr;
           qcd_systDn += thisFractionsErr*thisFractionsErr;

//	   ////// Suicidal attempt (to comment) 
//	   datacard << "qcd_RPHIstat_" << htName << " lnN - - - " <<  1.+thisRqcdErr<<"/0.01" << std::endl;
//	   qcd_systUp += thisRqcdErr*thisRqcdErr;
//           qcd_systDn += thisRqcdErr*thisRqcdErr;
//	   //////

	   ////// Suicidal attempt (to uncomment)
	   datacard << "qcd_RPHIstat_" << htName << " lnN - - - " <<  1.+thisRqcdErr  << std::endl;
	   qcd_systUp += thisRqcdErr*thisRqcdErr;
	   qcd_systDn += thisRqcdErr*thisRqcdErr;
	 
	   datacard << "qcd_RPHIsyst_" << htName << " lnN - - - " <<  1.+thisRFitVarErr  << std::endl;
	   qcd_systUp += thisRFitVarErr*thisRFitVarErr;
	   qcd_systDn += thisRFitVarErr*thisRFitVarErr;
	   //////
	 }
	
	 qcd_nCR = NQCD_cr;
	 
       }


       datacard.close();

       std::cout << "-> Created template datacard: " << datacardName << std::endl;



       // make them absolute uncertainties

       zinv_systUp = yield_zinv*sqrt(zinv_systUp);
       zinv_systDn = yield_zinv*sqrt(zinv_systDn);

       llep_systUp = yield_llep*sqrt(llep_systUp);
       llep_systDn = yield_llep*sqrt(llep_systDn);

       qcd_systUp = yield_qcd*sqrt(qcd_systUp);
       qcd_systDn = yield_qcd*sqrt(qcd_systDn);

       // now print the table:
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
  //////  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( mc_fileName, "SMS", "" );

//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagHeavy = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1bbbb_1500_100_eth.root", "T1bbbb_1500_100", "btagsf_heavy" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagLight = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1bbbb_1500_100_eth.root", "T1bbbb_1500_100", "btagsf_light" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_isr       = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1bbbb_1500_100_eth.root", "T1bbbb_1500_100", "isr" );

//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagHeavy = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1bbbb_1000_900_eth.root", "T1bbbb_1000_900", "btagsf_heavy" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagLight = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1bbbb_1000_900_eth.root", "T1bbbb_1000_900", "btagsf_light" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_isr       = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1bbbb_1000_900_eth.root", "T1bbbb_1000_900", "isr" );

//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagHeavy = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1qqqq_1000_800_eth.root", "T1qqqq_1000_800", "btagsf_heavy" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagLight = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1qqqq_1000_800_eth.root", "T1qqqq_1000_800", "btagsf_light" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_isr       = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1qqqq_1000_800_eth.root", "T1qqqq_1000_800", "isr" );

//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagHeavy = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1qqqq_1400_100_eth.root", "T1qqqq_1400_100", "btagsf_heavy" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagLight = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1qqqq_1400_100_eth.root", "T1qqqq_1400_100", "btagsf_light" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_isr       = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1qqqq_1400_100_eth.root", "T1qqqq_1400_100", "isr" );

//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagHeavy = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1tttt_1200_800_eth.root", "T1tttt_1200_800", "btagsf_heavy" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagLight = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1tttt_1200_800_eth.root", "T1tttt_1200_800", "btagsf_light" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_isr       = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1tttt_1200_800_eth.root", "T1tttt_1200_800", "isr" );
 
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagHeavy = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1tttt_1500_100_eth.root", "T1tttt_1500_100", "btagsf_heavy" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagLight = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1tttt_1500_100_eth.root", "T1tttt_1500_100", "btagsf_light" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_isr       = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1tttt_1500_100_eth.root", "T1tttt_1500_100", "isr" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals           = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/scratch/mmasciov/T1tttt_1500_100_eth.root", "T1tttt_1500_100", "isr" );

//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagHeavy = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/shome/mmasciov/signalScansFromDominick/T1bbbb_eth.root", "T1bbbb", "btagsf_heavy" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagLight = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/shome/mmasciov/signalScansFromDominick/T1bbbb_eth.root", "T1bbbb", "btagsf_light" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_isr       = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/shome/mmasciov/signalScansFromDominick/T1bbbb_eth.root", "T1bbbb", "isr" );
//  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals           = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "/shome/mmasciov/signalScansFromDominick/T1bbbb_eth.root", "T1bbbb", "isr" );
  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagHeavy = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "./signalScansFromDominick/T2cc_eth.root", "T2cc", "btagsf_heavy" );
  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_bTagLight = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "./signalScansFromDominick/T2cc_eth.root", "T2cc", "btagsf_light" );
  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals_isr       = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "./signalScansFromDominick/T2cc_eth.root", "T2cc", "isr" );
  std::vector<MT2Analysis<MT2EstimateSigSyst>*> signals           = MT2Analysis<MT2EstimateSigSyst>::readAllSystFromFile( "./signalScansFromDominick/T2cc_eth.root", "T2cc", "isr" );

 
  //std::vector<MT2Analysis<MT2Estimate>*> signals = MT2Analysis<MT2Estimate>::readAllFromFile( mc_fileName, "SMS_T1bbbb_fullScan" );
  //std::vector<MT2Analysis<MT2Estimate>*> signals = MT2Analysis<MT2Estimate>::readAllFromFile( mc_fileName, "DMS" );
  //std::vector<MT2Analysis<MT2Estimate>*> signals = MT2Analysis<MT2Estimate>::readAllFromFile( mc_fileName, "DMV" );
  
  for( unsigned  isig=0; isig<signals.size(); ++isig ) {

    std::string sigName;
    if( signals[isig]->getName().find("fullScan") != std::string::npos )
      sigName = signals[isig]->getName();
    else if( signals[isig]->getName().find("DarkMatter") != std::string::npos || signals[isig]->getName().find("prime") != std::string::npos || signals[isig]->getName().find("DM") != std::string::npos)
      sigName = signals[isig]->getName();
    else
      //      sigName = getSimpleSignalName( signals[isig]->getName() );
      sigName = signals[isig]->getName();

    //std::string sigName = getSimpleSignalName( signals[isig]->getName() );
    
    std::string path = dir + "/datacards_" + sigName;
    system(Form("mkdir -p %s", path.c_str()));
    
    std::string pathSE = "datacards_" + sigName;
    
    std::string path_mass = path;
    float xs_norm=1.;
    if( signals[isig]->getName().find("T2qq") != std::string::npos )
      xs_norm=8./10.;

    //    if( signals[isig]->getName().find("fullScan") != std::string::npos )
      for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
	
	TH3D* this_signal3d;
	TH1D* this_signalParent;

	MT2EstimateSigSyst* thisSigSystCentral = signals[isig]->get(*iR);
	if( thisSigSystCentral->yield3d!=0 && thisSigSystCentral->yield3d_systUp!=0){
	  
	  this_signal3d = signals[isig]->get(*iR)->yield3d;

	}
	else continue;
	
	
	this_signalParent = this_signal3d->ProjectionY("mParent");

	TH3D* this_signal3d_central;
	TH3D* this_signal3d_bTagHeavy_Up;
	TH3D* this_signal3d_bTagLight_Up;
	TH3D* this_signal3d_isr_Up;
	
	MT2EstimateSigSyst* thisSigSyst = signals_bTagHeavy[0]->get(*iR);

	if( thisSigSyst->yield3d!=0 && thisSigSyst->yield3d_systUp!=0){
	  
	  this_signal3d_central      = (TH3D*) signals_bTagHeavy[isig]->get(*iR)->yield3d->Clone();
	  this_signal3d_bTagHeavy_Up = (TH3D*) signals_bTagHeavy[isig]->get(*iR)->yield3d_systUp->Clone();
	  this_signal3d_bTagLight_Up = (TH3D*) signals_bTagLight[isig]->get(*iR)->yield3d_systUp->Clone();
	  this_signal3d_isr_Up       = (TH3D*) signals_isr[isig]->get(*iR)->yield3d_systDown->Clone();

	}
	else continue;
//	else{
//	  
//	  this_signal3d_central      = (TH3D*) signals[isig]->get(*iR)->yield3d->Clone();
//	  this_signal3d_bTagHeavy_Up = (TH3D*) signals[isig]->get(*iR)->yield3d->Clone();
//	  this_signal3d_bTagLight_Up = (TH3D*) signals[isig]->get(*iR)->yield3d->Clone();
//	  this_signal3d_isr_Up       = (TH3D*) signals[isig]->get(*iR)->yield3d->Clone();
//	
//	}

	for( int iBinY=1; iBinY<this_signalParent->GetNbinsX()+1; ++iBinY ){

	  float mParent = this_signalParent->GetBinLowEdge(iBinY);
	  if( !(mParent >= m1-1 && mParent < m2-1) ) continue;
	  
	  TH1D* this_signalLSP = this_signal3d_central->ProjectionZ("mLSP", 0, -1, iBinY, iBinY);
	  
	  int iBinYforZ = this_signalLSP->GetXaxis()->FindBin(mParent);
	  
	  for( int iBinZ=1; iBinZ < iBinYforZ; ++iBinZ ) {
	  	  
	    float mLSP = this_signalLSP->GetBinLowEdge(iBinZ);
	    if( !(mLSP >= m11-1 && mLSP < m22-1) ) continue;
	    
	    TH1D* this_signal = this_signal3d_central->ProjectionX("mt2", iBinY, iBinY, iBinZ, iBinZ);
	    
	    if( this_signal->Integral() <=0 ) continue;
	    
	    for( int iBin=1; iBin<this_signal->GetNbinsX()+1; ++iBin ) {
	      
	      if( this_signal->GetBinLowEdge( iBin ) > iR->htMax() && iR->htMax()>0 ) continue;
	      
	      float mt2Min = this_signal->GetBinLowEdge( iBin );
	      float mt2Max = (iBin==this_signal->GetNbinsX()) ?  -1. : this_signal->GetBinLowEdge( iBin+1 );
	      
	      
	      if( this_signal->GetBinContent(iBin) <= 0 );
	      else{
		
		std::string binName;
		if( mt2Max>=0. )
		  binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
		else
		  binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );
		
		Long_t id, size, flags, modtime;
		std::string fullPathSE = Form( "/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/datacards_T2cc_04Feb/%s_%.0f_%.0f/datacard_%s_%s_%.0f_%.0f.txt", pathSE.c_str(), mParent, mLSP, binName.c_str(), sigName.c_str(), mParent, mLSP);
		int checkFileSE = (int) gSystem->GetPathInfo(fullPathSE.c_str(),&id, &size, &flags, &modtime);
		
		std::string rmOnSE( Form("env --unset=LD_LIBRARY_PATH gfal-rm srm://t3se01.psi.ch/%s", fullPathSE.c_str()) );
		//std::string rmOnSE( Form("gfal-rm srm://t3se01.psi.ch/%s", fullPathSE.c_str()) );
		
		if( checkFileSE==0 && size==0 ){
		  
		  std::cout << "Removing. File " << fullPathSE << " exists and has zero size " << size << ". Removing." << std::endl;
		  system( rmOnSE.c_str() );
		  
		}
		else if ( checkFileSE==0 && size>0 ){
		  
		  std::cout << "Skipping. File " << fullPathSE << " exists and has non-zero size  " << size << ". Skipping." << std::endl;
		  
		  continue;
		  
		}
		
		
		std::string templateDatacard( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
		
		std::string newDatacard( Form("%s/datacard_%s_%s_%.0f_%.0f.txt", path_mass.c_str(), binName.c_str(), sigName.c_str(), mParent, mLSP) );
		std::string helpDatacard( Form("%s/datacard_%s_%s_%.0f_%.0f_forSed.txt", path_mass.c_str(), binName.c_str(), sigName.c_str(), mParent, mLSP) );

		std::ifstream thisNewDatacard( newDatacard.c_str() );
		if( thisNewDatacard.good() ) continue;

//////		std::string newDatacard( Form("%s/datacard_%s_%s.txt", path_mass.c_str(), binName.c_str(), sigName.c_str()) );
//////		std::string helpDatacard( Form("%s/datacard_%s_%s_forSed.txt", path_mass.c_str(), binName.c_str(), sigName.c_str()) );
		
		float sig = this_signal->GetBinContent(iBin);
		float sigErr = this_signal->GetBinError(iBin)/sig;
		//sig*=xs_norm;
		
		float isrErr = this_signal3d_isr_Up->GetBinContent(iBin, iBinY, iBinZ);
		isrErr = 2 - isrErr/sig;
		
		float bTagErr_heavy = this_signal3d_bTagHeavy_Up->GetBinContent(iBin, iBinY, iBinZ);
		bTagErr_heavy = bTagErr_heavy/sig;

		float bTagErr_light = this_signal3d_bTagLight_Up->GetBinContent(iBin, iBinY, iBinZ);
		bTagErr_light = bTagErr_light/sig;
		
		float totUncorrErr = 1.+sqrt(sigErr*sigErr+2*0.05*0.05); // MC stat + PDF + scales + JEC

		std::string mvCommand( Form("mv %s %s", newDatacard.c_str(), helpDatacard.c_str()) );
		std::string rmCommand( Form("rm -f %s", newDatacard.c_str()) );

		sig*=2.26355/2.26;
		
		std::string sedCommand( Form("sed 's/XXX/%.3f/' %s > %s", sig, templateDatacard.c_str(), newDatacard.c_str()) );
                system( sedCommand.c_str() );

		std::string sedCommand_uncErr( Form("sed -i 's/UUU/%.3f/' %s", totUncorrErr, newDatacard.c_str()) );
                system( sedCommand_uncErr.c_str() );

		std::string sedCommand_isrErr( Form("sed -i 's/III/%.3f/' %s", isrErr, newDatacard.c_str()) );
                system( sedCommand_isrErr.c_str() );

		std::string sedCommand_bTagHErr( Form("sed -i 's/HHH/%.3f/' %s", bTagErr_heavy, newDatacard.c_str()) );
                system( sedCommand_bTagHErr.c_str() );

		std::string sedCommand_bTagLErr( Form("sed -i 's/LLL/%.3f/' %s", bTagErr_light, newDatacard.c_str()) );
                system( sedCommand_bTagLErr.c_str() );

		

		std::string mkdirOnSE( Form("env --unset=LD_LIBRARY_PATH gfal-mkdir -p srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/datacards_T2cc_04Feb/%s_%.0f_%.0f", pathSE.c_str(), mParent, mLSP) );
		//std::string mkdirOnSE( Form("gfal-mkdir -p srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/datacards_T2cc_04Feb/%s_%.0f_%.0f", pathSE.c_str(), mParent, mLSP) );
		std::string copyOnSE( Form("xrdcp -v %s root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/mmasciov/datacards_T2cc_04Feb/%s_%.0f_%.0f/datacard_%s_%s_%.0f_%.0f.txt", newDatacard.c_str(), pathSE.c_str(), mParent, mLSP, binName.c_str(), sigName.c_str(), mParent, mLSP) );
		system( mkdirOnSE.c_str() );
		system( copyOnSE.c_str() );
		
		for(int c=0; c<3; ++c){
		  
		  checkFileSE = (int) gSystem->GetPathInfo(fullPathSE.c_str(),&id, &size, &flags, &modtime);
		  
		  if( checkFileSE==0 && size==0 ){
		    
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

//		std::string sedCommand( Form("sed 's/XXX/%.3f/g' %s > %s", sig, templateDatacard.c_str(), newDatacard.c_str()) );
//		system( sedCommand.c_str() );
//		
//		system( mvCommand.c_str() );
//
//		std::string sedCommand_uncErr( Form("sed 's/UUU/%.3f/g' %s > %s", totUncorrErr, helpDatacard.c_str(), newDatacard.c_str()) );
//		system( sedCommand_uncErr.c_str() );
//
//		system( mvCommand.c_str() );
//
//		std::string sedCommand_isrErr( Form("sed 's/III/%.3f/g' %s > %s", isrErr, helpDatacard.c_str(), newDatacard.c_str()) );
//		system( sedCommand_isrErr.c_str() );
//
//		system( mvCommand.c_str() );
//
//		std::string sedCommand_bTagHErr( Form("sed 's/HHH/%.3f/g' %s > %s", bTagErr_heavy, helpDatacard.c_str(), newDatacard.c_str()) );
//		system( sedCommand_bTagHErr.c_str() );
//
//		system( mvCommand.c_str() );
//
//		std::string sedCommand_bTagLErr( Form("sed 's/LLL/%.3f/g' %s > %s", bTagErr_light, helpDatacard.c_str(), newDatacard.c_str()) );
//		system( sedCommand_bTagLErr.c_str() );
//
//		system( rmCommand.c_str() );
		
	      }
	      
	    } // for bins X (MT2)
	  } // for bins Z (mLSP)
	}// for bins Y (mParent)      
      } // for regions
    
//    else
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




MT2Analysis<MT2Estimate>* get( const std::string& name, std::vector< MT2Analysis<MT2Estimate>* > analyses, const std::string& name1, const std::string& name2, const std::string& name3, const std::string& name4 ) {


  std::cout << "Looking for: " << name << std::endl;
  MT2Analysis<MT2Estimate>* returnAnalysis = new MT2Analysis<MT2Estimate>( name, analyses[0]->getHTRegions(), analyses[0]->getSignalRegions() );

  for( unsigned i=0; i<analyses.size(); ++i ) {

    if( analyses[i]->getName() == name1 || analyses[i]->getName() == name2 || analyses[i]->getName() == name3 || analyses[i]->getName() == name4 ) {
      std::cout << "  added: " << analyses[i]->getName() << std::endl;
      (*returnAnalysis) += (*analyses[i]);
    }

  }

  return returnAnalysis;

}




void writeToTemplateFile( TFile* file, MT2Analysis<MT2Estimate>* analysis, float err_uncorr ) {

  // err_uncorr is the bin-by-bin error
  // if it's zero, no error will be assigned
  // if it's > 0., it needs to be set as a fractional error (eg. 0.03 will give a 3% error)
  // if it's negative (-1), the histogram bin error will be used

  file->cd();

  TString analysisName(analysis->getName());

  std::set<MT2Region> regions = analysis->getRegions();
  
  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

      TH1D* h1 = analysis->get( *iR )->yield;

      if(h1->Integral() == 0.){
        for( int b=1; b < h1->GetNbinsX()+1; ++b)
          h1->SetBinContent(b, analysisName.Contains("SMS") ? 1e-4 : 1e-2);
      }
      
      h1->Write();

      if( err_uncorr==0. ) continue;

      for( int iBin=1; iBin<h1->GetNbinsX()+1; ++iBin ) {

        float binContent = h1->GetBinContent(iBin);

        float thisErrUncorr = (err_uncorr>0.) ? err_uncorr : h1->GetBinError(iBin)/binContent;

        TH1D* h1_binUp = new TH1D(*h1);
        h1_binUp->SetName(Form("%s_bin_%dUp", h1->GetName(), iBin));
        h1_binUp->SetBinContent( iBin, binContent*( 1. + thisErrUncorr ) );
        h1_binUp->Write();

        TH1D* h1_binDown = new TH1D(*h1);
        h1_binDown->SetName(Form("%s_bin_%dDown", h1->GetName(), iBin));
        h1_binDown->SetBinContent( iBin, binContent/( 1. + thisErrUncorr ) );
        h1_binDown->Write();

      } // for bins

  } // for regions

}



void writeToTemplateFile_poisson( TFile* file, MT2Analysis<MT2Estimate>* analysis, const std::string& name ) {

  file->cd();

  std::set<MT2Region> regions = analysis->getRegions();
  
  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
  

      TH1D* h1 = (TH1D*) (analysis->get( *iR )->yield->Clone());
      std::string oldName(h1->GetName());
      h1->SetName(Form("%s_%s", oldName.c_str(), name.c_str()));

      h1->Write();

      int nBJetsMin = iR->nBJetsMin();
      if( nBJetsMin>=2 ) continue;

      //////CHANGE HERE for ITERATION 1
      //Fake uncertainty
      //float k = (nBJetsMin==0) ? 2. : 20.;

      //Real uncertainty
      float k = 1.;

      for( int iBin=1; iBin<h1->GetNbinsX()+1; ++iBin ) {

        float binContent = h1->GetBinContent(iBin);
        int N_zinv = (int)binContent;
        float error = (N_zinv>0) ? 1./sqrt(k*N_zinv) : 0.;

        TH1D* h1_binUp = new TH1D(*h1);
        h1_binUp->SetName(Form("%s_bin_%dUp", h1->GetName(), iBin));
        h1_binUp->SetBinContent( iBin, binContent*( 1. + error ) );
        h1_binUp->SetLineColor(kGreen);
        h1_binUp->Write();

        TH1D* h1_binDown = new TH1D(*h1);
        h1_binDown->SetName(Form("%s_bin_%dDown", h1->GetName(), iBin));
        h1_binDown->SetBinContent( iBin, binContent/( 1. + error ) );
        h1_binDown->SetLineColor(kRed);
        h1_binDown->Write();

      } // for bins


      h1->SetName(oldName.c_str());

  } // for regions

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



void getQCDestimate( float htMin, float mt2Min, float mt2Max, float nB, int& NQCD_cr, float& r ){


  TF1* QCDpow = new TF1("QCDpow","[0]*TMath::Power(x,[1])");

  Int_t NQCD[6];

  int eLow = (mt2Min-200.)/100.;
  int eHigh = (mt2Max>0) ? (mt2Max-200.)/100. : 5;
  if( eLow  > 5 ) eLow =5;
  if( eHigh > 5 ) eHigh=5;

  if(htMin==450){

    QCDpow->SetParameter(0,73792.1);  QCDpow->SetParameter(1,-2.53584);
    Int_t NQCD_[6]={131, 9, 1, 0, 0, 0};
    for(int n=0; n<6; n++)
      NQCD[n]=NQCD_[n];

  }
  else if(htMin==575){

    QCDpow->SetParameter(0,57785.8);  QCDpow->SetParameter(1,-2.49123);
    Int_t NQCD_[6]={449, 49, 6, 1, 0, 0};
    for(int n=0; n<6; n++)
      NQCD[n]=NQCD_[n];


  }
  else if(htMin==1000){

    QCDpow->SetParameter(0,1958.03);  QCDpow->SetParameter(1,-1.70589);
    Int_t NQCD_[6]={192, 34, 7, 2, 0, 0};
    for(int n=0; n<6; n++)
      NQCD[n]=NQCD_[n];

  }
  else if(htMin==1500){
 
    QCDpow->SetParameter(0,385.588);  QCDpow->SetParameter(1,-1.30116);
    Int_t NQCD_[6]={93, 22, 7, 2, 1, 0};
    for(int n=0; n<6; n++)
      NQCD[n]=NQCD_[n];
    
  }
  
  NQCD_cr=0;
  for(int e=eLow; e<eHigh; ++e)
    NQCD_cr+=NQCD[e];

  r = QCDpow->Eval(mt2Min);

  if (nB==1)
    r*=0.16;
  else if (nB>=2)
    r*=0.03;
    

}
