#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip> 

#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TList.h"
#include "TObject.h"
#include "TString.h"
#include "RooHistError.h"

#include "interface/MT2Config.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateSyst.h"


bool use_gamma = true;
bool use_purity = true;


int round(float d) {
  return (int)(floor(d + 0.5));
}

//double lumi = 4;


void writeToTemplateFile( TFile* file, MT2Analysis<MT2Estimate>* analysis, float err_uncorr );
void writeToTemplateFile_poisson( TFile* file, MT2Analysis<MT2Estimate>* analysis, const std::string& name="stat" );
MT2Analysis<MT2Estimate>* get( const std::string& name, std::vector< MT2Analysis<MT2Estimate>* > analyses, const std::string& name1, const std::string& name2="", const std::string& name3="", const std::string& name4="" );
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


  if( argc!=2 ) {
    std::cout << "USAGE: ./createDatacards [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);



  std::string dir = cfg.getEventYieldDir();
  std::string mc_fileName = dir + "/analyses.root";

  std::string samplesName = "PHYS14_v5_skimprune";
  std::string regionsName = "zurich";


  bool useMC_qcd  = true;
  bool useMC_zinv = false;
  bool useMC_llep = true;

  float err_qcd_corr    = 0.0;
  float err_qcd_uncorr  = 1.0; // 100% of QCD MC yield
  float err_llep_corr   = 0.;
  //float err_llep_uncorr = 0.075;
  float err_llep_shape = 0.075;
  float err_llep_lepEff = 0.15;

  // float err_zinv_corr   = 0.0671; //  added in quadrature syst on templates (2%) and on f (4%) and MC stat on Rzg (5%) -> sqrt( 2*2 + 4*4 +5*5 ) = 6.708
  float err_zinv_corr   = 0.05; //  added in quadrature syst on templates (2%) and on f (4%) -> sqrt( 2*2 + 4*4  ) approx = 5

  //float err_zinv_corr   = 0.21; // 20% on Z/gamma ratio plus added in quadrature syst on templates (2%) and on f (4%) and MC stat on Rzg (5%) -> sqrt( 20*20 + 2*2 + 4*4 +5*5 ) = 21
  float err_zinv_uncorr = -1.; // will take histogram bin error
  float err_zinv_alpha_extra  = 0.2; // 20% extra uncertainty on alpha if using lower MT2 as CR
  float err_zinv_uncorr_2b = 1.0;
  float err_sig_corr    = 0.1;
  float err_sig_uncorr  = 0.;

  float llep_weighted   = 0.0314378*cfg.lumi();
  float zinv_weighted   = 0.0011701*cfg.lumi();

  MT2Analysis<MT2Estimate>* data  = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "data" );
  MT2Analysis<MT2Estimate>* qcd;
  if( useMC_qcd )
    qcd = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "QCD"  );
  else
    qcd = MT2Analysis<MT2Estimate>::readFromFile( "MT2QCDEstimate.root" );
  qcd->setName("qcd");


  
  MT2Analysis<MT2Estimate>* zinv;
  MT2Analysis<MT2Estimate>* zinvCR;
  MT2Analysis<MT2Estimate>* zinv_ratio;
  MT2Analysis<MT2EstimateSyst>* purity;


  MT2Analysis<MT2Estimate>* zll;
  MT2Analysis<MT2Estimate>* zll_mt2;
  MT2Analysis<MT2Estimate>* zll_yield;
  MT2Analysis<MT2Estimate>* zll_ht;
  MT2Analysis<MT2Estimate>* zll_nJets;
  MT2Analysis<MT2Estimate>* zll_nBJets;


  if( useMC_zinv )
    zinv = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "ZJets");
  else {
    zinvCR      = MT2Analysis<MT2Estimate>    ::readFromFile( dir + "/gammaControlRegion/data.root", "gammaCR");
    if( use_purity ){
      zinvCR      = MT2Analysis<MT2Estimate>    ::readFromFile( Form("GammaControlRegion_%s_%s_%.0ffb/data.root", samplesName.c_str(), regionsName.c_str(), cfg.lumi()), "gammaCR");
      zinv        = MT2Analysis<MT2Estimate>    ::readFromFile( Form("ZinvEstimateFromGamma_%s_%s_%.0ffb_type1/MT2ZinvEstimate.root", samplesName.c_str(), regionsName.c_str(), cfg.lumi()), "ZinvEstimate");
      zinv_ratio  = MT2Analysis<MT2Estimate>    ::readFromFile( Form("ZinvEstimateFromGamma_%s_%s_%.0ffb_type1/MT2ZinvEstimate.root", samplesName.c_str(), regionsName.c_str(), cfg.lumi()), "ZgammaRatio");
      purity      = MT2Analysis<MT2EstimateSyst>::readFromFile( Form("ZinvEstimateFromGamma_%s_%s_%.0ffb_type1/MT2ZinvEstimate.root", samplesName.c_str(), regionsName.c_str(), cfg.lumi()), "purity");
 
      zll      = MT2Analysis<MT2Estimate>::readFromFile( Form("ZllGamma_Ratio_%s_%.0ffb/zll_ratio.root", samplesName.c_str(), cfg.lumi()), "zllY_mt2");
      zll_mt2  = MT2Analysis<MT2Estimate>    ::readFromFile( Form("ZllGamma_Ratio_%s_%.0ffb/zll_ratio.root", samplesName.c_str(), cfg.lumi()), "zllG_mt2");
      zll_yield  = MT2Analysis<MT2Estimate>    ::readFromFile( Form("ZllGamma_Ratio_%s_%.0ffb/zll_ratio.root", samplesName.c_str(), cfg.lumi()), "zllY_mt2");
      zll_ht  = MT2Analysis<MT2Estimate>    ::readFromFile( Form("ZllGamma_Ratio_%s_%.0ffb/zll_ratio.root", samplesName.c_str(), cfg.lumi()), "zllG_ht");
      zll_nJets  = MT2Analysis<MT2Estimate>    ::readFromFile( Form("ZllGamma_Ratio_%s_%.0ffb/zll_ratio.root", samplesName.c_str(), cfg.lumi()), "zllG_nJets");
      zll_nBJets  = MT2Analysis<MT2Estimate>    ::readFromFile( Form("ZllGamma_Ratio_%s_%.0ffb/zll_ratio.root", samplesName.c_str(), cfg.lumi()), "zllG_nBJets");

   }   else{
      zinv        = MT2Analysis<MT2Estimate>    ::readFromFile( dir + "/zinvFromGamma_noPurity.root", "ZinvEstimate");
      zinv_ratio  = MT2Analysis<MT2Estimate>    ::readFromFile( dir + "/zinvFromGamma_noPurity.root", "ZgammaRatio");
    }
  }
  zinv->setName("zinv");
  zinv->addToFile( mc_fileName, true );


  MT2Analysis<MT2Estimate>* llep;
  if( useMC_llep ) {
    MT2Analysis<MT2Estimate>* wjets = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "WJets");
    MT2Analysis<MT2Estimate>* top   = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "Top");
    llep = new MT2Analysis<MT2Estimate>( (*wjets) + (*top) );
  } else {
    llep = MT2Analysis<MT2Estimate>::readFromFile( cfg.getEventYieldDir() + "/llepEstimate.root" );
  }
  llep->setName( "llep" );
  llep->addToFile( mc_fileName, true );

  //MT2Analysis<MT2Estimate>* llepCR = MT2Analysis<MT2Estimate>::readFromFile( Form("llep_%s_%s_llep_%.0ffb.root", samplesName.c_str(), regionsName.c_str(), lumi) );
  MT2Analysis<MT2Estimate>* llepCR = MT2Analysis<MT2Estimate>::readFromFile( cfg.getEventYieldDir() + "/llepEstimate.root" );



  std::set<MT2Region> regions = data->getRegions();


  //Zll yields (inclusive)
  std::set<MT2Region> inclRegions=  zll_ht->getRegions();
  MT2Region inclusiveRegion( (*inclRegions.begin() ) );

  TH1D* this_zll_ht = zll_ht->get(inclusiveRegion)->yield;
  TH1D* this_zll_yield = zll_yield->get(inclusiveRegion)->yield;
  TH1D* this_zll_nJets = zll_nJets->get(inclusiveRegion)->yield;
  TH1D* this_zll_nBJets = zll_nBJets->get(inclusiveRegion)->yield;
  TH1D* this_zll_mt2 = zll_mt2->get(inclusiveRegion)->yield;

  //Zll yield
  TH1D* this_zll = zll->get(inclusiveRegion)->yield;
 


  // first create template datacards

  std::string path_templ = dir + "/datacard_templates";
  system(Form("mkdir -p %s", path_templ.c_str()));

  int emptyZllBins = 0; //"empty" == (yield <5)

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

  

    TH1D* this_data = data->get(*iR)->yield;
    TH1D* this_qcd  = qcd ->get(*iR)->yield;
    TH1D* this_zinv = zinv->get(*iR)->yield;
    TH1D* this_zinvCR     = (use_gamma) ? zinvCR->get(*iR)->yield : 0;
    TH1D* this_zinv_ratio     = (use_gamma) ? zinv_ratio->get(*iR)->yield : 0;
    TH1D* this_llep = llep->get(*iR)->yield;
    TH1D* this_llepCR;
    //    if(iR->nJetsMin()>=7 && iR->nBJetsMin()>=1)
    //      this_llepCR = llepCR->get(MT2Region(iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), 1, 2))->yield;
    //    else  //commented out for mario
    this_llepCR = llepCR->get(*iR)->yield;
     
    TGraphAsymmErrors* this_zinv_purity;
    if ( use_purity ) this_zinv_purity = (use_gamma) ? purity->get(*iR)->getGraph() : 0;
    //TGraphAsymmErrors* this_zinv_purity = (use_gamma) ? zinv_purity->get(*iR)->getGraph() : 0;

    float N_llep_CR = this_llepCR->Integral();
    std::string llepCR_name;
    //     if(iR->nJetsMin()>=7 && iR->nBJetsMin()>=1){
    //       MT2Region* thisCR = new MT2Region(iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), 1, 2);
    //       llepCR_name = thisCR->getName();
    //     }
    //     else
    llepCR_name = iR->getName();

    //     if( iR->mtCut()!="" ) { 
    //       std::string choppedName = llepCR_name.substr(0, llepCR_name.size()-5);
    //       llepCR_name = choppedName;
    //       if( iR->mtCut()=="loMT" ) {
    //         N_llep_CR += llepCR->get(MT2Region(iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), iR->nBJetsMin(), iR->nBJetsMax(), "hiMT"))->yield->Integral();
    //       } else {
    //         N_llep_CR += llepCR->get(MT2Region(iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), iR->nBJetsMin(), iR->nBJetsMax(), "loMT"))->yield->Integral();
    //       }
    //     }

  
    
       
     unsigned iEmptyZinvBin=this_data->GetNbinsX()+1;
     int nEmptyCR=0;

     for( int iBin=1; iBin<this_data->GetNbinsX()+1; ++iBin ) {
       
       if(this_data->GetBinLowEdge( iBin ) > iR->htMax() && iR->htMax()>0 ) continue;

       float mt2Min = this_data->GetBinLowEdge( iBin );
       float mt2Max = (iBin==this_data->GetNbinsX()) ?  -1. : this_data->GetBinLowEdge( iBin+1 );

       std::string binName;
       if( mt2Max>=0. )
         binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
       else
         binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );


       std::string datacardName( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
       ofstream datacard( datacardName.c_str() );

       
       std::string tableName( Form("%s/table_%s.txt", path_templ.c_str(), binName.c_str()) );
       ofstream table( tableName.c_str() );
       table << std::setprecision(3);
       

       datacard << "imax 1" << std::endl;
       datacard << "jmax 3" << std::endl;
       datacard << "kmax *" << std::endl;
       //datacard << "shapes * * FAKE" << std::endl; // To fix instabilities of RooFit in case of too many bins
       datacard << "-------------" << std::endl;
       datacard << std::endl << std::endl;


       datacard << std::fixed;
       datacard << std::setprecision(3) << std::endl << std::endl;
       datacard << "bin  " << binName<< std::endl;
       datacard << "observation  " << round(this_data->GetBinContent(iBin)) << std::endl;
       datacard << "-------------" << std::endl;
       datacard << std::endl << std::endl;


       float yield_llep = this_llep->GetBinContent(iBin);
       float yield_qcd  = this_qcd->GetBinContent(iBin);
       float yield_zinv = this_zinv->GetBinContent(iBin);

       if( use_gamma ) {
             
         //if( yield_zinv<0.001 ) yield_zinv = 0.;
         if( yield_zinv<0.001 ) yield_zinv = zinv_weighted;
         if( yield_qcd <0.001 ) yield_qcd  = 0.;
         //if( yield_llep<0.001 ) yield_llep = 0.;
         if( yield_llep<0.001 ) yield_llep = llep_weighted;
         if( yield_llep<0.001 && yield_zinv<0.001 && yield_qcd<0.001 ) {
	   yield_qcd  = 0.01;
         }

       } else {
 
         if( yield_zinv<0.001 ) yield_zinv = 0.01;
         if( yield_qcd <0.001 ) yield_qcd  = 0.01;
         if( yield_llep<0.001 ) yield_llep = 0.01;

       }



       // sig qcd zinv llep
       datacard << "bin \t" << binName << "\t" << binName << "\t" << binName << "\t" << binName << std::endl;
       datacard << "process \t sig \t zinv \t llep \t qcd" << std::endl;
       datacard << "process \t 0 \t 1 \t 2 \t 3" << std::endl;
       datacard << "rate \t XXX";
       datacard << " \t " << yield_zinv << " \t " << yield_llep << " \t " << yield_qcd << std::endl;
       datacard << "-------------" << std::endl;

       datacard << "sig_syst    lnN    " << 1.+err_sig_corr << " - - -" << std::endl;



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






       // Z INVISIBLE SYSTEMATICS:

       if( yield_zinv>0. || use_gamma ) {

         if( iR->nBJetsMin()<2 ) { // 0 and 1 btag

           // correlated:
           datacard << "zinv_ZGratio lnN   - " << 1.+err_zinv_corr << " - -" << std::endl;
	   zinv_systUp += err_zinv_corr*err_zinv_corr;
           zinv_systDn += err_zinv_corr*err_zinv_corr;


	   //PROJECTION METHOD//////////////////////////////////////
	   int bin_mt2_zll = mt2Min/100 - 1;
	   if(bin_mt2_zll == 7) bin_mt2_zll = 6; //change this if you ever change the binning T.T
	   if(bin_mt2_zll == 9) bin_mt2_zll = 7;
	   if(bin_mt2_zll < 1) bin_mt2_zll = 7;

	   if(this_zll_yield->GetBinContent(bin_mt2_zll) < 5.){
	     //Low yield -> do uncorrelated uncertainties
	     datacard << "zinv_mt2_" << binName << " lnN  - " << 2 << " - -" << std::endl;
	     std::cout << "EMPTY BIN at " << binName << std::endl;
	     emptyZllBins+=1;
	   }else{
	     float zll_mt2 = 1+ this_zll_mt2->GetBinError(bin_mt2_zll) ; //change back here
	     datacard << "zll_mt2_"<< int(mt2Min)  << " lnN - " << zll_mt2 << " - -" << std::endl; 
	     //	     datacard << "zll_mt2_"<< bin_mt2_zll << " lnN - " << zll_mt2 << " - -" << std::endl; 
	   }
	   



	   if(this_zll_yield->GetBinContent(bin_mt2_zll) < 5.) {
	     //doing nothing, already taken care of by 100% uncorr uncertainty (nope I don't code elegantly
	   }else{
	     if(iR->htMax()<600){ 
	       float zll_ht = 1+ this_zll_ht->GetBinError(1) ;
	       datacard << "zll_ht_"<< int(iR->htMin()) << " lnN - " << zll_ht << " - -" << std::endl;
	  
	     }else if(iR->htMax()<1001){
	       float zll_ht = 1+ this_zll_ht->GetBinError(2) ;
	       datacard << "zll_ht_"<< int(iR->htMin())  << " lnN - " << zll_ht << " - -" << std::endl;
	 
	     }else if(iR->htMax()<1501 ){
	       float zll_ht = 1+ this_zll_ht->GetBinError(3);
	       datacard << "zll_ht_"<<  int(iR->htMin()) << " lnN - " << zll_ht << " - -" << std::endl;
	     } else{
	       float zll_ht = 1+ this_zll_ht->GetBinError(4) ;
	       datacard << "zll_ht_"<< int(iR->htMin())  << " lnN - " << zll_ht << " - -" << std::endl;	  
	     }

	     if(iR->nJetsMax()==3){
	       float zll_nJets = 1+ this_zll_nJets->GetBinError(1) ;
	       datacard << "zll_nJets_"<< iR->nJetsMin() << " lnN - " << zll_nJets << " - -" << std::endl;
	   
	     }else if(iR->nJetsMax()==6){
	       float zll_nJets = 1+ this_zll_nJets->GetBinError(2) ;
	       datacard << "zll_nJets_"<< iR->nJetsMin() << " lnN - " << zll_nJets << " - -" << std::endl;
	 
	     }else{
	       float zll_nJets = 1+ this_zll_nJets->GetBinError(3);
	       datacard << "zll_nJets_"<< iR->nJetsMin() << " lnN - " << zll_nJets << " - -" << std::endl;  
	     }

	     if(iR->nBJetsMax()==0){
	       float zll_nBJets = 1+ this_zll_nBJets->GetBinError(1);
	       datacard << "zll_nBJets_"<< iR->nBJetsMin() << " lnN - " << zll_nBJets << " - -" << std::endl;
	  
	     }else if(iR->nBJetsMax()==1){
	       float zll_nBJets = 1+ this_zll_nBJets->GetBinError(2) ;
	       datacard << "zll_nBJets_"<< iR->nBJetsMin()  << " lnN - " << zll_nBJets << " - -" << std::endl;
	  
	     }else if(iR->nBJetsMax()==2){
	       float zll_nBJets = 1+ this_zll_nBJets->GetBinError(3);
	       datacard << "zll_nBJets_"<< 3 << " lnN - " << zll_nBJets << " - -" << std::endl;
	     }else{
	       float zll_nBJets = 1+ this_zll_nBJets->GetBinError(4);
	       datacard << "zll_nBJets_"<< 4 << " lnN - " << zll_nBJets << " - -" << std::endl;	   
	     }
	   }//end of if yield <5 statement


         }

         // uncorrelated:
	 //  float thisError_zinv_uncorr = 1. + this_zinv->GetBinError(iBin)/yield_zinv;
        float thisError_zinv_uncorr_rel = this_zinv->GetBinError(iBin)/yield_zinv;
       
         if( !use_gamma ) {

           std::string iname = (iR->nBJetsMin()<2) ? "CRstat" : "MC";
           datacard << "zinv_" << iname << "_" << binName << " lnN - " << thisError_zinv_uncorr_rel << " - -" << std::endl;
	   zinv_systUp += thisError_zinv_uncorr_rel*thisError_zinv_uncorr_rel;
           zinv_systDn += thisError_zinv_uncorr_rel*thisError_zinv_uncorr_rel;


         } else {

           if( iR->nBJetsMin()>=2 ) {

	     // if( yield_zinv>0. )
	     // datacard << "zinv_MC_" << binName << " lnN - " << thisError_zinv_uncorr << " - -" << std::endl;
	     // else
	     datacard << "zinv_MC_" << binName << " lnN - " << 1.+err_zinv_uncorr_2b << " - -" << std::endl;
             zinv_systUp += err_zinv_uncorr_2b*err_zinv_uncorr_2b;
             zinv_systDn += err_zinv_uncorr_2b*err_zinv_uncorr_2b;

           } else {

             int Ngamma = round(this_zinvCR->GetBinContent(iBin));

             Double_t x_tmp, p, p_errUp, p_errDown;
	     if( use_purity ){

	       this_zinv_purity->GetPoint( iBin-1, x_tmp, p);
	       p_errUp   = this_zinv_purity->GetErrorYhigh( iBin -1 );
	       p_errDown = this_zinv_purity->GetErrorYlow ( iBin -1 ); 

	       if( Ngamma>0 ) {
		 datacard << "zinv_purity_" << binName << " lnN  - " << 1.+p_errUp/p << "/" << 1.-p_errDown/p << " - -" << std::endl;
		 zinv_systUp += (p_errUp/p)*(p_errUp/p);
		 zinv_systDn += (p_errDown/p)*(p_errDown/p);
	       }
	       
	     }//end of if use purity
             


	     float R = this_zinv_ratio->GetBinContent(iBin);
	     if( use_purity ) {
	       datacard << "zinv_CRstat_" << binName << " gmN " << Ngamma << " - " << R*p*0.92 << " - -" << std::endl;
	       double yield_zinv_up, yield_zinv_dn;
	       RooHistError::instance().getPoissonInterval(Ngamma,yield_zinv_dn,yield_zinv_up,1.);
	       yield_zinv_up *= R*p*0.92;
	       yield_zinv_dn *= R*p*0.92;
	       zinv_statUp = yield_zinv_up-yield_zinv;
	       zinv_statDn = yield_zinv-yield_zinv_dn;
	     } else {
	       datacard << "zinv_CRstat_" << binName << " gmN " << Ngamma << " - " << R << " - -" << std::endl;
	     }

	     float alphaErr = this_zinv_ratio->GetBinError(iBin)/R;
	     datacard << "zinv_alphaErr_" << binName << " lnN  - " << 1.+alphaErr << " - -" << std::endl;
	     zinv_systUp += alphaErr*alphaErr;
	     zinv_systDn += alphaErr*alphaErr;



/*
	     bool isEmptyCR=false;
             int Ngamma = round(this_zinvCR->GetBinContent(iBin));
	     if( Ngamma == 0 ){
	       if( nEmptyCR == 0 )
		 iEmptyZinvBin=iBin;
	       ++nEmptyCR;
	       Ngamma = round(this_zinvCR->GetBinContent(iEmptyZinvBin-1));
	       
	       isEmptyCR=true;
	       
	     }
             datacard << "zinv_CRstat_" << gammaConvention( yield_zinv, Ngamma, 1, binName ) << std::endl;
             if(  yield_zinv>0. && this_zinv_ratio->GetBinContent(iBin) > 0.) {
               float alphaErr = 1. + this_zinv_ratio->GetBinError(iBin)/this_zinv_ratio->GetBinContent(iBin); 
               datacard << "zinv_alphaErr_" << binName << " lnN  - " << alphaErr << " - -" << std::endl;
             }
	     if ( isEmptyCR ){
	       float alphaErr_extra = 1. + err_zinv_alpha_extra;
	       datacard << "zinv_alphaErr_extra lnN  - " << alphaErr_extra << " - -" << std::endl;
	     }
	     else ;
	     
*/
           } // if nbjets >= 2

         } // if use gamma

       } // if zinv




       // LOST LEPTON SYSTEMATICS:

       if( yield_llep>0. || use_gamma ) {

         // correlated within the SR (stat-like):
         float llep_stat_err = (N_llep_CR>0) ? 1./sqrt((float)N_llep_CR) : 0.;
         float llep_tot_err = sqrt( llep_stat_err*llep_stat_err + err_llep_lepEff*err_llep_lepEff );
         llep_tot_err+=1.;


         if( !use_gamma ) {

           datacard << "llep_CRstat_" << llepCR_name << "  lnN   - - " << llep_tot_err << " -" << std::endl;
           //datacard << "llep_shape_" << binName << " lnN - - " << 1.+err_llep_uncorr << " - " << std::endl;
	   datacard << "llep_shape_" << llepCR_name << " lnN - - " << 1.+err_llep_shape << " - " << std::endl;

         } else {

           datacard << "llep_lepeff_" << llepCR_name << "  lnN  - - " << 1.+err_llep_lepEff << " -" << std::endl;
	   //datacard << "llep_CRstat_" << gammaConvention( yield_llep, round(N_llep_CR), 2, llepCR_name, binName ) << std::endl;
           datacard << "llep_CRstat_" << gammaConvention( yield_llep, round(N_llep_CR), 2, llepCR_name ) << std::endl;
           if( yield_llep>0. ) {
             datacard << "llep_MCstat_" << binName << " lnN  - - " << 1.+this_llep->GetBinError(iBin)/yield_llep << " -" << std::endl;
             //datacard << "llep_shape_" << binName << " lnN - - " << 1.+err_llep_uncorr << " - " << std::endl;
	     datacard << "llep_shape_" << llepCR_name << " lnN - - " << 1.+err_llep_shape << " - " << std::endl;
           }

         }
           


       }



       // QCD SYSTEMATICS:

       if( yield_qcd>0. ) {
         datacard << "qcd_syst_" << binName << " lnN - - - " << 1.+err_qcd_uncorr << std::endl;
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
       table.close();

       std::cout << "-> Created BG table: " << tableName << std::endl;
       

    } // for bins

  } // for regions

  std::cout << "EMPTY ZLL BINS AMOUNT TOOO" << emptyZllBins << std::endl;




  // now create datacards for all signals
  //std::vector<MT2Analysis<MT2Estimate>*> signals = MT2Analysis<MT2Estimate>::readAllFromFile( mc_fileName, "SMS" );
  std::vector<MT2Analysis<MT2Estimate>*> signals = MT2Analysis<MT2Estimate>::readAllFromFile( mc_fileName, "DarkMatter" );

  for( unsigned  isig=0; isig<signals.size(); ++isig ) { 

    std::string sigName;
    if( signals[isig]->getName().find("fullScan") != std::string::npos )
      sigName = signals[isig]->getName();
    else if( signals[isig]->getName().find("DarkMatter") != std::string::npos )
      sigName = signals[isig]->getName();
    else
      sigName = getSimpleSignalName( signals[isig]->getName() );

    //    std::string sigName = getSimpleSignalName( signals[isig]->getName() );
    
    std::string path = dir + "/datacards_" + sigName;
    system(Form("mkdir -p %s", path.c_str()));

    std::string path_mass = path;
    float xs_norm=1.;
    if( signals[isig]->getName().find("T2qq") != std::string::npos )
      xs_norm=8./10.;

    if( signals[isig]->getName().find("fullScan") != std::string::npos )
      for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
	
	TH3D* this_signal3d = signals[isig]->get(*iR)->yield3d;
	
	TH1D* this_signalParent = this_signal3d->ProjectionY("mParent");
	
	for( int iBinY=1; iBinY<this_signalParent->GetNbinsX()+1; ++iBinY ){
	  
	  TH1D* this_signalLSP = this_signal3d->ProjectionZ("mLSP", 0, -1, iBinY, iBinY);
	  
	  for( int iBinZ=1; iBinZ < iBinY; ++iBinZ ) {
	  	  
	    float mParent = this_signalParent->GetBinLowEdge(iBinY);
	    float mLSP = this_signalLSP->GetBinLowEdge(iBinZ);
	    
	    TH1D* this_signal = this_signal3d->ProjectionX("mt2", iBinY, iBinY, iBinZ, iBinZ);
	    
	    if( this_signal->Integral() < 1e-3 ) continue;
	    
	    for( int iBin=1; iBin<this_signal->GetNbinsX()+1; ++iBin ) {
	      
	      if( this_signal->GetBinLowEdge( iBin ) > iR->htMax() && iR->htMax()>0 ) continue;
	      
	      float mt2Min = this_signal->GetBinLowEdge( iBin );
	      float mt2Max = (iBin==this_signal->GetNbinsX()) ?  -1. : this_signal->GetBinLowEdge( iBin+1 );
	      
	      if( this_signal->GetBinContent(iBin) < 1e-3 );
	      else{
		
		std::string binName;
		if( mt2Max>=0. )
		  binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
		else
		  binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );
		
		std::string templateDatacard( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
		
		std::string newDatacard( Form("%s/datacard_%s_%s_%.0f_%.0f.txt", path_mass.c_str(), binName.c_str(), sigName.c_str(), mParent, mLSP) );
		
		float sig = this_signal->GetBinContent(iBin);
		sig*=xs_norm;
		
		std::string sedCommand( Form("sed 's/XXX/%.3f/g' %s > %s", sig, templateDatacard.c_str(), newDatacard.c_str()) );
		system( sedCommand.c_str() );
		
	      }
	      
	    } // for bins X (MT2)
	  } // for bins Z (mLSP)
	}// for bins Y (mParent)      
      } // for regions
    
    else
      for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
	
	TH1D* this_signal = signals[isig]->get(*iR)->yield;
	
	for( int iBin=1; iBin<this_signal->GetNbinsX()+1; ++iBin ) {
	  
	  if( this_signal->GetBinLowEdge( iBin ) > iR->htMax() && iR->htMax()>0 ) continue;
	  
	  float mt2Min = this_signal->GetBinLowEdge( iBin );
	  float mt2Max = (iBin==this_signal->GetNbinsX()) ?  -1. : this_signal->GetBinLowEdge( iBin+1 );
	  
	  if( this_signal->GetBinContent(iBin) < 1e-3 );
	  else{
	    
	    std::string binName;
	    if( mt2Max>=0. )
	      binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
	    else
	      binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );
	    
	    std::string templateDatacard( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
	    
	    std::string newDatacard( Form("%s/datacard_%s_%s.txt", path.c_str(), binName.c_str(), sigName.c_str()) );
	    
	    float sig = this_signal->GetBinContent(iBin);
	    sig*=xs_norm;

	    std::string sedCommand( Form("sed 's/XXX/%.3f/g' %s > %s", sig, templateDatacard.c_str(), newDatacard.c_str()) );
	    system( sedCommand.c_str() );
	    
	  }
	  
	} // for bins X (MT2)
      } // for regions
    
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
  longName_tstr.ReplaceAll( "mGl", " " );
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
  // [2]: 2J
  // [3]: parent mass
  // [4]: lsp mass


  std::string simpleName = parts[1] + "_" + parts[3] + "_" + parts[4];

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
