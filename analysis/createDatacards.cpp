#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip> 

#include "TFile.h"
#include "TH1D.h"
#include "TList.h"
#include "TObject.h"
#include "TString.h"

#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"


//bool use_gamma = false;
bool use_gamma = true;


int round(float d) {
  return (int)(floor(d + 0.5));
}


void writeToTemplateFile( TFile* file, MT2Analysis<MT2Estimate>* analysis, float err_uncorr );
void writeToTemplateFile_poisson( TFile* file, MT2Analysis<MT2Estimate>* analysis, const std::string& name="stat" );
MT2Analysis<MT2Estimate>* get( const std::string& name, std::vector< MT2Analysis<MT2Estimate>* > analyses, const std::string& name1, const std::string& name2="", const std::string& name3="", const std::string& name4="" );
std::string getSimpleSignalName( const std::string& longName );
std::string gammaConvention( float yieldSR, int yieldCR, int position, const std::string& corrName, const std::string& uncorrName="", float testAlpha=1. );



int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./createDatacards [dir]" << std::endl;
    exit(113);
  } 


  std::string dir( argv[1] );
  std::string mc_fileName = dir + "/analyses.root";


  bool useMC_qcd  = true;
  bool useMC_zinv = false;
  bool useMC_llep = true;

  float err_qcd_corr    = 0.0;
  float err_qcd_uncorr  = 1.0; // 100% of QCD MC yield
  float err_llep_corr   = 0.;
  float err_llep_uncorr = 0.075;
  float err_llep_lepEff = 0.15;
  float err_zinv_corr   = 0.2; // 20% on Z/gamma ratio
  float err_zinv_uncorr = -1.; // will take histogram bin error
  float err_zinv_uncorr_2b = 1.0;
  float err_sig_corr    = 0.1;
  float err_sig_uncorr  = 0.;


  MT2Analysis<MT2Estimate>* data  = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "data" );
  MT2Analysis<MT2Estimate>* qcd;
  if( useMC_qcd )
    qcd = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "QCD"  );
  else
    qcd = MT2Analysis<MT2Estimate>::readFromFile( "MT2QCDEstimate.root" );
  qcd->setName("qcd");

  //////CHANGE HERE for ITERATION 1
  //MT2Analysis<MT2Estimate>* zinv  = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "ZJets");
  //zinv->setName("zinv");
  //MT2Analysis<MT2Estimate>* wjets = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "WJets");
  //MT2Analysis<MT2Estimate>* top   = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "Top"  );
  //MT2Analysis<MT2Estimate>* llep = new MT2Analysis<MT2Estimate>( *top + *wjets );
  //llep->setName( "llep" );
  ////*llep += *wjets;
  
  MT2Analysis<MT2Estimate>* zinv;
  MT2Analysis<MT2Estimate>* zinvCR;
  MT2Analysis<MT2Estimate>* zinv_ratio;
  if( useMC_zinv )
    zinv = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "ZJets");
  else {
    //zinv       = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "ZJets");
    //zinv       = MT2Analysis<MT2Estimate>::readFromFile( "ZinvEstimateFromGamma_PHYS14_v3_Zinv_13TeV_CSA14_4fb/MT2ZinvEstimate.root", "ZinvEstimate");
    zinv       = MT2Analysis<MT2Estimate>::readFromFile( "ZinvEstimateFromGamma_PHYS14_v3_Zinv_13TeV_PHYS14_hiHT_4fb/mc.root", "Zinv");
    zinvCR     = MT2Analysis<MT2Estimate>::readFromFile( "ZinvEstimateFromGamma_PHYS14_v3_Zinv_13TeV_PHYS14_hiHT_4fb/mc.root", "gammaJet");
    zinv_ratio = MT2Analysis<MT2Estimate>::readFromFile( "ZinvEstimateFromGamma_PHYS14_v3_Zinv_13TeV_PHYS14_hiHT_4fb/mc.root", "ZgammaRatio");
  }
  zinv->setName("zinv");
  zinv->addToFile( mc_fileName, true );


  MT2Analysis<MT2Estimate>* llep;
  if( useMC_llep ) {
    MT2Analysis<MT2Estimate>* wjets = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "WJets");
    MT2Analysis<MT2Estimate>* top   = MT2Analysis<MT2Estimate>::readFromFile( mc_fileName, "Top");
    llep = new MT2Analysis<MT2Estimate>( (*wjets) + (*top) );
  } else {
    //llep = MT2Analysis<MT2Estimate>::readFromFile( "llep_PHYS14_v3_CSA14.root" );
    llep = MT2Analysis<MT2Estimate>::readFromFile( "llep_PHYS14_v3_13TeV_PHYS14_hiHT.root" );
  }
  llep->setName( "llep" );
  llep->addToFile( mc_fileName, true );


  //MT2Analysis<MT2Estimate>* llepCR = MT2Analysis<MT2Estimate>::readFromFile( "llep_PHYS14_v3_CSA14.root" );
  MT2Analysis<MT2Estimate>* llepCR = MT2Analysis<MT2Estimate>::readFromFile( "llep_PHYS14_v3_13TeV_PHYS14_hiHT.root" );
  //MT2Analysis<MT2Estimate>* llepCR = llep;


  std::set<MT2Region> regions = data->getRegions();




  // first create template datacards

  std::string path_templ = dir + "/datacard_templates";
  system(Form("mkdir -p %s", path_templ.c_str()));

  
  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {


     TH1D* this_data = data->get(*iR)->yield;
     TH1D* this_qcd  = qcd ->get(*iR)->yield;
     TH1D* this_zinv = zinv->get(*iR)->yield;
     TH1D* this_llep = llep->get(*iR)->yield;
     TH1D* this_llepCR = llepCR->get(*iR)->yield;
     TH1D* this_zinvCR     = (use_gamma) ? zinvCR->get(*iR)->yield : 0;
     TH1D* this_zinv_ratio = (use_gamma) ? zinv_ratio->get(*iR)->yield : 0;

     float N_llep_CR = this_llepCR->Integral();
     std::string llepCR_name = iR->getName();
     if( iR->mtCut()!="" ) { 
       std::string choppedName = llepCR_name.substr(0, llepCR_name.size()-5);
       llepCR_name = choppedName;
       if( iR->mtCut()=="loMT" ) {
         N_llep_CR += llepCR->get(MT2Region(iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), iR->nBJetsMin(), iR->nBJetsMax(), "hiMT"))->yield->Integral();
       } else {
         N_llep_CR += llepCR->get(MT2Region(iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), iR->nBJetsMin(), iR->nBJetsMax(), "loMT"))->yield->Integral();
       }
     }
         


     for( unsigned iBin=1; iBin<this_data->GetNbinsX()+1; ++iBin ) {

       float mt2Min = this_data->GetBinLowEdge( iBin );
       float mt2Max = (iBin==this_data->GetNbinsX()) ?  -1. : this_data->GetBinLowEdge( iBin+1 );

       std::string binName;
       if( mt2Max>=0. )
         binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
       else
         binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );


       std::string datacardName( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
       ofstream datacard( datacardName.c_str() );


       datacard << "imax 1" << std::endl;
       datacard << "jmax 3" << std::endl;
       datacard << "kmax *" << std::endl;
       datacard << "-------------" << std::endl;
       datacard << std::endl << std::endl;


       datacard << std::fixed;
       datacard << std::setprecision(3) << std::endl << std::endl;
       datacard << "bin  " << binName<< std::endl;
       datacard << "observation  " << this_data->GetBinContent(iBin) << std::endl;
       datacard << "-------------" << std::endl;
       datacard << std::endl << std::endl;


       float yield_llep = this_llep->GetBinContent(iBin);
       float yield_qcd  = this_qcd->GetBinContent(iBin);
       float yield_zinv = this_zinv->GetBinContent(iBin);

       if( use_gamma ) {
             
         if( yield_zinv<0.001 ) yield_zinv = 0.;
         if( yield_qcd <0.001 ) yield_qcd  = 0.;
         if( yield_llep<0.001 ) yield_llep = 0.;
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






       // Z INVISIBLE SYSTEMATICS:

       if( yield_zinv>0. || use_gamma ) {

         if( iR->nBJetsMin()<2 ) { // 0 and 1 btag

           // correlated:
           datacard << "zinv_ZGratio lnN   - " << 1.+err_zinv_corr << " - -" << std::endl;

         }

         // uncorrelated:
         float thisError_zinv_uncorr = 1. + this_zinv->GetBinError(iBin)/yield_zinv;
         if( !use_gamma ) {

           std::string iname = (iR->nBJetsMin()<2) ? "CRstat" : "MC";
           datacard << "zinv_" << iname << "_" << binName << " lnN - " << thisError_zinv_uncorr << " - -" << std::endl;

         } else {

           if( iR->nBJetsMin()>=2 ) {

//             if( yield_zinv>0. )
//               datacard << "zinv_MC_" << binName << " lnN - " << thisError_zinv_uncorr << " - -" << std::endl;
//	     else
	     datacard << "zinv_MC_" << binName << " lnN - " << 1.+err_zinv_uncorr_2b << " - -" << std::endl;

           } else {

             int Ngamma = round(this_zinvCR->GetBinContent(iBin));
             datacard << "zinv_CRstat_" << gammaConvention( yield_zinv, Ngamma, 1, binName ) << std::endl;
             if(  yield_zinv>0. && this_zinv_ratio->GetBinContent(iBin) > 0.) {
               float alphaErr = 1. + this_zinv_ratio->GetBinError(iBin)/this_zinv_ratio->GetBinContent(iBin); 
               datacard << "zinv_alphaErr_" << binName << " lnN  - " << alphaErr << " - -" << std::endl;
             }

           }

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
           datacard << "llep_shape_" << binName << " lnN - - " << 1.+err_llep_uncorr << " - " << std::endl;

         } else {

           datacard << "llep_lepeff_" << llepCR_name << "  lnN  - - " << 1.+err_llep_lepEff << " -" << std::endl;
	   datacard << "llep_CRstat_" << gammaConvention( yield_llep, round(N_llep_CR), 2, llepCR_name, binName ) << std::endl;
           if( yield_llep>0. ) {
             datacard << "llep_MCstat_" << binName << " lnN  - - " << 1.+this_llep->GetBinError(iBin)/yield_llep << " -" << std::endl;
             datacard << "llep_shape_" << binName << " lnN - - " << 1.+err_llep_uncorr << " - " << std::endl;
           }

         }
           


       }



       // QCD SYSTEMATICS:

       if( yield_qcd>0. ) {
         datacard << "qcd_syst_" << binName << " lnN - - - " << 1.+err_qcd_uncorr << std::endl;
       }



       datacard.close();

       std::cout << "-> Created template datacard: " << datacardName << std::endl;

    } // for bins

  } // for regions





  // now create datacards for all signals
  std::vector<MT2Analysis<MT2Estimate>*> signals = MT2Analysis<MT2Estimate>::readAllFromFile( mc_fileName, "SMS" );


  for( unsigned  isig=0; isig<signals.size(); ++isig ) { 

    std::string sigName = getSimpleSignalName( signals[isig]->getName() );
    std::string path = dir + "/datacards_" + sigName;
    system(Form("mkdir -p %s", path.c_str()));


    for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

      TH1D* this_signal = signals[isig]->get(*iR)->yield;


      for( unsigned iBin=1; iBin<this_signal->GetNbinsX()+1; ++iBin ) {

        float mt2Min = this_signal->GetBinLowEdge( iBin );
        float mt2Max = (iBin==this_signal->GetNbinsX()) ?  -1. : this_signal->GetBinLowEdge( iBin+1 );

	//if( this_signal->GetBinContent(iBin) < 1e-3 );
	//if( this_signal->GetBinContent(iBin) < 1e-2 );
	//if( this_signal->GetBinContent(iBin) < 1e-3 || this_signal->GetBinContent(iBin) > 1e-2 );
	if( this_signal->GetBinContent(iBin) < 1e-1 );
	else{
	  
	std::string binName;
        if( mt2Max>=0. )
          binName = std::string( Form("%s_m%.0fto%.0f", iR->getName().c_str(), mt2Min, mt2Max) );
        else
          binName = std::string( Form("%s_m%.0ftoInf", iR->getName().c_str(), mt2Min) );
        
        std::string templateDatacard( Form("%s/datacard_%s.txt", path_templ.c_str(), binName.c_str()) );
        
        std::string newDatacard( Form("%s/datacard_%s_%s.txt", path.c_str(), binName.c_str(), sigName.c_str()) );
        
        
        float sig = this_signal->GetBinContent(iBin);
        
        std::string sedCommand( Form("sed 's/XXX/%.3f/g' %s > %s", sig, templateDatacard.c_str(), newDatacard.c_str()) );
        system( sedCommand.c_str() );

      }

      } // for bins

    } // for regions

    std::cout << "-> Created datacards in " << path << std::endl;
       
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

      for( unsigned iBin=1; iBin<h1->GetNbinsX()+1; ++iBin ) {

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

      for( unsigned iBin=1; iBin<h1->GetNbinsX()+1; ++iBin ) {

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

  for( unsigned i=0; i<position; ++i )
    line << " - ";

  line << syst;

  for( unsigned i=position+1; i<4; ++i )
    line << " - ";

  std::string line_str = line.str();
  return line_str;

}
