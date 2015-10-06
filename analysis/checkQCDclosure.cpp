#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateQCD.h"
#include "../interface/MT2DrawTools.h"





int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|             Running checkQCDclosure                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  MT2DrawTools::setStyle();

  if( argc<2 ) {
    std::cout << "USAGE: ./checkQCDclosure [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  std::string mc_or_data = cfg.gammaTemplateType();
  if( argc>2 ) {

    mc_or_data = std::string(argv[2]); 
    std::cout << std::endl;
    std::cout << "-> Will disobey the cfg and use mc_or_data = " << argv[2] << std::endl;
    std::cout << std::endl;

  } 

  if( mc_or_data=="mc" ) mc_or_data="MC";


  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion";
  MT2Analysis<MT2EstimateQCD>* qcd = MT2Analysis<MT2EstimateQCD>::readFromFile( qcdCRdir + "/mcFits.root", "mc" );

  std::set<MT2Region> regions = qcd->getRegions();

  std::ofstream ofs_problem("qcdProblematic.txt");
  std::ofstream ofs_ok     ("qcdOK.txt");

  TH1D* h1_pull = new TH1D("pull", "", 100, -10., 10.);


  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateQCD* thisQCD = qcd->get( *iR );
   
    TH1D* thisRatio = thisQCD->ratio;    
    TF1* thisFit = thisQCD->exp;

    Double_t expPar = thisFit->GetParameter(0);
    Double_t expParErr = thisFit->GetParError(0);

    TF1* thisFit_errUp = new TF1(*thisFit);
    thisFit_errUp->SetParameter(0, expPar + expParErr);
    TF1* thisFit_errDown = new TF1(*thisFit);
    thisFit_errDown->SetParameter(0, expPar - expParErr);

    Double_t xMin = 150.;
    Double_t xMax = thisRatio->GetXaxis()->GetXmax();

    int binMin = thisRatio->FindBin(xMin);
    int binMax = thisRatio->FindBin(xMax);

    Double_t intErr_ratio;
    Double_t int_ratio = thisRatio->IntegralAndError( binMin, binMax, intErr_ratio );
    Double_t int_f1    = thisFit->Integral( xMin, xMax );
    Double_t int_f1_up = thisFit_errUp->Integral( xMin, xMax );
    Double_t int_f1_down = thisFit_errDown->Integral( xMin, xMax );
    Double_t int_f1_errUp = int_f1_up-int_f1;
    Double_t int_f1_errDown = int_f1-int_f1_down;
    Double_t intErr_f1 = TMath::Max( int_f1_errUp, int_f1_errDown );
 


    float zTest = (int_ratio-int_f1) / ( sqrt( intErr_ratio*intErr_ratio + intErr_f1*intErr_f1 ) );
    if( fabs(zTest) > 3. || int_f1/int_ratio > 10. ) {
      ofs_problem << iR->getName() << " ratio: " << int_ratio << " +- " << intErr_ratio << "     f1: " << int_f1 << " +- " << intErr_f1 << "  (z=" << zTest << ")" << std::endl;
    } else {
      ofs_ok      << iR->getName() << " ratio: " << int_ratio << " +- " << intErr_ratio << "     f1: " << int_f1 << " +- " << intErr_f1 << "  (z=" << zTest << ")" << std::endl;
    }

    h1_pull->Fill( zTest );

  }

  TCanvas* c1 = new TCanvas("c3", "", 600, 600 );
  c1->cd();
  h1_pull->SetXTitle("(MC - fit) / (#sigma_{MC} #oplus #sigma_{fit})");
  h1_pull->Draw();
  c1->SaveAs("pull.eps");

  std::cout << h1_pull->GetMean() << " +- " << h1_pull->GetMeanError() << std::endl;
  std::cout << h1_pull->GetRMS()  << " +- " << h1_pull->GetRMSError() << std::endl;

  ofs_problem.close();
  ofs_ok     .close();

  return 0;

}


