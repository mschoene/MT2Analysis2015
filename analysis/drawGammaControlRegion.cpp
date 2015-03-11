#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2Estimate.h"
#include "../interface/MT2DrawTools.h"




void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2Estimate>* analysis );


int main( int argc, char* argv[] ) {


  std::string samples = "PHYS14_v2_Zinv";
  if( argc>1 ) {
    std::string samples_tmp(argv[1]); 
    samples = samples_tmp;
  }


  std::string regionsSet = "13TeV_CSA14";

  MT2DrawTools::setStyle();

  std::string outputdir = "Plots_GammaControlRegion_" + samples + "_" + regionsSet;
  system( Form("mkdir -p %s", outputdir.c_str()) );

  MT2Analysis<MT2Estimate>* purity_isoCut = MT2Analysis<MT2Estimate>::readFromFile("GammaControlRegion_" + samples + "_13TeV_CSA14/purityMC.root", "purity_isoCut");
  MT2Analysis<MT2Estimate>* eff_isoCut    = MT2Analysis<MT2Estimate>::readFromFile("GammaControlRegion_" + samples + "_13TeV_CSA14/purityMC.root", "eff_isoCut");

  purity_isoCut->setFullName( "Purity" );
  eff_isoCut->setFullName( "Efficiency" );


  std::vector<MT2Region> r_lowHT_vs_njet;
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2,  3, 0, 0 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, -1, 0, 0 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2,  3, 1, 1 ) );
  r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, -1, 1, 1 ) );

  compareRegions( outputdir, r_lowHT_vs_njet, purity_isoCut );
  compareRegions( outputdir, r_lowHT_vs_njet, eff_isoCut );

  std::vector<MT2Region> r_medHT_vs_njet;
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2,  3, 0, 0 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, -1, 0, 0 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2,  3, 1, 1 ) );
  r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, -1, 1, 1 ) );

  compareRegions( outputdir, r_medHT_vs_njet, purity_isoCut );
  compareRegions( outputdir, r_medHT_vs_njet, eff_isoCut );

  std::vector<MT2Region> r_highHT_vs_njet;
  r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 2,  3, 0, 0 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 4, -1, 0, 0 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 2,  3, 1, 1 ) );
  r_highHT_vs_njet.push_back( MT2Region( 1000., -1., 4, -1, 1, 1 ) );

  compareRegions( outputdir, r_highHT_vs_njet, purity_isoCut );
  compareRegions( outputdir, r_highHT_vs_njet, eff_isoCut );




  //std::vector<MT2Region> r_vsHT;
  //r_vsHT.push_back( MT2Region( 450., 575. , 2, 3, 0, 0 ) );
  //r_vsHT.push_back( MT2Region( 575., 1000., 2, 3, 0, 0 ) );
  //r_vsHT.push_back( MT2Region( 1000., -1. , 2, 3, 0, 0 ) );

  //compareRegions( outputdir, r_vsHT, purity_isoCut, true );
  //compareRegions( outputdir, r_vsHT, eff_isoCut , false );



  //std::vector<MT2Region> r_vsHT2;
  //r_vsHT2.push_back( MT2Region( 450., 575. , 4, -1, 0, 0 ) );
  //r_vsHT2.push_back( MT2Region( 575., 1000., 4, -1, 0, 0 ) );
  //r_vsHT2.push_back( MT2Region( 1000., -1. , 4, -1, 0, 0 ) );

  //compareRegions( outputdir, r_vsHT2, purity_isoCut, true );
  //compareRegions( outputdir, r_vsHT2, eff_isoCut , false );

  return 0;

}





void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2Estimate>* analysis ) {


  bool loopOnHT=true;

  if( regions.size()>1 ) 
    if( (*(regions[0].htRegion())) == (*(regions[1].htRegion())) ) loopOnHT=false;


  std::vector<int> colors;
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 38 );
  colors.push_back( 42 );
  colors.push_back( kRed );
  colors.push_back( kBlack );
  
  std::vector<int> markers;
  markers.push_back( 21 );
  markers.push_back( 20 );
  markers.push_back( 23 );
  markers.push_back( 24 );
  markers.push_back( 25 );
  markers.push_back( 26 );
  


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );


  TH2D* h2_axes = new TH2D( "axes", "", 10, 200., 1500., 10, 0.7, 1.0001 );
  h2_axes->SetYTitle( analysis->getFullName().c_str() );
  h2_axes->SetXTitle( "M_{T2} [GeV]" );
  c1->cd();
  h2_axes->Draw("");


  std::string legendTitle = (loopOnHT) ? regions[0].sigRegion()->getNiceName() : regions[0].htRegion()->getNiceName();

  TLegend* legend = new TLegend( 0.45, 0.2, 0.9, 0.2+0.08*regions.size(), legendTitle.c_str() );
  legend->SetTextSize(0.038);
  legend->SetFillColor(0);
  

  for( unsigned i=0; i<regions.size(); ++i ) {

    MT2Estimate* thisEstimate = analysis->get( regions[i] );
    if( thisEstimate==0 ) {
      std::cout << "ERROR! Didn't find estimate for region: " << regions[i].getName() << " ! Exiting." << std::endl;
      exit(119);
    }

    TH1D* thisYield = (TH1D*)(thisEstimate->yield->Clone());
    thisYield->SetLineColor(colors[i]);
    thisYield->SetLineWidth(2);

    thisYield->Draw("same" );

    if( loopOnHT )
      legend->AddEntry( thisYield, regions[i].htRegion()->getNiceName().c_str(), "L" );
    else
      legend->AddEntry( thisYield, regions[i].sigRegion()->getNiceName().c_str(), "L" );

  }



  TPaveText* labelTop = MT2DrawTools::getLabelTop();


  c1->cd();
  labelTop->Draw("same");
  legend->Draw("same");
  gPad->RedrawAxis();



  std::string saveName = (loopOnHT) ? regions[0].sigRegion()->getName() : regions[0].htRegion()->getName();

  c1->SaveAs( Form("%s/%s_%s.eps", outputdir.c_str(), analysis->getName().c_str(), saveName.c_str()) );
  c1->SaveAs( Form("%s/%s_%s.pdf", outputdir.c_str(), analysis->getName().c_str(), saveName.c_str()) );
  c1->SaveAs( Form("%s/%s_%s.png", outputdir.c_str(), analysis->getName().c_str(), saveName.c_str()) );


  delete c1;
  delete h2_axes;

}
