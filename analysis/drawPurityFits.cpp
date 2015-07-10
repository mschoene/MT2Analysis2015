#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateZinvGamma.h"
#include "../interface/MT2EstimateSyst.h"
#include "../interface/MT2DrawTools.h"







class PurityFit {

 public:
  PurityFit( const std::string& aNiceName, const std::string& aRegionsSet, MT2Analysis<MT2EstimateSyst>* apurity, int amarker, int acolor ) {
    niceName = aNiceName;
    regionsSet = aRegionsSet;
    purity = apurity;
    marker = amarker;
    color = acolor;
  }

  std::string niceName;
  std::string regionsSet;
  MT2Analysis<MT2EstimateSyst>* purity;
  int marker;
  int color;

 private:

};




void doAllPurityPlots( const MT2Config& cfg, const std::string& mc_or_data, const std::string& purityName );
void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2EstimateSyst>* analysis, const std::string& suffix="" );




int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|              Running drawPurityFits                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  MT2DrawTools::setStyle();

  if( argc<2 ) {
    std::cout << "USAGE: ./drawPurityFits [configFileName]" << std::endl;
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

  if( mc_or_data=="data" ) mc_or_data="DataRC";
  if( mc_or_data=="dataRC" ) mc_or_data="DataRC";
  if( mc_or_data=="mc" ) mc_or_data="MC";




  doAllPurityPlots( cfg, mc_or_data, "purityLoose" ); 
  doAllPurityPlots( cfg, mc_or_data, "purity" ); 

  return 0;

}





void doAllPurityPlots( const MT2Config& cfg, const std::string& mc_or_data, const std::string& purityName ) {

  std::string gammaCRdir = cfg.getEventYieldDir() + "/gammaControlRegion";

  MT2Analysis<MT2EstimateSyst>* purityMC = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaCRdir + "/purityMC.root", purityName );

  std::string outputdir = gammaCRdir + "/PurityFits" + mc_or_data;
  system( Form("mkdir -p %s", outputdir.c_str() ));


  std::vector< PurityFit > fits;
  if( mc_or_data=="MC" ) {
    //fits.push_back( PurityFit( "All Bins"  , "13TeV_CSA14"     , MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC_" + samples + "_13TeV_CSA14/purityFit_"     + samples + "_13TeV_CSA14.root"    , purityName), 20, kRed+2 ));
    //fits.push_back( PurityFit( "HT Bins"   , "13TeV_onlyHT"    , MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC_" + samples + "_13TeV_onlyHT/purityFit_"    + samples + "_13TeV_onlyHT.root"   , purityName), 21, 29 ));
    //fits.push_back( PurityFit( "Jet Bins"  , "13TeV_onlyJet"   , MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC_" + samples + "_13TeV_onlyJets/purityFit_"  + samples + "_13TeV_onlyJets.root" , purityName), 24, kAzure ));
    fits.push_back( PurityFit( "Template Fit" , cfg.gammaTemplateRegions(), MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC/purityFit.root", purityName), 25, kOrange+1 ));
  } else {
    if( cfg.lumi()<=1. ) {
      fits.push_back( PurityFit( "Template Fit", cfg.gammaTemplateRegions(), MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFits" + mc_or_data + "/purityFit.root", purityName), 20, kOrange+1 ));
    } else {
      fits.push_back( PurityFit( "Template Fit (MC)"  , cfg.gammaTemplateRegions(), MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC/purityFit.root", purityName), 21, 29 ));
      fits.push_back( PurityFit( "Template Fit (Data)", cfg.gammaTemplateRegions(), MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFits" + mc_or_data + "/purityFit.root", purityName), 20, kOrange+1 ));
    }
  }


  std::set<MT2Region> regions = purityMC->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()>1 ) continue;

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();


    MT2EstimateSyst* thisPurityMC = purityMC->get( *iR );

    TGraphAsymmErrors* gr_purityMC = thisPurityMC->getGraph();
    gr_purityMC->SetLineColor( kBlack );
    gr_purityMC->SetMarkerColor( kBlack );
    gr_purityMC->SetMarkerStyle( 24 );
    gr_purityMC->SetMarkerSize( 2 );
    gr_purityMC->SetLineWidth( 2 );


    float yMin = (purityName=="purity") ? 0.7 : 0.;


    TH2D* axes = new TH2D( "axes", "", 10, thisPurityMC->yield->GetXaxis()->GetXmin(), thisPurityMC->yield->GetXaxis()->GetXmax(), 10, yMin, 1.0001 );
    axes->SetXTitle( "M_{T2} [GeV]");
    axes->SetYTitle( "Photon Purity" );
    axes->Draw("");

    TPaveText* labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
    labelTop->Draw("same");


    gr_purityMC->Draw("p same");

    //float xMin_legend = 0.6;
    float xMin_legend = (mc_or_data=="MC" || fits.size()==1) ? 0.6 : 0.52;
    TLegend* legend = new TLegend( xMin_legend, 0.2, 0.9, 0.2+0.06*(fits.size()+1.) );
    legend->SetTextSize(0.038); 
    legend->SetFillColor(0);
    legend->AddEntry( gr_purityMC, "MC Purity", "PL" );
   
    for( unsigned i=0; i<fits.size(); ++i ) {

      MT2EstimateSyst* thisPurityFit = fits[i].purity->get( *iR );
      TGraphAsymmErrors* graph = thisPurityFit->getGraph();
      graph->SetMarkerStyle( fits[i].marker );
      graph->SetMarkerColor( fits[i].color );
      graph->SetLineColor( fits[i].color );
      graph->SetMarkerSize( 1.3 );
      graph->Draw("Psame");

      legend->AddEntry( graph, fits[i].niceName.c_str(), "P" );

    }

    legend->Draw("same");

    std::vector<std::string> regionNames = iR->getNiceNames();
    TPaveText* labelRegion = new TPaveText( 0.23, 0.18, 0.48, 0.29, "brNDC" );
    labelRegion->SetTextSize(0.034); 
    labelRegion->SetFillColor(0);
    for( unsigned i=0; i<regionNames.size(); ++i ) labelRegion->AddText( regionNames[i].c_str() );
    labelRegion->Draw("same");

    gPad->RedrawAxis();

    c1->SaveAs( Form("%s/fits_%s_%s.eps", outputdir.c_str(), purityName.c_str(), iR->getName().c_str()) );
    c1->SaveAs( Form("%s/fits_%s_%s.png", outputdir.c_str(), purityName.c_str(), iR->getName().c_str()) );
    c1->SaveAs( Form("%s/fits_%s_%s.pdf", outputdir.c_str(), purityName.c_str(), iR->getName().c_str()) );

    delete c1;
    delete axes;

  } // for regions



  TString gammaCRdir_tstr(gammaCRdir);
  
  if( gammaCRdir_tstr.Contains("zurich_onlyHT") ) {


    std::vector<MT2Region> r_vsHT;
    r_vsHT.push_back( MT2Region( 450., 575.  ) );
    r_vsHT.push_back( MT2Region( 575., 1000. ) );
    r_vsHT.push_back( MT2Region( 1000.,1500. ) );
    r_vsHT.push_back( MT2Region( 1500., -1.  ) );

    compareRegions( outputdir, r_vsHT, fits[0].purity );


  } else if( gammaCRdir_tstr.Contains("zurich_onlyJets") ) {


    std::vector<MT2Region> r_vs_njet;
    r_vs_njet.push_back( MT2Region( 450., -1., 2, 3, 0, 0 ) );
    r_vs_njet.push_back( MT2Region( 450., -1., 4, 6, 0, 0 ) );
    r_vs_njet.push_back( MT2Region( 450., -1., 2, 3, 1, 1 ) );
    r_vs_njet.push_back( MT2Region( 450., -1., 4, 6, 1, 1 ) );

    compareRegions( outputdir, r_vs_njet, fits[0].purity );



  } else if( gammaCRdir_tstr.Contains("zurich") ) {


    std::vector<MT2Region> r_lowHT_vs_njet;
    r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2, 3, 0, 0 ) );
    r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, 6, 0, 0 ) );
    r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2, 3, 1, 1 ) );
    r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, 6, 1, 1 ) );

    compareRegions( outputdir, r_lowHT_vs_njet, fits[0].purity );


    std::vector<MT2Region> r_medHT_vs_njet;
    r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2, 3, 0, 0 ) );
    r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, 6, 0, 0 ) );
    r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2, 3, 1, 1 ) );
    r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, 6, 1, 1 ) );

    compareRegions( outputdir, r_medHT_vs_njet, fits[0].purity );


    std::vector<MT2Region> r_highHT_vs_njet;
    r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 2, 3, 0, 0 ) );
    r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 4, 6, 0, 0 ) );
    r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 2, 3, 1, 1 ) );
    r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 4, 6, 1, 1 ) );

    compareRegions( outputdir, r_highHT_vs_njet, fits[0].purity );


    std::vector<MT2Region> r_veryhighHT_vs_njet;
    r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 2, 3, 0, 0 ) );
    r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 4, 6, 0, 0 ) );
    r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 2, 3, 1, 1 ) );
    r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 4, 6, 1, 1 ) );

    compareRegions( outputdir, r_veryhighHT_vs_njet, fits[0].purity );



    std::vector<MT2Region> r_vsHT;
    r_vsHT.push_back( MT2Region( 450., 575. , 2, 3, 0, 0 ) );
    r_vsHT.push_back( MT2Region( 575., 1000., 2, 3, 0, 0 ) );
    r_vsHT.push_back( MT2Region( 1000.,1500., 2, 3, 0, 0 ) );
    r_vsHT.push_back( MT2Region( 1500., -1. , 2, 3, 0, 0 ) );

    compareRegions( outputdir, r_vsHT, fits[0].purity );



    std::vector<MT2Region> r_vsHT_j46;
    r_vsHT_j46.push_back( MT2Region( 450., 575. , 4, 6, 0, 0 ) );
    r_vsHT_j46.push_back( MT2Region( 575., 1000., 4, 6, 0, 0 ) );
    r_vsHT_j46.push_back( MT2Region( 1000.,1500., 4, 6, 0, 0 ) );
    r_vsHT_j46.push_back( MT2Region( 1500., -1. , 4, 6, 0, 0 ) );

    compareRegions( outputdir, r_vsHT_j46, fits[0].purity );



    std::vector<MT2Region> r_vsHT_b1;
    r_vsHT_b1.push_back( MT2Region( 450., 575. , 2, 3, 1, 1 ) );
    r_vsHT_b1.push_back( MT2Region( 575., 1000., 2, 3, 1, 1 ) );
    r_vsHT_b1.push_back( MT2Region( 1000.,1500., 2, 3, 1, 1 ) );
    r_vsHT_b1.push_back( MT2Region( 1500., -1. , 2, 3, 1, 1 ) );

    compareRegions( outputdir, r_vsHT_b1, fits[0].purity );



    std::vector<MT2Region> r_vsHT_j46_b1;
    r_vsHT_j46_b1.push_back( MT2Region( 450., 575. , 4, 6, 1, 1 ) );
    r_vsHT_j46_b1.push_back( MT2Region( 575., 1000., 4, 6, 1, 1 ) );
    r_vsHT_j46_b1.push_back( MT2Region( 1000.,1500., 4, 6, 1, 1 ) );
    r_vsHT_j46_b1.push_back( MT2Region( 1500., -1. , 4, 6, 1, 1 ) );

    compareRegions( outputdir, r_vsHT_j46_b1, fits[0].purity );

  }

}




void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2EstimateSyst>* analysis, const std::string& suffix ) {

  bool loopOnHT=true;

  if( regions.size()>1 ) 
    if( (*(regions[0].htRegion())) == (*(regions[1].htRegion())) ) loopOnHT=false;



  std::vector<int> colors;
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 38 );
  colors.push_back( 42 );
  colors.push_back( kGray+1 );
  colors.push_back( kRed );
  
  std::vector<int> markers;
  markers.push_back( 21 );
  markers.push_back( 20 );
  markers.push_back( 25 );
  markers.push_back( 24 );
  markers.push_back( 26 );
  
  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D( "axes", "", 10, 200., 1500., 10, 0., 1.0001 );
  h2_axes->SetXTitle( "M_{T2} [GeV]" );
  h2_axes->SetYTitle( "Purity (MC Fit)");
  h2_axes->Draw("");
  

  std::string legendTitle = (loopOnHT) ? regions[0].sigRegion()->getNiceName() : regions[0].htRegion()->getNiceName();

  TLegend* legend = new TLegend( 0.5, 0.5-0.07*regions.size(), 0.9, 0.5, legendTitle.c_str() );
  legend->SetTextSize(0.038);
  legend->SetFillColor(0);
  

  for( unsigned i=0; i<regions.size(); ++i ) {

    MT2EstimateSyst* thisEstimate = analysis->get( regions[i] );
    if( thisEstimate==0 ) {
      std::cout << "-> Didn't find estimate for region: " << regions[i].getName() << " ! Skipping." << std::endl;
      delete c1;
      delete h2_axes;
      return;
    }

    TGraphAsymmErrors* thisGraph = thisEstimate->getGraph();

    thisGraph->SetMarkerSize(1.);
    thisGraph->SetMarkerStyle(markers[i]);
    thisGraph->SetMarkerColor(colors[i]);
    thisGraph->SetLineColor(kBlack);

    thisGraph->SetLineWidth(2);
    thisGraph->SetLineColor(colors[i]);

    thisGraph->Draw("p same" );

    if( loopOnHT )
      legend->AddEntry( thisGraph, regions[i].htRegion()->getNiceName().c_str(), "PL" );
    else
      legend->AddEntry( thisGraph, regions[i].sigRegion()->getNiceName().c_str(), "PL" );

  }


  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  std::string saveName = (loopOnHT) ? regions[0].sigRegion()->getName() : regions[0].htRegion()->getName();

  c1->SaveAs( Form("%s/%s_%s%s.eps", outputdir.c_str(), analysis->getName().c_str(), saveName.c_str(), suffix.c_str()) );
  c1->SaveAs( Form("%s/%s_%s%s.pdf", outputdir.c_str(), analysis->getName().c_str(), saveName.c_str(), suffix.c_str()) );
  c1->SaveAs( Form("%s/%s_%s%s.png", outputdir.c_str(), analysis->getName().c_str(), saveName.c_str(), suffix.c_str()) );

  delete c1;
  delete h2_axes;

}



