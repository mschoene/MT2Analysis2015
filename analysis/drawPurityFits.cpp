#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

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




void doAllPurityPlots( const std::string& gammaCRdir, const std::string& samples, const std::string& mc_or_data, const std::string& purityName );



int main( int argc, char* argv[] ) {


  MT2DrawTools::setStyle();

  std::string mc_or_data = "MC";
  if( argc>1 ) {
    mc_or_data = std::string(argv[1]);
    if( mc_or_data=="data" ) mc_or_data = "Data";
  }

  //std::string samples = "CSA14_Zinv";
  std::string samples = "PHYS14_v2_Zinv";

  std::string gammaCRdir = "GammaControlRegion_" + samples + "_zurich";
  
  doAllPurityPlots( gammaCRdir, samples, mc_or_data, "purityLoose" ); 
  doAllPurityPlots( gammaCRdir, samples, mc_or_data, "purity" ); 

  return 0;

}





void doAllPurityPlots( const std::string& gammaCRdir, const std::string& samples, const std::string& mc_or_data, const std::string& purityName ) {

  MT2Analysis<MT2EstimateSyst>* purityMC = MT2Analysis<MT2EstimateSyst>::readFromFile( gammaCRdir + "/purityMC.root", purityName );

  std::string outputdir = gammaCRdir + "/PurityFits" + mc_or_data;
  system( Form("mkdir -p %s", outputdir.c_str() ));


  std::vector< PurityFit > fits;
  if( mc_or_data=="MC" ) {
    //fits.push_back( PurityFit( "All Bins"  , "13TeV_CSA14"     , MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC_" + samples + "_13TeV_CSA14/purityFit_"     + samples + "_13TeV_CSA14.root"    , purityName), 20, kRed+2 ));
    //fits.push_back( PurityFit( "HT Bins"   , "13TeV_onlyHT"    , MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC_" + samples + "_13TeV_onlyHT/purityFit_"    + samples + "_13TeV_onlyHT.root"   , purityName), 21, 29 ));
    //fits.push_back( PurityFit( "Jet Bins"  , "13TeV_onlyJet"   , MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC_" + samples + "_13TeV_onlyJets/purityFit_"  + samples + "_13TeV_onlyJets.root" , purityName), 24, kAzure ));
    fits.push_back( PurityFit( "Inclusive" , "13TeV_inclusive" , MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC/purityFit.root", purityName), 25, kOrange+1 ));
  } else {
    fits.push_back( PurityFit( "Template Fit (MC)"  , "13TeV_inclusive" , MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsMC/purityFit.root", purityName), 21, 29 ));
    fits.push_back( PurityFit( "Template Fit (Data)", "13TeV_inclusive" , MT2Analysis<MT2EstimateSyst>::readFromFile(gammaCRdir+"/PurityFitsData/purityFit.root", purityName), 20, kOrange+1 ));
  }


  std::set<MT2Region> regions = purityMC->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()>1 ) continue;

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();


    MT2EstimateSyst* thisPurityMC = purityMC->get( *iR );

    TGraphAsymmErrors* gr_purityMC = thisPurityMC->getGraph();
    gr_purityMC->SetLineColor( kBlack );
    gr_purityMC->SetLineWidth( 2 );


    float yMin = (purityName=="purity") ? 0.7 : 0.;


    TH2D* axes = new TH2D( "axes", "", 10, thisPurityMC->yield->GetXaxis()->GetXmin(), thisPurityMC->yield->GetXaxis()->GetXmax(), 10, yMin, 1.0001 );
    axes->SetXTitle( "M_{T2} [GeV]");
    axes->SetYTitle( "Photon Purity" );
    axes->Draw("");

    TPaveText* labelTop = MT2DrawTools::getLabelTop();
    labelTop->Draw("same");


    gr_purityMC->Draw("p same");

    float xMin_legend = (mc_or_data=="MC") ? 0.65 : 0.52;
    TLegend* legend = new TLegend( xMin_legend, 0.2, 0.9, 0.2+0.06*(fits.size()+1.) );
    legend->SetTextSize(0.038); 
    legend->SetFillColor(0);
    legend->AddEntry( gr_purityMC, "MC Purity", "L" );
   
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


}
