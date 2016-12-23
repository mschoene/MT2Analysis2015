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




void doAllShapePlots( const MT2Config& cfg, const std::string& mc_or_data, const std::string& shapeName, const std::string& var , const std::string& label= "M_{T2} (Photon Removed) [GeV]", const std::string& topoRegion = "");

void compareRegions( const std::string& outputdir, std::vector<MT2Region> regions, MT2Analysis<MT2EstimateSyst>* analysis, const std::string& suffix=""  );




int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|              Running drawShape                     |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  MT2DrawTools::setStyle();

  if( argc<2 ) {
    std::cout << "USAGE: ./drawShapeFits [configFileName] data/MC (optional: axes)" << std::endl;
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

  bool doAxes = 0;
  if( argc>3 ) {
    doAxes = true; 
    std::cout << std::endl;
    std::cout << "-> Will also draw the shape along the other binning axes == " << doAxes << std::endl;
    std::cout << std::endl;
  } 


  if( mc_or_data=="mc" ) mc_or_data="MC";


  // doAllShapePlots( cfg, mc_or_data, "shapeLoose", "","" ); 
  doAllShapePlots( cfg, mc_or_data, "shape", "", "" ); 


  //    doAllShapePlots( cfg, mc_or_data, "shape", "incl_nbjets", "b-Jet Multiplicity", "#geq1j, #geq0b" );


  return 0;

}





void doAllShapePlots( const MT2Config& cfg, const std::string& mc_or_data, const std::string& shapeName, const std::string& var, const std::string& label , const std::string& topoRegion ) {

  std::string score = "_";
  std::string variable = var;
  if(var != "")
    variable = score + var;


  std::string zllCRdir = cfg.getEventYieldDir();

  //MT2Analysis<MT2EstimateSyst>* data_shape = MT2Analysis<MT2EstimateSyst>::readFromFile( zllCRdir + "/zinvFromZll.root", "shape" );
  //MT2Analysis<MT2EstimateSyst>* MC_shape = MT2Analysis<MT2EstimateSyst>::readFromFile( zllCRdir + "/zinvFromZll.root", "zllMC_shape" );

  MT2Analysis<MT2EstimateSyst>* data_shape = MT2Analysis<MT2EstimateSyst>::readFromFile( zllCRdir + "/zinvFromZll.root", "shape_data_forPR" );
  // MT2Analysis<MT2EstimateSyst>* MC_shape = MT2Analysis<MT2EstimateSyst>::readFromFile( zllCRdir + "/zinvFromZll.root", "zllMC_shape_TR" );
  MT2Analysis<MT2EstimateSyst>* hybrid_shape = MT2Analysis<MT2EstimateSyst>::readFromFile( zllCRdir + "/zinvFromZll.root", "shape_hybrid_forPR" );

 MT2Analysis<MT2Estimate>* extrapol = MT2Analysis<MT2Estimate>::readFromFile( zllCRdir + "/zinvFromZll.root", "bin_extrapol_forPR");


  std::string outputdir = zllCRdir + "/zllControlRegion/shapes_forPR";
  system( Form("mkdir -p %s", outputdir.c_str() ));
 
  std::set<MT2Region> regions = data_shape->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    // if( iR->nBJetsMin()>1 ) continue;

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
   
    MT2EstimateSyst* thisShapeMC = data_shape->get( *iR );
 
    TH1D* h_extrapol = extrapol->get( *iR )->yield;
    int extrapolBin = h_extrapol->GetBinContent(1);

    TH1D* h_data_shape = data_shape->get( *iR )->yield;
    h_data_shape->SetLineColor( kBlack );
    h_data_shape->SetMarkerColor( kBlack );
    h_data_shape->SetMarkerStyle( 20 );
    h_data_shape->SetMarkerSize( 2 );
    h_data_shape->SetLineWidth( 2 );
    
    // TH1D* h_MC_shape = MC_shape->get( *iR )->yield;
    // h_MC_shape->SetLineColor( 38 );
    // h_MC_shape->SetMarkerColor( 38 );
    // h_MC_shape->SetMarkerStyle( 24 );
    // h_MC_shape->SetMarkerSize( 2 );
    // h_MC_shape->SetLineWidth( 2 );


    TH1D* h_hybrid_shape = hybrid_shape->get( *iR )->yield;
    h_hybrid_shape->SetLineColor( kBlack );
    //    h_hybrid_shape->SetFillStyle( 3244 );
    h_hybrid_shape->SetFillColor( 430 );
    h_hybrid_shape->SetMarkerColor( 430 );
    //   h_hybrid_shape->SetMarkerStyle( 25 );
    h_hybrid_shape->SetMarkerSize( 2 );
    h_hybrid_shape->SetLineWidth( 2 );

 

    TH1D* hybrid_shaded = (TH1D*)h_hybrid_shape->Clone("hybrid_shaded");
    hybrid_shaded->SetFillStyle( 3244 );
    hybrid_shaded->SetFillColor( kWhite );
    if(extrapolBin>1)
      hybrid_shaded->SetAxisRange( h_hybrid_shape->GetBinLowEdge( extrapolBin ), hybrid_shaded->GetBinLowEdge( -1 ),"X");
    else
          hybrid_shaded->SetFillColor( kWhite );
   
    // if(extrapolBin>1)
    //   hybrid_shaded->SetAxisRange( h_hybrid_shape->GetBinLowEdge(1), hybrid_shaded->GetBinLowEdge( extrapolBin-1 ),"X");
    // else
    //       hybrid_shaded->SetFillColor( kWhite );
   
    std::cout << hybrid_shaded->GetBinLowEdge( extrapolBin ) << std::endl;


    TH1D* gr_ratio = (TH1D*)h_data_shape->Clone("gr_ratio");
    gr_ratio->SetLineColor( kBlack );
    gr_ratio->SetMarkerColor( kBlack );
    gr_ratio->SetMarkerStyle( 20 );
    gr_ratio->SetMarkerSize( 1 );
    gr_ratio->SetLineWidth( 2 );


    std::string thisName_Band = "band";
    TH1D* h_band = (TH1D*)h_data_shape->Clone("band");
    h_band->SetMarkerSize(0);
    h_band->SetFillColor (kGray+2);
    h_band->SetFillStyle (3244);
    for ( int iBin=1; iBin <= h_hybrid_shape->GetNbinsX(); iBin++){
      h_band->SetBinContent(iBin,1);
      double error=0;
      if( h_hybrid_shape->GetBinContent(iBin)>0)
	error = h_hybrid_shape->GetBinError(iBin)/h_hybrid_shape->GetBinContent(iBin);
      h_band->SetBinError(iBin, error);
    }

    float yMin = 0.01;
 
    TPad* pad1 = MT2DrawTools::getCanvasMainPad();
    pad1->Draw();
    pad1->cd();

    TH2D* axes = new TH2D( "axes", "", 10, h_data_shape->GetXaxis()->GetXmin(), h_data_shape->GetXaxis()->GetXmax(), 10, yMin, h_data_shape->GetMaximum()*10. );
    // TH2D* axes = new TH2D( "axes", "", 10, h_data_shape->GetXaxis()->GetXmin(), h_data_shape->GetXaxis()->GetXmax(), 10, yMin, 0.02 );
    //    TH2D* axes = new TH2D( "axes", "", 10, thisShapeMC->yield->GetXaxis()->GetXmin(), thisShapeMC->yield->GetXaxis()->GetXmax(), 10, yMin, 0.02 );
    axes->SetXTitle( label.c_str());
    axes->SetYTitle( "Entries" );
    axes->Draw("");

    gPad->SetLogy();


   std::cout << extrapolBin << std::endl;

    for( int iBin=1; iBin<= h_hybrid_shape->GetNbinsX(); iBin++){
      std::cout << h_hybrid_shape->GetBinContent(iBin) << std::endl;
      std::cout << h_hybrid_shape->GetBinError(iBin) << std::endl;
      std::cout << h_hybrid_shape->GetBinError(iBin)/ h_hybrid_shape->GetBinContent(iBin) << std::endl;
      std::cout << std::endl;
    }




    //    MT2DrawTools::addLabels( c1, 40, "CMS Simulation" );
    MT2DrawTools::addLabels( (TCanvas*)pad1, cfg.lumi(), "CMS Preliminary" );

 
    //   h_MC_shape->Draw("p same");
    //    h_MC_shape->Draw("p same");
    h_hybrid_shape->Draw("histo E same");

    //    hybrid_shaded->Draw("hist same");
  
    h_data_shape->Draw("p same");

    //float xMin_legend = 0.6;
    float xMin_legend = 0.71;
    //   float xMin_legend = (mc_or_data=="MC" || fits.size()==1) ? 0.6 : 0.52;
    TLegend* legend = new TLegend( xMin_legend, 0.2+0.6, 0.89, 0.2+0.6+0.06*( 2.) );
    legend->SetTextSize(0.038); 
    legend->SetFillColor(0);
    legend->AddEntry( h_data_shape, "Data", "PL" );
    //  legend->AddEntry( h_MC_shape, "DY MC", "PL" );
    // legend->AddEntry( h_MC_shape, "DY MC", "PL" );
    legend->AddEntry( h_hybrid_shape, "Hybrid", "F" );
    
    std::vector<TGraphAsymmErrors*> graphs;

    legend->Draw("same");
    std::vector<std::string> regionNames = iR->getNiceNames();
    TPaveText* labelRegion = new TPaveText( 0.43, 0.79, 0.58, 0.89, "brNDC" );
    labelRegion->SetTextSize(0.034); 
    labelRegion->SetFillColor(0);
    for( unsigned i=0; i<regionNames.size(); ++i ) {
      //    for( unsigned i=0; i<regionNames.size(); ++i ) {
      if( topoRegion!="" && i==1){
	labelRegion->AddText( topoRegion.c_str() );
      }else{
	labelRegion->AddText( regionNames[i].c_str() );  }
    }


    //float x, y;
    double y= h_data_shape->GetMaximum();
    //    double y= gr_data_shape->GetY()[0];

    if(y>0.001)
      std::cout << iR->getName().c_str() << " " <<  y << std::endl;

    /*    for( unsigned i=0; i<regionNames.size(); ++i ) {
    //    for( unsigned i=0; i<regionNames.size(); ++i ) {
    if( topoRegion!="" && i==1){
    labelRegion->AddText( topoRegion.c_str() );
    }else{
    ; }//	labelRegion->AddText( regionNames[i].c_str() );  }
    }*/

    labelRegion->Draw("same");

    gPad->RedrawAxis();

    c1->cd();

    TPad* pad2 =  MT2DrawTools::getCanvasRatioPad();
    pad2->Draw();
    pad2->cd();
   

    TH2D* axes_ratio = MT2DrawTools::getRatioAxes( h_data_shape->GetXaxis()->GetXmin(), h_data_shape->GetXaxis()->GetXmax(), 0.0 , 2.0 );
    //   axes_ratio->SetXTitle( label.c_str() );
    axes_ratio->SetYTitle( "Data/Hybrid" );
    axes_ratio->Draw("");

    TLine* line_one = new TLine(h_data_shape->GetXaxis()->GetXmin(), 1, h_data_shape->GetXaxis()->GetXmax(), 1);

    line_one->Draw("same");
    //   gPad->SetLogy();
    h_band->Draw("E2same");

    gr_ratio->Divide( h_hybrid_shape );

    gr_ratio->Draw("p same");
    gPad->RedrawAxis();


    c1->SaveAs( Form("%s/Shape%s_%s.eps", outputdir.c_str(), variable.c_str(), iR->getName().c_str()) );
    c1->SaveAs( Form("%s/Shape%s_%s.png", outputdir.c_str(), variable.c_str(), iR->getName().c_str()) );
    c1->SaveAs( Form("%s/Shape%s_%s.pdf", outputdir.c_str(), variable.c_str(), iR->getName().c_str()) );

    delete c1;
    delete axes;

  } // for regions



  // TString zllCRdir_tstr(zllCRdir);
  
  // if( zllCRdir_tstr.Contains("zurich_onlyHT") ) {


  //   std::vector<MT2Region> r_vsHT;
  //   r_vsHT.push_back( MT2Region( 450., 575.  ) );
  //   r_vsHT.push_back( MT2Region( 575., 1000. ) );
  //   r_vsHT.push_back( MT2Region( 1000.,1500. ) );
  //   r_vsHT.push_back( MT2Region( 1500., -1.  ) );

  //   compareRegions( outputdir, r_vsHT, fits[0].shape );


  // } else if( zllCRdir_tstr.Contains("zurich_onlyJets") ) {


  //   std::vector<MT2Region> r_vs_njet;
  //   r_vs_njet.push_back( MT2Region( 450., -1., 2, 3, 0, 0 ) );
  //   r_vs_njet.push_back( MT2Region( 450., -1., 4, 6, 0, 0 ) );
  //   r_vs_njet.push_back( MT2Region( 450., -1., 2, 3, 1, 1 ) );
  //   r_vs_njet.push_back( MT2Region( 450., -1., 4, 6, 1, 1 ) );

  //   compareRegions( outputdir, r_vs_njet, fits[0].shape );



  // } else if( zllCRdir_tstr.Contains("zurichPlus") ) { //to be adjusted for new zurich2016 regions


  //   std::vector<MT2Region> r_lowHT_vs_njet;
  //   r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2, 3, 0, 0 ) );
  //   r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, 6, 0, 0 ) );
  //   r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 2, 3, 1, 1 ) );
  //   r_lowHT_vs_njet.push_back( MT2Region( 450., 575., 4, 6, 1, 1 ) );

  //   compareRegions( outputdir, r_lowHT_vs_njet, fits[0].shape );


  //   std::vector<MT2Region> r_medHT_vs_njet;
  //   r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2, 3, 0, 0 ) );
  //   r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, 6, 0, 0 ) );
  //   r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 2, 3, 1, 1 ) );
  //   r_medHT_vs_njet.push_back( MT2Region( 575., 1000., 4, 6, 1, 1 ) );

  //   compareRegions( outputdir, r_medHT_vs_njet, fits[0].shape );


  //   std::vector<MT2Region> r_highHT_vs_njet;
  //   r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 2, 3, 0, 0 ) );
  //   r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 4, 6, 0, 0 ) );
  //   r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 2, 3, 1, 1 ) );
  //   r_highHT_vs_njet.push_back( MT2Region( 1000., 1500., 4, 6, 1, 1 ) );

  //   compareRegions( outputdir, r_highHT_vs_njet, fits[0].shape );


  //   std::vector<MT2Region> r_veryhighHT_vs_njet;
  //   r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 2, 3, 0, 0 ) );
  //   r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 4, 6, 0, 0 ) );
  //   r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 2, 3, 1, 1 ) );
  //   r_veryhighHT_vs_njet.push_back( MT2Region( 1500., -1., 4, 6, 1, 1 ) );

  //   compareRegions( outputdir, r_veryhighHT_vs_njet, fits[0].shape );



  //   std::vector<MT2Region> r_vsHT;
  //   r_vsHT.push_back( MT2Region( 450., 575. , 2, 3, 0, 0 ) );
  //   r_vsHT.push_back( MT2Region( 575., 1000., 2, 3, 0, 0 ) );
  //   r_vsHT.push_back( MT2Region( 1000.,1500., 2, 3, 0, 0 ) );
  //   r_vsHT.push_back( MT2Region( 1500., -1. , 2, 3, 0, 0 ) );

  //   compareRegions( outputdir, r_vsHT, fits[0].shape );



  //   std::vector<MT2Region> r_vsHT_j46;
  //   r_vsHT_j46.push_back( MT2Region( 450., 575. , 4, 6, 0, 0 ) );
  //   r_vsHT_j46.push_back( MT2Region( 575., 1000., 4, 6, 0, 0 ) );
  //   r_vsHT_j46.push_back( MT2Region( 1000.,1500., 4, 6, 0, 0 ) );
  //   r_vsHT_j46.push_back( MT2Region( 1500., -1. , 4, 6, 0, 0 ) );

  //   compareRegions( outputdir, r_vsHT_j46, fits[0].shape );



  //   std::vector<MT2Region> r_vsHT_b1;
  //   r_vsHT_b1.push_back( MT2Region( 450., 575. , 2, 3, 1, 1 ) );
  //   r_vsHT_b1.push_back( MT2Region( 575., 1000., 2, 3, 1, 1 ) );
  //   r_vsHT_b1.push_back( MT2Region( 1000.,1500., 2, 3, 1, 1 ) );
  //   r_vsHT_b1.push_back( MT2Region( 1500., -1. , 2, 3, 1, 1 ) );

  //   compareRegions( outputdir, r_vsHT_b1, fits[0].shape );



  //   std::vector<MT2Region> r_vsHT_j46_b1;
  //   r_vsHT_j46_b1.push_back( MT2Region( 450., 575. , 4, 6, 1, 1 ) );
  //   r_vsHT_j46_b1.push_back( MT2Region( 575., 1000., 4, 6, 1, 1 ) );
  //   r_vsHT_j46_b1.push_back( MT2Region( 1000.,1500., 4, 6, 1, 1 ) );
  //   r_vsHT_j46_b1.push_back( MT2Region( 1500., -1. , 4, 6, 1, 1 ) );

  //   compareRegions( outputdir, r_vsHT_j46_b1, fits[0].shape );

  // }

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
  h2_axes->SetYTitle( "Shape (MC Fit)");
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
  
  TPaveText* labelTop = MT2DrawTools::getLabelTop("40 fb^{-1} (13TeV)");
  labelTop->Draw("same");

  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Simulation");
  labelCMS->Draw("same");
  
  gPad->RedrawAxis();

  std::string saveName = (loopOnHT) ? regions[0].sigRegion()->getName() : regions[0].htRegion()->getName();

  c1->SaveAs( Form("%s/%s_%s%s.eps", outputdir.c_str(), analysis->getName().c_str(), saveName.c_str(), suffix.c_str()) );
  c1->SaveAs( Form("%s/%s_%s%s.pdf", outputdir.c_str(), analysis->getName().c_str(), saveName.c_str(), suffix.c_str()) );
  c1->SaveAs( Form("%s/%s_%s%s.png", outputdir.c_str(), analysis->getName().c_str(), saveName.c_str(), suffix.c_str()) );

  delete c1;
  delete h2_axes;

}



