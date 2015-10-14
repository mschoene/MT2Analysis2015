#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TVirtualFitter.h"
#include "TLegend.h"


#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateQCD.h"
#include "../interface/MT2EstimateTree.h"
#include "../interface/MT2DrawTools.h"







void fillFromTreeAndRatio( MT2Estimate* estimate, MT2Estimate* r_effective, TTree* tree, TF1* f1_ratio );
void drawSingleFit( const std::string& outdir, TF1* thisFitQCD, MT2EstimateQCD* qcd, MT2EstimateQCD* all, float xMin_fit, float xMax_fit );



int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|              Running fitDeltaPhiQCD                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./fitDeltaPhiQCD [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  MT2DrawTools::setStyle();


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool useMC = true;

  if( argc>2 ) {

    std::string mc_or_data_templates = std::string(argv[2]); 
    if( mc_or_data_templates=="mc" ) mc_or_data_templates="MC";
    std::cout << std::endl;
    std::cout << "-> Will disobey the cfg and use mc_or_data_templates = " << mc_or_data_templates << std::endl;
    cfg.set_gammaTemplateType(mc_or_data_templates);
    if( mc_or_data_templates=="MC" ) useMC = true;
    else useMC=false;
    std::cout << std::endl;

  } 




  // always start from inclusive qcd tree:
  MT2Analysis<MT2EstimateTree>* qcdTree_mc = MT2Analysis<MT2EstimateTree>::readFromFile( "EventYields_data_Run2015D_25nsGolden_v4/qcdControlRegion/mc.root", "qcdCRtree" );
  //MT2Analysis<MT2EstimateTree>* qcdTree_mc = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/mc.root", "qcdCRtree" );


  MT2Analysis<MT2Estimate>* estimate     = new MT2Analysis<MT2Estimate>("qcdEstimate", cfg.regionsSet());
  MT2Analysis<MT2Estimate>* r_effective  = new MT2Analysis<MT2Estimate>("r_effective", cfg.regionsSet());
  //MT2Analysis<MT2Estimate>* r_hat        = new MT2Analysis<MT2Estimate>("r_hat"      , cfg.regionsSet());
  //MT2Analysis<MT2Estimate>* f_jets       = new MT2Analysis<MT2Estimate>("f_jets"     , cfg.regionsSet());


  TH1D::AddDirectory(kTRUE);

  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion/";
  std::string outputdir = qcdCRdir;
  std::string fitsDir = qcdCRdir + "/fits";
  system( Form("mkdir -p %s", fitsDir.c_str() ));

  MT2Analysis<MT2EstimateQCD>* mc_all  = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "mc"     , cfg.qcdRegionsSet(), qcdTree_mc, "" );
  MT2Analysis<MT2EstimateQCD>* qcdOnly = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "qcdOnly", cfg.qcdRegionsSet(), qcdTree_mc, "id>=100 && id<200" );



  std::set<MT2Region> regions = estimate->getRegions();


  std::string lastRegionName = "";

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* this_estimate     = estimate    ->get( *iR );
    MT2Estimate* this_r_effective  = r_effective ->get( *iR );
    //MT2Estimate* this_r_hat        = r_hat       ->get( *iR );
    //MT2Estimate* this_f_jets       = f_jets      ->get( *iR );

    MT2Region* matchedRegion = qcdOnly->matchRegion(*iR);
    MT2EstimateQCD* matchedEstimate_all = mc_all->get( *matchedRegion );
    MT2EstimateQCD* matchedEstimate_qcd = qcdOnly->get( *matchedRegion );


    float xMin_fit = (iR->htMin()==1500.) ? 70. : 60.;
    float xMax_fit = 100.;
    TF1* f1_ratio = matchedEstimate_qcd->getFit( "pow", xMin_fit, xMax_fit );

    if( matchedRegion->getName()!=lastRegionName ) // draw only one per HT region
      drawSingleFit( fitsDir, f1_ratio, matchedEstimate_qcd, matchedEstimate_all, xMin_fit, xMax_fit );

    lastRegionName = matchedRegion->getName();

    fillFromTreeAndRatio( this_estimate, this_r_effective, matchedEstimate_qcd->tree, f1_ratio );

  }  // for regions
      
    

  mc_all ->writeToFile( outputdir + "/mcFits.root", "recreate" );
  qcdOnly->writeToFile( outputdir + "/mcFits.root" );


  estimate ->writeToFile( outputdir + "/qcdEstimate.root", "recreate" );
  r_effective ->writeToFile( outputdir + "/qcdEstimate.root" );
  //r_hat ->writeToFile( outputdir + "/qcdEstimate.root" );
  //f_jets ->writeToFile( outputdir + "/qcdEstimate.root" );

  return 0;

}




void fillFromTreeAndRatio( MT2Estimate* estimate, MT2Estimate* r_effective, TTree* tree, TF1* f1_ratio ) {

  int nBins;
  double* bins;
  estimate->region->getBins(nBins, bins);

  TProfile* rprof = new TProfile( "rprof", "", nBins, bins );

  float weight;
  tree->SetBranchAddress( "weight", &weight );
  float mt2;
  tree->SetBranchAddress( "mt2", &mt2 );

  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    float r = f1_ratio->Eval( mt2 );
    estimate->yield->Fill( mt2, weight*r );

    rprof->Fill( mt2, r, weight );

  } // for entries

  for( int iBin=1; iBin<nBins+1; ++iBin ) {
    r_effective->yield->SetBinContent( iBin, rprof->GetBinContent(iBin) );
    r_effective->yield->SetBinError( iBin, rprof->GetBinError(iBin) );
  }

  delete rprof;

}



void drawSingleFit( const std::string& outdir, TF1* thisFitQCD, MT2EstimateQCD* qcd, MT2EstimateQCD* all, float xMin_fit, float xMax_fit ) {


  TH1D* thisRatioAll = all->getRatio();
  TH1D* thisRatioQCD = qcd->getRatio();

  TCanvas* c1 = new TCanvas( "c2", "", 600, 600 );
  c1->cd();
  c1->SetLogx();
  c1->SetLogy();

  
  float xMin = thisRatioAll->GetXaxis()->GetXmin();
  float xMax = thisRatioAll->GetXaxis()->GetXmax();

  float yMax    = thisRatioAll->GetMaximum()*5.;
  float yMinAll = thisRatioQCD->GetMinimum()/2.;
  float yMin    = thisRatioQCD->GetMinimum()/2.;
  if( yMin < 0.03 ) yMin = 0.03;
  if( yMin > yMinAll && yMinAll>0.001 ) yMin = yMinAll;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax );
  h2_axes->SetXTitle( "M_{T2} [GeV]" );
  h2_axes->SetYTitle( "Ratio");
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->GetXaxis()->SetMoreLogLabels();
  h2_axes->GetYaxis()->SetNoExponent();
  h2_axes->Draw();

  std::vector<std::string> regionNiceNames = qcd->region->getNiceNames();

  TPaveText* regionName = new TPaveText( 0.4, 0.81, 0.9, 0.9, "brNDC" );
  regionName->SetTextAlign( 11 );
  regionName->SetTextSize( 0.035 );
  regionName->SetFillColor( 0 );
  regionName->AddText( regionNiceNames[0].c_str() );
  regionName->AddText( regionNiceNames[1].c_str() );
  regionName->Draw("same");

  TLegend* legend = new TLegend( 0.4, 0.65, 0.8, 0.82 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( thisRatioAll, "MC (all)", "P" );
  legend->AddEntry( thisRatioQCD, "MC (QCD Only)", "P" );
  legend->AddEntry( thisFitQCD, "Fit", "L" );
  legend->Draw("same");

  TLine* lineLeft = new TLine( xMin_fit, yMin, xMin_fit, yMax );
  lineLeft->SetLineStyle(2);
  lineLeft->Draw("same");

  TLine* lineRight = new TLine( 100., yMin, 100., yMax );
  lineRight->SetLineStyle(2);
  lineRight->Draw("same");

  TH1D* h_band = new TH1D(Form("band_%s", thisFitQCD->GetName()) , "", 500, xMin, xMax);
  h_band->SetMarkerSize(0);
  h_band->SetFillColor(18); 
  h_band->SetFillStyle(3001);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(h_band, 0.68);

  h_band->Draw("C E3 same");

  thisFitQCD->SetLineColor(46); 
  thisFitQCD->SetLineWidth(2); 
  thisFitQCD->Draw("L same");

  thisRatioAll->SetMarkerStyle(20);
  thisRatioAll->SetMarkerSize(1.3);
  thisRatioAll->SetLineColor(kBlack);
  thisRatioAll->Draw("P same");

  thisRatioQCD->SetMarkerStyle(24);
  thisRatioQCD->SetMarkerSize(1.3);
  thisRatioQCD->SetLineColor(kBlack);
  thisRatioQCD->Draw("P same");

  gPad->RedrawAxis();

  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
  labelTop->Draw("same");

  c1->SaveAs( Form("%s/ratio_%s.eps", outdir.c_str(), qcd->region->getName().c_str()) );
  c1->SaveAs( Form("%s/ratio_%s.pdf", outdir.c_str(), qcd->region->getName().c_str()) );
  

  delete c1;
  delete h2_axes;
  delete h_band;

}
