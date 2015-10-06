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
  MT2Analysis<MT2EstimateQCD>* qcd = MT2Analysis<MT2EstimateQCD>::readFromFile( qcdCRdir + "/mcFits.root", "qcdOnly" );
  MT2Analysis<MT2EstimateQCD>* all = MT2Analysis<MT2EstimateQCD>::readFromFile( qcdCRdir + "/mcFits.root", "mc" );

  std::string fitsdir = qcdCRdir + "/fits";
  system( Form("mkdir -p %s", fitsdir.c_str()) );

  std::set<MT2Region> regions = qcd->getRegions();

  std::ofstream ofs_problem("qcdProblematic.txt");
  std::ofstream ofs_ok     ("qcdOK.txt");

  TH1D* h1_pull = new TH1D("pull", "", 100, -10., 10.);


  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateQCD* thisQCD = qcd->get( *iR );
    MT2EstimateQCD* thisAll = all->get( *iR );

    thisQCD->doFit();
   
    TH1D* thisRatio = thisAll->ratio;    
    TH1D* thisRatioQCD = thisQCD->ratio;    
    //TF1* thisFit = thisAll->exp;
    TF1* thisFit = thisQCD->exp; // this should be true closure but doesn't work!

    Double_t xMin = thisRatio->GetBinLowEdge(thisRatio->FindBin(150.));
    Double_t xMax = thisRatio->GetXaxis()->GetXmax();

    int binMin = thisRatio->FindBin(xMin);
    int binMax = thisRatio->FindBin(xMax);

    Double_t intErr_ratio;
    Double_t int_ratio = thisRatio->IntegralAndError( binMin, binMax, intErr_ratio );
    Double_t int_f1    = thisFit->Integral( xMin, xMax );
    Double_t intErr_f1 = 0.; //thisFit->IntegralError( xMin, xMax );


    float zTest = (int_ratio-int_f1) / ( sqrt( intErr_ratio*intErr_ratio + intErr_f1*intErr_f1 ) );
    if( fabs(zTest) > 3. || int_f1/int_ratio > 10. ) {
      ofs_problem << iR->getName() << " ratio: " << int_ratio << " +- " << intErr_ratio << "     f1: " << int_f1 << " +- " << intErr_f1 << "  (z=" << zTest << ")" << std::endl;
    } else {
      ofs_ok      << iR->getName() << " ratio: " << int_ratio << " +- " << intErr_ratio << "     f1: " << int_f1 << " +- " << intErr_f1 << "  (z=" << zTest << ")" << std::endl;
    }

    h1_pull->Fill( zTest );

    
    TCanvas* c1 = new TCanvas( "c2", "", 600, 600 );
    c1->cd();
    c1->SetLogx();
    c1->SetLogy();


    float yMax = thisRatio->GetMaximum()*5.;
    float yMinAll = thisRatio->GetMinimum()/2.;
    float yMin = thisRatioQCD->GetMinimum()/2.;
    if( yMin < 0.03 ) yMin = 0.03;
    if( yMin > yMinAll && yMinAll>0.001 ) yMin = yMinAll;

    TH2D* h2_axes = new TH2D("axes", "", 10, 40., xMax, 10, yMin, yMax );
    h2_axes->SetXTitle( "M_{T2} [GeV]" );
    h2_axes->SetYTitle( "Ratio");
    h2_axes->GetXaxis()->SetNoExponent();
    h2_axes->GetXaxis()->SetMoreLogLabels();
    h2_axes->GetYaxis()->SetNoExponent();
    //h2_axes->GetYaxis()->SetMoreLogLabels();
    h2_axes->Draw();

    std::vector<std::string> regionNiceNames = iR->getNiceNames();

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
    legend->AddEntry( thisRatio, "MC (all)", "P" );
    legend->AddEntry( thisRatioQCD, "MC (QCD Only)", "P" );
    legend->AddEntry( thisFit, "Exp. Fit", "L" );
    legend->Draw("same");

    TLine* lineLeft = new TLine( 60., yMin, 60., yMax );
    lineLeft->SetLineStyle(2);
    lineLeft->Draw("same");

    TLine* lineRight = new TLine( 100., yMin, 100., yMax );
    lineRight->SetLineStyle(2);
    lineRight->Draw("same");

    TH1D* fit_band = MT2DrawTools::getBand(thisFit, Form("band_%s", thisFit->GetName()) );
    fit_band->Draw("C E3 same");

    thisFit->SetLineColor(46); 
    thisFit->SetLineWidth(2); 
    thisFit->Draw("L same");
 
    thisRatio->SetMarkerStyle(20);
    thisRatio->SetMarkerSize(1.3);
    thisRatio->Draw("P same");
 
    thisRatioQCD->SetMarkerStyle(24);
    thisRatioQCD->SetMarkerSize(1.3);
    thisRatioQCD->Draw("P same");
 
    gPad->RedrawAxis();

    TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
    labelTop->Draw("same");

    c1->SaveAs( Form("%s/ratio_%s.eps", fitsdir.c_str(), iR->getName().c_str()) );
    c1->SaveAs( Form("%s/ratio_%s.pdf", fitsdir.c_str(), iR->getName().c_str()) );

    delete c1;
    delete h2_axes;

  }

  TCanvas* c1 = new TCanvas("c3", "", 600, 600 );
  c1->cd();
  h1_pull->SetXTitle("(MC - fit) / (#sigma_{MC} #oplus #sigma_{fit})");
  h1_pull->Draw();
  c1->SaveAs( Form("%s/pull.eps", qcdCRdir.c_str()) );
  c1->SaveAs( Form("%s/pull.pdf", qcdCRdir.c_str()) );

  std::cout << h1_pull->GetMean() << " +- " << h1_pull->GetMeanError() << std::endl;
  std::cout << h1_pull->GetRMS()  << " +- " << h1_pull->GetRMSError() << std::endl;

  ofs_problem.close();
  ofs_ok     .close();

  return 0;

}


