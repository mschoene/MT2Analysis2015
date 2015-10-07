#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateQCD.h"
#include "../interface/MT2DrawTools.h"



//void randomizePoisson( MT2EstimateQCD* thisEstimate, float lumi, TRandom3 rand );
float drawSinglePlot( const std::string& dir, const MT2Region& region, float lumi, TH1D* histo, float mcValue, float mcErr=0. );


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
    std::cout << "USAGE: ./checkQCDclosure [configFileName] [toyLumi=0.]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  float lumiToCheck = 0.; //fb
  if( argc>2 ) {
    lumiToCheck = std::atof(argv[2]); 
  } 

  if( lumiToCheck>0. ) {
    std::cout << Form("-> Will run toyMC for lumi = %.0f fb-1", lumiToCheck) << std::endl;
  } 

  int niter = (lumiToCheck>0.) ? 301 : 1;
  //TRandom3 rand(13);


  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion";
  MT2Analysis<MT2EstimateQCD>* qcd = MT2Analysis<MT2EstimateQCD>::readFromFile( qcdCRdir + "/mcFits.root", "qcdOnly" );
  MT2Analysis<MT2EstimateQCD>* all = MT2Analysis<MT2EstimateQCD>::readFromFile( qcdCRdir + "/mcFits.root", "mc" );

  std::string fitsdir = qcdCRdir + "/fits";
  system( Form("mkdir -p %s", fitsdir.c_str()) );

  std::string toymcdir(Form("%s/toyMC_lumi%.0f_iter%d", qcdCRdir.c_str(), lumiToCheck*1000., niter-1));
  system( Form("mkdir -p %s", toymcdir.c_str()) );


  std::set<MT2Region> regions = qcd->getRegions();

  std::string fileName_problem(Form("%s/qcdProblematic.txt", toymcdir.c_str()) );
  std::string fileName_ok     (Form("%s/qcdOK.txt", toymcdir.c_str()) );
  std::ofstream ofs_problem(fileName_problem);
  std::ofstream ofs_ok     (fileName_ok);

  TH1D* h1_pull = new TH1D("pull", "", 100, -10., 10.);



  int iRegion = 0;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()>=2 ) continue; // who cares

    float prescale = 1.;
    if( iR->htMin() < 500. ) prescale = 180.;  // HLT_HT350 for HT>=450
    else if( iR->htMin() < 600. ) prescale = 60.; // HLT_HT_475 for HT >=575

    TH1D* h1_thisPull = new TH1D( Form("pull_%s", iR->getName().c_str()), "", 100, -10., 10.);
    TH1D* h1_fitPar0 = new TH1D( Form("fitPar0_%s", iR->getName().c_str()), "", 300, 0., 20.);
    TH1D* h1_fitPar1 = new TH1D( Form("fitPar1_%s", iR->getName().c_str()), "", 300, -0.14, 0.01 );
    TH1D* h1_fitInt  = new TH1D( Form("fitInt_%s", iR->getName().c_str()), "", 300, 0., 200.);
    h1_fitPar0 ->SetXTitle( "Fit Parameter 0" );
    h1_fitPar1 ->SetXTitle( "Fit Parameter 1" );
    h1_fitInt  ->SetXTitle( "Fit Integral (M_{T2} > 150)" );


    float fitPar0_MC    = -999.;
    float fitPar0Err_MC = -999.;
    float fitPar1_MC    = -999.;
    float fitPar1Err_MC = -999.;
    float fitInt_MC     = -999.;


    for( int iter=0; iter<niter; ++iter ) {

     
      TH1D::AddDirectory(kTRUE);

      MT2EstimateQCD* thisQCD = new MT2EstimateQCD( *(qcd->get( *iR )) );
      MT2EstimateQCD* thisAll = new MT2EstimateQCD( *(all->get( *iR )) );


      if( lumiToCheck>0. && iter>0 ) {
        thisQCD->randomizePoisson(lumiToCheck/prescale, 13+iRegion+iter);
        //thisAll->randomizePoisson(lumiToCheck, 17+iRegion+iter);
      }

      //thisAll->finalize();
      thisQCD->finalize();

     
      TH1D* thisRatioAll = thisAll->ratio;    
      TH1D* thisRatioQCD = thisQCD->ratio;    
      //Th2_axesF1* thisFit = thisAll->exp;
      TF1* thisFitQCD = thisQCD->exp; 

      Double_t xMin = thisRatioAll->GetBinLowEdge(thisRatioAll->FindBin(150.));
      Double_t xMax = thisRatioAll->GetXaxis()->GetXmax();

      int binMin = thisRatioAll->FindBin(xMin);
      int binMax = thisRatioAll->FindBin(xMax);

      Double_t intErr_ratio;
      Double_t int_ratio = thisRatioQCD->IntegralAndError( binMin, binMax, intErr_ratio );
      Double_t int_f1    = thisFitQCD->Integral( xMin, xMax );
      Double_t intErr_f1 = 0.; //thisFit->IntegralError( xMin, xMax );


      float zTest = (int_ratio-int_f1) / ( sqrt( intErr_ratio*intErr_ratio + intErr_f1*intErr_f1 ) );


      if( iter>0 ) {

        h1_thisPull->Fill( zTest );
   
        h1_fitPar0->Fill( thisFitQCD->GetParameter(0) );
        h1_fitPar1->Fill( thisFitQCD->GetParameter(1) );
        h1_fitInt->Fill( int_f1 );

      } else {

        h1_pull->Fill( zTest );
        fitPar0_MC    = thisFitQCD->GetParameter(0);
        fitPar0Err_MC = thisFitQCD->GetParError(0);
        fitPar1_MC    = thisFitQCD->GetParameter(1);
        fitPar1Err_MC = thisFitQCD->GetParError(1);
        fitInt_MC     = int_f1;

      }

      
      if( lumiToCheck<=0. ) {

        TCanvas* c1 = new TCanvas( "c2", "", 600, 600 );
        c1->cd();
        c1->SetLogx();
        c1->SetLogy();


        float yMax    = thisRatioAll->GetMaximum()*5.;
        float yMinAll = thisRatioQCD->GetMinimum()/2.;
        float yMin    = thisRatioQCD->GetMinimum()/2.;
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
        legend->AddEntry( thisRatioAll, "MC (all)", "P" );
        legend->AddEntry( thisRatioQCD, "MC (QCD Only)", "P" );
        legend->AddEntry( thisFitQCD, "Exp. Fit", "L" );
        legend->Draw("same");

        TLine* lineLeft = new TLine( 60., yMin, 60., yMax );
        lineLeft->SetLineStyle(2);
        lineLeft->Draw("same");

        TLine* lineRight = new TLine( 100., yMin, 100., yMax );
        lineRight->SetLineStyle(2);
        lineRight->Draw("same");

        TH1D* fit_band = MT2DrawTools::getBand(thisFitQCD, Form("band_%s", thisFitQCD->GetName()) );
        fit_band->Draw("C E3 same");

        thisFitQCD->SetLineColor(46); 
        thisFitQCD->SetLineWidth(2); 
        thisFitQCD->Draw("L same");
     
        thisRatioAll->SetMarkerStyle(20);
        thisRatioAll->SetMarkerSize(1.3);
        thisRatioAll->Draw("P same");
     
        thisRatioQCD->SetMarkerStyle(24);
        thisRatioQCD->SetMarkerSize(1.3);
        thisRatioQCD->Draw("P same");
     
        gPad->RedrawAxis();

        TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
        labelTop->Draw("same");

        if( lumiToCheck>0. ) {
          c1->SaveAs( Form("%s/ratio_%s_iter%d.eps", toymcdir.c_str(), iR->getName().c_str(), iter) );
          c1->SaveAs( Form("%s/ratio_%s_iter%d.pdf", toymcdir.c_str(), iR->getName().c_str(), iter) );
        } else {
          c1->SaveAs( Form("%s/ratio_%s.eps", fitsdir.c_str(), iR->getName().c_str()) );
          c1->SaveAs( Form("%s/ratio_%s.pdf", fitsdir.c_str(), iR->getName().c_str()) );
        }

        delete c1;
        delete h2_axes;
        delete fit_band;

      } // draw only if lumi <=0 


      iRegion++;

    } // for iter

    if( lumiToCheck>0. ) {

//    TCanvas* c3 = new TCanvas("c3", "", 600, 600 );
//    c3->cd();
//    h1_thisPull->SetXTitle("(MC - fit) / (#sigma_{MC} #oplus #sigma_{fit})");
//    h1_thisPull->Draw();
//    TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation( lumiToCheck );
//    labelTop->Draw("same");
//    c3->SaveAs( Form("%s/pull_%s.eps", toymcdir.c_str(), iR->getName().c_str()) );
//    c3->SaveAs( Form("%s/pull_%s.pdf", toymcdir.c_str(), iR->getName().c_str()) );

      float par0_reso = drawSinglePlot( toymcdir, *iR, lumiToCheck, h1_fitPar0, fitPar0_MC, fitPar0Err_MC );
      float par1_reso = drawSinglePlot( toymcdir, *iR, lumiToCheck, h1_fitPar1, fitPar1_MC, fitPar1Err_MC );
      drawSinglePlot( toymcdir, *iR, lumiToCheck, h1_fitInt , fitInt_MC );

      bool golden = (par0_reso<0.2) && (par1_reso<0.2);
      bool good   = (par0_reso<0.5) && (par1_reso<0.5);

      
      if( good ) {
        ofs_ok       << iR->getName();
        if( golden ) ofs_ok << " ***";
        ofs_ok << std::endl;
      } else {
        ofs_problem  << iR->getName();
        //if( !par0_ok ) ofs_problem << " (Par0 notOK)";
        //if( !par1_ok ) ofs_problem << " (Par1 notOK)";
        ofs_problem << std::endl;
      }


    } // for iter


  } // for regions


  if( lumiToCheck>0. ) {
    ofs_problem.close();
    ofs_ok     .close();
    std::cout << "-> Wrote OK regions to: " << fileName_ok << std::endl;
    std::cout << "-> Wrote problematic regions to: " << fileName_problem << std::endl;
  }

  return 0;

}



float drawSinglePlot( const std::string& dir, const MT2Region& region, float lumi, TH1D* histo, float mcValue, float mcErr ) {

  TCanvas* c3 = new TCanvas("c3", "", 600, 600 );
  c3->cd();

  float xMin = histo->GetXaxis()->GetXmin();
  float xMax = histo->GetXaxis()->GetXmax();
  float yMax = histo->GetMaximum()*1.2;

  TH2D* h2_axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0., yMax );
  h2_axes->SetXTitle( histo->GetXaxis()->GetTitle() );
  h2_axes->SetYTitle( "Number of Toys" );
  h2_axes->Draw();


  TLine* line = new TLine( mcValue, 0., mcValue, yMax );
  line->SetLineColor(46);
  line->SetLineWidth(2);

  TLine* lineUp = new TLine( mcValue+mcErr, 0., mcValue+mcErr, yMax );
  lineUp->SetLineColor(46);
  lineUp->SetLineWidth(2);
  lineUp->SetLineStyle(2);

  TLine* lineDown = new TLine( mcValue-mcErr, 0., mcValue-mcErr, yMax );
  lineDown->SetLineColor(46);
  lineDown->SetLineWidth(2);
  lineDown->SetLineStyle(2);

  line->Draw("same");
  lineUp->Draw("same");
  lineDown->Draw("same");
  histo->SetFillColor(kGray);
  histo->Draw("same");

  std::vector<std::string> regionNiceNames = region.getNiceNames();

  float xMin_label = (xMin<0.) ? 0.2 : 0.6;
  float xMax_label = (xMin<0.) ? 0.5 : 0.9;

  TPaveText* regionName = new TPaveText( xMin_label, 0.8, xMax_label, 0.9, "brNDC" );
  regionName->SetTextAlign( 11 );
  regionName->SetTextSize( 0.035 );
  regionName->SetFillColor( 0 );
  regionName->AddText( regionNiceNames[0].c_str() );
  regionName->AddText( regionNiceNames[1].c_str() );
  regionName->Draw("same");



  float mean = histo->GetMean();
  float rms  = histo->GetRMS();
  
  TPaveText* meanText = new TPaveText( xMin_label, 0.6, xMax_label, 0.75, "brNDC" );
  meanText->SetTextAlign( 11 );
  meanText->SetFillColor(0);
  meanText->SetTextSize(0.035);
  meanText->AddText( Form("%.0lf Toys", histo->GetEntries()) );
  if( mean > 0.01 ) {
    meanText->AddText( Form("Mean = %.3f", mean ) );
    meanText->AddText( Form("RMS  = %.3f", rms ) );
  } else {
    meanText->AddText( Form("Mean = %.5f", mean ) );
    meanText->AddText( Form("RMS  = %.5f", rms ) );
  }
  meanText->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation( lumi );
  labelTop->Draw("same");

  gPad->RedrawAxis();

  c3->SaveAs( Form("%s/%s.eps", dir.c_str(), histo->GetName()) );
  c3->SaveAs( Form("%s/%s.pdf", dir.c_str(), histo->GetName()) );

  delete c3; 
  delete h2_axes; 
  delete line; 
  delete lineDown; 
  delete lineUp; 

  //bool returnBool = fabs(mean) > rms;

  return rms/fabs(mean);

}
