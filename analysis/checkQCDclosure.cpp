#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TFitter.h"
#include "TMath.h"

#include "../interface/MT2Config.h"
#include "../interface/MT2Analysis.h"
#include "../interface/MT2Region.h"
#include "../interface/MT2EstimateQCD.h"
#include "../interface/MT2DrawTools.h"




bool doLumiScan=false;
bool use_powerlaw=true;




Double_t powerlaw(Double_t *x,Double_t *par)
{
  //Double_t fitval = TMath::Exp( par[0] + par[1]*TMath::Log( x[0] ) );
  Double_t fitval = par[0]*TMath::Power( x[0], par[1] );
  return fitval;
}



//void randomizePoisson( MT2EstimateQCD* thisEstimate, float lumi, TRandom3 rand );
void drawSinglePlot( const std::string& dir, const MT2Region& region, float lumi, TH1D* histo, float mcValue, float mcErr=0. );
float chiSquareProb( TF1* fit, TH1D* histo, float xMin, float xMax );
MT2EstimateQCD* randomizePoisson( const std::string& name, MT2EstimateQCD* est, float scale, TRandom3* rand );


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


  float lumiScenario = 0.; //fb
  if( argc>2 ) {
    lumiScenario = std::atof(argv[2]); 
  } 

  if( lumiScenario>0. && !doLumiScan ) {
    std::cout << Form("-> Will run toyMC for lumi = %.0f fb-1", lumiScenario) << std::endl;
  } 


  if( doLumiScan ) {

    std::cout << std::endl << std::endl;
    std::cout << "*********** Will perform a lumi scan! This might take a while!!!!" << std::endl;
    std::cout << std::endl << std::endl;
    lumiScenario=3.; // starting point of the scan

  }


  int niter = (lumiScenario>0.) ? 301 : 1;
  if( doLumiScan ) niter = 201;


  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion";
  MT2Analysis<MT2EstimateQCD>* qcd = MT2Analysis<MT2EstimateQCD>::readFromFile( qcdCRdir + "/mcFits.root", "qcdOnly" );
  MT2Analysis<MT2EstimateQCD>* all = MT2Analysis<MT2EstimateQCD>::readFromFile( qcdCRdir + "/mcFits.root", "mc" );

  std::string fitsdir = qcdCRdir + "/fits";
  if( use_powerlaw ) fitsdir += "_pow";
  system( Form("mkdir -p %s", fitsdir.c_str()) );

  std::string toymcdir(Form("%s/toyMC_lumi%.0f_iter%d", qcdCRdir.c_str(), lumiScenario*1000., niter-1));
  if( use_powerlaw ) toymcdir += "_pow";
  system( Form("mkdir -p %s", toymcdir.c_str()) );


  std::set<MT2Region> regions = qcd->getRegions();

  std::string fileName_problem(Form("%s/qcdProblematic.txt", toymcdir.c_str()) );
  std::string fileName_ok     (Form("%s/qcdOK.txt", toymcdir.c_str()) );
  std::string fileName_fitProblem(Form("%s/fitProblem.txt", toymcdir.c_str()) );
  std::string fileName_lumiScan(Form("%s/lumiScan.txt", toymcdir.c_str()) );
  std::ofstream ofs_problem(fileName_problem);
  std::ofstream ofs_ok     (fileName_ok);
  std::ofstream ofs_fitProblem(fileName_fitProblem);
  std::ofstream ofs_lumiScan(fileName_lumiScan);

  TH1D* h1_pull = new TH1D("pull", "", 100, -10., 10.);

  TH1D* h1_chiSquare_fit  = new TH1D("chiSquare_fit" , "", 25, 0., 1.0001 );
  TH1D* h1_chiSquare_tail = new TH1D("chiSquare_tail", "", 25, 0., 1.0001 );



  int iRegion = 0;
  TRandom3* rand = new TRandom3(13);

  int nRegions = regions.size();
  int binMax_regions = 0;
  TH1D* h1_yield_mc = new TH1D("yield_mc", "", nRegions, 0., nRegions );
  TH1D* h1_yield_est = new TH1D("yield_est", "", nRegions, 0., nRegions );


  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {


    if( iR->nBJetsMin()>=2 ) continue; // who cares
    //if( iR->htMin()<1500. ) continue; // who cares
    //if( iR->getName()=="HT450to575_j7toInf_b0" ) continue; // this fit fails!

    if( doLumiScan )
      std::cout << "region: " << iR->getName() << std::endl;

    float prescale = 1.;
    if( iR->htMin() < 500. ) prescale = 180.;  // HLT_HT350 for HT>=450
    else if( iR->htMin() < 600. ) prescale = 60.; // HLT_HT_475 for HT >=575


    bool stopLumiScan = false;
    float lumiStep = 0.5;

    MT2EstimateQCD* thisQCD_old = new MT2EstimateQCD( *(qcd->get( *iR )) );
    MT2EstimateQCD* thisAll = new MT2EstimateQCD( *(all->get( *iR )) );



    for( float lumiToCheck=lumiScenario; (lumiToCheck<=100. || !doLumiScan) && !stopLumiScan; lumiToCheck+=lumiStep ) {

      if( lumiToCheck>=10. ) lumiStep = 1.;
      if( lumiToCheck>=20. ) lumiStep = 5.;
      if( lumiToCheck>=30. ) lumiStep = 10.;
      if( lumiToCheck>=60. ) lumiStep = 20.;

      if( doLumiScan )
        std::cout << "  lumi: " << lumiToCheck << std::endl;

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


      TH1D::AddDirectory(kTRUE);


      for( int iter=0; iter<niter; ++iter ) {


        MT2EstimateQCD* thisQCD = new MT2EstimateQCD(*thisQCD_old);
        if( lumiToCheck>0. && iter>0 ) {
          //int seed = 13+iRegion+(int)100000.*iter+(int)1000.*lumiToCheck;
          std::string newname(Form("%s_%d_%.0f", iR->getName().c_str(), iter, 1000.*lumiToCheck));
          thisQCD = randomizePoisson( newname, thisQCD_old, lumiToCheck/prescale, rand );
          //thisQCD->randomizePoisson(lumiToCheck/prescale, seed );
          //thisAll->randomizePoisson(lumiToCheck, 17+iRegion+iter);

        }

        //thisAll->finalize();
        thisQCD->finalize();

       
        TH1D* thisRatioAll = thisAll->ratio;    
        TH1D* thisRatioQCD = thisQCD->ratio;    

        Double_t xMin = thisRatioAll->GetXaxis()->GetXmin();
        Double_t xMax = thisRatioAll->GetXaxis()->GetXmax();

        Double_t xMin_fit = 70;
        

        TF1* thisFitQCD;
        if( use_powerlaw ) {
          TF1* thisFit = new TF1( "fit", powerlaw, xMin_fit, 600., 2 );
          thisRatioQCD->Fit( thisFit, "QR0" );
          thisFitQCD = new TF1( "fit_draw", powerlaw, xMin, xMax, 2 );
          thisFitQCD->SetParameter( 0, thisFit->GetParameter(0) );
          thisFitQCD->SetParameter( 1, thisFit->GetParameter(1) );
        } else {
          thisFitQCD = thisQCD->exp; 
        }


        h1_yield_mc ->GetXaxis()->SetBinLabel(iRegion+1, iR->getName().c_str());
        h1_yield_est->GetXaxis()->SetBinLabel(iRegion+1, iR->getName().c_str());

        Double_t xMin_int = thisRatioAll->GetBinLowEdge(thisRatioAll->FindBin(150.));
        int binMin = thisRatioAll->FindBin(xMin_int);
        int binMax = thisRatioAll->FindBin(xMax);

        Double_t intErr_mc;
        Double_t integral_mc = thisQCD->hDphi->IntegralAndError( binMin, binMax, intErr_mc );
        h1_yield_mc->SetBinContent(iRegion+1, integral_mc);
        h1_yield_mc->SetBinError(iRegion+1, intErr_mc );

        Double_t integral_est = 0.;
        for( int iBin=binMin; iBin<=binMax; ++iBin ) {
          float x = thisQCD->hDphi->GetBinCenter(iBin);
          //float x = thisQCD->hDphi->GetBinLowEdge(iBin);
          integral_est += thisFitQCD->Eval(x)*thisQCD->lDphi->GetBinContent(iBin);
        }
        h1_yield_est->SetBinContent(iRegion+1, integral_est);
        h1_yield_est->SetBinError(iRegion+1, sqrt(integral_est) );
        binMax_regions = iRegion+1;


        Double_t intErr_ratio;
        Double_t int_ratio = thisRatioQCD->IntegralAndError( binMin, binMax, intErr_ratio );
        Double_t int_f1    = thisFitQCD->Integral( xMin_int, xMax );
        Double_t intErr_f1 = 0.; //thisFit->IntegralError( xMin, xMax );


        float zTest = (int_ratio-int_f1) / ( sqrt( intErr_ratio*intErr_ratio + intErr_f1*intErr_f1 ) );


        if( iter>0 ) {

          h1_thisPull->Fill( zTest );
     
          h1_fitPar0->Fill( thisFitQCD->GetParameter(0) );
          h1_fitPar1->Fill( thisFitQCD->GetParameter(1) );
          h1_fitInt->Fill( int_f1 );

        } else {

          // compute chisquare:
          float chisquare_fit  = chiSquareProb( thisFitQCD, thisRatioQCD, xMin_fit, 100. );
          float chisquare_tail = chiSquareProb( thisFitQCD, thisRatioQCD, 100., xMax );
          h1_chiSquare_fit ->Fill( chisquare_fit  );
          h1_chiSquare_tail->Fill( chisquare_tail );

          if( chisquare_tail<0.05 ) ofs_fitProblem << iR->getName() << std::endl;

          h1_pull->Fill( zTest );
          fitPar0_MC    = thisFitQCD->GetParameter(0);
          fitPar0Err_MC = thisFitQCD->GetParError(0);
          fitPar1_MC    = thisFitQCD->GetParameter(1);
          fitPar1Err_MC = thisFitQCD->GetParError(1);
          fitInt_MC     = int_f1;

        }

        
        if( lumiToCheck<=0. && !doLumiScan ) {

          TCanvas* c1 = new TCanvas( "c2", "", 600, 600 );
          c1->cd();
          c1->SetLogx();
          c1->SetLogy();


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
          delete h_band;

        } // draw only if lumi <=0 


        iRegion++;
        if( thisQCD!=0 ) delete thisQCD;
        if( thisFitQCD!=0 ) delete thisFitQCD;

      } // for iter



//      if( lumiToCheck>0. ) {


        if( !doLumiScan && lumiToCheck>0. ) {
          drawSinglePlot( toymcdir, *iR, lumiToCheck, h1_fitPar0, fitPar0_MC, fitPar0Err_MC );
          drawSinglePlot( toymcdir, *iR, lumiToCheck, h1_fitPar1, fitPar1_MC, fitPar1Err_MC );
          drawSinglePlot( toymcdir, *iR, lumiToCheck, h1_fitInt , fitInt_MC );
        }

        float par0_reso = h1_fitPar0->GetRMS()/fabs(h1_fitPar0->GetMean());
        float par1_reso = h1_fitPar1->GetRMS()/fabs(h1_fitPar1->GetMean());

        bool golden = (par0_reso<0.2) && (par1_reso<0.2);
        bool good   = (par0_reso<0.5) && (par1_reso<0.5);
        //std::cout << "par0_reso: " << par0_reso << " par1_reso: " << par1_reso << std::endl;

        
        if( doLumiScan ) {

          if( good ) {
            ofs_lumiScan << iR->getName() << " " << lumiToCheck << std::endl;
            stopLumiScan = true;
          }

        } else {

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

          break; // dont want to increase lumiToCheck

        } // if/else doLumiScan

 //     } // if lumiToCheck>0

      delete  h1_thisPull;
      delete  h1_fitPar0;
      delete  h1_fitPar1;
      delete  h1_fitInt;

      //if( thisFitQCD!=0 ) delete thisFitQCD;

    } // lumiScan loop

    delete thisQCD_old;
    delete thisAll;

  } // for regions


  if( lumiScenario>0. ) {

    if( doLumiScan ) {
      ofs_lumiScan.close();
      std::cout << "-> Wrote lumi scan to: " << fileName_lumiScan << std::endl;
    } else {
      ofs_problem.close();
      ofs_ok     .close();
      std::cout << "-> Wrote OK regions to: " << fileName_ok << std::endl;
      std::cout << "-> Wrote problematic regions to: " << fileName_problem << std::endl;
    }

  } else {

    ofs_fitProblem.close();
    std::cout << "-> Wrote fit problems to: " << fileName_fitProblem << std::endl;

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();

    h1_chiSquare_fit->SetFillColor(46);
    h1_chiSquare_fit->SetFillStyle(3004);
    h1_chiSquare_fit->SetLineColor(46);
    h1_chiSquare_fit->SetLineWidth(2);
    h1_chiSquare_fit->SetXTitle("#chi^{2} Probability");
    h1_chiSquare_fit->SetYTitle("Number of Regions" );

    h1_chiSquare_tail->SetFillColor(38);
    h1_chiSquare_tail->SetFillStyle(3004);
    h1_chiSquare_tail->SetLineColor(38);
    h1_chiSquare_tail->SetLineWidth(2);
    h1_chiSquare_tail->SetXTitle("#chi^{2} Probability");
    h1_chiSquare_tail->SetYTitle("Number of Regions" );

    h1_chiSquare_tail->Draw();
    h1_chiSquare_fit->Draw("same");

    TLegend* legend = new TLegend( 0.35, 0.7, 0.9, 0.9 );
    legend->SetTextSize(0.035);
    legend->SetFillColor(0);
    legend->AddEntry( h1_chiSquare_fit , Form("Fit Range (average=%.2f)", h1_chiSquare_fit->GetMean()), "F" );
    legend->AddEntry( h1_chiSquare_tail, Form("Tail (average=%.2f)", h1_chiSquare_tail->GetMean()), "F" );
    legend->Draw("same");

    TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
    labelTop->Draw("same");

    gPad->RedrawAxis();

    c1->SaveAs( Form("%s/chiSquareProb.eps", qcdCRdir.c_str()) );
    c1->SaveAs( Form("%s/chiSquareProb.pdf", qcdCRdir.c_str()) );



    c1 = new TCanvas( "c1", "", 1200, 600 );
    c1->SetLogy();

    TH2D* h2_axes = new TH2D("axes", "", binMax_regions, 0., binMax_regions, 10, 0.2, h1_yield_mc->GetMaximum()*4.);
    h2_axes->SetYTitle("Events");
    for( unsigned iBinx=1; iBinx<binMax_regions+1; ++iBinx ) 
      h2_axes->GetXaxis()->SetBinLabel(iBinx, h1_yield_mc->GetXaxis()->GetBinLabel(iBinx));
    h2_axes->Draw();

    h1_yield_mc->SetFillColor(kQCD);
    h1_yield_mc->SetLineColor(kBlack);

    h1_yield_est->SetMarkerStyle(20);
    h1_yield_est->SetMarkerSize(1.3);
    h1_yield_est->SetMarkerColor(kBlack);

    labelTop->Draw("same");

    TLegend* legend1 = new TLegend( 0.6, 0.8, 0.9, 0.9 );
    legend1->SetFillColor(0);
    legend1->SetTextSize(0.04);
    legend1->AddEntry( h1_yield_mc, "MC", "F" );
    legend1->Draw("same");

    TLegend* legend2 = new TLegend( 0.75, 0.8, 0.9, 0.9 );
    legend2->SetFillColor(0);
    legend2->SetTextSize(0.04);
    legend2->AddEntry( h1_yield_est, "Estimate", "PL" );
    legend2->Draw("same");

    h1_yield_mc->Draw("histo same");
    h1_yield_est->Draw("p same");

    gPad->RedrawAxis();

    c1->SaveAs(Form("%s/mc_closure.eps", qcdCRdir.c_str()) );
    c1->SaveAs(Form("%s/mc_closure.pdf", qcdCRdir.c_str()) );


    TFile* outfile = TFile::Open("prova.root", "recreate");
    outfile->cd();
    h1_pull->Write();
    h1_chiSquare_fit->Write();
    h1_chiSquare_tail->Write();
    h1_yield_est->Write();
    h1_yield_mc->Write();
    outfile->Close();
    
    delete c1;
    delete legend;

  }

  return 0;

}



void drawSinglePlot( const std::string& dir, const MT2Region& region, float lumi, TH1D* histo, float mcValue, float mcErr ) {

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

}



float chiSquareProb( TF1* fit, TH1D* histo, float xMin, float xMax ) {

  int binMin = histo->FindBin( xMin );
  int binMax = histo->FindBin( xMax );

  int nPoints = 0;
  int nPars = fit->GetNpar();
  float chiSquare = 0.;

  for( int iBin=binMin; iBin<=binMax; ++iBin ) {

    float x = histo->GetBinCenter(iBin); 
    float fValue = fit->Eval( x );
    float hValue = histo->GetBinContent(iBin);
    float hError = histo->GetBinError(iBin);

    if( hValue<=0. ) continue;

    float thisChi = (fValue-hValue)/hError;
    chiSquare += thisChi*thisChi;

    nPoints += 1;

  } // for points

  int NDF = nPoints - nPars;

  //std::cout << "chi2/NDF: " << chiSquare << "/" << NDF << std::endl;

  float chiSquareProb = TMath::Prob( chiSquare, NDF );
  //float chiSquareNorm = (NDF>0) ? chiSquare/((float)NDF) : 0.;

  return chiSquareProb;

}




MT2EstimateQCD* randomizePoisson( const std::string& name, MT2EstimateQCD* est, float scale, TRandom3* rand ) {

  MT2EstimateQCD* returnEst = new MT2EstimateQCD( name, *(est->region) );
  //MT2EstimateQCD* returnEst = new MT2EstimateQCD( *est );

  for( int ibin=1; ibin<est->lDphi->GetXaxis()->GetNbins()+1; ++ibin ) {

    int poisson_data = rand->Poisson(scale * est->lDphi->GetBinContent(ibin));
    returnEst->lDphi->SetBinContent(ibin, poisson_data);
    returnEst->lDphi->SetBinError( ibin, TMath::Sqrt(poisson_data) ); // here i want an approximation of the Poisson error
  }

  for( int ibin=1; ibin<est->hDphi->GetXaxis()->GetNbins()+1; ++ibin ) {

    int poisson_data = rand->Poisson(scale * est->hDphi->GetBinContent(ibin));
    returnEst->hDphi->SetBinContent(ibin, poisson_data);
    returnEst->hDphi->SetBinError( ibin, TMath::Sqrt(poisson_data) ); // here i want an approximation of the Poisson error
  }

  //returnEst->finalize();

  return returnEst;

}
