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







void projectFromInclusive( MT2Analysis<MT2Estimate>* analysis, MT2Analysis<MT2EstimateTree>* ana_inclusive, const std::string& selection );
void fillFromTreeAndRatio( MT2Estimate* estimate, MT2Estimate* r_effective, TTree* tree, TF1* f1_ratio );
void get_rHat( MT2Analysis<MT2Estimate>* rHat, MT2Analysis<MT2EstimateTree>* analysis );
void get_fJets( MT2Analysis<MT2Estimate>* fJets, MT2Analysis<MT2EstimateTree>* analysis );
void drawSingleFit( const std::string& outdir, TF1* thisFitQCD, MT2EstimateQCD* qcd, MT2EstimateQCD* all, float xMin_fit, float xMax_fit );
void drawClosure( const std::string& outputdir, MT2Analysis<MT2Estimate>* estimate, MT2Analysis<MT2Estimate>* mcTruth );



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

  TH1D::AddDirectory(kTRUE);



  // always start from inclusive qcd tree:
  MT2Analysis<MT2EstimateTree>* qcdTree_mc = MT2Analysis<MT2EstimateTree>::readFromFile( "EventYields_data_Run2015D_25nsGolden_v4/qcdControlRegion/mc.root", "qcdCRtree" );
  //MT2Analysis<MT2EstimateTree>* qcdTree_mc = MT2Analysis<MT2EstimateTree>::readFromFile( qcdCRdir + "/mc.root", "qcdCRtree" );


  std::string regionsSet = cfg.regionsSet();

  MT2Analysis<MT2Estimate>*  mcTruth = new MT2Analysis<MT2Estimate>( "mcTruth", regionsSet );
  //projectFromInclusive( mcTruth, qcdTree_mc, "id>=154 && id<=157 && deltaPhiMin>0.3 && mt2>200." );
  projectFromInclusive( mcTruth, qcdTree_mc, "id>=100 && id<200 && deltaPhiMin>0.3 && mt2>200." );
  //projectFromInclusive( mcTruth, qcdTree_mc, "((id>=154 && id<=157 && ht>=1500.) || ( id>=100 && id<200 && ht<1500.)) && deltaPhiMin>0.3 && mt2>200." );


  MT2Analysis<MT2Estimate>* estimate     = new MT2Analysis<MT2Estimate>("qcdEstimate", regionsSet);
  MT2Analysis<MT2Estimate>* r_effective  = new MT2Analysis<MT2Estimate>("r_effective", regionsSet);
  MT2Analysis<MT2Estimate>* r_hat        = new MT2Analysis<MT2Estimate>("r_hat"      , "zurich_onlyJets_noB");
  MT2Analysis<MT2Estimate>* f_jets       = new MT2Analysis<MT2Estimate>("f_jets"     , "zurich_onlyHT");


  MT2Analysis<MT2EstimateTree>* mcTruth_ht   = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mcTruth_ht"  , "zurich_onlyHT"      , qcdTree_mc, "id>=154 && id<=157 && deltaPhiMin>0.3 && mt2>200." );
  MT2Analysis<MT2EstimateTree>* mcTruth_jets = MT2EstimateTree::makeAnalysisFromInclusiveTree( "mcTruth_jets", "zurich_onlyJets_noB", qcdTree_mc, "id>=154 && id<=157 && deltaPhiMin>0.3 && mt2>200." );

  get_rHat ( r_hat , mcTruth_jets );
  get_fJets( f_jets, mcTruth_ht   );


  MT2Analysis<MT2EstimateQCD>* mc_all  = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "mc"     , cfg.qcdRegionsSet(), qcdTree_mc, "" );
  MT2Analysis<MT2EstimateQCD>* qcdOnly = MT2EstimateQCD::makeAnalysisFromInclusiveTree( "qcdOnly", cfg.qcdRegionsSet(), qcdTree_mc, "id>=100 && id<200" );



  std::string qcdCRdir = cfg.getEventYieldDir() + "/qcdControlRegion/";
  std::string outputdir = qcdCRdir;
  std::string fitsDir = qcdCRdir + "/fits";
  system( Form("mkdir -p %s", fitsDir.c_str() ));


  std::set<MT2Region> regions = estimate->getRegions();


  std::string lastRegionName = "";

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* this_estimate     = estimate    ->get( *iR );
    MT2Estimate* this_r_effective  = r_effective ->get( *iR );

    MT2Region* regionToMatch;
    if( iR->nBJetsMin()==3 && iR->nJetsMin()==2 ) 
      regionToMatch = new MT2Region( iR->htMin(), iR->htMax(), 4, 6, iR->nBJetsMin(), iR->nBJetsMax() );
    else
      regionToMatch = new MT2Region( *iR );

    MT2Estimate* this_r_hat = r_hat->getWithMatch( *regionToMatch );
    MT2Estimate* this_f_jets = f_jets->getWithMatch( *regionToMatch );


    MT2EstimateQCD* matchedEstimate_all = mc_all->getWithMatch( *iR );
    MT2EstimateQCD* matchedEstimate_qcd = qcdOnly->getWithMatch( *iR );


    float xMin_fit = (iR->htMin()>=1000.) ? 70. : 60.;
    float xMax_fit = 100.;
    TF1* f1_ratio = matchedEstimate_qcd->getFit( "pow", xMin_fit, xMax_fit );

    if( matchedEstimate_qcd->region->getName()!=lastRegionName ) // draw only one per HT region
      drawSingleFit( fitsDir, f1_ratio, matchedEstimate_qcd, matchedEstimate_all, xMin_fit, xMax_fit );

    lastRegionName = matchedEstimate_qcd->region->getName();

    fillFromTreeAndRatio( this_estimate, this_r_effective, matchedEstimate_qcd->tree, f1_ratio );

    int bin_bJets = this_r_hat ->yield->FindBin(iR->nBJetsMin());
    float thisRhatValue = this_r_hat ->yield->GetBinContent( bin_bJets );
    this_estimate->yield->Scale( thisRhatValue );

    int bin_jets = this_f_jets->yield->FindBin(iR->nJetsMin() );
    float thisFjetsValue = this_f_jets ->yield->GetBinContent( bin_jets );
    if( iR->nBJetsMin()==3 && iR->nJetsMin()==2 ) {
      int bin_jets_2 = this_f_jets->yield->FindBin(4);
      thisFjetsValue += this_f_jets ->yield->GetBinContent( bin_jets_2 );
    }
    this_estimate->yield->Scale( this_f_jets->yield->GetBinContent( bin_jets  ) );


  }  // for regions
      
    
  estimate ->writeToFile( outputdir + "/qcdEstimate.root", "recreate" );
  mcTruth->writeToFile( outputdir + "/qcdEstimate.root" );
  r_effective ->writeToFile( outputdir + "/qcdEstimate.root" );
  r_hat ->writeToFile( outputdir + "/qcdEstimate.root" );
  f_jets ->writeToFile( outputdir + "/qcdEstimate.root" );


  mcTruth->setColor(kQCD);
  estimate->setColor(kBlack);
  drawClosure( qcdCRdir, estimate, mcTruth );

  mc_all ->writeToFile( outputdir + "/mcFits.root", "recreate" );
  qcdOnly->writeToFile( outputdir + "/mcFits.root" );



  return 0;

}





void projectFromInclusive( MT2Analysis<MT2Estimate>* analysis, MT2Analysis<MT2EstimateTree>* ana_inclusive, const std::string& selection ) {

  MT2EstimateTree* treeInclusive = ana_inclusive->get( MT2Region("HT450toInf_j2toInf_b0toInf") );
  if( treeInclusive==0 ) {
    std::cout << "ERROR!! You need to pass an inclusive MT2EstimateTree Analysis to use this function!" << std::endl;
    exit(19191);
  }

  std::set<MT2Region> regions = analysis->getRegions();


  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* thisEstimate= analysis->get( *iR );

    std::string fullSelection = iR->getRegionCuts();
    if( selection!="" ) fullSelection = fullSelection + " && " + selection;
    treeInclusive->tree->Project( thisEstimate->yield->GetName(), "mt2", Form("weight*(%s)", fullSelection.c_str()) );

  } // for regions


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

  // this tree has only lowDeltaPhi events
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

  c1->SaveAs( Form("%s/ratio_%s.pdf", outdir.c_str(), qcd->region->getName().c_str()) );
  c1->SaveAs( Form("%s/ratio_%s.eps", outdir.c_str(), qcd->region->getName().c_str()) );
  

  delete c1;
  delete h2_axes;
  delete h_band;

}






void get_rHat( MT2Analysis<MT2Estimate>* rHat, MT2Analysis<MT2EstimateTree>* analysis ) {

  std::set<MT2Region> regions = rHat->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* thisEst = rHat->get(*iR);

    std::string name(thisEst->yield->GetName());
    int nBins = 4;
    Double_t bins[nBins+1];
    bins[0] = 0.;
    bins[1] = 1.;
    bins[2] = 2.;
    bins[3] = 3.;
    bins[4] = 10.;

    delete thisEst->yield;
    thisEst->yield = new TH1D( name.c_str(), "", nBins, bins );
    thisEst->yield->Sumw2();

    MT2EstimateTree* thisTree = analysis->get(*iR);
    thisTree->tree->Project( name.c_str(), "nBJets", "weight" );

    thisEst->yield->Scale( 1./thisEst->yield->Integral(1, nBins+1) );

  } // for regions

    
}


void get_fJets( MT2Analysis<MT2Estimate>* fJets, MT2Analysis<MT2EstimateTree>* analysis ) {

  std::set<MT2Region> regions = fJets->getRegions();

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* thisEst = fJets->get(*iR);

    std::string name(thisEst->yield->GetName());
    int nBins = 3;
    Double_t bins[nBins+1];
    bins[0] = 2.;
    bins[1] = 4.;
    bins[2] = 7.;
    bins[3] = 25.;

    delete thisEst->yield;
    thisEst->yield = new TH1D( name.c_str(), "", nBins, bins );
    thisEst->yield->Sumw2();

    MT2EstimateTree* thisTree = analysis->get(*iR);
    thisTree->tree->Project( name.c_str(), "nJets", "weight" );

    thisEst->yield->Scale( 1./thisEst->yield->Integral(1, nBins+1) );

  } // for regions

}



void drawClosure( const std::string& outputdir, MT2Analysis<MT2Estimate>* estimate, MT2Analysis<MT2Estimate>* mcTruth ) {


  system(Form("mkdir -p %s", outputdir.c_str()));

  
  std::set<MT2Region> MT2Regions = estimate->getRegions();
  
  TH1D* hdata = new TH1D("hdata", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  hdata->Sumw2();
  hdata->GetYaxis()->SetTitle("Entries");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( kBlack );
  hdata->SetMarkerColor( kBlack );
  
  TH1D* hestimate = new TH1D("hestimate", "", (int) MT2Regions.size(), 0, (int) MT2Regions.size());
  hestimate->Sumw2();
  hestimate->GetYaxis()->SetTitle("Entries");
  hestimate->SetFillColor(0);
  hestimate->SetLineColor( mcTruth->getColor() );
  hestimate->SetMarkerColor( mcTruth->getColor() );
  hestimate->SetMarkerStyle(20);
  hestimate->SetMarkerSize(1.6);
  
  TH1D* hPull = new TH1D("hPull", "", 20, -5, 5);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Data Driven - MC)/#sigma");
  hPull->GetYaxis()->SetTitle("Entries");
  
  
  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::string fullPath = outputdir;

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      TH1D* h_first = estimate->get(*iMT2)->yield;
      h_first->SetMarkerStyle(20);
      h_first->SetMarkerSize(1.6);
      h_first->SetLineColor( estimate->getColor() );
      h_first->SetMarkerColor( estimate->getColor() );


      int nBins = h_first->GetXaxis()->GetNbins();
      double err_data;
      double int_data = h_first->IntegralAndError(1, nBins+1, err_data);
 
      hdata->SetBinContent(iRegion, int_data);
      hdata->SetBinError(iRegion, err_data);
      

      hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
      
      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.15);
      pad1->Draw();
      pad1->cd();

      TH1D* h_second = mcTruth->get(*iMT2)->yield;

      h_second->SetLineColor( mcTruth->getColor() );
      h_second->SetMarkerColor( mcTruth->getColor() );
      h_second->SetMarkerStyle(20);
      h_second->SetMarkerSize(1.6);
      
   
      double err_int;
      double integral = h_second->IntegralAndError(1, nBins+1, err_int);
      hestimate->SetBinContent(iRegion, h_second->Integral());
      hestimate->SetBinError(iRegion, err_int);
      hestimate->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );

      if(int_data>0 && integral>0)
      hPull->Fill((int_data-integral)/sqrt(err_data*err_data+err_int*err_int));


      float xMin = h_first->GetXaxis()->GetXmin();
      float xMax = h_first->GetXaxis()->GetXmax();
      float yMax_1 = h_first->GetMaximum()*1.5;
      float yMax_2 = 1.2*(h_first->GetMaximum() + h_first->GetBinError(h_first->GetMaximumBin()));
      float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
      float yMax_3 = h_second->GetMaximum()*1.5;
      float yMax_4 = 1.2*(h_second->GetMaximum() + h_second->GetBinError(h_second->GetMaximumBin()));
      float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
      float yMax = (yMax1>yMax2) ? yMax1 : yMax2;

      TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
      h2_axes->SetXTitle("M_{T2} [GeV]");
      h2_axes->SetYTitle("Entries");

      h2_axes->Draw();
 
      //std::vector<std::string> niceNames = iMT2->getNiceNames();
      for( unsigned i=0; i<niceNames.size(); ++i ) {

        float yMax = 0.9-(float)i*0.05;
        float yMin = yMax - 0.05;
        TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
        regionText->SetTextSize(0.035);
        regionText->SetTextFont(42);
        regionText->SetFillColor(0);
        regionText->SetTextAlign(11);
        regionText->AddText( niceNames[i].c_str() );
        regionText->Draw("same");

      }


      TLegend* legend = new TLegend( 0.6, 0.9-2.*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
      legend->AddEntry( h_first, "data-driven", "P" );
      legend->AddEntry( h_second, "MC QCD", "P" );

      legend->Draw("same");

      h_first->Draw("P same");
      h_second->Draw("P same");
      //      bgStack.Draw("histoE, same");

      TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
      labelTop->Draw("same");

      gPad->RedrawAxis();

      c1->cd();
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.1);
      pad2->Draw();
      pad2->cd();

      std::string thisName = Form("%s_ratio", h_first->GetName());
      TH1D* h_ratio = (TH1D*) h_first->Clone(thisName.c_str());
      h_ratio->Divide(h_second);
      h_ratio->SetStats(0);	    
      h_ratio->SetMarkerStyle(20);
      h_ratio->SetLineColor(1);
      //      h_ratio->SetMarkerSize(0.02);
      h_ratio->GetXaxis()->SetLabelSize(0.00);
      h_ratio->GetXaxis()->SetTickLength(0.09);
      h_ratio->GetYaxis()->SetNdivisions(5,5,0);
      h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
      h_ratio->GetYaxis()->SetTitleSize(0.17);
      h_ratio->GetYaxis()->SetTitleOffset(0.4);
      h_ratio->GetYaxis()->SetLabelSize(0.17);
      h_ratio->GetYaxis()->SetTitle("Ratio");
      
      

      h_ratio->SetLineWidth(2);
      //h_ratio->Draw("PE");
      
      TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, 0.0, 2.0 );
      
      TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
      lineCentral->SetLineColor(1);
      

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      h_ratio->Draw("pe,same");

      gPad->RedrawAxis();
      
      c1->cd();

      c1->SaveAs( Form("%s/closure_%s.eps", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/closure_%s.pdf", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;
      delete h2_axes;
      delete h2_axes_ratio;
      delete h_ratio;
      delete h_first;
      delete h_second;
      
      ++iRegion;

  } // for MT2 regions

  
  TCanvas* c2 = new TCanvas("c2", "", 1200, 600);
  c2->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  pad1->SetLogy();
  
  float yMax_1 = hdata->GetMaximum()*1.5;
  float yMax_2 = 1.2*(hdata->GetMaximum() + hdata->GetBinError(hestimate->GetMaximumBin()));
  float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
  float yMax_3 = hestimate->GetMaximum()*1.5;
  float yMax_4 = 1.2*(hestimate->GetMaximum() + hestimate->GetBinError(hestimate->GetMaximumBin()));
  float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
  float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
  
  float yMin = 1e-1;
  yMax*=20.;
  
  hestimate->GetXaxis()->SetRangeUser(0, (int) MT2Regions.size());
  hestimate->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate->GetXaxis()->LabelsOption("v");
  hestimate->Draw("PE");


  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( estimate->getColor() );
  hdata->SetMarkerColor( estimate->getColor() );

  hdata->Draw("pe,same");

  TLegend* legend = new TLegend( 0.8, 0.9-0.06-0.06-0.06, 0.93, 0.9-0.06 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( hdata, "data-driven", "PL" );
  legend->AddEntry( hestimate, "QCD MC", "PL" );

  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
  labelTop->Draw("same");
  
  TLine* lHT[3];
  for( int iHT=1; iHT < 4; iHT++ ){
    lHT[iHT-1] = new TLine(11*iHT, 0.0, 11*iHT, yMax );
    lHT[iHT-1]->SetLineColor(kBlack);
    lHT[iHT-1]->SetLineStyle(3);
    lHT[iHT-1]->SetLineWidth(2);

    lHT[iHT-1]->Draw("same");
  }

  int nHTRegions = 4;
  std::vector< std::string > htRegions;
  htRegions.push_back("low H_{T}");
  htRegions.push_back("medium H_{T}");
  htRegions.push_back("high H_{T}");
  htRegions.push_back("extreme H_{T}");
  
  TPaveText* htBox[nHTRegions];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){
    
    htBox[iHT] = new TPaveText(0.16+0.2*iHT, 0.9-0.06, 0.34+0.2*iHT, 0.9, "brNDC");
    htBox[iHT]->AddText( htRegions[iHT].c_str() );
    
    htBox[iHT]->SetBorderSize(0);
    htBox[iHT]->SetFillColor(kWhite);
    htBox[iHT]->SetTextSize(0.038);
    htBox[iHT]->SetTextAlign(21); // align centered
    htBox[iHT]->SetTextFont(62);
    htBox[iHT]->Draw("same");

  }

  gPad->RedrawAxis();
  
  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();

  std::string thisName = Form("%s_ratio", hdata->GetName());
  TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
  h_Ratio->Divide(hestimate);
  h_Ratio->SetStats(0);
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->GetXaxis()->SetLabelSize(0.00);
  h_Ratio->GetXaxis()->SetTickLength(0.09);
  h_Ratio->GetYaxis()->SetNdivisions(5,5,0);
  h_Ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  h_Ratio->GetYaxis()->SetTitleSize(0.17);
  h_Ratio->GetYaxis()->SetTitleOffset(0.4);
  h_Ratio->GetYaxis()->SetLabelSize(0.17);
  h_Ratio->GetYaxis()->SetTitle("Ratio");

  h_Ratio->SetLineWidth(2);
  //h_Ratio->Draw("P,E");
  

  TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( 0, MT2Regions.size(), 0.0, 2.0 );

  TLine* LineCentral = new TLine(0, 1.0, MT2Regions.size(), 1.0);
  LineCentral->SetLineColor(1);


  h2_axes_ratio->Draw("");
  LineCentral->Draw("same");
  h_Ratio->Draw("pe,same");

  for( int iHT=1; iHT < 4; iHT++ ){
    lHT[iHT-1]->Draw("same");
  }


  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/closure_allRegions.pdf", outputdir.c_str()) );
  c2->SaveAs( Form("%s/closure_allRegions.eps", outputdir.c_str()) );


  gStyle->SetOptStat(1110);
  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  hPull->SetStats(1110);
  hPull->Draw("hist");
  c3->SaveAs( Form("%s/closure_pull.pdf", outputdir.c_str()) );
  c3->SaveAs( Form("%s/closure_pull.eps", outputdir.c_str()) );

  delete c2;
  delete c3;
  delete h2_axes_ratio;
  
}

