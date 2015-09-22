#include <iostream>



#include "interface/MT2Config.h"
#include "interface/MT2DrawTools.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2EstimateSyst.h"

#define mt2_cxx
#include "../interface/mt2.h"


#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMath.h"






void drawMt2VsB( MT2Config cfg, TTree* tree, const std::string& suffix, const std::string& legendTitle, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax, const std::string& additionalSel="");
float getP( MT2Config cfg, TTree* tree, const std::string& name, TTree* tree_data=0, const std::string& name_data="Data" );
TH1D* getRatioHisto( MT2Config cfg, TTree* tree, const std::string& name, const std::string& niceName, const std::string& selection = "ht>450. && mt2 > 200." );
MT2Analysis<MT2Estimate>* compute2bFromRatio( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, TF1* func );
void fillFromTreeRatio( TTree* tree, TH1D* yield_2b_extrapMC, TF1* func );
void compareHistos( MT2Config cfg, const std::string& saveName, TH1D* histo1, TH1D* histo2, TH1D* histo3=0, TH1D* histo4=0 );



int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|             Running computeZinv2B                  |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc!=2 ) {
    std::cout << "USAGE: ./computeZinv2B [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  MT2DrawTools::setStyle();


  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  
  std::string mainDir  = cfg.getEventYieldDir();
  std::string gammaDir = mainDir + "/gammaControlRegion";
  std::string zllDir = mainDir + "/zllControlRegion";

  std::string mcFile = mainDir + "/analyses.root";

  std::cout << "-> Running on : " << mcFile << std::endl;

  TH1D* histo_mc;
  TH1D* histo_mc_zll;
  TH1D* histo_data;
  TH1D* histo_mc_gjet;

  if( cfg.regionsSet()=="13TeV_inclusive") {

    TFile* file_Zinv = TFile::Open( mcFile.c_str() );
    TTree* tree_Zinv = (TTree*)file_Zinv->Get("ZJets/HT450toInf_j2toInf_b0toInf/tree_ZJets_HT450toInf_j2toInf_b0toInf");

    histo_mc = getRatioHisto( cfg, tree_Zinv, "zinv", "Z #rightarrow #nu#nu MC" );
    //histo_mc = getPHisto2( cfg, tree_Zinv, "zinv", "Z #rightarrow #nu#nu MC" );

    drawMt2VsB( cfg, tree_Zinv, "zinv", "Z #rightarrow #nu#nu", "mt2", "M_{T2} [GeV]", 100, 0., 1450. );
    //drawMt2VsB( cfg, tree_Zinv, "zinv", "Z #rightarrow #nu#nu", "nJets", "Number of Jets (p_{T}>30 GeV)", 8, 1.5, 9.5 );

    //TFile* file_gjet = TFile::Open( Form("%s/mc.root"  , gammaDir.c_str()) );
    //TFile* file_data = TFile::Open( Form("%s/data.root", gammaDir.c_str()) );

    //TTree* tree_gjet = (TTree*)file_gjet->Get("gammaCRtree/HT450toInf_j2toInf_b0toInf/tree_gammaCRtree_HT450toInf_j2toInf_b0toInf");
    //TTree* tree_data = (TTree*)file_data->Get("gammaCRtree/HT450toInf_j2toInf_b0toInf/tree_gammaCRtree_HT450toInf_j2toInf_b0toInf");

    //histo_mc_gjet = getPHisto2( cfg, tree_gjet, "gjetMC", "#gamma + Jets MC" );


    TFile* file_zll = TFile::Open( Form("%s/mc.root"  , zllDir.c_str()) );
    TTree* tree_zll = (TTree*)file_zll->Get("zllCR/HT450toInf_j2toInf_b0toInf/tree_zllCR_HT450toInf_j2toInf_b0toInf");
    TFile* file_zllData = TFile::Open( Form("%s/data.root"  , zllDir.c_str()) );
    TTree* tree_zllData = (TTree*)file_zllData->Get("data/HT450toInf_j2toInf_b0toInf/tree_data_HT450toInf_j2toInf_b0toInf");

    histo_mc_zll = getRatioHisto( cfg, tree_zll, "zllMC", "Z #rightarrow ll MC", "mt2>200. && ht>450. && Z_mass > 70. && Z_mass<110." );
    TH1D* histo_mc_zll_loose = getRatioHisto( cfg, tree_zll, "zllMC_loose", "Z #rightarrow ll MC (loose)", "ht>450. && Z_mass > 70. && Z_mass<110." );

    histo_data = getRatioHisto( cfg, tree_zllData, "zllData", "Z #rightarrow ll Data (loose)", "ht>450." );

    compareHistos( cfg, "compare_mc", histo_mc, histo_mc_zll, histo_mc_zll_loose );

    TFile* file = TFile::Open("prova.root", "RECREATE");
    file->cd();
    histo_mc->SetName("p");
    histo_mc->Write();
    file->Close();

  } else {

    TFile* file = TFile::Open("prova.root");
    histo_mc = (TH1D*)file->Get("p");
 
  }


  MT2Analysis<MT2EstimateTree>* zinv = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "ZJets");
  zinv->setFullName("Z + Jets");
  compute2bFromRatio( cfg, zinv, histo_mc->GetFunction("line") );


  return 0;

}




TH1D* getRatioHisto( MT2Config cfg, TTree* tree, const std::string& name, const std::string& niceName, const std::string& selection ) {

  float njetMin = 2;
  float njetMax = 12;
  int nbins_histo = njetMax-njetMin + 1;

  TH1D* h1_r_vs_nj = new TH1D(Form("r_vs_nj_%s", name.c_str()), "", nbins_histo, njetMin-0.5, njetMax+0.5 );
  

  std::string outdir = cfg.getEventYieldDir() + "/fits2b/" + name;
  system( Form("mkdir -p %s", outdir.c_str()) );

  //float yMin = 0.;
  //float yMax = 0.5;
  ////float yMax = 0.09;

  

  for( int njet=(int)njetMin; njet<=(int)njetMax; ++njet ) {

    std::string name_bjets(Form("nbjets_%d", njet));


    int nbjet_min = 0;
    int nbjet_max = 6;
    int nbjetbins = nbjet_max-nbjet_min;

    TH1D* h1_nbjets = new TH1D( name_bjets.c_str(), "", nbjetbins, nbjet_min-0.5, nbjet_max-0.5 );
    h1_nbjets->Sumw2();
    h1_nbjets->SetYTitle("Events");
    h1_nbjets->SetXTitle("Number of b-Jets");

    tree->Project( name_bjets.c_str(), "nBJets", Form("weight*(nJets==%d && %s )", njet, selection.c_str() ) );

    float n1 = h1_nbjets->GetBinContent( h1_nbjets->FindBin( 1 ) );
    float n1_err = h1_nbjets->GetBinError( h1_nbjets->FindBin( 1 ) );
    float n2 = h1_nbjets->GetBinContent( h1_nbjets->FindBin( 2 ) );
    float n2_err = h1_nbjets->GetBinError( h1_nbjets->FindBin( 2 ) );

    if( n1>0.0001 ) {

      float ratio = n2/n1;
      float ratio_err = sqrt( n2_err*n2_err/(n1*n1) + n1_err*n1_err*n2*n2/(n1*n1*n1*n1) );

      int bin = h1_r_vs_nj->FindBin( njet );
      h1_r_vs_nj->SetBinContent( bin, ratio );
      h1_r_vs_nj->SetBinError( bin, ratio_err );

    }

    delete h1_nbjets;

  }



  TF1* line = new TF1("line", "[0] + [1]*x", njetMin-0.5, njetMax+0.5 );
  line->SetLineColor(46);
  h1_r_vs_nj->Fit( line, "RQ+" );



  TCanvas* c1 = new TCanvas("c1_", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, njetMin-0.5, njetMax+0.5, 10, 0., 0.5 );
  h2_axes->SetXTitle("Number of Jets");
  h2_axes->SetYTitle("2b/1b Ratio");
  h2_axes->Draw();

  h1_r_vs_nj->SetMarkerStyle(20);
  h1_r_vs_nj->SetMarkerSize(1.6);


  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  h1_r_vs_nj->Draw("p same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/ratio.eps", outdir.c_str()) );
  c1->SaveAs( Form("%s/ratio.pdf", outdir.c_str()) );

  delete c1;
  delete h2_axes;

  h1_r_vs_nj->SetTitle( niceName.c_str() );

  return h1_r_vs_nj;

}

  

 

void drawMt2VsB( MT2Config cfg, TTree* tree, const std::string& suffix, const std::string& legendTitle, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax, const std::string& additionalSel) {


  std::string outputdir = cfg.getEventYieldDir() + "/gammaControlRegion/plots2b";
  system(Form("mkdir -p %s", outputdir.c_str()) );


  TCanvas* c1 = new TCanvas("c1_", "", 600, 600);

  TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600);
  c1_log->SetLogy();


  TH1D* h1_0b = new TH1D("0b", "", nBins, xMin, xMax); 
  TH1D* h1_1b = new TH1D("1b", "", nBins, xMin, xMax); 
  TH1D* h1_2b = new TH1D("2b", "", nBins, xMin, xMax); 

  h1_0b->Sumw2();
  h1_1b->Sumw2();
  h1_2b->Sumw2();

  if( additionalSel=="" ) {
    tree->Project("0b", varName.c_str(), "weight*(nBJets==0)");
    tree->Project("1b", varName.c_str(), "weight*(nBJets==1)");
    tree->Project("2b", varName.c_str(), "weight*(nBJets==2)");
  } else {
    tree->Project("0b", varName.c_str(), Form("weight*(nBJets==0 && %s)", additionalSel.c_str()) );
    tree->Project("1b", varName.c_str(), Form("weight*(nBJets==1 && %s)", additionalSel.c_str()) );
    tree->Project("2b", varName.c_str(), Form("weight*(nBJets==2 && %s)", additionalSel.c_str()) );
  }

  h1_0b->SetLineColor(kBlack);
  h1_1b->SetLineColor(kRed);
  h1_2b->SetLineColor(kBlue);

  h1_0b->SetLineWidth(2);
  h1_1b->SetLineWidth(2);
  h1_2b->SetLineWidth(2);

  h1_0b->SetLineStyle(1);
  h1_1b->SetLineStyle(2);
  h1_2b->SetLineStyle(3);


  float yMax_0b = h1_0b->GetMaximum()/h1_0b->Integral();
  float yMax_1b = h1_1b->GetMaximum()/h1_1b->Integral();
  float yMax_2b = h1_2b->GetMaximum()/h1_2b->Integral();

  float yMax = yMax_0b;
  if( yMax_1b>yMax ) yMax = yMax_1b;
  if( yMax_2b>yMax ) yMax = yMax_2b;
  yMax *= 1.1;

  TH2D* h2_axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0., yMax );
  h2_axes->SetXTitle(axisName.c_str());
  h2_axes->SetYTitle("Normalized to Unity");
  c1->cd();
  h2_axes->Draw();

  TH2D* h2_axes_log = new TH2D( "axes_log", "", 10, xMin, xMax, 10, 0.000001, yMax*2. );
  h2_axes_log->SetXTitle(axisName.c_str());
  h2_axes_log->SetYTitle("Normalized to Unity");
  c1_log->cd();
  h2_axes_log->Draw();


  TLegend* legend = new TLegend( 0.66, 0.66, 0.9, 0.9, legendTitle.c_str() );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_0b, "b = 0", "L" );
  legend->AddEntry( h1_1b, "b = 1", "L" );
  legend->AddEntry( h1_2b, "b = 2", "L" );
  c1->cd();
  legend->Draw("same");
  c1_log->cd();
  legend->Draw("same");

  c1->cd();
  h1_0b->DrawNormalized("histo same");
  h1_1b->DrawNormalized("histo same");
  h1_2b->DrawNormalized("histo same");

  c1_log->cd();
  h1_0b->DrawNormalized("histo same");
  h1_1b->DrawNormalized("histo same");
  h1_2b->DrawNormalized("histo same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop("CMS Simulation, #sqrt{s} = 13 TeV");

  c1->cd();
  labelTop->Draw("Same");
  gPad->RedrawAxis();

  c1_log->cd();
  labelTop->Draw("Same");
  gPad->RedrawAxis();

  std::string saveName = varName + "vsB_" + suffix;

  c1->SaveAs(Form("%s/%s.eps", outputdir.c_str(), saveName.c_str()));
  c1->SaveAs(Form("%s/%s.pdf", outputdir.c_str(), saveName.c_str()));
  c1->SaveAs(Form("%s/%s.png", outputdir.c_str(), saveName.c_str()));

  c1_log->SaveAs(Form("%s/%s_log.eps", outputdir.c_str(), saveName.c_str()));
  c1_log->SaveAs(Form("%s/%s_log.pdf", outputdir.c_str(), saveName.c_str()));
  c1_log->SaveAs(Form("%s/%s_log.png", outputdir.c_str(), saveName.c_str()));

  delete h2_axes;
  delete h2_axes_log;
  delete c1;
  delete c1_log;
  delete h1_0b;
  delete h1_1b;
  delete h1_2b;

}



MT2Analysis<MT2Estimate>* compute2bFromRatio( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, TF1* func ) {



  std::set<MT2Region> regions = mc->getRegions();

  MT2Analysis<MT2Estimate>* extrap = new MT2Analysis<MT2Estimate>( mc->getName() + "_extrap", regions );

  std::string outdir = cfg.getEventYieldDir() + "/fits2b";
  system(Form("mkdir -p %s", outdir.c_str()));


  if( cfg.regionsSet()=="13TeV_inclusive" ) { // simple closure test

    TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

    MT2EstimateTree* thisEst = mc->get( *(regions.begin()) );

    TTree* tree = thisEst->tree;

    TH1D* h1_mcTruth = new TH1D( "mcTruth", "", 50, 200., 1700. );
    TH1D* h1_mcFakes = new TH1D( "mcFakes", "", 50, 200., 1700. );
    TH1D* h1_closure = new TH1D( "closure", "", 50, 200., 1700. );

    h1_mcTruth->Sumw2();
    h1_mcFakes->Sumw2();
    h1_closure->Sumw2();

    tree->Project( "mcTruth", "mt2", "weight*( mt2>200. && ht > 450. && nBJets==2 )" );
    tree->Project( "mcFakes", "mt2", "weight*( mt2>200. && ht > 450. && nBJets==2 && nTrueB==0. && nTrueC==0. )" );

    fillFromTreeRatio( tree, h1_closure, func );

    TCanvas* c1 = new TCanvas( "c1_", "", 600, 600 );
    c1->SetLogy();
    c1->cd();

    TH2D* h2_axes = new TH2D( "axes", "", 10, 200., 1200., 10, 0.0001, 0.5 );
    h2_axes->SetXTitle( "M_{T2} [GeV]");
    h2_axes->SetYTitle( "Events" );
    h2_axes->Draw();

    h1_mcTruth->SetMarkerStyle(20);
    h1_mcTruth->SetMarkerSize(1.1);
    h1_mcTruth->SetMarkerColor(kBlack);
    h1_mcTruth->SetLineColor(kBlack);

    h1_mcFakes->SetLineStyle(2);
    h1_mcFakes->SetLineColor(46);
    //h1_mcFakes->SetLineWidth(2);

    h1_closure->SetMarkerStyle(24);
    h1_closure->SetMarkerSize(1.1);
    h1_closure->SetMarkerColor(kBlack);
    h1_closure->SetLineColor(kBlack);

    h1_mcTruth->Draw("P same");
    h1_mcFakes->Draw("histo same");
    h1_closure->Draw("P same");

    TLegend* legend = new TLegend( 0.53, 0.73, 0.9, 0.9 );
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry( h1_mcTruth, "MC b=2 (all)", "P" );
    legend->AddEntry( h1_mcFakes, "MC b=2 (fakes)", "L" );
    legend->AddEntry( h1_closure, "Extrap. from b=1", "P" );
    legend->Draw("same");

    TPaveText* labelTop = MT2DrawTools::getLabelTop();
    labelTop->Draw("same");

    gPad->RedrawAxis();

    c1->SaveAs(Form("%s/closure.eps", outdir.c_str()));
    c1->SaveAs(Form("%s/closure.pdf", outdir.c_str()));


  } else { // more complex region set


    for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

      if( iR->nBJetsMin()!=1 ) continue; //reweight Nb=1 to get Nb=2
      if( iR->nBJetsMax()!=1 ) continue; //reweight Nb=1 to get Nb=2

      MT2EstimateTree* thisEst_1b = mc->get( *iR );
      TTree* tree_1b = thisEst_1b->tree;

      MT2Region region_2b( iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), 2, 2 );
      MT2EstimateTree* thisEst_2b = mc->get( region_2b );
      if( thisEst_2b==0 ) continue;

      TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

      TH1D* yield_2b = thisEst_2b->yield;

      TH1D* yield_2b_extrapMC = new TH1D( *yield_2b ); // same binning
      yield_2b_extrapMC->Reset();
      std::string newNameMC = "extrap2bMC_" + iR->getName();
      yield_2b_extrapMC->SetName(newNameMC.c_str());

      //TH1D* yield_2b_extrapData = new TH1D( *yield_2b ); // same binning
      //yield_2b_extrapData->Reset();
      //std::string newNameData = "extrap2bData_" + iR->getName();
      //yield_2b_extrapData->SetName(newNameData.c_str());


      fillFromTreeRatio( tree_1b, yield_2b_extrapMC, func );


      TCanvas* c1 = new TCanvas("c1", "", 600, 600);
      c1->cd();

      float xMin = yield_2b->GetXaxis()->GetXmin();
      float xMax = yield_2b->GetXaxis()->GetXmax();
      float yMax_mc = yield_2b->GetMaximum();
      float yMax_extrap = yield_2b_extrapMC->GetMaximum();
      float yMax = (yMax_mc>yMax_extrap) ? yMax_mc: yMax_extrap;
      yMax *= 1.2;

      TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
      h2_axes->SetXTitle( "M_{T2} [GeV]");
      h2_axes->SetYTitle( "Events" );
      h2_axes->Draw();

      yield_2b->SetLineColor(kGray+3);
      yield_2b->SetLineWidth(2.);

      yield_2b_extrapMC->SetLineColor(46);
      yield_2b_extrapMC->SetFillColor(46);
      yield_2b_extrapMC->SetFillStyle(3004);
      yield_2b_extrapMC->SetLineWidth(2);

      yield_2b_extrapMC->SetMarkerColor(46);
      yield_2b_extrapMC->SetMarkerSize(2.);
      yield_2b_extrapMC->SetMarkerStyle(24);

      //yield_2b_extrapData->SetMarkerColor(kBlack);
      //yield_2b_extrapData->SetMarkerSize(1.6);
      //yield_2b_extrapData->SetMarkerStyle(20);

      TPaveText* regionName = new TPaveText(0.5, 0.78, 0.9, 0.88, "brNDC");
      regionName->SetFillColor(0);
      regionName->SetTextAlign(11);
      regionName->SetTextSize(0.035);
      regionName->AddText( region_2b.getNiceNames()[0].c_str() );
      regionName->AddText( region_2b.getNiceNames()[1].c_str() );
      regionName->Draw("same");

      yield_2b->Draw("L same");
      //if( histo_p_data!=0 )
      //  yield_2b_extrapData->Draw("P same");
      yield_2b_extrapMC->Draw("P same");



      TLegend* legend = new TLegend( 0.5, 0.63, 0.9, 0.78);
      legend->SetFillColor(0);
      legend->SetTextSize(0.035);
      legend->AddEntry( yield_2b, "MC Truth", "L" );
      legend->AddEntry( yield_2b_extrapMC, "Extrap. from 1b", "PL");
      legend->Draw("same");

      TPaveText* labelTop = MT2DrawTools::getLabelTop("CMS Simulation, #sqrt{s} = 13 TeV");
      labelTop->Draw("same");

      gPad->RedrawAxis();

      c1->SaveAs(Form("%s/fits2b/closure_%s.eps", cfg.getEventYieldDir().c_str(), region_2b.getName().c_str()) );
      c1->SaveAs(Form("%s/fits2b/closure_%s.pdf", cfg.getEventYieldDir().c_str(), region_2b.getName().c_str()) );
      c1->SaveAs(Form("%s/fits2b/closure_%s.png", cfg.getEventYieldDir().c_str(), region_2b.getName().c_str()) );

      delete c1;
      delete h2_axes;

    } // for regions

    extrap->finalize();

  } // if regions


  return extrap;

}






void fillFromTreeRatio( TTree* tree, TH1D* yield_2b_extrapMC, TF1* func ) {

  float weight;
  tree->SetBranchAddress( "weight", &weight );
  int njets;
  tree->SetBranchAddress( "nJets", &njets );
  int nbjets;
  tree->SetBranchAddress( "nBJets", &nbjets );
  float mt2;
  tree->SetBranchAddress( "mt2", &mt2 );
  float ht;
  tree->SetBranchAddress( "ht", &ht );

  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( njets<2 ) continue;
    if( nbjets!=1 ) continue;  // now correcting 1b
    if( mt2<200. ) continue;
    if( ht<450. ) continue;

    float corr_mc = TMath::Max( func->Eval(njets), 0. );

    yield_2b_extrapMC->Fill( mt2, weight*corr_mc );

  }

}







void compareHistos( MT2Config cfg, const std::string& saveName, TH1D* histo1, TH1D* histo2, TH1D* histo3, TH1D* histo4 ) {


  std::string outdir = cfg.getEventYieldDir() + "/fits2b";

  float xMin = histo1->GetXaxis()->GetXmin();
  float xMax = histo1->GetXaxis()->GetXmax();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();


  float yMax = 0.7;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
  h2_axes->SetXTitle("Number of Jets");
  h2_axes->SetYTitle("2b/1b Ratio");
  h2_axes->Draw();

  histo1->SetMarkerStyle(20);
  histo1->SetMarkerSize(1.6);
  histo1->GetFunction( "line" )->SetLineStyle(1);

  histo2->SetMarkerStyle(24);
  histo2->SetMarkerSize(1.6);

  if( histo3!=0 ) {

    histo1->GetFunction("line")->SetLineColor(kBlack);

    histo2->SetMarkerStyle(20);
    histo2->SetLineColor(46);
    histo2->SetMarkerColor(46);
    histo2->GetFunction("line")->SetLineColor(46);

    histo3->SetMarkerStyle( 21 );
    histo3->SetMarkerSize( 1.6 );
    histo3->SetMarkerColor( 38 );
    histo3->SetLineColor( 38 );
    histo3->GetFunction( "line" )->SetLineColor(38);
  }

  if( histo4!=0 ) {

    histo1->SetMarkerStyle(24);

    histo4->SetMarkerStyle(20);
    histo4->SetLineColor(kBlack);
    histo4->SetMarkerColor(kBlack);
    histo4->GetFunction("line")->SetLineColor(kBlack);

  }

  float yMin_legend = 0.75;
  if( histo3!=0 ) yMin_legend = 0.7;
  if( histo4!=0 ) yMin_legend = 0.65;
  TLegend* legend = new TLegend( 0.2, yMin_legend, 0.6, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( histo1, histo1->GetTitle(), "P" );
  legend->AddEntry( histo2, histo2->GetTitle(), "P" );
  if( histo3!=0 )
    legend->AddEntry( histo3, histo3->GetTitle(), "P" );
  if( histo4!=0 )
    legend->AddEntry( histo4, histo4->GetTitle(), "P" );
  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  if( histo3!=0 )
    histo3->Draw("p same");
  histo1->Draw("p same");
  histo2->Draw("p same");
  if( histo4!=0 )
    histo4->Draw("p same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s.eps", outdir.c_str(), saveName.c_str()) );
  c1->SaveAs( Form("%s/%s.pdf", outdir.c_str(), saveName.c_str()) );

  delete c1;
  delete h2_axes;






}




