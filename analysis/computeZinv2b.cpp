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
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"




void drawMt2VsB( MT2Config cfg, TTree* tree, const std::string& suffix, const std::string& legendTitle, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax, const std::string& additionalSel="");
float getP( MT2Config cfg, TTree* tree, const std::string& name, TTree* tree_data=0, const std::string& name_data="Data" );
TH1D* getHisto( MT2Config cfg, TTree* tree, const std::string& name, TGraphErrors* gr_compare=0, const std::string& name_compare="" );
//TGraphErrors* getGraph( MT2Config cfg, TTree* tree, const std::string& name, TGraphErrors* gr_compare=0, const std::string& name_compare="" );
void getPFromHisto( TH1D* gr, float& p, float& p_err );
MT2Analysis<MT2Estimate>* compute2bFrom01b( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, TH1D* histo_p ); //float p_mc, float p_data=0., MT2Analysis<MT2EstimateTree>* data=0 );
//MT2Analysis<MT2Estimate>* compute2bFrom01b( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, float p_mc, float p_data=0., MT2Analysis<MT2EstimateTree>* data=0 );
void fillFromTree( TTree* tree, TH1D* yield_2b_extrapMC, TH1D* yield_2b_extrapData, TH1D* histo_p );
//void fillFromTree( TTree* tree, TH1D* yield_2b_extrapMC, TH1D* yield_2b_extrapData, float p_mc, float p_data );



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


  //TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  
  std::string mainDir  = cfg.getEventYieldDir();
  std::string gammaDir = mainDir + "/gammaControlRegion";

  std::string mcFile = mainDir + "/analyses.root";


  float p_mc   = 0.036; // default
  //float p_mc   = 0.019; // default
  float p_data = 0.; // default
  TH1D* histo_mc;

  if( cfg.regionsSet()=="13TeV_inclusive") {

    TFile* file_Zinv = TFile::Open( mcFile.c_str() );
    TTree* tree_Zinv = (TTree*)file_Zinv->Get("ZJets/HT450toInf_j2toInf_b0toInf/tree_ZJets_HT450toInf_j2toInf_b0toInf");

    //TGraphErrors* gr_mc = getGraph( cfg, tree_Zinv, "mc" );
    histo_mc = getHisto( cfg, tree_Zinv, "mc" );
    float p_mc, p_mc_err;
    getPFromHisto(histo_mc, p_mc, p_mc_err);

    drawMt2VsB( cfg, tree_Zinv, "zinv", "Z #rightarrow #nu#nu", "mt2", "M_{T2} [GeV]", 100, 0., 1450. );
    //drawMt2VsB( cfg, tree_Zinv, "zinv", "Z #rightarrow #nu#nu", "nJets", "Number of Jets (p_{T}>30 GeV)", 8, 1.5, 9.5 );

    TFile* file_gjet = TFile::Open( Form("%s/mc.root"  , gammaDir.c_str()) );
    TFile* file_data = TFile::Open( Form("%s/data.root", gammaDir.c_str()) );

    TTree* tree_gjet = (TTree*)file_gjet->Get("gammaCRtree_loose/HT450toInf_j2toInf_b0toInf/tree_gammaCRtree_loose_HT450toInf_j2toInf_b0toInf");
    TTree* tree_data = (TTree*)file_data->Get("gammaCRtree_loose/HT450toInf_j2toInf_b0toInf/tree_gammaCRtree_loose_HT450toInf_j2toInf_b0toInf");

    //TGraphErrors* gr_data = getGraph( cfg, tree_data, "Data", gr_mc, "MC" );
    //float p_data, p_data_err;
    //getPFromGraph(gr_data, p_data, p_data_err);

    drawMt2VsB( cfg, tree_gjet, "gjet", "Prompt Photons", "mt2", "M_{T2} [GeV]", 100, 0., 1450., "prompt==2" );
    //drawMt2VsB( cfg, tree_gjet, "gjet", "Prompt Photons", "nJets", "Number of Jets (p_{T}>30 GeV)", 8, 1.5, 9.5, "prompt==2" );

    //drawMt2VsB( cfg, tree_gjet  , "2b_nip"   , "Fragm. Photons", "prompt==1" );
    //drawMt2VsB( cfg, tree_gjet  , "2b_fake"  , "Fake   Photons", "prompt==0" );

    TFile* file = TFile::Open("prova.root", "RECREATE");
    file->cd();
    histo_mc->Write();
    file->Close();

  } else {

    TFile* file = TFile::Open("prova.root");
    histo_mc = (TH1D*)file->Get("histo_mc");
 
  }

  MT2Analysis<MT2EstimateTree>* zinv = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "ZJets");
  zinv->setFullName("Z + Jets");
  compute2bFrom01b( cfg, zinv, histo_mc ); //p_mc, p_data );


  return 0;

}




TH1D* getHisto( MT2Config cfg, TTree* tree, const std::string& name, TGraphErrors* gr_compare, const std::string& name_compare ) {


  //int i=0;
  //TGraphErrors* graph = new TGraphErrors(0);
  TH1D* histo = new TH1D(Form("histo_%s", name.c_str()), "", 7, 1.5, 8.5);
  

  float yMin = 0.1;
  float yMax = 0.2;

  float mean_w = 0.;
  float sumw = 0.;

  //TH1D* h1 = new TH1D("prova", "", 100, 0., 0.5);
  

  for( unsigned njet=2; njet<=8; ++njet ) {

    std::string name_0b(Form("%s_%dj_0b", name.c_str(), njet));
    std::string name_1b(Form("%s_%dj_1b", name.c_str(), njet));

    int nBins = 100;
    TH1D* h1_0b = new TH1D( name_0b.c_str(), "", nBins, 0., 20000. );
    TH1D* h1_1b = new TH1D( name_1b.c_str(), "", nBins, 0., 20000. );

    h1_0b->Sumw2();
    h1_1b->Sumw2();

    tree->Project( name_0b.c_str(), "mt2", Form("weight*(nJets==%d && nBJets==0 && mt2>200. && ht>450.)", njet) );
    tree->Project( name_1b.c_str(), "mt2", Form("weight*(nJets==%d && nBJets==1 && mt2>200. && ht>450.)", njet) );

    //float int_0b = h1_0b->Integral();
    //float int_1b = h1_1b->Integral();
    Double_t int_0b_err;
    Double_t int_0b = h1_0b->IntegralAndError(1, nBins, int_0b_err);
    Double_t int_1b_err;
    Double_t int_1b = h1_1b->IntegralAndError(1, nBins, int_1b_err);

    if( int_1b==0 ) continue;

    float integralRatio = int_1b/int_0b;
    float integralRatioErr = (1./int_0b)*(1./int_0b)*int_1b_err*int_1b_err + int_1b*int_1b/(int_0b*int_0b*int_0b*int_0b)*int_0b_err*int_0b_err;
    integralRatioErr = sqrt(integralRatioErr);

    if( integralRatioErr/integralRatio > 0.5 ) continue;

    int bin = histo->FindBin(njet);
    histo->SetBinContent( bin, integralRatio/(float)njet );
    histo->SetBinError( bin, integralRatioErr );

    //graph->SetPoint( i, njet, integralRatio/(float)njet );
    //graph->SetPointError( i, 0., integralRatioErr );
    //i++;

    //h1->Fill(integralRatio/(float)njet);
    float w = 1/integralRatioErr;
    mean_w += integralRatio/(float)njet*w;
    sumw += w;

    if( integralRatio+integralRatioErr>yMax ) yMax = integralRatio+integralRatioErr;
    if( integralRatio-integralRatioErr<yMin ) yMin = integralRatio-integralRatioErr;

    delete h1_0b;
    delete h1_1b;

  }

  mean_w /= sumw;

std::cout << "**** mean_w: " << mean_w << std::endl;
//std::cout << "**** h1->GetMean(): " <<  h1->GetMean() << std::endl;

  yMin *= 0.9;
  yMax *= 1.1;

yMin = 0.;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 1.5, 8.5, 10, yMin, yMax );
  h2_axes->SetXTitle("Number of Jets");
  h2_axes->SetYTitle("Events(b=1) / Events(b=0)");
  h2_axes->Draw();

  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(1.6);
  //graph->SetMarkerStyle(20);
  //graph->SetMarkerSize(1.6);

  //if( gr_compare!=0 ) {
  //  gr_compare->SetMarkerStyle(24);
  //  gr_compare->SetMarkerSize(1.6);
  //  TF1* line_compare = gr_compare->GetFunction("line");
  //  line_compare->SetLineStyle(2);
  //}

  TF1* line = new TF1("line", "[0] + [1]*x", 1.5, 8.5 );
  line->SetLineColor(kRed);
  histo->Fit( line, "RQ" );
  //graph->Fit( line, "RQ" );

  float p, p_err;
  getPFromHisto( histo, p, p_err );
  //getPFromGraph( graph, p, p_err );


  //if( gr_compare!=0 ) {

  //  float p_compare, p_compare_err;
  //  getPFromGraph( gr_compare, p_compare, p_compare_err );

  //  TLegend* legend = new TLegend( 0.17, 0.75, 0.5, 0.9 );
  //  legend->SetFillColor(0);
  //  legend->SetTextSize(0.035);
  //  legend->AddEntry( gr_compare, Form("%s (p=%.1f#pm%.1f%%)", name_compare.c_str(), 100.*p_compare, 100.*p_compare_err), "P");
  //  legend->AddEntry( graph, Form("%s (p=%.2f#pm%.2f%%)", name.c_str(), 100.*p, 100.*p_err), "P");
  //  legend->Draw("same");

  //} else {

    TPaveText* pText = new TPaveText( 0.55, 0.2, 0.9, 0.25, "brNDC" );
    pText->SetTextSize(0.038);
    pText->SetFillColor(0);
    pText->AddText( Form("p = %.2f #pm %.2f%%", 100.*p, 100.*p_err) );
    pText->Draw("P same");

  //}


  histo->Draw("P same");
  //graph->Draw("P same");
  //if( gr_compare!=0 ) {
  //  gr_compare->Draw("P same");
  //}

  TPaveText* labelTop;
  if( name=="data" || name=="Data" ) 
    labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
  else
    labelTop = MT2DrawTools::getLabelTop("CMS Simulation, #sqrt{s} = 13 TeV");
  labelTop->Draw("same");
  

  gPad->RedrawAxis();

  c1->SaveAs(Form("%s/gammaControlRegion/plots2b/p_%s.eps", cfg.getEventYieldDir().c_str(), name.c_str()));
  c1->SaveAs(Form("%s/gammaControlRegion/plots2b/p_%s.pdf", cfg.getEventYieldDir().c_str(), name.c_str()));
  c1->SaveAs(Form("%s/gammaControlRegion/plots2b/p_%s.png", cfg.getEventYieldDir().c_str(), name.c_str()));

  std::cout << "q: " << line->GetParameter(0) << std::endl;
  std::cout << "m: " << line->GetParameter(1) << std::endl;
  std::cout << "p(" << name << "): " << p << " +- " << p_err << std::endl;

  delete c1;
  delete h2_axes;

  return histo;
  //return graph;

}


void getPFromHisto( TH1D* histo, float& p, float& p_err ) {
//void getPFromGraph( TGraphErrors* gr, float& p, float& p_err ) {

  TF1* line = histo->GetFunction("line");

  float m = line->GetParameter(1);
  float m_err = line->GetParError(1);

  p = m/(1.+m);
  p_err = 1./( (1.+m)*(1.+m) ) * m_err;

}

 

void drawMt2VsB( MT2Config cfg, TTree* tree, const std::string& suffix, const std::string& legendTitle, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax, const std::string& additionalSel) {


  std::string outputdir = cfg.getEventYieldDir() + "/gammaControlRegion/plots2b";
  system(Form("mkdir -p %s", outputdir.c_str()) );


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);

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



MT2Analysis<MT2Estimate>* compute2bFrom01b( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, TH1D* histo_p ) { //float p_mc, float p_data, MT2Analysis<MT2EstimateTree>* data ) {


  std::set<MT2Region> regions = mc->getRegions();

  MT2Analysis<MT2Estimate>* extrap = new MT2Analysis<MT2Estimate>( mc->getName() + "_extrap", regions );


  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()!=0 ) continue; //reweight Nb=0 to get Nb=2

    MT2EstimateTree* thisEst_0b = mc->get( *iR );
    TTree* tree_0b = thisEst_0b->tree;

    MT2Region region_1b( iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), 1, 1 );
    MT2EstimateTree* thisEst_1b = mc->get( region_1b );
    if( thisEst_1b==0 ) continue;

    MT2Region region_2b( iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), 2, 2 );
    MT2EstimateTree* thisEst_2b = mc->get( region_2b );
    if( thisEst_2b==0 ) continue;

    TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

    TH1D* yield_2b = thisEst_2b->yield;

    //TH1D* yield_2b_extrap = extrap->get( *iR )->yield;
    //TH1D* yield_2b_extrap = new TH1D( *(extrap->get( *iR )->yield) );
    //yield_2b_extrap->SetName("tt");
    TH1D* yield_2b_extrapMC = new TH1D( *yield_2b ); // same binning
    yield_2b_extrapMC->Reset();
    std::string newNameMC = "extrap2bMC_" + iR->getName();
    yield_2b_extrapMC->SetName(newNameMC.c_str());

    TH1D* yield_2b_extrapData = new TH1D( *yield_2b ); // same binning
    yield_2b_extrapData->Reset();
    std::string newNameData = "extrap2bData_" + iR->getName();
    yield_2b_extrapData->SetName(newNameData.c_str());

    TH1D* yield_2b_MC = new TH1D( *yield_2b ); // same binning
    yield_2b_MC->Reset();
    std::string newNameMCMC = "MC_" + iR->getName();
    yield_2b_MC->SetName(newNameMCMC.c_str());
    thisEst_2b->tree->Project( newNameMCMC.c_str(), "mt2", "weight");

    fillFromTree( tree_0b, yield_2b_extrapMC, yield_2b_extrapData, histo_p );
    //fillFromTree( tree_0b, yield_2b_extrapMC, yield_2b_extrapData, p_mc, p_data );


    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();

    float xMin = yield_2b->GetXaxis()->GetXmin();
    float xMax = yield_2b->GetXaxis()->GetXmax();
    float yMax = 1.15*yield_2b->GetMaximum();

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle( "M_{T2} [GeV]");
    h2_axes->SetYTitle( "Events" );
    h2_axes->Draw();

    //yield_2b->SetLineColor(kRed);
    //yield_2b->SetLineWidth(2.);

    yield_2b_MC->SetLineColor(kRed);
    yield_2b_MC->SetLineWidth(2.);

    yield_2b_extrapMC->SetMarkerColor(kBlue);
    yield_2b_extrapMC->SetMarkerSize(1.6);
    yield_2b_extrapMC->SetMarkerStyle(24);

    yield_2b_extrapData->SetMarkerColor(kBlack);
    yield_2b_extrapData->SetMarkerSize(1.6);
    yield_2b_extrapData->SetMarkerStyle(20);

    TPaveText* regionName = new TPaveText(0.5, 0.78, 0.9, 0.88, "brNDC");
    regionName->SetFillColor(0);
    regionName->SetTextAlign(11);
    regionName->SetTextSize(0.035);
    regionName->AddText( region_2b.getNiceNames()[0].c_str() );
    regionName->AddText( region_2b.getNiceNames()[1].c_str() );
    regionName->Draw("same");

    yield_2b_MC->Draw("L same");
    //yield_2b->Draw("L same");
    yield_2b_extrapMC->Draw("P same");
    //if( p_data>0. )
    //  yield_2b_extrapData->Draw("P same");



    TLegend* legend = new TLegend( 0.5, 0.6, 0.9, 0.75);
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry( yield_2b, "MC Truth", "L" );
    legend->AddEntry( yield_2b_extrapMC, "Extrapolation", "P");
    //legend->AddEntry( yield_2b_extrapMC, Form("Extrap. (p=%.2f%%)", 100.*p_mc), "P" );
    legend->Draw("same");

    TPaveText* labelTop = MT2DrawTools::getLabelTop("CMS Simulation, #sqrt{s} = 13 TeV");
    labelTop->Draw("same");

    gPad->RedrawAxis();

    c1->SaveAs(Form("%s/closure_%s.eps", cfg.getEventYieldDir().c_str(), region_2b.getName().c_str()) );
    c1->SaveAs(Form("%s/closure_%s.pdf", cfg.getEventYieldDir().c_str(), region_2b.getName().c_str()) );
    c1->SaveAs(Form("%s/closure_%s.png", cfg.getEventYieldDir().c_str(), region_2b.getName().c_str()) );

    delete c1;
    delete h2_axes;

  } // for regions


  extrap->finalize();
  
  return extrap;

}





void fillFromTree( TTree* tree, TH1D* yield_2b_extrapMC, TH1D* yield_2b_extrapData, TH1D* histo_p ) {
//void fillFromTree( TTree* tree, TH1D* yield_2b_extrapMC, TH1D* yield_2b_extrapData, float p_mc, float p_data ) {

  float weight;
  tree->SetBranchAddress( "weight", &weight );
  int njets;
  tree->SetBranchAddress( "nJets", &njets );
  int nbjets;
  tree->SetBranchAddress( "nBJets", &nbjets );
  float mt2;
  tree->SetBranchAddress( "mt2", &mt2 );

  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    float p_mc = histo_p->GetBinContent( histo_p->FindBin(njets) );

    float correction_mc = 0.5*njets*(njets-1.)*p_mc*p_mc*pow(1.-p_mc, njets-2);
    float correction_data = correction_mc;
    //float correction_data = 0.5*njets*(njets-1.)*p_data*p_data*pow(1.-p_data, njets-2);

    //float correction_mc = p_mc*njets;

    yield_2b_extrapMC  ->Fill( mt2, weight*correction_mc );
    yield_2b_extrapData->Fill( mt2, weight*correction_data );

  }

}


/*
MT2Analysis<MT2Estimate>* compute2bFrom01b( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, MT2Analysis<MT2EstimateTree>* data ) {


  std::set<MT2Region> regions = mc->getRegions();

  MT2Analysis<MT2Estimate>* extrap = new MT2Analysis<MT2Estimate>( mc->getName() + "_extrap", regions );


  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    if( iR->nBJetsMin()!=2 ) continue;

    MT2EstimateTree* thisEst_2b = mc->get( *iR );

    MT2Region region_0b( iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), 0, 0 );
    MT2Region region_1b( iR->htMin(), iR->htMax(), iR->nJetsMin(), iR->nJetsMax(), 1, 1 );

    MT2EstimateTree* thisEst_0b = mc->get( region_0b );
    MT2EstimateTree* thisEst_1b = mc->get( region_1b );
    MT2EstimateTree* thisEstData_1b = (data!=0) ? data->get( region_1b ) : 0;

    TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

    TH1D* yield_2b = thisEst_2b->yield;

    //TH1D* yield_2b_extrap = extrap->get( *iR )->yield;
    //TH1D* yield_2b_extrap = new TH1D( *(extrap->get( *iR )->yield) );
    //yield_2b_extrap->SetName("tt");
    TH1D* yield_2b_extrap = new TH1D( *yield_2b ); // same binning
    std::string newName = "extrap2b_" + iR->getName();
    yield_2b_extrap->SetName(newName.c_str());
    
    float transferFactor = thisEst_1b->yield->Integral()/thisEst_0b->yield->Integral(); 

    if( thisEstData_1b!=0 ) 
      thisEstData_1b->tree->Project(yield_2b_extrap->GetName(), "mt2", Form("%f*weight", transferFactor) );
    else {
      thisEst_1b    ->tree->Project(newName.c_str(), "mt2", Form("%f*weight", transferFactor) );
    }


    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();

    yield_2b->SetLineColor(kBlack);
    yield_2b_extrap->SetLineColor(kRed);

    yield_2b->Draw("histo");
    yield_2b_extrap->Draw("histo same");

    c1->SaveAs(Form("%s/closure_%s.eps", cfg.getEventYieldDir().c_str(), iR->getName().c_str()) );
    //c1->SaveAs(Form("%s/closure_%s.eps", cfg.getEventYieldDir().c_str(), iR->getName().c_str()) );
    //c1->SaveAs(Form("%s/closure_%s.eps", cfg.getEventYieldDir().c_str(), iR->getName().c_str()) );

    delete c1;

  } // for regions


  extrap->finalize();
  
  return extrap;

}*/
