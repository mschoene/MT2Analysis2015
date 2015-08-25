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



Double_t bern( Double_t *x, Double_t *par) {

  // par[0] = normalization
  // par[1] = p
  // par[2] = N
  // x[0] = n
  return par[0]*TMath::BinomialI( par[1], par[2], floor(x[0]+0.5) );

}


//Double_t func( Double_t *x, Double_t *par ) {
//
//  return 0.5*par[0]/(1.-par[0])*(x[0]-1.);
//
//}



void drawMt2VsB( MT2Config cfg, TTree* tree, const std::string& suffix, const std::string& legendTitle, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax, const std::string& additionalSel="");
float getP( MT2Config cfg, TTree* tree, const std::string& name, TTree* tree_data=0, const std::string& name_data="Data" );
TH1D* getPHisto( MT2Config cfg, TTree* tree, const std::string& name, TH1D* histo_compare=0, const std::string& name_compare="" );
TH1D* getPHisto2( MT2Config cfg, TTree* tree, const std::string& name, const std::string& niceName );
void getPfromFunc( TF1* line, float& p, float& p_err );
void getPfromFunc21( TF1* line, float& p, float& p_err );
MT2Analysis<MT2Estimate>* compute2bFrom01b( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, float p, float p_err );
MT2Analysis<MT2Estimate>* compute2bFrom01b_2( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, TF1* func );
void fillFromTree( TTree* tree, TH1D* yield_2b_extrapMC, float p, float p_err ); 
void fillFromTree_2( TTree* tree, TH1D* yield_2b_extrapMC, TF1* func );
float getCorrection( int njets, float p );
float getCorrection_1b( int njets, float p );



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

  TH1D* histo_mc, *histo_mc_gjet, *histo_data;

  if( cfg.regionsSet()=="13TeV_inclusive") {

    TFile* file_Zinv = TFile::Open( mcFile.c_str() );
    TTree* tree_Zinv = (TTree*)file_Zinv->Get("ZJets/HT450toInf_j2toInf_b0toInf/tree_ZJets_HT450toInf_j2toInf_b0toInf");

    histo_mc = getPHisto2( cfg, tree_Zinv, "zinv", "Z #rightarrow #nu#nu MC" );

    drawMt2VsB( cfg, tree_Zinv, "zinv", "Z #rightarrow #nu#nu", "mt2", "M_{T2} [GeV]", 100, 0., 1450. );
    //drawMt2VsB( cfg, tree_Zinv, "zinv", "Z #rightarrow #nu#nu", "nJets", "Number of Jets (p_{T}>30 GeV)", 8, 1.5, 9.5 );

    TFile* file_gjet = TFile::Open( Form("%s/mc.root"  , gammaDir.c_str()) );
    TFile* file_data = TFile::Open( Form("%s/data.root", gammaDir.c_str()) );

    TTree* tree_gjet = (TTree*)file_gjet->Get("gammaCRtree/HT450toInf_j2toInf_b0toInf/tree_gammaCRtree_HT450toInf_j2toInf_b0toInf");
    TTree* tree_data = (TTree*)file_data->Get("gammaCRtree/HT450toInf_j2toInf_b0toInf/tree_gammaCRtree_HT450toInf_j2toInf_b0toInf");

    histo_mc_gjet = getPHisto2( cfg, tree_gjet, "gjetMC", "#gamma + Jets MC" );


    TFile* file_zll = TFile::Open( Form("%s/mc.root"  , zllDir.c_str()) );
    TTree* tree_zll = (TTree*)file_zll->Get("zllCR/HT450toInf_j2toInf_b0toInf/tree_zllCR_HT450toInf_j2toInf_b0toInf");
    TFile* file_zllData = TFile::Open( Form("%s/data.root"  , zllDir.c_str()) );
    TTree* tree_zllData = (TTree*)file_zllData->Get("data/HT450toInf_j2toInf_b0toInf/tree_data_HT450toInf_j2toInf_b0toInf");

    getPHisto2( cfg, tree_zll, "zllMC", "Z #rightarrow ll MC" );
    histo_data = getPHisto2( cfg, tree_zllData, "zllData", "Z #rightarrow ll Data" );
    //histo_data = getPHisto( cfg, tree_data, "Data", histo_mc, "MC" );

    //drawMt2VsB( cfg, tree_gjet, "gjet", "Prompt Photons", "mt2", "M_{T2} [GeV]", 100, 0., 1450., "prompt==2" );
    ////drawMt2VsB( cfg, tree_gjet, "gjet", "Prompt Photons", "nJets", "Number of Jets (p_{T}>30 GeV)", 8, 1.5, 9.5, "prompt==2" );

    ////drawMt2VsB( cfg, tree_gjet  , "2b_nip"   , "Fragm. Photons", "prompt==1" );
    ////drawMt2VsB( cfg, tree_gjet  , "2b_fake"  , "Fake   Photons", "prompt==0" );

    TFile* file = TFile::Open("prova.root", "RECREATE");
    file->cd();
    histo_mc->SetName("p");
    histo_mc->Write();
    file->Close();

  } else {

    TFile* file = TFile::Open("prova.root");
    histo_mc = (TH1D*)file->Get("p");
 
  }

  //float p_mc, p_mc_err;
  //getPfromFunc( histo_mc->GetFunction("line"), p_mc, p_mc_err );

  MT2Analysis<MT2EstimateTree>* zinv = MT2Analysis<MT2EstimateTree>::readFromFile(mcFile, "ZJets");
  zinv->setFullName("Z + Jets");
  compute2bFrom01b_2( cfg, zinv, histo_mc->GetFunction("line") );


  return 0;

}




TH1D* getPHisto2( MT2Config cfg, TTree* tree, const std::string& name, const std::string& niceName ) {

  float njetMin = 2;
  float njetMax = 12;
  int nbins_histo = njetMax-njetMin + 1;

  TH1D* h1_p_vs_nj_12 = new TH1D(Form("p_vs_nj_12_%s", name.c_str()), "", nbins_histo, njetMin-0.5, njetMax+0.5 );
  TH1D* h1_p_vs_nj_02 = new TH1D(Form("p_vs_nj_02_%s", name.c_str()), "", nbins_histo, njetMin-0.5, njetMax+0.5 );
  

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

    tree->Project( name_bjets.c_str(), "nBJets", Form("weight*(nJets==%d && mt2>200. && ht>450. )", njet) );

    TF1* bern_12 = new TF1(Form("bern_%d", njet), bern, 0.5, 2.5, 3 );
    bern_12->SetParameter( 0, h1_nbjets->Integral() );
    bern_12->SetParameter( 1, 0.06 );
    bern_12->FixParameter( 2, njet );
    
    TF1* bern_02 = new TF1(Form("bern_%d", njet), bern, -0.5, 2.5, 3 );
    bern_02->SetParameter( 0, h1_nbjets->Integral() );
    bern_02->SetParameter( 1, 0.06 );
    bern_02->FixParameter( 2, njet );
    
    h1_nbjets->Fit( bern_12, "RQN" ); 
    h1_nbjets->Fit( bern_02, "RQN" ); 

    int bin = h1_p_vs_nj_02->FindBin( njet );
    h1_p_vs_nj_12->SetBinContent( bin, bern_12->GetParameter(1) );
    h1_p_vs_nj_12->SetBinError  ( bin, bern_12->GetParError(1) );
    h1_p_vs_nj_02->SetBinContent( bin, bern_02->GetParameter(1) );
    h1_p_vs_nj_02->SetBinError  ( bin, bern_02->GetParError(1) );

    TF1* bern_12_draw = new TF1(Form("bern_%d", njet), bern, nbjet_min-0.5, nbjet_max-0.5, 3 );
    bern_12_draw->SetParameter( 0, bern_12->GetParameter(0) );
    bern_12_draw->SetParameter( 1, bern_12->GetParameter(1) );
    bern_12_draw->SetParameter( 2, bern_12->GetParameter(2) );

    TF1* bern_02_draw = new TF1(Form("bern_%d", njet), bern, nbjet_min-0.5, nbjet_max-0.5, 3 );
    bern_02_draw->SetParameter( 0, bern_02->GetParameter(0) );
    bern_02_draw->SetParameter( 1, bern_02->GetParameter(1) );
    bern_02_draw->SetParameter( 2, bern_02->GetParameter(2) );
    
    bern_12_draw->SetLineColor(46); 
    bern_12_draw->SetLineWidth(2); 

    bern_02_draw->SetLineColor(38); 
    bern_02_draw->SetLineWidth(2); 
    bern_02_draw->SetLineStyle(2); 

    h1_nbjets->SetMarkerStyle(20);
    h1_nbjets->SetMarkerSize(1.5);

    TLegend* legend = new TLegend( 0.6, 0.7, 0.9, 0.9, Form("N(jets) = %d", njet) );
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry( bern_02_draw, "Fit Range 0-2", "L" );
    legend->AddEntry( bern_12_draw, "Fit Range 1-2", "L" );

    TPaveText* labelTop = MT2DrawTools::getLabelTop();
    
    TCanvas* c1_fit = new TCanvas("c1_fit", "", 600, 600);
    c1_fit->cd();
    c1_fit->SetLogy();
    h1_nbjets->Draw("P");
    legend->Draw("same");
    labelTop->Draw("same");
    bern_12_draw->Draw("same");
    bern_02_draw->Draw("same");
    h1_nbjets->Draw("P same");
    gPad->RedrawAxis();
    c1_fit->SaveAs( Form("%s/fit_%d.eps", outdir.c_str(), njet) );
    c1_fit->SaveAs( Form("%s/fit_%d.pdf", outdir.c_str(), njet) );


    delete c1_fit;
    delete bern_12;
    delete bern_02;
    delete bern_12_draw;
    delete bern_02_draw;
    delete h1_nbjets;

  }



  TF1* line_02 = new TF1("line_02", "[0] + [1]*x", njetMin-0.5, njetMax+0.5 );
  line_02->SetLineColor(46);
  h1_p_vs_nj_02->Fit( line_02, "RQ+" );

  TF1* line_12 = new TF1("line", "[0] + [1]*x", njetMin-0.5, njetMax+0.5 );
  line_12->SetLineColor(46);
  line_12->SetLineStyle(2);
  h1_p_vs_nj_12->Fit( line_12, "RQ+" );



  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, njetMin-0.5, njetMax+0.5, 10, 0., 0.16 );
  h2_axes->SetXTitle("Number of Jets");
  h2_axes->SetYTitle("p");
  h2_axes->Draw();

  h1_p_vs_nj_02->SetMarkerStyle(20);
  h1_p_vs_nj_02->SetMarkerSize(1.6);

  h1_p_vs_nj_12->SetMarkerStyle(24);
  h1_p_vs_nj_12->SetMarkerSize(1.6);

  TLegend* legend = new TLegend( 0.5, 0.7, 0.9, 0.9, niceName.c_str() );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( h1_p_vs_nj_02, "Fit Range 0-2", "P" );
  legend->AddEntry( h1_p_vs_nj_12, "Fit Range 1-2", "P" );
  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  h1_p_vs_nj_12->Draw("p same");
  h1_p_vs_nj_02->Draw("p same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/p.eps", outdir.c_str()) );
  c1->SaveAs( Form("%s/p.pdf", outdir.c_str()) );

  delete c1;
  delete h2_axes;

  return h1_p_vs_nj_12;

}

  

TH1D* getPHisto( MT2Config cfg, TTree* tree, const std::string& name, TH1D* histo_compare, const std::string& name_compare ) {

  float njetMin = 2;
  float njetMax = 12;
  int nbins_histo = njetMax-njetMin + 1;
  TH1D* histo = new TH1D(Form("histo_%s", name.c_str()), "", nbins_histo, njetMin-0.5, njetMax+0.5 );
  TH1D* histo_2b1b = new TH1D(Form("histo_2b1b_%s", name.c_str()), "", nbins_histo, njetMin-0.5, njetMax+0.5 );


  float yMin = 0.;
  float yMax = 0.5;
  //float yMax = 0.09;

  

  for( int njet=(int)njetMin; njet<=(int)njetMax; ++njet ) {

    std::string name_0b(Form("%s_%dj_0b", name.c_str(), njet));
    std::string name_1b(Form("%s_%dj_1b", name.c_str(), njet));
    std::string name_2b(Form("%s_%dj_2b", name.c_str(), njet));

    int nBins = 100;
    TH1D* h1_0b = new TH1D( name_0b.c_str(), "", nBins, 0., 20000. );
    TH1D* h1_1b = new TH1D( name_1b.c_str(), "", nBins, 0., 20000. );
    TH1D* h1_2b = new TH1D( name_2b.c_str(), "", nBins, 0., 20000. );

    h1_0b->Sumw2();
    h1_1b->Sumw2();
    h1_2b->Sumw2();

    tree->Project( name_0b.c_str(), "mt2", Form("weight*(nJets==%d && nBJets==0 && mt2>200. && ht>450. )", njet) );
    tree->Project( name_1b.c_str(), "mt2", Form("weight*(nJets==%d && nBJets==1 && mt2>200. && ht>450. )", njet) );
    tree->Project( name_2b.c_str(), "mt2", Form("weight*(nJets==%d && nBJets==2 && mt2>200. && ht>450. )", njet) );
    //tree->Project( name_0b.c_str(), "mt2", Form("weight*(nJets==%d && nBJets==0 && mt2>200. && ht>450. && nTrueB==0 && nTrueC==0)", njet) );
    //tree->Project( name_1b.c_str(), "mt2", Form("weight*(nJets==%d && nBJets==1 && mt2>200. && ht>450. && nTrueB==0 && nTrueC==0)", njet) );


    Double_t int_0b_err;
    Double_t int_0b = h1_0b->IntegralAndError(1, nBins, int_0b_err);
    Double_t int_1b_err;
    Double_t int_1b = h1_1b->IntegralAndError(1, nBins, int_1b_err);
    Double_t int_2b_err;
    Double_t int_2b = h1_2b->IntegralAndError(1, nBins, int_2b_err);

    if( int_1b==0 ) continue;

    float integralRatio = int_1b/int_0b;
    float integralRatioErr = (1./int_0b)*(1./int_0b)*int_1b_err*int_1b_err + int_1b*int_1b/(int_0b*int_0b*int_0b*int_0b)*int_0b_err*int_0b_err;
    integralRatioErr = sqrt(integralRatioErr);

    float integralRatio_2b1b = int_2b/int_1b;
    float integralRatioErr_2b1b = (1./int_1b)*(1./int_1b)*int_2b_err*int_2b_err + int_2b*int_2b/(int_1b*int_1b*int_1b*int_1b)*int_1b_err*int_1b_err;
    integralRatioErr_2b1b = sqrt(integralRatioErr_2b1b);

    //if( integralRatioErr/integralRatio > 0.5 ) continue;

    //float m = integralRatio/(float)njet;
    //float m_err = integralRatioErr/(float)njet;
    //float p = m/(1.+m);
    //float p_err = 1./( (1.+m)*(1.+m) ) * m_err;


    //p *= 100.; // in percent
    //p_err *= 100.; // in percent

    int bin = histo->FindBin(njet);
    histo->SetBinContent( bin, integralRatio );
    histo->SetBinError( bin, integralRatioErr );

    histo_2b1b->SetBinContent( bin, integralRatio_2b1b );
    histo_2b1b->SetBinError( bin, integralRatioErr_2b1b );

    if( integralRatio+integralRatioErr>yMax ) yMax = integralRatio+integralRatioErr;

    //histo->SetBinContent( bin, p );
    //histo->SetBinError( bin, p_err );

    //if( p+p_err>yMax ) yMax = p+p_err;
    //if( p-p_err<yMin ) yMin = p-p_err;

    delete h1_0b;
    delete h1_1b;
    delete h1_2b;

  }

  yMax *= 1.1;

  if( yMax > 1. ) yMax = 1.;


  TF1* line = new TF1("line", "[0] + [1]*x", njetMin-0.5, njetMax+0.5 );
  line->SetLineColor(kRed);
  histo->Fit( line, "RQ" );

  //TF1* line21 = new TF1("line21", func,  njetMin-0.5, njetMax+0.5, 1 );
  TF1* line21 = new TF1("line21", "[0] + [1]*x", 1.5, 4.5 );
  line21->SetLineColor(kRed);
  line21->SetParameter(0, 0.06);
  histo_2b1b->Fit( line21, "R" );

  float p, p_err;
  getPfromFunc( line, p, p_err );

  float p21, p21_err;
  getPfromFunc21( line21, p21, p21_err );


  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, njetMin-0.5, njetMax+0.5, 10, yMin, yMax );
  h2_axes->SetXTitle("Number of Jets");
  //h2_axes->SetYTitle("Events(b=1) / Events(b=0)");
  h2_axes->SetYTitle("Event Ratio");
  h2_axes->Draw();

  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(1.6);

  histo_2b1b->SetMarkerStyle(24);
  histo_2b1b->SetMarkerSize(1.6);

  if( histo_compare!=0 ) {

    //histo_compare->Scale(100.);
    histo_compare->SetMarkerStyle(24);
    histo_compare->SetMarkerSize(1.6);

    TF1* line_compare = histo_compare->GetFunction("line");
    line_compare->SetLineWidth(1);
    line_compare->SetLineStyle(2);

    float p_compare, p_compare_err;
    getPfromFunc(line_compare, p_compare, p_compare_err);

    TLegend* legend = new TLegend( 0.2, 0.73, 0.6, 0.88 );
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry( histo_compare, Form("%s (p=%.2f#pm%.2f%%)", name_compare.c_str(), 100.*p_compare, 100.*p_compare_err), "P");
    legend->AddEntry( histo, Form("%s (p=%.2f#pm%.2f%%)", name.c_str(), 100.*p, 100.*p_err), "P");
    legend->Draw("same");

  } else {

    TPaveText* pLabel = new TPaveText( 0.2, 0.7, 0.55, 0.8, "brNDC" );
    pLabel->SetFillColor(0);
    pLabel->SetTextSize(0.038);
    pLabel->AddText( Form("p = %.2f #pm %.2f %%", p*100. , p_err*100.) );
    pLabel->Draw("same");

  }

  histo->Draw("P same");
  histo_2b1b->Draw("P same");
  if( histo_compare!=0 ) {
    histo_compare->Draw("P same");
  }

  TLegend* legend = new TLegend( 0.2, 0.73, 0.6, 0.88 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( histo     , Form("1b/0b (p=%.2f#pm%.2f%%)", 100.*p, 100.*p_err), "P");
  legend->AddEntry( histo_2b1b, Form("2b/1b (p=%.2f#pm%.2f%%)", 100.*p21, 100.*p21_err), "P");
  legend->Draw("same");

  TPaveText* labelTop;
  if( name=="data" || name=="Data" ) 
    labelTop = MT2DrawTools::getLabelTop(cfg.lumi());
  else
    labelTop = MT2DrawTools::getLabelTop("CMS Simulation, #sqrt{s} = 13 TeV");
  labelTop->Draw("same");
  

  gPad->RedrawAxis();

  //c1->SaveAs(Form("%s/gammaControlRegion/plots2b/p_%s.eps", cfg.getEventYieldDir().c_str(), name.c_str()));
  //c1->SaveAs(Form("%s/gammaControlRegion/plots2b/p_%s.pdf", cfg.getEventYieldDir().c_str(), name.c_str()));
  //c1->SaveAs(Form("%s/gammaControlRegion/plots2b/p_%s.png", cfg.getEventYieldDir().c_str(), name.c_str()));

  c1->SaveAs("p.eps");

  delete c1;
  delete h2_axes;

  //histo->Scale(1./100.);
  //if( histo_compare!=0 )
  //  histo_compare->Scale(1./100.);

  return histo;

}


void getPfromFunc( TF1* line, float& p, float& p_err ) {

  float m = line->GetParameter(1);
  float m_err = line->GetParError(1);

  p = m/(1.+m);
  p_err = 1./( (1.+m)*(1.+m) ) * m_err;

}
 


void getPfromFunc21( TF1* line, float& p, float& p_err ) {

  p = line->GetParameter(0);
  p_err = line->GetParError(0);

  //float m = line->GetParameter(1);
  //float m_err = line->GetParError(1);

  //p = m/(0.5+m);
  //p_err = 0.5/( (0.5+m)*(0.5+m) ) * m_err;

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



MT2Analysis<MT2Estimate>* compute2bFrom01b_2( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, TF1* func ) {



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

    fillFromTree_2( tree, h1_closure, func );

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
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


      fillFromTree_2( tree_1b, yield_2b_extrapMC, func );


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




MT2Analysis<MT2Estimate>* compute2bFrom01b( MT2Config cfg, MT2Analysis<MT2EstimateTree>* mc, float p, float p_err ) {  //TH1D* histo_p, TH1D* histo_p_data ) { 


  std::cout << std::endl;
  std::cout << "-> Starting computation from 0b. Using p = " << p << " +- " << p_err << std::endl;


  std::set<MT2Region> regions = mc->getRegions();

  MT2Analysis<MT2Estimate>* extrap = new MT2Analysis<MT2Estimate>( mc->getName() + "_extrap", regions );



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

    tree->Project( "mcFakes", "mt2", "weight*( mt2>200. && ht > 450. && nBJets==2 && nTrueB==0. && nTrueC==0. )" );
    tree->Project( "mcTruth", "mt2", "weight*( mt2>200. && ht > 450. && nBJets==2 )" );

    fillFromTree( tree, h1_closure, p, p_err );

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
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
    legend->AddEntry( h1_closure, "Extrap. from b=0", "P" );
    legend->Draw("same");

    TPaveText* labelTop = MT2DrawTools::getLabelTop();
    labelTop->Draw("same");

    gPad->RedrawAxis();

    c1->SaveAs("closure.eps");

    TFile* fileprova = TFile::Open("prova.root", "recreate");
    fileprova->cd();
    h1_mcTruth->Write();
    h1_closure->Write();
    fileprova->Close();



  } else { // more complex regions sets
     


    for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

      if( iR->nBJetsMin()!=0 ) continue; //reweight Nb=0 to get Nb=2
      if( iR->nBJetsMax()!=0 ) continue; //reweight Nb=0 to get Nb=2

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

      TH1D* yield_2b_extrapMC = new TH1D( *yield_2b ); // same binning
      yield_2b_extrapMC->Reset();
      std::string newNameMC = "extrap2bMC_" + iR->getName();
      yield_2b_extrapMC->SetName(newNameMC.c_str());

      TH1D* yield_2b_extrapData = new TH1D( *yield_2b ); // same binning
      yield_2b_extrapData->Reset();
      std::string newNameData = "extrap2bData_" + iR->getName();
      yield_2b_extrapData->SetName(newNameData.c_str());


      fillFromTree( tree_0b, yield_2b_extrapMC, p, p_err );  //histo_p, histo_p_data );


      TCanvas* c1 = new TCanvas("c1", "", 600, 600);
      c1->cd();

      float xMin = yield_2b->GetXaxis()->GetXmin();
      float xMax = yield_2b->GetXaxis()->GetXmax();
      float yMax = 1.15*yield_2b->GetMaximum();

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

      yield_2b->Draw("L same");
      //if( histo_p_data!=0 )
      //  yield_2b_extrapData->Draw("P same");
      yield_2b_extrapMC->Draw("P same");



      TLegend* legend = new TLegend( 0.5, 0.63, 0.9, 0.78);
      legend->SetFillColor(0);
      legend->SetTextSize(0.035);
      legend->AddEntry( yield_2b, "MC Truth", "L" );
      //if( histo_p_data!=0 ) {
      //  legend->AddEntry( yield_2b_extrapMC, "Extrap. (MC)", "P");
      //  legend->AddEntry( yield_2b_extrapData, "Extrap. (Data)", "P");
      //} else {
        legend->AddEntry( yield_2b_extrapMC, "Extrap. from 0b", "PL");
        //legend->AddEntry( yield_2b_extrapMC, "Extrap. from 0b", "FLP");
      //}
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

  } // if regions

  
  return extrap;

}


void fillFromTree_2( TTree* tree, TH1D* yield_2b_extrapMC, TF1* func ) {

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

    float corr_mc = getCorrection_1b( njets, func->Eval(njets) );

    yield_2b_extrapMC->Fill( mt2, weight*corr_mc );

  }

}




void fillFromTree( TTree* tree, TH1D* yield_2b_extrapMC, float p, float p_err ) {  //TH1D* histo_p, TH1D* histo_p_data ) {


  TH1D* yield_2b_extrapMCUp   = new TH1D(*yield_2b_extrapMC);
  TH1D* yield_2b_extrapMCDown = new TH1D(*yield_2b_extrapMC);

  yield_2b_extrapMCUp->SetName(Form("%sUp"  , yield_2b_extrapMC->GetName()));
  yield_2b_extrapMCDown->SetName(Form("%sDown", yield_2b_extrapMC->GetName()));

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
    if( nbjets!=0 ) continue;
    if( mt2<200. ) continue;
    if( ht<450. ) continue;

    float corr_mc       = getCorrection( njets, p );
    float corr_mcUp     = getCorrection( njets, p+p_err );
    float corr_mcDown   = getCorrection( njets, p-p_err );

    yield_2b_extrapMC    ->Fill( mt2, weight*corr_mc );
    yield_2b_extrapMCUp  ->Fill( mt2, weight*corr_mcUp );
    yield_2b_extrapMCDown->Fill( mt2, weight*corr_mcDown );

  }


  //for( int ibin=0; ibin<yield_2b_extrapMC->GetNbinsX(); ++ibin ) {

  //  float x     = yield_2b_extrapMC    ->GetBinContent(ibin+1);
  //  float xUp   = yield_2b_extrapMCUp  ->GetBinContent(ibin+1);
  //  float xDown = yield_2b_extrapMCDown->GetBinContent(ibin+1);
 
  //  float xSyst = (xUp-x > x-xDown) ? xUp-x : x-xDown;
  //  float xStat = yield_2b_extrapMC->GetBinError(ibin+1);
  //  float xErr = sqrt( xStat*xStat + xSyst*xSyst );

  //  yield_2b_extrapMC->SetBinError(ibin+1, xErr);

  //}
    

}


float getCorrection( int njets, float p ) {

  float corr = TMath::Binomial(njets,2) * p*p / ( (1.-p)*(1.-p) ); //TMath::Power( 1.-p, njets-2 ); 
  //float corr = (float)njets*(njets-1.)*p*p/((1.-p)*(1.-p)); 
  //float corr = 0.5*(float)njets*(njets-1.)*p*p/((1.-p)*(1.-p)); 

  return corr;

}



float getCorrection_1b( int njets, float p ) {

  float corr = 0.5*(njets-1.)*p/(1.-p);

  return corr;

}
