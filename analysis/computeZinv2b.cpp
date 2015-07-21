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




void drawMt2VsB( MT2Config cfg, TTree* tree, const std::string& saveName, const std::string& legendTitle, const std::string& additionalSel="");



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

  TFile* file_Zinv = TFile::Open( Form("%s/analyses.root"  , mainDir.c_str()) );
  TTree* tree_Zinv = (TTree*)file_Zinv->Get("ZJets/HT450toInf_j2toInf_b0toInf/tree_ZJets_HT450toInf_j2toInf_b0toInf");

  TFile* file_mc   = TFile::Open( Form("%s/mc.root"  , gammaDir.c_str()) );
  TFile* file_data = TFile::Open( Form("%s/data.root", gammaDir.c_str()) );

  TTree* tree_mc   = (TTree*)file_mc  ->Get("gammaCRtree_loose/HT450toInf_j2toInf_b0toInf/tree_gammaCRtree_loose_HT450toInf_j2toInf_b0toInf");
  TTree* tree_data = (TTree*)file_data->Get("gammaCRtree_loose/HT450toInf_j2toInf_b0toInf/tree_gammaCRtree_loose_HT450toInf_j2toInf_b0toInf");


  drawMt2VsB( cfg, tree_Zinv, "2b_zinv", "Z #rightarrow #nu#nu" );
  drawMt2VsB( cfg, tree_mc  , "2b_prompt", "Prompt Photons", "prompt==2" );
  drawMt2VsB( cfg, tree_mc  , "2b_nip"   , "Fragm. Photons", "prompt==1" );
  drawMt2VsB( cfg, tree_mc  , "2b_fake"  , "Fake   Photons", "prompt==0" );


  return 0;

}



void drawMt2VsB( MT2Config cfg, TTree* tree, const std::string& saveName, const std::string& legendTitle, const std::string& additionalSel) {


  std::string outputdir = cfg.getEventYieldDir() + "/gammaControlRegion/plots2b";
  system(Form("mkdir -p %s", outputdir.c_str()) );


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);

  TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600);
  c1_log->SetLogy();

  int nBins = 100;
  float xMin = 0.;
  float xMax = 1450.;

  TH1D* h1_0b = new TH1D("mt2_0b", "", nBins, xMin, xMax); 
  TH1D* h1_1b = new TH1D("mt2_1b", "", nBins, xMin, xMax); 
  TH1D* h1_2b = new TH1D("mt2_2b", "", nBins, xMin, xMax); 

  h1_0b->Sumw2();
  h1_1b->Sumw2();
  h1_2b->Sumw2();

  if( additionalSel=="" ) {
    tree->Project("mt2_0b", "mt2", "weight*(nBJets==0)");
    tree->Project("mt2_1b", "mt2", "weight*(nBJets==1)");
    tree->Project("mt2_2b", "mt2", "weight*(nBJets==2)");
  } else {
    tree->Project("mt2_0b", "mt2", Form("weight*(nBJets==0 && %s)", additionalSel.c_str()) );
    tree->Project("mt2_1b", "mt2", Form("weight*(nBJets==1 && %s)", additionalSel.c_str()) );
    tree->Project("mt2_2b", "mt2", Form("weight*(nBJets==2 && %s)", additionalSel.c_str()) );
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
  h2_axes->SetXTitle("M_{T2} [GeV]");
  h2_axes->SetYTitle("Normalized to Unity");
  c1->cd();
  h2_axes->Draw();

  TH2D* h2_axes_log = new TH2D( "axes_log", "", 10, xMin, xMax, 10, 0.000001, yMax*2. );
  h2_axes_log->SetXTitle("M_{T2} [GeV]");
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
