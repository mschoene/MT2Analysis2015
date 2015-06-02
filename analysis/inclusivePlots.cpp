#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"

#include "TRandom3.h"

#define mt2_cxx
#include "interface/mt2.h"

double lumi=4.; //fb-1

int main( int argc, char* argv[] ) {

  MT2DrawTools::setStyle();

  if( argc > 2 ) {
    std::cout << "USAGE: ./inclusivePlots [samplesFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  
  std::string sampleName;
  if( argc == 2)
    sampleName = argv[1];
  else 
    sampleName = "PHYS14_v5_plots";

  std::string samplesFileName = "../samples/samples_" + sampleName + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName);
  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }

  std::vector<MT2Sample> bSamples = MT2Sample::loadSamples(samplesFileName, 100, 699);
  std::vector<MT2Sample> sSamples = MT2Sample::loadSamples(samplesFileName, 1000);
  
  std::string outputdir = "Plots/"+sampleName;
  system(Form("mkdir -p %s", outputdir.c_str()));


  // DECLARING STACK HISTOGRAM - ONE PER VARIABLE
  THStack* h_njets_stack  = new THStack("h_njets_stack", "");
  THStack* h_nbjets_stack = new THStack("h_nbjets_stack", "");
  THStack* h_nleps_stack  = new THStack("h_nleps_stack", "");
  THStack* h_ht_stack     = new THStack("h_ht_stack", "");
  THStack* h_met_stack    = new THStack("h_met_stack", "");
  THStack* h_mt2_stack    = new THStack("h_mt2_stack", "");
  THStack* h_deltaPhiMin_stack  = new THStack("h_deltaPhiMin_stack", "");
  THStack* h_diffMetMht_stack  = new THStack("h_diffMetMht_stack", "");


  // Declaring histograms - one per sample
  int nSamples = fSamples.size(); 
  std::cout<< "nSamples " << nSamples << std::endl;

  TH1F* h_njets[nSamples+1];  
  TH1F* h_nbjets[nSamples+1];
  TH1F* h_nleps[nSamples+1];
  TH1F* h_ht[nSamples+1];
  TH1F* h_met[nSamples+1];
  TH1F* h_mt2[nSamples+1];
  TH1F* h_deltaPhiMin[nSamples+1];
  TH1F* h_diffMetMht[nSamples+1];
  
  for( unsigned i=0; i<fSamples.size(); ++i ){
    
    //Initializing histograms:    
    int nbins;
    float xmin =0., xmax;

    nbins = 12;
    xmax  = 12.;
    std::string name_njets = "njets"+fSamples[i].name;
    h_njets[i] = new TH1F( name_njets.c_str(), "number of jets", nbins, xmin, xmax );
    h_njets[i]->GetXaxis()->SetTitle("N(jets)");
    h_njets[i]->GetYaxis()->SetTitle("Events");
  

    nbins = 6;
    xmax  = 6.;
    std::string name_nbjets="nbjets"+fSamples[i].name;
    h_nbjets[i] = new TH1F(name_nbjets.c_str(), "number of b-jets", nbins, xmin, xmax );
    h_nbjets[i]->GetXaxis()->SetTitle("N(b-jets)");
    h_nbjets[i]->GetYaxis()->SetTitle("Events");
  
    nbins = 10;
    xmax  = 10.;
    std::string name_nleps="nleps"+fSamples[i].name;
    h_nleps[i] = new TH1F(name_nleps.c_str(), "number of leptons", nbins, xmin, xmax );
    h_nleps[i]->GetXaxis()->SetTitle("N(leptons)");
    h_nleps[i]->GetYaxis()->SetTitle("Events");
    
    nbins = 120;
    xmax  = 3000.;
    std::string name_ht="ht"+fSamples[i].name;
    h_ht[i] = new TH1F(name_ht.c_str(), "H_{T}", nbins, xmin, xmax );
    h_ht[i]->GetXaxis()->SetTitle("H_{T} [GeV]");
    h_ht[i]->GetYaxis()->SetTitle("Events/25 GeV");
    h_ht[i]->SetMinimum(1.e-1);
    
    nbins = 60;
    xmax  = 1500.;
    std::string name_met="met"+fSamples[i].name;
    h_met[i] = new TH1F(name_met.c_str(), "E_{T}^{miss}", nbins, xmin, xmax );
    h_met[i]->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    h_met[i]->GetYaxis()->SetTitle("Events/25 GeV");
    h_met[i]->SetMinimum(1.e-1);

    nbins = 60;
    xmax  = 1500.;
    std::string name_mt2="mt2"+fSamples[i].name;
    h_mt2[i] = new TH1F(name_mt2.c_str(), "M_{T2}", nbins, xmin, xmax );
    h_mt2[i]->GetXaxis()->SetTitle("M_{T2} [GeV]");
    h_mt2[i]->GetYaxis()->SetTitle("Events/25 GeV");
    h_mt2[i]->SetMinimum(1.e-1);

    nbins = 32;
    xmax  = 3.2;
    std::string name_deltaPhiMin="mindPhi"+fSamples[i].name;
    h_deltaPhiMin[i] = new TH1F(name_deltaPhiMin.c_str(), "min #delta #Phi (4 jets, E_{T}^{miss})", nbins, xmin, xmax );
    h_deltaPhiMin[i]->GetXaxis()->SetTitle("min #delta #Phi (4 jets, E_{T}^{miss})");
    h_deltaPhiMin[i]->GetYaxis()->SetTitle("Events");
    
    nbins = 20;
    xmax  = 1.;
    std::string name_diffMetMht="diffMetMht"+fSamples[i].name;
    h_diffMetMht[i] = new TH1F(name_diffMetMht.c_str(), "(E_{T}^{miss} - H_{T}^{miss})/E_{T}^{miss}", nbins, xmin, xmax );
    h_diffMetMht[i]->GetXaxis()->SetTitle("(E_{T}^{miss} - H_{T}^{miss})/E_{T}^{miss}");
    h_diffMetMht[i]->GetYaxis()->SetTitle("Events");

    if( fSamples[i].id >= 100 && fSamples[i].id < 200 ){
      h_njets[i]->SetLineColor(401);
      h_njets[i]->SetFillColor(401);
      h_nbjets[i]->SetLineColor(401);
      h_nbjets[i]->SetFillColor(401);
      h_nleps[i]->SetLineColor(401);
      h_nleps[i]->SetFillColor(401);
      h_ht[i]->SetLineColor(401);
      h_ht[i]->SetFillColor(401);
      h_met[i]->SetLineColor(401);
      h_met[i]->SetFillColor(401);
      h_mt2[i]->SetLineColor(401);
      h_mt2[i]->SetFillColor(401);
      h_deltaPhiMin[i]->SetLineColor(401);
      h_deltaPhiMin[i]->SetFillColor(401);
      h_diffMetMht[i]->SetLineColor(401);
      h_diffMetMht[i]->SetFillColor(401);
    }
    else if( fSamples[i].id >= 500 && fSamples[i].id < 600 ){
      h_njets[i]->SetLineColor(417);
      h_njets[i]->SetFillColor(417);
      h_nbjets[i]->SetLineColor(417);
      h_nbjets[i]->SetFillColor(417);
      h_nleps[i]->SetLineColor(417);
      h_nleps[i]->SetFillColor(417);
      h_ht[i]->SetLineColor(417);
      h_ht[i]->SetFillColor(417);
      h_met[i]->SetLineColor(417);
      h_met[i]->SetFillColor(417);
      h_mt2[i]->SetLineColor(417);
      h_mt2[i]->SetFillColor(417);
      h_deltaPhiMin[i]->SetLineColor(417);
      h_deltaPhiMin[i]->SetFillColor(417);
      h_diffMetMht[i]->SetLineColor(417);
      h_diffMetMht[i]->SetFillColor(417);
    }
    else if( fSamples[i].id >= 600 && fSamples[i].id < 700 ){
      h_njets[i]->SetLineColor(419);
      h_njets[i]->SetFillColor(419);
      h_nbjets[i]->SetLineColor(419);
      h_nbjets[i]->SetFillColor(419);
      h_nleps[i]->SetLineColor(419);
      h_nleps[i]->SetFillColor(419);
      h_ht[i]->SetLineColor(419);
      h_ht[i]->SetFillColor(419);
      h_met[i]->SetLineColor(419);
      h_met[i]->SetFillColor(419);
      h_mt2[i]->SetLineColor(419);
      h_mt2[i]->SetFillColor(419);
      h_deltaPhiMin[i]->SetLineColor(419);
      h_deltaPhiMin[i]->SetFillColor(419);
      h_diffMetMht[i]->SetLineColor(419);
      h_diffMetMht[i]->SetFillColor(419);
    }
    else if( fSamples[i].id >= 300 && fSamples[i].id < 500 ){
      h_njets[i]->SetLineColor(855);
      h_njets[i]->SetFillColor(855);
      h_nbjets[i]->SetLineColor(855);
      h_nbjets[i]->SetFillColor(855);
      h_nleps[i]->SetLineColor(855);
      h_nleps[i]->SetFillColor(855);
      h_ht[i]->SetLineColor(855);
      h_ht[i]->SetFillColor(855);
      h_met[i]->SetLineColor(855);
      h_met[i]->SetFillColor(855);
      h_mt2[i]->SetLineColor(855);
      h_mt2[i]->SetFillColor(855);
      h_deltaPhiMin[i]->SetLineColor(855);
      h_deltaPhiMin[i]->SetFillColor(855);
      h_diffMetMht[i]->SetLineColor(855);
      h_diffMetMht[i]->SetFillColor(855);
    }

    std::cout << std::endl << std::endl;
    std::cout << "-> Starting filling histograms for sample: " << fSamples[i].name << std::endl;

    TFile* file = TFile::Open(fSamples[i].file.c_str());
    TTree* tree = (TTree*)file->Get("mt2");

    MT2Tree myTree;
    myTree.Init(tree);
    
    int nentries;
    nentries = tree->GetEntries(); 
           
    std::cout << std::endl << std::endl;
    std::cout << "-> Starting loop over " << nentries << " entries..." << std::endl;

    for( long iEntry=0; iEntry<nentries; ++iEntry ) {

      if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

      myTree.GetEntry(iEntry);
    
      if( !myTree.passSelection() ) continue;
  
      float ht       = myTree.ht;
      float met      = myTree.met_pt;
      float mt2      = myTree.mt2;
      int njets      = myTree.nJet40;
      int nbjets     = myTree.nBJet40;
      int nleps      = myTree.nLepLowMT;
      
      double weight   = myTree.evt_scale1fb*lumi;

      if( mt2 < 200. ) continue;
      if( ht < 1000. && met < 200. ) continue;
      if( ht >= 1000. && met < 30. ) continue;

      h_deltaPhiMin[i] ->Fill(myTree.deltaPhiMin, weight);
      h_diffMetMht[i]  ->Fill(myTree.diffMetMht/myTree.met_pt, weight);
      h_njets[i]       ->Fill(njets, weight);
      h_nbjets[i]      ->Fill(nbjets, weight);
      h_nleps[i]       ->Fill(nleps, weight);
      h_ht[i]          ->Fill(ht, weight);
      h_met[i]         ->Fill(met, weight);
      h_mt2[i]         ->Fill(mt2, weight);      
      
    }// entries
    
    std::cout << "Done with sample " << fSamples[i].name << std::endl << std::endl;

    delete tree;
    
    file->Close();
    delete file;

  }// samples
  
  for( unsigned i=0; i < bSamples.size(); ++i ){
  
    h_njets_stack    ->Add(h_njets[i]);
    h_nbjets_stack   ->Add(h_nbjets[i]);
    h_nleps_stack    ->Add(h_nleps[i]);
    h_ht_stack       ->Add(h_ht[i]);
    h_met_stack      ->Add(h_met[i]);
    h_mt2_stack      ->Add(h_mt2[i]);
    h_deltaPhiMin_stack ->Add(h_deltaPhiMin[i]);
    h_diffMetMht_stack ->Add(h_diffMetMht[i]);
    
  }// BG samples

  int sigColors[]={1, 2, 6, 5, 7, 9};
  for( unsigned i=0; i < sSamples.size(); ++i ){
    
    h_njets[i+bSamples.size()]      ->SetLineColor(sigColors[i]);
    h_nbjets[i+bSamples.size()]     ->SetLineColor(sigColors[i]);
    h_nleps[i+bSamples.size()]      ->SetLineColor(sigColors[i]);
    h_ht[i+bSamples.size()]         ->SetLineColor(sigColors[i]);
    h_met[i+bSamples.size()]        ->SetLineColor(sigColors[i]);
    h_mt2[i+bSamples.size()]        ->SetLineColor(sigColors[i]);
    h_deltaPhiMin[i+bSamples.size()]->SetLineColor(sigColors[i]);
    h_diffMetMht[i+bSamples.size()] ->SetLineColor(sigColors[i]);

    h_njets[i+bSamples.size()]      ->SetLineWidth(2);
    h_nbjets[i+bSamples.size()]     ->SetLineWidth(2);
    h_nleps[i+bSamples.size()]      ->SetLineWidth(2);
    h_ht[i+bSamples.size()]         ->SetLineWidth(2);
    h_met[i+bSamples.size()]        ->SetLineWidth(2);
    h_mt2[i+bSamples.size()]        ->SetLineWidth(2);
    h_deltaPhiMin[i+bSamples.size()]->SetLineWidth(2);
    h_diffMetMht[i+bSamples.size()] ->SetLineWidth(2);    
    
    h_njets[i+bSamples.size()]      ->Scale(50);
    h_nbjets[i+bSamples.size()]     ->Scale(50);
    h_nleps[i+bSamples.size()]      ->Scale(50);
    h_ht[i+bSamples.size()]         ->Scale(50);
    h_met[i+bSamples.size()]        ->Scale(50);
    h_mt2[i+bSamples.size()]        ->Scale(50);
    h_deltaPhiMin[i+bSamples.size()]->Scale(50);
    h_diffMetMht[i+bSamples.size()] ->Scale(50);
    
  }


  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  
  //For inclusive plots:
  std::vector< std::string > bgYields;
  //bgYields.push_back("QCD");
  bgYields.push_back("Z+jets");
  bgYields.push_back("W+jets");
  bgYields.push_back("Top");
  
  std::vector< int > colors;
  //colors.push_back(401);
  colors.push_back(419);
  colors.push_back(417);
  colors.push_back(855);

  std::vector< std::string > sigYields;
  sigYields.push_back("T1tttt 1500, 100 x50");
  sigYields.push_back("T1bbbb 1500, 100 x50");
  sigYields.push_back("T1qqqq 1400, 100 x50");
  sigYields.push_back("T2tt 850, 100 x50");
  sigYields.push_back("T2bb 900, 100 x50");
  sigYields.push_back("T2qq 1200, 100 x50");

  std::vector< std::string > niceNames;
  niceNames.push_back("H_{T} > 450 GeV");
  niceNames.push_back("#geq 2j");
  niceNames.push_back("M_{T2} > 200 GeV");
  TPaveText* regionText[niceNames.size()];
  for( unsigned i=0; i<niceNames.size(); ++i ) {
    
    float yMax = 0.9-(float)i*0.04;
    float yMin = yMax - 0.04;
    regionText[i] = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
    //regionText[i]->SetTextSize(0.035);
    regionText[i]->SetTextFont(42);
    regionText[i]->SetFillColor(0);
    regionText[i]->SetTextAlign(11);
    regionText[i]->AddText( niceNames[i].c_str() );

  }


  TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+sigYields.size()+1)*0.03, 0.93, 0.9 );
  //  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  
  // Declaring empty histos for legend definition
  TH1F* hLegend[(int)bgYields.size()];
  for( unsigned b=0; b<bgYields.size(); ++b){ 
    
    hLegend[b] = new TH1F( bgYields[b].c_str(), "", 1, 0, 1 );
    hLegend[b]->SetFillColor(colors[b]);

    legend->AddEntry( hLegend[b], bgYields[b].c_str(), "F" );                                                                                                            	

  }
  for( unsigned s=0; s<sigYields.size(); ++s )
    legend->AddEntry( h_njets[s + bSamples.size()], sigYields[s].c_str(), "l" );
  

  float YMax;
  
  YMax = 1.5*(h_njets_stack->GetMaximum());
  TCanvas* c1=new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  
  h_njets_stack ->Draw("hist");
  h_njets_stack ->GetXaxis()->SetTitle("N(jets)");
  h_njets_stack ->GetYaxis()->SetTitle("Events");
  h_njets_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_njets_stack ->SetMaximum(YMax);
  h_njets_stack ->Draw("hist");

  for( unsigned s=0; s<sigYields.size(); ++s )
    h_njets[s+bSamples.size()]->Draw("same");

  legend->Draw("same");
  for( unsigned i=0; i<niceNames.size(); ++i )
    regionText[i]->Draw("same");
  labelTop->Draw("same");
  //gPad->RedrawAxis();

  c1->Update();

  c1->SaveAs(Form("%s/inclusiveNJets.pdf", outputdir.c_str()));
  
  YMax = 1.5*(h_nbjets_stack->GetMaximum());
  TCanvas* c2=new TCanvas("c2", "c2", 600, 600);
  c2->cd();
  h_nbjets_stack ->Draw("hist");
  h_nbjets_stack ->GetXaxis()->SetTitle("N(b-jets)");
  h_nbjets_stack ->GetYaxis()->SetTitle("Events");
  h_nbjets_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_nbjets_stack ->SetMaximum(YMax);
  h_nbjets_stack ->Draw("hist");
  
  for( unsigned s=0; s<sigYields.size(); ++s )
    h_nbjets[s+bSamples.size()]->Draw("same");
  
  legend->Draw("same");
  for( unsigned i=0; i<niceNames.size(); ++i )
    regionText[i]->Draw("same");
  labelTop->Draw("same");
  //gPad->RedrawAxis();
  

  c2->Update();
  c2->SaveAs(Form("%s/inclusiveNbjets.pdf", outputdir.c_str()));
  
  YMax = 1.5*(h_nleps_stack->GetMaximum());
  TCanvas* c3=new TCanvas("c3", "c3", 600, 600);
  c3->cd();
  h_nleps_stack ->Draw("hist");
  h_nleps_stack ->GetXaxis()->SetTitle("N(leptons)");
  h_nleps_stack ->GetYaxis()->SetTitle("Events");
  h_nleps_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_nleps_stack ->SetMaximum(YMax);
  h_nleps_stack ->Draw("hist");
  
  for( unsigned s=0; s<sigYields.size(); ++s )
    h_nleps[s+bSamples.size()]->Draw("same");
  
  legend->Draw("same");
  for( unsigned i=0; i<niceNames.size(); ++i )
    regionText[i]->Draw("same");
  labelTop->Draw("same");
  //gPad->RedrawAxis();
  
  c3->Update();
  c3->SaveAs(Form("%s/inclusiveNleps.pdf", outputdir.c_str()));

  YMax = 15*(h_ht_stack->GetMaximum());
  TCanvas* c4=new TCanvas("c4", "c4", 600, 600);
  c4->cd();
  gPad->SetLogy();

  h_ht_stack ->Draw("hist");
  h_ht_stack ->GetXaxis()->SetTitle("H_{T} [GeV]");
  h_ht_stack ->GetYaxis()->SetTitle("Events/25 GeV");
  h_ht_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_ht_stack ->SetMinimum(1.e-1);
  h_ht_stack ->SetMaximum(YMax);
  h_ht_stack ->GetYaxis()->SetRangeUser(1.e-1, YMax);
  h_ht_stack ->Draw("hist");
  
  for( unsigned s=0; s<sigYields.size(); ++s )
    h_ht[s+bSamples.size()]->Draw("same");

  legend->Draw("same");
  for( unsigned i=0; i<niceNames.size(); ++i )
    regionText[i]->Draw("same");
  labelTop->Draw("same");
  //gPad->RedrawAxis();


  c4->Update();
  c4->SaveAs(Form("%s/HTinclusive.pdf", outputdir.c_str()));

  YMax = 15*(h_met_stack->GetMaximum());
  TCanvas* c5=new TCanvas("c5", "c5", 600, 600);
  c5->cd();
  gPad->SetLogy();

  h_met_stack ->Draw("hist");
  h_met_stack ->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  h_met_stack ->GetYaxis()->SetTitle("Events/25 GeV");
  h_met_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_met_stack ->SetMinimum(1.e-1);
  h_met_stack ->SetMaximum(YMax);
  h_met_stack ->GetYaxis()->SetRangeUser(1.e-1, YMax);
  h_met_stack ->Draw("hist");
  
  for( unsigned s=0; s<sigYields.size(); ++s )
    h_met[s+bSamples.size()]->Draw("same");
    
  legend->Draw("same");
  for( unsigned i=0; i<niceNames.size(); ++i )
    regionText[i]->Draw("same");
  labelTop->Draw("same");
  //gPad->RedrawAxis();
  
  c5->Update();
  c5->SaveAs(Form("%s/METinclusive.pdf", outputdir.c_str()));

  YMax = 15*(h_mt2_stack->GetMaximum());
  TCanvas* c6=new TCanvas("c6", "c6", 600, 600);
  c6->cd();
  gPad->SetLogy();

  h_mt2_stack ->Draw("hist");
  h_mt2_stack ->GetXaxis()->SetTitle("M_{T2} [GeV]");
  h_mt2_stack ->GetYaxis()->SetTitle("Events/25 GeV");
  h_mt2_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_mt2_stack ->SetMinimum(1.e-1);
  h_mt2_stack ->SetMaximum(YMax);
  h_mt2_stack ->GetYaxis()->SetRangeUser(1.e-1, YMax);
  h_mt2_stack ->Draw("hist");

  for( unsigned s=0; s<sigYields.size(); ++s )
    h_mt2[s+bSamples.size()]->Draw("same");
  
  legend->Draw("same");
  for( unsigned i=0; i<niceNames.size(); ++i )
    regionText[i]->Draw("same");
  labelTop->Draw("same");
  //gPad->RedrawAxis();

  c6->Update();
  c6->SaveAs(Form("%s/MT2inclusive.pdf", outputdir.c_str()));

  YMax = 1.5*(h_deltaPhiMin_stack->GetMaximum());
  TCanvas* c7=new TCanvas("c7", "c7", 600, 600);
  c7->cd();
  //gPad->SetLogy();

  h_deltaPhiMin_stack ->Draw("hist");
  h_deltaPhiMin_stack ->GetXaxis()->SetTitle("min #Delta#Phi(4 jets, E_{T}^{miss})");
  h_deltaPhiMin_stack ->GetYaxis()->SetTitle("Events/0.1");
  h_deltaPhiMin_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_deltaPhiMin_stack ->SetMaximum(YMax);
  h_deltaPhiMin_stack ->Draw("hist");
 
  for( unsigned s=0; s<sigYields.size(); ++s )
    h_deltaPhiMin[s+bSamples.size()]->Draw("same");

  legend->Draw("same");
  for( unsigned i=0; i<niceNames.size(); ++i )
    regionText[i]->Draw("same");
  labelTop->Draw("same");
  //gPad->RedrawAxis();
  
  c7->Update();
  c7->SaveAs(Form("%s/inclusiveDeltaPhiMin.pdf", outputdir.c_str()));

  YMax = 15*(h_diffMetMht_stack->GetMaximum());
  TCanvas* c8=new TCanvas("c8", "c8", 600, 600);
  c8->cd();
  gPad->SetLogy();

  h_diffMetMht_stack ->Draw("hist");
  h_diffMetMht_stack ->GetXaxis()->SetTitle("|E_{T}^{miss}-H_{T}^{miss}|/E_{T}^{miss}");
  h_diffMetMht_stack ->GetYaxis()->SetTitle("Events");
  h_diffMetMht_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_diffMetMht_stack ->SetMinimum(0.1);
  h_diffMetMht_stack ->SetMaximum(YMax);
  h_diffMetMht_stack ->Draw("hist");

  for( unsigned s=0; s<sigYields.size(); ++s )
    h_diffMetMht[s+bSamples.size()]->Draw("same");

  legend->Draw("same");
  for( unsigned i=0; i<niceNames.size(); ++i )
    regionText[i]->Draw("same");
  labelTop->Draw("same");
  //gPad->RedrawAxis();

  c8->Update();
  c8->SaveAs(Form("%s/inclusivediffMetMht.pdf", outputdir.c_str()));

  return 0;

}// main
