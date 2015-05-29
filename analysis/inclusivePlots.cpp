#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>


#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
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


  if( argc!=2 ) {
    std::cout << "USAGE: ./inclusivePlots [samplesFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  
  bool isReduceTree = false;

  std::string sampleName(argv[1]);
  
  std::string samplesFileName = "../samples/samples_" + sampleName + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName);
  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }
  
  std::string outputdir = "inclusivePlots_"+sampleName;
  system(Form("mkdir -p %s", outputdir.c_str()));


  // DECLARING STACK HISTOGRAM - ONE PER VARIABLE
  THStack* h_njets_stack  = new THStack("h_njets_stack", "");
  THStack* h_nbjets_stack = new THStack("h_nbjets_stack", "");
  THStack* h_nleps_stack  = new THStack("h_nleps_stack", "");
  THStack* h_ht_stack     = new THStack("h_ht_stack", "");
  THStack* h_met_stack    = new THStack("h_met_stack", "");
  THStack* h_mt2_stack    = new THStack("h_mt2_stack", "");
  THStack* h_deltaPhiMin_stack  = new THStack("h_deltaPhiMin_stack", "");
  THStack* h_MetMinusMht_stack  = new THStack("h_MetMinusMht_stack", "");


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
  TH1F* h_MetMinusMht[nSamples+1];

  int S=-1;
  
  for( unsigned i=0; i<fSamples.size(); ++i ){
    
    //Initializing histograms:    
    std::string name_njets="njets"+fSamples[i].name;
    h_njets[i] = new TH1F(name_njets.c_str(), "number of jets", 12, 0, 12 );
    h_njets[i]->GetXaxis()->SetTitle("N(jets)");
    h_njets[i]->GetYaxis()->SetTitle("Events");
  
    std::string name_nbjets="nbjets"+fSamples[i].name;
    h_nbjets[i] = new TH1F(name_nbjets.c_str(), "number of b-jets", 6, 0, 6 );
    h_nbjets[i]->GetXaxis()->SetTitle("N(b-jets)");
    h_nbjets[i]->GetYaxis()->SetTitle("Events");
  
    std::string name_nleps="nleps"+fSamples[i].name;
    h_nleps[i] = new TH1F(name_nleps.c_str(), "number of leptons", 10, 0, 10 );
    h_nleps[i]->GetXaxis()->SetTitle("N(leptons)");
    h_nleps[i]->GetYaxis()->SetTitle("Events");
    
    std::string name_ht="ht"+fSamples[i].name;
    h_ht[i] = new TH1F(name_ht.c_str(), "H_{T}", 300, 0, 3000 );
    h_ht[i]->GetXaxis()->SetTitle("H_{T} [GeV]");
    h_ht[i]->GetYaxis()->SetTitle("Events/10 GeV");
  
    std::string name_met="met"+fSamples[i].name;
    h_met[i] = new TH1F(name_met.c_str(), "E_{T}^{miss}", 100, 0, 1000 );
    h_met[i]->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    h_met[i]->GetYaxis()->SetTitle("Events/10 GeV");
   
    std::string name_mt2="mt2"+fSamples[i].name;;
    h_mt2[i] = new TH1F(name_mt2.c_str(), "M_{T2}", 100, 0, 1000 );
    h_mt2[i]->GetXaxis()->SetTitle("M_{T2} [GeV]");
    h_mt2[i]->GetYaxis()->SetTitle("Events/10 GeV");
    
    std::string name_deltaPhiMin="mindPhi"+fSamples[i].name;;
    h_deltaPhiMin[i] = new TH1F(name_deltaPhiMin.c_str(), "min #delta #Phi (4 jets, E_{T}^{miss})", 32, 0, 3.2 );
    h_deltaPhiMin[i]->GetXaxis()->SetTitle("min #delta #Phi (4 jets, E_{T}^{miss})");
    h_deltaPhiMin[i]->GetYaxis()->SetTitle("Events");
    
    std::string name_MetMinusMht="MetMinusMht"+fSamples[i].name;;
    h_MetMinusMht[i] = new TH1F(name_MetMinusMht.c_str(), "(E_{T}^{miss} - H_{T}^{miss})/E_{T}^{miss}", 20, 0, 1 );
    h_MetMinusMht[i]->GetXaxis()->SetTitle("(E_{T}^{miss} - H_{T}^{miss})/E_{T}^{miss}");
    h_MetMinusMht[i]->GetYaxis()->SetTitle("Events");


    if( fSamples[i].name.find("QCD") != std::string::npos ){
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
      h_MetMinusMht[i]->SetLineColor(401);
      h_MetMinusMht[i]->SetFillColor(401);
    }
    else if( fSamples[i].name.find("WJets") != std::string::npos ){
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
      h_MetMinusMht[i]->SetLineColor(417);
      h_MetMinusMht[i]->SetFillColor(417);
    }
    else if(fSamples[i].name.find("ZJets") != std::string::npos ){
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
      h_MetMinusMht[i]->SetLineColor(419);
      h_MetMinusMht[i]->SetFillColor(419);
    }
    else if( fSamples[i].name.find("TTJets") != std::string::npos ){
      h_njets[i]->SetLineColor(4);
      h_njets[i]->SetFillColor(4);
      h_nbjets[i]->SetLineColor(4);
      h_nbjets[i]->SetFillColor(4);
      h_nleps[i]->SetLineColor(4);
      h_nleps[i]->SetFillColor(4);
      h_ht[i]->SetLineColor(4);
      h_ht[i]->SetFillColor(4);
      h_met[i]->SetLineColor(4);
      h_met[i]->SetFillColor(4);
      h_mt2[i]->SetLineColor(4);
      h_mt2[i]->SetFillColor(4);
      h_deltaPhiMin[i]->SetLineColor(4);
      h_deltaPhiMin[i]->SetFillColor(4);
      h_MetMinusMht[i]->SetLineColor(4);
      h_MetMinusMht[i]->SetFillColor(4);
    }
    else if( fSamples[i].name.find("T1") != std::string::npos ){
      h_njets[i]->SetLineColor(1);
      h_nbjets[i]->SetLineColor(1);
      h_nleps[i]->SetLineColor(1);
      h_ht[i]->SetLineColor(1);
      h_met[i]->SetLineColor(1);
      h_mt2[i]->SetLineColor(1);
      h_deltaPhiMin[i]->SetLineColor(1);
      h_MetMinusMht[i]->SetLineColor(1);
      S=i;
    }

    std::cout << std::endl << std::endl;
    std::cout << "-> Starting filling histograms for sample: " << fSamples[i].name << std::endl;

    TFile* file = TFile::Open(fSamples[i].file.c_str());
    TTree* tree = (TTree*)file->Get("mt2");

    TFile* tmpFile = TFile::Open("tmp.root", "recreate");
    TTree* tree_reduced;

    if(isReduceTree){
    
      // Define selection if reducing the tree before looping over entries:

      std::ostringstream preselectionStream;
      preselectionStream << " "
			 << "(nTaus20==0 && nMuons10==0 && nElectrons10==0)"                   << " && "
			 << "(nVert > 0)"                      << " && "
			 << "(nJet40 > 1)"                     << " && "
			 << "(ht > 450)"                       << " && "
			 << "(mt2 > 50)";
        //		       << "(jet_pt[1] > 100)"                << " && "
        //		       << "(deltaPhiMin > 0.3)"              << " && "
        //		       << "(diffMetMht < 70)";
      
      TString preselection = preselectionStream.str().c_str();
      TString cuts = preselection;
      
      /*
	TFile* tmpFile = TFile::Open("tmp.root", "recreate");
	tmpFile->cd();
	TTree* tree_reduced = tree->CopyTree(cuts);
      */

      tmpFile->cd();
      tree_reduced = tree->CopyTree(cuts);
      
    }// reduceTree

    MT2Tree myTree;
    if(isReduceTree)
      myTree.Init(tree_reduced);
    else
      myTree.Init(tree);
    
    int nentries;

    
    if(isReduceTree){      

      nentries = tree_reduced->GetEntries();
         
    }
    else {

      nentries = tree->GetEntries(); 
       
    }
    
    std::cout << std::endl << std::endl;
    std::cout << "-> Starting loop over " << nentries << " entries..." << std::endl;

    for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

      if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

      myTree.GetEntry(iEntry);
      
      if(!isReduceTree){
	
	//Define event selection if NOT reducing the input tree:
	if( myTree.nMuons10 > 0) continue;
	if( myTree.nElectrons10 > 0 ) continue;
	if( myTree.nPFLep5LowMT > 0) continue;
	if( myTree.nPFHad10LowMT > 0) continue;
	if( myTree.nVert==0 ) continue;
	if( myTree.nJet40<2 ) continue;
	if( myTree.njet<2 ) continue;
	if( myTree.mt2<200 ) continue; 
	if( myTree.ht<450 ) continue;
	if( myTree.ht<1000 && myTree.met_pt<200 ) continue;
	if( myTree.ht>1000 && myTree.met_pt<30 ) continue;
	if( myTree.jet1_pt < 40. &&  myTree.jet2_pt < 40.) continue;

//	float jetCentral_pt[2];
//	int njetsCentral = 0;
//	for(int j=0; j<myTree.njet; ++j){
//	  if( fabs( myTree.jet_eta[j] ) < 2.5 ) {
//	    jetCentral_pt[njetsCentral] = myTree.jet_pt[j];
//	    ++njetsCentral;
//	  }
//	  if( njetsCentral >= 2 ) break;
//	}
//	if (jetCentral_pt[1] < 100. ) continue;
//
      }
         
      float ht       = myTree.ht;
      float met      = myTree.met_pt;
      float mt2      = myTree.mt2;
      int njets      = myTree.nJet40;
      int nbjets     = myTree.nBJet40;
      int nmuons     = myTree.nMuons10;
      int nelectrons = myTree.nElectrons10;
      int ntaus      = myTree.nTaus20;
      int nleps   = nmuons+nelectrons+ntaus;
      int njet_all   = myTree.njet;
      
      double weight   = myTree.evt_scale1fb*lumi;


      //if( myTree.diffMetMht < 0.5*myTree.met_pt )
      h_deltaPhiMin[i] ->Fill(myTree.deltaPhiMin);

      if( myTree.deltaPhiMin<0.3 ) continue;

      h_MetMinusMht[i] ->Fill(myTree.diffMetMht/myTree.met_pt);
      
      if( myTree.diffMetMht>0.5*myTree.met_pt ) continue;
      
      h_njets[i]    ->Fill(njets, weight);
      h_nbjets[i]   ->Fill(nbjets, weight);
      h_nleps[i]    ->Fill(nleps, weight);
      h_ht[i]       ->Fill(ht, weight);
      h_met[i]      ->Fill(met, weight);
      h_mt2[i]      ->Fill(mt2, weight);      
      
    }// entries
    
    std::cout << "Done with sample " << fSamples[i].name << std::endl << std::endl; 
 
    if(i!=S){
      h_njets_stack    ->Add(h_njets[i]);
      h_nbjets_stack   ->Add(h_nbjets[i]);
      h_nleps_stack    ->Add(h_nleps[i]);
      h_ht_stack       ->Add(h_ht[i]);
      h_met_stack      ->Add(h_met[i]);
      h_mt2_stack      ->Add(h_mt2[i]);
      h_deltaPhiMin_stack ->Add(h_deltaPhiMin[i]);
      h_MetMinusMht_stack ->Add(h_MetMinusMht[i]);
    }

    if(isReduceTree)
      delete tree_reduced;
    delete tree;
    
    tmpFile->Close();
    delete tmpFile;
       
    file->Close();
    delete file;

  }// samples

  system( "rm tmp.root" );

  // Declaring empty histos for legend definition
  TH1F* hQCD = new TH1F("QCD", "QCD", 1, 0, 1);
  hQCD->SetFillColor(401);
  
  TH1F* hTop = new TH1F("Top", "Top", 1, 0, 1);
  hTop->SetFillColor(4);

  TH1F* hWJets = new TH1F("WJets", "WJets", 1, 0, 1);
  hWJets->SetFillColor(417);

  TH1F* hZJets = new TH1F("ZJets", "ZJets", 1, 0, 1);
  hZJets->SetFillColor(419);

  TH1F* hT1 = new TH1F("T1", "T1", 1, 0, 1);
  hT1->SetLineColor(1);


  TLegend* Leg = new TLegend(.65,.65,.85,.85, "");
  Leg->SetFillColor(0);
  Leg->AddEntry(hQCD, "QCD", "F");
  //Leg->AddEntry(hTop, "Top", "F");
  //Leg->AddEntry(hWJets, "W+jets", "F");
  //Leg->AddEntry(hZJets, "Z+jets", "F");
  //Leg->AddEntry(hT1, "T1bbbb (1500, 100)", "l");

  TCanvas* c1=new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  h_njets_stack ->Draw("hist");
  h_njets_stack ->GetXaxis()->SetTitle("N(jets)");
  h_njets_stack ->GetYaxis()->SetTitle("Events");
  h_njets_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_njets_stack ->Draw("hist");
  //  h_njets[S]    ->Draw("same");
  Leg->Draw("same");
  
  c1->Update();

  c1->SaveAs(Form("%s/inclusiveNJets.eps", outputdir.c_str()));
  
  TCanvas* c2=new TCanvas("c2", "c2", 600, 600);
  c2->cd();
  h_nbjets_stack ->Draw("hist");
  h_nbjets_stack ->GetXaxis()->SetTitle("N(b-jets)");
  h_nbjets_stack ->GetYaxis()->SetTitle("Events");
  h_nbjets_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_nbjets_stack ->Draw("hist");
  //  h_nbjets[S]    ->Draw("same");

  Leg->Draw("same");
  
  c2->Update();
  c2->SaveAs(Form("%s/inclusiveNbjets.eps", outputdir.c_str()));

  TCanvas* c3=new TCanvas("c3", "c3", 600, 600);
  c3->cd();
  h_nleps_stack ->Draw("hist");
  h_nleps_stack ->GetXaxis()->SetTitle("N(leptons)");
  h_nleps_stack ->GetYaxis()->SetTitle("Events");
  h_nleps_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_nleps_stack ->Draw("hist");
  //  h_nleps[S]    ->Draw("hist");

  Leg->Draw("same");

  c3->Update();
  c3->SaveAs(Form("%s/inclusiveNleps.eps", outputdir.c_str()));

  TCanvas* c4=new TCanvas("c4", "c4", 600, 600);
  c4->cd();
  gPad->SetLogy();

  h_ht_stack ->Draw("hist");
  h_ht_stack ->GetXaxis()->SetTitle("H_{T} [GeV]");
  h_ht_stack ->GetYaxis()->SetTitle("Events/10 GeV");
  h_ht_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_ht_stack ->SetMinimum(0.1);
  h_ht_stack ->Draw("hist");
  h_ht_stack ->Draw("same");
  //  h_ht[S]    ->Draw("same");

  Leg->Draw("same");

  c4->Update();
  c4->SaveAs(Form("%s/HTinclusive.eps", outputdir.c_str()));

  TCanvas* c5=new TCanvas("c5", "c5", 600, 600);
  c5->cd();
  gPad->SetLogy();

  h_met_stack ->Draw("hist");
  h_met_stack ->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  h_met_stack ->GetYaxis()->SetTitle("Events/10 GeV");
  h_met_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_met_stack ->SetMinimum(0.1);
  h_met_stack ->Draw("hist");
  //  h_met[S]    ->Draw("same");

  Leg->Draw("same");
  
  c5->Update();
  c5->SaveAs(Form("%s/METinclusive.eps", outputdir.c_str()));

  TCanvas* c6=new TCanvas("c6", "c6", 600, 600);
  c6->cd();
  gPad->SetLogy();

  h_mt2_stack ->Draw("hist");
  h_mt2_stack ->GetXaxis()->SetTitle("M_{T2} [GeV]");
  h_mt2_stack ->GetYaxis()->SetTitle("Events/10 GeV");
  h_mt2_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_mt2_stack ->SetMinimum(0.1);
  h_mt2_stack ->Draw("hist");
  //  h_mt2[S]    ->Draw("same");

  Leg->Draw("same");

  c6->Update();
  c6->SaveAs(Form("%s/MT2inclusive.eps", outputdir.c_str()));

  TCanvas* c7=new TCanvas("c7", "c7", 600, 600);
  c7->cd();
  //gPad->SetLogy();

  h_deltaPhiMin_stack ->Draw("hist");
  h_deltaPhiMin_stack ->GetXaxis()->SetTitle("min #Delta#Phi(4 jets, E_{T}^{miss})");
  h_deltaPhiMin_stack ->GetYaxis()->SetTitle("Events/0.1");
  h_deltaPhiMin_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_deltaPhiMin_stack ->SetMinimum(0.1);
  h_deltaPhiMin_stack ->Draw("hist");
  //  h_deltaPhiMin[S]    ->Draw("same");

  Leg->Draw("same");

  c7->Update();
  c7->SaveAs(Form("%s/inclusiveDeltaPhiMin.eps", outputdir.c_str()));

  TCanvas* c8=new TCanvas("c8", "c8", 600, 600);
  c8->cd();
  gPad->SetLogy();

  h_MetMinusMht_stack ->Draw("hist");
  h_MetMinusMht_stack ->GetXaxis()->SetTitle("|E_{T}^{miss}-H_{T}^{miss}|/E_{T}^{miss}");
  h_MetMinusMht_stack ->GetYaxis()->SetTitle("Events");
  h_MetMinusMht_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_MetMinusMht_stack ->SetMinimum(0.1);
  h_MetMinusMht_stack ->Draw("hist");
  //  h_MetMinusMht[S]    ->Draw("same");

  Leg->Draw("same");

  c8->Update();
  c8->SaveAs(Form("%s/inclusiveMetMinusMht.eps", outputdir.c_str()));

  

  return 0;

}// main
