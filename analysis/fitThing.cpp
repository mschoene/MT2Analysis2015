#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooPlot.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"

#include "../interface/MT2DrawTools.h"

#define mt2_cxx
#include "../interface/mt2.h"



using namespace RooFit;


RooRealVar* x = new RooRealVar( "x", "", 400., 1000. ); //global var - fuck it


void fitDataset( RooDataSet* data );
RooDataSet* getDataset( TFile* file, const std::string& name );


int main() {

  MT2DrawTools::setStyle();

  TFile* fileEE = new TFile("selEventsEE.root");
  TFile* fileMM = new TFile("selEventsMM.root");
  TFile* fileEEMM = new TFile("selEventsEEMM.root");

  RooDataSet* dataEE = getDataset(fileEE, "ee");
  RooDataSet* dataMM = getDataset(fileMM, "mm");
  RooDataSet* dataEEMM = getDataset(fileEEMM, "eemm");

  fitDataset( dataEE );
  fitDataset( dataMM );
  fitDataset( dataEEMM );


  return 0;

}



RooDataSet* getDataset( TFile* file, const std::string& name ) {



  TTree* tree = (TTree*)file->Get("mt2");

  MT2Tree myTree;
  myTree.Init(tree);


  RooDataSet* data = new RooDataSet( Form("data_%s", name.c_str()), "", RooArgSet(*x) );
  

  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( myTree.nlep<2 ) continue;

    TLorentzVector lep0;
    lep0.SetPtEtaPhiM( myTree.lep_pt[0], myTree.lep_eta[0], myTree.lep_phi[0], myTree.lep_mass[0] );
    TLorentzVector lep1;
    lep1.SetPtEtaPhiM( myTree.lep_pt[1], myTree.lep_eta[1], myTree.lep_phi[1], myTree.lep_mass[1] );

    if( myTree.lep_pt[0]<25. ) continue;
    if( myTree.lep_pt[1]<25. ) continue;
    if( abs(myTree.lep_pdgId[0])==11 ) 
      if( !(myTree.lep_tightId[0]>0) ) continue;
    if( abs(myTree.lep_pdgId[1])==11 ) 
      if( !(myTree.lep_tightId[1]>0) ) continue;

    TLorentzVector dilep = lep0 + lep1;
    float mll = dilep.M();

    std::vector<TLorentzVector> jets;
    std::vector<TLorentzVector> bjets;
    for( int i=0; i<myTree.njet; ++i ) {
      TLorentzVector thisJet; 
      thisJet.SetPtEtaPhiM( myTree.jet_pt[i], myTree.jet_eta[i], myTree.jet_phi[i], myTree.jet_mass[i] );
      if( thisJet.Pt()<30. ) continue;
      //if( fabs(thisJet.Eta())>2.5 ) continue;
      float deltaR0 = thisJet.DeltaR( lep0 );
      if( deltaR0<0.4 ) continue;
      float deltaR1 = thisJet.DeltaR( lep1 );
      if( deltaR1<0.4 ) continue;
      jets.push_back(thisJet);
      if( myTree.jet_btagCSV[i]>0.89 )
        bjets.push_back(thisJet);
    }

    int nJets = jets.size();
    int nBJets = bjets.size();


    if( nJets!=2 ) continue;
    if( nBJets!=0 ) continue;
    if( myTree.met_pt>60. ) continue;
    
    x->setVal(mll);
    data->add( RooArgList(*x) ); 

  }


  return data;

}



void fitDataset( RooDataSet* data ) {


  std::string name(data->GetName());


  RooRealVar sigmean ("sigmean","Mthing", 760., 700., 800.);
  RooRealVar sigwidth("sigwidth","Wthing", 20., 10., 100.); 
  RooGaussian gauss("gauss","gaussian PDF", *x, sigmean, sigwidth);

  RooRealVar expoSlope("expoSlope","exponential slope parameter", -0.01, -0.1, 0.) ; 
  RooExponential expo("expoBG","",*x,expoSlope);

  // --- Construct signal+background PDF ---
  RooRealVar nsig("nsig","#signal events", 15, 0.,100);
  RooRealVar nbkg("nbkg","#background events", 200, 0., 10000) ;
  RooAddPdf sum("sum","g+a",RooArgList(gauss,expo),RooArgList(nsig,nbkg)) ;

  // --- Perform extended ML fit of composite PDF to toy data ---
  sum.fitTo(*data,Extended()) ;


  RooPlot* xframe = x->frame();
  data->plotOn(xframe, Binning(60, 400, 1000));
  sum.plotOn(xframe);

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  gPad->SetLeftMargin(0.15);
  TH2D* h2_axes = new TH2D("axes", "", 10, 400., 1000, 10, 0., xframe->GetMaximum()*1.1 );
  h2_axes->SetXTitle("Dilepton Invariant Mass [GeV]");
  h2_axes->SetYTitle("Events");
  h2_axes->Draw();
  xframe->GetYaxis()->SetTitleOffset(1.4); 
  xframe->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop(2.1);
  labelTop->Draw("same");

  c1->SaveAs(Form("thingFit_%s.eps", data->GetName()));
  c1->SaveAs(Form("thingFit_%s.pdf", data->GetName()));

  delete c1;
  delete xframe;


}
