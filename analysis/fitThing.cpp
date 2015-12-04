#include <iostream>
#include <fstream>

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

  RooDataSet* dataEE_j23  = getDataset(fileEE  , "ee_j23"  , 2, 3);
  RooDataSet* dataMM_j23  = getDataset(fileMM  , "mm_j23"  , 2, 3);
  RooDataSet* dataEEMM_j23= getDataset(fileEEMM, "eemm_j23", 2, 3);

  RooDataSet* dataEE_j2   = getDataset(fileEE  , "ee_j2"  , 2, 2);
  RooDataSet* dataMM_j2   = getDataset(fileMM  , "mm_j2"  , 2, 2);
  RooDataSet* dataEEMM_j2 = getDataset(fileEEMM, "eemm_j2", 2, 2);

  ofstream ofs("logThing.log");

  fitDataset( ofs, dataEE_j23  );
  fitDataset( ofs, dataMM_j23  );
  fitDataset( ofs, dataEEMM_j23);

  fitDataset( ofs, dataEE_j2  );
  fitDataset( ofs, dataMM_j2  );
  fitDataset( ofs, dataEEMM_j2);

  ofs.close();

  return 0;

}



RooDataSet* getDataset( TFile* file, const std::string& name, int nJetMin, int nJetMax ) {



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

    bool jetsOK = nJets>=nJetMin && nJets<=nJetMax;

    if( !jetsOK ) continue;
    if( nBJets!=0 ) continue;
    if( myTree.met_pt>60. ) continue;
    
    x->setVal(mll);
    data->add( RooArgList(*x) ); 

  }


  return data;

}



void fitDataset( ofstream ofs, RooDataSet* data ) {


  ofs << "fit: " << data->GetName() << std::endl;
  std::string name(data->GetName());


  RooRealVar sigmean ("sigmean","Mthing", 770., 750., 800.);
  RooRealVar sigwidth("sigwidth","Wthing", 15., 5., 40.); 
  RooGaussian gauss("gauss","gaussian PDF", *x, sigmean, sigwidth);

  RooRealVar expoSlope("expoSlope","exponential slope parameter", -0.01, -0.1, 0.) ; 
  RooExponential expo("expoBG","",*x,expoSlope);

  // --- Construct signal+background PDF ---
  RooRealVar nsig("nsig","#signal events", 15, 0.,100);
  RooRealVar nbkg("nbkg","#background events", 200, 0., 10000) ;
  RooAddPdf sum("sum","g+a",RooArgList(gauss,expo),RooArgList(nsig,nbkg)) ;

  // --- Perform extended ML fit of composite PDF to toy data ---
  sum.fitTo(*data,Extended()) ;


  int nBins = 60;
  float xMin = 400.;
  float xMax = 1000.;
  float binWidth = (xMax-xMin)/((float)nBins);

  RooPlot* xframe = x->frame();
  data->plotOn(xframe, Binning(nBins, xMin, xMax));
  sum.plotOn(xframe);
  sum.plotOn(xframe,Components(expo),LineStyle(kDashed)) ;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  gPad->SetLeftMargin(0.15);
  TH2D* h2_axes = new TH2D("axes", "", 10, 400., 1000, 10, 0., xframe->GetMaximum() );
  h2_axes->SetXTitle("Dilepton Invariant Mass [GeV]");
  h2_axes->SetYTitle(Form("Events / %.0f GeV", binWidth));
  h2_axes->Draw();
  xframe->GetYaxis()->SetTitleOffset(1.4); 
  xframe->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop(2.1);
  labelTop->Draw("same");

  c1->SaveAs(Form("thingFit_%s.eps", data->GetName()));
  c1->SaveAs(Form("thingFit_%s.pdf", data->GetName()));



  float xMin_int = 750.;
  float xMax_int = 790.;
  x->setRange("sigregion", xMin_int, xMax_int);
  //x->setRange("sigregion", sigmean.getVal()-2.*sigwidth.getVal(), sigmean.getVal()+2.*sigwidth.getVal());
  
  RooAbsReal* intSig = gauss.createIntegral(*x,NormSet(*x),Range("sigregion")) ;
  RooAbsReal* intBkg = expo .createIntegral(*x,NormSet(*x),Range("sigregion")) ;
  ofs << "  signal  = " << intSig->getVal()*nsig.getVal() << std::endl; // " (+" << intSig->getErrorHi()*nsig.getVal() << ")(-" << intSig->getErrorLo()*nsig.getVal() << std::endl ;
  ofs << "  bkg     = " << intBkg->getVal()*nbkg.getVal() << std::endl; // " (+" << intBkg->getErrorHi()*nbkg.getVal() << ")(-" << intBkg->getErrorLo()*nbkg.getVal() << std::endl ;
  ofs << "  sigmean = " << sigmean->getVal() << " +- " << sigmean.getError()  << std::endl
  ofs << "  sigwidth = " << sigwidth->getVal() << " +- " << sigwidth.getError()  << std::endl



  delete c1;
  delete xframe;
  delete h2_axes;


}
