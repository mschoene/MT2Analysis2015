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
#include "TRandom3.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"

#include "../interface/MT2DrawTools.h"

#define mt2_cxx
#include "../interface/mt2.h"



using namespace RooFit;


RooRealVar* x = new RooRealVar( "x", "", 400., 1000. ); //global var - fuck it


RooDataSet* getDataset( TTree* tree, const std::string& name, int nJetMin, int nJetMax );
void fitDataset( const std::string& outdir, ofstream& ofs, RooDataSet* data );



int main( int argc, char* argv[] ) {


  int nJetMin=2;
  int nJetMax=100;
  if( argc>1 ) {
    std::string nJetMin_str(argv[1]);
    nJetMin = atoi(nJetMin_str.c_str());
  }
  if( argc>2 ) {
    std::string nJetMax_str(argv[2]);
    nJetMax = atoi(nJetMax_str.c_str());
  }



  MT2DrawTools::setStyle();

  TFile* fileEE = new TFile("selEventsEE.root");
  TFile* fileMM = new TFile("selEventsMM.root");
  TFile* fileEEMM = new TFile("selEventsEEMM.root");

  TTree* treeEE   = (TTree*)fileEE  ->Get("mt2");
  TTree* treeMM   = (TTree*)fileMM  ->Get("mt2");
  TTree* treeEEMM = (TTree*)fileEEMM->Get("mt2");

  RooDataSet* dataEE  = getDataset(treeEE  , "ee"  , nJetMin, nJetMax);
  RooDataSet* dataMM  = getDataset(treeMM  , "mm"  , nJetMin, nJetMax);
  RooDataSet* dataEEMM= getDataset(treeEEMM, "eemm", nJetMin, nJetMax);


  std::string outdir;
  if( nJetMax>=0 ) 
    outdir = std::string(Form("thing_nJ%d%d", nJetMin, nJetMax));
  else
    outdir = std::string(Form("thing_nJ%dInf", nJetMin));
  system( Form("mkdir -p %s", outdir.c_str()) );

  std::ofstream ofs(Form("%s/log.log", outdir.c_str()) );

  fitDataset( outdir, ofs, dataEE  );
  fitDataset( outdir, ofs, dataMM  );
  fitDataset( outdir, ofs, dataEEMM);

  //fitDataset( dataEE_j2  );
  //fitDataset( dataMM_j2  );
  //fitDataset( dataEEMM_j2);


  return 0;

}



RooDataSet* getDataset( TTree* tree, const std::string& name, int nJetMin, int nJetMax ) {


  //TTree* tree = (TTree*)file->Get("mt2");

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

    //if( myTree.lep_pdgId[0]!=myTree.lep_pdgId[1] ) continue;

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

    bool jetsOK = nJets>=nJetMin;
    if( nJetMax>=nJetMin ) jetsOK = jetsOK && nJets<=nJetMax;

    if( !jetsOK ) continue;
    if( nBJets!=0 ) continue;
    if( myTree.met_pt>60. ) continue;
    
    x->setVal(mll);
    data->add( RooArgList(*x) ); 

  }


  return data;

}



void fitDataset( const std::string& outdir, ofstream& ofs, RooDataSet* data ) {


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

  c1->SaveAs(Form("%s/fit_%s.eps", outdir.c_str(), data->GetName()));
  c1->SaveAs(Form("%s/fit_%s.pdf", outdir.c_str(), data->GetName()));



  float xMin_int = sigmean.getVal()-sigwidth.getVal();
  float xMax_int = sigmean.getVal()+sigwidth.getVal();
  //float xMin_int = 750.;
  //float xMax_int = 790.;
  x->setRange("sigregion", xMin_int, xMax_int);
  ofs << "  range  : " << xMin_int << "-" << xMax_int << " GeV" << std::endl;
  
  Double_t obs    = data->sumEntries(Form("x>%f && x<%f", xMin_int, xMax_int));
  RooAbsReal* intSig = gauss.createIntegral(*x,NormSet(*x),Range("sigregion")) ;
  RooAbsReal* intBkg = expo .createIntegral(*x,NormSet(*x),Range("sigregion")) ;
  ofs << "  signal  = " << intSig->getVal()*nsig.getVal() << std::endl; // " (+" << intSig->getErrorHi()*nsig.getVal() << ")(-" << intSig->getErrorLo()*nsig.getVal() << std::endl ;

  double bgMean = intBkg->getVal()*nbkg.getVal();
  double bgErr = intBkg->getVal()*nbkg.getError();
  ofs << "  bkg     = " << bgMean << " +- " << bgErr << std::endl;
  ofs << "  obs     = " << obs << std::endl;

  ofs << "  sigmean = " << sigmean.getVal() << " +- " << sigmean.getError()  << std::endl;
  ofs << "  sigwidth = " << sigwidth.getVal() << " +- " << sigwidth.getError()  << std::endl;

  
  TRandom3 randGen;//(int(obs));

  unsigned int counterN(0);
  unsigned int counterD(0);


  bool doLeft = obs<=bgMean;

  for(int i=0; i<1e6; i++){
    double meanShift = randGen.Gaus(0.,bgErr);
    //cout << "meanShift: " << meanShift << endl;
    double rand = randGen.Poisson(bgMean+meanShift); counterD++;

    if(doLeft) 
      {if(rand<=obs) counterN++;}
    else{
      if(rand>=obs) counterN++;
    }

    //cout << "rand " << i << " : " << rand.Poisson(bgMean) << endl;

  }
  
  double prob=1.0*counterN/counterD;
  double significance  = TMath::NormQuantile(1-prob);

  ofs << "  probability: " << prob  << std::endl;
  ofs << "  significance: " << significance << std::endl;


  delete c1;
  delete xframe;
  delete h2_axes;


}
