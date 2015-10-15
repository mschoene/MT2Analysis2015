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
#include "interface/MT2EstimateQCD.h"
#include "interface/MT2DrawTools.h"

#define mt2_cxx
#include "../interface/mt2.h"

#include "TEventList.h"

// this analysis:
// - gets the ratio function from subtracted data
// - uses it to estimate QCD
// - draws MT2 with QCDfromData and QCDfromMC
// note: need to rerun over samples to have an event by event estimation
// todo: extend MT2EstimateQCD to maybe inherit from MT2EstimateTree so no need to rerun over samples

float lumi = 0.454; //fb-1


// todo: make optional arguments for these guys
bool doDiffMetMht = true;
bool doHFjetVeto = false;

bool doRand = false;

bool recompute     = true; // don't loop over samples if already done and only drawing modifications are desired
bool recomputeData = true; // force recompute data. For instance if only data or fit to data has changed. 
                            // only makes a difference if recompute==false
                            // NOTE: to use with caution. if fit changed then subtraction from non-QCD will be different so things must be recomputed

MT2Analysis<MT2Estimate>* merge ( std::vector<MT2Analysis<MT2Estimate> *> anas, const std::string& regionsSet, const std::string& name, int id_min, int id_max );
void computeYield( const MT2Sample& sample, const std::string& regionsSet,
		   MT2Analysis<MT2Estimate>*  analysis,
		   MT2Analysis<MT2Estimate>*  analysis_QCDfromData, 
		   MT2Analysis<MT2EstimateQCD>* dataQCD, MT2Analysis<MT2EstimateQCD>* dataQCDpsLow, MT2Analysis<MT2EstimateQCD>* dataQCDpsMed, 
		   float lumi=1. );

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::vector<MT2Analysis<MT2Estimate>* > bgYields , const std::string& prefix);


int main( int argc, char* argv[] ) {


  if( argc>4 ) {
    std::cout << "USAGE: ./qcdControlRegion [samplesFileName] [regionsSet] [postfix]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  //std::string samplesFileName = "74X_jecV4_MET30_QCD";
  std::string samplesFileName = "Spring15_25ns_qcdSkim";
  //std::string regionsSet = "zurich_HTtriggers2";
  //std::string regionsSet = "zurich_HTtriggers";
  std::string regionsSet = "zurich_onlyHT";
  std::string postfix = "";

  if( argc>1 ) {
    std::string samplesFileName_tmp(argv[1]); 
    samplesFileName = samplesFileName_tmp;
    if( argc>2 ) {
      std::string regionsSet_tmp(argv[2]); 
      regionsSet = regionsSet_tmp; 
      if( argc>3 ) {
	std::string postfix_tmp(argv[3]); 
	postfix = postfix_tmp; 
      }
    }
  }


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string dir( Form("QCDcontrolRegion_%s_%s_%.3ffb%s%s%s", samplesFileName.c_str(), regionsSet.c_str(), lumi,(doDiffMetMht ? "" : "_noDiffMetMht"), (!doHFjetVeto ? "" : "_HFjetVeto"), postfix.c_str()) );


  // get QCD fits to subtracted data
  std::string rand_text = doRand ? "_rand" : "";
  MT2Analysis<MT2EstimateQCD>* dataQCD      = MT2Analysis<MT2EstimateQCD>::readFromFile(dir + "/qcdCRmod"+rand_text+".root", "dataSubCR");
  MT2Analysis<MT2EstimateQCD>* dataQCDps180 = MT2Analysis<MT2EstimateQCD>::readFromFile(dir + "/qcdCRmod"+rand_text+".root", "dataSubCRps180");
  MT2Analysis<MT2EstimateQCD>* dataQCDps60  = MT2Analysis<MT2EstimateQCD>::readFromFile(dir + "/qcdCRmod"+rand_text+".root", "dataSubCRps60");

  // get samples
  std::string samplesFile = "../samples/samples_" + samplesFileName + ".dat";
  std::vector<MT2Sample> samples = MT2Sample::loadSamples(samplesFile, 1, 999); // not interested in signal here (see later)
  if( samples.size()==0 ) {
    std::cout << "There must be an error: didn't find any qcd files in " << samplesFile << "!" << std::endl;
    exit(1209);
  }

  MT2Analysis<MT2Estimate>* allBkg = new MT2Analysis<MT2Estimate>("allBkg", regionsSet);
  std::vector< MT2Analysis<MT2Estimate>* > bgYields;
  std::vector< MT2Analysis<MT2Estimate>* > bgYields_QCDfromData;
  std::vector< MT2Analysis<MT2Estimate>* > EventYield;
  std::vector< MT2Analysis<MT2Estimate>* > EventYield_QCDfromData;

  MT2Analysis<MT2Estimate>* data;   
  MT2Analysis<MT2Estimate>* qcd;   
  MT2Analysis<MT2Estimate>* top;
  MT2Analysis<MT2Estimate>* wjets;
  MT2Analysis<MT2Estimate>* zjets;
  MT2Analysis<MT2Estimate>* qcd_QCDfromData;

  if(recompute){
     for( unsigned i=0; i<samples.size(); ++i ) {
       MT2Analysis<MT2Estimate>* analysis = new MT2Analysis<MT2Estimate>( samples[i].sname, regionsSet, samples[i].id );
       MT2Analysis<MT2Estimate>* analysis_QCDfromData = new MT2Analysis<MT2Estimate>( samples[i].sname, regionsSet, samples[i].id );
       computeYield( samples[i], regionsSet, analysis, analysis_QCDfromData, dataQCD, dataQCDps180, dataQCDps60, lumi);
       EventYield.push_back(analysis);
       EventYield_QCDfromData.push_back(analysis_QCDfromData);
     }
       
     data  = merge( EventYield, regionsSet, "data" ,   1,  99 );
     qcd   = merge( EventYield, regionsSet, "qcd"  , 100, 199 );
     top   = merge( EventYield, regionsSet, "top"  , 300, 499 );
     wjets = merge( EventYield, regionsSet, "wjets", 500, 599 );
     zjets = merge( EventYield, regionsSet, "zjets", 600, 699 );
     data ->finalize();
     qcd  ->finalize();
     top  ->finalize();
     wjets->finalize();
     zjets->finalize();
     
     //qcd_QCDfromData = merge( EventYield_QCDfromData, regionsSet, "qcd_QCDfromData"  , 100, 199 );
     qcd_QCDfromData = merge( EventYield_QCDfromData, regionsSet, "qcd_QCDfromData"  , 1, 99 );
     qcd_QCDfromData->finalize();
  }
  else{

    qcd   = MT2Analysis<MT2Estimate>::readFromFile(dir + "/yields.root", "qcd");
    top   = MT2Analysis<MT2Estimate>::readFromFile(dir + "/yields.root", "top");
    wjets = MT2Analysis<MT2Estimate>::readFromFile(dir + "/yields.root", "wjets");
    zjets = MT2Analysis<MT2Estimate>::readFromFile(dir + "/yields.root", "zjets");
    
    if (!recomputeData ) {
      qcd_QCDfromData = MT2Analysis<MT2Estimate>::readFromFile(dir + "/yields.root", "qcd_QCDfromData");
      data  = MT2Analysis<MT2Estimate>::readFromFile(dir + "/yields.root", "data");
    }
    else {
      for( unsigned i=0; i<samples.size(); ++i ) {
        if (samples[i].id>99) continue;
        MT2Analysis<MT2Estimate>* analysis = new MT2Analysis<MT2Estimate>( samples[i].sname, regionsSet, samples[i].id );
        MT2Analysis<MT2Estimate>* analysis_QCDfromData = new MT2Analysis<MT2Estimate>( samples[i].sname, regionsSet, samples[i].id );
        computeYield( samples[i], regionsSet, analysis, analysis_QCDfromData, dataQCD, dataQCDps180, dataQCDps60, lumi);
        EventYield.push_back(analysis);
        EventYield_QCDfromData.push_back(analysis_QCDfromData);
      }
      data  = merge( EventYield, regionsSet, "data" ,   1,  99 );
      data ->finalize();
      qcd_QCDfromData = merge( EventYield_QCDfromData, regionsSet, "qcd_QCDfromData"  , 1, 99 );
      qcd_QCDfromData->finalize();
    }

  }

  *allBkg = *qcd + *top + *wjets + *zjets;
  allBkg->finalize();

  bgYields.push_back( qcd );
  bgYields.push_back( wjets );
  bgYields.push_back( zjets );
  bgYields.push_back( top );

  bgYields_QCDfromData.push_back( qcd_QCDfromData );
  bgYields_QCDfromData.push_back( wjets           );
  bgYields_QCDfromData.push_back( zjets           );
  bgYields_QCDfromData.push_back( top             );

  drawYields( dir, data, bgYields, "_QCDfromMC"   );
  drawYields( dir, data, bgYields_QCDfromData, "_QCDfromData" );

  if (recompute){
    allBkg->writeToFile( dir + "/yields.root" );
    qcd_QCDfromData->addToFile( dir + "/yields.root" );
    qcd  ->addToFile( dir + "/yields.root" );
    top  ->addToFile( dir + "/yields.root" );
    wjets->addToFile( dir + "/yields.root" );
    zjets->addToFile( dir + "/yields.root" );
    data->addToFile ( dir + "/yields.root" );
  }



}


MT2Analysis<MT2Estimate>* merge( std::vector<MT2Analysis<MT2Estimate> *> anas, const std::string& regionsSet, const std::string& name, int id_min, int id_max ) {

  if( id_max<0 ) id_max=id_min;

  MT2Analysis<MT2Estimate>* ana = new MT2Analysis<MT2Estimate>(name, regionsSet);

  for( unsigned i=0; i<anas.size(); ++i ) {
    if( anas[i]->getId() >= id_min && anas[i]->getId() <= id_max )
      *(ana) += *(anas[i]);
  } 

  return ana;

}


void computeYield( const MT2Sample& sample, const std::string& regionsSet,
		   MT2Analysis<MT2Estimate>* analysis,
		   MT2Analysis<MT2Estimate>* analysis_QCDfromData, 
		   MT2Analysis<MT2EstimateQCD>* dataQCD, MT2Analysis<MT2EstimateQCD>* dataQCDpsLow, MT2Analysis<MT2EstimateQCD>* dataQCDpsMed,
		   float lumi ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
  
  // In absence of skimmed ntuples let's filter to gain some speed
  TString filter = "met_pt>30&&ht>450";
  tree->Draw(">>selList", filter);
  TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
  tree->SetEventList(myEvtList);

  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  //analysis             = new MT2Analysis<MT2Estimate>( sample.sname, regionsSet, sample.id );
  //analysis_QCDfromData = new MT2Analysis<MT2Estimate>( sample.sname, regionsSet, sample.id );

  // int nentries = tree->GetEntries();
  // for( int iEntry=0; iEntry<nentries; ++iEntry ) {

  //   if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
  int nentries = myEvtList->GetN();
  for( int jEntry=0; jEntry<nentries; ++jEntry ) {
    int iEntry = myEvtList->GetEntry(jEntry);

    if( jEntry % 50000 == 0 ) std::cout << "    Entry: " << jEntry << " / " << nentries << std::endl;
    myTree.GetEntry(iEntry);

    // no dphi cut
    if ( !(myTree.met_pt>30 && myTree.nVert>0 && myTree.nJet30>=2)   ) continue;
    if ( !(myTree.passLeptonVeto() && myTree.passIsoTrackVeto())     ) continue;
    if ( doDiffMetMht && !(myTree.diffMetMht < 0.5*myTree.met_pt)    ) continue;
    if ( myTree.isData==1 && ( !myTree.isGolden || 
			       !(myTree.Flag_CSCTightHaloFilter==1 && 
				 myTree.Flag_HBHENoiseFilter   ==1 &&
				 myTree.Flag_eeBadScFilter     ==1) ) )  continue;
    

    if (doHFjetVeto) {
      int nJetHF30 = 0;
      for(int j=0; j<myTree.njet; ++j){
	if( myTree.jet_pt[j] < 30. || fabs(myTree.jet_eta[j]) < 3.0 ) continue;
	else ++nJetHF30;
      }
      if ( !(nJetHF30==0))  continue;
    }

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = myTree.mt2;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20;    
    float dphi = myTree.deltaPhiMin;
    int isData = myTree.isData;

    if (isData && ( (ht<575  && myTree.HLT_ht350prescale==0) ||
		    (ht>575 && ht<1000 && myTree.HLT_ht475prescale==0) ) )  continue;


    if (met>30) met = 200.;  // don't do met>200 for low HT regions (do met>30)
    //if (ht<1000 && met>30) met = 200.;  // don't do met>200 for low HT regions (do met>30)

    Double_t weight = isData ? 1.0 : myTree.evt_scale1fb*lumi;

    if (dphi > 0.3){ 
      MT2Estimate* thisEstimate = analysis->get( ht, njets, nbjets,  minMTBmet, mt2 );
      if( thisEstimate==0 ) continue;
      thisEstimate->yield->Fill(mt2, weight );
    }
    else{
      MT2Estimate* thisEstimate = analysis_QCDfromData->get( ht, njets, nbjets, minMTBmet, mt2 );
      if( thisEstimate==0 ) continue;
      float r;
      if (ht>1000)
	r = dataQCD->get( ht, njets, nbjets, minMTBmet, mt2 )->exp->Eval(mt2);
      else if (ht>575)
	r = dataQCDpsMed->get( ht, njets, nbjets, minMTBmet, mt2 )->exp->Eval(mt2);
      else
	r = dataQCDpsLow->get( ht, njets, nbjets, minMTBmet, mt2 )->exp->Eval(mt2);

      thisEstimate->yield->Fill(mt2, r*weight );
    }
    
  } // for entries

  analysis->finalize();
  analysis_QCDfromData->finalize();
  
  delete tree;

  file->Close();
  delete file;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::vector< MT2Analysis<MT2Estimate> *> bgYields,  const std::string& prefix) {


  MT2DrawTools::setStyle();

  std::vector<int> colors;
  if( bgYields.size()==3 ) { // estimates
    colors.push_back(402); 
    colors.push_back(430); 
    colors.push_back(418); 
  } else { // mc
    colors.push_back(401); // qcd
    colors.push_back(417); // w+jets
    colors.push_back(419); // z+jets
    colors.push_back(855); // top
    //colors.push_back(); // other
  }



  std::string fullPath = outputdir;
  std::string fullPathPlots = outputdir + "/plots";
  system( Form("mkdir -p %s", fullPathPlots.c_str()) );

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  

    //std::string fullPath = outputdir + "/" + iHT->getName() + "/" + iSR->getName();
    //std::string mkdircommand = "mkdir -p " + fullPath;
    //system( mkdircommand.c_str() );


    MT2Region thisRegion( (*iMT2) );


    TH1D* h1_data = (TH1D*)data->get(thisRegion)->yield->Clone();
    h1_data->SetName(Form("%s%s",h1_data->GetName(),prefix.c_str()));
    // blinded data
    for (int b = 0; b<=h1_data->GetNbinsX(); b++) {
      if (h1_data->GetBinLowEdge(b)+h1_data->GetBinWidth(b)>200)
	h1_data->SetBinContent(b,0);
    }

    TFile* histoFile = TFile::Open( Form("%s/histograms%s_%s.root", fullPath.c_str(), prefix.c_str(), thisRegion.getName().c_str()), "recreate" );
    histoFile->cd();
    h1_data->Write();

    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.2);

    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
      //for( int i=bgYields.size()-1; i>=0;  --i ) { 
      int index = bgYields.size() - i - 1;
      TH1D* h1_bg = (TH1D*) bgYields[index]->get(thisRegion)->yield->Clone();
      h1_bg->SetName(Form("%s%s",h1_bg->GetName(),prefix.c_str()));
      h1_bg->SetFillColor( colors[index] );
      h1_bg->SetLineColor( kBlack );
      
      if (prefix.find("fromData") != std::string::npos ){
      	for (int b = 0; b<=6; b++) {
      	  h1_bg->SetBinContent(b,0);
      	}
      }
      
      // apply prescales to mc components; ps values hardcoded for now
      std::cout << "bkg is: " << bgYields[index]->getFullName() << std::endl;
      if ( !(((TString)bgYields[index]->getFullName()).Contains("QCDfromData")) ) {
	std::cout << "entering mc component: " << bgYields[index]->getFullName() << std::endl;
	if ( ((TString)thisRegion.getName()).Contains("HT450to575") ){
	  std::cout << "prescale 180: " << thisRegion.getName() << std::endl;
	  std::cout << "integral before prescale: " << h1_bg->Integral() << std::endl;
	  h1_bg->Scale(1./180.);
	  std::cout << "integral after prescale: " << h1_bg->Integral() << std::endl;
	}
	else if ( ((TString)thisRegion.getName()).Contains("HT575to1000") ) {
	  std::cout << "prescale 60: " << thisRegion.getName() << std::endl;
	  std::cout << "integral before prescale: " << h1_bg->Integral() << std::endl;
	  h1_bg->Scale(1./60.);
	  std::cout << "integral after prescale: " << h1_bg->Integral() << std::endl;
	}
      }	

      bgStack.Add(h1_bg);
    }

    TH1D* h1_ratio = (TH1D*)h1_data->Clone("h1_ratio");
    h1_ratio->Sumw2();
    for( int iBin=1; iBin<(h1_ratio->GetXaxis()->GetNbins()+1); ++iBin ){
      h1_ratio->SetBinContent(iBin, TMath::Nint(h1_ratio->GetBinContent(iBin)));
      h1_ratio->SetBinError(iBin, TMath::Sqrt(TMath::Nint(h1_ratio->GetBinContent(iBin))));
    }
    h1_ratio->Divide((TH1D*)bgStack.GetStack()->Last());
    //h1_ratio->SetBinErrorOption(TH1::kPoisson);
    

    TCanvas* c1 = new TCanvas( "canny", "", 600, 700 );
    //TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    //c1->SetLogy();
    c1->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
    pad1->SetBottomMargin(0.15);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();


    float xMin = h1_data->GetXaxis()->GetXmin();
    float xMax = h1_data->GetXaxis()->GetXmax();
    xMax = 300;
    float multiplier = pad1->GetLogy() ? 10 : 1.2;
    float yMax1 = h1_data->GetMaximum()*multiplier;
    float yMax2 = multiplier*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
    float yMax3 = multiplier*(bgStack.GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    //float yMax = TMath::Max( h1_data->GetMaximum()*1.5, (h1_data->GetMaximum() + h1_data->GetBinError(h1_data->GetMaximumBin()))*1.2);
    //float yMax = h1_data->GetMaximum()*1.5;
    if( h1_data->GetNbinsX()<2 ) yMax *=3.;

    //TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax ); \\ TH2D is not a good option for logy
    TH1D* h2_axes = new TH1D("axes", "", 10, xMin, xMax );
    h2_axes->SetXTitle("M_{T2} [GeV]");
    h2_axes->SetYTitle("Entries");
    if (pad1->GetLogy()) h2_axes->GetYaxis()->SetRangeUser(0.1, yMax);
    else                 h2_axes->GetYaxis()->SetRangeUser(0.0, yMax);

    h2_axes->Draw();
   


    std::vector<std::string> niceNames = thisRegion.getNiceNames();

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
    

    TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( gr_data, "Data", "PE" );
    histoFile->cd();
    for( unsigned i=0; i<bgYields.size(); ++i ) {  
      TH1D* h1_bg = (TH1D*)bgYields[i]->get(thisRegion)->yield->Clone();
      h1_bg->SetName(Form("leg%s%s",h1_bg->GetName(),prefix.c_str()));
      h1_bg->SetFillColor( colors[i] );
      h1_bg->SetLineColor( kBlack );
      TString fullName = bgYields[i]->getFullName().c_str();
      if (fullName.Contains("QCDfromData"))
	fullName = "QCD from data";
      else if (fullName.Contains("qcd"))
	fullName = "QCD from MC";	
      fullName = fullName.Contains("QCDfromData") ? "QCD from data" : (fullName.Contains("qcd") ? "qcd" : (fullName.Contains("zjets") ? "Z+jets" : (fullName.Contains("wjets") ? "W+jets" : (fullName.Contains("top") ? "Top" : fullName ))));
      legend->AddEntry( h1_bg, fullName, "F" );
      h1_bg->Write();
    }

    histoFile->Close();

    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");

    TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
    labelTop->Draw("same");

    gPad->RedrawAxis();

    c1->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.1);
    pad2->Draw();
    pad2->cd();


    TH2D* h2_axes_rat = new TH2D("axes_rat", "", 10, xMin, xMax, 5 , 0.,2  );
    h2_axes_rat->SetYTitle("Data / MC");

    h2_axes_rat->GetXaxis()->SetTitleSize(0.2);
    h2_axes_rat->GetXaxis()->SetTitleOffset(5);
    h2_axes_rat->GetXaxis()->SetLabelSize(0.00);
    h2_axes_rat->GetXaxis()->SetTickLength(0.09);
    h2_axes_rat->GetYaxis()->SetNdivisions(5,5,0);
    h2_axes_rat->GetYaxis()->SetTitleSize(0.2);
    h2_axes_rat->GetYaxis()->SetTitleOffset(0.34);
    h2_axes_rat->GetYaxis()->SetLabelSize(0.17); 

    h1_ratio->SetLineColor(1);
    h1_ratio->SetMarkerStyle(20);
    h1_ratio->SetMarkerSize(1);

    h2_axes_rat->Draw();
    h1_ratio->Draw("pesame");

    TLine *line1 = new TLine(xMin, 1.0, xMax, 1.0);
    line1->SetLineStyle(2);
    line1->Draw("same");

    gPad->RedrawAxis();

    c1->SaveAs( Form("%s/mt2%s_%s.eps", fullPathPlots.c_str(), prefix.c_str() , thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/mt2%s_%s.png", fullPathPlots.c_str(), prefix.c_str() , thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/mt2%s_%s.pdf", fullPathPlots.c_str(), prefix.c_str() , thisRegion.getName().c_str()) );

    delete c1;
    delete h2_axes;
    delete h1_ratio;
  }// for MT2 regions

}
