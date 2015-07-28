#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGaxis.h"
#include "THStack.h"
#include "TLorentzVector.h"

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2Config.h"


#include "../interface/MT2DrawTools.h"


#include <iostream>


#define mt2_cxx
#include "interface/mt2.h"


//float lumi = 0.1;
//float lumi = 1.;


//This file creates the Zll trees used to estimate the backgrounds in the 
//Zll control region.
//Then run zllPurity to get the nice plots and figures.

void drawYields( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data_zll, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields, float lumi );

MT2Analysis<MT2EstimateTree>* computeYield( const MT2Sample& sample, const MT2Config& cfg, float lumi=1., bool doSameFlavor=1 );

MT2Analysis<MT2EstimateTree>* mergeYields( std::vector< MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max=-1, const std::string& legendName="" );


int main(int argc, char* argv[]){

  std::string regionsSet = "zurich";
  if( argc>2 ) {
    regionsSet = std::string(argv[2]);
  }

  if( argc<2 ) {
    std::cout << "USAGE: ./coputeZllGammaRatio [configFileName] regionSet" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
 

  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);
  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat"; 
  std::string samples = cfg.mcSamples();

  regionsSet = cfg.regionsSet();


  std::string outputdir( Form("ZllPurity_%s", configFileName.c_str() ) );
  std::string outputdir_of( Form("ZllPurity_OF_%s", configFileName.c_str()) );

 
  std::cout << "-> Using regions: " << regionsSet << std::endl;


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  //std::string outputdir = "ZllGammaRatio_"+ cfg.mcSamples + regionsSet ;
  //  std::string outputdir = "Zll_" + configFileName;
  double intpart;
  double fracpart = modf(cfg.lumi(), &intpart);
  std::string suffix;
  if( fracpart>0. )
    suffix = std::string( Form("_%.0fp%.0ffb", intpart, 10.*fracpart ) );
  else
    suffix = std::string( Form("_%.0ffb", intpart ) );
  //  outputdir += suffix;
  //  outputdir_of += suffix;
  
  system(Form("mkdir -p %s", outputdir.c_str()));


 
  std::cout << "-> Using regions: " << regionsSet << std::endl;
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;



  //DATA
  std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;
  std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data, "Double");
  // std::vector<MT2Sample> samples_data_of = MT2Sample::loadSamples(samplesFile_data, "MuonEG");

  std::vector< MT2Analysis<MT2EstimateTree>* > dataTree;
  std::vector< MT2Analysis<MT2EstimateTree>* > dataTree_of;

  if( samples_data.size()==0 ) {
    std::cout << std::endl;
    std::cout << "-> WARNING!! Didn't find any data in file: " << samplesFile_data << "!" << std::endl;
    std::cout << "-> Exiting." << std::endl;
    std::cout << std::endl;
  } else {
 
    // = new MT2Analysis<MT2EstimateTree>( "zllCRtree", cfg.regionsSet() );   
    for( unsigned i=0; i<samples_data.size(); ++i ) {
      dataTree.push_back( computeYield( samples_data[i], cfg, cfg.lumi(),1 ));
      dataTree_of.push_back( computeYield( samples_data[i], cfg, cfg.lumi(),0 ));
    }
  }

  MT2Analysis<MT2EstimateTree>* EventYield_data = mergeYields( dataTree, cfg.regionsSet(), "data", 0, 2000, "" );
  MT2Analysis<MT2EstimateTree>* EventYield_data_of = mergeYields( dataTree_of, cfg.regionsSet(), "data_of", 0, 2000, "" );


  //MC
  std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 1, 999); // not interested in signal here
  if( fSamples.size()==0 ) {
    std::cout << "There must be an error: samples is empty!" << std::endl;
    exit(1209);
  }
  
  std::vector< MT2Analysis<MT2EstimateTree>* > EventYield;
  for( unsigned i=0; i<fSamples.size(); ++i ) 
    EventYield.push_back( computeYield( fSamples[i], cfg, cfg.lumi(), 1 ) );
    
 
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields; 

  MT2Analysis<MT2EstimateTree>* EventYield_zll = mergeYields( EventYield, cfg.regionsSet(), "DYJets", 700, 799, "DYJets" );
 

  MT2Analysis<MT2EstimateTree>* EventYield_top   = mergeYields( EventYield, cfg.regionsSet(), "Top", 300, 499 ); // ttbar, single top, ttW, ttZ...
  MT2Analysis<MT2EstimateTree>* EventYield_qcd   = mergeYields( EventYield, cfg.regionsSet(), "QCD", 100, 199 );
  MT2Analysis<MT2EstimateTree>* EventYield_wjets = mergeYields( EventYield, cfg.regionsSet(), "WJets", 500, 599, "W+jets" );
  MT2Analysis<MT2EstimateTree>* EventYield_zjets = mergeYields( EventYield, cfg.regionsSet(), "ZJets", 600, 699, "Z+jets" );

  bgYields.push_back( EventYield_qcd );
  bgYields.push_back( EventYield_wjets );
  bgYields.push_back( EventYield_zjets );
  bgYields.push_back( EventYield_top );
 

 drawYields( outputdir, EventYield_zll, bgYields, cfg.lumi() );
 
  std::string outFile = outputdir + "/ZllPurityTrees.root";

  EventYield_zll->writeToFile( outFile );
  EventYield_top->addToFile( outFile );
  EventYield_qcd->addToFile( outFile );
  EventYield_wjets->addToFile( outFile );
  EventYield_zjets->addToFile( outFile );

  std::string outFile_data = outputdir + "/ZllPurityTrees_data.root";
  EventYield_data->writeToFile(outFile_data);

  
 
  
  std::vector< MT2Analysis<MT2EstimateTree>* > EventYield_of;
  for( unsigned i=0; i<fSamples.size(); ++i ) 
    EventYield_of.push_back( computeYield( fSamples[i], cfg, cfg.lumi(), 0 ) );
    
  MT2Analysis<MT2EstimateTree>* EventYield_zll_of = mergeYields( EventYield_of, cfg.regionsSet(), "DYJets", 700, 799, "DYJets" );

  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields_of; 

  MT2Analysis<MT2EstimateTree>* EventYield_top_of   = mergeYields( EventYield_of, cfg.regionsSet(), "Top", 300, 499 ); // ttbar, single top, ttW, ttZ...
  MT2Analysis<MT2EstimateTree>* EventYield_qcd_of   = mergeYields( EventYield_of, cfg.regionsSet(), "QCD", 100, 199 );
  MT2Analysis<MT2EstimateTree>* EventYield_wjets_of = mergeYields( EventYield_of, cfg.regionsSet(), "WJets", 500, 599, "W+jets" );
  MT2Analysis<MT2EstimateTree>* EventYield_zjets_of = mergeYields( EventYield_of, cfg.regionsSet(), "ZJets", 600, 699, "Z+jets" );
 
  bgYields_of.push_back( EventYield_qcd_of );
  bgYields_of.push_back( EventYield_wjets_of );
  bgYields_of.push_back( EventYield_zjets_of );
  bgYields_of.push_back( EventYield_top_of );

  drawYields( outputdir_of, EventYield_zll_of, bgYields_of, cfg.lumi() );
 
 

 
  
  std::string outFile_of = outputdir_of + "/ZllPurityTrees_of.root";

  EventYield_zll_of->writeToFile( outFile_of );
  EventYield_top_of->addToFile( outFile_of );
  EventYield_qcd_of->addToFile( outFile_of );
  EventYield_wjets_of->addToFile( outFile_of );
  EventYield_zjets_of->addToFile( outFile_of );

  std::string outFile_data_of = outputdir_of + "/ZllPurityTrees_data_of.root";
  EventYield_data_of->writeToFile(outFile_data_of);
  

  return 0;
}

























void drawYields( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields, float lumi ) {


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
    MT2Region thisRegion( (*iMT2) );

    TH1D* h1_data = data->get(thisRegion)->yield;

    TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), thisRegion.getName().c_str()), "recreate" );
    histoFile->cd();
    h1_data->Write();


    /*
    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.6);
    */
    //   h1_data->SetMarkerStyle(20);
    //   h1_data->SetMarkerSize(1.6);


    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
      int index = bgYields.size() - i - 1;
      TH1D* h1_bg = bgYields[index]->get(thisRegion)->yield;
      h1_bg->SetFillColor( colors[index] );
      h1_bg->SetLineColor( kBlack );
      bgStack.Add(h1_bg);
    }

    h1_data->SetLineColor(kBlack);
    h1_data->SetFillColor(430);
    bgStack.Add(h1_data);

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();

    float xMin = h1_data->GetXaxis()->GetXmin();
    float xMax = h1_data->GetXaxis()->GetXmax();
    float yMax1 = h1_data->GetMaximum()*1.5;
    float yMax2 = 1.5*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
    float yMax3 = 1.5*(bgStack.GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    if( h1_data->GetNbinsX()<2 ) yMax *=3.;

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle("M_{T2} [GeV]");
    h2_axes->SetYTitle("Entries");

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
    

    TLegend* legend = new TLegend( 0.7, 0.9-(bgYields.size()+1)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);

    legend->AddEntry( h1_data, "Zll", "F" );
    //      legend->AddEntry( gr_data, "Zll", "P" );
    histoFile->cd();
    for( unsigned i=0; i<bgYields.size(); ++i ) {  
      TH1D* h1_bg = bgYields[i]->get(thisRegion)->yield;
      legend->AddEntry( h1_bg, bgYields[i]->getFullName().c_str(), "F" );
      h1_bg->Write();
    }


    histoFile->Close();



    legend->Draw("same");
    bgStack.Draw("histo same");
    // h1_data->Draw("p same");
    //    gr_data->Draw("p same");


    TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
    labelTop->Draw("same");


    gPad->RedrawAxis();

    c1->SaveAs( Form("%s/mt2_%s.eps", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/mt2_%s.png", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/mt2_%s.pdf", fullPathPlots.c_str(), thisRegion.getName().c_str()) );
 

    delete c1;
    delete h2_axes;


  }// for MT2 regions


}























MT2Analysis<MT2EstimateTree>* computeYield( const MT2Sample& sample, const MT2Config& cfg, float lumi, bool doSameFlavor ) {

  std::string regionsSet = cfg.regionsSet();

  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  std::cout << "-> Getting mt2 tree from file: " << sample.file << std::endl;

  TTree* tree = (TTree*)file->Get("mt2");
  

  MT2Tree myTree;
  if( cfg.additionalStuff()=="qgVars" ) {
     myTree.loadGenStuff = true;
  } else {
    myTree.loadGenStuff = false;
  }
  myTree.Init(tree);



  std::cout << "-> Setting up MT2Analysis with name: " << sample.sname << std::endl;
  MT2Analysis<MT2EstimateTree>* analysis = new MT2Analysis<MT2EstimateTree>( sample.sname, regionsSet, sample.id );


  MT2EstimateTree::addVar( analysis, "Z_pt" );
  MT2EstimateTree::addVar( analysis, "Z_phi" );
  MT2EstimateTree::addVar( analysis, "Z_mass" );
  MT2EstimateTree::addVar( analysis, "Z_lepId" );
  MT2EstimateTree::addVar( analysis, "nLep" );
  MT2EstimateTree::addVar( analysis, "sample_Id");
  MT2EstimateTree::addVar( analysis, "lep_pt0");
  MT2EstimateTree::addVar( analysis, "lep_pt1");
  MT2EstimateTree::addVar( analysis, "lep_eta0");
  MT2EstimateTree::addVar( analysis, "lep_eta1");
  MT2EstimateTree::addVar( analysis, "raw_mt2");
  
  MT2EstimateTree::addVar( analysis, "HLT_DoubleMu");
  MT2EstimateTree::addVar( analysis, "HLT_DoubleEl");
  


  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "  Entry: " << iEntry << " / " << nentries << std::endl;
    myTree.GetEntry(iEntry);

    // if(!(myTree.passSelection("zll"))) continue; 
    //this the the baseline selection:
    if(myTree.nVert < 0) continue;
    if(myTree.nJet30 < 2  ) continue;
    if(myTree.zll_deltaPhiMin < 0.3) continue;
    if(myTree.zll_diffMetMht > 0.5*myTree.zll_met_pt) continue;

    if(!(myTree.nlep==2)) continue; 
    
    if(myTree.mt2>200) continue; //change once we have more data

    //Sample  are the Z leptons
    //and thus that if they don't have the same flavor they are rejected
 
    if(doSameFlavor==1 && !(myTree.lep_pdgId[0] == -myTree.lep_pdgId[1]) )   continue;
   
    if(doSameFlavor==0 && (myTree.lep_pdgId[0] == -myTree.lep_pdgId[1])) continue;

    if(( myTree.lep_pdgId[0]*myTree.lep_pdgId[1])>0 )   continue;
    
    if(myTree.lep_pt[0]<25) continue;
    if(myTree.lep_pt[1]<20) continue;

    if(  doSameFlavor==1 && !(myTree.HLT_DoubleMu || myTree.HLT_DoubleEl) ) continue;
    //  if( myTree.isData && doSameFlavor==0 && !(myTree.HLT_MuEl) ) continue;
    if( doSameFlavor==0 && !(myTree.HLT_DoubleMu || myTree.HLT_DoubleEl ) ) continue;

    //Need the lorentz vectors of the leptons first
    TLorentzVector *LVec = new TLorentzVector[5];
    for(int i=0; i< 2; i++){
      LVec[i].SetPtEtaPhiM(myTree.lep_pt[i], myTree.lep_eta[i],myTree.lep_phi[i], myTree.lep_mass[i]);
    }

    double Z_invM_true = 91.19;
    TLorentzVector z = LVec[0] + LVec[1]; //leptons invariant mass
    double M_ll = z.M(); //Z mass

    //   if( abs(M_ll - Z_invM_true)>20.) continue;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = myTree.mt2;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet30;
    int nbjets = myTree.nBJet20;

    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi(); 

    MT2EstimateTree* thisEstimate = analysis->get( myTree.zll_ht, njets, nbjets, myTree.zll_met_pt, minMTBmet, myTree.zll_mt2 );
    if( thisEstimate==0 ) continue; 

    //initialize
    thisEstimate->assignVar("Z_pt", z.Perp() );
    thisEstimate->assignVar("Z_phi", z.Phi() );
    thisEstimate->assignVar("Z_mass", z.M() );
    thisEstimate->assignVar("Z_lepId", abs(myTree.lep_pdgId[0])  );
    thisEstimate->assignVar("sample_Id", myTree.evt_id);
    thisEstimate->assignVar("nLep", myTree.nlep );
    thisEstimate->assignVar("lep_pt0", myTree.lep_pt[0] );
    thisEstimate->assignVar("lep_pt1", myTree.lep_pt[1] );
    thisEstimate->assignVar("lep_eta0", myTree.lep_eta[0] );
    thisEstimate->assignVar("lep_eta1", myTree.lep_eta[1] );
    thisEstimate->assignVar("raw_mt2", myTree.mt2 );

    thisEstimate->assignVar("HLT_DoubleMu", myTree.HLT_DoubleMu );
    thisEstimate->assignVar("HLT_DoubleEl", myTree.HLT_DoubleEl );

    thisEstimate->fillTree_zll(myTree, weight);

    thisEstimate->yield->Fill(myTree.zll_mt2, weight );
  
  } // for entries

  //ofs.close();

  analysis->finalize();
  
  delete tree;

  file->Close();
  delete file;
  
  return analysis;
}



























MT2Analysis<MT2EstimateTree>* mergeYields( std::vector<MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max, const std::string& legendName ) {

  if( id_max<0 ) id_max=id_min;
  MT2Analysis<MT2EstimateTree>* return_EventYield = new MT2Analysis<MT2EstimateTree>(name, regionsSet, id_min, legendName);

  for( unsigned i=0; i<EventYield.size(); ++i ) {
    if( EventYield[i]->getId() >= id_min && EventYield[i]->getId() <= id_max ) {
       *(return_EventYield) += *(EventYield[i]);
     }

  } // for EventYield

  return return_EventYield;
}
