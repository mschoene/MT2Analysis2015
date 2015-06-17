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
#include "interface/MT2EstimateSigTree.h"
#include "interface/MT2DrawTools.h"

#include "TRandom3.h"


#define mt2_cxx
#include "interface/mt2.h"



bool dummyAnalysis;
double lumi = 4; // in fb-1



class MT2Config {

 public:

  MT2Config( const std::string& configFileName );

  std::string regionsSet()      const { return regionsSet_; };
  std::string mcSamples()       const { return mcSamples_; };
  std::string sigSamples()      const { return sigSamples_; };
  std::string dataSamples()     const { return dataSamples_; };
  std::string additionalStuff() const { return additionalStuff_; };

  bool useMC() {
    return mcSamples()!="";
  }

 private:

  std::string regionsSet_;
  std::string mcSamples_;
  std::string sigSamples_;
  std::string dataSamples_;
  std::string additionalStuff_;

};





void randomizePoisson( MT2Analysis<MT2EstimateTree>* data );
MT2Analysis<MT2EstimateTree>* computeYield( const MT2Sample& sample, const MT2Config& cfg, float lumi=1. );
MT2Analysis<MT2EstimateSigTree>* computeSigYield( const MT2Sample& sample, const MT2Config& cfg, float lumi=1. );
MT2Analysis<MT2EstimateTree>* mergeYields( std::vector< MT2Analysis<MT2EstimateTree> *> EventYield, const std::string& regionsSet, const std::string& name, int id_min, int id_max=-1, const std::string& legendName="" );
int matchPartonToJet( int index, MT2Tree* myTree );


void drawYields( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields );



int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./regionEventYields [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }




  std::string configFileName(argv[1]);

  MT2Config cfg("cfgs/" + configFileName + ".txt");
  dummyAnalysis = cfg.dataSamples()=="datatest";

  std::string outputdir = "EventYields_" + configFileName;
  if( dummyAnalysis ) {
    double intpart;
    double fracpart = modf(lumi, &intpart);
    std::string suffix;
    if( fracpart>0. )
      suffix = std::string( Form("_dummy_%.0fp%.0ffb", intpart, 10.*fracpart ) );
    else
      suffix = std::string( Form("_dummy_%.0ffb", intpart ) );
    outputdir += suffix;
  }

    
  system(Form("mkdir -p %s", outputdir.c_str()));



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields;



  if( cfg.useMC() ) { // use MC BG estimates


    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 1, 999); // not interested in signal here (see later)
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(1209);
    }


    
    std::vector< MT2Analysis<MT2EstimateTree>* > EventYield;
    for( unsigned i=0; i<fSamples.size(); ++i ) {
      int this_id = fSamples[i].id;
      if( this_id>=200 && this_id<300 ) continue; // skip GJets
      EventYield.push_back( computeYield( fSamples[i], cfg, lumi ) );
    }
    


    std::cout << "-> Done looping on samples. Start merging." << std::endl;

    MT2Analysis<MT2EstimateTree>* EventYield_top   = mergeYields( EventYield, cfg.regionsSet(), "Top", 300, 499 ); // ttbar, single top, ttW, ttZ...
    MT2Analysis<MT2EstimateTree>* EventYield_qcd   = mergeYields( EventYield, cfg.regionsSet(), "QCD", 100, 199 );
    MT2Analysis<MT2EstimateTree>* EventYield_wjets = mergeYields( EventYield, cfg.regionsSet(), "WJets", 500, 599, "W+jets" );
    MT2Analysis<MT2EstimateTree>* EventYield_zjets = mergeYields( EventYield, cfg.regionsSet(), "ZJets", 600, 699, "Z+jets" );
    //MT2Analysis<MT2EstimateTree>* EventYield_other = mergeYields( EventYield, cfg.regionsSet(), "Diboson", 700, 899, "Other" );

    bgYields.push_back( EventYield_qcd );
    bgYields.push_back( EventYield_wjets );
    bgYields.push_back( EventYield_zjets );
    bgYields.push_back( EventYield_top );
    //bgYields.push_back( EventYield_other );

  }

  // load signal samples, if any
  std::vector< MT2Analysis< MT2EstimateSigTree>* > signals;
  if( cfg.mcSamples()!="" ) {

    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading signal samples from file: " << samplesFileName << std::endl;

    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 1000); // only signal (id>=1000)


    if( fSamples.size()==0 ) {

      std::cout << "No signal samples found, skipping." << std::endl;

    } else {
    
      for( unsigned i=0; i<fSamples.size(); ++i ) 
        signals.push_back( computeSigYield( fSamples[i], cfg, lumi ) );
    
    } // if samples != 0

  } // if mc samples

  else if ( cfg.sigSamples()!="" ){

    std::string samplesFileName = "../samples/samples_" + cfg.sigSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading signal samples from file: " << samplesFileName << std::endl;

    std::vector<MT2Sample> fSamples = MT2Sample::loadSamples(samplesFileName, 1000); // only signal (id>=1000)

    if( fSamples.size()==0 ) {

      std::cout << "No signal samples found, skipping." << std::endl;

    } else {

      for( unsigned i=0; i<fSamples.size(); ++i )
        signals.push_back( computeSigYield( fSamples[i], cfg, lumi ) );

    } // if samples != 0
    
  } // if sig samples
  

  //MT2Analysis<MT2EstimateTree>* data = new MT2Analysis<MT2EstimateTree>( "data", cfg.regionsSet() );
  MT2Analysis<MT2EstimateTree>* data;
 
  if( dummyAnalysis ) { // use same as MC

    data = new MT2Analysis<MT2EstimateTree>( *(bgYields[0]) );
    data->setName("data");
    for( unsigned i=1; i < bgYields.size(); ++i ) (*data) += *(bgYields[i]);
    //randomizePoisson( data );

  } else {

    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    std::vector<MT2Sample> samples_data = MT2Sample::loadSamples(samplesFile_data);
    if( samples_data.size()==0 ) {
      std::cout << "There must be an error: samples_data is empty!" << std::endl;
      exit(1209);
    }

    data = new MT2Analysis<MT2EstimateTree>( "data", cfg.regionsSet() );
    for( unsigned i=0; i<samples_data.size(); ++i ) (*data) += *(computeYield( samples_data[i], cfg ));

  }


  drawYields( outputdir, data, bgYields );


  // save MT2Analyses:
  data->writeToFile(outputdir + "/analyses.root");
  for( unsigned i=0; i<bgYields.size(); ++i )
    bgYields[i]->writeToFile(outputdir + "/analyses.root", "UPDATE");
  for( unsigned i=0; i<signals.size(); ++i )
    signals[i]->writeToFile(outputdir + "/analyses.root", "UPDATE");

  return 0;

}





MT2Analysis<MT2EstimateTree>* computeYield( const MT2Sample& sample, const MT2Config& cfg, float lumi ) {


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
  
 
  if( cfg.additionalStuff()=="qgVars" ) {
    MT2EstimateTree::addVar( analysis, "partId0" );
    MT2EstimateTree::addVar( analysis, "partId1" );
    MT2EstimateTree::addVar( analysis, "partId2" );
    MT2EstimateTree::addVar( analysis, "partId3" );
    MT2EstimateTree::addVar( analysis, "qgl0" );
    MT2EstimateTree::addVar( analysis, "qgl1" );
    MT2EstimateTree::addVar( analysis, "qgl2" );
    MT2EstimateTree::addVar( analysis, "qgl3" );
    MT2EstimateTree::addVar( analysis, "qglProd" );
    MT2EstimateTree::addVar( analysis, "qglAve" );
  }
  

  bool isData = sample.id<100 && sample.id>0;



  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( !myTree.passSelection() ) continue;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = myTree.mt2;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet40;
    int nbjets = myTree.nBJet20;    
    
    Double_t weight = (isData) ? 1. : myTree.evt_scale1fb*lumi;
    //weight *= myTree.weight_lepsf;
   
    MT2EstimateTree* thisEstimate = analysis->get( ht, njets, nbjets, met, minMTBmet, mt2 );
    if( thisEstimate==0 ) continue;


    if( cfg.additionalStuff()=="qgVars" ) {

      // initialize
      thisEstimate->assignVar( "qgl0", -1. );
      thisEstimate->assignVar( "qgl1", -1. );
      thisEstimate->assignVar( "qgl2", -1. );
      thisEstimate->assignVar( "qgl3", -1. );
      thisEstimate->assignVar( "partId0", 0 );
      thisEstimate->assignVar( "partId1", 0 );
      thisEstimate->assignVar( "partId2", 0 );
      thisEstimate->assignVar( "partId3", 0 );

      float qglProd = 1.;
      float qglAve = 0.;
      int denom = 0;


      if( njets>0 && fabs(myTree.jet_eta[0])<2.5 ) {

        float qgl0 = myTree.jet_qgl[0];
        thisEstimate->assignVar( "qgl0", qgl0 );
        qglProd *= qgl0;
        qglAve += qgl0;
        denom++;
        thisEstimate->assignVar( "partId0", matchPartonToJet( 0, &myTree ) );
        //thisEstimate->assignVar( "partId0", myTree.jet_mcFlavour[0] );

      }


      if( njets>1 && fabs(myTree.jet_eta[1])<2.5 ) {

        float qgl1 = myTree.jet_qgl[1];
        thisEstimate->assignVar( "qgl1", qgl1 );
        qglProd *= qgl1;
        qglAve += qgl1;
        denom++;

        thisEstimate->assignVar( "partId1", matchPartonToJet( 1, &myTree ) );
        //thisEstimate->assignVar( "partId1", myTree.jet_mcFlavour[1] );

      }
        
      if( njets>2 && fabs(myTree.jet_eta[2])<2.5 ) {

        float qgl2 = myTree.jet_qgl[2];
        thisEstimate->assignVar( "qgl2", qgl2 );
        qglProd *= qgl2;
        qglAve += qgl2;
        denom++;

        thisEstimate->assignVar( "partId2", matchPartonToJet( 2, &myTree ) );
        //thisEstimate->assignVar( "partId2", myTree.jet_mcFlavour[2] );

      }
        

      if( njets>3 && fabs(myTree.jet_eta[3])<2.5 ) {

        float qgl3 = myTree.jet_qgl[3];
        thisEstimate->assignVar( "qgl3", qgl3 );
        qglProd *= qgl3;
        qglAve += qgl3;
        denom++;

        thisEstimate->assignVar( "partId3", matchPartonToJet( 3, &myTree ) );
        //thisEstimate->assignVar( "partId3", myTree.jet_mcFlavour[3] );

      }

      qglAve /= (float)denom;
        
      thisEstimate->assignVar( "qglProd", qglProd );
      thisEstimate->assignVar( "qglAve", qglAve );

      thisEstimate->assignTree(myTree, weight );
      thisEstimate->tree->Fill();

    } else {

      thisEstimate->fillTree(myTree, weight );

    }

    thisEstimate->yield->Fill(mt2, weight );
    
  } // for entries
    
  //ofs.close();

  analysis->finalize();
  
  delete tree;

  file->Close();
  delete file;
  
  return analysis;

}



MT2Analysis<MT2EstimateSigTree>* computeSigYield( const MT2Sample& sample, const MT2Config& cfg, float lumi ) {


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
  MT2Analysis<MT2EstimateSigTree>* analysis = new MT2Analysis<MT2EstimateSigTree>( sample.sname, regionsSet, sample.id );
  
 
  if( cfg.additionalStuff()=="qgVars" ) {
    MT2EstimateSigTree::addVar( analysis, "partId0" );
    MT2EstimateSigTree::addVar( analysis, "partId1" );
    MT2EstimateSigTree::addVar( analysis, "partId2" );
    MT2EstimateSigTree::addVar( analysis, "partId3" );
    MT2EstimateSigTree::addVar( analysis, "qgl0" );
    MT2EstimateSigTree::addVar( analysis, "qgl1" );
    MT2EstimateSigTree::addVar( analysis, "qgl2" );
    MT2EstimateSigTree::addVar( analysis, "qgl3" );
    MT2EstimateSigTree::addVar( analysis, "qglProd" );
    MT2EstimateSigTree::addVar( analysis, "qglAve" );
  }
  

  bool isData = sample.id<100 && sample.id>0;



  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( !myTree.passSelection() ) continue;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = myTree.mt2;
    float minMTBmet = myTree.minMTBMet;
    int njets  = myTree.nJet40;
    int nbjets = myTree.nBJet20;    
    
    float GenSusyMScan1 = myTree.GenSusyMScan1;
    float GenSusyMScan2 = myTree.GenSusyMScan2;

    Double_t weight = (isData) ? 1. : myTree.evt_scale1fb*lumi;
    //weight *= myTree.weight_lepsf;
   
    MT2EstimateSigTree* thisEstimate = analysis->get( ht, njets, nbjets, met, minMTBmet, mt2 );
    if( thisEstimate==0 ) continue;


    if( cfg.additionalStuff()=="qgVars" ) {

      // initialize
      thisEstimate->assignVar( "qgl0", -1. );
      thisEstimate->assignVar( "qgl1", -1. );
      thisEstimate->assignVar( "qgl2", -1. );
      thisEstimate->assignVar( "qgl3", -1. );
      thisEstimate->assignVar( "partId0", 0 );
      thisEstimate->assignVar( "partId1", 0 );
      thisEstimate->assignVar( "partId2", 0 );
      thisEstimate->assignVar( "partId3", 0 );

      float qglProd = 1.;
      float qglAve = 0.;
      int denom = 0;


      if( njets>0 && fabs(myTree.jet_eta[0])<2.5 ) {

        float qgl0 = myTree.jet_qgl[0];
        thisEstimate->assignVar( "qgl0", qgl0 );
        qglProd *= qgl0;
        qglAve += qgl0;
        denom++;
        thisEstimate->assignVar( "partId0", matchPartonToJet( 0, &myTree ) );
        //thisEstimate->assignVar( "partId0", myTree.jet_mcFlavour[0] );

      }


      if( njets>1 && fabs(myTree.jet_eta[1])<2.5 ) {

        float qgl1 = myTree.jet_qgl[1];
        thisEstimate->assignVar( "qgl1", qgl1 );
        qglProd *= qgl1;
        qglAve += qgl1;
        denom++;

        thisEstimate->assignVar( "partId1", matchPartonToJet( 1, &myTree ) );
        //thisEstimate->assignVar( "partId1", myTree.jet_mcFlavour[1] );

      }
        
      if( njets>2 && fabs(myTree.jet_eta[2])<2.5 ) {

        float qgl2 = myTree.jet_qgl[2];
        thisEstimate->assignVar( "qgl2", qgl2 );
        qglProd *= qgl2;
        qglAve += qgl2;
        denom++;

        thisEstimate->assignVar( "partId2", matchPartonToJet( 2, &myTree ) );
        //thisEstimate->assignVar( "partId2", myTree.jet_mcFlavour[2] );

      }
        

      if( njets>3 && fabs(myTree.jet_eta[3])<2.5 ) {

        float qgl3 = myTree.jet_qgl[3];
        thisEstimate->assignVar( "qgl3", qgl3 );
        qglProd *= qgl3;
        qglAve += qgl3;
        denom++;

        thisEstimate->assignVar( "partId3", matchPartonToJet( 3, &myTree ) );
        //thisEstimate->assignVar( "partId3", myTree.jet_mcFlavour[3] );

      }

      qglAve /= (float)denom;
        
      thisEstimate->assignVar( "qglProd", qglProd );
      thisEstimate->assignVar( "qglAve", qglAve );

      thisEstimate->assignTree(myTree, weight );
      thisEstimate->tree->Fill();

    } else {

      thisEstimate->fillTree(myTree, weight );

    }
    
    thisEstimate->yield3d->Fill(mt2, GenSusyMScan1, GenSusyMScan2, weight );
    thisEstimate->yield->Fill(mt2, weight );
    
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

    if( EventYield[i]->id >= id_min && EventYield[i]->id <= id_max ) {

       *(return_EventYield) += *(EventYield[i]);

    }

  } // for EventYield


  return return_EventYield;

}




void drawYields( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields ) {


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


    TH1D* h1_data = data->get(thisRegion)->yield;

    TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), thisRegion.getName().c_str()), "recreate" );
    histoFile->cd();
    h1_data->Write();

    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h1_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.6);


    THStack bgStack("bgStack", "");
    for( unsigned i=0; i<bgYields.size(); ++i ) { // reverse ordered stack is prettier
      int index = bgYields.size() - i - 1;
      TH1D* h1_bg = bgYields[index]->get(thisRegion)->yield;
      h1_bg->SetFillColor( colors[index] );
      h1_bg->SetLineColor( kBlack );
      bgStack.Add(h1_bg);
    }


    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();

    float xMin = h1_data->GetXaxis()->GetXmin();
    float xMax = h1_data->GetXaxis()->GetXmax();
    float yMax1 = h1_data->GetMaximum()*1.5;
    float yMax2 = 1.5*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum()));
    float yMax3 = 1.5*(bgStack.GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    //float yMax = TMath::Max( h1_data->GetMaximum()*1.5, (h1_data->GetMaximum() + h1_data->GetBinError(h1_data->GetMaximumBin()))*1.2);
    //float yMax = h1_data->GetMaximum()*1.5;
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
    if( dummyAnalysis )
      legend->AddEntry( gr_data, "Dummy", "P" );
    else
      legend->AddEntry( gr_data, "Data", "P" );
    histoFile->cd();
    for( unsigned i=0; i<bgYields.size(); ++i ) {  
      TH1D* h1_bg = bgYields[i]->get(thisRegion)->yield;
      legend->AddEntry( h1_bg, bgYields[i]->getFullName().c_str(), "F" );
      h1_bg->Write();
    }

    histoFile->Close();

    legend->Draw("same");
    bgStack.Draw("histo same");
    gr_data->Draw("p same");

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


MT2Config::MT2Config( const std::string& configFileName ) {

  std::cout << std::endl;
  std::cout << "-> Reading config file: " << configFileName << std::endl;
  std::cout << std::endl;

  regionsSet_ = "";
  mcSamples_ = "";
  sigSamples_ = "";
  dataSamples_ = "";
  additionalStuff_ = "";

  ifstream IN(configFileName.c_str());
  char buffer[200];
  char StringValue[1000];


  while( IN.getline(buffer, 200, '\n') ) {

    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'                                                                                                                                                                                 
    }

    std::cout << buffer << std::endl;

    char name_c[200];
    sscanf(buffer, "%s %s", name_c, StringValue);
    std::string name(name_c);

    if( name=="regionsSet" )
      regionsSet_ = std::string(StringValue);
    else if( name=="mcSamples" )
      mcSamples_ = std::string(StringValue);
    else if( name=="sigSamples" )
      sigSamples_ = std::string(StringValue);
    else if( name=="dataSamples" )
      dataSamples_ = std::string(StringValue);
    else if( name=="additionalStuff" )
      additionalStuff_ = std::string(StringValue);

  } // while getline



  std::cout << std::endl;
     
}



void randomizePoisson( MT2Analysis<MT2EstimateTree>* data ) {

  TRandom3 rand(13);


//  std::set<MT2HTRegion> HTRegions = data->getHTRegions();
//  std::set<MT2SignalRegion> signalRegions = data->getSignalRegions();

  std::set<MT2Region> MT2Regions = data->getRegions();

//  for( std::set<MT2HTRegion>::iterator iHT = HTRegions.begin(); iHT!=HTRegions.end(); ++iHT ) {
//    for( std::set<MT2SignalRegion>::iterator iSR = signalRegions.begin(); iSR!=signalRegions.end(); ++iSR ) {

      for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

	MT2Region thisRegion( (*iMT2) );
	
	TH1D* h1_data = data->get(thisRegion)->yield;

	for( int ibin=1; ibin<h1_data->GetXaxis()->GetNbins()+1; ++ibin ) {

	  int poisson_data = rand.Poisson(h1_data->GetBinContent(ibin, 0, 0));
	  h1_data->SetBinContent(ibin, 0, 0, poisson_data);
	  h1_data->SetBinError(ibin, 0, 0,  0.);
	  
	}  // for bins

      }// for MT2 regions

//    }// for signal regions
//  }// for HT regions

}


int matchPartonToJet( int index, MT2Tree* myTree ) {

  float deltaRMin = 0.3; // at least 0.3
  int foundId = 0;

  TLorentzVector jet;
  jet.SetPtEtaPhiM( myTree->jet_pt[index], myTree->jet_eta[index], myTree->jet_phi[index], myTree->jet_mass[index] );


  for( int i=0; i<myTree->ngenPart; ++i ) {

    if( myTree->genPart_status[i]!=23 ) continue;
    if( !(myTree->genPart_pdgId[i]==21 || abs(myTree->genPart_pdgId[i]<6)) ) continue;

    TLorentzVector thisPart;
    thisPart.SetPtEtaPhiM( myTree->genPart_pt[i], myTree->genPart_eta[i], myTree->genPart_phi[i], myTree->genPart_mass[i] );

    float thisDeltaR = jet.DeltaR(thisPart);
    if( thisDeltaR<deltaRMin ) {
      deltaRMin = thisDeltaR;
      foundId = myTree->genPart_pdgId[i];
    }
  }

  return foundId;

}
