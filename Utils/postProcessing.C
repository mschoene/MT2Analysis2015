/*
How to run:
root -l
.L postProcessing.C+
run()
*/

#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TFile.h"
#include "TFileMerger.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TBranch.h"
#include "TString.h"
#include "TH1D.h"

using namespace std;


int postProcessing(string inputString="input",
		   string inputFolder="./",
		   string outputFile="output.root",
		   string treeName="tree",
		   float filter=1.0, float kfactor=1.0, float xsec=-1.0, int id=1,
		   string crabExt="");


int run(string cfg="postProcessing_74X_50ns.cfg",
	string treeName="tree", 
	string inputFolder = "/pnfs/psi.ch/cms/trivcat/store/user/casal/babies/PHYS14_Production_QCDpt_noSietaieta/", 
	string outputFolder = "./test/",  
	string fileExtension = "_post.root",
        string crabExt = ""){
  
  // for measuring timing
  time_t start = time(0);

  cout<<"Configuration file is: "<<cfg.c_str()<<endl;
  
  ifstream configuration(cfg.c_str());
  string line;
  
  while(std::getline(configuration,line)){
    istringstream ss(line);
    if((ss.str()[0])==(string)"#") continue;
    if(ss.str().empty()) continue;


    string idS,name,xsecS,filterS,kfactorS;
    ss >> idS;
    ss >> name;
    ss >> xsecS;
    ss >> filterS;
    ss >> kfactorS;

    int id;
    float filter,kfactor,xsec;
    filter = (float)atof(filterS.c_str());
    kfactor = (float)atof(kfactorS.c_str());
    xsec = (float)atof(xsecS.c_str());
    id = (int)atof(idS.c_str());

    int debug=0;
    if(debug){
      cout << "id,name,x,f,k: " 
	   << id << " , " 
	   << name << " , "
	   << xsec << " , " 
	   << filter << " , " 
	   << kfactor << endl;
    }  

    string outputFile = outputFolder + "/" + name + fileExtension;
    postProcessing(name, inputFolder, outputFile, treeName, filter, kfactor, xsec, id, crabExt);
  }
  
  // for printing measured timing
  time_t stop = time(0);
  double time = difftime(stop, start);
  cout << "real time in seconds: " << time << endl;
  
  return 0;
}

int postProcessing(string inputString,
		   string inputFolder,
		   string outputFile,
		   string treeName,
		   float filter, float kfactor, float xsec, int id,
		   string crabExt)
{
  TChain* chain = new TChain(treeName.c_str());
  // Add all files in the input folder
  string dcap = inputFolder.find("pnfs")!=std::string::npos ? "dcap://t3se01.psi.ch:22125/" : "";
  string fullInputString;
  if(crabExt!="")
    fullInputString = dcap + inputFolder + "/" + inputString + "/"+ crabExt +"/0000/mt2*.root";
  else
    fullInputString = dcap + inputFolder + "/" + inputString + "/mt2*.root";
  //string fullInputString = dcap + inputFolder + "/" + inputString + "/mt2_14.root";
  int chainReturn = chain->Add( fullInputString.c_str() );
  if (chainReturn < 1) {
    cout << "ERROR: input folder/fileName is not well defined. Exit!" << endl;
    cout << "fullInputString: " << fullInputString << endl;
    return 1;
  }
  
  // here I set the "Count" histograms
  TIter nextfile(chain->GetListOfFiles());
  TChainElement *elem;
  bool isFirst=true;
  TH1D* newH=0;
  unsigned long int allHistoEntries=0;
  TH1D* newSumW=0;
  unsigned long int allSumGenWeight=0;
  while ((elem = (TChainElement*)nextfile())) {
    TFile *f; f = TFile::Open(elem->GetTitle(),"READ");
    TH1D *countH = (TH1D*)f->Get("Count");
    TH1D *sumW = (TH1D*)f->Get("SumGenWeights");
    if(isFirst){
      newH = (TH1D*) countH->Clone();
      newH->SetDirectory(0);  
      isFirst=false;
      
      newSumW = (TH1D*) sumW->Clone();
      newSumW->SetDirectory(0);  
      isFirst=false;      
    }
    allHistoEntries += countH->GetBinContent(1);
    allSumGenWeight += sumW->GetBinContent(1);
    f->Close();
    delete f;      
  }
  newH->SetBinContent(1,allHistoEntries);
  newSumW->SetBinContent(1,allSumGenWeight);

  // This line should be uncommented for all the branches that we want to overwrite.
  // If the branch is not in the input tree, we don't need this.
  //
  //t->SetBranchStatus("scale1fb", 0);

  TFile *puIN = TFile::Open("/shome/mmasciov/PUProfile_JetHT.root");
  TH1D* hPU_data = (TH1D*) puIN->Get("hPU_data");
  
  TFile *out = TFile::Open(outputFile.c_str(), "RECREATE");
  TTree *clone = new TTree("mt2", "post processed baby tree for mt2 analysis");

  clone = chain->CloneTree(-1, "fast"); 
  clone->SetName("mt2");
  
  //if(SortBasketsByEntry)
  // clone = t->CloneTree(-1, "fastSortBasketsByEntry");
  //else 
  // clone = t->CloneTree(-1, "fast");
  
  //-------------------------------------------------------------

  //Calculate scaling factor and put variables into tree 
  ULong64_t nEventsTree = clone->GetEntries();
  ULong64_t nEventsHisto = (ULong64_t) newH->GetBinContent(1); 					 
//  Int_t nEventsTree = clone->GetEntries();
//  Int_t nEventsHisto = (Int_t) newH->GetBinContent(1); 					 

  float genWeight_=1.0;
  float genWeight=0.;
  chain->SetBranchAddress("genWeight", &genWeight);
  
  int isData=0; 
  chain->SetBranchAddress("isData",&isData);
  
  int nVert=0;
  chain->SetBranchAddress("nVert", &nVert);

  if( nEventsTree > 0 ){

    chain->GetEntry(0);

    if(isData)
      genWeight_ = 0.;
    else
      genWeight_ = fabs(genWeight);

  }



  TH1D* hPU = new TH1D("hPU", "", 100, 0, 100);
  hPU->Sumw2();
  
  chain->Project("hPU", "nVert");
  hPU->Scale(1.0/nEventsTree);
  
  TH1D* hPU_r = (TH1D*) hPU_data->Clone("hPU_r");
  hPU_r->Divide(hPU);


  //  float scale1fb = xsec*kfactor*1000*filter/(Float_t)nEventsHisto;
  ULong64_t sumGenWeightsHisto; 
  ULong64_t nEffEventsHisto; 
  float scale1fb_noGenWeight;
  float scale1fb_sumGenWeights;
  float scale1fb;
  float puWeight;

  if(isData){
    
    sumGenWeightsHisto = (ULong64_t) 0;
    nEffEventsHisto = nEventsHisto;
    
    scale1fb_noGenWeight = 0.0;
    scale1fb_sumGenWeights = 0.0;
    
  }
  else{
    
    sumGenWeightsHisto = (ULong64_t)  newSumW->GetBinContent(1);
    nEffEventsHisto = ( (double) 1.0*sumGenWeightsHisto/genWeight_  > (ULong64_t) (1.0*sumGenWeightsHisto/genWeight_ + 0.5) ) ? (ULong64_t) (1.0*sumGenWeightsHisto/genWeight_)+1 : (ULong64_t) (1.0*sumGenWeightsHisto/genWeight_);
  
    scale1fb_noGenWeight = xsec*kfactor*1000*filter/(Float_t)nEffEventsHisto;
    scale1fb_sumGenWeights = xsec*kfactor*1000*filter/(Float_t)sumGenWeightsHisto;

  }
  
  if (nEventsHisto < nEventsTree) // this should not happen
    cout << "ERROR: histogram count has less events than tree. This indicates something went wrong" << endl
	 << "#events histo: "  << nEventsHisto << endl
	 << "#events tree: "  << nEventsTree << endl;
  else if (nEventsHisto < nEventsTree) // this can happen
    cout << "WARNING: histogram count has more events than tree. This should only happen if tree was skimmed" << endl
	 << "#events histo: "  << nEventsHisto << endl
	 << "#events tree: "  << nEventsTree << endl;
    
  TBranch* b1 = clone->Branch("evt_scale1fb", &scale1fb, "evt_scale1fb/F");
  TBranch* b2 = clone->Branch("evt_scale1fb_noGenWeight", &scale1fb_noGenWeight, "evt_scale1fb_noGenWeights/F");
  TBranch* b3 = clone->Branch("evt_scale1fb_sumGenWeights", &scale1fb_sumGenWeights, "evt_scale1fb_sumGenWeights/F");
  TBranch* b4 = clone->Branch("evt_xsec", &xsec, "evt_xsec/F");
  TBranch* b5 = clone->Branch("evt_kfactor", &kfactor, "evt_kfactor/F");
  TBranch* b6 = clone->Branch("evt_filter", &filter, "evt_filter/F");
  TBranch* b7 = clone->Branch("evt_nEvts", &nEventsHisto, "evt_nEvts/l");
  TBranch* b8 = clone->Branch("evt_nEffectiveEvts", &nEffEventsHisto, "evt_nEffectiveEvts/l");
  TBranch* b9 = clone->Branch("evt_sumGenWeights", &sumGenWeightsHisto, "evt_sumGenWeights/l");
  TBranch* b10 = clone->Branch("evt_id", &id, "evt_id/I");
  TBranch* b11 = clone->Branch("evt_puWeight", &puWeight, "evt_puWeight/F");

  for(Long64_t i = 0; i < (Long64_t) nEventsTree; i++) {
    
    chain->GetEntry(i);
    if(isData){
      //data should never be rescaled with scale1fb when making plots
      //set scaler to zero so that user cannot miss the mistake
      scale1fb = 0.0; 
      puWeight = 1.0;
    }
    else{
      scale1fb= xsec*kfactor*1000*filter*genWeight/(Float_t)sumGenWeightsHisto;
      
      int puBin = (int) hPU_r->GetXaxis()->FindBin(nVert);
      puWeight = hPU_r->GetBinContent(puBin);

    }
    
    b1->Fill();
    b2->Fill();
    b3->Fill();
    b4->Fill();
    b5->Fill();
    b6->Fill();
    b7->Fill();
    b8->Fill();
    b9->Fill();
    b10->Fill();
    b11->Fill();
    
  }
  //-------------------------------------------------------------

  delete chain; 


  hPU_r->Write();
  delete hPU_r;
  delete hPU;
  delete hPU_data;
  newH->Write();
  delete newH;
  newSumW->Write();
  delete newSumW;
  clone->Write();
  delete clone;
  out->Close();
  return 0;
  
}
