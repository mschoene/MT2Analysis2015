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
		   float filter=1.0, float kfactor=1.0, float xsec=-1.0, int id=1 );


int run(string cfg="postProcessing.cfg",
	string treeName="tree", 
	string inputFolder = "/pnfs/psi.ch/cms/trivcat/store/user/casal/babies/PHYS14_Production_QCDpt_noSietaieta/", 
	string outputFolder = "./test/",  
	string fileExtension = "_post.root"){
  
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
    postProcessing(name, inputFolder, outputFile, treeName, filter, kfactor, xsec, id);
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
		   float filter, float kfactor, float xsec, int id)
{
  TChain* chain = new TChain(treeName.c_str());
  // Add all files in the input folder 
  string dcap = inputFolder.find("pnfs")!=std::string::npos ? "dcap://t3se01.psi.ch:22125/" : "";
  string fullInputString = dcap + inputFolder + "/" + inputString + "/*.root";
  //string fullInputString = dcap + inputFolder + "/" + inputString + "/mt2_14.root";
  int chainReturn = chain->Add(fullInputString.c_str()  );
  if (chainReturn < 1) {
    cout << "ERROR: input folder/fileName is not well defined. Exit!" << endl;
    cout << "fullInputString: " << fullInputString << endl;
    return 1;
  }
  

  // Merge (add) Count histograms from tchain files
  if(chainReturn>1){
    TFileMerger merger = TFileMerger(kFALSE);
    TObjArray *fileElements = chain->GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl=0;
    while (( chEl=(TChainElement*)next() )){    
      merger.AddFile(chEl->GetTitle());
    }
    merger.SetNotrees(kTRUE);
    merger.OutputFile(outputFile.c_str(), "RECREATE");
    bool returnCode=merger.Merge();
    if(!returnCode) {
      cout << "ERROR: merger.Merge() returned " << returnCode << endl;
      return 1;
    }
  }else{
    // have to do this because the merger doesn't work as expected when the input is one single file
    TFile* fileIn = TFile::Open(fullInputString.c_str());    
    TH1* histo = (TH1D*)fileIn->Get("Count");
    TFile output(outputFile.c_str(), "RECREATE");
    histo->Write();
    fileIn->Close();
    delete fileIn;
  }


  // This line should be uncommented for all the branches that we want to overwrite.
  // If the branch is not in the input tree, we don't need this.
  //
  //t->SetBranchStatus("scale1fb", 0);


  
  TFile *out = TFile::Open(outputFile.c_str(), "UPDATE");
  TTree *clone = new TTree("mt2", "post processed baby tree for mt2 analysis");

  clone = chain->CloneTree(-1, "fast"); 
  clone->SetName("mt2");

  
  //if(SortBasketsByEntry)
  // clone = t->CloneTree(-1, "fastSortBasketsByEntry");
  //else 
  // clone = t->CloneTree(-1, "fast");
  

  //-------------------------------------------------------------

  //Calculate scaling factor and put variables into tree 
  int nEventsTree = clone->GetEntries();
  int nEventsHisto = (int)((TH1D*)out->Get("Count"))->GetBinContent(1);
  float scale1fb = xsec*kfactor*1000*filter/(Float_t)nEventsHisto;
 
  if (nEventsHisto < nEventsTree) // this should not happen
    cout << "ERROR: histogram count has less events than tree. This indicates something went wrong" << endl
	 << "#events histo: "  << nEventsHisto << endl
	 << "#events tree: "  << nEventsTree << endl;
  else if (nEventsHisto < nEventsTree) // this can happen
    cout << "WARNING: histogram count has more events than tree. This should only happen if tree was skimmed" << endl
	 << "#events histo: "  << nEventsHisto << endl
	 << "#events tree: "  << nEventsTree << endl;
    

  int isData; 
  clone->SetBranchAddress("isData",&isData);
  clone->GetEntry(0); 
  if(isData){
    //data should never be rescaled with scale1fb when making plots
    //set scaler to zero so that user cannot miss the mistake
    scale1fb = 0.0; 
  }
  
  TBranch* b1 = clone->Branch("evt_scale1fb", &scale1fb, "evt_scale1fb/F");
  TBranch* b2 = clone->Branch("evt_xsec", &xsec, "evt_xsec/F");
  TBranch* b3 = clone->Branch("evt_kfactor", &kfactor, "evt_kfactor/F");
  TBranch* b4 = clone->Branch("evt_filter", &filter, "evt_filter/F");
  TBranch* b5 = clone->Branch("evt_nEvts", &nEventsHisto, "evt_nEvts/I");
  TBranch* b6 = clone->Branch("evt_id", &id, "evt_id/I");

  for(Int_t i = 0; i < nEventsTree; i++) {
    b1->Fill();
    b2->Fill();
    b3->Fill();
    b4->Fill();
    b5->Fill();
    b6->Fill();
  }
  //-------------------------------------------------------------

  delete chain; 

  clone->Write();
  delete clone;
  out->Close();
  return 0;
  
}
