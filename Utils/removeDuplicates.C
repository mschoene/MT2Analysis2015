/*
How to run:
root -l
.L removeDuplicates.C+
removeDuplicates(inputfile.root, false, outputfile.root, treeName)

Note: you should run on each baby-tree with the option 'false' as 2nd argument the first time. Only if duplicates will be found, run again with option 'true' to fill a new tree without duplicates. This will save time! 
*/


#include "TTree.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TChain.h"

#include <iostream>
#include <set>
#include <ctime>

using namespace std;

class EventKey {
public:
  //  EventKey(unsigned int input_run=0, unsigned int input_lumi=0, unsigned long long input_evt=0) : 
  EventKey(int input_run=0, int input_lumi=0, unsigned long long input_evt=0) : 
    run_(input_run), lumi_(input_lumi), evt_(input_evt){;}

//  unsigned int run() const {return run_;}
//  unsigned int lumi() const {return lumi_;}
  int run() const {return run_;}
  int lumi() const {return lumi_;}
  unsigned long long evt() const {return evt_;}

  bool operator<(EventKey const& right) const{
    if (run_ == right.run()) {
      if (lumi_ == right.lumi()) {
	return evt_ < right.evt();
      }
      return lumi_ < right.lumi();
    }
    return run_ < right.run();
  }

private:
//  unsigned int run_;
//  unsigned int lumi_;
  int run_;
  int lumi_;
  unsigned long long evt_;

};



void removeDuplicates(string inputFile="duplicates.root",
		      bool fillNewTree=false,
		      string outputFile="output.root",
		      string treeName="mt2"){

  int start_s=clock();

  //Get input tree
  TFile *oldfile = new TFile(inputFile.c_str());
  TTree *oldtree = (TTree*)oldfile->Get(treeName.c_str());
  Long64_t nentries = oldtree->GetEntries();
  
  cout << "In input tree, nentries = " << nentries << endl;

  //  unsigned int run,lumi;  //CMG
  int run,lumi;         // americans
  unsigned long long evt;
  oldtree->SetBranchAddress("run",&run);
  oldtree->SetBranchAddress("lumi",&lumi);
  oldtree->SetBranchAddress("evt",&evt);
  
  //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile(outputFile.c_str(),"recreate");
  TTree *newtree = oldtree->CloneTree(0);
  

  //Create set where we store list of event keys
  std::set<EventKey> previousEvents;

  int nDuplicates = 0;

  for (Long64_t i=0;i<nentries; i++) {
    //for (Long64_t i=0;i<1000; i++) {
    oldtree->GetEntry(i);

    EventKey newEvent(run,lumi,evt);
    bool isDuplicate = !previousEvents.insert(newEvent).second;

    if(i%100000==0) {
      time_t t = time(0);   // get time now
      tm * now = localtime( & t );
      cout << "Processing event: " << i << " at time " 
	   << now->tm_hour << ":"
	   << now->tm_min << ":"
	   << now->tm_sec 
	   << endl;
    }
    
    if(!isDuplicate) {
      if(fillNewTree) newtree->Fill();
    }else{
      nDuplicates++;
      //cout << "Found duplicate! run,lumi,evt: " 
      //     << run << " , " << lumi << " , " << evt <<endl;
    }

  }

  cout << "Number of duplicates found: " << nDuplicates << endl;


  //newtree->Print();
  newtree->AutoSave();
  delete oldfile;
  delete newfile;

  int stop_s=clock();
  cout << "elapsed time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << endl;
}

int removeDuplicatesFromChain(string inputFilesList,
			      bool fillNewTree=false,
			      string outputFile="output.root",
			      string treeName="mt2"){
  int start_s=clock();

  // -------- Same code for TChain on remote files as in postProcessing.C
  TChain* chain = new TChain(treeName.c_str());

  //std::cout << inputFilesList.c_str() << std::endl;
  TFileCollection *filelist = new TFileCollection("listOfFiles");
  
  std::cout << "Adding files list... " << std::endl;
  filelist->AddFromFile( inputFilesList.c_str() );
  
  std::cout << "Adding files to chain... " << std::endl;  
  int chainReturn = chain->AddFileInfoList((TCollection*) filelist->GetList());
 
  if (chainReturn < 1) {
    std::cout << "ERROR: input folder/fileName is not well defined. Exit!" << std::endl;
    std::cout << "InputFileList: " << inputFilesList << std::endl;
    return 1;
  }
  // ------------------------------------------------------

  cout << "Chain created" << endl;
  Long64_t nentries = chain->GetEntries();
  
  cout << "In input tree, nentries = " << nentries << endl;

  //  unsigned int run,lumi;  //CMG
  int run,lumi;         // americans
  unsigned long long evt;
  chain->SetBranchAddress("run",&run);
  chain->SetBranchAddress("lumi",&lumi);
  chain->SetBranchAddress("evt",&evt);
  
  //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile(outputFile.c_str(),"recreate");
  TTree *newtree = chain->CloneTree(0);

  cout << "Tree cloned" << endl;

  //Create set where we store list of event keys
  std::set<EventKey> previousEvents;

  int nDuplicates = 0;

  for (Long64_t i=0;i<nentries; i++) {
    chain->GetEntry(i);

    EventKey newEvent(run,lumi,evt);
    bool isDuplicate = !previousEvents.insert(newEvent).second;

    if(i%100000==0) {
      time_t t = time(0);   // get time now
      tm * now = localtime( & t );
      cout << "Processing event: " << i << " at time " 
	   << now->tm_hour << ":"
	   << now->tm_min << ":"
	   << now->tm_sec 
	   << endl;
    }
    
    if(!isDuplicate) {
      if(fillNewTree) newtree->Fill();
    }else{
      nDuplicates++;
      //cout << "Found duplicate! run,lumi,evt: " 
      //     << run << " , " << lumi << " , " << evt <<endl;
    }

  }

  cout << "Number of duplicates found: " << nDuplicates << endl;


  //newtree->Print();
  newtree->AutoSave();
  delete chain;
  delete newfile;

  int stop_s=clock();
  cout << "elapsed time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << endl;

  return 0;
}

