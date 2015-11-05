/*
How to run:
root -l
.L filterFromTxt.C
filterFromTxt(filterList.txt, inputfile.root, false, outputfile.root, treeName, flagName)

Note: if flagname is specified, a branch with flag information is added, otherwise events are filtered

Note: you should run on each baby-tree with the option 'false' as 3nd argument the first time. Only if duplicates will be found, run again with option 'true' to fill a new tree without duplicates. This will save time! 
*/


#include "TTree.h"
#include "TFile.h"

#include <fstream>
#include <iostream>
#include <set>
#include <ctime>
#include <string>

using namespace std;

class EventKey {
public:
  EventKey(unsigned int input_run=0, unsigned int input_lumi=0, unsigned long long input_evt=0) : 
    run_(input_run), lumi_(input_lumi), evt_(input_evt){;}

  unsigned int  run () const {return run_; }
  unsigned int  lumi() const {return lumi_;}
  unsigned long long evt () const {return evt_; }

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
  unsigned int run_;
  unsigned int lumi_;
  unsigned long long evt_;

};

// overload >> operator
std::ifstream& operator>>(std::ifstream& is, EventKey& event) {
  
    std::string run, lumi, evt;
    std::getline(is, run , ':');
    std::getline(is, lumi, ':');
    std::getline(is, evt , '\n');
    if (!run.empty() && !lumi.empty() && !evt.empty())
      event = EventKey(stoi(run),stoi(lumi),stol(evt));
  return is;
  
};

void filterFromTxt(string filterList="eventlist_JetHT_csc2015.txt",
		   string inputFile="duplicates.root",
		   bool fillNewTree=false,
		   string outputFile="output.root",
		   string treeName="mt2",
		   string flagName=""){

  int start_s=clock();

  //Get input tree
  TFile *oldfile = new TFile(inputFile.c_str());
  TTree *oldtree = (TTree*)oldfile->Get(treeName.c_str());
  Long64_t nentries = oldtree->GetEntries();
  
  cout << "In input tree, nentries = " << nentries << endl;

  unsigned int run,lumi;
  unsigned long long evt;
  oldtree->SetBranchAddress("run" ,&run );
  oldtree->SetBranchAddress("lumi",&lumi);
  oldtree->SetBranchAddress("evt" ,&evt );
  
  //Create a new file + a clone of old tree in new file
  TFile *newfile;
  TTree *newtree;
  float flag=1.0; // our current Flag are float, keep it float for consistency
  if ( fillNewTree ) {
    newfile = new TFile(outputFile.c_str(),"recreate");
    newtree = oldtree->CloneTree(0);
    if ( flagName != "" )
      TBranch* b_flag = newtree->Branch(flagName.c_str(), &flag, (flagName+"/F").c_str());
  }

  //Create set where we store list of event keys
  std::set<EventKey> listFromTxt;

  ifstream infile(filterList.c_str());
  while (!infile.eof()) {
    EventKey aEvent;
    infile >> aEvent;
    listFromTxt.insert( aEvent );
  }

  int nRemoved=0;

  for (Long64_t i=0;i<nentries; i++) {
    //for (Long64_t i=0;i<1000; i++) {
    oldtree->GetEntry(i);

    EventKey newEvent(run,lumi,evt);
    bool inList = listFromTxt.find(newEvent)!=listFromTxt.end();

    if(i%100000==0) {
      time_t t = time(0);   // get time now
      tm * now = localtime( & t );
      cout << "Processing event: " << i << " at time " 
	   << now->tm_hour << ":"
	   << now->tm_min << ":"
	   << now->tm_sec 
	   << endl;
    }
    

    if ( fillNewTree && flagName != "" ){
      flag = inList ? 0.0 : 1.0;
      newtree->Fill();
      
    }else if ( !inList ) {
      if(fillNewTree) newtree->Fill();
    }

    if ( inList ){
      nRemoved++;
      // cout << "Found event to filter! run,lumi,evt: " 
      // 	   << run << " , " << lumi << " , " << evt <<endl;
    }

  }


  cout << "Number of events found to be contained in the list = " << nRemoved << endl;

  //newtree->Print();
  if(fillNewTree) {
    newtree->AutoSave();
    delete newfile;
  }
  delete oldfile;

  int stop_s=clock();
  cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;

}

