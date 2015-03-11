/*
How to run:
root -l
.L postProcessing.C+
run("samples_50ns_miniaod.txt", "treeProducerSusyFullHad", "/scratch/mmasciov/ToPostProcess/", "/scratch/mmasciov/PostProcessED/", "_babytree.root")
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
#include "TH1F.h"

using namespace std;


int postProcessing(string inputFile="input.root",
		   string outputFile="output.root",
		   string treeName="tree",
		   float filter=1.0, float kfactor=1.0, float xsec=-1.0, int id=1 );


int run(string sampleFileName="samples_50ns_miniaod.txt",
	string treeName="tree", 
	string inputPath = "/scratch/mmasciov/CSA14v2/ToPostProcess/", 
	string outputPath = "/scratch/mmasciov/CSA14v2/PostProcess/",  
	string fileExtension = "_post.root"){
  
//int run(string sampleFileName="samples_50ns_miniaod.txt",
//	string treeName="treeProducerSusyFullHad",
//	string inputPath = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/",
//	string outputPath = "/scratch/mmasciov/12Dec/PostProcess/",
//	string fileExtension = "_PU_S14_POSTLS170_babytree_post.root"){
//  
  int start_s=clock();

  ifstream IN(sampleFileName.c_str());

  if(!IN)
    system(Form("wget https://mangano.web.cern.ch/mangano/public/MECCA/%s", sampleFileName.c_str()));

  char buffer[512];

  std::cout << "------------------------------------" << std::endl;
  std::cout << "Sample File  " << sampleFileName << std::endl;

  int counter(0);

  while( IN.getline(buffer, 512, '\n') ) {

    if (buffer[0] == '#' || buffer[0] == '\n' || buffer[0] == ' ' || buffer[0] == '\0' || strlen(buffer)==1)
      continue;

    ++counter;
    std::cout << "# sample " << counter << std::endl;
    std::cout << buffer << std::endl;

    int id;
    char datasetName[512];
    float xs, filter, kfactor;
    sscanf(buffer, "%d\t/%s\t%f\t%f\t%f\n", &id, datasetName, &xs, &filter, &kfactor);
    std::cout << id << "\t" << xs << "\t" << filter << "\t" << "\t" << kfactor << std::endl;

    std::string dataset(datasetName);
    std::size_t length = dataset.find("_Tune");
    if (length==std::string::npos) // if "_Tune" doesn't exist try "_13TeV" (some QCD samples)
      length = dataset.find("_13TeV");

    char fileName[512];
    std::size_t nameLength = dataset.copy(fileName, length, 0);
    fileName[nameLength]='\0';

    for (int c = 0; c <= nameLength; ++c) {
      
      if (fileName[c] == '-' || fileName[c] == '_')
	fileName[c] = '?';
      
      if (fileName[c]=='\0')
	break;
      
    }

    std::string file(fileName);
    TString wildcard;
    if (inputPath.find("/pnfs")==std::string::npos)
      wildcard = TString::Format("%s%s*.root", inputPath.c_str(), file.c_str());
    else // mind the directory structure in the SE 
      wildcard = TString::Format("%s%s*/*.root", inputPath.c_str(), file.c_str());

    system(Form("ls %s > inputFILE.txt", wildcard.Data()));
      
    
    ifstream inputFILE("inputFILE.txt");
    if(!inputFILE) {

      std::cout << "ATTENTION: no input file, looking for next" << std::endl;
      continue;
    
    }

    char inputBuffer[512];

    inputFILE.getline(inputBuffer, 512, '\n'); 
      
    if (inputBuffer[0] == '\0' || strlen(inputBuffer)==1){
      
      std::cout << "ATTENTION: no input file, looking for next" << std::endl;
      continue;
    
    }

    std::string inputFile(inputBuffer);

    char unused[512];
    inputFILE.getline(unused, 512, '\n'); 
    if (unused[0] == '\0' || strlen(unused)==1){      
      std::cout << "single input file: " << inputFile << std::endl;
    }
    else{ // in case of more than one file use wildcard
      if (inputPath.find("/pnfs")==std::string::npos) // if not in SE use directly wildcar
	inputFile = wildcard.Data();
      else {// double wildcard does not work in SE. dataset name in directory not in filename
	std::size_t len = inputFile.rfind("/");
	inputFile = inputFile.substr(0,len)+"/*.root";
      }
      std::cout << "multiple input files: " << inputFile << std::endl;
    }

    std::string outputFileName(inputBuffer);    
    std::size_t outlength = outputFileName.find(inputPath);
    outputFileName.erase(outlength, inputPath.length());
    std::string extension = ".root";
    std::size_t extensionlength = outputFileName.find(".root");
    outputFileName.erase(extensionlength, extension.length());
    
    if (inputPath.find("/pnfs")!=std::string::npos) // in SE output name is input dir
      outputFileName.erase(outputFileName.find("/"), outputFileName.length());

    std::string outputFile = outputPath + outputFileName + fileExtension;
    std::cout << inputFile << std::endl;
    std::cout << outputFile << std::endl;

    FILE* rootfile;
    rootfile = fopen(outputFile.c_str(), "r");
    if (rootfile == NULL);
    else {

      std::cout << "ATTENTION: Target file already exists!!! I will skip this sample..." << std::endl;
      continue;

    }

    postProcessing(inputFile, outputFile, treeName, filter, kfactor, xs, id);

  }
  
  int stop_s=clock();
  std::cout << std::endl << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

  system(Form("rm -f inputFILE.txt"));
  
  return 0;
}

int postProcessing(string inputFile,
		   string outputFile,
		   string treeName,
		   float filter, float kfactor, float xsec, int id){

  //TString infname, TString outfile, unsigned int events, 
  //Float_t xsec, Float_t kfactor,
  //Float_t filt_eff, bool SortBasketsByEntry ) {
  
  cout << "Processing File " << inputFile << endl;

  // TFile *f = TFile::Open(inputFile.c_str(), "READ");
  // if (! f || f->IsZombie()) {
  //   cout << "File does not exist!" << endl;
  //   return 1;
  // }
  
  // TTree* t = (TTree*)f->Get(treeName.c_str());
  // if (! t || t->IsZombie()) {
  //   cout << "Tree does not exist!" << endl;
  //   return 2;
  // }

  TChain *c = new TChain(treeName.c_str());
  TString dcap = inputFile.find("pnfs")!=std::string::npos ? "dcap://t3se01.psi.ch:22125/" : "";
  int chainReturn = c->Add(dcap + inputFile.c_str());
  if (chainReturn < 1) {
    cout << "File does not exist!" << endl;
    return 1;
  }

  // merge (add) Count histograms from tchain files
  TFileMerger *merger = new TFileMerger(kFALSE);
  TObjArray *fileElements = c->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() ))
    merger->AddFile(chEl->GetTitle());
  merger->OutputFile(outputFile.c_str(), "RECREATE");
  merger->SetNotrees(kTRUE);
  merger->Merge();
  delete merger;
  
  // This line should be uncommented for all the branches that we want to overwrite.
  // If the branch is not in the input tree, we don't need this.
  //
  //t->SetBranchStatus("scale1fb", 0);

  TFile *out = TFile::Open(outputFile.c_str(), "UPDATE");
  TTree *clone = new TTree("mt2", "post processed baby tree for mt2 analysis");

  //clone = t->CloneTree(-1, "fast");
  clone = c->CloneTree(-1, "fast");
  clone->SetName("mt2");

  /*
  if(SortBasketsByEntry)
    clone = t->CloneTree(-1, "fastSortBasketsByEntry");
  else 
    clone = t->CloneTree(-1, "fast");
  */

  //-------------------------------------------------------------

  //Calculate scaling factor and put variables into tree 
  //int events = t->GetEntries();
  int events = clone->GetEntries();
  // get Nevents from Count histogram
  int count = (int)((TH1F*)out->Get("Count"))->GetBinContent(1);
  float scale1fb = xsec*kfactor*1000*filter/(Float_t)count;
 
  if (count < events) // this should not happen
    cout << "ERROR: histogram count has less events than tree" << endl
	 << "count: "  << count << endl
	 << "events: "  << events << endl;
  else if (count > events) // this can happen
    cout << "WARNING: histogram count has more events than tree" << endl
	 << "count: "  << count << endl
	 << "events: "  << events << endl;
    
  /*
  if(isdata){
	scale1fb = 1.0;
	cout<< "Data file. scale1fb: " << scale1fb << endl;
	
  }else{
	cout << "scale1fb: " << scale1fb << endl; 
  }
  */
  
  TBranch* b1 = clone->Branch("evt_scale1fb", &scale1fb, "evt_scale1fb/F");
  TBranch* b2 = clone->Branch("evt_xsec", &xsec, "evt_xsec/F");
  TBranch* b3 = clone->Branch("evt_kfactor", &kfactor, "evt_kfactor/F");
  TBranch* b4 = clone->Branch("evt_filter", &filter, "evt_filter/F");
  TBranch* b5 = clone->Branch("evt_nEvts", &count, "evt_nEvts/I");
  TBranch* b6 = clone->Branch("evt_id", &id, "evt_id/I");

  for(Int_t i = 0; i < events; i++) {
    b1->Fill();
    b2->Fill();
    b3->Fill();
    b4->Fill();
    b5->Fill();
    b6->Fill();
  }
  //-------------------------------------------------------------

  // delete t;
  // f->Close();
  delete c;

  //out->cd();
  clone->Write();
  delete clone;
  out->Close();

  
  return 0;
  
}
