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
#include <vector>

#include <string>
#include <iostream>
#include "TFile.h"
#include "TFileMerger.h"
#include "TFileCollection.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TBranch.h"
#include "TString.h"
#include "TH1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH2.h"
#include "TLorentzVector.h"



using namespace std;

int preProcessing(std::string inputString="input",
		  std::string inputFolder="./",
		  std::string outputFile="output.txt",
		  std::string treeName="tree",
		  int id=1,
		  std::string txtFileList="");


int preProcessing(std::string inputString,
		  std::string inputFolder,
		  std::string outputFile,
		  std::string treeName,
		  int id,
		  std::string txtFileList){
  

  TChain* chain = new TChain(treeName.c_str());
  // Add all files in the input folder

  std::cout << txtFileList << std::endl;

  TFileCollection *filelist = new TFileCollection("listOfFiles");
  
  std::cout << "Adding file list... " << std::endl;
  filelist->AddFromFile( txtFileList.c_str() );
  
  std::cout << "Adding files to chain... " << std::endl;  
  int chainReturn = chain->AddFileInfoList((TCollection*) filelist->GetList());
   
  if (chainReturn < 1) {
    std::cout << "ERROR: input folder/fileName is not well defined. Exit!" << std::endl;
    // std::cout << "fullInputString: " << fullInputString << std::endl;
    return 1;
  }


  // TFile *out = TFile::Open(outputFile.c_str(), "RECREATE");

  // here I set the "Count" histograms
  TIter nextfile(chain->GetListOfFiles());
  TChainElement *elem;
  bool isFirst=true;
  TH1D* newH=0;
  unsigned long int allHistoEntries=0;
  TH1D* newSumW=0;
  unsigned long int allSumGenWeight=0;
  
  std::cout << "Entering loop over chain elements..." << std::endl;
  
  int nFiles=0;

  while ((elem = (TChainElement*)nextfile())) {
    TFile *f; f = TFile::Open(elem->GetTitle(),"READ");
    TH1D *countH = (TH1D*)f->Get("Count");
    
    // ++nFiles;
    // std::cout << elem->GetTitle() << std::endl;
    // std::cout << "Read Count histogram for file "<< nFiles <<": " << countH->GetEntries() << std::endl;

    TH1D *sumW = (TH1D*)f->Get("SumGenWeights");
    //
    //std::cout << "Read SumGenWeights histogram for file "<< nFiles <<": " << sumW->GetEntries() << std::endl;
    //

    if(isFirst){
      newH = (TH1D*) countH->Clone();
      newH->SetDirectory(0);  
      
      if(sumW)
	newSumW = (TH1D*) sumW->Clone();
      else
	newSumW = (TH1D*) countH->Clone();
      newSumW->SetDirectory(0);  
      
      isFirst=false;
      
    }
    allHistoEntries += countH->GetBinContent(1);
    if(sumW)
      allSumGenWeight += sumW->GetBinContent(1);
    else
      allSumGenWeight += countH->GetBinContent(1);
    f->Close();
    delete f;      
  }
  newH->SetBinContent(1,allHistoEntries);
  newSumW->SetBinContent(1,allSumGenWeight);
  
  std::cout << "Read gen weights and count histogram..." << std::endl;


  //-------------------------------------------------------------

  //Calculate scaling factor and put variables into tree 
  ULong64_t nEventsTree = chain->GetEntries();
  ULong64_t nEventsHisto = (ULong64_t)newH->GetBinContent(1);		 


  float genWeight_=1.0;
  float genWeight=0.;
  TBranch* thisGenWeight = (TBranch*) chain->GetListOfBranches()->FindObject("genWeight");
  if (thisGenWeight) chain->SetBranchAddress("genWeight", &genWeight);


  int isData=0; 
  chain->SetBranchAddress("isData",&isData);
  chain->GetEntry(0);

  if(isData) return 1;

  UInt_t run=0;
  chain->SetBranchAddress("run", &run);
  UInt_t lumi=0;
  chain->SetBranchAddress("lumi", &lumi);

  int nVert=0;
  chain->SetBranchAddress("nVert", &nVert);



  ///float nTrueInt=0;
  int nTrueInt=0;
  TBranch* thisNTrueInt = (TBranch*) chain->GetListOfBranches()->FindObject("nTrueInt");
  if (thisNTrueInt) chain->SetBranchAddress("nTrueInt", &nTrueInt);
  
  if( nEventsTree > 0 )
    genWeight_ = fabs(genWeight);


  

  //  float scale1fb = xsec*kfactor*1000*filter/(Float_t)nEventsHisto;
  ULong64_t sumGenWeightsHisto; 
  ULong64_t nEffEventsHisto; 


  float weight_toppt;


  ////// isr re-weight
  Int_t nGenPart;
  Float_t GenPart_pt[100];
  Int_t GenPart_status[100];
  Int_t GenPart_pdgId[100];

  if (!isData) {
    std::cout << "Loading the generator parts " << std::endl;
    chain->SetBranchAddress("nGenPart", &nGenPart);
    chain->SetBranchAddress("GenPart_pt", GenPart_pt);
    chain->SetBranchAddress("GenPart_status", GenPart_status);
    chain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
  }

 


 
  sumGenWeightsHisto = (ULong64_t) newSumW->GetBinContent(1);
  nEffEventsHisto = ( (double) 1.0*sumGenWeightsHisto/genWeight_  > (ULong64_t) (1.0*sumGenWeightsHisto/genWeight_ + 0.5) ) ? (ULong64_t) (1.0*sumGenWeightsHisto/genWeight_ + 1.0) : (ULong64_t) (1.0*sumGenWeightsHisto/genWeight_);
  
 
  if (nEventsHisto < nEventsTree) // this should not happen
    std::cout << "ERROR: histogram count has less events than tree. This indicates something went wrong" << std::endl
	      << "#events histo: "  << nEventsHisto << std::endl
	      << "#events tree: "  << nEventsTree << std::endl;
  else if (nEventsHisto > nEventsTree) // this can happen
    std::cout << "WARNING: histogram count has more events than tree. This should only happen if tree was skimmed" << std::endl
	      << "#events histo: "  << nEventsHisto << std::endl
	      << "#events tree: "  << nEventsTree << std::endl;
    

  double average = 0;
  average = 1.0; //dont' average for now
  /*
  //Top pt reweighting 
  //Values from Run1 still valid for now, see here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  float a = 0.156;
  float b = -0.00137;
  //First loop over events for the normalization of the toppt reweighting///////

  if( id>300 && id<400 ){
    for( ULong64_t k = 0; k < nEventsTree; k++) {

      if( (k% 1000000) == 0)
	std::cout << "Finished event " << k << " / " << nEventsTree << std::endl;
      
      weight_toppt=1;

      int foundTop=0;
      int foundAntiTop=0;

      chain->GetEntry(k);
      if(nGenPart>0)
	for(int o=0; o<nGenPart; ++o){

	  if(GenPart_status[o] != 62) continue;
	  //Only apply weight to top or antitop
	  if( abs(GenPart_pdgId[o])==6 ){
	    weight_toppt *= sqrt(exp( a + b * GenPart_pt[o] ) );

	    if( GenPart_pdgId[o]== 6 )  foundTop++;
	    if( GenPart_pdgId[o]==-6 )  foundAntiTop++;
	  }
	  else continue;

	}//end loop over objects  

      if( foundTop==1 && foundAntiTop==1 )
	average += weight_toppt;  
      else 
	average += 1.0;
    }
    average /= (double) nEventsTree;
  }//end or first loop for toppt average

  */

  std::cout << "average = " << average << std::endl;

  std::cout << "sumGenWeightsHisto = " << sumGenWeightsHisto << std::endl;

  std::ofstream ofs (Form("%s", outputFile.c_str()), std::ofstream::out);
  ofs << sumGenWeightsHisto << "    ";
  ofs << average;
  ofs.close();



  delete chain; 
  //newH->Write();
  delete newH;
  // newSumW->Write();
  delete newSumW;
  return 0;
  
  }


