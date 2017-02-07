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

#include "btagSF.C"

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
  
  BTagSFHelper* bTagSFHelper =  new BTagSFHelper();

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



  float nTrueInt=0;
  //  int nTrueInt=0;
  TBranch* thisNTrueInt = (TBranch*) chain->GetListOfBranches()->FindObject("nTrueInt");
  if (thisNTrueInt) chain->SetBranchAddress("nTrueInt", &nTrueInt);
  
  if( nEventsTree > 0 )
    genWeight_ = fabs(genWeight);


  

  //  float scale1fb = xsec*kfactor*1000*filter/(Float_t)nEventsHisto;
  ULong64_t sumGenWeightsHisto; 
  ULong64_t nEffEventsHisto; 


  float weight_toppt;

  ////// b-tag SF
  Int_t njet;
  chain->SetBranchAddress("njet", &njet);
  Float_t jet_pt[100];
  chain->SetBranchAddress("jet_pt", jet_pt);
  Float_t jet_eta[100];
  chain->SetBranchAddress("jet_eta", jet_eta);
  Int_t jet_mcFlavour[100];
  chain->SetBranchAddress( "jet_hadronFlavour", jet_mcFlavour);
  Float_t jet_btagCSV[100];
  chain->SetBranchAddress("jet_btagCSV", jet_btagCSV);
  
  ////// isr re-weight
  Int_t nGenPart;
  Float_t GenPart_pt[100];
  Int_t GenPart_status[100];
  Int_t GenPart_pdgId[100];

  Int_t nisrMatch;

  Int_t GenSusyMNeutralino;
  Int_t GenSusyMGluino;

  bool isFastSim = 0;
  if( id>1000) isFastSim = 1;

  if( (id>1000 ) ){
    chain->SetBranchAddress("GenSusyMNeutralino", &GenSusyMNeutralino);

    if(  txtFileList.find("T2bb") != std::string::npos ){
      std::cout << "FOUND T2bb " << std::endl;  
      chain->SetBranchAddress("GenSusyMSbottom", &GenSusyMGluino);

    }else if(  txtFileList.find("T2qq") != std::string::npos ){
      std::cout << "FOUND T2qq " << std::endl;  
      chain->SetBranchAddress("GenSusyMSquark", &GenSusyMGluino);

    }else if(  txtFileList.find("T2tt") != std::string::npos ){
      std::cout << "FOUND T2tt " << std::endl;  
      chain->SetBranchAddress("GenSusyMStop", &GenSusyMGluino);

    }else if(  txtFileList.find("T1") != std::string::npos ){
      std::cout << "FOUND T1* " << std::endl;  
      chain->SetBranchAddress("GenSusyMGluino", &GenSusyMGluino);

    }	else std::cout << "Fuck" << std::endl;

  }   

  if (!isData) {
    std::cout << "Loading the generator parts " << std::endl;
    chain->SetBranchAddress("nisrMatch", &nisrMatch);
    // chain->SetBranchAddress("nGenPart", &nGenPart);
    // chain->SetBranchAddress("GenPart_pt", GenPart_pt);
    // chain->SetBranchAddress("GenPart_status", GenPart_status);
    // chain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
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
    

  //ISR /// bins of size 5GeV for T2cc
  Int_t nBins = 120*5;
  Float_t min = 12.5;
  Float_t max = 3012.5;

  //Need number of events per mass point for correct point
  TH2F* h_weight = new TH2F("h_weight", "", nBins, min, max, nBins, -min, max-min); h_weight->Sumw2();
  //Counter histogram for averaging
  TH2F* h_totalSumGenWeightsHisto = new TH2F("h_totalSumGenWeightsHisto", "", nBins, min, max, nBins, -min, max-min); h_totalSumGenWeightsHisto->Sumw2();

  //For SIGNAL (need average per mass point)
  TH2F* h_isr = new TH2F("h_isr", "", nBins, min, max, nBins, -min, max-min); h_isr->Sumw2();
  TH2F* h_isr_UP = new TH2F("h_isr_UP", "", nBins, min, max, nBins, -min, max-min); h_isr_UP->Sumw2();
  TH2F* h_isr_DN = new TH2F("h_isr_DN", "", nBins, min, max, nBins, -min, max-min); h_isr_DN->Sumw2();

  //For SIGNAL (need average per mass point)
  TH2F* h_btag = new TH2F("h_btag", "", nBins, min, max, nBins, -min, max-min); h_btag->Sumw2();
  TH2F* h_btag_heavy_UP = new TH2F("h_btag_heavy_UP", "", nBins, min, max, nBins, -min, max-min); h_btag_heavy_UP->Sumw2();
  TH2F* h_btag_heavy_DN = new TH2F("h_btag_heavy_DN", "", nBins, min, max, nBins, -min, max-min); h_btag_heavy_DN->Sumw2();
  TH2F* h_btag_light_UP = new TH2F("h_btag_light_UP", "", nBins, min, max, nBins, -min, max-min); h_btag_light_UP->Sumw2();
  TH2F* h_btag_light_DN = new TH2F("h_btag_light_DN", "", nBins, min, max, nBins, -min, max-min); h_btag_light_DN->Sumw2();

  //For MC 
  double isr_average = 0.;
  double weight_isr_av;

  double average = 0;
  average = 1.0; //dont' average for now
  
  //Btag weights, will get them damn averages
  Float_t weight_btagsf;
  Float_t weight_btagsf_heavy_UP;
  Float_t weight_btagsf_heavy_DN;
  Float_t weight_btagsf_light_UP;
  Float_t weight_btagsf_light_DN;

  double weight_btagsf_av=0.;
  double weight_btagsf_heavy_UP_av=0.;
  double weight_btagsf_heavy_DN_av=0.;
  double weight_btagsf_light_UP_av=0.;
  double weight_btagsf_light_DN_av=0.;

  //Top pt reweighting 
  //Values from Run1 still valid for now, see here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  float a = 0.156;
  float b = -0.00137;
  //First loop over events for the normalization of the toppt reweighting///////

  //if( id>300 && id<400 ){
  std::cout << "starting loop over events" << std::endl;

  for( ULong64_t k = 0; k <= nEventsTree; k++) {
    //   if( (k% 1000000) == 0)
    //    std::cout << "Finished event " << k << " / " << nEventsTree << std::endl;
     
    chain->GetEntry(k);

    weight_btagsf          = 1.;
    weight_btagsf_heavy_UP = 1.;
    weight_btagsf_heavy_DN = 1.;
    weight_btagsf_light_UP = 1.;
    weight_btagsf_light_DN = 1.;

 
    bTagSFHelper->get_weight_btag(njet, jet_pt, jet_eta, jet_mcFlavour, jet_btagCSV, weight_btagsf, weight_btagsf_heavy_UP, weight_btagsf_heavy_DN, weight_btagsf_light_UP, weight_btagsf_light_DN, isFastSim);

    weight_btagsf_av          += weight_btagsf/(double) nEventsTree ;
    weight_btagsf_heavy_UP_av += weight_btagsf_heavy_UP/(double) nEventsTree ;
    weight_btagsf_heavy_DN_av += weight_btagsf_heavy_DN/(double) nEventsTree ;
    weight_btagsf_light_UP_av += weight_btagsf_light_UP/(double) nEventsTree ;
    weight_btagsf_light_DN_av += weight_btagsf_light_DN/(double) nEventsTree ;

    weight_isr_av = 1.0;

    //NEW method with number of ISR jets
    float isr_err = 0.;
    if( nisrMatch == 1 ){
      weight_isr_av = 0.882;	    isr_err = 0.059;
    }else if( nisrMatch == 2 ){
      weight_isr_av = 0.792;	    isr_err = 0.104;
    }else if( nisrMatch == 3 ){ 
      weight_isr_av = 0.702;	    isr_err = 0.149;
    }else if( nisrMatch == 4 ){
      weight_isr_av = 0.648;	    isr_err = 0.176;
    }else if( nisrMatch == 5 ){
      weight_isr_av = 0.601;	    isr_err = 0.199;
    }else if( nisrMatch > 5 ){  
      weight_isr_av = 0.515;	    isr_err = 0.242;
    }

    if(id>10 && id<1000)
      isr_average += weight_isr_av /(double) nEventsTree ;
    else if( id>=1000){
      h_isr->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino,  weight_isr_av);

      h_isr_UP->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino,  weight_isr_av + isr_err);
      h_isr_DN->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino,  weight_isr_av - isr_err);

      h_totalSumGenWeightsHisto->Fill((float)GenSusyMGluino,(float)GenSusyMNeutralino,  1.0);

      h_btag         ->Fill((float)GenSusyMGluino, (float)GenSusyMNeutralino, weight_btagsf );
      h_btag_heavy_UP->Fill((float)GenSusyMGluino, (float)GenSusyMNeutralino, weight_btagsf_heavy_UP );
      h_btag_heavy_DN->Fill((float)GenSusyMGluino, (float)GenSusyMNeutralino, weight_btagsf_heavy_DN );
      h_btag_light_UP->Fill((float)GenSusyMGluino, (float)GenSusyMNeutralino, weight_btagsf_light_UP );
      h_btag_light_DN->Fill((float)GenSusyMGluino, (float)GenSusyMNeutralino, weight_btagsf_light_DN );

    }
      
    //   weight_toppt=1;
    //   int foundTop=0;   int foundAntiTop=0;
    //   chain->GetEntry(k);
    //   if(nGenPart>0)
    // 	for(int o=0; o<nGenPart; ++o){
    // 	  if(GenPart_status[o] != 62) continue;
    // 	  //Only apply weight to top or antitop
    // 	  if( abs(GenPart_pdgId[o])==6 ){
    // 	    weight_toppt *= sqrt(exp( a + b * GenPart_pt[o] ) );
    // 	    if( GenPart_pdgId[o]== 6 )  foundTop++;
    // 	    if( GenPart_pdgId[o]==-6 )  foundAntiTop++;
    // 	  }
    // 	  else continue;
    // 	}//end loop over objects  
    //   if( foundTop==1 && foundAntiTop==1 )
    // 	average += weight_toppt;  
    //   else 
    // 	average += 1.0;
    // }


  } // end of loop over events


  // average /= (double) nEventsTree;
  //  isr_average /= (double) nEventsTree;
  h_isr->Divide(h_totalSumGenWeightsHisto);
  h_isr_UP->Divide(h_totalSumGenWeightsHisto);
  h_isr_DN->Divide(h_totalSumGenWeightsHisto);
  
  h_btag->Divide(h_totalSumGenWeightsHisto);
  h_btag_heavy_UP->Divide(h_totalSumGenWeightsHisto);
  h_btag_heavy_DN->Divide(h_totalSumGenWeightsHisto);
  h_btag_light_UP->Divide(h_totalSumGenWeightsHisto);
  h_btag_light_DN->Divide(h_totalSumGenWeightsHisto);

  std::cout << "average = " << average << std::endl;

  std::cout << "sumGenWeightsHisto = " << sumGenWeightsHisto << std::endl;

  std::ofstream ofs (Form("%s.cfg", outputFile.c_str()), std::ofstream::out);
  ofs << sumGenWeightsHisto << " ";
  ofs << isr_average << " ";
  ofs << weight_btagsf_av << " ";
  ofs << weight_btagsf_heavy_UP_av << " ";
  ofs << weight_btagsf_heavy_DN_av << " ";
  ofs << weight_btagsf_light_UP_av << " ";
  ofs << weight_btagsf_light_DN_av << " ";

  // ofs << weight_btagsf_av /(double) nEventsTree << " ";
  // ofs << weight_btagsf_heavy_UP_av /(double) nEventsTree << " ";
  // ofs << weight_btagsf_heavy_DN_av /(double) nEventsTree << " ";
  // ofs << weight_btagsf_light_UP_av /(double) nEventsTree << " ";
  // ofs << weight_btagsf_light_DN_av /(double) nEventsTree << " ";
  // //  ofs << average;

  ofs.close();


  TFile* file_out = new TFile( Form("%s.root", outputFile.c_str()) ,"RECREATE");
  h_isr->Write();
  h_isr_UP->Write();
  h_isr_DN->Write();
  h_totalSumGenWeightsHisto->Write();

  h_btag->Write();
  h_btag_heavy_UP->Write();
  h_btag_heavy_DN->Write();
  h_btag_light_UP->Write();
  h_btag_light_DN->Write();

  file_out->Close();

  delete chain; 
  //newH->Write();
  delete newH;
  // newSumW->Write();
  delete newSumW;
  return 0;
  
  }


