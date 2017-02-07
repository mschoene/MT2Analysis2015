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
#include "TH3.h"
#include "TLorentzVector.h"

#include "goodrunClass.h"
#include "btagSF.C"

using namespace std;

int postProcessing(std::string inputString="input",
		   std::string inputFileList="", //		   std::string inputFolder="./",
		   std::string outputFile="output.root",
		   std::string treeName="tree",
		   float filter=1.0, float kfactor=1.0, float xsec=-1.0, int id=1,
		   std::string crabExt="",
		   std::string inputPU="",
		   std::string PUvar="nTrueInt",
		   bool applyJSON=true,
		   bool applySF=true,
		   bool doSilver=false,
		   std::string normFile="");


int run(std::string cfg="postProcessing.cfg",
	std::string treeName="tree", 
	std::string inputFolder = "pnfs/psi.ch/cms/trivcat/store/user/casal/babies/PHYS14_Production_QCDpt_noSietaieta/", 
	std::string outputFolder = "./test/",  
	std::string fileExtension = "_post.root",
        std::string crabExt = "",
	std::string inputPU = "",
	std::string PUvar = "nTrueInt",
	bool applyJSON=true){
  
  // for measuring timing
  time_t start = time(0);

  std::cout<<"Configuration file is: "<< cfg.c_str()<<std::endl;
  
  ifstream configuration(cfg.c_str());
  std::string line;
  
  while(std::getline(configuration,line)){
    istringstream ss(line);
    if((ss.str()[0])==(std::string)"#") continue;
    if(ss.str().empty()) continue;

    std::string idS,name,xsecS,filterS,kfactorS;
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
      std::cout << "id,name,x,f,k: " 
	   << id << " , " 
	   << name << " , "
	   << xsec << " , " 
	   << filter << " , " 
	   << kfactor << std::endl;
    }  

    std::string outputFile = outputFolder + "/" + name + fileExtension;
    postProcessing(name, inputFolder, outputFile, treeName, filter, kfactor, xsec, id, crabExt);
  }
  
  // for printing measured timing
  time_t stop = time(0);
  double time = difftime(stop, start);
  std::cout << "real time in seconds: " << time << std::endl;
  
  return 0;
}

int postProcessing(std::string inputString,
		   std::string inputFileList,
		   std::string outputFile,
		   std::string treeName,
		   float filter, float kfactor, float xsec, int id,
		   std::string crabExt,
		   std::string inputPU,
		   std::string PUvar,
		   bool applyJSON,
		   bool applySF,
		   bool doSilver,
		   std::string normFile){

  double totalSumGenWeightsHisto, topAverageWeight, isrAverageWeight, weight_btag_average, weight_btag_heavy_UP_average, weight_btag_heavy_DN_average, weight_btag_light_UP_average, weight_btag_light_DN_average;
 
  TH2F* h_isr;
  TH2F* h_isr_UP;
  TH2F* h_isr_DN;
  TH2F* h_totalSumGenWeightsHisto;

  TH2F* h_btag;
  TH2F* h_btag_heavy_UP;
  TH2F* h_btag_heavy_DN;
  TH2F* h_btag_light_UP;
  TH2F* h_btag_light_DN;


  if( normFile!="" ){
    ifstream configuration( Form("%s.cfg", normFile.c_str()) );
    std::string line;
    while(std::getline(configuration,line)){
      istringstream ss(line);
      if((ss.str()[0])==(std::string)"#") continue;
      if(ss.str().empty()) continue;

      //      std::string normName;
      //      ss >> normName;
      ss >> totalSumGenWeightsHisto;
      ss >> isrAverageWeight;
     
      ss >> weight_btag_average; 
      ss >> weight_btag_heavy_UP_average; 
      ss >> weight_btag_heavy_DN_average;
      ss >> weight_btag_light_UP_average;
      ss >> weight_btag_light_DN_average;

      //ss >> topAverageWeight;

      std::cout<< "Got a line with " << std::endl;
      std::cout<< "sumGenWeightsHisto = " << totalSumGenWeightsHisto << std::endl;
      std::cout<< "weight_isr_av    = " << isrAverageWeight << std::endl;
      std::cout<< "weight_btag_av    = " << weight_btag_average << std::endl;
      //      std::cout<< "weight_topPt_av    = " << topAverageWeight << std::endl;

      TFile* f_temp = new TFile(Form("%s.root", normFile.c_str() ));
      h_isr = (TH2F*) f_temp->Get("h_isr");
      h_isr_UP = (TH2F*) f_temp->Get("h_isr_UP");
      h_isr_DN = (TH2F*) f_temp->Get("h_isr_DN");

      h_btag = (TH2F*) f_temp->Get("h_btag");
      h_btag_heavy_UP = (TH2F*) f_temp->Get("h_btag_heavy_UP");
      h_btag_heavy_DN = (TH2F*) f_temp->Get("h_btag_heavy_DN");
      h_btag_light_UP = (TH2F*) f_temp->Get("h_btag_light_UP");
      h_btag_light_DN = (TH2F*) f_temp->Get("h_btag_light_DN");

      h_totalSumGenWeightsHisto = (TH2F*) f_temp->Get("h_totalSumGenWeightsHisto");

      if (!h_totalSumGenWeightsHisto || !h_isr || !h_btag ) std::cout << "ERROR: Could not find averaging histogram"<< std::endl;
    }
  }


  //bool applyJSON=true;
  const char* goldenjson_file = "goodruns_golden.txt";
  const char* silverjson_file = "goodruns_silver.txt";

  GoodRun golden;
  GoodRun silver;

  if (applyJSON) {
    //std::cout << gSystem->pwd() << std::endl;
    //gSystem->Load("goodrun_cc");
    std::cout << "Loading golden json file: " << goldenjson_file << std::endl;
    golden.set_goodrun_file(goldenjson_file);
    if (doSilver) {
      std::cout << "Loading silver json file: " << silverjson_file << std::endl;
      silver.set_goodrun_file(silverjson_file);
    }
  }

  TH2D* h_id = 0;
  TH2D* h_iso = 0;
  TH2D* h_elSF = 0;
  TH2D* h_id_mu = 0;
  TH2D* h_iso_mu = 0;
  TH2D* h_dxyz_mu = 0;
  TH2D* h_muSF = 0;

  TH2D* h_elTrk = 0;
  TH1D* h_muTrk_hi = 0;
  TH1D* h_muTrk_lo = 0;
 
  TH2D* h_fast_muSF = 0;
  TH2D* h_fast_elSF = 0;
  TH2D* h_eff_full_mu = 0;
  TH2D* h_eff_full_el = 0;
  TH2D* h_eff_fast_mu =  0;
  TH2D* h_eff_fast_el =  0;

  BTagSFHelper* bTagSFHelper = 0;
  if( applySF ){
    bTagSFHelper =  new BTagSFHelper();

    //Getting the lepton scale factor histograms/////////////////
    //Electrons//
    std::string filename = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/scaleFactors.root";
    TFile * f_ele = new TFile(filename.c_str() );
    if (!f_ele->IsOpen()) std::cout << " ERROR: Could not find scale factor file " << filename << std::endl; 
    //Uncomment for loose Id
    //TH2D* h_id = (TH2D*) f_ele->Get("CutBasedLoose");
    //(TH2D*) f_ele->Get("CutBasedVeto");
    h_id = (TH2D*) f_ele->Get("GsfElectronToVeto");
    h_iso = (TH2D*) f_ele->Get("MVAVLooseElectronToMini");
    if (!h_id || !h_iso) std::cout << "ERROR: Could not find scale factor histogram"<< std::endl;
    h_elSF = (TH2D*) h_id->Clone("h_elSF");
    h_elSF->SetDirectory(0);
    h_elSF->Multiply(h_iso);

    std::string filenameElTrk = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/egammaEffi_SF2D.root";
    TFile * f_eleTrk = new TFile(filenameElTrk.c_str() );
    if (!f_eleTrk->IsOpen()) std::cout << " ERROR: Could not find scale factor file " << filenameElTrk << std::endl; 
    h_elTrk = (TH2D*) f_eleTrk->Get("EGamma_SF2D");
    h_elTrk->SetDirectory(0);
    f_eleTrk->Close(); delete f_eleTrk; 


    //Muons//
    std::string filenameID = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/TnP_MuonID_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root";
    std::string filenameISO = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/TnP_MuonID_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root";
    std::string filenamedxyz = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/TnP_MuonID_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root";
    TFile * f1 = new TFile(filenameID.c_str() );
    TFile * f2 = new TFile(filenameISO.c_str() );
    TFile * f3 = new TFile(filenamedxyz.c_str() );
    if (!f1->IsOpen()) { std::cout<<" ERROR: Could not find ID scale factor file "<<filenameID<<std::endl; return 0;}
    if (!f2->IsOpen()) { std::cout<<"ERROR: Could not find ISO scale factor file "<<filenameISO<<std::endl; return 0;}
    if (!f3->IsOpen()) { std::cout<<"ERROR: Could not find dxy dz scale factor file "<<filenamedxyz<<std::endl; return 0;}
    h_id_mu = (TH2D*) f1->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0");
    h_iso_mu = (TH2D*) f2->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_PF_pass");
    h_dxyz_mu = (TH2D*) f3->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_PF_pass");
    if (!h_id_mu || !h_iso_mu  || !h_dxyz_mu) { std::cout<<"ERROR: Could not find scale factor histogram"<<std::endl; return 0;}
    h_muSF = (TH2D*) h_id_mu->Clone("h_muSF");
    h_muSF->SetDirectory(0);
    h_muSF->Multiply(h_iso_mu);
    h_muSF->Multiply(h_dxyz_mu);

    f_ele->Close();  f1->Close(); f2->Close();
    delete f_ele; delete f1; delete f2;

    TH1D* h_trk_mu_hi = 0;
    TH1D* h_trk_mu_lo = 0;

    std::string filenameTrk = "/mnt/t3nfs01/data01/shome/mschoene/lepSF/general_tracks_and_early_general_tracks_corr_ratio.root";
    TFile * fTrk = new TFile(filenameTrk.c_str() );
    if (!fTrk->IsOpen()) { std::cout<<" ERROR: Could not find track ineff scale factor file "<<filenameTrk<<std::endl; return 0;}
    h_trk_mu_hi = (TH1D*) fTrk->Get("mutrksfptg10");
    if (!h_trk_mu_hi) { std::cout<<"ERROR: Could not find trk sf histogram"<<std::endl; return 0;}
    h_muTrk_hi = (TH1D*) h_trk_mu_hi->Clone("h_muTrk_hi");
    h_muTrk_hi->SetDirectory(0);
    h_trk_mu_lo = (TH1D*) fTrk->Get("mutrksfptl10");
    if (!h_trk_mu_lo) { std::cout<<"ERROR: Could not find trk sf histogram"<<std::endl; return 0;}
    h_muTrk_lo = (TH1D*) h_trk_mu_lo->Clone("h_muTrk_lo");
    h_muTrk_lo->SetDirectory(0);
    fTrk->Close(); delete fTrk;


    std::cout << std::endl;
    std::cout << "Using Loose Muon ID, MiniIso 0.2 lepton scale factors" << std::endl;
    std::cout << "Using Veto Electrons ID, MiniIso 0.1 lepton scale factors" << std::endl;
    std::cout << "Be aware that Veto Electrons are not suited for selecting Electrons." << std::endl;
    std::cout << std::endl;

    std::cout << "Also loading the FullSim efficiency map" << std::endl;
    TFile * f_eff_full = new TFile("/mnt/t3nfs01/data01/shome/mschoene/CMSSW_7_4_12_pP/src/analysisCode/Utils/vetoeff_emu_etapt_lostlep.root" );
    if(!f_eff_full->IsOpen()) {std::cout<<" ERROR: Could not find muon Fullsim scale factor file" <<std::endl; return 0;}
    h_eff_full_mu = (TH2D*) f_eff_full->Get("h_mu_comb_eff");
    h_eff_full_el = (TH2D*) f_eff_full->Get("h_ele_comb_eff");
    if(!h_eff_full_mu || !h_eff_full_el ) {std::cout << " ERROR: Could not find the 2D histogram in your files " << std::endl; return 0;}
    h_eff_full_mu->SetDirectory(0);
    h_eff_full_el->SetDirectory(0);
    f_eff_full->Close();
    delete f_eff_full;


    if( id >= 1000 ){
      std::cout << "Also loading the FastSim/FullSim Lepton scale factors" << std::endl;
      TFile * f_mu = new TFile("/mnt/t3nfs01/data01/shome/mschoene/CMSSW_7_4_12_pP/src/analysisCode/Utils/sf_mu_looseID_mini02.root" );
      TFile * f_el = new TFile("/mnt/t3nfs01/data01/shome/mschoene/CMSSW_7_4_12_pP/src/analysisCode/Utils/sf_el_vetoCB_mini01.root" );
      if(!f_mu->IsOpen()) { std::cout<<" ERROR: Could not find muon Fastsim scale factor file " <<std::endl; return 0;}
      if(!f_el->IsOpen()) { std::cout<<" ERROR: Could not find electron Fastsim scale factor file " <<std::endl; return 0;}

      h_fast_muSF = (TH2D*) f_mu->Get("histo2D");
      h_fast_elSF = (TH2D*) f_el->Get("histo2D");
      if(!h_fast_muSF || !h_fast_elSF ) {std::cout << " ERROR: Could not find the 3D histogram in your files " << std::endl; return 0;}

      h_fast_muSF->SetDirectory(0);
      h_fast_elSF->SetDirectory(0);

      f_mu->Close(); f_el->Close();
      delete f_mu; delete f_el;


      std::cout << "Also loading the FastSim efficiency map" << std::endl; 
      TFile * f_eff_fast = new TFile("/mnt/t3nfs01/data01/shome/mschoene/CMSSW_7_4_12_pP/src/analysisCode/Utils/vetoeff_emu_etapt_T1tttt_mGluino-1500to1525.root" );
      if(!f_eff_fast->IsOpen()) {std::cout<<" ERROR: Could not find muon Fastsim scale factor file" <<std::endl; return 0;}
      h_eff_fast_mu = (TH2D*) f_eff_fast->Get("h_ele_comb_eff");
      h_eff_fast_el = (TH2D*) f_eff_fast->Get("h_mu_comb_eff");
      if(!h_eff_fast_mu || !h_eff_fast_el ) {std::cout << " ERROR: Could not find the 2D histogram in your files " << std::endl; return 0;}
      h_eff_fast_mu->SetDirectory(0);
      h_eff_fast_el->SetDirectory(0);
      f_eff_fast->Close();
      delete f_eff_fast;
    }
  }


  TChain* chain = new TChain(treeName.c_str());

  std::cout << inputFileList << std::endl;
  TFileCollection *filelist = new TFileCollection("listOfFiles");
  
  std::cout << "Adding file list... " << std::endl;
  filelist->AddFromFile( inputFileList.c_str() );
  
  std::cout << "Adding files to chain... " << std::endl;  
  int chainReturn = chain->AddFileInfoList((TCollection*) filelist->GetList());
 
  if (chainReturn < 1) {
    std::cout << "ERROR: input folder/fileName is not well defined. Exit!" << std::endl;
    std::cout << "InputFileList: " << inputFileList << std::endl;
    return 1;
  }
  

  TFile * f_puData = new TFile( inputPU.c_str() );
  if (!f_puData->IsOpen()) { std::cout<<"ERROR: Could not find data pile up histogram"<<std::endl; return 0;}
  TH2D* hPU_data = (TH2D*) f_puData->Get("pileup");

  ULong64_t nData = hPU_data->Integral();
  hPU_data->Scale(1.0/nData);

  std::cout << "Initialized PU file." << std::endl;
  
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

  // This line should be uncommented for all the branches that we want to overwrite.
  // If the branch is not in the input tree, we don't need this.
  //
  // t->SetBranchStatus("scale1fb", 0);


  TFile *out = TFile::Open(outputFile.c_str(), "RECREATE");
  TTree *clone = new TTree("mt2", "post processed baby tree for mt2 analysis");

  // clone->AutoSave();
  //  clone->SetAutoSave( - 10000000 );

  TBranch* thisPUWeight = (TBranch*) chain->GetListOfBranches()->FindObject("puWeight");
  // if (thisPUWeight) chain->SetBranchStatus("puWeight", 0);

  std::cout << "Cloning tree..." << std::endl;

  clone = chain->CloneTree(-1, "fast"); 
  clone->SetName("mt2");
  
  std::cout << "Cloned tree." << std::endl;
 
  //if(SortBasketsByEntry)
  // clone = t->CloneTree(-1, "fastSortBasketsByEntry");
  //else 
  // clone = t->CloneTree(-1, "fast");
  
  //-------------------------------------------------------------

  //Calculate scaling factor and put variables into tree 
  ULong64_t nEventsTree = clone->GetEntries();
  ULong64_t nEventsHisto = (ULong64_t)newH->GetBinContent(1);		 
  //  Int_t nEventsTree = clone->GetEntries();
  //  Int_t nEventsHisto = (Int_t) newH->GetBinContent(1); 

  float genWeight_=1.0;
  float genWeight=0.;
  TBranch* thisGenWeight = (TBranch*) chain->GetListOfBranches()->FindObject("genWeight");
  if (thisGenWeight) chain->SetBranchAddress("genWeight", &genWeight);


  int isData=0; 
  chain->SetBranchAddress("isData",&isData);
  chain->GetEntry(0);

  UInt_t run=0;
  chain->SetBranchAddress("run", &run);
  UInt_t lumi=0;
  chain->SetBranchAddress("lumi", &lumi);

  int nVert=0;
  chain->SetBranchAddress("nVert", &nVert);



  float nTrueInt=0;
  //int nTrueInt=0;
  TBranch* thisNTrueInt = (TBranch*) chain->GetListOfBranches()->FindObject("nTrueInt");
  if (thisNTrueInt) chain->SetBranchAddress("nTrueInt", &nTrueInt);
  
  if( nEventsTree > 0 ){
    chain->GetEntry(0);
    if(isData){
      genWeight_ = 1.;
      applySF = false; //to make absolutely sure no plonker tries to fill them for data
    }
    else
      genWeight_ = fabs(genWeight);
  }



  TH1D* hPU = (TH1D*) hPU_data->Clone("hPU");
 
  
  // if(PUvar == "nVert")
  //   chain->Project("hPU", "nVert");
  // else 
  if(!isData){
    hPU->Reset();
    chain->Project("hPU", "nTrueInt");
  }

  if( hPU->Integral() >0.0 )
    hPU->Scale(1.0/hPU->Integral() );
  
  TH1D* hPU_r = (TH1D*) hPU_data->Clone("hPU_r");
  hPU_r->Divide(hPU);
  //  std::cout << hPU_r->GetBinContent(1)<< std::endl;

  //  float scale1fb = xsec*kfactor*1000*filter/(Float_t)nEventsHisto;
  ULong64_t sumGenWeightsHisto; 
  ULong64_t nEffEventsHisto; 
  float scale1fb_noGenWeight;
  float scale1fb_sumGenWeights;
  float scale1fb;
  float puWeight;
  int isGolden;
  int isSilver;

  Float_t weight_btagsf;
  Float_t weight_btagsf_heavy_UP;
  Float_t weight_btagsf_heavy_DN;
  Float_t weight_btagsf_light_UP;
  Float_t weight_btagsf_light_DN;
  Float_t weight_lepsf;
  Float_t weight_lepsf_UP;
  Float_t weight_lepsf_DN;
  Float_t weight_toppt;
  Float_t weight_isr;
  Float_t weight_isr_UP;
  Float_t weight_isr_DN;
  Float_t weight_isr_av;
  Float_t weight_isr_UP_av;
  Float_t weight_isr_DN_av;
  Float_t weight_scales[110];
  Float_t weight_scales_av[110];

  Float_t weight_lepsf_0l;
  Float_t weight_lepsf_0l_UP;
  Float_t weight_lepsf_0l_DN;


  ////// Lepton Efficiency SF
  Int_t nlep;
  chain->SetBranchAddress("nlep", &nlep);
  Float_t lep_pt[100];
  chain->SetBranchAddress("lep_pt", lep_pt);
  Float_t lep_eta[100];
  chain->SetBranchAddress("lep_eta", lep_eta);
  Int_t lep_pdgId[100];
  chain->SetBranchAddress("lep_pdgId", lep_pdgId);
  

  Int_t ngenLep;
  Float_t genLep_pt[100];
  Float_t genLep_eta[100];  
  Int_t genLep_pdgId[100];
  Int_t ngenLepFromTau;
  Float_t genLepFromTau_pt[100];
  Float_t genLepFromTau_eta[100];  
  Int_t genLepFromTau_pdgId[100];
  Int_t nMuons10;
  Int_t nElectrons10;
  Int_t nPFLep5LowMT;
  Int_t nPFHad10LowMT;

  Int_t nisrMatch;

  if(!isData){
    chain->SetBranchAddress("nisrMatch", &nisrMatch);
    
    chain->SetBranchAddress("ngenLep", &ngenLep);
    chain->SetBranchAddress("genLep_pt", genLep_pt);
    chain->SetBranchAddress("genLep_eta", genLep_eta);
    chain->SetBranchAddress("genLep_pdgId", genLep_pdgId);

    chain->SetBranchAddress("ngenLepFromTau", &ngenLepFromTau);
    chain->SetBranchAddress("genLepFromTau_pt", genLepFromTau_pt);
    chain->SetBranchAddress("genLepFromTau_eta", genLepFromTau_eta);
    chain->SetBranchAddress("genLepFromTau_pdgId", genLepFromTau_pdgId);

    chain->SetBranchAddress("nMuons10", &nMuons10);
    chain->SetBranchAddress("nElectrons10", &nElectrons10);
    chain->SetBranchAddress("nPFLep5LowMT", &nPFLep5LowMT);
    chain->SetBranchAddress("nPFHad10LowMT", &nPFHad10LowMT);
  }


  ////// b-tag SF
  Int_t njet;
  chain->SetBranchAddress("njet", &njet);
  Float_t jet_pt[100];
  chain->SetBranchAddress("jet_pt", jet_pt);
  Float_t jet_eta[100];
  chain->SetBranchAddress("jet_eta", jet_eta);
  Int_t jet_mcFlavour[100];
  if(!isData)
    chain->SetBranchAddress( "jet_hadronFlavour", jet_mcFlavour);
  //chain->SetBranchAddress("jet_mcFlavour", jet_mcFlavour);
  Float_t jet_btagCSV[100];
  chain->SetBranchAddress("jet_btagCSV", jet_btagCSV);
  
  ////// isr re-weight
  Int_t nGenPart;
  Float_t GenPart_pt[100];
  Float_t GenPart_eta[100];
  Float_t GenPart_phi[100];
  Float_t GenPart_mass[100];
  Int_t GenPart_status[100];
  Int_t GenPart_pdgId[100];

  if (!isData) {
    std::cout << "Loading the generator parts " << std::endl;
    chain->SetBranchAddress("nGenPart", &nGenPart);
    chain->SetBranchAddress("GenPart_pt", GenPart_pt);
    chain->SetBranchAddress("GenPart_eta", GenPart_eta);
    chain->SetBranchAddress("GenPart_phi", GenPart_phi);
    chain->SetBranchAddress("GenPart_mass", GenPart_mass);
    chain->SetBranchAddress("GenPart_status", GenPart_status);
    chain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
  }

 
  Float_t LHEweight_original;
  Float_t LHEweight_wgt[510];
  if (id>1000 && id <2000) {
    chain->SetBranchAddress("LHEweight_original", &LHEweight_original);
    chain->SetBranchAddress("LHEweight_wgt", LHEweight_wgt);
  }


  if( isData ){
    sumGenWeightsHisto = (ULong64_t) 0.0;
    nEffEventsHisto    = (ULong64_t) nEventsHisto;
    
    scale1fb_noGenWeight   = (Float_t) 0.0;
    scale1fb_sumGenWeights = (Float_t) 0.0;
  }else{
 
    sumGenWeightsHisto = (ULong64_t) newSumW->GetBinContent(1);
    if( normFile!="" )
      sumGenWeightsHisto =  totalSumGenWeightsHisto; 

    nEffEventsHisto = ( (double) 1.0*sumGenWeightsHisto/genWeight_  > (ULong64_t) (1.0*sumGenWeightsHisto/genWeight_ + 0.5) ) ? (ULong64_t) (1.0*sumGenWeightsHisto/genWeight_ + 1.0) : (ULong64_t) (1.0*sumGenWeightsHisto/genWeight_);
  
    scale1fb_noGenWeight = xsec*kfactor*1000*filter/(Float_t)nEffEventsHisto;
    scale1fb_sumGenWeights = xsec*kfactor*1000*filter/(Float_t)sumGenWeightsHisto;
  }
  
  if (nEventsHisto < nEventsTree) // this should not happen
    std::cout << "ERROR: histogram count has less events than tree. This indicates something went wrong" << std::endl
	      << "#events histo: "  << nEventsHisto << std::endl
	      << "#events tree: "  << nEventsTree << std::endl;
  else if (nEventsHisto > nEventsTree) // this can happen
    std::cout << "WARNING: histogram count has more events than tree. This should only happen if tree was skimmed" << std::endl
	      << "#events histo: "  << nEventsHisto << std::endl
	      << "#events tree: "  << nEventsTree << std::endl;
    
  bool isFastSim=0;
  if( id>=1000 && id<=2000)
    isFastSim=1;


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
  TBranch* b11 = clone->Branch("puWeight", &puWeight, "puWeight/F");
  TBranch* b12 = clone->Branch("isGolden", &isGolden, "isGolden/I");
  TBranch* b13 = clone->Branch("isSilver", &isSilver, "isSilver/I");
 
  TBranch* b14 = 0;
  TBranch* b15 = 0; 
  TBranch* b16 = 0;
  TBranch* b17 = 0;
  TBranch* b18 = 0;
  TBranch* b19 = 0; 
  TBranch* b20 = 0; 
  TBranch* b21 = 0; 
  TBranch* b22 = 0;
  TBranch* b23 = 0; 
  TBranch* b24 = 0; 
  TBranch* b25 = 0; 
  TBranch* b26 = 0;
  TBranch* b27 = 0;
  TBranch* b28 = 0; 
  TBranch* b29 = 0; 
  TBranch* b30 = 0; 
  TBranch* b31 = 0; 

  TBranch* b32 = 0; 
  TBranch* b33 = 0; 
  TBranch* b34 = 0; 
  TBranch* b35 = 0; 
  TBranch* b36 = 0; 

  TBranch* b37 = 0; 
  TBranch* b38 = 0; 

  bool doIsrAv = 1;
  
  if( applySF ){
    b14 = clone->Branch("weight_btagsf"         , &weight_btagsf         , "weight_btagsf/F"         );
    b15 = clone->Branch("weight_btagsf_heavy_UP", &weight_btagsf_heavy_UP, "weight_btagsf_heavy_UP/F");
    b16 = clone->Branch("weight_btagsf_heavy_DN", &weight_btagsf_heavy_DN, "weight_btagsf_heavy_DN/F");
    b17 = clone->Branch("weight_btagsf_light_UP", &weight_btagsf_light_UP, "weight_btagsf_light_UP/F");
    b18 = clone->Branch("weight_btagsf_light_DN", &weight_btagsf_light_DN, "weight_btagsf_light_DN/F");
    b19 = clone->Branch("weight_lepsf",    &weight_lepsf,    "weight_lepsf/F");
    b20 = clone->Branch("weight_lepsf_UP", &weight_lepsf_UP, "weight_lepsf_UP/F");
    b21 = clone->Branch("weight_lepsf_DN", &weight_lepsf_DN, "weight_lepsf_DN/F");

    b22 = clone->Branch("weight_toppt", &weight_toppt, "weight_toppt/F");

    b23 = clone->Branch("weight_isr",      &weight_isr,      "weight_isr/F");
    b24 = clone->Branch("weight_isr_av", &weight_isr_av, "weight_isr_av/F");
    b25 = clone->Branch("weight_isr_UP",   &weight_isr_UP,   "weight_isr_UP/F");
    b26 = clone->Branch("weight_isr_DN",   &weight_isr_DN,   "weight_isr_DN/F");

    b27 = clone->Branch("weight_isr_UP_av",   &weight_isr_UP_av,   "weight_isr_UP_av/F");
    b28 = clone->Branch("weight_isr_DN_av",   &weight_isr_DN_av,   "weight_isr_DN_av/F");

    b29 = clone->Branch("weight_scales", &weight_scales, "weight_scales[110]/F");
    b30 = clone->Branch("weight_scales_av", &weight_scales_av, "weight_scales_av[110]/F");

    b31 = clone->Branch("weight_lepsf_0l", &weight_lepsf_0l, "weight_lepsf_0l/F");
    b32 = clone->Branch("weight_lepsf_0l_UP", &weight_lepsf_0l_UP, "weight_lepsf_0l_UP/F");
    b33 = clone->Branch("weight_lepsf_0l_DN", &weight_lepsf_0l_DN, "weight_lepsf_0l_DN/F");

    b34 = clone->Branch("weight_btagsf_av"         , &weight_btag_average         , "weight_btagsf_av/F"         );
    b35 = clone->Branch("weight_btagsf_heavy_UP_av", &weight_btag_heavy_UP_average, "weight_btagsf_heavy_UP_av/F");
    b36 = clone->Branch("weight_btagsf_heavy_DN_av", &weight_btag_heavy_DN_average, "weight_btagsf_heavy_DN_av/F");
    b37 = clone->Branch("weight_btagsf_light_UP_av", &weight_btag_light_UP_average, "weight_btagsf_light_UP_av/F");
    b38 = clone->Branch("weight_btagsf_light_DN_av", &weight_btag_light_DN_average, "weight_btagsf_light_DN_av/F");
  }


  vector<TH2F*> h_scales;
  vector<TH2F*> h_evt;

  //ISR /// bins of size 5GeV for T2cc
  Int_t nBins = 120*5;
  Float_t min = 12.5;
  Float_t max = 3012.5;

  Float_t weight_isr_av_firstLoop;
  if(normFile=="" && id>1000){
    h_isr = new TH2F("h_isr", "", nBins, min, max, nBins, -min, max-min); h_isr->Sumw2();
  
    h_btag = new TH2F("h_btag", "", nBins, min, max, nBins, -min, max-min); h_btag->Sumw2();
    h_btag_heavy_UP = new TH2F("h_btag_heavy_UP", "", nBins, min, max, nBins, -min, max-min); h_btag_heavy_UP->Sumw2();
    h_btag_heavy_DN = new TH2F("h_btag_heavy_DN", "", nBins, min, max, nBins, -min, max-min); h_btag_heavy_DN->Sumw2();
    h_btag_light_UP = new TH2F("h_btag_light_UP", "", nBins, min, max, nBins, -min, max-min); h_btag_light_UP->Sumw2();
    h_btag_light_DN = new TH2F("h_btag_light_DN", "", nBins, min, max, nBins, -min, max-min); h_btag_light_DN->Sumw2();
  }

  TH2F* h_counter = new TH2F("h_counter", "", nBins, min, max, nBins, -min, max-min); h_counter->Sumw2();
  Int_t GenSusyMNeutralino;
  Int_t GenSusyMGluino;

  if( (id>1000) ){
    chain->SetBranchAddress("GenSusyMNeutralino", &GenSusyMNeutralino);

    if(  normFile.find("T2bb") != std::string::npos ){
      std::cout << "FOUND T2bb " << std::endl;  
      chain->SetBranchAddress("GenSusyMSbottom", &GenSusyMGluino);

    }else if(  normFile.find("T2qq") != std::string::npos ){
      std::cout << "FOUND T2qq " << std::endl;  
      chain->SetBranchAddress("GenSusyMSquark", &GenSusyMGluino);

    }else if(  normFile.find("T2tt") != std::string::npos ){
      std::cout << "FOUND T2tt " << std::endl;  
      chain->SetBranchAddress("GenSusyMStop", &GenSusyMGluino);

    }else if(  normFile.find("T1") != std::string::npos ){
      std::cout << "FOUND T1* " << std::endl;  
      chain->SetBranchAddress("GenSusyMGluino", &GenSusyMGluino);

    }	else std::cout << "Fuck" << std::endl;
  }

  double isr_average = 1.;
   
  if( applySF  ){

    if( id>1000 ) {
      for(int k=0;k<110;k++){
	std::string name = "var"+  std::to_string(k);
	TH2F* h_scales_temp = new TH2F(name.c_str(), "", nBins, min, max, nBins, -min, max-min); h_scales_temp->Sumw2();
	h_scales.push_back(h_scales_temp);
 
	std::string name_evt = "evt"+  std::to_string(k); 
	TH2F* h_evt_temp = new TH2F(name_evt.c_str(), "", nBins, min, max, nBins, -min, max-min); h_evt_temp->Sumw2();
	h_evt.push_back(h_evt_temp);
      }
    }


    if( normFile=="" ){

      std::cout << "Entering first loop to determine the average ISR weights" << std::endl;
      for(  ULong64_t j = 0; j < nEventsTree; j++) {
    
	chain->GetEntry(j);
 
	if(normFile==""){
	  weight_isr_av_firstLoop = 1.0;
     
	  //NEW method with number of ISR jets

	  if( nisrMatch == 1 )      weight_isr_av_firstLoop = 0.882;
	  else if( nisrMatch == 2 ) weight_isr_av_firstLoop = 0.792;
	  else if( nisrMatch == 3 ) weight_isr_av_firstLoop = 0.702;
	  else if( nisrMatch == 4 ) weight_isr_av_firstLoop = 0.648;
	  else if( nisrMatch == 5 ) weight_isr_av_firstLoop = 0.601;
	  else if( nisrMatch > 5 )  weight_isr_av_firstLoop = 0.515;
 
	  if(id>10 && id<1000)
	    isr_average += weight_isr_av_firstLoop;
	  else{
	    h_isr->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino,  weight_isr_av_firstLoop);
	    h_counter->Fill((float)GenSusyMGluino,(float)GenSusyMNeutralino,  1.0);
	  }

	  if( id>1000){
	    Float_t weight_btagsf = 1.0;
	    Float_t weight_btagsf_heavy_UP = 1.0;
	    Float_t weight_btagsf_heavy_DN = 1.0;
	    Float_t weight_btagsf_light_UP = 1.0;
	    Float_t weight_btagsf_light_DN = 1.0;
	  
	    bTagSFHelper->get_weight_btag(njet, jet_pt, jet_eta, jet_mcFlavour, jet_btagCSV, weight_btagsf, weight_btagsf_heavy_UP, weight_btagsf_heavy_DN, weight_btagsf_light_UP, weight_btagsf_light_DN, isFastSim);

	    h_btag->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino, weight_btagsf );
	    h_btag_heavy_UP->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino, weight_btagsf_heavy_UP );
	    h_btag_heavy_DN->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino, weight_btagsf_heavy_DN );
	    h_btag_light_UP->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino, weight_btagsf_light_UP );
	    h_btag_light_DN->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino, weight_btagsf_light_DN );
	  }

	}
      
	if( id < 2000. && id > 1000.){
	  for(int v=0; v<110; ++v){
	    weight_scales[v] = 1.;
	    weight_scales[v] = LHEweight_wgt[v]/LHEweight_original;

	    h_scales[v]->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino,  weight_scales[v]);
	    h_evt[v]->Fill((float)GenSusyMGluino,(float)GenSusyMNeutralino,  1.0);
	  }
	}
    
      }//end or first loop over events for calculating the average ISR


      if(normFile==""){
	h_isr->Divide(h_counter);

  	h_btag->Divide(h_counter);
	h_btag_heavy_UP->Divide(h_counter);
	h_btag_heavy_DN->Divide(h_counter);
	h_btag_light_UP->Divide(h_counter);
	h_btag_light_DN->Divide(h_counter);
      }
    
      if( id<2000 && id>1000. ){
	for(int v=0; v<110; ++v){
	  h_scales[v]->Divide(h_evt[v]);
	}
      }
    } //end of if normfile, ie this is only done if not split
    
    /* 
       for(int l = 0; l<120*5; l++){
       for(int m = 0; m<120*5; m++){
	if( h_scales[5]->GetBinContent(l,m) > 0.5)
	  std::cout << l<< ", " << m << " has content " <<h_scales[5]->GetBinContent(l,m)  << std::endl;
      }
    }    
    */
    std::cout << "Finished first loop over the events to get the average weight of ISR" << std::endl;
  }//end of ISR norm histo calculation

  if( normFile=="" && id>10 )
    isr_average /= (double) nEventsTree;
  
  //Top pt reweighting 
  //Values from Run1 still valid for now, see here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  float a = 0.156;
  float b = -0.00137;

  //NOT averaging the top pt weight 
  /*
  //First loop over events for the normalization of the toppt reweighting///////
  double average = 0;
  if( applySF && normFile==""&& id>300 && id<400 ){
    for( ULong64_t k = 0; k < nEventsTree; k++) {
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

  if(normFile!="")
    average = topAverageWeight;

  */



  std::cout << "Entering the final loop over the events" << std::endl;

  for( Long64_t i = 0; i < (Long64_t)nEventsTree; i++) {
  
    if( i==0)
      std::cout << "Start of loop over tree entries" << std::endl;

    
    weight_btag_average = 1.;
    weight_btag_heavy_UP_average=  1.;
    weight_btag_heavy_DN_average=  1.;
    weight_btag_light_UP_average=  1.;
    weight_btag_light_DN_average=  1.;
  
    weight_btagsf          = 1.;
    weight_btagsf_heavy_UP = 1.;
    weight_btagsf_heavy_DN = 1.;
    weight_btagsf_light_UP = 1.;
    weight_btagsf_light_DN = 1.;
    weight_lepsf           = 1.;
    weight_lepsf_UP        = 1.;
    weight_lepsf_DN        = 1.;
    weight_toppt           = 1.;
    weight_isr             = 1.;
    weight_isr_UP          = 1.;
    weight_isr_DN          = 1.;
    weight_lepsf_0l        = 1.;
    weight_lepsf_0l_UP     = 1.;
    weight_lepsf_0l_DN     = 1.;

    for(int v=0; v<110; ++v){
      weight_scales[v]=1.;
      weight_scales_av[v]=1.;
    }
    

    chain->GetEntry(i);

    if(isData){ 
      //data should never be rescaled with scale1fb when making plots
      //set scaler to zero so that user cannot miss the mistake
      scale1fb = 0.0; 
      puWeight = 1.0;
    }else{
      if(id >= 1000 && id<2000. && normFile!=""){
	Int_t binx = h_totalSumGenWeightsHisto->GetXaxis()->FindBin( (float)GenSusyMGluino );
	Int_t biny = h_totalSumGenWeightsHisto->GetYaxis()->FindBin( (float)GenSusyMNeutralino );
	sumGenWeightsHisto = h_totalSumGenWeightsHisto->GetBinContent(binx,biny); // weight per mass point
      }
      scale1fb = xsec*kfactor*1000.*filter*genWeight/(Float_t)sumGenWeightsHisto;
 
      int nPU = ( PUvar == "nVert") ? nVert : nTrueInt;
      int puBin = (int) hPU_r->GetXaxis()->FindBin(nPU);
      puWeight = hPU_r->GetBinContent(puBin);
    }

    if( applyJSON && isData && !golden.goodrun(run, lumi) ) isGolden=0;
    else isGolden=1;
    if( applyJSON && doSilver && isData && !silver.goodrun(run, lumi) ) isSilver=0;
    else isSilver=1;



    if( applySF ){
      /////////Add ISR scale factor////////// 
      /////////ISR for all MC now//////////////
      //NEW method with number of ISR jets
      float isr_err = 0.;
      if( nisrMatch == 1 ){
	weight_isr = 0.882;	    isr_err = 0.059;
      }else if( nisrMatch == 2 ){
	weight_isr = 0.792;	    isr_err = 0.104;
      }else if( nisrMatch == 3 ){ 
	weight_isr = 0.702;	    isr_err = 0.149;
      }else if( nisrMatch == 4 ){
	weight_isr = 0.648;	    isr_err = 0.176;
      }else if( nisrMatch == 5 ){
	weight_isr = 0.601;	    isr_err = 0.199;
      }else if( nisrMatch > 5 ){  
	weight_isr = 0.515;	    isr_err = 0.242;
      }

      double central = 1.;

      double central_UP = 1.;
      double central_DN = 1.;

      if( id<1000 )
	if( normFile == "" )
	  central = isr_average;
	else
	  central = isrAverageWeight;
      else{
	Int_t binx = h_isr->GetXaxis()->FindBin( (float)GenSusyMGluino );
	Int_t biny = h_isr->GetYaxis()->FindBin( (float)GenSusyMNeutralino );
	central = h_isr->GetBinContent(binx,biny);
	central_UP = h_isr_UP->GetBinContent(binx,biny);
	central_DN = h_isr_DN->GetBinContent(binx,biny);

	weight_btag_average          = h_btag->GetBinContent(binx,biny);        
	weight_btag_heavy_UP_average = h_btag_heavy_UP->GetBinContent(binx,biny);
	weight_btag_heavy_DN_average = h_btag_heavy_DN->GetBinContent(binx,biny);
	weight_btag_light_UP_average = h_btag_light_UP->GetBinContent(binx,biny);
	weight_btag_light_DN_average = h_btag_light_DN->GetBinContent(binx,biny);

      }

      //weight_isr /= central;
      weight_isr_av = central;
      weight_isr_UP_av = central_UP; 
      weight_isr_DN_av = central_DN;

      weight_isr_UP = weight_isr + isr_err;
      weight_isr_DN = weight_isr - isr_err;

      //Scale variations
      if( id>=1000 && id<2000 ){
	for(int v=0; v<110; ++v){
	  weight_scales[v] = 1.;
	  weight_scales_av[v] = 1.;
	  weight_scales[v] = LHEweight_wgt[v]/LHEweight_original;
	  //bins & binning are the same as for isr
	  Int_t binx = h_scales[v]->GetXaxis()->FindBin( (float)GenSusyMGluino );
	  Int_t biny = h_scales[v]->GetYaxis()->FindBin( (float)GenSusyMNeutralino );
	  double scale_av = h_scales[v]->GetBinContent(binx,biny);
	  weight_scales_av[v] = weight_scales[v] / scale_av;
	}
      }//finished scale weights


      int foundTop=0;
      int foundAntiTop=0;

      /////////Add Top pt scale factor//////////
      if( id >399 || id < 300 ) ;
      else if(nGenPart>0){
	for(int o=0; o<nGenPart; ++o){
	  if(GenPart_status[o] != 62) continue;
	  //Only apply weight to top or antitop
	  if( abs(GenPart_pdgId[o])==6 ){

	    float pt_withCutOff = GenPart_pt[o];
	    if(pt_withCutOff>400.0)
	      pt_withCutOff=400.0;

	    weight_toppt *= sqrt(exp( a + b * pt_withCutOff ) );

	    if( GenPart_pdgId[o]== 6 )  foundTop++;
	    if( GenPart_pdgId[o]==-6 )  foundAntiTop++;
	  }
	  else continue;
	}//end loop over objects

	//	if( foundTop==1 && foundAntiTop==1 ) //not averaging till we see that it has a big effect
	// weight_toppt /= average;
	//	else
	//  weight_toppt=1.0;
       
	//Only apply the weight if you find a pair
	if( !(foundTop==1 && foundAntiTop==1) )
	  weight_toppt=1.0;
      }



      /////////Add b-tagging scale factor//////////     
      bTagSFHelper->get_weight_btag(njet, jet_pt, jet_eta, jet_mcFlavour, jet_btagCSV, weight_btagsf, weight_btagsf_heavy_UP, weight_btagsf_heavy_DN, weight_btagsf_light_UP, weight_btagsf_light_DN, isFastSim);



      /////////Add lepton scale factor//////////
      if(nlep>0){
	Float_t uncert = 0.; //Place holder for total uncertainty
	Float_t central = 1; 
	Float_t uncert_UP = 1.; 	Float_t uncert_DN = 1.; 

	Float_t fast_central = 1; 
	Float_t fast_err = 0;


	for(int o=0; o<nlep; ++o){

	  Float_t err = 0.;

	  float pt = lep_pt[o];
	  float eta = fabs( lep_eta[o] );
	  int pdgId = abs( lep_pdgId[o] );

	  float pt_cutoff = std::max( 10.1, std::min( 100., double(pt) ) );
	  float pt_cutoff_eff = std::max( 5.1, std::min( 100., double(pt) ) );

	  //Electrons
	  if( pdgId == 11) {
	    Int_t binx = h_elSF->GetXaxis()->FindBin(pt_cutoff);
	    Int_t biny = h_elSF->GetYaxis()->FindBin(eta);
	    central = h_elSF->GetBinContent(binx,biny);
	    err  = h_elSF->GetBinError(binx,biny);
	    if (central > 1.2 || central < 0.8) 
	      std::cout<<"STRANGE: Electron with pT/eta of "<< pt <<"/"<< eta <<". SF is "<< central <<std::endl;

	    float central_trk = 1.;
	    Int_t binx_trk = h_elTrk->GetXaxis()->FindBin(  lep_eta[o] );
	    if( binx_trk>28 ) binx_trk = 28; //though we shouldn't really get many electrons with eta = 2.5
	    else if( binx_trk<1 ) binx_trk = 1;
	    central_trk = h_elTrk->GetBinContent( binx_trk,1  ); // y=pT range of 20-200

	    central *= central_trk;
	     
	    float trk_err = h_elTrk->GetBinError( binx_trk,1  ); // y=pT range of 20-200 

	    if( pt_cutoff < 20. )
	      err = sqrt( trk_err*trk_err + err*err + 0.03*0.03 );
	    else
	      err = sqrt( trk_err*trk_err + err*err  );

	    // uncert_UP +=  err;
	    // uncert_DN -=  err;
	    uncert_UP = central + err;
	    uncert_DN = central - err;

	    if( id>999 ){//FASTSIM SCALEFACTORS
	      Int_t fast_binx = h_fast_elSF->GetXaxis()->FindBin(pt_cutoff);
	      Int_t fast_biny = h_fast_elSF->GetYaxis()->FindBin(eta);
	      fast_central = h_fast_elSF->GetBinContent(fast_binx,fast_biny);
	      fast_err = h_fast_elSF->GetBinError(fast_binx,fast_biny);
	      fast_err = sqrt(fast_err*fast_err + 0.05*0.05); // 5% systematic uncertainty

	      if( fast_central > 1.3 || fast_central < 0.7 )
		std::cout << "Strange FastSim Electron with pT/eta of" << pt <<"/"<< eta <<". SF is "<< fast_central <<std::endl;

	      central   *= fast_central;
	      uncert_UP *= ( fast_central + fast_err );
	      uncert_DN *= ( fast_central - fast_err );
	      // uncert_UP += fast_err ;
	      //uncert_DN -= fast_err ;

	    }

	  } //else Muons
	  else if ( pdgId == 13) {
	    Int_t binx = h_muSF->GetXaxis()->FindBin(pt);
	    Int_t biny = h_muSF->GetYaxis()->FindBin(fabs(eta));

	    float central_trk = 1;
	    Int_t binx_trk = h_muTrk_hi->GetXaxis()->FindBin(  lep_eta[o] );
	    if( binx_trk>10 ) binx_trk = 10;
	    else if( binx_trk<1 ) binx_trk = 1;
	    central_trk = h_muTrk_hi->GetBinContent( binx_trk );


	    if ( binx >7 ) binx = 7; //overflow bin empty for the muons...
	    central = h_muSF->GetBinContent(binx,biny);

	    central *= central_trk;

	    err  = 0.03; //current recomendation is 3% //   err  = 0.014; // adding in quadrature 1% unc. on ID and 1% unc. on ISO

	    if (central > 1.3 || central < 0.7) 
	      std::cout<<"STRANGE: Muon with pT/eta of "<<pt<<"/"<< fabs(eta) <<". SF is "<< central <<std::endl;
	    //uncert_UP +=  err;
	    //uncert_DN -=  err;
	    uncert_UP = central + err;
	    uncert_DN = central - err;

	    if( id>999 ){ //FASTSIM SCALE FACTORS
	      Int_t fast_binx = h_fast_muSF->GetXaxis()->FindBin(pt);
	      Int_t fast_biny = h_fast_muSF->GetYaxis()->FindBin(fabs(eta));
	      fast_central = h_fast_muSF->GetBinContent(fast_binx,fast_biny);
	      fast_err  = h_fast_muSF->GetBinError(fast_binx,fast_biny);
	      if(lep_pt[0] < 20)
		fast_err= sqrt(fast_err*fast_err+ 0.03*0.03); // 3% systematic uncertainty
	      else
		fast_err= sqrt(fast_err*fast_err+ 0.01*0.01); // 1% systematic uncertainty

	      if( fast_central > 1.3 || fast_central < 0.7 )
		std::cout << "Strange FastSim Muon with pT/eta of" <<pt<<"/"<<eta<<". SF is "<< fast_central <<std::endl;

	      central   *= fast_central;

	      uncert_UP *= ( fast_central + fast_err );
	      uncert_DN *= ( fast_central - fast_err );

	      //uncert_UP += fast_err ;
	    //uncert_DN -= fast_err ;
	    }
	  }//done with one  electron/muon 
	  weight_lepsf    *= central;

	  weight_lepsf_UP *= uncert_UP;
	  weight_lepsf_DN *= uncert_DN;	  

 	  //weight_lepsf_UP += uncert_UP;
	  //weight_lepsf_DN += uncert_DN;	  
  	
	}//end of loop over objects

	//  weight_lepsf_UP += weight_lepsf;
	//  weight_lepsf_DN += weight_lepsf;	
	  
      }//end of lepton sf




      //0lepton efficiency correction
      if( nMuons10+nElectrons10+nPFLep5LowMT+nPFHad10LowMT == 0 ){ //otherwise we have caught the lep

	for(int o=0; o < ngenLep+ngenLepFromTau; ++o){

	  float central=1.;
	  float err=0;
	  float uncert=0;
	  float fast_central=1.;
	  float fast_err=0.;
	  float uncert_fast=0;

	  float pt, eta;
	  int  pdgId;
	  if(o < ngenLep){
	    pt = genLep_pt[o];
	    eta = genLep_eta[o];
	    pdgId = genLep_pdgId[o];
	  }else{
	    pt = genLepFromTau_pt[o-ngenLep];
	    eta = genLepFromTau_eta[o-ngenLep];
	    pdgId = genLepFromTau_pdgId[o-ngenLep] ;
	  }
	  // acceptance of lepton veto: pt>5, |eta| < 2.4
	  if( pt < 5.) continue;
	  if( fabs(eta) > 2.4) continue; 

	  float pt_cutoff = std::max( 10.1, std::min( 100., double(pt) ) );
	  float pt_cutoff_eff = std::max( 5.1, std::min( 100., double(pt) ) );
	  float veto_eff = 1.; //just need one veto eff, not one for full and fast sim

	  if ( abs( pdgId ) == 11) { 
	    Int_t binx = h_elSF->GetXaxis()->FindBin(pt_cutoff);
	    Int_t biny = h_elSF->GetYaxis()->FindBin( fabs(eta) );
	    central = h_elSF->GetBinContent(binx,biny);
	    uncert  = h_elSF->GetBinError(binx,biny);

	    Int_t binx_eff = h_eff_full_el->GetXaxis()->FindBin(pt_cutoff_eff);
	    Int_t biny_eff = h_eff_full_el->GetYaxis()->FindBin( fabs(eta) );
	    veto_eff = h_eff_full_el->GetBinContent(binx_eff,biny_eff);

	  } else if( abs( pdgId ) == 13) {
	    Int_t binx = h_muSF->GetXaxis()->FindBin(pt_cutoff);
	    Int_t biny = h_muSF->GetYaxis()->FindBin( fabs(eta) );
	    central = h_muSF->GetBinContent(binx,biny);
	    uncert = 0.03;// uncert  = 0.014;   // adding in quadrature 1% unc. on ID and 1% unc. on ISO

	    float central_trk = 1;
	    if( pt >= 10.){
	      Int_t binx_trk = h_muTrk_hi->GetXaxis()->FindBin(  lep_eta[o] );
	      if( binx_trk>10 ) binx_trk = 10;
	      else if( binx_trk<1 ) binx_trk = 1;
	      central_trk = h_muTrk_hi->GetBinContent( binx_trk );
	    }else if( pt < 10 ){
	      Int_t binx_trk = h_muTrk_lo->GetXaxis()->FindBin(  lep_eta[o] );
	      if( binx_trk>10 ) binx_trk = 10;
	      else if( binx_trk<1 ) binx_trk = 1;
	      central_trk = h_muTrk_lo->GetBinContent( binx_trk );
	    }
	    central *= central_trk;

	    Int_t binx_eff = h_eff_full_mu->GetXaxis()->FindBin( pt_cutoff_eff );
	    Int_t biny_eff = h_eff_full_mu->GetYaxis()->FindBin( fabs(eta) );
	    veto_eff = h_eff_full_mu->GetBinContent(binx_eff,biny_eff);
	  }

	  //Fast sim correction included, change weight and uncertainty
	  if( id>999 ){
	    if ( abs( pdgId ) == 11) {
	      Int_t fast_binx = h_fast_elSF->GetXaxis()->FindBin(pt_cutoff);
	      Int_t fast_biny = h_fast_elSF->GetYaxis()->FindBin( fabs(eta) );
	      fast_central = h_fast_elSF->GetBinContent(fast_binx,fast_biny);
	      //	      fast_err  = h_fast_elSF->GetBinError(fast_binx,fast_biny);
	      //	      uncert_fast = sqrt(fast_err*fast_err+ 0.05*0.05); // 5% systematic uncertainty
	      uncert_fast = 0.05; // 5% systematic uncertainty

	      veto_eff = h_eff_fast_el->GetBinContent(fast_binx,fast_biny);

	    }else if( abs( pdgId ) == 13) {

	      Int_t fast_binx = h_fast_muSF->GetXaxis()->FindBin(pt_cutoff);
	      Int_t fast_biny = h_fast_muSF->GetYaxis()->FindBin( fabs(eta) );
	      fast_central = h_fast_muSF->GetBinContent(fast_binx,fast_biny);
	      fast_err  = h_fast_muSF->GetBinError(fast_binx,fast_biny);
	      uncert_fast = 0.01;
	      if( pt_cutoff < 20.)
		uncert_fast= sqrt(fast_err*fast_err+ 0.03*0.03); // 3% systematic uncertainty
	      else
		uncert_fast= sqrt(fast_err*fast_err+ 0.01*0.01); // 1% systematic uncertainty
	      
	      veto_eff = h_eff_fast_mu->GetBinContent(fast_binx,fast_biny);

	    }//got the sf
	    float sf = central*fast_central;

	    float veto_eff_corr = veto_eff* sf; //corrected veto eff with the

	    weight_lepsf_0l *= (1. - veto_eff) / (1. - veto_eff_corr);

	    float unc = uncert + uncert_fast;
	    float veto_eff_cor_unc_UP = veto_eff_corr * ( 1. + unc);
	    float unc_UP_0l = (  1. - veto_eff_cor_unc_UP) / (1. - veto_eff_corr) -1. ;

	    weight_lepsf_0l_UP *= ( ( (1. - veto_eff) / (1. - veto_eff_corr)  + unc_UP_0l ));
	    weight_lepsf_0l_DN *= ( ( (1. - veto_eff) / (1. - veto_eff_corr)  - unc_UP_0l ));
	  
	  }else{//Only full sim correction, only uncertainty, not weight
	    
	    float sf = central;
	    // float veto_eff_corr = veto_eff* sf; //corrected veto eff with the 

	    float unc = uncert;
	    float veto_eff_unc_UP = veto_eff* (1. + unc);

	    float unc_UP_0l = (( 1. - veto_eff_unc_UP) / (1. - veto_eff))  -1.;
	    weight_lepsf_0l_UP *= ( 1. + unc_UP_0l);
	    weight_lepsf_0l_DN *= ( 1. - unc_UP_0l);
	  }


	}//end loop over genLep

      }//end 0Lep region
 

    }//end of if applySF







    if( i==(Long64_t)(nEventsTree-1))
      std::cout << "End of loop over tree entries" << std::endl;

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
    b12->Fill();
    b13->Fill();

    if( applySF ){   
      b14->Fill();
      b15->Fill();
      b16->Fill();
      b17->Fill();
      b18->Fill();
      b19->Fill();
      b20->Fill();
      b21->Fill();
      b22->Fill();
      b23->Fill();
      b24->Fill();
      b25->Fill();
      b26->Fill();
      b27->Fill();
      b28->Fill();
      b29->Fill();

      b30->Fill();
      b31->Fill();
      b32->Fill();
      b33->Fill();
      b34->Fill();
      b35->Fill();
      b36->Fill();
    }
    
  }//end loop over events
  //-------------------------------------------------------------


  if(applySF &&id > 1000 ){  
    for(int k=0;k<110;k++){
      delete h_scales[k];
      delete h_evt[k];
    }
    if( h_isr )
      delete h_isr; 
    delete h_counter;
  }
  //  if(bTagSFHelper) delete bTagSFHelper;
  delete chain; 
  hPU->Write();
  hPU_data->Write();
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
  delete out;
  return 0;
  
}

