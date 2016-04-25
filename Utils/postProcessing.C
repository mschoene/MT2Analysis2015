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

#include "BTagCalibrationStandalone.cc"
#include "goodrunClass.cc"

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

  double totalSumGenWeightsHisto, topAverageWeight;
 
  if( normFile!="" ){
    ifstream configuration( Form("%s", normFile.c_str()) );
    std::string line;
    while(std::getline(configuration,line)){
      istringstream ss(line);
      if((ss.str()[0])==(std::string)"#") continue;
      if(ss.str().empty()) continue;

      //      std::string normName;
      //      ss >> normName;
      ss >> totalSumGenWeightsHisto;
      ss >> topAverageWeight;

      std::cout<< "Got a line with " << std::endl;
      std::cout<< "sumGenWeightsHisto = " << totalSumGenWeightsHisto << std::endl;
      std::cout<< "weight_topPt_av    = " << topAverageWeight << std::endl;

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
  //Getting the lepton scale factor histograms/////////////////
  //Electrons//
  std::string filename = "kinematicBinSFele.root";
  TFile * f_ele = new TFile(filename.c_str() );
  if (!f_ele->IsOpen()) std::cout << " ERROR: Could not find scale factor file " << filename << std::endl; 
  //Uncomment for loose Id
  //TH2D* h_id = (TH2D*) f_ele->Get("CutBasedLoose");
  TH2D* h_id = (TH2D*) f_ele->Get("CutBasedVeto");
  TH2D* h_iso = (TH2D*) f_ele->Get("MiniIso0p1_vs_AbsEta");
  if (!h_id || !h_iso) std::cout << "ERROR: Could not find scale factor histogram"<< std::endl;
  TH2D* h_elSF = (TH2D*) h_id->Clone("h_elSF");
  h_elSF->SetDirectory(0);
  h_elSF->Multiply(h_iso);

  //Muons//
  std::string filenameID = "TnP_MuonID_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root";
  std::string filenameISO = "TnP_MuonID_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root";
  TFile * f1 = new TFile(filenameID.c_str() );
  TFile * f2 = new TFile(filenameISO.c_str() );
  if (!f1->IsOpen()) { std::cout<<" ERROR: Could not find ID scale factor file "<<filenameID<<std::endl; return 0;}
  if (!f2->IsOpen()) { std::cout<<"ERROR: Could not find ISO scale factor file "<<filenameISO<<std::endl; return 0;}
  TH2D* h_id_mu = (TH2D*) f1->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_tag_IsoMu20_pass");
  TH2D* h_iso_mu = (TH2D*) f2->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_PF_pass_&_tag_IsoMu20_pass");
  if (!h_id_mu || !h_iso_mu) { std::cout<<"ERROR: Could not find scale factor histogram"<<std::endl; return 0;}
  TH2D* h_muSF = (TH2D*) h_id_mu->Clone("h_muSF");
  h_muSF->SetDirectory(0);
  h_muSF->Multiply(h_iso_mu);

  f_ele->Close();  f1->Close(); f2->Close();
  delete f_ele; delete f1; delete f2;

  std::cout << std::endl;
  std::cout << "Using Loose Muon ID, MiniIso 0.2 lepton scale factors" << std::endl;
  std::cout << "Using Veto Electrons ID, MiniIso 0.1 lepton scale factors" << std::endl;
  std::cout << "Be aware that Veto Electrons are not suited for selecting Electrons." << std::endl;
  std::cout << std::endl;


  // if(id>=1000 && id <2000){
  std::cout << "Also loading the FastSim/FullSim Lepton scale factors" << std::endl;

  TFile * f_mu = new TFile("sf_mu_looseID_mini02.root" );
  TFile * f_el = new TFile("sf_el_vetoCB_mini01.root" );
  if(!f_mu->IsOpen()) { std::cout<<" ERROR: Could not find muon Fastsim scale factor file " <<std::endl; return 0;}
  if(!f_el->IsOpen()) { std::cout<<" ERROR: Could not find electron Fastsim scale factor file " <<std::endl; return 0;}

  TH3D* h_fast_muSF = (TH3D*) f_mu->Get("histo2D");
  TH3D* h_fast_elSF = (TH3D*) f_el->Get("histo2D");
  if(!h_fast_muSF || !h_fast_elSF ) {std::cout << " ERROR: Could not find the 3D histogram in your files " << std::endl; return 0;}

  h_fast_muSF->SetDirectory(0);
  h_fast_elSF->SetDirectory(0);

  f_mu->Close(); f_el->Close();
  delete f_mu; delete f_el;
  // }


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
  if (thisPUWeight) chain->SetBranchStatus("puWeight", 0);

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



  //float nTrueInt=0;
  int nTrueInt=0;
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
  Float_t weight_scales[110];
  Float_t weight_scales_av[110];


  ////// Lepton Efficiency SF
  Int_t nlep;
  chain->SetBranchAddress("nlep", &nlep);
  Float_t lep_pt[100];
  chain->SetBranchAddress("lep_pt", lep_pt);
  Float_t lep_eta[100];
  chain->SetBranchAddress("lep_eta", lep_eta);
  Int_t lep_pdgId[100];
  chain->SetBranchAddress("lep_pdgId", lep_pdgId);
  
  ////// b-tag SF
  Int_t njet;
  chain->SetBranchAddress("njet", &njet);
  Float_t jet_pt[100];
  chain->SetBranchAddress("jet_pt", jet_pt);
  Float_t jet_eta[100];
  chain->SetBranchAddress("jet_eta", jet_eta);
  Int_t jet_mcFlavour[100];
  if(!isData)
    chain->SetBranchAddress("jet_mcFlavour", jet_mcFlavour);
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
 
  TBranch* b14;
  TBranch* b15; 
  TBranch* b16;
  TBranch* b17;
  TBranch* b18;
  TBranch* b19; 
  TBranch* b20; 
  TBranch* b21; 
  TBranch* b22;
  TBranch* b23; 
  TBranch* b24; 
  TBranch* b25; 
  
  if( applySF ){
    b14 = clone->Branch("weight_btagsf"         , &weight_btagsf         , "weight_btagsf/F"         );
    b15 = clone->Branch("weight_btagsf_heavy_UP", &weight_btagsf_heavy_UP, "weight_btagsf_heavy_UP/F");
    b16 = clone->Branch("weight_btagsf_heavy_DN", &weight_btagsf_heavy_DN, "weight_btagsf_heavy_DN/F");
    b17 = clone->Branch("weight_btagsf_light_UP", &weight_btagsf_light_UP, "weight_btagsf_light_UP/F");
    b18 = clone->Branch("weight_btagsf_light_DN", &weight_btagsf_light_DN, "weight_btagsf_light_DN/F");
    b19 = clone->Branch("weight_lepsf", &weight_lepsf, "weight_lepsf/F");
    b20 = clone->Branch("weight_lepsf_UP", &weight_lepsf_UP, "weight_lepsf_UP/F");
    b21 = clone->Branch("weight_lepsf_DN", &weight_lepsf_DN, "weight_lepsf_DN/F");
    b22 = clone->Branch("weight_toppt", &weight_toppt, "weight_toppt/F");
    b23 = clone->Branch("weight_isr", &weight_isr, "weight_isr/F");
    b24 = clone->Branch("weight_scales", &weight_scales, "weight_scales[110]/F");
    b25 = clone->Branch("weight_scales_av", &weight_scales_av, "weight_scales_av[110]/F");
  }




  vector<TH2F*> h_scales;
  vector<TH2F*> h_evt;



  //ISR /// bins of size 5GeV for T2cc
  Float_t weight_isr_av;
  TH2F* h_isr = new TH2F("h_isr", "", 120*5, 12.5, 3012.5, 120*5, -12.5, 2987.5); h_isr->Sumw2();
  TH2F* h_counter = new TH2F("h_counter", "", 120*5, 12.5, 3012.5, 120*5, -12.5, 2987.5); h_counter->Sumw2();
  Int_t GenSusyMNeutralino;
  Int_t GenSusyMGluino;

  if( id>1000 && applySF ){

    for(int k=0;k<110;k++){
      std::string name = "var"+  std::to_string(k);
      TH2F* h_scales_temp = new TH2F(name.c_str(), "", 120*5, 12.5, 3012.5, 120*5, -12.5, 2987.5); h_scales_temp->Sumw2();
      h_scales.push_back(h_scales_temp);
 
      std::string name_evt = "evt"+  std::to_string(k); 
      TH2F* h_evt_temp = new TH2F(name_evt.c_str(), "", 120*5, 12.5, 3012.5, 120*5, -12.5, 2987.5); h_evt_temp->Sumw2();
      h_evt.push_back(h_evt_temp);
    }

    std::cout << "Entering first loop to determine the average ISR weights" << std::endl;

    chain->SetBranchAddress("GenSusyMNeutralino", &GenSusyMNeutralino);
    chain->SetBranchAddress("GenSusyMGluino", &GenSusyMGluino);
   
    
    for(  ULong64_t j = 0; j < nEventsTree; j++) {
    
      chain->GetEntry(j);
      weight_isr_av = 1.0;

      TLorentzVector s;
      if(nGenPart>0){
 	for(int o=0; o<nGenPart; ++o){
	  if(GenPart_status[o] != 62)   continue;
	  TLorentzVector s_;
	  s_.SetPtEtaPhiM(GenPart_pt[o], GenPart_eta[o], GenPart_phi[o], GenPart_mass[o]);
	  s+=s_;
	}  
	double pt_hard = s.Pt();
	if( pt_hard < 0. || pt_hard > 99999 )         continue;
	else if( pt_hard < 400. )  weight_isr_av = 1.0;
	else if( pt_hard < 600. )  weight_isr_av = 1.0-0.15;
	else if( pt_hard >= 600. ) weight_isr_av = 1.0-0.30;
      }

      h_isr->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino,  weight_isr_av);
      h_counter->Fill((float)GenSusyMGluino,(float)GenSusyMNeutralino,  1.0);

      
      if( id < 2000. ){
	for(int v=0; v<110; ++v){
	  weight_scales[v] = 1.;
	  weight_scales[v] = LHEweight_wgt[v]/LHEweight_original;

	  h_scales[v]->Fill((float)GenSusyMGluino , (float)GenSusyMNeutralino,  weight_scales[v]);
	  h_evt[v]->Fill((float)GenSusyMGluino,(float)GenSusyMNeutralino,  1.0);
	}
      }
    
    }//end or first loop over events for calculating the average ISR


    h_isr->Divide(h_counter);

    
    if( id<2000 ){
      for(int v=0; v<110; ++v){
	h_scales[v]->Divide(h_evt[v]);
      }
    }
    
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


  //Top pt reweighting 
  //Values from Run1 still valid for now, see here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  float a = 0.156;
  float b = -0.00137;
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




  std::cout << "Entering the final loop over the events" << std::endl;

  for( Long64_t i = 0; i < (Long64_t)nEventsTree; i++) {
  
    if( i==0)
      std::cout << "Start of loop over tree entries" << std::endl;


    weight_btagsf = 1.;
    weight_btagsf_heavy_UP = 1.;
    weight_btagsf_heavy_DN = 1.;
    weight_btagsf_light_UP = 1.;
    weight_btagsf_light_DN = 1.;
    weight_lepsf = 1.;
    weight_lepsf_UP = 1.;
    weight_lepsf_DN = 1.;
    weight_toppt = 1.;
    weight_isr=1.;
    
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
      scale1fb= xsec*kfactor*1000.*filter*genWeight/(Float_t)sumGenWeightsHisto;
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
      if( id < 1000 ) ;
      else{
	TLorentzVector s;
	if(nGenPart>0)
	  for(int o=0; o<nGenPart; ++o){
	    if(GenPart_status[o] != 62) continue;
	    TLorentzVector s_;
	    s_.SetPtEtaPhiM(GenPart_pt[o], GenPart_eta[o], GenPart_phi[o], GenPart_mass[o]);
	    s+=s_;	  
	  }
	float pt_hard = s.Pt();
	if( pt_hard < 0. )         continue;
	else if( pt_hard < 400. )  weight_isr = 1.0;
	else if( pt_hard < 600. )  weight_isr = 1.0-0.15;
	else if( pt_hard >= 600. ) weight_isr = 1.0-0.30;

	Int_t binx = h_isr->GetXaxis()->FindBin( (float)GenSusyMGluino );
	Int_t biny = h_isr->GetYaxis()->FindBin( (float)GenSusyMNeutralino );
	double central = h_isr->GetBinContent(binx,biny);

	weight_isr /= central;
	
	//Scale variations
	if( id<2000 ){

	  for(int v=0; v<110; ++v){
	    weight_scales[v] = 1.;
	    weight_scales_av[v] = 1.;
	    weight_scales[v] = LHEweight_wgt[v]/LHEweight_original;
	  
	    //bins & binning are the same as for isr
	    double scale_av = h_scales[v]->GetBinContent(binx,biny);
	    weight_scales_av[v] = weight_scales[v] / scale_av;

	  }
	}
	
      }//finished isr (and scale) weights


      int foundTop=0;
      int foundAntiTop=0;

      /////////Add Top pt scale factor//////////
      if( id >399 || id < 300 ) ;
      else if(nGenPart>0){
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
	  weight_toppt /= average;
	else
	  weight_toppt=1.0;

      }


      /////////Add b-tagging scale factor//////////     
      get_weight_btag(njet, jet_pt, jet_eta, jet_mcFlavour, jet_btagCSV, weight_btagsf, weight_btagsf_heavy_UP, weight_btagsf_heavy_DN, weight_btagsf_light_UP, weight_btagsf_light_DN, isFastSim);

      /////////Add lepton scale factor//////////
      if(nlep>0){
	Float_t uncert = 0; //Place holder for total uncertainty
	Float_t central = 1; 
	Float_t err = 0;

	Float_t fast_central = 1; 
	Float_t fast_err = 0;


	for(int o=0; o<nlep; ++o){
	  //Electrons
	  if (abs( lep_pdgId[o]) == 11) {
	    Int_t binx = h_elSF->GetXaxis()->FindBin(lep_pt[o]);
	    Int_t biny = h_elSF->GetYaxis()->FindBin(fabs(lep_eta[o]));
	    central = h_elSF->GetBinContent(binx,biny);
	    err  = h_elSF->GetBinError(binx,biny);
	    if (central > 1.2 || central < 0.8) 
	      std::cout<<"STRANGE: Electron with pT/eta of "<<lep_pt[o]<<"/"<<lep_eta[o]<<". SF is "<< central <<std::endl;

	    if( id>999 ){//FASTSIM SCALEFACTORS
	    Int_t fast_binx = h_fast_elSF->GetXaxis()->FindBin(lep_pt[o]);
	    Int_t fast_biny = h_fast_elSF->GetYaxis()->FindBin(fabs(lep_eta[o]));
	    fast_central = h_fast_elSF->GetBinContent(fast_binx,fast_biny);
	    fast_err  = h_fast_elSF->GetBinError(fast_binx,fast_biny);
	    fast_err= sqrt(fast_err*fast_err+ 0.05*0.05); // 5% systematic uncertainty

	    if( fast_central > 1.2 || fast_central < 0.8 )
	      std::cout << "Strange FastSim Electron with pT/eta of" <<lep_pt[o]<<"/"<<lep_eta[o]<<". SF is "<< fast_central <<std::endl;
	    central *= fast_central;
	    err = sqrt( fast_err*fast_err + err*err );
	    }

	  } //else Muons
	  else if (abs( lep_pdgId[o]) == 13) {
	    Int_t binx = h_muSF->GetXaxis()->FindBin(lep_pt[o]);
	    Int_t biny = h_muSF->GetYaxis()->FindBin(fabs(lep_eta[o]));
	    if ( binx >7 ) binx = 7; //overflow bin empty for the muons...
	    central = h_muSF->GetBinContent(binx,biny);
	    err  = 0.014; // adding in quadrature 1% unc. on ID and 1% unc. on ISO
	    if (central > 1.3 || central < 0.7) 
	      std::cout<<"STRANGE: Muon with pT/eta of "<<lep_pt[o]<<"/"<< fabs(lep_eta[o]) <<". SF is "<< central <<std::endl;
	  
	    if( id>999 ){ //FASTSIM SCALE FACTORS
	    Int_t fast_binx = h_fast_muSF->GetXaxis()->FindBin(lep_pt[o]);
	    Int_t fast_biny = h_fast_muSF->GetYaxis()->FindBin(fabs(lep_eta[o]));
	    fast_central = h_fast_muSF->GetBinContent(fast_binx,fast_biny);
	    fast_err  = h_fast_muSF->GetBinError(fast_binx,fast_biny);
	    if(lep_pt[0]>20)
	      fast_err= sqrt(fast_err*fast_err+ 0.03*0.03); // 5% systematic uncertainty
	    else
	      fast_err= sqrt(fast_err*fast_err+ 0.01*0.01); // 5% systematic uncertainty

	    if( fast_central > 1.3 || fast_central < 0.7 )
	      std::cout << "Strange FastSim Muon with pT/eta of" <<lep_pt[o]<<"/"<<lep_eta[o]<<". SF is "<< fast_central <<std::endl;

	    central *= fast_central;
	    err = sqrt( fast_err*fast_err+ err*err );
	    }
	  } 
	  weight_lepsf *= central;
	  //uncertainties are supposed to be summed up linearly for the lepton SF
	  uncert += err; 	
	}//end of loop over objects
    
	weight_lepsf_UP = central + uncert; 
	weight_lepsf_DN = central - uncert;
	 
      }//end of lepton sf
    }//end of if applySF

    if( i==(nEventsTree-1))
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
    }
    
  }//end loop over events
  //-------------------------------------------------------------


  if(applySF &&id > 1000 ){  
    for(int k=0;k<110;k++){
      delete h_scales[k];
      delete h_evt[k];
    }
    delete h_isr; delete h_counter;
  }

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


