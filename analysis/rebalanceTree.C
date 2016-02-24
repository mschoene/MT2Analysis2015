#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMinuit.h"
#include "TMath.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "rooDoubleCB.h"
#include <vector>

Double_t sigmaSoft = 20.; 

const Int_t NJ=50;
Int_t Nj =6;
Double_t jPtReco [NJ];
Double_t jPhiReco[NJ];
Double_t jEtaReco[NJ];
Double_t met, metphi;
Double_t softPxReco, softPyReco;

Int_t jPtIndex[NJ], jEtaIndex[NJ];

Int_t NPUj=6;
Double_t jPtPU [NJ];
Double_t jPhiPU[NJ];
Double_t jEtaPU[NJ];


// open file with response functions
//PtBinEdges = 0, 20, 30, 50, 80, 120, 170, 230, 300, 380, 470, 570, 680, 800, 1000, 1300, 1700, 2200, 2800, 3500, 4300, 5200, 6500 GeV
//EtaBinEdges = 0, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.3, 2.8, 3.2, 4.1, 5.0
double PtBinEdges[23] = {0, 20, 30, 50, 80, 120, 170, 230, 300, 380, 470, 570, 680, 800, 1000, 1300, 1700, 2200, 2800, 3500, 4300, 5200, 6500};
double EtaBinEdges[12] = {0, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.3, 2.8, 3.2, 4.1, 5.0};
RooWorkspace *w[22][12];
TFile *fTemplates = new TFile("responseTemplates.root");

class EventKey {
public:
  EventKey(unsigned int input_run=0, unsigned int input_lumi=0, unsigned long long input_evt=0) : 
    run_(input_run), lumi_(input_lumi), evt_(input_evt){;}

  unsigned int run() const {return run_;}
  unsigned int lumi() const {return lumi_;}
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
  unsigned int run_;
  unsigned int lumi_;
  unsigned long long evt_;

};


void getTemplates(){
  for (int pt=0; pt<22; pt++){
    for (int eta=0; eta<12; eta++){
      TString temp = TString::Format("wsp_h_tot_JetAll_ResponsePt_Pt%d_Eta%d", pt, eta);
      if ( fTemplates->GetListOfKeys()->Contains(temp) ){
  	w[pt][eta] = (RooWorkspace*)fTemplates->Get(temp);
      }
    }
  }
}

void updateIndex(int j, double x){
  int p=0, e=0;
  double pt = jPtReco[j]/x;
  //for (p=0; p<22; p++){
  for (p=0; p<16; p++){  // only up to bin Pt16, empty for higher pt
    if (PtBinEdges[p+1]>pt)
      break;
  }	  

  if(pt < 300.0){ // all eta bins
    for (e=0; e<11; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReco[j]))
	break;
    }
  }
  else if(pt < 570.0){ // up to bin Eta9
    for (e=0; e<9; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReco[j]))
	break;
    }
  }
  else if(pt < 1000.0){ // up to bin Eta8
    for (e=0; e<8; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReco[j]))
	break;
    }
  }
  else if(pt < 1300.0){ // up to bin Eta7
    for (e=0; e<7; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReco[j]))
  	  break;
    }
  }
  else if(pt < 1700.0){ // up to bin Eta6
    for (e=0; e<6; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReco[j]))
	break;
    }
  }
  else{ // up to bin Eta2
    for (e=0; e<2; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReco[j]))
	break;
    }
  }
  jPtIndex [j]=p;
  jEtaIndex[j]=e;
}

Double_t getResponseProb(Double_t x, Int_t p, Int_t e){
  w[p][e]->var("x")->setVal(x);
  return w[p][e]->pdf("cb")->getVal();
}

Double_t getResponseProbGaus(Double_t x, Int_t p, Int_t e){
  double mean  = w[p][e]->var("mean" )->getVal();
  double sigma = w[p][e]->var("sigma")->getVal();
  return TMath::Gaus(x, mean, sigma);
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  // par[j] = 1/c[j], ptrue = 1/c*ptreco
  Double_t logl = 0;
  Double_t softPxTrue=0, softPyTrue=0;
  for (int j=0; j<Nj; j++) {
    //logl += TMath::Log(TMath::Gaus(par[j], 1.0, 0.1)); // example response function for trial
    updateIndex(j,par[j]); // the binning is gen Pt, find the corresponding bin
    //logl += TMath::Log(getResponseProb(par[j], jPtIndex[j], jEtaIndex[j])); // take appropriate response function
    logl += TMath::Log(getResponseProbGaus(par[j], jPtIndex[j], jEtaIndex[j])); // take appropriate response function (only gaussian)
    softPxTrue -= jPtReco[j]*TMath::Cos(jPhiReco[j])/par[j];
    softPyTrue -= jPtReco[j]*TMath::Sin(jPhiReco[j])/par[j];
  }
  //logl +=  TMath::Log(TMath::Gaus((softPxReco-softPxTrue)/sigmaSoft , 0.0, 1.0));
  //logl +=  TMath::Log(TMath::Gaus((softPyReco-softPyTrue)/sigmaSoft , 0.0, 1.0));
  logl +=  -TMath::Power((softPxReco-softPxTrue)/sigmaSoft,2)/2; // log(gaus)
  logl +=  -TMath::Power((softPyReco-softPyTrue)/sigmaSoft,2)/2;

  f = -logl;
}


void rebalanceTree(){

  TStopwatch tt;
  tt.Start();
  cout << "starting" << endl;

  getTemplates();

  //Get input tree
  TFile *oldfile = new TFile("/shome/pandolf/QCDtreesRS/smearTree_MG_QCD_ht1000to1500.root");
  TTree *oldtree = (TTree*)oldfile->Get("SmearTree");
  Long64_t nentries = oldtree->GetEntries();
  
  cout << "In input tree, nentries = " << nentries << endl;

   Int_t           run;
   Int_t           lumi;
   ULong64_t       evt;
   oldtree->SetBranchAddress("evt_run",&run);
   oldtree->SetBranchAddress("evt_lumi",&lumi);
   oldtree->SetBranchAddress("evt_event",&evt);

  // vector<float>   *jet_puId = 0;
  // vector<float>   *jet_pt   = 0;
  // vector<float>   *jet_eta  = 0;
  // vector<float>   *jet_phi  = 0;
  vector<float>   *rsjet_pt = new vector<float>();
  vector<float>   *jet_puId = new vector<float>();
  vector<float>   *jet_pt   = new vector<float>();
  vector<float>   *jet_eta  = new vector<float>();
  vector<float>   *jet_phi  = new vector<float>();
  Float_t         met_pt;
  Float_t         met_phi;

  TBranch        *b_before_jet_pt;   //!
  TBranch        *b_rs_jet_pt;   //!
  TBranch        *b_rs_jet_puId;     //!
  TBranch        *b_rs_jet_eta;      //!
  TBranch        *b_rs_jet_phi;      //!
  TBranch        *b_before_met_pt;   //!
  TBranch        *b_before_met_phi;  //!

  oldtree->SetBranchAddress("before_jet_pt" , &jet_pt  , &b_before_jet_pt );
  oldtree->SetBranchAddress("rs_jet_puId"   , &jet_puId, &b_rs_jet_puId   );
  oldtree->SetBranchAddress("rs_jet_pt"     , &rsjet_pt, &b_rs_jet_pt    );
  oldtree->SetBranchAddress("rs_jet_eta"    , &jet_eta , &b_rs_jet_eta    );
  oldtree->SetBranchAddress("rs_jet_phi"    , &jet_phi , &b_rs_jet_phi    );
  oldtree->SetBranchAddress("before_met_pt" , &met_pt  , &b_before_met_pt );
  oldtree->SetBranchAddress("before_met_phi", &met_phi , &b_before_met_phi);

 
  //Create a new file + a clone of old tree in new file
  TFile *newfile;
  TTree *newtree;
  newfile = new TFile("rebalancedTree.root","recreate");
  newtree = oldtree->CloneTree(0);

  //vector<float>   *new_jet_pt  = 0;
  vector<float>   *new_jet_pt  =  new vector<float>();
  Float_t          new_met_pt;
  Float_t          new_met_phi;
  Float_t          new2_met_pt;
  Float_t          new2_met_phi;
  Float_t          new_true_softpt;
  Int_t            errflag;
  Int_t            nTries;
  Float_t          realTime;
  Float_t          cpuTime;

  newtree->Branch("rebalanced_jet_pt" , &new_jet_pt     );
  newtree->Branch("reb_brunos_met_pt" , &new_met_pt     );
  newtree->Branch("reb_brunos_met_phi", &new_met_phi    );
  newtree->Branch("reb_jasons_met_pt" , &new2_met_pt    );
  newtree->Branch("reb_jasons_met_phi", &new2_met_phi   );
  newtree->Branch("rebalanced_softpt" , &new_true_softpt);
  newtree->Branch("errflag"           , &errflag        );
  newtree->Branch("nTries"            , &nTries         );
  newtree->Branch("realTime"          , &realTime       );
  newtree->Branch("cpuTime"           , &cpuTime        );

  oldfile->cd();
  
  //Create set where we store list of event keys
  std::set<EventKey> previousEvents;
  
  cout << "starting loop over tree events" << endl;
  

  //nentries = 5000;
  for (int i=0;i<nentries; i++) {

    TStopwatch t;
    t.Start();
    
    oldtree->GetEntry(i);
    
    if(i%100000==0) {
      cout << "Processing event: " << i  << endl;
      newtree->AutoSave();
    }
    
    EventKey newEvent(run,lumi,evt);
    bool wasRebalanced = !previousEvents.insert(newEvent).second;

    if (wasRebalanced) continue;
    

    // Nj=0;
    // NPUj=0;
    // for (unsigned int j=0; j<jet_pt->size(); j++){
    //   if ( jet_pt->at(j)>100 || jet_puId->at(j)==1 ){
    // 	jPtReco [Nj] = jet_pt ->at(j);
    // 	jEtaReco[Nj] = jet_eta->at(j);
    // 	jPhiReco[Nj] = jet_phi->at(j);
    // 	Nj++;
    //   }
    //   else{
    // 	jPtPU [NPUj] = jet_pt ->at(j);
    // 	jEtaPU[NPUj] = jet_eta->at(j);
    // 	jPhiPU[NPUj] = jet_phi->at(j);
    // 	NPUj++;
    //   }
    // }

    // PU jets are demoted to last, order of rest is according to original pt
    std::vector<float> pujetspt;
    for (unsigned int j=jet_puId->size()-1; j>0; j--) { // first get a pt list of pu jets (reverse order)
      if ( jet_puId->at(j)==0 )
	pujetspt.push_back(rsjet_pt->at(j));
      else 
	break;
    }
    //cout << "pujets size = " << pujetspt.size() << endl;
    Nj=0;
    NPUj=0;
    int njets   = jet_pt->size();
    int npujets = pujetspt.size();
    for (unsigned int j=0; j<njets; j++){
      //cout << j << " " << Nj  << " " << NPUj << endl;
      if ( npujets>NPUj && jet_pt->at(j) == pujetspt.at(npujets-NPUj-1) ){
	//cout << "got pu jet " << j << " " << NPUj << " " << pujetspt.size()-NPUj-1 << endl;
    	jPtPU [NPUj] = jet_pt ->at(j);
    	jEtaPU[NPUj] = jet_eta->at(njets-npujets+NPUj); // pu jets are last in rs collection
    	jPhiPU[NPUj] = jet_phi->at(njets-npujets+NPUj);
    	NPUj++;
	//cout << "pu jet filled" << endl;
      }
      else{
	//cout << "got jet " << j << " " << Nj << endl;
    	jPtReco [Nj] = jet_pt ->at(j);
    	jEtaReco[Nj] = jet_eta->at(j-NPUj); // jet j-NPUj in rs collection
    	jPhiReco[Nj] = jet_phi->at(j-NPUj);
    	Nj++;
	//cout << "jet filled" << endl;
      }
      //cout << j << " " << Nj  << " " << NPUj << endl;
    }


    
    met = met_pt;
    metphi = met_phi;
    softPxReco=-met*TMath::Cos(metphi);
    softPyReco=-met*TMath::Sin(metphi);
    for (int j=0; j<Nj; j++) {
      softPxReco -= jPtReco[j]*TMath::Cos(jPhiReco[j]);
      softPyReco -= jPtReco[j]*TMath::Sin(jPhiReco[j]);
    }
    for (int j=0; j<NPUj; j++) {
      softPxReco -= jPtPU[j]*TMath::Cos(jPhiPU[j]);
      softPyReco -= jPtPU[j]*TMath::Sin(jPhiPU[j]);
    }

     // for (int j=0; j<Nj+NPUj; j++) {
     //   if (j<Nj)
     // 	cout << "jet " << j << " with reco pt = " << jPtReco[j] << " and  eta = " << jEtaReco[j] << endl;
     //   else 
     // 	cout << "PU jet " << j << " with reco pt = " << jPtPU[j-Nj] << " and  eta = " << jEtaPU[j-Nj] << endl;
     // }
     // cout << "met: "        << met        << "; metphi: "     << metphi     << endl
     //  	 << "softPxreco: " << softPxReco << "; softPyreco: " << softPyReco << endl;

    
    TMinuit *gMinuit = new TMinuit(Nj);
    gMinuit->SetFCN(fcn);
    
    Double_t arglist[10];
    Int_t ierflg = 0;
    
    //arglist[0] = 0; // quiet printout (minimum output)
    arglist[0] = -1; // quietest printout (no output)
    gMinuit->mnexcm("SET PRI", arglist, 1, ierflg);
    
    arglist[0] = 0.5; // for likelihood
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
    
    
    for (int j=0; j<Nj; j++) 
      gMinuit->mnparm(j, TString::Format("c%d",j), 1., 0.01, 0, 10, ierflg);
    
    
	 
    nTries = 0;
    arglist[0] = 10000; // max calls
    arglist[1] = 1.0;   // tolerance
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    nTries++;

    errflag = ierflg;
    
    if (ierflg != 0){
      cout << "*******ERROR: MIGRAD minimization failed: " << ierflg << endl;
      cout << "attempting relaxed tolerance threshold" << endl;
      arglist[1] = 10.0;   // tolerance
      gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
      nTries++;
      errflag = ierflg;
    }
      
    if (ierflg != 0){
      cout << "*******ERROR: MIGRAD minimization failed again: " << ierflg << endl;
      cout << "attempting MINOS minimization" << endl;
      arglist[0] = 5000; // max calls
      arglist[1] = 0;
      arglist[2] = 1;
      gMinuit->mnexcm("MINOS", arglist, 3, ierflg);
      nTries++;
      errflag = ierflg;
    }

    if (ierflg != 0)
      cout << "*******ERROR: MINOS minimization failed: " << ierflg << endl;
      

    Double_t ic [NJ];
    Double_t eic[NJ];
    
    Double_t softPxTrue=0, softPyTrue=0;
    for (int j=0; j<Nj; j++) {
      gMinuit->GetParameter(j, ic[j], eic[j]);
      //ic[j] = jasonsC[j];
      //cout << "1/c["<<j<<"] = " << ic[j] << endl;
    
      softPxTrue -= jPtReco[j]*TMath::Cos(jPhiReco[j])/ic[j];
      softPyTrue -= jPtReco[j]*TMath::Sin(jPhiReco[j])/ic[j];
      //cout << "jet " << j << " with reco pt = " << jPtReco[j] << "(" << jPtIndex[j] << ") and  eta = " << jEtaReco[j] << "(" << jEtaIndex[j] << ") --> true pt = " << jPtReco[j] << "/" << ic[j] << " = " << jPtReco[j]/ic[j] << endl;
    }
    // cout << "met: "        << met        << "; metphi: "     << metphi     << endl
    // 	 << "softPxreco: " << softPxReco << "; softPyreco: " << softPyReco << endl
    // 	 << "softPxTrue: " << softPxTrue << "; softPyTrue: " << softPyTrue << endl;

    // my rebalanced met. uses softPtTrue (balanced with rebalanced jets) and PU jets
    Double_t metx_re=0., mety_re=0., met_re, metphi_re;
    for (int j=0; j<NPUj; j++) {
      metx_re -= jPtPU[j]*TMath::Cos(jPhiPU[j]);
      mety_re -= jPtPU[j]*TMath::Sin(jPhiPU[j]);
    }
    met_re = sqrt(metx_re*metx_re+mety_re*mety_re);
    metphi_re = atan(mety_re/metx_re);
    //cout << "bruno's rebalanced met: " << met_re << "; metphi: " << metphi_re << endl;
    
    // jason's rebalanced met. uses softPtReco and ignores PU jets
    Double_t metx_re2 = softPxTrue - softPxReco;
    Double_t mety_re2 = softPyTrue - softPyReco;
    Double_t met_re2 = sqrt(metx_re2*metx_re2+mety_re2*mety_re2);
    Double_t metphi_re2 = atan(mety_re2/metx_re2);
    //cout << "jason's rebalanced met: " << met_re2 << "; metphi: " << metphi_re2 << endl;
    
    
    
    new_jet_pt ->clear();
    
    // int nj=0;
    // for (unsigned int j=0; j<jet_pt->size(); j++){
    //   if ( jet_pt->at(j)>100 || jet_puId->at(j)==1 ){
    // 	new_jet_pt ->push_back(jet_pt ->at(j)/ic[nj]);
    // 	nj++;
    //   }
    //   else{
    // 	new_jet_pt ->push_back(jet_pt ->at(j));
    //   }
    // }
    for (int j=0; j<Nj; j++)
      new_jet_pt ->push_back(jPtReco[j]/ic[j]);
    for (int j=0; j<NPUj; j++) 
      new_jet_pt ->push_back(jPtPU[j]);
   
    new_met_pt      = met_re;
    new_met_phi     = metphi_re;
    new2_met_pt     = met_re2;
    new2_met_phi    = metphi_re2;
    new_true_softpt = sqrt(softPxTrue*softPxTrue+softPyTrue*softPyTrue);
    
    t.Stop();
    realTime = t.RealTime();
    cpuTime  = t.CpuTime();

    newtree->Fill();

    gMinuit->Delete();

  }

  newfile->cd();
  newtree->AutoSave();

  tt.Stop();
  tt.Print();
}
