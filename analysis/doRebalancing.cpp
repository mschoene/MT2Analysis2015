#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"
#include "TMinuit.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "interface/rooDoubleCB.h"
#include "TStopwatch.h"

#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"

int ijob=0, Njobs=1;

Double_t sigmaSoft = 20.; 

const Int_t NJ=50;
Int_t Nj, NPUj;
Double_t jPtReco [NJ];
Double_t jPhiReco[NJ];
Double_t jEtaReco[NJ];
Double_t jPtPU [NJ];
Double_t jPhiPU[NJ];
Double_t jEtaPU[NJ];
Int_t jPtIndex[NJ], jEtaIndex[NJ];
Double_t softPxReco, softPyReco;


// open file with response functions
double PtBinEdges[23] = {0, 20, 30, 50, 80, 120, 170, 230, 300, 380, 470, 570, 680, 800, 1000, 1300, 1700, 2200, 2800, 3500, 4300, 5200, 6500};
double EtaBinEdges[12] = {0, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.3, 2.8, 3.2, 4.1, 5.0};
RooWorkspace *w[22][12];
TFile *fTemplates = new TFile("/shome/casal/templatesRandS/responseTemplates.root");

void rebalance( const MT2Sample& sample, MT2Analysis<MT2EstimateTree>* anaTree );
void getTemplates();
void updateIndex(int j, double x);
Double_t getResponseProb(Double_t x, Int_t p, Int_t e);
Double_t getResponseProbGaus(Double_t x, Int_t p, Int_t e);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

int main( int argc, char* argv[] ) {



  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|             rebalancing events            |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;

  if( argc<2 ) {
    std::cout << "USAGE: ./doRebalancing configFileName [data/MC/all] [sampleID] [job_i] [Njobs]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  
  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool onlyData = false;
  bool onlyMC   = false;
  if( argc > 2 ) {
    std::string dataMC(argv[2]);
    if( dataMC=="data" ) onlyData = true;
    else if( dataMC=="MC" || dataMC=="mc" ) onlyMC = true;
  }

  if( onlyData ) {
    std::cout << "-> Will run only on data." << std::endl;
  } else if( onlyMC ) {
    std::cout << "-> Will run only on MC." << std::endl;
  } else {
    std::cout << "-> Will run only on both data and MC." << std::endl;
  }
 
  int sampleID = -1;
  if( argc > 3 ) {
    sampleID = atoi(argv[3]);
    std::cout << "-> Will run over sample ID" << sampleID << std::endl;
  }

  if( argc > 5 ) {
    ijob  = atoi(argv[4]);
    Njobs = atoi(argv[5]);
    std::cout << "-> Will run job " << ijob << " out of " << Njobs << std::endl;
  }



  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  std::string outputdir = cfg.getEventYieldDir() + "/rebalancedTrees"; 
  system(Form("mkdir -p %s", outputdir.c_str()));


  // load response templates on memory
  getTemplates();

  if( cfg.useMC() && !onlyData ) { // run on MC

    std::string samplesFile = "../samples/samples_" + cfg.mcSamples() + ".dat";
    
    std::vector<MT2Sample> samples_qcd;

    if (sampleID==-1) // do all qcd samples
      samples_qcd = MT2Sample::loadSamples(samplesFile, 100, 199);
    else 
      samples_qcd = MT2Sample::loadSamples(samplesFile, sampleID, sampleID);

    MT2Analysis<MT2EstimateTree>* qcdCRtree = new MT2Analysis<MT2EstimateTree>( "qcdRebalancedTree", "13TeV_inclusive" );
    MT2EstimateTree::addVar( qcdCRtree, "jet1_pt" );
    MT2EstimateTree::addVar( qcdCRtree, "jet2_pt" );
    MT2EstimateTree::addVar( qcdCRtree, "metphi"  );
    MT2EstimateTree::addVar( qcdCRtree, "errflag" );
    MT2EstimateTree::addVar( qcdCRtree, "ntries"  );
    MT2EstimateTree::addVar( qcdCRtree, "cputime" );
    MT2EstimateTree::addVar( qcdCRtree, "realtime" );
    MT2EstimateTree::addVar( qcdCRtree, "reco_soft_pt" );
    MT2EstimateTree::addVar( qcdCRtree, "reco_soft_phi" );
    MT2EstimateTree::addVar( qcdCRtree, "true_soft_pt" );
    MT2EstimateTree::addVar( qcdCRtree, "true_soft_phi" );
    MT2EstimateTree::addVar( qcdCRtree, "reb_jasons_met_pt"  );
    MT2EstimateTree::addVar( qcdCRtree, "reb_jasons_met_phi" );
    MT2EstimateTree::addVar( qcdCRtree, "reb_brunos_met_pt"  );
    MT2EstimateTree::addVar( qcdCRtree, "reb_brunos_met_phi" );
    MT2EstimateTree::addVector( qcdCRtree, "before_jet_pt"  );
    MT2EstimateTree::addVector( qcdCRtree, "reb_jet_pt"  );
    MT2EstimateTree::addVector( qcdCRtree, "jet_eta" );
    MT2EstimateTree::addVector( qcdCRtree, "jet_phi" );
    MT2EstimateTree::addVector( qcdCRtree, "jet_mass" );
    MT2EstimateTree::addVector( qcdCRtree, "jet_puId" );
    MT2EstimateTree::addVector( qcdCRtree, "jet_id" );
    MT2EstimateTree::addVector( qcdCRtree, "jet_btagCSV" );
    
    
    for( unsigned i=0; i<samples_qcd.size(); ++i ) 
      rebalance( samples_qcd[i], qcdCRtree );
    

   
    TString mcFile = outputdir + "/mc";
    if (sampleID!=-1)
      mcFile += TString::Format("_id%d_job%dof%d", sampleID, ijob, Njobs);
    mcFile += ".root";
    qcdCRtree->writeToFile( mcFile.Data(), "RECREATE" );

  }


  if( !(cfg.dummyAnalysis()) && cfg.dataSamples()!="" && !onlyMC  ) {  // run on data

    //to be implemented

  }



  return 0;

}


void rebalance( const MT2Sample& sample, MT2Analysis<MT2EstimateTree>* anaTree ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  int nentries = tree->GetEntries();

  std::cout << "-> Loaded tree: it has " << nentries << " entries." << std::endl
	    << "   will run " << Njobs << " jobs (roughly " << nentries/Njobs << " events per job)" << std::endl;
  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  //nentries = 1000;

  //for( int iEntry=0; iEntry<nentries; ++iEntry ) {
  for( int iEntry=0+ijob; iEntry<nentries; iEntry+=Njobs ) {

    TStopwatch t;
    t.Start();

    if( iEntry/Njobs % 50000 == 0 ) std::cout << "    Entry: " << iEntry/Njobs << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    if( myTree.isData ) {
      //if (  myTree.isGolden == 0 ) continue;
      if ( !myTree.passFilters() ) continue;
    }
    

    // remove leptons
    if( (myTree.nElectrons10 + mt2tree.nMuons10 + mt2tree.nPFLep5LowMT + mt2tree.nPFHad10LowMT)>0 )
      continue;

    float minMTBmet = myTree.minMTBMet;
    float met_pt    = myTree.met_pt;
    float met_phi   = myTree.met_phi;
    int njets       = myTree.nJet30;
    int nbjets      = myTree.nBJet20;    
    float mt2       = (njets>1) ? myTree.mt2 : myTree.jet1_pt;
    float ht        = myTree.ht;

    Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;//*cfg.lumi(); 
    //Double_t weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*myTree.weight_lepsf*myTree.weight_btagsf*myTree.weight_toppt;
    
    
    //float myht = njets==1 ? 201. : ht; // let everything (mt2=ht>40) pass for monojet
    MT2EstimateTree* thisTree = anaTree->get( ht, njets, nbjets, minMTBmet, mt2 );
    if( thisTree==0 ) continue;


    std::vector<float> all_jet_pt     ;
    std::vector<float> all_jet_eta    ;
    std::vector<float> all_jet_phi    ;
    std::vector<float> all_jet_mass   ;
    std::vector<float> all_jet_puId   ;
    std::vector<float> all_jet_id     ;
    std::vector<float> all_jet_btagCSV;
    std::vector<float> reb_jet_pt     ;
    
    Nj=0;
    NPUj=0;
    for (int j=0; j< myTree.njet; j++) {
      all_jet_pt     .push_back(myTree.jet_pt     [j]);
      all_jet_eta    .push_back(myTree.jet_eta    [j]);
      all_jet_phi    .push_back(myTree.jet_phi    [j]);      
      all_jet_mass   .push_back(myTree.jet_mass   [j]);
      all_jet_puId   .push_back(myTree.jet_puId   [j]);
      all_jet_id     .push_back(myTree.jet_id     [j]);
      all_jet_btagCSV.push_back(myTree.jet_btagCSV[j]);
      if (myTree.jet_pt[j]>100 || myTree.jet_puId[j]==1){
    	jPtReco [Nj] = myTree.jet_pt [j];
    	jEtaReco[Nj] = myTree.jet_eta[j];
    	jPhiReco[Nj] = myTree.jet_phi[j];
    	Nj++;
      }
      else{
    	jPtPU [NPUj] = myTree.jet_pt [j];
    	jEtaPU[NPUj] = myTree.jet_eta[j];
    	jPhiPU[NPUj] = myTree.jet_phi[j];
    	NPUj++;
      }
    }
    
    softPxReco=-met_pt*TMath::Cos(met_phi);
    softPyReco=-met_pt*TMath::Sin(met_phi);
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
    // 	std::cout << "jet " << j << " with reco pt = " << jPtReco[j] << " and  eta = " << jEtaReco[j] << std::endl;
    //   else 
    // 	std::cout << "PU jet " << j << " with reco pt = " << jPtPU[j-Nj] << " and  eta = " << jEtaPU[j-Nj] << std::endl;
    // }
    // std::cout << "met: "        << met_pt        << "; metphi: "     << met_phi     << std::endl
    // 	      << "softPxreco: " << softPxReco << "; softPyreco: " << softPyReco << std::endl;

    TMinuit *gMinuit = new TMinuit(Nj);
    gMinuit->SetFCN(fcn);
    
    Double_t arglist[10];
    Int_t ierflg = 0;
    
    //arglist[0] = 0; // quiet printout (minimum output)
    arglist[0] = -1; // quietest printout (no output... [although that's not really true])
    gMinuit->mnexcm("SET PRI", arglist, 1, ierflg);
    
    arglist[0] = 0.5; // for likelihood
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
    
    
    for (int j=0; j<Nj; j++) 
      gMinuit->mnparm(j, TString::Format("c%d",j), 1., 0.01, 0, 10, ierflg);
    
	 
    int nTries = 0;
    arglist[0] = 10000; // max calls
    arglist[1] = 1.0;   // tolerance. minimization stops when estimated vertical distance to the minimum (EDM) is les than 0.001*[tolerance]*UP  (UP=ERRordef set to 0.5 above for LL)
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    nTries++;  // 1

    if (ierflg != 0){
      //std::cout << "*******ERROR: MIGRAD minimization failed: " << ierflg << std::endl;
      //std::cout << "attempting relaxed tolerance threshold" << std::endl;
      arglist[1] = 5.0;   // tolerance
      gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
      nTries++;  // 2
    }
      
    if (ierflg != 0){
      //std::cout << "*******ERROR: MIGRAD minimization failed again: " << ierflg << std::endl;
      //std::cout << "attempting SIMPLEX method" << std::endl;
      arglist[1] = 1.0;   // tolerance
      gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflg);
      nTries++;  // 3
    }
      
    if (ierflg != 0){
      //std::cout << "*******ERROR: SIMPLEX minimization failed again: " << ierflg << std::endl;
      //std::cout << "attempting SIMPLEX method with relaxed threshold" << std::endl;
      arglist[1] = 5.0;   // tolerance
      gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflg);
      nTries++;  // 4
    }
      
    if (ierflg != 0){
      //std::cout << "*******ERROR: SIMPLEX minimization failed again: " << ierflg << std::endl;
      //std::cout << "attempting even more relaxed MIGRAD/SIMPLEX tolerance threshold" << std::endl;
      arglist[1] = 10.0;   // tolerance
      gMinuit->mnexcm("MINI", arglist, 2, ierflg); // first tries MIGRAD, switches to SIMPLEX if MIGRAD fails
      nTries++;  // 5
    }
      

    // no clue how MINOS works, apparently with the options below only two parameters are used for error analysis
    // just took it from Jason as the last option (i never found it failing)
    if (ierflg != 0){
      //std::cout << "*******ERROR: MIGRAD minimization failed again: " << ierflg << std::endl;
      //std::cout << "attempting MINOS minimization" << std::endl;
      arglist[0] = 5000; // max calls
      arglist[1] = 0;
      arglist[2] = 1;
      gMinuit->mnexcm("MINOS", arglist, 3, ierflg);
      nTries++;
    }

    //if (ierflg != 0)
    //  std::cout << "*******ERROR: MINOS minimization failed: " << ierflg << std::endl;
      

    Double_t ic [NJ];
    Double_t eic[NJ];
    
    Double_t softPxTrue=0, softPyTrue=0;
    for (int j=0; j<Nj; j++) {
      gMinuit->GetParameter(j, ic[j], eic[j]);
      softPxTrue -= jPtReco[j]*TMath::Cos(jPhiReco[j])/ic[j];
      softPyTrue -= jPtReco[j]*TMath::Sin(jPhiReco[j])/ic[j];
    }
    // std::cout << "met: "        << met        << "; metphi: "     << metphi     << std::endl
    // 	 << "softPxreco: " << softPxReco << "; softPyreco: " << softPyReco << std::endl
    // 	 << "softPxTrue: " << softPxTrue << "; softPyTrue: " << softPyTrue << std::endl;

    // my rebalanced met. uses softPtTrue (completely balanced with rebalanced jets) and PU jets
    Double_t metx_re=0., mety_re=0., met_re=0., metphi_re=0.;
    for (int j=0; j<NPUj; j++) {
      metx_re -= jPtPU[j]*TMath::Cos(jPhiPU[j]);
      mety_re -= jPtPU[j]*TMath::Sin(jPhiPU[j]);
    }
    
    if (NPUj>0){
      met_re = sqrt(metx_re*metx_re+mety_re*mety_re);
      metphi_re = mety_re==0.0 && metx_re==0.0 ? 0.0 : TMath::ATan2(mety_re,metx_re);
    }
    //std::cout << "bruno's rebalanced met: " << met_re << "; metphi: " << metphi_re << std::endl;
    
    // jason's rebalanced met. uses softPtReco and ignores PU jets
    Double_t metx_re2 = softPxTrue - softPxReco;
    Double_t mety_re2 = softPyTrue - softPyReco;
    Double_t met_re2 = sqrt(metx_re2*metx_re2+mety_re2*mety_re2);
    Double_t metphi_re2 =  mety_re2==0.0 && metx_re2==0.0 ? 0.0 : TMath::ATan2(mety_re2,metx_re2);
    //cout << "jason's rebalanced met: " << met_re2 << "; metphi: " << metphi_re2 << std::endl;
    
    int nj=0; // keep original ordering for rebalanced jets
    for (unsigned int j=0; j<all_jet_pt.size(); j++){
      if ( all_jet_pt.at(j)>100 || all_jet_puId.at(j)==1 ){
     	reb_jet_pt.push_back(all_jet_pt.at(j)/ic[nj]);
     	nj++;
      }
      else{
     	reb_jet_pt.push_back(all_jet_pt.at(j));
      }
    }

    float true_softpt  = sqrt(softPxTrue*softPxTrue+softPyTrue*softPyTrue);
    float true_softphi = softPyTrue==0.0 && softPxTrue==0.0 ? 0.0 : TMath::ATan2(softPyTrue,softPxTrue);
    float reco_softpt  = sqrt(softPxReco*softPxReco+softPyReco*softPyReco);
    float reco_softphi = softPyReco==0.0 && softPxReco==0.0 ? 0.0 : TMath::ATan2(softPyReco,softPxReco);
    
    //thisTree->yield->Fill( mt2, weight );
    thisTree->assignVar( "jet1_pt"           , myTree.jet1_pt );
    thisTree->assignVar( "jet2_pt"           , myTree.jet2_pt );
    thisTree->assignVar( "metphi"            , myTree.met_phi );
    thisTree->assignVar( "errflag"           , ierflg         );
    thisTree->assignVar( "ntries"            , nTries         );
    thisTree->assignVar( "reco_soft_pt"      , reco_softpt    );
    thisTree->assignVar( "reco_soft_phi"     , reco_softphi   );
    thisTree->assignVar( "true_soft_pt"      , true_softpt    );
    thisTree->assignVar( "true_soft_phi"     , true_softphi   );
    thisTree->assignVar( "reb_jasons_met_pt" , met_re2        );
    thisTree->assignVar( "reb_jasons_met_phi", metphi_re2     );
    thisTree->assignVar( "reb_brunos_met_pt" , met_re         );
    thisTree->assignVar( "reb_brunos_met_phi", metphi_re      );
    thisTree->assignVector( "before_jet_pt", all_jet_pt   );
    thisTree->assignVector( "reb_jet_pt"   , reb_jet_pt   );
    thisTree->assignVector( "jet_eta"      , all_jet_eta  );
    thisTree->assignVector( "jet_phi"      , all_jet_phi  );
    thisTree->assignVector( "jet_mass"     , all_jet_mass );
    thisTree->assignVector( "jet_puId"     , all_jet_puId );
    thisTree->assignVector( "jet_id"       , all_jet_id );
    thisTree->assignVector( "jet_btagCSV"  , all_jet_btagCSV );


    t.Stop();
    float realTime = t.RealTime();
    float cpuTime  = t.CpuTime();
    thisTree->assignVar( "cputime"           , cpuTime        );
    thisTree->assignVar( "realtime"          , realTime       );

    thisTree->fillTree( myTree, weight );

    delete gMinuit;
    
  } // for entries


  anaTree->finalize();


  delete tree;


  file->Close();
  delete file;
  

  std::cout << " ---> Rebalancing finished" << std::endl;

}

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
    updateIndex(j,par[j]); // the binning is gen Pt, find the corresponding bin
    logl += TMath::Log(getResponseProb(par[j], jPtIndex[j], jEtaIndex[j])); // take appropriate response function
    //logl += TMath::Log(getResponseProbGaus(par[j], jPtIndex[j], jEtaIndex[j])); // take appropriate response function (only gaussian)
    softPxTrue -= jPtReco[j]*TMath::Cos(jPhiReco[j])/par[j];
    softPyTrue -= jPtReco[j]*TMath::Sin(jPhiReco[j])/par[j];
  }
  //logl +=  TMath::Log(TMath::Gaus((softPxReco-softPxTrue)/sigmaSoft , 0.0, 1.0));
  //logl +=  TMath::Log(TMath::Gaus((softPyReco-softPyTrue)/sigmaSoft , 0.0, 1.0));
  logl +=  -TMath::Power((softPxReco-softPxTrue)/sigmaSoft,2)/2; // = log(gaus); just faster
  logl +=  -TMath::Power((softPyReco-softPyTrue)/sigmaSoft,2)/2;

  f = -logl;
}
