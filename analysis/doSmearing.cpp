#include "interface/MT2Config.h"
#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"
#include "interface/Hemisphere.h"
#include "interface/Davismt2.h"
#include "TMinuit.h"
#include "TSystem.h"
#include "TRandom.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "interface/rooDoubleCB.h"
#include "TStopwatch.h"

#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1F.h"
#include "TF1.h"

int ijob=0, Njobs=1;

Double_t sigmaSoft = 20.; 

int nSmearings = 10;

bool smearPUjets = false; // smear also PU jets? (Jason didn't)

const Int_t NJ=50;
Int_t Nj, NPUj;
Double_t jPtReb  [NJ];
Double_t jPhiReb [NJ];
Double_t jEtaReb [NJ];
Double_t jMassReb[NJ];
Int_t jPtIndex[NJ], jEtaIndex[NJ];


// open file with response functions
double PtBinEdges[23] = {0, 20, 30, 50, 80, 120, 170, 230, 300, 380, 470, 570, 680, 800, 1000, 1300, 1700, 2200, 2800, 3500, 4300, 5200, 6500};
double EtaBinEdges[12] = {0, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.3, 2.8, 3.2, 4.1, 5.0};
//RooWorkspace *w[22][12];
TF1 *ftemplate[22][12];
TFile *fTemplates = new TFile("/shome/casal/templatesRandS/responseTemplates.root");

void smear( MT2Analysis<MT2EstimateTree>* inputAnaTree,  MT2Analysis<MT2EstimateTree>* outputAnaTree );
void getTemplates();
void updateIndex(int j, double x);
Double_t getRandomFromTemplate(int j);
std::vector<float> smearJets(std::vector<float> *before_jet_pt, std::vector<float> *jet_puID);
void smearSoftPt(float &pt, float &phi);
void recalculateVars(MT2EstimateTree *tree, std::vector<float> smear_jet_pt, float softpt, float softphi, float &jet1_pt, float &jet2_pt, float &ht, float &met, float &mt2, float &dPhiMin, float &diffMetMht, int &njets, int &nbjets);
Float_t calcMT2(std::vector<float> smear_jet_pt, float met, float metphi);
void getHemispheres(std::vector<float> smear_jet_pt, TLorentzVector *v1, TLorentzVector *v2);
void fillTree(MT2Analysis<MT2EstimateTree>* outputAnaTree, MT2EstimateTree* inputTree, std::vector<float> smear_jet_pt, float softpt, float softphi, float jet1_pt, float jet2_pt, float ht, float met, float mt2, float dPhiMin, float diffMetMht, int njets, int nbjets, int iSmear);
void addVars(MT2Analysis<MT2EstimateTree>* anaTree, bool smearTree);


int main( int argc, char* argv[] ) {



  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|              smearing events              |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  if( argc<2 ) {
    std::cout << "USAGE: ./doSmearing configFileName [data/MC/all] [sampleID] [job_i] [Njobs]" << std::endl;
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

  std::string inputdir  = cfg.getEventYieldDir() + "/rebalancedTrees"; 
  std::string outputdir = cfg.getEventYieldDir() + "/smearedTrees"; 
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

    TString imcFile = inputdir + "/mc";
    if (sampleID!=-1)
      imcFile += TString::Format("_id%d_job%dof%d", sampleID, ijob, Njobs);
    imcFile += ".root";

    MT2Analysis<MT2EstimateTree>* rebalanceTree = MT2Analysis<MT2EstimateTree>::readFromFile( imcFile.Data(), "qcdRebalancedTree" );
    addVars(rebalanceTree, false); // add extra vars for rebalance tree
    

    MT2Analysis<MT2EstimateTree>* smearTree = new MT2Analysis<MT2EstimateTree>( "qcdSmearedTree", "13TeV_inclusive" );
    addVars(smearTree, true); // add extra vars for smear tree

    
    smear( rebalanceTree, smearTree );
    

    TString omcFile = outputdir + "/mc";
    if (sampleID!=-1)
      omcFile += TString::Format("_id%d_job%dof%d", sampleID, ijob, Njobs);
    omcFile += ".root";
    smearTree->writeToFile( omcFile.Data(), "RECREATE" );

  }


  if( !(cfg.dummyAnalysis()) && cfg.dataSamples()!="" && !onlyMC  ) {  // run on data

    //to be implemented

  }



  return 0;

}


void smear( MT2Analysis<MT2EstimateTree>* inputAnaTree, MT2Analysis<MT2EstimateTree>* outputAnaTree ) {


  // loop over regions
  std::set<MT2Region> regions = inputAnaTree->getRegions();
  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateTree* estimateTree  = inputAnaTree->get( *iR );
    TTree* tree = estimateTree->tree;

    int nentries = tree->GetEntries();

    std::cout << "-> Loaded tree: it has " << nentries << " entries." << std::endl;


    // std::vector<float> *before_jet_pt  = new std::vector<float>();
    // std::vector<float> *reb_jet_pt     = new std::vector<float>();
    // std::vector<float> *jet_eta  = new std::vector<float>();
    // std::vector<float> *jet_phi  = new std::vector<float>();
    // std::vector<float> *jet_puId = new std::vector<float>();
    // tree->SetBranchAddress( "before_jet_pt", &(estimateTree->extraVectors["before_jet_pt"]));//before_jet_pt );
    // tree->SetBranchAddress( "reb_jet_pt"   , &(estimateTree->extraVectors["reb_jet_pt"   ]));//reb_jet_pt    );
    // tree->SetBranchAddress( "jet_eta"      , &(estimateTree->extraVectors["jet_eta"      ]));//jet_eta    );
    // tree->SetBranchAddress( "jet_phi"      , &(estimateTree->extraVectors["jet_phi"      ]));//jet_phi    );
    // tree->SetBranchAddress( "jet_puId"     , &(estimateTree->extraVectors["jet_puId"     ]));//jet_puId   );

    estimateTree->initTree4read();
  
    nentries = 100;

    for( int iEntry=0; iEntry<nentries; ++iEntry ) {
      
      
      if( iEntry % 25000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
      
      tree->GetEntry(iEntry);


      std::vector<float> *before_jet_pt  = estimateTree->extraVectors["before_jet_pt"];
      std::vector<float> *reb_jet_pt     = estimateTree->extraVectors["reb_jet_pt"   ];
      std::vector<float> *jet_eta        = estimateTree->extraVectors["jet_eta"      ];
      std::vector<float> *jet_phi        = estimateTree->extraVectors["jet_phi"      ];
      std::vector<float> *jet_puId       = estimateTree->extraVectors["jet_puId"     ];
      std::vector<float> *jet_mass       = estimateTree->extraVectors["jet_mass"     ];

      // update rebalanced jets and response templates for the current event
      Nj=reb_jet_pt->size();
      for (int j=0; j<Nj; j++) {
	jPtReb  [j] = reb_jet_pt ->at(j);
	jEtaReb [j] = jet_eta ->at(j);
	jPhiReb [j] = jet_phi ->at(j);
	jMassReb[j] = jet_mass->at(j);
	updateIndex(j, 1.0); // 1.0 means this jet acts as genJet
      }
      
      for ( int iSmear=0; iSmear<nSmearings; iSmear++) {
	
	std::vector<float> smear_jet_pt = smearJets(before_jet_pt, jet_puId);
	
	// jason's choice is recoSoftPt w/o smearing
	float softPt  = *(estimateTree->extraVars["true_soft_pt" ]);
	float softPhi = *(estimateTree->extraVars["true_soft_phi"]);
	smearSoftPt(softPt, softPhi);
	
	
	float j1pt, j2pt, ht, met, mt2, dPhiMin, diffMetMht;
	int njets, nbjets;
	recalculateVars(estimateTree, smear_jet_pt, softPt, softPhi, j1pt, j2pt,ht, met, mt2, dPhiMin, diffMetMht, njets, nbjets);
	
	fillTree(outputAnaTree,estimateTree, smear_jet_pt, softPt, softPhi, j1pt, j2pt,ht, met, mt2, dPhiMin, diffMetMht, njets, nbjets, iSmear+1); 
	
      }
      
    } // for entries
   
  } 
  
  outputAnaTree->finalize();

  std::cout << " ---> Smearing finished" << std::endl;
  
}
  
  
void getTemplates(){
  for (int pt=0; pt<22; pt++){
    for (int eta=0; eta<12; eta++){
      TString temp = TString::Format("wsp_h_tot_JetAll_ResponsePt_Pt%d_Eta%d", pt, eta);
      if ( fTemplates->GetListOfKeys()->Contains(temp) ){
  	RooWorkspace* w = (RooWorkspace*)fTemplates->Get(temp);
	ftemplate[pt][eta] = (TF1*)w->pdf("cb")->asTF( RooArgList(*(w->var("x"))), RooArgList(*(w->var("mean")),*(w->var("sigma")),*(w->var("a")),*(w->var("n")),*(w->var("aDx")),*(w->var("nDx")))); // "arbitrary" normalization
      }
    }
  }
}
 
void updateIndex(int j, double x){
  int p=0, e=0;
  double pt = jPtReb[j]/x;
  //for (p=0; p<22; p++){
  for (p=0; p<16; p++){  // only up to bin Pt16, empty for higher pt
    if (PtBinEdges[p+1]>pt)
      break;
  }	  

  if(pt < 300.0){ // all eta bins
    for (e=0; e<11; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReb[j]))
	break;
    }
  }
  else if(pt < 570.0){ // up to bin Eta9
    for (e=0; e<9; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReb[j]))
	break;
    }
  }
  else if(pt < 1000.0){ // up to bin Eta8
    for (e=0; e<8; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReb[j]))
	break;
    }
  }
  else if(pt < 1300.0){ // up to bin Eta7
    for (e=0; e<7; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReb[j]))
  	  break;
    }
  }
  else if(pt < 1700.0){ // up to bin Eta6
    for (e=0; e<6; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReb[j]))
	break;
    }
  }
  else{ // up to bin Eta2
    for (e=0; e<2; e++){
      if (EtaBinEdges[e+1]>fabs(jEtaReb[j]))
	break;
    }
  }
  jPtIndex [j]=p;
  jEtaIndex[j]=e;
}

Double_t getRandomFromTemplate(int j){
  return ftemplate[jPtIndex[j]][jEtaIndex[j]]->GetRandom();
}

std::vector<float> smearJets(std::vector<float> *before_jet_pt, std::vector<float> *jet_puId){
  std::vector<float> smearedJets_pt;
  // if smearPUjets, use response template assuming their reco pt is the true pt
  for (int j=0; j<Nj; j++){
    if ( smearPUjets || before_jet_pt->at(j)>100 || jet_puId->at(j)==1 )
      smearedJets_pt.push_back(jPtReb[j]*getRandomFromTemplate(j));
    else
      smearedJets_pt.push_back(jPtReb[j]);
  }
  
  return smearedJets_pt;
}

void smearSoftPt(float &pt, float &phi){
  float px = gRandom->Gaus(pt*TMath::Cos(phi), sigmaSoft);
  float py = gRandom->Gaus(pt*TMath::Sin(phi), sigmaSoft);
  pt = sqrt(px*px+py*py);
  phi = py==0.0 && px==0.0 ? 0.0 : TMath::ATan2(py,px);
}


void recalculateVars(MT2EstimateTree *tree, std::vector<float> smear_jet_pt, float softpt, float softphi, float &jet1_pt, float &jet2_pt, float &ht, float &met, float &mt2, float &dPhiMin, float &diffMetMht, int &njets, int &nbjets){

  // recalculate met ( = -ptsoft -PUjetspt (=reb_brunos_met) - smearjetspt
  float rebMetPt  = *(tree->extraVars["reb_brunos_met_pt" ]);
  float rebMetPhi = *(tree->extraVars["reb_brunos_met_phi"]);
  float metx = -softpt*TMath::Cos(softphi) + rebMetPt*TMath::Cos(rebMetPhi);
  float mety = -softpt*TMath::Sin(softphi) + rebMetPt*TMath::Sin(rebMetPhi);
  float diffmetmhtx=metx, diffmetmhty=mety;
  for (unsigned int j = 0; j<smear_jet_pt.size(); j++){
    metx -= smear_jet_pt.at(j)*TMath::Cos(jPhiReb[j]);
    mety -= smear_jet_pt.at(j)*TMath::Sin(jPhiReb[j]);
    if ( !(smear_jet_pt.at(j)>30 && fabs(jEtaReb[j])<2.5) ){
      diffmetmhtx -= smear_jet_pt.at(j)*TMath::Cos(jPhiReb[j]);
      diffmetmhty -= smear_jet_pt.at(j)*TMath::Sin(jPhiReb[j]);
    }
  }
  met    = sqrt(metx*metx+mety*mety);
  float metphi = mety==0.0 && metx==0.0 ? 0.0 : TMath::ATan2(mety,metx);
  diffMetMht = sqrt(diffmetmhtx*diffmetmhtx+diffmetmhty*diffmetmhty);


  ht = 0.; jet1_pt = 0.;  jet2_pt = 0.; dPhiMin = 3.5;
  njets = 0; nbjets = 0;
  for (unsigned int j = 0; j<smear_jet_pt.size(); j++){
    if (smear_jet_pt.at(j) > jet1_pt){
      jet2_pt = jet1_pt;
      jet1_pt = smear_jet_pt.at(j);
    }
    else if (smear_jet_pt.at(j) > jet2_pt)
      jet2_pt = smear_jet_pt.at(j);
    if (smear_jet_pt.at(j)>30 && fabs(jEtaReb[j])<2.5){
      ht += smear_jet_pt.at(j);
      njets++;
    }
    float dphi = fabs(TVector2::Phi_mpi_pi(metphi-jPhiReb[j]));
    if ( smear_jet_pt.at(j)>30 && fabs(jEtaReb[j])<4.7 && dphi < dPhiMin)
      dPhiMin = dphi;
  }

  //recalculate mt2
  mt2 = njets<2 ? -9.9 : calcMT2(smear_jet_pt, met, metphi);

}

Float_t calcMT2(std::vector<float> smear_jet_pt, float met, float metphi){

  TLorentzVector *visible1 = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *visible2 = new TLorentzVector(0.,0.,0.,0.);

  getHemispheres(smear_jet_pt, visible1, visible2);

  double pa[3];
  double pb[3];
  double pmiss[3];

  pmiss[0] = 0;
  pmiss[1] = static_cast<double> (met*TMath::Cos(metphi));
  pmiss[2] = static_cast<double> (met*TMath::Sin(metphi));

  pa[0] = 0;
  pa[1] = static_cast<double> (visible1->Px());
  pa[2] = static_cast<double> (visible1->Py());
  
  pb[0] = 0;
  pb[1] = static_cast<double> (visible2->Px());
  pb[2] = static_cast<double> (visible2->Py());

  Davismt2 *mt2 = new Davismt2();
  mt2->set_momenta(pa, pb, pmiss);
  mt2->set_mn(0);
  Float_t MT2=mt2->get_mt2();
  delete mt2;
  delete visible1;
  delete visible2;
  return MT2;

}

void getHemispheres(std::vector<float> smear_jet_pt, TLorentzVector *v1, TLorentzVector *v2){
  
  std::vector<float> px, py, pz, E;
  int NJets = smear_jet_pt.size();
  for(int j=0; j<NJets; ++j){
    TLorentzVector jet;
    if (smear_jet_pt.at(j) > 30 && fabs(jEtaReb[j])<2.5)
      jet.SetPtEtaPhiM(smear_jet_pt.at(j), jEtaReb[j], jPhiReb[j], jMassReb[j]);
    else
      continue;
    px.push_back(jet.Px());
    py.push_back(jet.Py());
    pz.push_back(jet.Pz());
    E .push_back(jet.E ());
  }

  Hemisphere* hemisp = new Hemisphere(px, py, pz, E, 2, 3);
  std::vector<int> grouping = hemisp->getGrouping();

  v1->SetPxPyPzE(0.,0.,0.,0.);
  v2->SetPxPyPzE(0.,0.,0.,0.);
  for(unsigned int i=0; i<px.size(); ++i){
	if(grouping[i]==1){
		v1->SetPx(v1->Px() + px[i]);
		v1->SetPy(v1->Py() + py[i]);
		v1->SetPz(v1->Pz() + pz[i]);
		v1->SetE (v1->E () + E [i]);	
	}else if(grouping[i] == 2){
		v2->SetPx(v2->Px() + px[i]);
		v2->SetPy(v2->Py() + py[i]);
		v2->SetPz(v2->Pz() + pz[i]);
		v2->SetE (v2->E () + E [i]);
	}
  }
  delete hemisp;

}

void fillTree(MT2Analysis<MT2EstimateTree>* outputAnaTree, MT2EstimateTree* inputTree, std::vector<float> smear_jet_pt, float softpt, float softphi, float jet1_pt, float jet2_pt, float ht, float met, float mt2, float dPhiMin, float diffMetMht, int njets, int nbjets, int iSmear){

  MT2EstimateTree* thisTree = outputAnaTree->get( ht, njets, nbjets, -1, mt2 );
  if( thisTree==0 ) return;

  thisTree->run         = inputTree->run        ;
  thisTree->lumi        = inputTree->lumi       ;
  thisTree->evt         = inputTree->evt        ;
  thisTree->weight      = inputTree->weight     ;
  thisTree->puWeight    = inputTree->puWeight   ;
  thisTree->id          = inputTree->id         ;
  thisTree->nVert       = inputTree->nVert      ;  
  thisTree->nElectrons  = inputTree->nElectrons ;
  thisTree->nMuons      = inputTree->nMuons     ;
  thisTree->nPFLep      = inputTree->nPFLep     ;
  thisTree->nPFHad      = inputTree->nPFHad     ;
  thisTree->nJetHF      = inputTree->nJetHF     ;

  thisTree->mt2         = mt2;
  thisTree->ht          = ht;
  thisTree->met         = met;
  thisTree->deltaPhiMin = dPhiMin;
  thisTree->diffMetMht  = diffMetMht;
  thisTree->nJets       = njets;
  thisTree->nBJets      = nbjets;

  
  thisTree->assignVar( "jet1_pt"           , jet1_pt );
  thisTree->assignVar( "jet2_pt"           , jet2_pt );
  //thisTree->assignVar( "reco_soft_pt"      ,  );
  //thisTree->assignVar( "reco_soft_phi"     ,  );
  thisTree->assignVar( "true_soft_pt"      ,  *(inputTree->extraVars["true_soft_pt" ]));
  thisTree->assignVar( "true_soft_phi"     ,  *(inputTree->extraVars["true_soft_phi"]));
  //thisTree->assignVar( "reb_jasons_met_pt" ,  );
  //thisTree->assignVar( "reb_jasons_met_phi",  );
  thisTree->assignVar( "reb_brunos_met_pt" , *(inputTree->extraVars["reb_brunos_met_pt" ]) ); // contribution from PU jets
  thisTree->assignVar( "reb_brunos_met_phi", *(inputTree->extraVars["reb_brunos_met_phi"]) );
  thisTree->assignVector( "jet_puId"      , *(inputTree->extraVectors["jet_puId"     ]));
  thisTree->assignVector( "before_jet_pt" , *(inputTree->extraVectors["before_jet_pt"]));
  thisTree->assignVector( "reb_jet_pt"    , *(inputTree->extraVectors["reb_jet_pt"   ]));
  thisTree->assignVector( "jet_eta"       , *(inputTree->extraVectors["jet_eta"      ]));
  thisTree->assignVector( "jet_phi"       , *(inputTree->extraVectors["jet_phi"      ]));
  thisTree->assignVar( "before_jet1_pt"    , *(inputTree->extraVars["jet1_pt"]) );
  thisTree->assignVar( "before_jet2_pt"    , *(inputTree->extraVars["jet2_pt"]) );
  thisTree->assignVar( "before_met"        , inputTree->met        );
  thisTree->assignVar( "before_mt2"        , inputTree->mt2        );
  thisTree->assignVar( "before_ht"         , inputTree->ht         );
  thisTree->assignVar( "before_nJets"      , inputTree->nJets      );
  thisTree->assignVar( "before_nBJets"     , inputTree->nBJets     );
  thisTree->assignVar( "before_deltaPhiMin", inputTree->deltaPhiMin);
  thisTree->assignVar( "before_diffMetMht" , inputTree->diffMetMht );
  thisTree->assignVar( "smear_soft_pt"  , softpt  );
  thisTree->assignVar( "smear_soft_phi" , softphi );
  thisTree->assignVar( "iSmear"         , iSmear  );
  thisTree->assignVector( "after_jet_pt", smear_jet_pt );
  
  thisTree->tree->Fill();
  

}

void addVars(MT2Analysis<MT2EstimateTree>* anaTree, bool smearTree){

    MT2EstimateTree::addVar( anaTree, "jet1_pt"      );
    MT2EstimateTree::addVar( anaTree, "jet2_pt"      );
    //MT2EstimateTree::addVar( anaTree, "reco_soft_pt" );
    //MT2EstimateTree::addVar( anaTree, "reco_soft_phi" );
    MT2EstimateTree::addVar( anaTree, "true_soft_pt" );
    MT2EstimateTree::addVar( anaTree, "true_soft_phi");
    //MT2EstimateTree::addVar( anaTree, "reb_jasons_met_pt"  );
    //MT2EstimateTree::addVar( anaTree, "reb_jasons_met_phi" );
    MT2EstimateTree::addVar( anaTree, "reb_brunos_met_pt"  ); // contribution from PU jets
    MT2EstimateTree::addVar( anaTree, "reb_brunos_met_phi" );
    MT2EstimateTree::addVector( anaTree, "before_jet_pt" );
    MT2EstimateTree::addVector( anaTree, "reb_jet_pt"    );
    MT2EstimateTree::addVector( anaTree, "jet_eta"       );
    MT2EstimateTree::addVector( anaTree, "jet_phi"       );
    MT2EstimateTree::addVector( anaTree, "jet_puId"      );
    MT2EstimateTree::addVector( anaTree, "jet_mass"      );
    MT2EstimateTree::addVector( anaTree, "jet_id"        );
    MT2EstimateTree::addVector( anaTree, "jet_btagCSV"   );

    if(smearTree){
      MT2EstimateTree::addVar( anaTree, "before_jet1_pt"     );
      MT2EstimateTree::addVar( anaTree, "before_jet2_pt"     );
      MT2EstimateTree::addVar( anaTree, "before_met"         );
      MT2EstimateTree::addVar( anaTree, "before_mt2"         );
      MT2EstimateTree::addVar( anaTree, "before_ht"          );
      MT2EstimateTree::addVar( anaTree, "before_nJets"       );
      MT2EstimateTree::addVar( anaTree, "before_nBJets"      );
      MT2EstimateTree::addVar( anaTree, "before_deltaPhiMin" );
      MT2EstimateTree::addVar( anaTree, "before_diffMetMht"  );
      MT2EstimateTree::addVar( anaTree, "smear_soft_pt"      );
      MT2EstimateTree::addVar( anaTree, "smear_soft_phi"     );
      MT2EstimateTree::addVar( anaTree, "iSmear"             );
      MT2EstimateTree::addVector( anaTree, "after_jet_pt"  );
    }
}
