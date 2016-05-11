#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TLorentzVector.h"


#define mt2_cxx
#include "interface/mt2.h"


//double PtBinEdges[23] = {0, 20, 30, 50, 80, 120, 170, 230, 300, 380, 470, 570, 680, 800, 1000, 1300, 1700, 2200, 2800, 3500, 4300, 5200, 6500};
//double EtaBinEdges[12] = {0, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.3, 2.8, 3.2, 4.1, 5.0};

std::vector< std::vector< TH1D* > > allocateTemplates( std::vector<double> ptBins, std::vector<double> etaBins );
void fillTemplates( float pt, float ptGen, float eta, float weight, std::vector< std::vector< TH1D*> > templates, std::vector<double> ptBins, std::vector<double> etaBins );
int getBin( float var, std::vector<double> v );



int main( int argc, char* argv[] ) {

  if( argc==1 ) {
    std::cout << "-> Usage: ./makeJetResponseTemplates [dataset]" << std::endl;
    exit(11);
  }


  std::vector<double> ptBins;
  ptBins.push_back( 0. );
  ptBins.push_back( 20. );
  ptBins.push_back( 30. );
  ptBins.push_back( 50. );
  ptBins.push_back( 80. );
  ptBins.push_back( 120. );
  ptBins.push_back( 170. );
  ptBins.push_back( 230. );
  ptBins.push_back( 300. );
  ptBins.push_back( 380. );
  ptBins.push_back( 470. );
  ptBins.push_back( 570. );
  ptBins.push_back( 680. );
  ptBins.push_back( 800. );
  ptBins.push_back( 1000. );
  ptBins.push_back( 1300. );
  ptBins.push_back( 1700. );
  ptBins.push_back( 2200. );
  ptBins.push_back( 2800. );
  ptBins.push_back( 3500. );
  ptBins.push_back( 4300. );
  ptBins.push_back( 5200. );
  ptBins.push_back( 6500. );

//double EtaBinEdges[12] = {0, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.3, 2.8, 3.2, 4.1, 5.0};
  std::vector<double> etaBins;
  etaBins.push_back( 0. );
  etaBins.push_back( 0.3 );
  etaBins.push_back( 0.5 );
  etaBins.push_back( 0.8 );
  etaBins.push_back( 1.1 );
  etaBins.push_back( 1.4 );
  etaBins.push_back( 1.7 );
  etaBins.push_back( 2.3 );
  etaBins.push_back( 2.8 );
  etaBins.push_back( 3.2 );
  etaBins.push_back( 4.1 );
  etaBins.push_back( 5.0 );


  std::vector< std::vector< TH1D* > > templates = allocateTemplates( ptBins, etaBins );

  std::string path = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/casal/MT2production/74X/Spring15/PostProcessed/QCDpt10_gen/";
  std::string dataset(argv[1]);
  std::string fileName(Form("%s/%s_post.root", path.c_str(), dataset.c_str()));
  TFile* file = TFile::Open(fileName.c_str());
  TTree* tree = (TTree*)file->Get("mt2");

  TString dataset_tstr(dataset);
  TString dataset_tstr2 = dataset_tstr.ReplaceAll("QCD_HT", "");
  TString dataset_tstr3 = dataset_tstr.ReplaceAll("to", " ");


  std::vector<std::string> internal;
  std::istringstream iss(dataset_tstr3.Data()); 
  std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                      std::istream_iterator<std::string>{}};
  
  float htMin = atof(tokens[0].c_str());
  float htMax = atof(tokens[1].c_str());
  

  MT2Tree myTree;
  myTree.Init(tree);

  int nentries = tree->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);

    if( myTree.njet<2 ) continue;
    if( myTree.ngenJet<2 ) continue;

    TLorentzVector jet0, jet1;
    jet0.SetPtEtaPhiM( myTree.jet_pt[0], myTree.jet_eta[0], myTree.jet_phi[0], myTree.jet_mass[0] );
    jet1.SetPtEtaPhiM( myTree.jet_pt[1], myTree.jet_eta[1], myTree.jet_phi[1], myTree.jet_mass[1] );
    if( fabs( jet0.DeltaPhi(jet1) ) < 3.14159 - 0.5 ) continue;

    bool dijetEvent = myTree.ngenJet==2 || ( myTree.ngenJet>2 && myTree.genJet_pt[2] < 0.3*0.5*(myTree.genJet_pt[0] + myTree.genJet_pt[1]) );
    if( !dijetEvent ) continue;

    if( myTree.jet_mcPt[0]>htMin/2. && myTree.jet_mcPt[0]<htMax/2. )
      fillTemplates( myTree.jet_pt[0], myTree.jet_mcPt[0], myTree.jet_eta[0], myTree.evt_scale1fb, templates, ptBins, etaBins ); 
    if( myTree.jet_mcPt[1]>htMin/2. && myTree.jet_mcPt[1]<htMax/2. )
      fillTemplates( myTree.jet_pt[1], myTree.jet_mcPt[1], myTree.jet_eta[1], myTree.evt_scale1fb, templates, ptBins, etaBins ); 

  }

  TFile* outfile = TFile::Open(Form("templates_%s.root", dataset.c_str()), "recreate");
  outfile->cd();
 
  for( unsigned i=0; i<templates.size(); ++i ) {
    for( unsigned j=0; j<templates[i].size(); ++j ) {
      templates[i][j]->Write();
    }
  }

  outfile->Close();

  std::cout << "-> Find your templates in: " << outfile->GetName() << std::endl;

  return 0;

}





std::vector< std::vector< TH1D* > > allocateTemplates( std::vector<double> ptBins, std::vector<double> etaBins ) {

  std::vector< std::vector< TH1D* > > templates;


  for( unsigned ipt=0; ipt<ptBins.size(); ++ipt ) {

    std::vector< TH1D* > thisEtaVec;

    for( unsigned ieta=0; ieta<etaBins.size(); ++ieta ) {

      std::string thisName(Form("h_tot_JetAll_ResponsePt_Pt%d_Eta%d", ipt, ieta));
      TH1D* thisTempl = new TH1D( thisName.c_str(), "", 150, 0., 3.); 
      thisTempl->Sumw2();
      thisEtaVec.push_back(thisTempl);

    } // for eta

    templates.push_back( thisEtaVec );

  } // for pt

  return templates;

}



void fillTemplates( float pt, float ptGen, float eta, float weight, std::vector< std::vector< TH1D*> > templates, std::vector<double> ptBins, std::vector<double> etaBins ) {

  int iPtBin  = getBin( ptGen    , ptBins  );
  int iEtaBin = getBin( fabs(eta), etaBins );


  if( iPtBin>=0 && iEtaBin>=0 )
    templates[iPtBin][iEtaBin]->Fill( pt/ptGen, weight );
  else {
    std::cout << "pt: " << pt << " corresponds to iPtBin: " << iPtBin << std::endl;
    std::cout << "eta: " << eta << " corresponds to iEtaBin: " << iEtaBin << std::endl;
  }

}


int getBin( float var, std::vector<double> v ) {

  if( var>v[v.size()-1] ) return v.size()-1;

  int returnBin = -1;
  for( unsigned i=0; i<v.size()-1; ++i ) {
    if( var>=v[i] && var<v[i+1] ) {
      returnBin = i;
      break;
    }
  }

  return returnBin;

}
