#include <iostream>


#include "TFile.h"
#include "TTree.h"

#define  MT2Tree8TeV_cxx
#include "MT2Tree8TeV.h"



int main( int argc, char* argv[] ) {


  std::string dir(argv[1]);
  std::string dataset(argv[2]);
  std::string fileName(argv[3]);
  

  std::cout << "-> Translating: " << dir << "/" << dataset << "/" << fileName << std::endl;
  //std::string pnfsDir = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV";


  TFile* file = TFile::Open(Form("dcap://t3se01.psi.ch:22125/%s/%s/%s", dir.c_str(), dataset.c_str(), fileName.c_str()));
  TTree* tree = (TTree*)file->Get("MassTree");



  MT2Tree8TeV mt2tree;
  mt2tree.Init(tree);


  TFile* outfile = TFile::Open(fileName.c_str(), "recreate");
  TTree* newtree = new TTree("mt2", "");
  
  int id = 1;
  newtree->Branch( "evt_id",        &id,                              "id/I" );
  newtree->Branch( "evt",           &mt2tree.misc_Event,              "mt2tree.misc_Event/I" );
  newtree->Branch( "run",           &mt2tree.misc_Run,                "mt2tree.misc_Run/I" );
  newtree->Branch( "lumi",          &mt2tree.misc_LumiSection,        "mt2tree.misc_LumiSection/I" );
  newtree->Branch( "nJet40",        &mt2tree.NJetsIDLoose40,          "mt2tree.NJetsIDLoose40/I" );
  newtree->Branch( "nBJet40",       &mt2tree.NBJets40CSVM,            "mt2tree.NBJets40CSVM/I" );
  newtree->Branch( "nBJet20",       &mt2tree.NBJets40CSVM,            "mt2tree.NBJets40CSVM/I" );
  newtree->Branch( "nMuons10",      &mt2tree.NMuons,                  "mt2tree.NMuons/I" );
  newtree->Branch( "nElectrons10",  &mt2tree.NEles,                   "mt2tree.NEles/I" );
  newtree->Branch( "nPFLep5LowMT",  &mt2tree.NTausIDLoose3Hits,       "mt2tree.NTausIDLoose3Hits/I" );
  newtree->Branch( "nPFHad10LowMT", &mt2tree.NTausIDLoose3Hits,       "mt2tree.NTausIDLoose3Hits/I" );
  newtree->Branch( "nVert",         &mt2tree.pileUp_NVertices,        "mt2tree.pileUp_NVertices/I" );
  newtree->Branch( "ht",            &mt2tree.misc_HT,                 "mt2tree.misc_HT/F" );
  newtree->Branch( "jet1_pt",       &mt2tree.misc_LeadingJPt,         "mt2tree.misc_LeadingJPt/F" );
  newtree->Branch( "jet2_pt",       &mt2tree.misc_SecondJPt,          "mt2tree.misc_SecondJPt/F" );
  newtree->Branch( "deltaPhiMin",   &mt2tree.misc_MinMetJetDPhi4Pt40, "mt2tree.misc_MinMetJetDPhi4Pt40/F" );
  newtree->Branch( "diffMetMht",    &mt2tree.misc_Vectorsumpt,        "mt2tree.misc_Vectorsumpt/F" );
  newtree->Branch( "mt2",           &mt2tree.misc_MT2,                "mt2tree.misc_MT2/F" );
  newtree->Branch( "met_pt",        &mt2tree.misc_MET,                "mt2tree.misc_MET/F" );


  int nentries = tree->GetEntries();

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry%10000 == 0) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;
    mt2tree.GetEntry(iEntry);
    newtree->Fill();

  }

  outfile->cd();
  newtree->Write();
  outfile->Close();

  return 0;

}
