#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"


#define mt2_cxx
#include "../interface/mt2.h"


#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"


float lumi = 4.; //fb-1
bool dummyData = true;



MT2Analysis<MT2EstimateTree> computeYield( const MT2Sample& sample, const std::string& regionsSet, const std::string& prompt_fake="all");
TProfile* fitSingleVariable( const std::string& outputdir, TTree* tree, const std::string& varName, const std::string& axisName, float yMax );
void compareTemplates( const std::string& outputdir, TTree* tree, TF1* corr );



int main( int argc, char* argv[] ) {


  std::string regionsSet = "13TeV_inclusive";

  std::string samplesFileName = "PHYS14_v2_Zinv_noSietaieta";
  //std::string samplesFileName = "CSA14_Zinv";

  std::string samplesFile = "../samples/samples_" + samplesFileName + ".dat";
  
  std::vector<MT2Sample> samples = MT2Sample::loadSamples(samplesFile, 100, 199); // QCD only
  if( samples.size()==0 ) {
    std::cout << "There must be an error: didn't find any files in " << samplesFile << "!" << std::endl;
    exit(1209);
  }


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  
  MT2Analysis<MT2EstimateTree>* templatesFake = new MT2Analysis<MT2EstimateTree>( "templatesFake", regionsSet );


  for( unsigned i=0; i<samples.size(); ++i ) {
    (*templatesFake) += (computeYield( samples[i], regionsSet, "fake" ));
  }


  std::string outputdir = "SigmaFits_" + samplesFileName + "_" + regionsSet;
  system(Form("mkdir -p %s", outputdir.c_str()));

  std::string outfilename = outputdir + "/fits_" + samplesFileName + "_" + regionsSet + ".root";
  templatesFake->writeToFile(outfilename);
  
  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this
  
  TTree* tree = templatesFake->get(*(templatesFake->getRegions().begin()))->tree;

  MT2DrawTools::setStyle();


  TProfile* hp_iso_vs_sigma     = fitSingleVariable( outputdir, tree, "iso", "Photon Relative Charged Isolation", 0.1 );
  TProfile* hp_ptGamma_vs_sigma = fitSingleVariable( outputdir, tree, "ptGamma", "Photon p_{T} [GeV]", 500. );
  TProfile* hp_isoAbs_vs_sigma  = fitSingleVariable( outputdir, tree, "isoAbs", "Photon Charged Isolation [GeV]", 30. );

  // add to file
  TFile* file = TFile::Open(outfilename.c_str(), "update");
  hp_iso_vs_sigma->Write();
  hp_isoAbs_vs_sigma->Write();
  hp_ptGamma_vs_sigma->Write();
  file->Close();


  compareTemplates( outputdir, tree, hp_iso_vs_sigma->GetFunction("f1") );


  return 0;

}




TProfile* fitSingleVariable( const std::string& outputdir, TTree* tree, const std::string& varName, const std::string& axisName, float yMax ) {

  std::string profileName = varName + "vs_sigma";

  TProfile* hp_var_vs_sigma = new TProfile(profileName.c_str(), "", 19, 0., 0.019);
  hp_var_vs_sigma->Sumw2();
  tree->Project( profileName.c_str(), Form("%s:sietaieta", varName.c_str()), "weight", "profile" );

  TF1* f1 = new TF1("f1", "[0] + [1]*x", 0.0085, 0.018);
  //TF1* f1 = new TF1("f1", "[0] + [1]*x + [2]*x*x", 0.0085, 0.018);
  f1->SetLineColor(kRed);
  f1->SetParameter( 0, 0.024 );
  f1->SetParameter( 1, 0.1 );
  hp_var_vs_sigma->Fit( f1, "R+" );

  
  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 0.008, 0.019, 10, 0., yMax );
  h2_axes->SetXTitle( "Photon #sigma_{i#eta i#eta}" );
  h2_axes->SetYTitle( axisName.c_str() );
  h2_axes->Draw();

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  TPaveText* fakeLabel = new TPaveText( 0.7, 0.18, 0.9, 0.28, "brNDC");
  fakeLabel->SetTextAlign(31);
  fakeLabel->SetFillColor(0);
  fakeLabel->SetTextSize(0.035);
  fakeLabel->AddText("Fake Photons");
  fakeLabel->AddText("Barrel Only");
  fakeLabel->Draw("same");

  TLine* lineCut = new TLine( 0.01, 0., 0.01, yMax );
  lineCut->SetLineStyle(2);
  lineCut->SetLineColor(kGray);
  lineCut->Draw("same");

  hp_var_vs_sigma->SetMarkerStyle(20);
  hp_var_vs_sigma->SetMarkerSize(1.3);
  hp_var_vs_sigma->SetLineColor(kBlack);
  hp_var_vs_sigma->SetMarkerColor(kBlack);
  hp_var_vs_sigma->Draw("P same");

  gPad->RedrawAxis();

  c1->SaveAs(Form("%s/%s.eps", outputdir.c_str(), varName.c_str()));
  c1->SaveAs(Form("%s/%s.pdf", outputdir.c_str(), varName.c_str()));
  c1->SaveAs(Form("%s/%s.png", outputdir.c_str(), varName.c_str()));
  
  delete c1;
  delete h2_axes;

  return hp_var_vs_sigma;

}




void compareTemplates( const std::string& outputdir, TTree* tree, TF1* corr ) {

  int nbins = 12;
  Double_t bins[nbins];
  bins[0] = 0.;
  bins[1] = 0.005;
  bins[2] = 0.01;
  bins[3] = 0.02;
  bins[4] = 0.03;
  bins[5] = 0.04;
  bins[6] = 0.05;
  bins[7] = 0.06;
  bins[8] = 0.07;
  bins[9] = 0.08;
  bins[10] = 0.09;
  bins[11] = 0.1;

  TH1D* templ_sr = new TH1D("templ_sr", "", nbins-1, bins);
  templ_sr->Sumw2();
  TH1D* templ_side = new TH1D("templ_side", "", nbins-1, bins);
  templ_side->Sumw2();
  TH1D* templ_sideCorr = new TH1D("templ_sideCorr", "", nbins-1, bins);
  templ_sideCorr->Sumw2();

  float weight;
  tree->SetBranchAddress("weight", &weight ); 
  float iso;
  tree->SetBranchAddress("iso", &iso ); 
  float sigma;
  tree->SetBranchAddress("sietaieta", &sigma ); 

  int nentries = tree->GetEntries();

  for( unsigned iEntry = 0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( sigma<0.01 ) {
      templ_sr->Fill( iso, weight );
    } else if( sigma < 0.014 ) {
      templ_side->Fill( iso, weight );
      //float f_cut  = 0.045;
      float f_cut  = corr->Eval(0.0092);
      float f_this  = corr->Eval(sigma);
      float k = f_cut/f_this;
      templ_sideCorr->Fill( iso*k, weight );
    }

  }

 
  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 0.1, 10, 0., 0.25);
  h2_axes->SetXTitle("Photon Relative Charged Isolation");      
  h2_axes->SetYTitle("Normalized to Unity");
  h2_axes->Draw();

  templ_sr->SetLineColor(46);
  templ_sr->SetLineWidth(2);
  templ_sr->DrawNormalized("same");

  templ_side->SetLineColor(kBlack);
  templ_side->SetLineWidth(2);
  templ_side->DrawNormalized("same");

  //templ_sideCorr->SetLineColor(38);
  //templ_sideCorr->SetLineWidth(2);
  //templ_sideCorr->DrawNormalized("histo same");

  TLegend* legend = new TLegend( 0.2, 0.72, 0.5, 0.9 );
  legend->SetTextSize( 0.035 );
  legend->SetFillColor( 0 );
  legend->AddEntry( templ_sr,   "#sigma < 0.01", "L" );
  legend->AddEntry( templ_side, "#sigma > 0.01", "L" );
  //legend->AddEntry( templ_sideCorr, "#sigma > 0.01 + Corr.", "L" );
  legend->Draw("same");

  TPaveText* labelTop = MT2DrawTools::getLabelTop();
  labelTop->Draw("same");

  c1->SaveAs( Form("%s/closure.eps", outputdir.c_str()) );
  c1->SaveAs( Form("%s/closure.pdf", outputdir.c_str()) );
  c1->SaveAs( Form("%s/closure.png", outputdir.c_str()) );

  delete c1;
  delete h2_axes;

}




MT2Analysis<MT2EstimateTree> computeYield( const MT2Sample& sample, const std::string& regionsSet, const std::string& prompt_fake ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("mt2");
  
  std::cout << "-> Loaded tree: it has " << tree->GetEntries() << " entries." << std::endl;



  MT2Analysis<MT2EstimateTree> analysis( sample.sname, regionsSet, sample.id );
  MT2EstimateTree::addVar( &analysis, "iso" );
  MT2EstimateTree::addVar( &analysis, "isoAbs" );
  MT2EstimateTree::addVar( &analysis, "ptGamma" );
  MT2EstimateTree::addVar( &analysis, "sietaieta" );

  
  MT2Tree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);

  int nentries = tree->GetEntries();




  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);

    //// template control region (use low-significance regions):
    //if( myTree.gamma_ht > 1000. ) continue;
    if( myTree.gamma_mt2 < 200.) continue;
    //if( myTree.gamma_mt2 > 300. ) continue;
    if( myTree.mt2 > 200.) continue;
    //if( myTree.gamma_nJet40 > 4) continue;
    if( myTree.gamma_nBJet40 > 1) continue;
    //if( myTree.gamma_ht>1000. && sample.id==204 ) continue; // remove high-weight spikes (remove GJet_400to600 leaking into HT>1000)



    if( myTree.nMuons10 > 0) continue;
    if( myTree.nElectrons10 > 0 ) continue;
    if( myTree.nPFLep5LowMT > 0) continue;
    if( myTree.nPFHad10LowMT > 0) continue;

    if( myTree.gamma_deltaPhiMin<0.3 ) continue;
    if( myTree.gamma_diffMetMht>0.5*myTree.gamma_met_pt ) continue;
  
    if( myTree.nVert==0 ) continue;

    if( myTree.gamma_nJet40<2 ) continue;

    if( myTree.ngamma==0 ) continue;
    if( myTree.gamma_pt[0]<160. ) continue;




    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );
    int closestJet = -1;
    float deltaRmin = 0.4;
    for( unsigned i=0; i<myTree.njet; ++i ) {
      if( fabs(myTree.jet_eta[i])>2.5 ) continue;
      if( myTree.jet_pt[i]<40. ) continue;
      TLorentzVector thisjet;
      thisjet.SetPtEtaPhiM( myTree.jet_pt[i], myTree.jet_eta[i], myTree.jet_phi[i], myTree.jet_mass[i] );
      float thisDeltaR = gamma.DeltaR(thisjet);
      if( thisDeltaR<deltaRmin ) {
        deltaRmin = thisDeltaR;
        closestJet = i;
      }
    }
    float found_pt = 0.;
    int jet_counter = 0;
    for( unsigned i=0; i<myTree.njet; ++i ) {
      if( i==closestJet ) continue;
      if( fabs(myTree.jet_eta[i])>2.5 ) continue;
      if( myTree.jet_pt[i]<40. ) continue;
      jet_counter++;
      if( jet_counter==2 ) {
        found_pt = myTree.jet_pt[i];
        break;
      }
    }

    if( found_pt<100. ) continue;


    float sietaieta = myTree.gamma_sigmaIetaIeta[0];

    // for now only EB:
    if( fabs(gamma.Eta())>1.479 ) continue;
    // safety upper thresh on sietaieta:
    if( sietaieta>0.02 ) continue;
    float iso = myTree.gamma_chHadIso[0]/myTree.gamma_pt[0];
    if( iso > 0.1 ) continue;


    // remove prompt photons from QCD (remove double counting)
    if( sample.id>=100 && sample.id<199 ) {
      int mcMatchId = myTree.gamma_mcMatchId[0];
      float genIso = myTree.gamma_genIso[0];
      if((mcMatchId==22 || mcMatchId==7) && genIso<5.) continue;
    }
    


    Double_t weight = myTree.evt_scale1fb*lumi; 

    MT2EstimateTree* thisEstimate = analysis.get( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt );
    if( thisEstimate==0 ) continue;

    thisEstimate->yield->Fill(myTree.gamma_mt2, weight );


    thisEstimate->assignTree(myTree, lumi*myTree.evt_scale1fb);
    thisEstimate->assignVars( myTree.gamma_ht, myTree.gamma_nJet40, myTree.gamma_nBJet40, myTree.gamma_met_pt, myTree.gamma_mt2 );
    thisEstimate->assignVar( "sietaieta", sietaieta );
    thisEstimate->assignVar( "iso", iso );
    thisEstimate->assignVar( "isoAbs", iso*gamma.Pt() );
    thisEstimate->assignVar( "ptGamma", gamma.Pt() );
    thisEstimate->tree->Fill();

    
  } // for entries


  analysis.finalize();
  

  delete tree;


  file->Close();
  delete file;
  
  return analysis;

}





