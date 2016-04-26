#include "interface/MT2DrawTools.h"
#include "interface/MT2Config.h"

#include "TH1F.h"
#include "TChain.h"
#include "TLegend.h"

void setBranches(TChain *c, std::vector<std::string> vars, std::map<std::string, float> &fvars, int &njets, int &nbjets);
void setHistos  (std::vector<std::string> vars, std::map<std::string, TH1F*> &hvars, std::string str="");
void fillHistos (float w, std::vector<std::string> vars, std::map< std::string, float > fVars, std::map< std::string, TH1F* > &hVars, std::string prefix);
void drawHistos (std::vector<std::string> vars, std::map< std::string, TH1F* > hVars, std::string str="");
void saveHistos (TFile *file, std::map< std::string, TH1F* > hVars);
int   getNbins (std::string var);
float getBinMax(std::string var);
float getBinMin(std::string var);


int main( int argc, char* argv[] ) {



  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|               closure of R&S              |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "|                                           |" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;

  if( argc<2 ) {
    std::cout << "USAGE: ./closureRS configFileName [data/MC/all] [sampleID] [job_i] [Njobs]" << std::endl
	      << "Exiting." << std::endl;
    exit(11);
  }

  
  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  bool onlyData = false;
  bool onlyMC   = false;
  int ijob=0, Njobs=1;
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

  MT2DrawTools::setStyle();

  int nSmearings = 100;

  bool smearGen = false;

  TString dcapPath = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/20apr" + std::string(smearGen ? "_gen/" : "/") + cfg.getEventYieldDir() + "/smearedTrees/";
  TString treename = "qcdSmearedTree/HT0toInf_j1toInf_b0toInf/tree_qcdSmearedTree_HT0toInf_j1toInf_b0toInf";

  std::vector<std::string> vars = {"ht", "met", "mht", "mt2", "nJets", "nBJets", "jet1_pt", "jet2_pt", "deltaPhiMin", "diffMetMht"};
  std::map< std::string, float > fVars;
  std::map< std::string, TH1F* > hVars;
  std::map< std::string, TH1F* > hVars_hHT;
  std::map< std::string, TH1F* > hVars_vhHT;
  std::map< std::string, TH1F* > hVars_lHT;
  setHistos(vars, hVars);
  setHistos(vars, hVars_hHT ,"_hHT" );
  setHistos(vars, hVars_vhHT,"_vhHT");
  setHistos(vars, hVars_lHT ,"_lHT" );

  TChain *ch = new TChain(treename);
  ch->Add(dcapPath+TString::Format("mc_id%d_job%dof%d.root", sampleID, ijob, Njobs));
  // for (int sampleID=156; sampleID<=156; sampleID++)
  //   for (int j=0; j<Njobs; j++){
  //     int job = j;
  //     // this is bullshit, be careful
  //     //if(sampleID==153 && j==10)  job=j-1;
  //     //if(j==10)  continue;
  //     if(j==4||j==5||j==14)  continue;
  //     //ch->Add(dcapPath+TString::Format("mc_id%d_job%dof25.root", sampleID, job));
  //     ch->Add(dcapPath+TString::Format("mc_id%d_job%dof30.root", sampleID, job));
  //   }
  // //ch->Add(dcapPath+TString::Format("mc_id%d_job%dof%d.root", sampleID, j, Njobs));
  
  int njets, nbjets; // the ones w/o "before_" are ints!!
  setBranches(ch,vars,fVars, njets, nbjets);


  ch->SetBranchStatus("iSmear"    ,1);
  ch->SetBranchStatus("weight"    ,1);
  ch->SetBranchStatus("nElectrons",1);
  ch->SetBranchStatus("nMuons"    ,1);
  ch->SetBranchStatus("nPFLep"    ,1);
  ch->SetBranchStatus("nPFHad"    ,1);
  float iSmear;  // it's a float!!
  ch->SetBranchAddress("iSmear", &iSmear);
  float weight;
  ch->SetBranchAddress("weight", &weight);
  int nele, nmuo, npfhad, npflep;
  ch->SetBranchAddress( "nElectrons"   , &nele   );
  ch->SetBranchAddress( "nMuons"       , &nmuo   );
  ch->SetBranchAddress( "nPFLep"       , &npfhad );
  ch->SetBranchAddress( "nPFHad"       , &npflep );
  
	     
  unsigned long long nentries = ch->GetEntries();

  std::cout << "nentries: " << nentries << std::endl;

  //nentries = 100000;

  for (unsigned long long ientry=0; ientry<nentries; ientry++){

    if( ientry % 1000000 == 0 ) std::cout << "    Entry: " << ientry << " / " << nentries << std::endl;
    ch->GetEntry(ientry);
    fVars["nJets" ] = (float)njets ;  // because they are ints
    fVars["nBJets"] = (float)nbjets;
    
    
    if (nele+nmuo+npfhad+npflep>0)  continue;


    if (iSmear==1){
      if (fVars["before_ht"]>200)  {
	fillHistos(weight, vars, fVars, hVars, "before_");
	if (fVars["before_ht"]>3000)
	  fillHistos(weight, vars, fVars, hVars_vhHT, "before_");
	if (fVars["before_ht"]>1000)
	  fillHistos(weight, vars, fVars, hVars_hHT, "before_");
	else
	  fillHistos(weight, vars, fVars, hVars_lHT, "before_");
      }
    }

    if (iSmear<=nSmearings){
      if (fVars["ht"]>200) { 
	fillHistos(weight/nSmearings, vars, fVars, hVars,"");
	if (fVars["ht"]>3000)  
	  fillHistos(weight/nSmearings, vars, fVars, hVars_vhHT, "");
	if (fVars["ht"]>1000)  
	  fillHistos(weight/nSmearings, vars, fVars, hVars_hHT, "");
	else
	  fillHistos(weight/nSmearings, vars, fVars, hVars_lHT, "");
      }
    }

  }
  std::cout << " ---> Loop finished" << std::endl;

  // drawHistos(vars, hVars);
  // drawHistos(vars, hVars_hHT,"_hHT");
  // drawHistos(vars, hVars_lHT,"_lHT");

  TString outputdir =  cfg.getEventYieldDir() + "/closureRS" + (smearGen ? "_gen/" : "/"); 
  system(Form("mkdir -p %s", outputdir.Data()));
  TFile *outFile = new TFile(outputdir+TString::Format("histos_id%d_job%dof%d.root", sampleID, ijob, Njobs),"RECREATE");

  saveHistos(outFile, hVars);
  saveHistos(outFile, hVars_vhHT);
  saveHistos(outFile, hVars_hHT);
  saveHistos(outFile, hVars_lHT);
  outFile->Close();

  return 0;

}

void setHistos(std::vector<std::string> vars, std::map<std::string, TH1F* > &hvars, std::string str){

  for (unsigned int v=0; v<vars.size();v++){
    hvars[          vars.at(v)] = new TH1F((          vars.at(v)+str.c_str()) .c_str(),"", getNbins(vars.at(v)), getBinMin(vars.at(v)), getBinMax(vars.at(v)));
    hvars["before_"+vars.at(v)] = new TH1F(("before_"+vars.at(v)+str.c_str()).c_str(),"", getNbins(vars.at(v)), getBinMin(vars.at(v)), getBinMax(vars.at(v)));
    hvars[          vars.at(v)]->Sumw2();
    hvars["before_"+vars.at(v)]->Sumw2();
    hvars[          vars.at(v)]->SetLineWidth(2);    hvars[          vars.at(v)]->SetLineColor(2);
    hvars["before_"+vars.at(v)]->SetLineWidth(2);    hvars["before_"+vars.at(v)]->SetLineColor(4);
    hvars[          vars.at(v)]->SetXTitle(vars.at(v).c_str());
    hvars["before_"+vars.at(v)]->SetXTitle(vars.at(v).c_str());
  }

}

void setBranches(TChain *c, std::vector<std::string> vars, std::map<std::string, float> &fvars, int &njets, int &nbjets){
  c->SetBranchStatus("*",0);
  for (unsigned int v=0; v<vars.size();v++){
    c->SetBranchStatus(           vars.at(v) .c_str() , 1 );
    c->SetBranchStatus(("before_"+vars.at(v)).c_str() , 1 );
    
    fvars[          vars.at(v)] = 0.0;
    fvars["before_"+vars.at(v)] = 0.0;

    c->SetBranchAddress(("before_"+vars.at(v)).c_str(), &(fvars["before_"+vars.at(v)]) );
    if (vars.at(v) == "nJets")
      c->SetBranchAddress( vars.at(v).c_str(), &njets );
      else if (vars.at(v) == "nBJets")
      c->SetBranchAddress( vars.at(v).c_str(), &nbjets );
    else
      c->SetBranchAddress( vars.at(v).c_str(), &(fvars[vars.at(v)]) );
  }

}

void fillHistos(float w, std::vector<std::string> vars, std::map< std::string, float > fVars, std::map< std::string, TH1F* > &hVars, std::string prefix){
   
  for (unsigned int v=0; v<vars.size();v++){
    std::string var = vars.at(v);
    if (var=="nJets" || fVars[prefix+"nJets"] >=2){
      float val = fVars[prefix+var];
      //if(prefix=="before_" && var=="ht") std::cout << val << std::endl;
      if (var=="diffMetMht") val /= fVars[prefix+"met"];
      hVars[prefix+var]->Fill(val, w);
    }
  }
  
}

// void drawHistos(std::vector<std::string> vars, std::map< std::string, TH1F* > hVars){
  

//   TCanvas *c = new TCanvas("c", "", 600, 600);
//   c->SetLogy();
//   for (unsigned int v=0; v<vars.size();v++){
//     std::string var = vars.at(v);
//     hVars["before_"+var]->Draw();
//     hVars[          var]->Draw("same");
//     c->SaveAs(("smearPlots/" + var + ".pdf").c_str());
//   }
  
// }


void drawHistos(std::vector<std::string> vars, std::map< std::string, TH1F* > hVars, std::string str){
  TCanvas* c = new TCanvas( "c1", "", 600, 700 );
  c->cd(); 
  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->SetLogy();
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();

  TLegend *leg = new TLegend(0.75, 0.8, 0.9, 0.9);
  leg->SetTextSize(0.038);
  leg->SetTextFont(42);
  leg->SetFillColor(0);
  leg->AddEntry(hVars["before_nJets"], "QCD MC", "PL");
  leg->AddEntry(hVars[       "nJets"], "QCD RS", "PL");
  
  for (unsigned int v=0; v<vars.size();v++){
    pad1->cd();
    std::string var = vars.at(v);
    float xmin = hVars["before_"+var]->GetXaxis()->GetXmin();
    float xmax = hVars["before_"+var]->GetXaxis()->GetXmax();
    hVars["before_"+var]->Draw();
    hVars[          var]->Draw("same");
    leg->Draw("same");

    pad2->cd();
    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xmin, xmax, 0.0, 2.0 );
    h2_axes_ratio->SetYTitle("MC/RS");
    h2_axes_ratio->Draw();
    TH2F* h_ratio = (TH2F*)hVars["before_"+var]->Clone((var+"_ratio").c_str());
    h_ratio->Divide(hVars[var]);
    h_ratio->Draw("same");

    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    line->Draw("same");

    c->SaveAs(("smearPlots/" + var + str + ".pdf").c_str());
    c->SaveAs(("smearPlots/" + var + str + ".C").c_str());
    c->SaveAs(("smearPlots/" + var + str + ".root").c_str());
    delete h2_axes_ratio;
    delete line;
  }

}

void saveHistos(TFile *file, std::map< std::string, TH1F* > hVars){

  file->cd();
  for (std::map< std::string, TH1F* >::iterator iH=hVars.begin(); iH!=hVars.end(); ++iH)
    iH->second->Write();

}


int   getNbins (std::string var){
  if (var.find("nJets") != std::string::npos )
    return 15;
  if (var.find("nBJets") != std::string::npos )
    return 8;
  if (var.find("delta") != std::string::npos )
    return 64;
  else 
    return 100;
}

float getBinMax(std::string var){
  if (var.find("nJets") != std::string::npos )
    return 15.;
  if (var.find("nBJets") != std::string::npos )
    return 8.;
  if (var.find("jet") != std::string::npos )
    return 2000.;
  if (var.find("diff") != std::string::npos )
    return 10.;
  if (var.find("met") != std::string::npos )
    return 1000.;
  if (var.find("mht") != std::string::npos )
    return 1000.;
  if (var.find("mt2") != std::string::npos )
    return 600.;
  if (var.find("delta") != std::string::npos )
    return 3.2;
  if (var.find("ht") != std::string::npos )
    return 6000.;
  else 
    return 100;

}

float getBinMin(std::string var){
  if (var == "ht" )
    return 200.;
  else 
    return 0.;
}
