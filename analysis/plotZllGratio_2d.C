{

  //READ in all the trees
  std::cout << "Reading in all input trees..." << std::endl;
  
  // Photon data
  TFile* fGdata = TFile::Open("EventYields_data_Run2016_12p9ifb/gammaControlRegion/data.root");
  TTree* tGdata = (TTree*) fGdata->Get("gammaCRtree/HT200toInf_j1toInf_b0toInf/tree_gammaCRtree_HT200toInf_j1toInf_b0toInf");

  // Photon MC
  TFile* fGmc = TFile::Open("EventYields_data_Run2016_12p9ifb/gammaControlRegion/mc.root");
  TTree* tGmc = (TTree*) fGmc->Get("gammaCRtree/HT200toInf_j1toInf_b0toInf/tree_gammaCRtree_HT200toInf_j1toInf_b0toInf");

  // Zll data
  TFile* fZdata = TFile::Open("EventYields_data_Run2016_12p9ifb/zllControlRegion/data.root");
  TTree* tZdata = (TTree*) fZdata->Get("data/HT200toInf_j1toInf_b0toInf/tree_data_HT200toInf_j1toInf_b0toInf");

  // Zll MC
  TFile* fZmc = TFile::Open("EventYields_data_Run2016_12p9ifb/zllControlRegion/mc.root");
  TTree* tZmc = (TTree*) fZmc->Get("zllCR/HT200toInf_j1toInf_b0toInf/tree_zllCR_HT200toInf_j1toInf_b0toInf");
  
  // Zll Purity (top)
  TFile* fZtop = TFile::Open("EventYields_data_Run2016_12p9ifb/zllControlRegion/ZllPurityTrees.root");
  TTree* tZtop = (TTree*) fZtop->Get("Top/HT200toInf_j1toInf_b0toInf/tree_Top_HT200toInf_j1toInf_b0toInf");

  std::cout << "All the trees have been read." << std::endl;

  std::cout << "Now initializing histograms..." << std::endl;
  int nBinsHT = 5;
  double binsHT[]={200, 450, 575, 1000, 1500, 3000};
  
  int nBinsJ = 4;
  double binsJ[]={1, 2, 4, 7, 12};

  int nBinsB = 4;
  double binsB[]={0, 1, 2, 3, 6};
  
  
  TH2D* hgammaData_jht = new TH2D("hgammaData_jht", "", nBinsHT, binsHT, nBinsJ, binsJ);
  TH2D* hgammaData_bj  = new TH2D("hgammaData_bj" , "", nBinsJ, binsJ, nBinsB, binsB);
  
  TH2D* hgammaMC_jht   = new TH2D("hgammaMC_jht", "", nBinsHT, binsHT, nBinsJ, binsJ);
  TH2D* hgammaMC_bj    = new TH2D("hgammaMC_bj" , "", nBinsJ, binsJ, nBinsB, binsB);
  
  TH2D* hgammaMCP_jht  = new TH2D("hgammaMCP_jht", "", nBinsHT, binsHT, nBinsJ, binsJ);
  TH2D* hgammaMCP_bj   = new TH2D("hgammaMCP_bj" , "", nBinsJ, binsJ, nBinsB, binsB);

  TH2D* hgammaMCL_jht  = new TH2D("hgammaMCL_jht", "", nBinsHT, binsHT, nBinsJ, binsJ);
  TH2D* hgammaMCL_bj   = new TH2D("hgammaMCL_bj" , "", nBinsJ, binsJ, nBinsB, binsB);
  

  TH2D* hzllData_jht = new TH2D("hzllData_jht", "", nBinsHT, binsHT, nBinsJ, binsJ);
  TH2D* hzllData_bj  = new TH2D("hzllData_bj" , "", nBinsJ, binsJ, nBinsB, binsB);
  
  TH2D* hzllMC_jht   = new TH2D("hzllMC_jht", "", nBinsHT, binsHT, nBinsJ, binsJ);
  TH2D* hzllMC_bj    = new TH2D("hzllMC_bj" , "", nBinsJ, binsJ, nBinsB, binsB);
  
  TH2D* hzllTop_jht  = new TH2D("hzllTop_jht", "", nBinsHT, binsHT, nBinsJ, binsJ);
  TH2D* hzllTop_bj   = new TH2D("hzllTop_bj" , "", nBinsJ, binsJ, nBinsB, binsB);

  
  std::cout << "Projecting gamma..." << std::endl;

  tGdata->Project("hgammaData_jht", "nJets:ht"    , "iso<2.5  && ((ht>200 && met>200)||(ht>1000. && met>30.)) && ptGamma>180 && nJets>0 && mt2>200 && ht>200");
  tGdata->Project("hgammaData_bj" , "nBJets:nJets", "iso<2.5  && ((ht>200 && met>200)||(ht>1000. && met>30.)) && ptGamma>180 && nJets>0 && mt2>200 && ht>200");
  
  tGmc->Project("hgammaMC_jht", "nJets:ht"    , "weight*(prompt==2 && iso<2.5  && ((ht>200 && met>200)||(ht>1000. && met>30.)) &&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23");
  tGmc->Project("hgammaMC_bj" , "nBJets:nJets", "weight*(prompt==2 && iso<2.5  && ((ht>200 && met>200)||(ht>1000. && met>30.)) &&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23");
  
  tGmc->Project("hgammaMCP_jht", "nJets:ht"    , "weight*(prompt>0 && iso<2.5  && ((ht>200 && met>200)||(ht>1000. && met>30.)) &&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23");
  tGmc->Project("hgammaMCP_bj" , "nBJets:nJets", "weight*(prompt>0 && iso<2.5  && ((ht>200 && met>200)||(ht>1000. && met>30.)) &&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23");
  
  tGmc->Project("hgammaMCL_jht", "nJets:ht"    , "weight*(iso<2.5  && ((ht>200 && met>200)||(ht>1000. && met>30.)) &&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23");
  tGmc->Project("hgammaMCL_bj" , "nBJets:nJets", "weight*(iso<2.5  && ((ht>200 && met>200)||(ht>1000. && met>30.)) &&  ptGamma>180 && nJets>0 && mt2>200 && ht>200 )*1.23");
  
  hgammaMCP_jht->Divide(hgammaMCL_jht);
  hgammaMCP_bj ->Divide(hgammaMCL_bj);
  
  hgammaData_jht->Multiply(hgammaMCP_jht);
  hgammaData_bj ->Multiply(hgammaMCP_bj);

  TH2D* purity_jht = (TH2D*) hgammaMCP_jht->Clone("purity_jht");
  TH2D* purity_bj = (TH2D*) hgammaMCP_bj->Clone("purity_bj");

  hgammaData_jht->Scale(0.92);
  hgammaData_bj ->Scale(0.92);

  std::cout << "Projected gamma." << std::endl;
  
  
  std::cout << "Projecting Z(ll)..." << std::endl;
  
  tZdata->Project("hzllData_jht", "nJets:ht"    , "(abs(Z_mass-91.19)<10 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && mt2>200 && ht>200 && nJets>0 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)))");
  tZdata->Project("hzllData_bj" , "nBJets:nJets", "(abs(Z_mass-91.19)<10 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && mt2>200 && ht>200 && nJets>0 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)))");
  
  tZmc->Project("hzllMC_jht", "nJets:ht"    , "weight*HLT_weight*weight_lep0*(abs(Z_mass-91.19)<10 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && mt2>200 && ht>200 && nJets>0 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)))");
  tZmc->Project("hzllMC_bj" , "nBJets:nJets", "weight**HLT_weight*weight_lep0*(abs(Z_mass-91.19)<10 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && mt2>200 && ht>200 && nJets>0 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)))");
  
  tZtop->Project("hzllTop_jht", "nJets:ht"    , "weight*(abs(Z_mass-91.19)<10 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && mt2>200 && ht>200 && nJets>0 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)))");
  tZtop->Project("hzllTop_bj" , "nBJets:nJets", "weight*(abs(Z_mass-91.19)<10 && ((ht>200 && met>200)||(ht>1000. && met>30.)) && mt2>200 && ht>200 && nJets>0 && (Z_lepId==13 || (Z_lepId==11 && lep_tightId0>0 && lep_tightId1>0)))");
  
  TH2D* zllPurity_jht = (TH2D*) hzllMC_jht->Clone("zllPurity_jht");
  zllPurity_jht->Add(hzllTop_jht);
  zllPurity_jht->Divide(hzllMC_jht);  

  TH2D* zllPurity_bj = (TH2D*) hzllMC_bj->Clone("zllPurity_bj");
  zllPurity_bj->Add(hzllTop_bj);
  zllPurity_bj->Divide(hzllMC_bj);  

  hzllData_jht->Divide(zllPurity_jht);
  hzllData_bj ->Divide(zllPurity_bj );

  std::cout << "Projeced Z(ll)." << std::endl;

  std::cout << "Performing Z/G ratio in data..." << std::endl;

  hzllData_jht->Divide(hgammaData_jht);
  hzllData_bj ->Divide(hgammaData_bj);

  TH2D* zllG_data_jht = (TH2D*) hzllData_jht->Clone("zllG_data_jht");
  TH2D* zllG_data_bj  = (TH2D*) hzllData_bj ->Clone("zllG_data_bj");

  std::cout << "Performing Z/G ratio in MC..." << std::endl;
  
  hzllMC_jht->Divide(hgammaMC_jht);
  hzllMC_bj ->Divide(hgammaMC_bj);

  TH2D* zllG_mc_jht = (TH2D*) hzllMC_jht->Clone("zllG_mc_jht");
  TH2D* zllG_mc_bj  = (TH2D*) hzllMC_bj ->Clone("zllG_mc_bj");

  std::cout << "Performing double ratio..." << std::endl;
  
  hzllData_jht->Divide(hzllMC_jht);
  hzllData_bj ->Divide(hzllMC_bj);

  std::cout << "All ratios have been calculated. Drawing..." << std::endl;
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(51);

  TCanvas* cp0 = new TCanvas("cp0", "", 600, 600);
  purity_jht->SetTitle("Photon Purity");
  purity_jht->GetXaxis()->SetTitle("H_{T} [GeV]");
  purity_jht->GetYaxis()->SetTitle("Jet Multiplicity");
  purity_jht->GetZaxis()->SetRangeUser(0.5, 1.0);
  purity_jht->GetZaxis()->SetLabelSize(0.02);
  purity_jht->Draw("colz");

  TCanvas* cp1 = new TCanvas("cp1", "", 600, 600);
  purity_bj->SetTitle("Photon Purity");
  purity_bj->GetXaxis()->SetTitle("Jet Multiplicity");
  purity_bj->GetYaxis()->SetTitle("b-jet Multiplicity");
  purity_bj->GetZaxis()->SetRangeUser(0.5, 1.0);
  purity_bj->GetZaxis()->SetLabelSize(0.02);
  purity_bj->Draw("colz");

  TCanvas* c0 = new TCanvas("c0", "", 600, 600);
  hzllData_jht->SetTitle("Z(l^{+}l^{-}) / #gamma ratio");
  hzllData_jht->GetXaxis()->SetTitle("H_{T} [GeV]");
  hzllData_jht->GetYaxis()->SetTitle("Jet Multiplicity");
  hzllData_jht->GetZaxis()->SetTitle("Data / MC");
  hzllData_jht->GetZaxis()->SetRangeUser(0.1,1.25);
  hzllData_jht->GetZaxis()->SetTitleOffset(0.6);
  hzllData_jht->GetZaxis()->SetLabelSize(0.02);
  hzllData_jht->Draw("colz");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  hzllData_bj->SetTitle("Z(l^{+}l^{-}) / #gamma ratio");
  hzllData_bj->GetXaxis()->SetTitle("Jet Multiplicity");
  hzllData_bj->GetYaxis()->SetTitle("b-jet Multiplicity");
  hzllData_bj->GetZaxis()->SetTitle("Data / MC");
  hzllData_bj->GetZaxis()->SetRangeUser(0.1,1.25);
  hzllData_bj->GetZaxis()->SetTitleOffset(0.6);
  hzllData_bj->GetZaxis()->SetLabelSize(0.02);
  hzllData_bj->Draw("colz");

}
