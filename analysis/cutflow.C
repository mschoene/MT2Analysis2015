{

  TChain* chain = new TChain("mt2"); 
  //chain->Add("");

  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/QCD_Pt1000to1400_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/QCD_Pt1400to1800_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/QCD_Pt1800to2400_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/QCD_Pt2400to3200_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/QCD_Pt300to470_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/QCD_Pt3200_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/QCD_Pt470to600_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/QCD_Pt600to800_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/QCD_Pt800to1000_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/ZJetsToNuNu_HT100to200_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/ZJetsToNuNu_HT200to400_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/ZJetsToNuNu_HT400to600_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/ZJetsToNuNu_HT600toInf_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/WJetsToLNu_HT100to200_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/WJetsToLNu_HT200to400_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/WJetsToLNu_HT400to600_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/WJetsToLNu_HT600toInf_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/TBarToLeptons_sch_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/TBarToLeptons_tch_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/TBar_tWch_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/TTH_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/TTJets_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/TTWJets_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/TTZJets_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/TToLeptons_sch_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/TToLeptons_tch_post_skim_prune.root");
  //chain->Add("$DCAPUSER//pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/T_tWch_post_skim_prune.root");
  //////
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T1bbbb_2J_mGl1000_mLSP900_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T1bbbb_2J_mGl1500_mLSP100_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T1qqqq_2J_mGl1000_mLSP800_post_skim_prune.root");
  //chain->Add("$DCAPUSER/mmasciov/MT2production/PostProcessed/01Jul2015_Signal/skimAndPrune/SMS_T1qqqq_2J_mGl1400_mLSP100_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T1tttt_2J_mGl1200_mLSP800_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T1tttt_2J_mGl1500_mLSP100_post_skim_prune.root");
  //chain->Add("$DCAPUSER/mmasciov/MT2production/PostProcessed/01Jul2015_Signal/skimAndPrune/SMS_T2bb_2J_mStop600_mLSP580_post_skim_prune.root");
  //chain->Add("$DCAPUSER/mmasciov/MT2production/PostProcessed/01Jul2015_Signal/skimAndPrune/SMS_T2bb_2J_mStop900_mLSP100_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T2qq_2J_mStop1200_mLSP100_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T2qq_2J_mStop600_mLSP550_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T2tt_2J_mStop425_mLSP325_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T2tt_2J_mStop500_mLSP325_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T2tt_2J_mStop650_mLSP325_post_skim_prune.root");
  //chain->Add("$DCAPUSER/pandolf/MT2production/PostProcessed/June15_jet30_v4/skimAndPrune/SMS_T2tt_2J_mStop850_mLSP100_post_skim_prune.root");
  //
  //chain->Add("/scratch/pandolf/data2015_50ns_DCSonly/JetHT_Run2015B.root");
  //
  chain->Add("/scratch/mmasciov/data2015_HBHEfilter_final/JetHT_Run2015B.root");

  float jet_pt[100];
  float jet_eta[100];
  float jet_phi[100];
  float jet_mass[100];
  float met_pt;
  float met_phi;
  float ht;
  float mt2;
  float diffMetMht;
  float deltaPhiMin;

  int nVert;
  int njet;
  int nJet30;
  int nElectrons10;
  int nMuons10;
  int nPFLep5LowMT;
  int nPFHad10LowMT;
  
  float evt_scale1fb=1;

//  float evt_scale1fb;
//  chain->SetBranchAddress("evt_scale1fb", &evt_scale1fb);
  
  chain->SetBranchAddress("nVert", &nVert);
  chain->SetBranchAddress("njet", &njet);
  chain->SetBranchAddress("nJet30", &nJet30);
  chain->SetBranchAddress("nElectrons10", &nElectrons10);
  chain->SetBranchAddress("nMuons10", &nMuons10);
  chain->SetBranchAddress("nPFLep5LowMT", &nPFLep5LowMT);
  chain->SetBranchAddress("nPFHad10LowMT", &nPFHad10LowMT);

  chain->SetBranchAddress("ht", &ht);
  chain->SetBranchAddress("mt2", &mt2);
  chain->SetBranchAddress("met_pt", &met_pt);
  chain->SetBranchAddress("met_phi", &met_phi);
  chain->SetBranchAddress("deltaPhiMin", &deltaPhiMin);
  chain->SetBranchAddress("diffMetMht", &diffMetMht);

  chain->SetBranchAddress("jet_pt", jet_pt);
  chain->SetBranchAddress("jet_eta", jet_eta);
  chain->SetBranchAddress("jet_phi", jet_phi);
  chain->SetBranchAddress("jet_mass", jet_mass);

  TH1F* hdeltaPhiMax[2];

  hdeltaPhiMax[0]= new TH1F("hdeltaPhiMax_all", "", 320, TMath::Pi()-3.2, TMath::Pi());
  hdeltaPhiMax[0]->GetXaxis()->SetTitle("max #Delta #phi");
  hdeltaPhiMax[0]->SetLineColor(1);
  hdeltaPhiMax[0]->SetMarkerColor(1); 

  hdeltaPhiMax[1]= new TH1F("hdeltaPhiMax_pass", "", 320, TMath::Pi()-3.2, TMath::Pi());
  hdeltaPhiMax[1]->GetXaxis()->SetTitle("max #Delta #phi");
  hdeltaPhiMax[1]->SetLineColor(2);
  hdeltaPhiMax[1]->SetMarkerColor(2);

  
  float nTotal=0;
  float nPass=0;
  float nPass_=0;
  
  std::cout << "Starting loop over " << chain->GetEntries() << " entries..." << std::endl;
  for( unsigned int i=0; i < chain->GetEntries(); ++i ){
      
    chain->GetEntry(i);
    
    nTotal+=evt_scale1fb;

    float deltaPhiMax=0.;
       
    TLorentzVector met;
    met.SetPtEtaPhiM(met_pt, 0.0, met_phi, 0.0);
    
    int maxNJ = 0;
    for( int j=0; j < njet; ++j ){

      if( maxNJ > 4 ) break;
      
      if( jet_pt[j] < 30. || jet_eta[j] > 4.7 ) continue;
      
      TLorentzVector jet;
      jet.SetPtEtaPhiM(jet_pt[j], jet_eta[j], jet_phi[j], jet_mass[j]);
            
      float thisDeltaPhi = jet.DeltaPhi(met);
      deltaPhiMax = (fabs(thisDeltaPhi) > deltaPhiMax) ? fabs(thisDeltaPhi) : deltaPhiMax;

      ++maxNJ;

    }
    
    hdeltaPhiMax[0]->Fill(deltaPhiMax, evt_scale1fb);
    
    // Standard Selection
    if( ht < 450. ) continue;
    //if( ht < 900. ) continue;
    if( (ht > 450. && met_pt < 200.) || (ht > 900. && met_pt < 30.) ) continue;
    if( nJet30 < 2 ) continue;
    if( nElectrons10+nMuons10+nPFLep5LowMT+nPFHad10LowMT > 0 ) continue;
    if( deltaPhiMin < 0.3 ) continue;
    if( diffMetMht > 0.5*met_pt ) continue; 
    //    if( mt2 < 200 ) continue;

    if( met_pt < 200 || nJet30 < 2 ) continue;

    nPass_+=evt_scale1fb;

    hdeltaPhiMax[1]->Fill(deltaPhiMax, evt_scale1fb);
	
    //    if( deltaPhiMax > TMath::Pi() - 0.3 ) continue;
    if( deltaPhiMax > TMath::Pi() - 0.015 ) continue;
    else nPass+=evt_scale1fb;

  }

  std::cout << "Standard Selection: " << 1.0*nPass_/nTotal*100.0 << "% of events" << std::endl;
  std::cout << "max deltaPhi(4 jets, MET) < #Pi - 0.3: " << 1.0*nPass/nPass_*100.0 << "% of events after standard-selection" << std::endl;
  std::cout << "i.e.: " << 1.0*nPass/nTotal*100.0 << "% of all events" << std::endl;

  TFile *f=TFile::Open("deltaPhiMax_data_HBHEfilter_all.root", "recreate");
  f->cd();
  hdeltaPhiMax[0]->Write();
  hdeltaPhiMax[1]->Write();
  f->Close();

}
