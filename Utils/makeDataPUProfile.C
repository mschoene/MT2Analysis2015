{

  TChain* chain_ = new TChain("mt2");
  chain_->Add("/scratch/pandolf/data2015_50ns_DCSonly/JetHT_Run2015B.root");
  
  TH1D* hPU_data = new TH1D("hPU_data", "", 100, 0, 100);
  chain_->Project("hPU_data", "nVert");
  
  ULong64_t nEntries = chain_->GetEntries();
  hPU_data->Scale(1.0/nEntries);
  
  TFile* output = TFile::Open("/shome/mmasciov/PUProfile_JetHT.root", "recreate");
  output->cd();
  hPU_data->Write();
  output->Close();

}
