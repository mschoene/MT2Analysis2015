#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2EstimateSigContSyst.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2DrawTools.h"
#include "interface/MT2Config.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>


#include "TMath.h"
#include "TRandom3.h"
#include "TTreeFormula.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"

#include "RooHistError.h"


struct BGTable {

  float zinv;
  float zinv_statUp;
  float zinv_statDn;
  float zinv_systUp;
  float zinv_systDn;

  float llep;
  float llep_statUp;
  float llep_statDn;
  float llep_systUp;
  float llep_systDn;

  float qcd;
  float qcd_statUp;
  float qcd_statDn;
  float qcd_systUp;
  float qcd_systDn;

};

bool drawSignals=false;
float lumi; //fb-1 

BGTable getTable( const std::string& tableFileName );
void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data, std::vector<MT2Analysis<MT2Estimate>*> analysesSignal, std::vector<MT2Analysis<MT2EstimateSigContSyst>*> analysesSignalCont, std::string dir );


int main( int argc, char* argv[] ) {
  
  
  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|           Running computeLostLepton                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;
  
  
  if( argc!=2 ) {
    std::cout << "USAGE: ./computeLostLepton [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  
  
  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);

  lumi = cfg.lumi();
  
  TH1::AddDirectory(kTRUE);
  
  std::string dir = cfg.getEventYieldDir();
  std::string outputdir = cfg.getEventYieldDir() + "/YieldComparison_dataMC_binned_post";
 
 
  MT2Analysis<MT2Estimate>* analysis = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "data" ); // any one is good, just need to know the regions                                                                    

  std::vector < MT2Analysis<MT2Estimate>* > analysesSignal;
  std::vector < MT2Analysis<MT2EstimateSigContSyst>* > analysesSignalCont;
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1tttt_1500_100") );
//  analysesSignal[0]->setName("T1tttt 1500,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1tttt_1200_800") );
//  analysesSignal[1]->setName("T1tttt 1200,800");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1bbbb_1500_100") );
//  analysesSignal[2]->setName("T1bbbb 1500,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1bbbb_1000_900") );
//  analysesSignal[3]->setName("T1bbbb 1000,900");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1qqqq_1400_100") );
//  analysesSignal[4]->setName("T1qqqq 1400,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "T1qqqq_1000_800") );
//  analysesSignal[5]->setName("T1qqqq 1000,800");


  std::string sigPath="./signalScansFromDominick";

//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1bbbb_eth.root", "T1bbbb") );
//  //  analysesSignal[0]->setName("T1bbbb 1500,100");
//  analysesSignal[0]->setName("pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow b#bar{b}#tilde{#chi}_{1}^{0}");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1bbbb_eth.root", "T1bbbb") );
//  analysesSignal[1]->setName("T1bbbb 700,600");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1qqqq_eth.root", "T1qqqq") );
//  analysesSignal[2]->setName("T1qqqq 1300,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1qqqq_eth.root", "T1qqqq") );
//  analysesSignal[3]->setName("T1qqqq 700, 600");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1tttt_eth.root", "T1tttt") );
//  //  analysesSignal[4]->setName("T1tttt 1200, 100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1tttt_eth.root", "T1tttt") );
//  //  analysesSignal[5]->setName("T1tttt 700,400");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2tt_eth.root", "T2tt") );
//  analysesSignal[6]->setName("T2tt 650, 0");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2tt_eth.root", "T2tt") );
//  analysesSignal[7]->setName("T2tt 600,200");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2tt_eth.root", "T2tt") );
//  analysesSignal[8]->setName("T2tt 200, 100");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1bbbb_eth.root", "T1bbbb") );
  //  analysesSignal[0]->setName("T1bbbb 1500, 100");
  analysesSignal[0]->setName("pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow b#bar{b}#tilde{#chi}_{1}^{0}");
  //  (*analysesSignal[0]) *= 2.26355/2.155;

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1bbbb_eth.root", "T1bbbb") );
  analysesSignal[1]->setName("T1bbbb 700, 600");
  //  (*analysesSignal[1]) *= 2.26355/2.155;

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1qqqq_eth.root", "T1qqqq") );
  analysesSignal[2]->setName("pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q}#tilde{#chi}_{1}^{0}");
  //  (*analysesSignal[2]) *= 2.26355/2.155;

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1qqqq_eth.root", "T1qqqq") );
  analysesSignal[3]->setName("T1qqqq 700, 600");
  //  (*analysesSignal[3]) *= 2.26355/2.155;

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2bb_eth.root", "T2bb") );
  analysesSignal[4]->setName("pp #rightarrow #tilde{b}#bar{#tilde{b}}, #tilde{b} #rightarrow b#tilde{#chi}_{1}^{0}");
  //  (*analysesSignal[4]) *= 2.26355/2.26;

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2bb_eth.root", "T2bb") );
  analysesSignal[5]->setName("T2bb 400, 200");
  //  (*analysesSignal[5]) *= 2.26355/2.26;

//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2qq_eth.root", "T2qq") );
//  analysesSignal[6]->setName("pp #rightarrow #tilde{q}#bar{#tilde{q}}, #tilde{q} #rightarrow q#tilde{#chi}_{1}^{0}");
//  (*analysesSignal[6]) *= 2.26355/2.26;
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2qq_eth.root", "T2qq") );
//  analysesSignal[7]->setName("T2qq 600, 0");
//  (*analysesSignal[7]) *= 2.26355/2.26;
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2qq_eth.root", "T2qq") );
//  analysesSignal[8]->setName("T2qq 500, 300");
//  (*analysesSignal[8]) *= 2.26355/2.26;


  analysesSignalCont.push_back( MT2Analysis<MT2EstimateSigContSyst>::readSystFromFile( sigPath + "/T1tttt_sigcontam_eth.root", "T1tttt_sigcontam", "isr") );
 
  analysesSignalCont.push_back( MT2Analysis<MT2EstimateSigContSyst>::readSystFromFile( sigPath + "/T1tttt_sigcontam_eth.root", "T1tttt_sigcontam", "isr") );
 
  analysesSignalCont.push_back( MT2Analysis<MT2EstimateSigContSyst>::readSystFromFile( sigPath + "/T2tt_sigcontam_eth.root", "T2tt_sigcontam", "isr") );
  analysesSignalCont[2]->setName("pp #rightarrow #tilde{t}#bar{#tilde{t}}, #tilde{t} #rightarrow t#tilde{#chi}_{1}^{0}");
 
  analysesSignalCont.push_back( MT2Analysis<MT2EstimateSigContSyst>::readSystFromFile( sigPath + "/T2tt_sigcontam_eth.root", "T2tt_sigcontam", "isr") );
  analysesSignalCont[3]->setName("T2tt 600, 200");
 
  analysesSignalCont.push_back( MT2Analysis<MT2EstimateSigContSyst>::readSystFromFile( sigPath + "/T2tt_sigcontam_eth.root", "T2tt_sigcontam", "isr") );
  analysesSignalCont[4]->setName("T2tt 200, 100");
  
  std::set<MT2Region> regions = analysis->getRegions();

  MT2Analysis<MT2Estimate>* data = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "data" );
  
  drawYields( outputdir.c_str(), data, analysesSignal, analysesSignalCont, dir );

  return 0;

}

void drawYields( const std::string& outputdir, MT2Analysis<MT2Estimate>* data,  std::vector < MT2Analysis<MT2Estimate>* > analysesSignal, std::vector < MT2Analysis<MT2EstimateSigContSyst>* > analysesSignalCont, std::string dir ) {

  
  MT2DrawTools::setStyle();

  system(Form("mkdir -p %s", outputdir.c_str()));

  std::vector<std::string> sigName;
  sigName.push_back("T1bbbb_1500_100");
  sigName.push_back("T1bbbb_700_600");
  sigName.push_back("T1qqqq_1300_100");
  sigName.push_back("T1qqqq_700_600");
  sigName.push_back("T2bb_700_100");
  sigName.push_back("T2bb_400_200");
//  sigName.push_back("T2qq_1000_100");
//  sigName.push_back("T2qq_600_0");
//  sigName.push_back("T2qq_500_300");
  sigName.push_back("T1tttt_1200_100");
  sigName.push_back("T1tttt_700_400");
  sigName.push_back("T2tt_650_100");
  sigName.push_back("T2tt_600_200");
  sigName.push_back("T2tt_200_100");

  std::vector<int> colors;
//  colors.push_back( 402 );
//  colors.push_back( 430 );
//  colors.push_back( 418 );
  
  colors.push_back( kYellow+1 );
  colors.push_back( kAzure+4 );
  colors.push_back( kGreen+2 );

  std::vector<int> colorsSig;
  colorsSig.push_back(2);
  colorsSig.push_back(2);
  colorsSig.push_back(kAzure+10);
  colorsSig.push_back(kAzure+10);
  colorsSig.push_back(2);
  colorsSig.push_back(2);
//  colorsSig.push_back(kAzure+10);
//  colorsSig.push_back(kAzure+10);
//  colorsSig.push_back(kAzure+10);
  colorsSig.push_back(6);
  colorsSig.push_back(6);
  colorsSig.push_back(6);
  colorsSig.push_back(6);
  colorsSig.push_back(6);

  std::vector<int> styleSig;
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
//  styleSig.push_back(1);
//  styleSig.push_back(1);
//  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);

  int bgSize = 3;
  int sigSize = 6;//9;
  int sigContSize = 5;

  int S=0;

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  TH1D* hdata = new TH1D("hdata", "", 174, 0, 174);
  hdata->Sumw2();
  hdata->GetYaxis()->SetTitle("Entries");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( 1 );
  hdata->SetMarkerColor( 1 );
  
  TH1D* hsig[sigSize+sigContSize];
  for(int s=0; s<sigSize+sigContSize; ++s){
    std::string thisNameS( Form("hsig_%d", s) );
    hsig[s] = new TH1D(thisNameS.c_str(), "", 174, 0, 174);
    hsig[s]->Sumw2();
    hsig[s]->GetYaxis()->SetTitle("Entries");
    hsig[s]->SetLineColor( colorsSig[s] );
    hsig[s]->SetFillColor( colorsSig[s] );
    hsig[s]->SetLineStyle( styleSig[s] );
  }

  TH1D* hestimate_all = new TH1D(Form("hestimate_all"), "", 174, 0, 174);
  hestimate_all->Sumw2();
  hestimate_all->GetYaxis()->SetTitle("Entries");
  
  TH1D* hestimate[bgSize];

  TH1D* hestimate_all_forRatio = new TH1D(Form("hestimate_all_forRatio"), "", 174, 0, 174);
  hestimate_all_forRatio->Sumw2();

  TH1D* hestimate_forRatio[bgSize];
  
  for(int b=0; b<bgSize; ++b){
  
    hestimate[b]= new TH1D(Form("hestimate_%d", b), "", 174, 0, 174);
    hestimate[b]->Sumw2();
    hestimate[b]->GetYaxis()->SetTitle("Entries");
    hestimate[b]->SetFillColor(colors[b]);
    hestimate[b]->SetLineColor(1);

    hestimate_forRatio[b]= new TH1D(Form("hestimate_forRatio%d", b), "", 174, 0, 174);
    hestimate_forRatio[b]->Sumw2();
    hestimate_forRatio[b]->GetYaxis()->SetTitle("Entries");
    hestimate_forRatio[b]->SetFillColor(colors[b]);
    hestimate_forRatio[b]->SetLineColor(1);
    
  }

  THStack bgStack("bgStack", "");

  TH1D* hPull = new TH1D("hPull", "", 101, -5.05, 5.05);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Est. - Obs.)/#sigma");
  hPull->GetYaxis()->SetTitle("Entries");

  TH1D* hPvalue = new TH1D("hPvalue", "", 14, 0, 1.05);
  hPvalue->Sumw2();
  hPvalue->GetXaxis()->SetTitle("p-value");
  hPvalue->GetYaxis()->SetTitle("Entries");

  // int Nobs08[100];
  // for (int t=0; t<100; ++t)
  //   Nobs08[t]=0;
  
  std::string fullPath = outputdir;

  std::string labelsMono[12]={"[200,250]","[250,350]","[350,450]","[450,575]","[575,700]","[700,1000]",">1000", "[200,250]","[250,350]","[350,450]","[450,575]",">575"};

  TFile* fmono=TFile::Open("mlfit_monojet.root");
  TFile* fvlht=TFile::Open("mlfit_veryLowHT.root");
  TFile* flht =TFile::Open("mlfit_LowHT.root");
  TFile* fmht =TFile::Open("mlfit_MediumHT.root");
  TFile* fhht =TFile::Open("mlfit_HighHT.root");
  TFile* feht =TFile::Open("mlfit_ExtremeHT.root");
  
  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      int nBins;
      double *bins;
      iMT2->getBins(nBins, bins);
      
      TH1D* h_first = data->get(*iMT2)->yield;
      TGraphAsymmErrors* g_first = MT2DrawTools::getPoissonGraph(h_first);      
      
      TH1D* h_sig[sigSize+sigContSize];
      TH3D* h_sig3d[sigSize+sigContSize];
      for(int s=0; s<sigSize; ++s)
        h_sig3d[s] = analysesSignal[s]->get(*iMT2)->yield3d;
      for(int s=0; s<sigContSize; ++s)
        h_sig3d[s+sigSize] = analysesSignalCont[s]->get(*iMT2)->yield3d;

      int nBinsMT2;
      double* binsMT2;
      iMT2->getBins(nBinsMT2, binsMT2);

      int nBinsM=81;
      double binWidthM=25.;
      double binsM[nBinsM+1];
      for (int b=0; b<=nBinsM; ++b)
        binsM[b]=b*binWidthM;

      for(int s=0; s<sigSize+sigContSize; ++s)
        if( h_sig3d[s] == 0 )
          h_sig3d[s] = new TH3D("emptyHisto", "", nBinsMT2, binsMT2, nBinsM, binsM, nBinsM, binsM);

      float m1[]={1800., 700., 1300., 700., 700., 400., 1000., 600., 500., 1200., 700., 650., 600., 200.};
      float m2[]={100., 600., 100., 600., 100., 200., 100., 100., 300., 100., 400., 100., 200., 100.};

      int binY, binZ;

      for(int s = 0; s<sigSize; ++s){

        binY = h_sig3d[s]->GetYaxis()->FindBin(m1[s]);
        binZ = h_sig3d[s]->GetZaxis()->FindBin(m2[s]);
        h_sig[s] = h_sig3d[s]->ProjectionX(Form("mt2_%d", s), binY, binY, binZ, binZ);

      }

      for(int s = 0; s<sigContSize; ++s){

	binY = h_sig3d[sigSize+s]->GetYaxis()->FindBin(m1[sigSize+s]);
        binZ = h_sig3d[sigSize+s]->GetZaxis()->FindBin(m2[sigSize+s]);
        h_sig[sigSize+s] = h_sig3d[sigSize+s]->ProjectionX(Form("mt2_%d", sigSize+s), binY, binY, binZ, binZ);

	TH3D* this_signal3d_crsl = (TH3D*) analysesSignalCont[s]->get(*iMT2)->yield3d_crsl->Clone();
	
	if( this_signal3d_crsl == 0 ){
	  continue;
	}
	else{
	  
	  TH1D* this_signal_crsl   = this_signal3d_crsl->ProjectionX("mt2_crsl", binY, binY, binZ, binZ);
	  
	  TH1D* this_signal_alpha  = (TH1D*) analysesSignalCont[s]->get(*iMT2)->yield_alpha->Clone();
	  
	  TH1D* this_signalContamination = (TH1D*) this_signal_crsl->Clone();
	  this_signalContamination->Multiply(this_signal_alpha);
	
	  h_sig[sigSize+s]->Add(this_signalContamination, -1.0);
	  //	  h_sig[sigSize+s]->Scale(2.26355/2.26);

	}
      
      }


      TFile* histoFile = TFile::Open( Form("%s/histograms_%s.root", fullPath.c_str(), iMT2->getName().c_str()), "recreate" );
      histoFile->cd();
      //      h_first->Write();
      
      g_first->SetMarkerStyle(20);
      g_first->SetMarkerSize(1.6);
      g_first->SetLineColor( 1 );
      g_first->SetMarkerColor( 1 );

      THStack bgStack_region("bgStack_region", "");
      
      TH1D* h_second_all;
      TH1D* h_second_forRatio_all;
      TH1D* h_second[bgSize];
      TH1D* h_second_forRatio[bgSize];

      for(int b=0; b< bgSize; ++b){
	
	h_second[b] = new TH1D(Form("h_second_%d", b), "", nBins, bins);
	
	h_second_forRatio[b] = new TH1D(Form("h_second_%d", b), "", nBins, bins);
	
	h_second[b]->SetFillColor( colors[b] );
	h_second[b]->SetLineColor( 1 );
	
      }

      h_second_all = new TH1D(Form("h_second_all"), "", nBins, bins);
      h_second_all->Sumw2();

      h_second_forRatio_all = new TH1D(Form("h_second_forRatio_all"), "", nBins, bins);
      h_second_forRatio_all->Sumw2();
      
      for( int iBin=0; iBin<nBins; ++iBin ) {

	std::string tableName;
	if( iBin < nBins-1 )
	  tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lfto%.0lf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin], bins[iBin+1]) );
	else
	  tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lftoInf.txt", dir.c_str(), iMT2->getName().c_str(), bins[iBin] ));

	//BGTable thisTable = getTable(tableName);
	getTable(tableName);
	

	float totalPost_llep;
	float totalPost_zinv;
	float totalPost_qcd;
	// float totalPost_Err_llep;
	// float totalPost_Err_zinv;
	// float totalPost_Err_qcd;

	float totalPost;
	float totalPost_Err;

	if(iRegion <=12){

	  int ch=iRegion+iBin;
	  fmono->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
	  
	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;

	  // totalPost_Err_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinError(1) : 0;
	  // totalPost_Err_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinError(1) : 0;
	  // totalPost_Err_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd"))  ? thisqcd ->GetBinError(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);

	  gDirectory->cd("..");
	      
	}
	else if(iRegion >=13 && iRegion <= 36){

	  int ch=iRegion-12+iBin;
	  
	  fvlht->cd();
	  gDirectory->cd("shapes_fit_b");
	    
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
	  
	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
	  gDirectory->cd("..");

	}
	else if(iRegion >=37 && iRegion <= 67){
	  
	  int ch=iRegion-36+iBin;
	  
	  flht->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
	  
	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
	  gDirectory->cd("..");
	
	}
	else if(iRegion >=68 && iRegion <= 109){

	  int ch=iRegion-67+iBin;

	  fmht->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");

	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
	  gDirectory->cd("..");
	  
	}
	else if(iRegion >=110 && iRegion <= 144){

	  int ch=iRegion-109+iBin;
	  
	  fhht->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");

	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);	
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//  	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);

	  gDirectory->cd("..");

	}
	else if(iRegion >=145){

	  int ch=iRegion-144+iBin;

	  feht->cd();
	  gDirectory->cd("shapes_fit_b");
	  
	  std::string thisCh = Form("ch%d", ch);
	  gDirectory->cd(thisCh.c_str());
	  
	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  totalPost = thisBG->GetBinContent(1);
	  totalPost_Err = thisBG->GetBinError(1);
//	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
//  	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
//	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
	  gDirectory->cd("..");

	}
	
	gDirectory->cd();
	
	std::cout << "Filling histograms" << std::endl;
	
	//QCD
	h_second[0]->SetBinContent(iBin+1, totalPost_qcd);
	h_second[0]->SetBinError(iBin+1, 0.);
	
	h_second_forRatio[0]->SetBinContent(iBin+1, totalPost_qcd);
	h_second_forRatio[0]->SetBinError(iBin+1, 0);
	
	//Lost Lepton
	h_second[1]->SetBinContent(iBin+1, totalPost_llep);
	h_second[1]->SetBinError(iBin+1, 0.);

	h_second_forRatio[1]->SetBinContent(iBin+1, totalPost_llep);
	h_second_forRatio[1]->SetBinError(iBin+1, 0.);

	//Invisible Z
	h_second[2]->SetBinContent(iBin+1, totalPost_zinv);
	h_second[2]->SetBinError(iBin+1, 0.);

	h_second_forRatio[2]->SetBinContent(iBin+1, totalPost_zinv);
	h_second_forRatio[2]->SetBinError(iBin+1, 0);

	h_second_all->SetBinContent(iBin+1, totalPost);
	h_second_all->SetBinError(iBin+1, totalPost_Err);

	h_second_forRatio_all->SetBinContent(iBin+1, totalPost);
	h_second_forRatio_all->SetBinError(iBin+1, 0);

      }	
	
      for(int b=0; b<bgSize; ++b){
      
	bgStack_region.Add(h_second[b]);
	      
      }
      
      for(int iBin=1; iBin<=nBins; ++iBin){

	float thisData    = h_first->GetBinContent(iBin);

	float thisEst     = h_second_all->GetBinContent(iBin);
	float thisEstErr  = h_second_all->GetBinError(iBin);
	
	int obs = thisData;
	float meanExp = thisEst;
	float meanUnc = thisEstErr;
	//	int N=100000;
	int N=1;

	TRandom3 randGen;

	unsigned int counterN(0);
	unsigned int counterD(0);

	//	bool doLeft = obs<=meanExp;
	
	//	TH1D* hrand = new TH1D("hrand", "", 100000, 0, 100000);

	for(int i=0; i<N; i++){
	  double meanShift = randGen.Gaus(0.,(double)meanUnc);
	  double rand = randGen.Poisson(double(meanExp)+meanShift); counterD++;

	  //	  hrand->Fill(rand);
	  
//	  if(doLeft)
//	    {
//	  if(rand<=obs) counterN++;
//	  }
//	  else{
//	    if(rand>=obs) counterN++;
//	  }

	  if(rand>=obs) counterN++;

 
	}

//	TFile* frand = new TFile("rand.root", "UPDATE");
//	frand->cd();
//	hrand->Write();
//	frand->Close();
      

	double prob=1.0*counterN/counterD;
	// double significance  = TMath::NormQuantile(1-prob);
	
//	std::cout << "probability: " << prob  << std::endl;
//	std::cout << "significance: " << significance << std::endl;
//	
//	std::cout << "Bin: " << iBin << std::endl;
//	std::cout << "Obs: " << obs << std::endl;
//	std::cout << "exp: " << meanExp << "\t" << meanUnc << std::endl;
	
	//	if(meanExp > 1.0 && obs > 0)
	//	if(obs > 0)
	hPvalue->Fill(prob);
	
//	for (int t=0; t<100; ++t){
//
//	  double estShift = randGen.Gaus(0.,(double)meanUnc);
//	  obs = randGen.Poisson(double(meanExp)+estShift);
////	  std::cout << "Obs: " << obs << std::endl;
////	  std::cout << "exp: " << hrand->GetMean() << std::endl;
//	  
//	  int tBin = hrand->GetXaxis()->FindBin(obs);
//	  prob = (hrand->Integral(tBin, -1))/(hrand->Integral(1, -1));
//	  //	  std::cout << "probability: " << prob  << std::endl;
//	  
//	  if( prob >= 0.5 ) ++Nobs08[t];
//	  
//  
//	}

      }
	
      
      double err_data;
      double int_data;
      for (int iBin=1; iBin<=nBins; ++iBin){
	
	int_data = h_first->GetBinContent(iBin);
	err_data = h_first->GetBinError(iBin);
	
	hdata->SetBinContent(iRegion, int_data);
	hdata->SetBinError(iRegion, err_data);

	hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
	
	double int_sig[sigSize+sigContSize];
	for(int s=0; s<sigSize+sigContSize; ++s){
	  
	  int_sig[s] = h_sig[s]->GetBinContent(iBin);
	  
	  if( int_sig[s] < 0. ) int_sig[s]=0.;
	  
	  hsig[s]->SetBinContent(iRegion, int_sig[s]);
	  hsig[s]->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );

	}

	for(int b=0; b<bgSize; ++b){

	  double err_int = fabs(h_second[b]->GetBinError(iBin));
	  double integral = fabs(h_second[b]->GetBinContent(iBin));
	  hestimate[b]->SetBinContent(iRegion, integral);
	  hestimate[b]->SetBinError(iRegion, err_int);

	  double integral_forRatio = fabs(h_second_forRatio[b]->GetBinContent(iBin));
	  hestimate_forRatio[b]->SetBinContent(iRegion, integral_forRatio);
	  hestimate_forRatio[b]->SetBinError(iRegion, 0);
	  
	  std::string thisLabel=Form("%s,%d", niceNames[1].c_str(), iBin); 
	  hestimate[b]->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
	  
	}

	hestimate_all->SetBinContent(iRegion, h_second_all->GetBinContent(iBin));
	hestimate_all->SetBinError(iRegion, h_second_all->GetBinError(iBin));
	
	if( iMT2->nJetsMax()==1 )
	  hestimate_all->GetXaxis()->SetBinLabel( iRegion, labelsMono[iRegion-1].c_str() );
	else{
	  std::string thisLabel;
	  if( iBin < nBins )
	    thisLabel=Form("[%.0lf,%.0lf]", bins[iBin-1], bins[iBin]);
	  else
	    thisLabel=Form(">%.0lf", bins[iBin-1]);
	  hestimate_all->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() ); 
	}

	std::string thisLabel=Form("%s,%d", niceNames[1].c_str(), iBin); 
//	hestimate_all->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
		
	hestimate_all_forRatio->SetBinContent(iRegion, h_second_all->GetBinContent(iBin));
	hestimate_all_forRatio->SetBinError(iRegion, 0);
	hestimate_all_forRatio->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
		
	std::cout <<"iRegion " << iRegion << std::endl;
	++iRegion;
       
      }

      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.15);
      pad1->Draw();
      pad1->cd();

      float xMin = h_first->GetXaxis()->GetXmin();
      float xMax = h_first->GetXaxis()->GetXmax();
      float yMax_1 = h_first->GetMaximum()*1.5;
      float yMax_2 = 1.2*(h_first->GetMaximum() + h_first->GetBinError(h_first->GetMaximumBin()));
      float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
      
      float yMax_3 = h_second_all->GetMaximum()*1.5;
      float yMax_4 = 1.2*(h_second_all->GetMaximum() + h_second_all->GetBinError(h_second_all->GetMaximumBin()));
      float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
      float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
      
      TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
      h2_axes->SetXTitle("M_{T2} [GeV]");
      h2_axes->SetYTitle("Entries");

      if(iMT2->nJetsMax()==1)
	h2_axes->SetXTitle("H_{T} [GeV]");

      h2_axes->Draw();
 
      for( unsigned i=0; i<niceNames.size(); ++i ) {

        float yMax = 0.9-(float)i*0.05;
        float yMin = yMax - 0.05;
        TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
        regionText->SetTextSize(0.035);
        regionText->SetTextFont(42);
        regionText->SetFillColor(0);
        regionText->SetTextAlign(11);
        regionText->AddText( niceNames[i].c_str() );
        regionText->Draw("same");
	
      }


      TLegend* legend = new TLegend( 0.6, 0.9-(bgSize)*0.06, 0.93, 0.9 );
      legend->SetTextSize(0.038);
      legend->SetTextFont(42);
      legend->SetFillColor(0);
      legend->AddEntry( g_first, "Data", "P" );

      legend->Draw("same");

      bgStack_region.Draw("histo, same");
      g_first->Draw("pe,same");
      
      TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
      labelTop->Draw("same");

      //      TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Supplementary");
      TPaveText* labelCMS = MT2DrawTools::getLabelCMS();
      labelCMS->Draw("same");

      gPad->RedrawAxis();

      c1->cd();
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.1);
      pad2->Draw();
      pad2->cd();

      std::string thisName = Form("%s_ratio", h_first->GetName());
      TH1D* h_ratio = (TH1D*) h_first->Clone(thisName.c_str());
      h_ratio->Divide(h_second_all);
      //      h_ratio->Write();
      h_ratio->SetStats(0);	    
      h_ratio->SetMarkerStyle(20);
      h_ratio->SetMarkerColor(1);
      h_ratio->SetLineColor(1);
      h_ratio->GetXaxis()->SetLabelSize(0.00);
      h_ratio->GetXaxis()->SetTickLength(0.09);
      h_ratio->GetYaxis()->SetNdivisions(5,5,0);
      h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
      h_ratio->GetYaxis()->SetTitleSize(0.17);
      h_ratio->GetYaxis()->SetTitleOffset(0.4);
      h_ratio->GetYaxis()->SetLabelSize(0.17);
      h_ratio->GetYaxis()->SetTitle("Ratio");
            
      TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, 0.0, 2.0 );
      h2_axes_ratio->GetYaxis()->SetTitle("Ratio");
      
      TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
      lineCentral->SetLineColor(1);

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      h_ratio->Draw("pe,same");
      
      gPad->RedrawAxis();
      
      c1->cd();

      c1->SaveAs( Form("%s/mt2_%s.pdf", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.png", fullPath.c_str(), iMT2->getName().c_str()) );
      c1->SaveAs( Form("%s/mt2_%s.eps", fullPath.c_str(), iMT2->getName().c_str()) );

      delete c1;

      delete h2_axes;
      delete h2_axes_ratio;

      delete h_ratio;

      delete h_first;

      delete h_second_all;

      for(int b=0; b<bgSize; ++b)
	delete h_second[b];
      
      //      ++iRegion;

  } // for MT2 regions


  for(int b=0; b<bgSize; ++b){

    hestimate[b]->SetLineWidth(0);
    bgStack.Add(hestimate[b]);
    //bgStack.Add(hestimate_forRatio[b]);
    
//    if(b==0) hestimate_all = (TH1D*) hestimate[b]->Clone("hestimate_all");
//    else hestimate_all->Add(hestimate[b]);
//
//    if(b==0) hestimate_all_forRatio = (TH1D*) hestimate_forRatio[b]->Clone("hestimate_all_forRatio");
//    else hestimate_all_forRatio->Add(hestimate_forRatio[b]);

  }

  
  TGraphAsymmErrors* g_Ratio = MT2DrawTools::getRatioGraph(hdata, hestimate_all_forRatio, "binWidth");
  g_Ratio->SetMarkerStyle(20);
  g_Ratio->SetMarkerSize(1.6);
  g_Ratio->SetMarkerColor( 1 );
  g_Ratio->SetLineColor(1);
  g_Ratio->GetXaxis()->SetLabelSize(0.00);
  g_Ratio->GetXaxis()->SetTickLength(0.09);
  g_Ratio->GetYaxis()->SetNdivisions(5,5,0);
  g_Ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  g_Ratio->GetYaxis()->SetTitleSize(0.17);
  g_Ratio->GetYaxis()->SetTitleOffset(0.4);
  g_Ratio->GetYaxis()->SetLabelSize(0.17);
  g_Ratio->GetYaxis()->SetTitle("Ratio");

  TGraphAsymmErrors* g_Ratio_zero = new TGraphAsymmErrors(*(g_Ratio));
  g_Ratio_zero->SetMarkerSize(0);
  g_Ratio_zero->SetLineColor( 1 );
  g_Ratio_zero->SetMarkerColor( 1 );

  TGraphAsymmErrors* gdata = MT2DrawTools::getPoissonGraph(hdata, true, "binWidth");
  gdata->GetYaxis()->SetTitle("Entries");
  gdata->SetMarkerStyle(20);
  gdata->SetMarkerSize(1.6);
  gdata->SetLineColor( 1 );
  gdata->SetMarkerColor( 1 );

  TGraphAsymmErrors* gdata_zero = new TGraphAsymmErrors(*(gdata));
  gdata_zero->SetMarkerSize(0);
  gdata_zero->SetLineColor( 1 );
  gdata_zero->SetMarkerColor( 1 );


  for(int iBin=1; iBin<=hestimate_all->GetNbinsX(); ++iBin){

    float thisData    = hdata->GetBinContent(iBin);
    float thisDataErr = gdata->GetErrorY(iBin-1);
    std::cout << "TEST!" << std::endl;
    std::cout << thisData << "\t" << gdata->GetY()[iBin-1] << "\t" << thisDataErr << std::endl;

    float thisEst     = hestimate_all->GetBinContent(iBin);
    float thisEstErr  = hestimate_all->GetBinError(iBin);


    hPull->Fill( (thisEst-thisData)/( TMath::Sqrt( thisDataErr*thisDataErr + thisEstErr*thisEstErr ) ) );

  
  }
  

  TCanvas* c2 = new TCanvas("c2", "", 1100, 600);
  c2->cd();
  

  std::string thisName = Form("%s_ratio", hdata->GetName());
  TH1D* h_Ratio = (TH1D*) hdata->Clone(thisName.c_str());
  h_Ratio->Divide(hestimate_all_forRatio);
  //  h_Ratio->Write();
  h_Ratio->SetStats(0);
  h_Ratio->SetMarkerStyle(20);
  h_Ratio->SetLineColor(1);
  h_Ratio->GetXaxis()->SetLabelSize(0.00);
  h_Ratio->GetXaxis()->SetTickLength(0.09);
  h_Ratio->GetYaxis()->SetNdivisions(5,5,0);
  h_Ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  h_Ratio->GetYaxis()->SetTitleSize(0.17);
  h_Ratio->GetYaxis()->SetTitleOffset(0.4);
  h_Ratio->GetYaxis()->SetLabelSize(0.17);
  h_Ratio->GetYaxis()->SetTitle("Ratio");

  TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  pad1->SetLogy();
  
  float yMax_1 = hdata->GetMaximum()*1.5;
  float yMax_2 = 1.2*(hdata->GetMaximum() + hdata->GetBinError(hestimate_all->GetMaximumBin()));
  float yMax1 = (yMax_1>yMax_2) ? yMax_1 : yMax_2;
  float yMax_3 = hestimate_all->GetMaximum()*1.5;
  float yMax_4 = 1.2*(hestimate_all->GetMaximum() + hestimate_all->GetBinError(hestimate_all->GetMaximumBin()));
  float yMax2 = (yMax_3>yMax_4) ? yMax_3 : yMax_4;
  float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
  
  float yMin = 1e-3;
  //  yMin=0;
  if(drawSignals)
    //    yMax*=50.;
    yMax*=100.;
  else
    yMax*=20;

  int thisBin=174;
  
  hestimate_all->GetXaxis()->SetRangeUser(0, thisBin);  
  gdata->GetXaxis()->SetRangeUser(0, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(0, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  for(int s=0; s<sigSize+sigContSize; ++s)
    hsig[s]->GetXaxis()->SetRangeUser(0, thisBin);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->GetXaxis()->SetLabelFont(62);
  hestimate_all->SetFillStyle(3244);
  hestimate_all->SetFillColor(kGray+2);

//  TGraphAsymmErrors* g_data = new TGraphAsymmErrors(0);
//  for( int iBin=1; iBin<(hdata->GetXaxis()->GetNbins()+1); ++iBin ) {
//
//    double y;
//    double x, xerr;
//
//    x = hdata->GetBinCenter(iBin);
//    xerr = hdata->GetBinWidth(iBin)/2.;
//
//    y = hdata->GetBinContent(iBin);
//    double yerr = hdata->GetBinError(iBin);
//
//    int thisPoint = g_data->GetN();
//    g_data->SetPoint( thisPoint, x, y );
//    g_data->SetPointError( thisPoint, xerr, xerr, yerr, yerr );
//
//  }
  

  if(drawSignals)
      bgStack.Add(hsig[S]);

//  hdata->Draw("pe");
//  bgStack.Draw("histo, same");
//  hestimate_all->Draw("E2,same");
//  hdata->Draw("pe,same");

  gdata->Draw("pe");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");
  
  TH1D* postfit=new TH1D("postfit", "", 1, 0, 1);
  postfit->SetLineColor(0);
  postfit->SetFillColor(0);
  

  TLegend* legend;// = new TLegend( 0.8, 0.9-(bgSize+1-1)*0.06-0.06+0.02, 0.93, 0.9-0.06+0.02 );
  if(drawSignals)
    legend = new TLegend( 0.7, 0.9-(bgSize+1-1)*0.06-0.06+0.02-0.06-0.01, 0.85, 0.9-0.06+0.02+0.02+0.04 );
  else
    legend = new TLegend( 0.8, 0.9-(bgSize+1-1)*0.06-0.06+0.02+0.02, 0.93, 0.9-0.06+0.02+0.02 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( hdata, "Data", "PL" );
  //  legend->AddEntry( postfit, "Post-fit SM", "F");
  legend->AddEntry( hestimate[0], "Multijet", "F");
  legend->AddEntry( hestimate[1], "Lost lepton", "F");
  legend->AddEntry( hestimate[2], "Z #rightarrow #nu#bar{#nu}", "F");

  legend->Draw("same");

  //  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
  labelTop->Draw("same");

  //  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Supplementary");
  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Preliminary");
  labelCMS->Draw("same");
  
  int nHTRegions = 6;
  std::vector< std::string > htRegions;
  htRegions.push_back("1 Jet");
  htRegions.push_back("H_{T} [200, 450] GeV");
  htRegions.push_back("H_{T} [450, 575] GeV");
  htRegions.push_back("H_{T} [575, 1000] GeV");
  htRegions.push_back("H_{T} [1000, 1500] GeV");
  htRegions.push_back("H_{T} > 1500 GeV");
//  htRegions.push_back("#it{H}_{T} [200, 450] GeV");
//  htRegions.push_back("#it{H}_{T} [450, 575] GeV");
//  htRegions.push_back("#it{H}_{T} [575, 1000] GeV");
//  htRegions.push_back("#it{H}_{T} [1000, 1500] GeV");
//  htRegions.push_back("#it{H}_{T} > 1500 GeV");
////  htRegions.push_back("very low H_{T}");
////  htRegions.push_back("1 Jet");
////  htRegions.push_back("low H_{T}");
////  htRegions.push_back("medium H_{T}");
////  htRegions.push_back("high H_{T}");
////  htRegions.push_back("extreme H_{T}");
  
  TPaveText* htBox[5];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){
    
    htBox[iHT] = new TPaveText(0.4, 0.9-0.06+0.02, 0.7, 0.85+0.02, "brNDC");
    htBox[iHT]->AddText( htRegions[iHT].c_str() );
    
    htBox[iHT]->SetBorderSize(0);
    htBox[iHT]->SetFillColor(kWhite);
    htBox[iHT]->SetTextSize(0.038);
    htBox[iHT]->SetTextAlign(21); // align centered
    htBox[iHT]->SetTextFont(62);
    //    htBox[iHT]->Draw("same");

  }
  //  htBox[0]->Draw("same");

  gPad->RedrawAxis();
  
  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();  

  TH2D* h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0., 3.0 );
  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.20);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.4);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitle("Ratio");
  
  TLine* LineCentral = new TLine(0, 1.0, thisBin, 1.0);
  LineCentral->SetLineColor(1);

  
  std::string thisName_Band =  Form("%s_band", hestimate_all->GetName());
  TH1D* h_band = (TH1D*)hestimate_all->Clone(thisName_Band.c_str());
  h_band->SetMarkerSize(0);
  h_band->SetFillColor (kGray+2);
  h_band->SetFillStyle (3244);
  for ( int iBin=1; iBin <= hestimate_all->GetNbinsX(); iBin++){
    
    h_band->SetBinContent(iBin,1);

    double error=0;

    if(hestimate_all_forRatio->GetBinContent(iBin)>0)
      error = hestimate_all->GetBinError(iBin)/hestimate_all->GetBinContent(iBin);
    //else error = hestimate_all->GetBinError(iBin);

    h_band->SetBinError(iBin, error);

  }


  h2_axes_ratio->Draw("");
  h_band->Draw("E2same");
  LineCentral->Draw("same");
  //h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2->cd();
  //  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.pdf", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.eps", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.png", fullPath.c_str()) );
  

  TCanvas* c2_0 = new TCanvas("c2_0", "", 1100, 600);
  c2_0->cd();
  
  TPad *pad1_0 = new TPad("pad1_0","pad1_0",0,0.3-0.1,1,1);
  pad1_0->SetBottomMargin(0.15);
  pad1_0->Draw();
  pad1_0->cd();


  int oldBin=0;
  pad1_0->SetLogy();
    
  thisBin=12;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  for(int s=0; s<sigSize+sigContSize; ++s)
    hsig[s]->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);

  hestimate_all->GetYaxis()->SetTitleOffset(0.95);
  hestimate_all->GetYaxis()->SetLabelSize(0.042);

  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");

  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");
  
  htBox[0]->Draw("same");

  gPad->RedrawAxis();
  
  c2_0->cd();
  TPad *pad2_0 = new TPad("pad2_0","pad2_0",0,0,1,0.21);
  pad2_0->SetTopMargin(0.05);
  pad2_0->SetBottomMargin(0.1);
  pad2_0->Draw();
  pad2_0->cd();

  bool doLogRatio=false;

  TH2D* h2_axes_ratio_0;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
    h_Ratio->GetYaxis()->SetRangeUser(0.1, 10.0);
    g_Ratio->GetYaxis()->SetRangeUser(0.1, 10.0);
  }
  else
    //    h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0., 2.0 );
    h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0., 1.5 );

  //  TH2D* h2_axes_ratio_0 = new TH2D("axes_ratio_0", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_0->SetStats(0);
  h2_axes_ratio_0->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_0->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_0->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_0->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_0->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_0->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_0->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_0 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_0->SetLineColor(1);

  h2_axes_ratio_0->Draw("");
  h_band->Draw("E2same");
  LineCentral_0->Draw("same");
  //h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_0->cd();
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.pdf", fullPath.c_str()) );
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate.png", fullPath.c_str()) );

  c2_0->Clear();

  c2_0->cd();

  gPad->SetLogy();

  hestimate_all->Draw("");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[0]->Draw("same");

  if( drawSignals ){
    if(S<sigSize)
      legend->AddEntry( hsig[S], analysesSignal[S]->getName().c_str(), "F" );
    else if(S==sigSize)
      legend->AddEntry( hsig[S], "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t}#tilde{#chi}_{1}^{0}", "F" );
    else if(S==sigSize+1)
      legend->AddEntry( hsig[S], "T1tttt 700, 400", "F" );
    else
      legend->AddEntry( hsig[S], analysesSignalCont[S-sigSize]->getName().c_str(), "F" );

    if(sigName[S]=="T1bbbb_1500_100"){
      legend->AddEntry( postfit, "m_{#tilde{g}}=1800 GeV", "F");
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}=100 GeV", "F");
    }
    else if(sigName[S]=="T1qqqq_1300_100"){
      legend->AddEntry( postfit, "m_{#tilde{g}}=1300 GeV", "F");
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}=100 GeV", "F");
    }
    else if(sigName[S]=="T1tttt_1200_100"){
      legend->AddEntry( postfit, "m_{#tilde{g}}=1200 GeV", "F");
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}=100 GeV", "F");
    }
    else if(sigName[S]=="T2qq_1000_100"){
      legend->AddEntry( postfit, "m_{#tilde{q}}=1000 GeV", "F");
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}=100 GeV", "F");
    }
    else if(sigName[S]=="T2bb_700_100"){
      legend->AddEntry( postfit, "m_{#tilde{b}}=700 GeV", "F");
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}=100 GeV", "F");
    }
    else if(sigName[S]=="T2tt_650_100"){
      legend->AddEntry( postfit, "m_{#tilde{t}}=650 GeV", "F");
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}=100 GeV", "F");
    }
  }
  //  legend->AddEntry( postfit, "(1500, 100)", "F");
  
  legend->Draw("same");

  gPad->RedrawAxis();

  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate_plus%s.pdf", fullPath.c_str(), sigName[S].c_str()) );
  c2_0->SaveAs( Form("%s/mt2_monojet_fullEstimate_plus%s.png", fullPath.c_str(), sigName[S].c_str()) );



  TCanvas* c2_1 = new TCanvas("c2_1", "", 1100, 600);
  c2_1->cd();
  
  TPad *pad1_1 = new TPad("pad1_1","pad1_1",0,0.3-0.1,1,1);
  pad1_1->SetBottomMargin(0.15);
  pad1_1->Draw();
  pad1_1->cd();

  pad1_1->SetLogy();
    
  oldBin=thisBin;
  thisBin=36;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  for(int s=0; s<sigSize+sigContSize; ++s)
    hsig[s]->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");  
  //  hdata->Draw("pe,same");
  gdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");

  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");
  
  htBox[1]->Draw("same");

  gPad->RedrawAxis();
  
  c2_1->cd();
  TPad *pad2_1 = new TPad("pad2_1","pad2_1",0,0,1,0.21);
  pad2_1->SetTopMargin(0.05);
  pad2_1->SetBottomMargin(0.1);
  pad2_1->Draw();
  pad2_1->cd();

  TH2D* h2_axes_ratio_1;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    //  h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0., 4.0 );
    h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0., 2.5 );
  //h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  //TH2D* h2_axes_ratio_1 = new TH2D("axes_ratio_1", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_1->SetStats(0);
  h2_axes_ratio_1->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_1->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_1->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_1->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_1->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_1->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_1->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_1 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_1->SetLineColor(1);

  h2_axes_ratio_1->Draw("");
  h_band->Draw("E2same");
  LineCentral_1->Draw("same");
  //  h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");

  gPad->RedrawAxis();

  c2_1->cd();
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate.png", fullPath.c_str()) );


  c2_1->Clear();

  c2_1->cd();

  gPad->SetLogy();

  hestimate_all->Draw("");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");

  labelCMS->Draw("same");

  htBox[1]->Draw("same");

  legend->Draw("same");

  gPad->RedrawAxis();

  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate_plus%s.pdf", fullPath.c_str(), sigName[S].c_str()) );
  c2_1->SaveAs( Form("%s/mt2_veryLowHT_fullEstimate_plus%s.png", fullPath.c_str(), sigName[S].c_str()) );


  TCanvas* c2_2 = new TCanvas("c2_2", "", 1100, 600);
  c2_2->cd();
  
  TPad *pad1_2 = new TPad("pad1_2","pad1_2",0,0.3-0.1,1,1);
  pad1_2->SetBottomMargin(0.15);
  pad1_2->Draw();
  pad1_2->cd();

  pad1_2->SetLogy();
    
  oldBin=thisBin;
  thisBin=67;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  for(int s=0; s<sigSize+sigContSize; ++s)
    hsig[s]->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");  
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");

  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");
  
  htBox[2]->Draw("same");

  gPad->RedrawAxis();

  c2_2->cd();
  TPad *pad2_2 = new TPad("pad2_2","pad2_2",0,0,1,0.21);
  pad2_2->SetTopMargin(0.05);
  pad2_2->SetBottomMargin(0.1);
  pad2_2->Draw();
  pad2_2->cd();

  TH2D* h2_axes_ratio_2;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    //  h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0., 4.0 );
    h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0., 2.5 );
  //    h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  //  TH2D* h2_axes_ratio_2 = new TH2D("axes_ratio_2", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_2->SetStats(0);
  h2_axes_ratio_2->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_2->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_2->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_2->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_2->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_2->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_2->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_2 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_2->SetLineColor(1);

  h2_axes_ratio_2->Draw("");
  h_band->Draw("E2same");
  LineCentral_2->Draw("same");
  //  h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_2->cd();
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate.png", fullPath.c_str()) );

  c2_2->Clear();

  c2_2->cd();

  gPad->SetLogy();

  hestimate_all->Draw("");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[2]->Draw("same");

  legend->Draw("same");

  gPad->RedrawAxis();

  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate_plus%s.pdf", fullPath.c_str(), sigName[S].c_str()) );
  c2_2->SaveAs( Form("%s/mt2_lowHT_fullEstimate_plus%s.png", fullPath.c_str(), sigName[S].c_str()) );


  TCanvas* c2_3 = new TCanvas("c2_3", "", 1100, 600);
  c2_3->cd();
  
  TPad *pad1_3 = new TPad("pad1_3","pad1_3",0,0.3-0.1,1,1);
  pad1_3->SetBottomMargin(0.15);
  pad1_3->Draw();
  pad1_3->cd();

  pad1_3->SetLogy();
    
  oldBin=thisBin;
  thisBin=109;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  for(int s=0; s<sigSize+sigContSize; ++s)
    hsig[s]->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  // hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[3]->Draw("same");

  gPad->RedrawAxis();
  
  c2_3->cd();
  TPad *pad2_3 = new TPad("pad2_3","pad2_3",0,0,1,0.21);
  pad2_3->SetTopMargin(0.05);
  pad2_3->SetBottomMargin(0.1);
  pad2_3->Draw();
  pad2_3->cd();

  TH2D* h2_axes_ratio_3;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    //    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0., 2.5 );
    h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0., 2.0 );
  //h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  //  TH2D* h2_axes_ratio_3 = new TH2D("axes_ratio_3", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_3->SetStats(0);
  h2_axes_ratio_3->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_3->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_3->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_3->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_3->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_3->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_3->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_3 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_3->SetLineColor(1);

  h2_axes_ratio_3->Draw("");
  h_band->Draw("E2same");
  LineCentral_3->Draw("same");
  //  h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_3->cd();
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate.png", fullPath.c_str()) );

  c2_3->Clear();

  c2_3->cd();

  gPad->SetLogy();

  hestimate_all->Draw("");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[3]->Draw("same");

  legend->Draw("same");

  gPad->RedrawAxis();

  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate_plus%s.pdf", fullPath.c_str(), sigName[S].c_str()) );
  c2_3->SaveAs( Form("%s/mt2_mediumHT_fullEstimate_plus%s.png", fullPath.c_str(), sigName[S].c_str()) );




  TCanvas* c2_4 = new TCanvas("c2_4", "", 1100, 600);
  c2_4->cd();
  
  TPad *pad1_4 = new TPad("pad1_4","pad1_4",0,0.3-0.1,1,1);
  pad1_4->SetBottomMargin(0.15);
  pad1_4->Draw();
  pad1_4->cd();

  pad1_4->SetLogy();
    
  oldBin=thisBin;
  thisBin=144;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  for(int s=0; s<sigSize+sigContSize; ++s)
    hsig[s]->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[4]->Draw("same");

  gPad->RedrawAxis();
  
  c2_4->cd();
  TPad *pad2_4 = new TPad("pad2_4","pad2_4",0,0,1,0.21);
  pad2_4->SetTopMargin(0.05);
  pad2_4->SetBottomMargin(0.1);
  pad2_4->Draw();
  pad2_4->cd();

  TH2D* h2_axes_ratio_4;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    //  h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0., 2.5 );
    h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0., 2.0 );
  //h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  //  TH2D* h2_axes_ratio_4 = new TH2D("axes_ratio_4", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_4->SetStats(0);
  h2_axes_ratio_4->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_4->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_4->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_4->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_4->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_4->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_4->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_4 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_4->SetLineColor(1);

  h2_axes_ratio_4->Draw("");
  h_band->Draw("E2same");
  LineCentral_4->Draw("same");
  //  h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_4->cd();
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate.png", fullPath.c_str()) );

  c2_4->Clear();

  c2_4->cd();

  gPad->SetLogy();

  hestimate_all->Draw("");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  //hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[4]->Draw("same");

  legend->Draw("same");

  gPad->RedrawAxis();

  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate_plus%s.pdf", fullPath.c_str(), sigName[S].c_str()) );
  c2_4->SaveAs( Form("%s/mt2_highHT_fullEstimate_plus%s.png", fullPath.c_str(), sigName[S].c_str()) );




  TCanvas* c2_5 = new TCanvas("c2_5", "", 1100, 600);
  //  TCanvas* c2_5 = new TCanvas("c2_5", "", 1300, 800);
  c2_5->cd();
  
  TPad *pad1_5 = new TPad("pad1_5","pad1_5",0,0.3-0.1,1,1);
  pad1_5->SetBottomMargin(0.15);
  pad1_5->Draw();
  pad1_5->cd();

  pad1_5->SetLogy();
    
  oldBin=thisBin;
  thisBin=174;
  hestimate_all->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata->GetXaxis()->SetRangeUser(oldBin, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(oldBin, thisBin);
  for(int s=0; s<sigSize+sigContSize; ++s)
    hsig[s]->GetXaxis()->SetRangeUser(oldBin, thisBin);
  h_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  g_Ratio->GetXaxis()->SetRangeUser(oldBin, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
  hestimate_all->GetXaxis()->LabelsOption("v");
  hestimate_all->GetXaxis()->SetLabelSize(0.042);
  hestimate_all->Draw("");  

  bgStack.Draw("histo,same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");
  
  legend->Draw("same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[5]->Draw("same");

  gPad->RedrawAxis();
  
  c2_5->cd();
  TPad *pad2_5 = new TPad("pad2_5","pad2_5",0,0,1,0.21);
  pad2_5->SetTopMargin(0.05);
  pad2_5->SetBottomMargin(0.1);
  pad2_5->Draw();
  pad2_5->cd();

  TH2D* h2_axes_ratio_5;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0.1, 10.0 );
  }
  else
    //    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0., 5.0 );
    h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0., 3.25 );
    //h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0., 2.0 );

  //  TH2D* h2_axes_ratio_5 = new TH2D("axes_ratio_5", "", 10, oldBin, thisBin, 10, 0., 3.0 );
  h2_axes_ratio_5->SetStats(0);
  h2_axes_ratio_5->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio_5->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio_5->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio_5->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio_5->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio_5->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio_5->GetYaxis()->SetTitle("Data/Est.");
  
  TLine* LineCentral_5 = new TLine(oldBin, 1.0, thisBin, 1.0);
  LineCentral_5->SetLineColor(1);

  h2_axes_ratio_5->Draw("");
  h_band->Draw("E2same");
  LineCentral_5->Draw("same");
  //  h_Ratio->Draw("pe,same");
  g_Ratio->Draw("pe,same");
  
  gPad->RedrawAxis();

  c2_5->cd();
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.pdf", fullPath.c_str()) );
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.png", fullPath.c_str()) );
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate.C", fullPath.c_str()) );

  c2_5->Clear();

  c2_5->cd();

  gPad->SetLogy();

  hestimate_all->Draw("");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");

  labelTop->Draw("same");
  labelCMS->Draw("same");

  htBox[5]->Draw("same");

  legend->Draw("same");

  gPad->RedrawAxis();

  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate_plus%s.pdf", fullPath.c_str(), sigName[S].c_str()) );
  c2_5->SaveAs( Form("%s/mt2_extremeHT_fullEstimate_plus%s.png", fullPath.c_str(), sigName[S].c_str()) );




  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);
  
  TF1* fgauss= new TF1("fgauss", "gaus", -5, 5);
  fgauss->SetLineColor(2);

  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  hPull->SetStats(1110);
  hPull->Draw("hist");
  hPull->Fit("fgauss");
  fgauss->Draw("same");
  c3->SaveAs( Form("%s/PullDistribution.pdf", fullPath.c_str()) );
  c3->SaveAs( Form("%s/PullDistribution.png", fullPath.c_str()) );

  TCanvas* c4 = new TCanvas("c4", "", 600, 600);
  c4->cd();
  hPvalue->SetStats(1110);
  hPvalue->Draw("hist");
  c4->SaveAs( Form("%s/PvalueDistribution.pdf", fullPath.c_str()) );
  c4->SaveAs( Form("%s/PvalueDistribution.png", fullPath.c_str()) );

//  TH1D* hNobs08 = new TH1D("hNobs08", "", 174, 0, 174);
//  for(int t=0; t<100; ++t)
//    hNobs08->Fill(Nobs08[t]);
//  
//  TLine* l08 = new TLine(106, 0, 106, 50);
//  l08->SetLineColor(2);
//
//  int NobsBin = hNobs08->GetXaxis()->FindBin(106);
//  std::cout << "Fraction of >0.84 = " << hNobs08->Integral(NobsBin, -1) << "/" << hNobs08->Integral(1, -1) << ": " << hNobs08->Integral(NobsBin, -1)/hNobs08->Integral(1, -1) << std::endl;
//
//  TCanvas* c5 = new TCanvas("c5", "", 600, 600);
//  c5->cd();
//  hNobs08->SetStats(1110);
//  hNobs08->Draw("hist");
//  l08->Draw("same");
//  c5->SaveAs( Form("%s/Nobs08.pdf", fullPath.c_str()) );
//  c5->SaveAs( Form("%s/Nobs08.png", fullPath.c_str()) );
    
}



BGTable getTable( const std::string& tableFileName ) {

  std::ifstream ifs( tableFileName.c_str() );

  BGTable table;

  while( ifs.good() ) {


    char thisLine[256];
    ifs.getline( thisLine, 256 );
    if( thisLine[0]=='#' ) continue;

    std::istringstream thisLine_iss(thisLine);

    std::string name;
    float yield, statUp, statDn, systUp, systDn;
    thisLine_iss >> name >> yield >> statUp >> statDn >> systUp >> systDn;

    if( name=="zinv" ) {
      table.zinv = yield;
      table.zinv_statUp = statUp;
      table.zinv_statDn = statDn;
      table.zinv_systUp = systUp;
      table.zinv_systDn = systDn;
    } else if( name=="llep" ) {
      table.llep = yield;
      table.llep_statUp = statUp;
      table.llep_statDn = statDn;
      table.llep_systUp = systUp;
      table.llep_systDn = systDn;
    } else if( name=="qcd" ) {
      table.qcd = yield;
      table.qcd_statUp = statUp;
      table.qcd_statDn = statDn;
      table.qcd_systUp = systUp;
      table.qcd_systDn = systDn;
    } else {
      continue;
    }

  }

  return table;

}
