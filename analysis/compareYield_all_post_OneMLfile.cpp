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
#include "TH3D.h"
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


bool drawSignals=true;
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

  //lumi = 18.1;
  lumi = cfg.lumi();
  
  TH1::AddDirectory(kTRUE);
  
  std::string dir = cfg.getEventYieldDir();
  //  std::string outputdir = cfg.getEventYieldDir() + "/YieldComparison_dataMC_post";
  std::string outputdir = cfg.getEventYieldDir() + "/YieldComparison_dataMC_post_oneMLfile";
  
 
  MT2Analysis<MT2Estimate>* analysis = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "data" ); // any one is good, just need to know the regions                                                                    


  //  std::string sigPath="./signalScansFromDominick";
  std::string sigPath="~mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/signalScansFromDominick/";


  std::vector < MT2Analysis<MT2Estimate>* > analysesSignal;
  std::vector < MT2Analysis<MT2EstimateSigContSyst>* > analysesSignalCont;
  
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1bbbb_eth.root", "T1bbbb") );
//  //  analysesSignal[0]->setName("T1bbbb 1500,100");
//  analysesSignal[0]->setName("pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow b#chi_{1}^{0}");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1bbbb_eth.root", "T1bbbb") );
//  //  analysesSignal[1]->setName("T1bbbb 700,600");
//  analysesSignal[1]->setName("pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow b#bar{b}#chi_{1}^{0}");
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
  
   analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1bbbb_eth.root", "T1bbbb") );
   analysesSignal[1]->setName("pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow b#bar{b}#tilde{#chi}_{1}^{0}");


   analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1qqqq_eth.root", "T1qqqq") );
   analysesSignal[2]->setName("pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q}#tilde{#chi}_{1}^{0}");

   analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T1qqqq_eth.root", "T1qqqq") );
   analysesSignal[3]->setName("pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q}#tilde{#chi}_{1}^{0}");

   analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2bb_eth.root", "T2bb") );
   analysesSignal[4]->setName("pp #rightarrow #tilde{b}#bar{#tilde{b}}, #tilde{b} #rightarrow b#tilde{#chi}_{1}^{0}");

   analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2bb_eth.root", "T2bb") );
   analysesSignal[5]->setName("T2bb 400, 200");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2qq_eth.root", "T2qq") );
  analysesSignal[6]->setName("pp #rightarrow #tilde{q}#bar{#tilde{q}}, #tilde{q} #rightarrow q#tilde{#chi}_{1}^{0}");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2qq_eth.root", "T2qq") );
  analysesSignal[7]->setName("T2qq 600, 0");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( sigPath + "/T2qq_eth.root", "T2qq") );
  analysesSignal[8]->setName("T2qq 500, 300");


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
  sigName.push_back("T1bbbb_2000_100");
  sigName.push_back("T1bbbb_700_600");
  sigName.push_back("T1qqqq_1300_100");
  sigName.push_back("T1qqqq_700_600");
  sigName.push_back("T2bb_700_0");
  sigName.push_back("T2bb_400_200");
//  sigName.push_back("T2qq_1000_0");
//  sigName.push_back("T2qq_600_0");
//  sigName.push_back("T2qq_500_300");
  sigName.push_back("T1tttt_1200_100");
  sigName.push_back("T1tttt_700_400");
  sigName.push_back("T2tt_650_0");
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
  styleSig.push_back(1);
//  styleSig.push_back(1);
//  styleSig.push_back(1);
//  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);
  styleSig.push_back(1);

  int bgSize = 3;
//  int sigSize = 0;//9;
//  int sigContSize = 0;

   int sigSize = 9;
   int sigContSize = 5;

  int S=1;

  std::set<MT2Region> MT2Regions = data->getRegions();
  
  TH1D* hdata = new TH1D("hdata", "", 63, 0, 63);
  hdata->Sumw2();
  hdata->GetYaxis()->SetTitle("Entries");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineColor( 1 );
  hdata->SetMarkerColor( 1 );

  TH1D* hsig[sigSize+sigContSize];
  for(int s=0; s<sigSize+sigContSize; ++s){
    std::string thisNameS( Form("hsig_%d", s) );
    hsig[s] = new TH1D(thisNameS.c_str(), "", 63, 0, 63);
    hsig[s]->Sumw2();
    hsig[s]->GetYaxis()->SetTitle("Entries");
    hsig[s]->SetLineColor( colorsSig[s] );
    hsig[s]->SetFillColor( colorsSig[s] );
    hsig[s]->SetLineStyle( styleSig[s] );
  }
  
  TH1D* hestimate_all = new TH1D(Form("hestimate_all"), "", 63, 0, 63);
  hestimate_all->Sumw2();
  hestimate_all->GetYaxis()->SetTitle("Entries");

  TH1D* hestimate[bgSize];

  TH1D* hestimate_all_forRatio = new TH1D(Form("hestimate_all_forRatio"), "", 63, 0, 63);
  hestimate_all_forRatio->Sumw2();

  TH1D* hestimate_forRatio[bgSize];
  
  for(int b=0; b<bgSize; ++b){
  
    hestimate[b]= new TH1D(Form("hestimate_%d", b), "", 63, 0, 63);
    hestimate[b]->Sumw2();
    hestimate[b]->GetYaxis()->SetTitle("Entries");
    hestimate[b]->SetFillColor(colors[b]);
    hestimate[b]->SetLineColor(1);

    hestimate_forRatio[b]= new TH1D(Form("hestimate_forRatio%d", b), "", 63, 0, 63);
    hestimate_forRatio[b]->Sumw2();
    hestimate_forRatio[b]->GetYaxis()->SetTitle("Entries");
    hestimate_forRatio[b]->SetFillColor(colors[b]);
    hestimate_forRatio[b]->SetLineColor(1);
    
  }

  THStack bgStack("bgStack", "");

  TH1D* hPull = new TH1D("hPull", "", 21, -5.25, 5.25);
  hPull->Sumw2();
  hPull->GetXaxis()->SetTitle("(Est. - Obs.)/#sigma");
  hPull->GetYaxis()->SetTitle("Entries");
  
  TH1D* hPvalue = new TH1D("hPvalue", "", 14, 0, 1.05);
  hPvalue->Sumw2();
  hPvalue->GetXaxis()->SetTitle("p-value");
  hPvalue->GetYaxis()->SetTitle("Entries");

  std::string fullPath = outputdir;

  std::string labelsMono[12]={"[250,350]","[350,450]","[450,575]","[575,700]","[700,1000]","[1000,1200]", ">1200","[250,350]","[350,450]","[450,575]","[575,700]", ">700"};


  TFile* fall=TFile::Open( Form( "%s/mlfit_all.root", dir.c_str() ) );

  // TFile* fmono=TFile::Open( Form( "%s/mlfit_monojetHT.root", dir.c_str() ) );
  // TFile* fvlht=TFile::Open( Form( "%s/mlfit_veryLowHT.root", dir.c_str() ) );
  // TFile* flht =TFile::Open( Form( "%s/mlfit_lowHT.root", dir.c_str() ) );
  // TFile* fmht =TFile::Open( Form( "%s/mlfit_mediumHT.root", dir.c_str() ) );
  // TFile* fhht =TFile::Open( Form( "%s/mlfit_highHT.root", dir.c_str() ) );
  // TFile* feht =TFile::Open( Form( "%s/mlfit_extremeHT.root", dir.c_str() ) );
  
  int iBinRegion = 1;
  int iRegion = 1;
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

      std::vector<std::string> niceNames = iMT2->getNiceNames();
      
      int nBins;
      double *bins;
      iMT2->getBins(nBins, bins);
      
      TH1D* h_first_forExtreme = data->get(*iMT2)->yield;
      TH1D* h_first;
      if( iMT2->htMin()==1500 && iMT2->nJetsMin()>1 ){
	double *binsExtreme = bins++;
	nBins--;
	h_first = new TH1D("h_first", "", nBins, binsExtreme);
	
	for( int iBin=0; iBin<=nBins; ++iBin )
	  h_first->SetBinContent( iBin, h_first_forExtreme->GetBinContent(iBin+1) );

      }else 
	h_first = (TH1D*)h_first_forExtreme->Clone("h_first");

      //TH1D* h_first = data->get(*iMT2)->yield;
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
      

      float m1[]={2000., 1000., 1300., 700., 700., 400., 1000., 600., 500., 1200., 700., 650., 600., 200.};
      float m2[]={100., 800., 100., 600., 0., 200., 0., 0., 300., 100., 400., 0., 200., 100.};

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
          h_sig[sigSize+s]->Scale(2.26355/2.26);

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

	fall->cd();
	gDirectory->cd("shapes_fit_b");

	//	if(iRegion <=12){

	  int ch=iBinRegion;//+iBin;
	  // fmono->cd();
	  // gDirectory->cd("shapes_fit_b");
	  
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
	      
// 	}
// 	else if(iRegion >=13 && iRegion <= (12+7)){
	 
// 	  int ch=iBinRegion-12;//+iBin;

// 	  fvlht->cd();
// 	  gDirectory->cd("shapes_fit_b");
	  
// 	  std::string thisCh = Form("ch%d", ch);
// 	  gDirectory->cd(thisCh.c_str());
	  
// 	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
// 	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
// 	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
// 	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
// 	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
// 	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
// 	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;

// 	  totalPost = thisBG->GetBinContent(1);
// 	  totalPost_Err = thisBG->GetBinError(1);

// //	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
// //	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
// //	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
// 	  gDirectory->cd("..");

// 	}
// 	else if(iRegion >= (12+7+1) && iRegion <= (12+7+11) ){
	  
// 	  int ch=iBinRegion-(12+21);//+iBin;
// 	  std::cout << "lowHT " << iBinRegion << "\t" << ch << "\t" << iRegion << std::endl;	 

// 	  flht->cd();
// 	  gDirectory->cd("shapes_fit_b");
	  
// 	  std::string thisCh = Form("ch%d", ch);
// 	  gDirectory->cd(thisCh.c_str());
	  
// 	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
// 	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
// 	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
// 	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
	  
// 	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
// 	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
// 	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
	  
// 	  totalPost = thisBG->GetBinContent(1);
// 	  totalPost_Err = thisBG->GetBinError(1);
// //	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
// //	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
// //	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
// 	  gDirectory->cd("..");
	
// 	}
// 	else if(iRegion >=(12+7+11+1) && iRegion <= (12+7+11+11) ){

// 	  int ch=iBinRegion-(12+21+40);//+iBin;
// 	  std::cout << "mediumHT " << iBinRegion << "\t" << ch << "\t" << iRegion << std::endl;	 

// 	  fmht->cd();
// 	  gDirectory->cd("shapes_fit_b");
	  
// 	  std::string thisCh = Form("ch%d", ch);
// 	  gDirectory->cd(thisCh.c_str());
	  
// 	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
// 	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
// 	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
// 	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");

// 	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
// 	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
// 	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
// 	  totalPost = thisBG->GetBinContent(1);
// 	  totalPost_Err = thisBG->GetBinError(1);
// //	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
// //	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
// //	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
// 	  gDirectory->cd("..");
	  
// 	}
// 	else if(iRegion >=(12+7+11+11+1 ) && iRegion <= (12+7+11+11+11)){

// 	  int ch=iBinRegion-(12+21+40+51);//+iBin;
// 	  std::cout << iBinRegion << "\t" << ch <<  "\t" << iRegion << std::endl;	  

// 	  fhht->cd();
// 	  gDirectory->cd("shapes_fit_b");
	  
// 	  std::string thisCh = Form("ch%d", ch);
// 	  gDirectory->cd(thisCh.c_str());
	  
// 	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
// 	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
// 	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
// 	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");

// 	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
// 	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
// 	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
// 	  totalPost = thisBG->GetBinContent(1);
// 	  totalPost_Err = thisBG->GetBinError(1);	
// //	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
// //  	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
// //	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
// 	  gDirectory->cd("..");

// 	}
// 	else if(iRegion >= (12+7+11+11+11+1) ){

// 	  int ch=iBinRegion-(12+21+40+51+53);//+iBin;
// 	  std::cout << iBinRegion << "\t" << ch <<  "\t" << iRegion << std::endl;	  

// 	  feht->cd();
// 	  gDirectory->cd("shapes_fit_b");
	  
// 	  std::string thisCh = Form("ch%d", ch);
// 	  gDirectory->cd(thisCh.c_str());
	  
// 	  TH1F* thisBG=(TH1F*)gDirectory->Get("total_background");
// 	  TH1F* thisllep=(TH1F*)gDirectory->Get("llep");
// 	  TH1F* thiszinv=(TH1F*)gDirectory->Get("zinv");
// 	  TH1F* thisqcd=(TH1F*)gDirectory->Get("qcd");
	  
// 	  totalPost_llep = (gDirectory->GetListOfKeys()->Contains("llep")) ? thisllep->GetBinContent(1) : 0;
// 	  totalPost_zinv = (gDirectory->GetListOfKeys()->Contains("zinv")) ? thiszinv->GetBinContent(1) : 0;
// 	  totalPost_qcd  = (gDirectory->GetListOfKeys()->Contains("qcd")) ? thisqcd ->GetBinContent(1) : 0;
	  
// 	  totalPost = thisBG->GetBinContent(1);
// 	  totalPost_Err = thisBG->GetBinError(1);
// //	  totalPost = totalPost_llep+totalPost_zinv+totalPost_qcd;
// //  	  totalPost_Err = TMath::Sqrt(totalPost_Err_llep*totalPost_Err_llep + totalPost_Err_zinv*totalPost_Err_zinv + totalPost_Err_qcd*totalPost_Err_qcd);
// //	  totalPost_Err = (totalPost_Err > thisBG->GetBinError(1)) ? totalPost_Err : thisBG->GetBinError(1);
	  
// 	  gDirectory->cd("..");

// 	}
	
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
	
	++iBinRegion;

      }
      
      
      for(int b=0; b<bgSize; ++b){
	
	bgStack_region.Add(h_second[b]);
	      
	double error;
	hestimate_all->SetBinContent(iRegion, h_second_all->IntegralAndError(1, nBins+1, error));
        hestimate_all->SetBinError(iRegion, error);

	std::string thisLabel=Form("%s", niceNames[1].c_str());
	if( iMT2->nJetsMax()==1 )
	  hestimate_all->GetXaxis()->SetBinLabel( iRegion, labelsMono[iRegion-1].c_str() );
	else
	  hestimate_all->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
	
	//hestimate_all->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );

        hestimate_all_forRatio->SetBinContent(iRegion, h_second_all->Integral());
        hestimate_all_forRatio->SetBinError(iRegion, 0);
	hestimate_all_forRatio->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
	
      }
      
//      for(int iBin=1; iBin<=nBins; ++iBin){
//
//	float thisData    = h_first->GetBinContent(iBin);
//	float thisDataErr = h_first->GetBinError(iBin);
//
//	float thisEst     = h_second_all->GetBinContent(iBin);
//	float thisEstErr  = h_second_all->GetBinError(iBin);
//	
//	
//	hPull->Fill( (thisEst-thisData)/( TMath::Sqrt( thisDataErr*thisDataErr + thisEstErr*thisEstErr ) ) );
//	
//      }

      
      int firstBin=1;
      double err_data;
      double int_data;

      if( iMT2->htMin()==1200 && iMT2->nJetsMin()==1 ){ 
	std::cout << "happening                here " << std::endl;
	int_data = h_first->IntegralAndError(firstBin, nBins, err_data);
      }else
	int_data = h_first->IntegralAndError(firstBin, nBins+1, err_data);
      
      hdata->SetBinContent(iRegion, int_data);
      hdata->SetBinError(iRegion, err_data);

      std::cout<<"Filled data  with: " << int_data << std::endl;
      
      hdata->GetXaxis()->SetBinLabel( iRegion, niceNames[1].c_str() );
      
      
      double int_sig[sigSize+sigContSize];
      for(int s=0; s<sigSize+sigContSize; ++s){

	int_sig[s]=0;

	for(int b=1; b<=h_sig[s]->GetXaxis()->GetNbins(); b++){

	  if(iMT2->htMin()>=1500 && b==1) continue;

	  if(h_sig[s]->GetBinContent(b) >= 0)
	    int_sig[s]+=h_sig[s]->GetBinContent(b);
	}

	hsig[s]->SetBinContent(iRegion, int_sig[s]);
	
	std::string thisLabel=Form("%s", niceNames[1].c_str());
	hsig[s]->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );

      }

      for(int b=0; b<bgSize; ++b){
	
	double err_int;
	double integral = h_second[b]->IntegralAndError(firstBin, nBins+1, err_int);
	hestimate[b]->SetBinContent(iRegion, integral);
	hestimate[b]->SetBinError(iRegion, err_int);
	
	double err_int_forRatio;
	double integral_forRatio = h_second_forRatio[b]->IntegralAndError(firstBin, nBins+1, err_int_forRatio);;
	hestimate_forRatio[b]->SetBinContent(iRegion, integral_forRatio);
	hestimate_forRatio[b]->SetBinError(iRegion, 0);
	
	std::cout<<"Filled estimate for background " << b << " with: " << integral << std::endl;
	
	std::string thisLabel=Form("%s", niceNames[1].c_str());
	hestimate[b]->GetXaxis()->SetBinLabel( iRegion, thisLabel.c_str() );
	
      }
      
      
      TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
      c1->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
      pad1->SetBottomMargin(0.18);
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
      
      ++iRegion;

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

  for(int iBin=1; iBin<=hestimate_all->GetNbinsX(); ++iBin){
    
    float thisData    = hdata->GetBinContent(iBin);
    
    float thisEst     = hestimate_all->GetBinContent(iBin);
    float thisEstErr  = hestimate_all->GetBinError(iBin);

    int obs = thisData;
    float meanExp = thisEst;
    float meanUnc = thisEstErr;    
    int N=100000;
    
    TRandom3 randGen;

    unsigned int counterN(0);
    unsigned int counterD(0);

    bool doLeft = obs<=meanExp;

    for(int i=0; i<N; i++){
      double meanShift = randGen.Gaus(0.,(double)meanUnc);
  
      double rand = randGen.Poisson(double(meanExp)+meanShift); counterD++;

      if(doLeft)
	{
	  if(rand<=obs) counterN++;
	}
      else{
	if(rand>=obs) counterN++;
      }

      //cout << "rand " << i << " : " << rand.Poisson(meanExp) << endl;                                                                                                                
    }


    double prob=1.0*counterN/counterD;
    // double significance  = TMath::NormQuantile(1-prob);

//    std::cout << "probability: " << prob  << std::endl;
//    std::cout << "significance: " << significance << std::endl;
//    
//    std::cout << "Bin: " << iBin << std::endl;
//    std::cout << "Obs: " << obs << std::endl;
//    std::cout << "exp: " << meanExp << "\t" << meanUnc << std::endl;
    
    hPvalue->Fill(prob);

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

  TCanvas* c2;
  if(drawSignals) c2 = new TCanvas("c2", "", 1200, 600);
  //if(drawSignals) c2 = new TCanvas("c2", "", 1300, 800);
  else c2 = new TCanvas("c2", "", 1100, 600);

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
  pad1->SetBottomMargin(0.18);
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
  
  float yMin;
  if(drawSignals){
    yMin = 1e-2;
    yMax*=1000;
  }
//  if(drawSignals){
//    yMin = 1e-1;
//    yMax*=10;
//  }
  else{
    yMin = 1e-1;
    yMax*=20.;
  }

  int thisBin=63;
  
  hestimate_all->GetXaxis()->SetRangeUser(0, thisBin);
  hdata->GetXaxis()->SetRangeUser(0, thisBin);
  gdata->GetXaxis()->SetRangeUser(0, thisBin);
  gdata_zero->GetXaxis()->SetRangeUser(0, thisBin);
  for(int s=0; s<sigSize+sigContSize; ++s)
    hsig[s]->GetXaxis()->SetRangeUser(0, thisBin);
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);
//  hestimate_all->GetXaxis()->LabelsOption("v");
//  hestimate_all->GetXaxis()->SetLabelSize(0.035);
  hestimate_all->GetXaxis()->SetLabelSize(0.043);
  hestimate_all->GetXaxis()->SetLabelFont(62);
  hestimate_all->SetFillStyle(3244);
  hestimate_all->SetFillColor(kGray+2);

  TGraphAsymmErrors* g_data = new TGraphAsymmErrors(0);
  for( int iBin=1; iBin<(hdata->GetXaxis()->GetNbins()+1); ++iBin ) {

    double y;
    double x, xerr;

    x = hdata->GetBinCenter(iBin);
    xerr = hdata->GetBinWidth(iBin)/2.;

    y = hdata->GetBinContent(iBin);
    double yerr = hdata->GetBinError(iBin);

    int thisPoint = g_data->GetN();
    g_data->SetPoint( thisPoint, x, y );
    g_data->SetPointError( thisPoint, xerr, xerr, yerr, yerr );

  }
  
  hestimate_all->GetYaxis()->SetTitleOffset(0.95);
  hestimate_all->GetYaxis()->SetLabelSize(0.042);

  if(drawSignals)
    bgStack.Add(hsig[S]);

  hestimate_all->Draw("");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  //  hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");

  //  TH1D* hsig5 = (TH1D*) hsig[5]->Clone("hsig5");

//  if(drawSignals)
//    for(int s=0; s<sigSize+sigContSize; ++s)
//      hsig[s]->Draw("hist,same");
  
  TH1D* postfit = new TH1D("postfit", "", 1, 0, 1);
  postfit->SetFillColor(0);
  postfit->SetLineColor(0);

  TLegend* legend;
  if(drawSignals)
    legend = new TLegend( 0.7, 0.9-(bgSize+1-1)*0.06-0.06-0.06-0.06-0.03, 0.85, 0.9-0.06 );
  else
    legend = new TLegend( 0.8, 0.9-(bgSize+1-1)*0.06-0.06-0.02, 0.93, 0.9-0.06 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( hdata, "Data", "PL" );
  //  legend->AddEntry( postfit, "Post-fit SM", "F" );
  legend->AddEntry( hestimate[0], "Multijet", "F");
  legend->AddEntry( hestimate[1], "Lost lepton", "F");
  legend->AddEntry( hestimate[2], "Z #rightarrow #nu#bar{#nu}", "F");

//  TLegend* legendA = new TLegend( 0.1, 0.9-(6)*0.06-0.06, 0.99, 0.99 );
//  legendA->SetTextSize(0.042);
//  legendA->SetTextFont(42);
//  legendA->SetFillColor(0);
//  legendA->SetBorderSize(0);
//  legendA->AddEntry( hdata, "Data", "PL" );
//  legendA->AddEntry( hestimate[0], "Multijet", "F");
//  legendA->AddEntry( hestimate[1], "Lost Lepton", "F");
//  legendA->AddEntry( hestimate[2], "Invisible Z", "F");
//  for(int s=0; s<sigSize+sigContSize; ++s)
//    legendA->AddEntry( hsig[s], analysesSignal[s]->getName().c_str(), "L" );


//  if( drawSignals ){
//    if(S<sigSize)
//      legend->AddEntry( hsig[S], analysesSignal[S]->getName().c_str(), "F" );
//    else if(S==sigSize)
//      legend->AddEntry( hsig[S], "pp #rightarrow #tilde{g}#bar{#tilde{g}}, #tilde{g} #rightarrow t#bar{t}#tilde{#chi}_{1}^{0}", "F" );
//    else if(S==sigSize+1)
//      legend->AddEntry( hsig[S], "pp #rightarrow #tilde{g}#bar{#tilde{g}}, #tilde{g} #rightarrow t#bar{t}#tilde{#chi}_{1}^{0}", "F" );
//    else
//      legend->AddEntry( hsig[S], analysesSignalCont[S-sigSize]->getName().c_str(), "F" );
//    
//    if(sigName[S]=="T1bbbb_700_600"){
//      legend->AddEntry( postfit, "m_{#tilde{g}} = 700 GeV", "F" );
//      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 600 GeV", "F" );
//    }
//    else if(sigName[S]=="T1qqqq_700_600"){
//      legend->AddEntry( postfit, "m_{#tilde{g}} = 700 GeV", "F" );
//      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 600 GeV", "F" );
//    }
//    else if(sigName[S]=="T1tttt_700_400"){
//      legend->AddEntry( postfit, "m_{#tilde{g}} = 700 GeV", "F" );
//      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 400 GeV", "F" );
//    }
//    else if(sigName[S]=="T2bb_400_200"){
//      legend->AddEntry( postfit, "m_{#tilde{b}} = 400 GeV", "F" );
//      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 200 GeV", "F" );
//    }
//    else if(sigName[S]=="T2qq_500_300"){
//      legend->AddEntry( postfit, "m_{#tilde{q}} = 500 GeV", "F" );
//      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 300 GeV", "F" );
//    }
//    else if(sigName[S]=="T2tt_200_100"){
//      legend->AddEntry( postfit, "m_{#tilde{t}} = 200 GeV", "F" );
//      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 100 GeV", "F" );
//    }
//
//  }
  //  legend->AddEntry( postfit, "(700, 600)", "F" );
  

  //  TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi);
  TPaveText* labelTop = MT2DrawTools::getLabelTop(lumi);
  labelTop->Draw("same");
  
  //  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Supplementary");
  TPaveText* labelCMS = MT2DrawTools::getLabelCMS("CMS Preliminary");
  labelCMS->Draw("same");
  

  int nHTRegions = 6;
  std::vector< std::string > htRegions;
  htRegions.push_back("1 Jet");
////  htRegions.push_back("very low H_{T}");
////  htRegions.push_back("low H_{T}");
////  htRegions.push_back("medium H_{T}");
////  htRegions.push_back("high H_{T}");
////  htRegions.push_back("extreme H_{T}");
//  htRegions.push_back("#it{H}_{T} [200,450] GeV");
//  htRegions.push_back("#it{H}_{T} [450,575] GeV");
//  htRegions.push_back("#it{H}_{T} [575,1000] GeV");
//  htRegions.push_back("#it{H}_{T} [1000,1500] GeV");
//  htRegions.push_back("#it{H}_{T} >1500 GeV");
  htRegions.push_back("H_{T} [250,450]");
  htRegions.push_back("H_{T} [450,575]");
  htRegions.push_back("H_{T} [575,1000]");
  htRegions.push_back("H_{T} [1000,1500]");
  htRegions.push_back("H_{T} > 1500");


  TPaveText* htBox[5];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){

    if (iHT==0) htBox[iHT] = new TPaveText(0.12+0.15*iHT, 0.9-0.06+0.02, 0.34+0.15*iHT, 0.85+0.02, "brNDC");
    else if (iHT==1) htBox[iHT] = new TPaveText(0.30, 0.9-0.06+0.02, 0.39, 0.85+0.02, "brNDC");
    else htBox[iHT] = new TPaveText(0.39+0.14*(iHT-2), 0.9-0.06+0.02, 0.39+0.14+0.14*(iHT-2), 0.85+0.02, "brNDC");
    htBox[iHT]->AddText( htRegions[iHT].c_str() );

    htBox[iHT]->SetBorderSize(0);
    htBox[iHT]->SetFillColor(kWhite);
    htBox[iHT]->SetTextSize(0.035);
    htBox[iHT]->SetTextAlign(21); // align centered
    htBox[iHT]->SetTextFont(62);
    htBox[iHT]->Draw("same");

  // TPaveText* htBox[6];
  // for( int iHT = 0; iHT < nHTRegions; ++iHT){

  //   if (iHT==0) htBox[iHT] = new TPaveText(0.12+0.15*iHT, 0.9-0.06+0.02, 0.34+0.15*iHT, 0.85+0.02, "brNDC");
  //   else htBox[iHT] = new TPaveText(0.13+0.13*iHT, 0.9-0.06+0.02, 0.34+0.13*iHT, 0.85+0.02, "brNDC");
  //   htBox[iHT]->AddText( htRegions[iHT].c_str() );

  //   htBox[iHT]->SetBorderSize(0);
  //   htBox[iHT]->SetFillColor(0);
  //   htBox[iHT]->SetTextSize(0.035);
  //   htBox[iHT]->SetTextAlign(21); // align centered
  //   //    htBox[iHT]->SetTextFont(62);
  //   htBox[iHT]->Draw("same");

  }
  

  TLine* lHT[5];
  for( int iHT=0; iHT < 5; iHT++ ){
    if( iHT==0)
      lHT[iHT-1] = new TLine(12+11*(iHT), 0.0, 12+11*(iHT), yMax );
    else if (iHT!=1)
      lHT[iHT-1] = new TLine(12-4+11*iHT, 0.0, 12-4+11*iHT, yMax );
    else
      lHT[iHT-1] = new TLine(12+7*iHT, 0.0, 12+7*iHT, yMax );
 
    lHT[iHT-1]->SetLineColor(kBlack);
    lHT[iHT-1]->SetLineStyle(3);
    lHT[iHT-1]->SetLineWidth(2);

    lHT[iHT-1]->Draw("same");
  }

  legend->Draw("same");

  float left = pad1->GetLeftMargin();
  float right = pad1->GetRightMargin();
  float bot = pad1->GetBottomMargin();
  float top = pad1->GetTopMargin();
  float binWidth = (1.0-right-left)/63;

  TLatex* text = new TLatex();
  text->SetNDC(1);

  text->SetTextAlign(23);
  text->SetTextFont(42);
  text->SetTextAngle(0);
  text->SetTextSize(0.05);
  text->DrawLatex((1-left-right)/2+0.15,1-top-0.12, "Post-fit background");
 
  gPad->RedrawAxis();
  
  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();
  pad2->cd();  

  bool doLogRatio = false;
  
  TH2D* h2_axes_ratio;
  if(doLogRatio){
    gPad->SetLogy();
    h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0.1, 10.0 );  
    h_Ratio->GetYaxis()->SetRangeUser(0.1, 10.0);
    g_Ratio->GetYaxis()->SetRangeUser(0.1, 10.0);
  }
  else{
    //    h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0.5, 1.5 );
    h2_axes_ratio = new TH2D("axes_ratio", "", 10, 0, thisBin, 10, 0.0, 2.0 );
  }

  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.18);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.26);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitle("Data/Est.");
  
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


  TLine* lHT_b[6];
  for( int iHT=1; iHT < 6; iHT++ ){
    //    lHT_b[iHT-1] = new TLine(12+11*(iHT-1), 0, 12+11*(iHT-1), 3.0 );
    if(doLogRatio){
      //	lHT_b[iHT-1] = new TLine(12+11*(iHT-1), 0.5, 12+11*(iHT-1), 1.5 );
      if (iHT!=2)
	lHT_b[iHT-1] = new TLine(12-4+11*(iHT-1), 0.5, 12-4+11*(iHT-1), 1.5 );
      else
	lHT_b[iHT-1] = new TLine(12+7*(iHT-1), 0.5, 12+7*(iHT-1), 1.5 );
      if( iHT==1)
	lHT_b[iHT-1] = new TLine(12+11*(iHT-1), 0.5, 12+11*(iHT-1), 1.5 );
    }
    else{
      //      lHT_b[iHT-1] = new TLine(12+11*(iHT-1), 0, 12+11*(iHT-1), 1.5 );
      if (iHT!=2)
	lHT_b[iHT-1] = new TLine(12-4+11*(iHT-1), 0, 12-4+11*(iHT-1), 2 );
      else
	lHT_b[iHT-1] = new TLine(12+7*(iHT-1), 0, 12+7*(iHT-1), 2 ); 
      if( iHT==1)
	lHT_b[iHT-1] = new TLine(12+11*(iHT-1), 0, 12+11*(iHT-1), 2 );  
    }

    lHT_b[iHT-1]->SetLineColor(kBlack);
    lHT_b[iHT-1]->SetLineStyle(3);
    lHT_b[iHT-1]->SetLineWidth(2);

    lHT_b[iHT-1]->Draw("same");
  }
  
  gPad->RedrawAxis();

  c2->cd();
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.pdf", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.png", fullPath.c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate.eps", fullPath.c_str()) );

  c2->Clear();

  c2->cd();  

  yMin=1e-1;
  hestimate_all->GetYaxis()->SetRangeUser(yMin, yMax);

  gPad->SetLogy();
  
  thisBin=63;

  bgStack.Add(hsig[S]);

  hestimate_all->Draw("");
  bgStack.Draw("histo, same");
  hestimate_all->Draw("E2,same");
  //hdata->Draw("pe,same");
  gdata_zero->Draw("pe,same");
  gdata->Draw("pe,same");
  
  labelTop->Draw("same");
  labelCMS->Draw("same");

  TPaveText* htBox2[5];
  for( int iHT = 0; iHT < nHTRegions; ++iHT){


    if (iHT==0) htBox2[iHT] = new TPaveText(0.12+0.15*iHT, 0.9-0.06+0.02, 0.34+0.15*iHT, 0.85+0.02, "brNDC");
    else if (iHT==1) htBox2[iHT] = new TPaveText(0.30, 0.9-0.06+0.02, 0.39, 0.85+0.02, "brNDC");
    else htBox2[iHT] = new TPaveText(0.39+0.14*(iHT-2), 0.9-0.06+0.02, 0.39+0.14+0.14*(iHT-2), 0.85+0.02, "brNDC");
    htBox2[iHT]->AddText( htRegions[iHT].c_str() );

//    if (iHT==0) htBox2[iHT] = new TPaveText(0.12+0.15*iHT, 0.9-0.04+0.01, 0.34+0.15*iHT, 0.87+0.01, "brNDC");
//    else htBox2[iHT] = new TPaveText(0.13+0.13*iHT, 0.9-0.04+0.01, 0.34+0.13*iHT, 0.87+0.01, "brNDC");
//    htBox2[iHT]->AddText( htRegions[iHT].c_str() );

    htBox2[iHT]->SetBorderSize(0);
    htBox2[iHT]->SetFillColor(kWhite);
    htBox2[iHT]->SetTextSize(0.03);
    htBox2[iHT]->SetTextAlign(21); // align centered
    htBox2[iHT]->SetTextFont(62);
    htBox2[iHT]->Draw("same");

  }
  
  for( int iHT=0; iHT < 5; iHT++ ){
    lHT[iHT-1]->Draw("same");
  }

  if( drawSignals ){
    if(S<sigSize)
      legend->AddEntry( hsig[S], analysesSignal[S]->getName().c_str(), "F" );
    else if(S==sigSize)
      legend->AddEntry( hsig[S], "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t}#tilde{#chi}_{1}^{0}", "F" );
    else if(S==sigSize+1)
      legend->AddEntry( hsig[S], "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t}#tilde{#chi}_{1}^{0}", "F" );
    else
      legend->AddEntry( hsig[S], analysesSignalCont[S-sigSize]->getName().c_str(), "F" );
    
    if(sigName[S]=="T1bbbb_700_600"){
      legend->AddEntry( postfit, "m_{#tilde{g}} = 1000 GeV", "F" );
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 800 GeV", "F" );
    }
    else if(sigName[S]=="T1qqqq_700_600"){
      legend->AddEntry( postfit, "m_{#tilde{g}} = 700 GeV", "F" );
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 600 GeV", "F" );
    }
    else if(sigName[S]=="T1tttt_700_400"){
      legend->AddEntry( postfit, "m_{#tilde{g}} = 700 GeV", "F" );
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 400 GeV", "F" );
    }
    else if(sigName[S]=="T2bb_400_200"){
      legend->AddEntry( postfit, "m_{#tilde{b}} = 400 GeV", "F" );
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 200 GeV", "F" );
    }
    else if(sigName[S]=="T2qq_500_300"){
      legend->AddEntry( postfit, "m_{#tilde{q}} = 500 GeV", "F" );
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 300 GeV", "F" );
    }
    else if(sigName[S]=="T2tt_200_100"){
      legend->AddEntry( postfit, "m_{#tilde{t}} = 200 GeV", "F" );
      legend->AddEntry( postfit, "m_{#tilde{#chi}_{1}^{0}}= 100 GeV", "F" );
    }

  }
  
  legend->Draw("same");
  
  left  = gPad->GetLeftMargin();
  right = gPad->GetRightMargin();
  bot   = gPad->GetBottomMargin();
  top   = gPad->GetTopMargin();
  binWidth = (1.0-right-left)/63;

  text->SetNDC(1);

  text->SetTextAlign(23);
  text->SetTextFont(42);
  text->SetTextAngle(0);
  text->SetTextSize(0.05);
  text->DrawLatex((1-left-right)/2+0.15,1-top-0.12, "Post-fit background");

  gPad->RedrawAxis();
  
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate_plus%s.pdf", fullPath.c_str(), sigName[S].c_str()) );
  c2->SaveAs( Form("%s/mt2_ALL_fullEstimate_plus%s.png", fullPath.c_str(), sigName[S].c_str()) );


//  TCanvas* cLegend = new TCanvas("cLegend", "", 200, 400);
//  cLegend->cd();
//  legendA->Draw();
//  cLegend->cd();
//  cLegend->SaveAs( Form("%s/LegendAll.pdf", fullPath.c_str()) );
//  cLegend->SaveAs( Form("%s/LegendAll.png", fullPath.c_str()) );
  

  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);
  
  TF1* fgauss= new TF1("fgauss", "gaus", -5, 5);
  fgauss->SetLineColor(2);

  TCanvas* c3 = new TCanvas("c3", "", 600, 600);
  c3->cd();
  hPull->SetStats(1110);
  //  hPull->GetYaxis()->SetRangeUser(0, 15);
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
