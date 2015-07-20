#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"

#include "TRandom3.h"

#define mt2_cxx
#include "interface/mt2.h"

double lumi=4.; //fb-1
bool doNminusOne=false;

void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* > bgYields,  std::string var, std::string label, std::string selection, bool logY );

int main( int argc, char* argv[] ) {

  MT2DrawTools::setStyle();

  if( argc!=2 ) {
    std::cout << "USAGE: ./scalesVariation [dir]" << std::endl;
    exit(113);
  }


  std::string dir( argv[1] );
  std::string mc_fileName = dir + "/analyses.root";

  std::string outputdir = dir + "/scalesVariations";

  MT2Analysis<MT2EstimateTree>* data  = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "data" );
  
    MT2Analysis<MT2EstimateTree>* top   = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "Top");
  //  MT2Analysis<MT2EstimateTree>* wjets   = MT2Analysis<MT2EstimateTree>::readFromFile( mc_fileName, "WJets");
  
  std::vector<MT2Analysis<MT2EstimateTree>* > bgYields;
  bgYields.push_back(top);
  //bgYields.push_back(wjets);
  
  drawHisto( outputdir, data, bgYields,  "mt2", "M_{T2} [GeV]", "", kFALSE );
//  drawHisto( outputdir, data, bgYields, "ht", 100, 0., 2500., "H_{T} [GeV]", "", kTRUE );
//  drawHisto( outputdir, data, bgYields, "nJets", 12, 0, 12, "N(jet)", "", kFALSE );
//  drawHisto( outputdir, data, bgYields, "nBJets", 6, 0, 6, "N(b-tag)", "", kFALSE );

  return 0;

}// main



void drawHisto( const std::string& outputdir, MT2Analysis<MT2EstimateTree>* data, std::vector< MT2Analysis<MT2EstimateTree> *> bgYields, std::string var, std::string label, std::string selection, bool logY ) {

  MT2DrawTools::setStyle();

  std::vector<int> colors; //mc
  colors.push_back(855); // top
  //colors.push_back(417); // wjets
    
  std::string fullPath = outputdir;
  system( Form("mkdir -p %s", fullPath.c_str()) );

  std::set<MT2Region> MT2Regions = data->getRegions();

  //TFile* fout = new TFile(Form("%s/MMHT2014.root", outputdir.c_str()), "recreate");  
  TFile* fout = new TFile(Form("%s/CT10.root", outputdir.c_str()), "recreate");  
  //TFile* fout = new TFile(Form("%s/NNPDF.root", outputdir.c_str()), "recreate");  
  TFile* fout2= new TFile(Form("%s/allHistos.root", outputdir.c_str()), "recreate");
    
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    MT2Region thisRegion( (*iMT2) );
   
    TH1F::AddDirectory(kTRUE);
    
    int nBins;
    double* bins;
    thisRegion.getBins(nBins, bins);

      
    TH1D* h_data = new TH1D("h_data", "", nBins, bins);
    
    // All selections:
    std::string MinusOneCut=selection;
    int nCuts = 8;
    std::string cut[nCuts];
    cut[0] = "((ht>1000. && met>30.) || (ht>450. && met>200.))";
    cut[1] = "nVert>0";
    cut[2] = "nJets>= 2";
    cut[3] = "nElectrons10==0 && nMuons10==0";
    cut[4] = "nPFLep5LowMT==0 && nPFHad10LowMT==0";
    cut[5] = "dPhiMin>0.3";
    cut[6] = "diffMetMht/met<0.5";
    cut[7] = "mt2>200.";

    std::string scale_=selection;
    if( !doNminusOne ){
      
      if( selection == "" ) selection = Form("(mt2>200.)*weight*%0.f", lumi);
      else selection = Form("(mt2>200. &&  %s)*weight*%0.f", selection.c_str(), lumi);
      
    }
    else{
      
      selection = "";
      for( int c=0; c < nCuts; ++c ){
	
	if( cut[c] == MinusOneCut );
	else selection += cut[c];
	
      }
      
      selection = Form("(%s)*weight*%0.f", selection.c_str(), lumi);
      
    }
    
    data->get( thisRegion )->tree->Project( h_data->GetName(), var.c_str(), selection.c_str() );
    
    TGraphAsymmErrors* gr_data = MT2DrawTools::getPoissonGraph(h_data);
    gr_data->SetMarkerStyle(20);
    gr_data->SetMarkerSize(1.6);
    
    THStack bgStack("bgStack", "");
    TH1D* h_bg[bgYields.size()];
    
//    int nVariations=50+1; //NNPDF3.0 LO
//    int firstId=10;
    int nVariations=25+1; //MMHT2014 LO
    int firstId=316;
//    int nVariations=1; //CTEQ61L
//    int firstId=315;
//    int nVariations=1; //NNPDF3.0 NLO
//    int firstId=391;
//    int nVariations=26+1; //CT10 NLO
//    int firstId=393;
//    int nVariations=1; //MMHT2014nlo68cl
//    int firstId=446;
    std::vector<std::string> variations;
    //    variations.push_back("NNPDF3.0 LO, #alpha_{s}=0.130");
    variations.push_back("MMHT2014 LO 68% CL");
    //    variations.push_back("CTEQ61L");
    //    variations.push_back("NNPDF3.0 NLO, #alpha_{s}=0.118");
    //    variations.push_back("CT10 NLO");
    //    variations.push_back("MMHT2014 NLO 68% CL");
    
    TH1D* h_bg_v[bgYields.size()][nVariations*2-1];
    TH1D* h_bg_r[bgYields.size()][nVariations*2-1];

    TGraphAsymmErrors* g_v[bgYields.size()][nVariations*2-1];
    TGraphAsymmErrors* g_r[bgYields.size()][nVariations*2-1];
    
    double normalization;
    float binContent[nVariations*2-1][nBins];

    float sigma2[nBins];
    float sigma2p[nBins];
    float sigma2m[nBins];

    float sigma[nBins];
    float sigmap[nBins];
    float sigmam[nBins];

    int Nplus[nBins];
    int Nminus[nBins];
    
    for( unsigned i=0; i<bgYields.size(); ++i ) {
      //for( unsigned i=1; i<bgYields.size(); ++i ) { //Not drawing QCD
      
      for(int b=0; b<nBins; ++b){
	
	for( int v=0; v<nVariations*2-1; ++v)
	  binContent[v][b]=0.;
	
	sigma2[b]=0.;
	sigma[b]=0.;

	sigma2p[b]=0.;
	sigmap[b]=0.;

	sigma2m[b]=0.;
	sigmam[b]=0.;
	
	Nplus[b]=0;
	Nminus[b]=0;

      }
      
      normalization=1.0;
      
      int index=i;
      
      //int index = bgYields.size() - i - 1; // reverse ordered stack is prettier with QCD
      
      h_bg[index] = new TH1D(Form("h_%s_%s", var.c_str(), bgYields[i]->getName().c_str()), "", nBins, bins);      
      h_bg[index]->Sumw2();

      bgYields[index]->get(thisRegion)->tree->Project( h_bg[index]->GetName(), var.c_str(), selection.c_str() );
            
      std::vector<int> colors_v;
      //colors_v.push_back(4); //NNPDF                                                                                                                                                                       
      colors_v.push_back(6); //CT10                                                                                                                                                                          
      //colors_v.push_back(3); //MMHT2014 

      h_bg_v[index][0] = new TH1D(Form("h_%s_%s_v%d", var.c_str(), bgYields[i]->getName().c_str(), firstId), "", nBins, bins);
      h_bg_v[index][0]->Sumw2();
      
      std::string scale = scale_;

      if( scale == "" ) scale = Form("(mt2>200. && LHEweight_id==%d)*LHEweight_wgt/LHE_originalWeight*weight*%0.f", firstId, lumi);
      else scale = Form("(mt2>200. &&  %s && LHEweight_id==%d)*LHEweight_wgt/LHE_originalWeight*weight*%0.f", scale.c_str(), firstId, lumi);
      bgYields[index]->get(thisRegion)->tree->Project( h_bg_v[index][0]->GetName(), var.c_str(), scale.c_str() );

      normalization=h_bg[index]->Integral();
      
      double error=0.;
      double thisIntegral = h_bg_v[index][0]->IntegralAndError(1, -1, error);
      double thisScale = normalization/thisIntegral;

      h_bg_v[index][0]->Scale(thisScale);

      h_bg_v[index][0]->SetLineColor(colors_v[0]);
      h_bg_v[index][0]->SetLineWidth(2);
      h_bg_v[index][0]->SetMarkerColor(colors_v[0]);
      
      h_bg_r[index][0]=(TH1D*) h_bg_v[index][0]->Clone(Form("h_%s_%s_r%d", var.c_str(), bgYields[i]->getName().c_str(), firstId));
      h_bg_r[index][0]->Sumw2();
      
      h_bg_r[index][0]->Divide(h_bg[index]);
      
      for( int b=0; b<nBins; ++b )
	binContent[0][b] = h_bg_v[index][0]->GetBinContent(b+1);

      
      for(int v=1; v < nVariations*2-1; ++v){
	
	colors_v.push_back(2);
	colors_v.push_back(2);

	scale = scale_;

	h_bg_v[index][v] = new TH1D(Form("h_%s_%s_v%d", var.c_str(), bgYields[i]->getName().c_str(), v+firstId), "", nBins, bins);
	h_bg_v[index][v]->Sumw2();
	
	if( scale == "" ) scale = Form("(mt2>200. && LHEweight_id==%d)*LHEweight_wgt/LHE_originalWeight*weight*%0.f", v+firstId, lumi);
	else scale = Form("(mt2>200. &&  %s && LHEweight_id==%d)*LHEweight_wgt/LHE_originalWeight*weight*%0.f", scale.c_str(), v+firstId, lumi);
	bgYields[index]->get(thisRegion)->tree->Project( h_bg_v[index][v]->GetName(), var.c_str(), scale.c_str() );

	//normalization=h_bg_v[index][0]->Integral();

        error=0.;
	thisIntegral = h_bg_v[index][v]->IntegralAndError(1, -1, error);
	thisScale = normalization/thisIntegral;
	
	h_bg_v[index][v]->Scale(thisScale);


	scale = scale_;

	h_bg_v[index][v+1] = new TH1D(Form("h_%s_%s_v%d", var.c_str(), bgYields[i]->getName().c_str(), v+1+firstId), "", nBins, bins);
	h_bg_v[index][v+1]->Sumw2();
	
	if( scale == "" ) scale = Form("(mt2>200. && LHEweight_id==%d)*LHEweight_wgt/LHE_originalWeight*weight*%0.f", v+1+firstId, lumi);
	else scale = Form("(mt2>200. &&  %s && LHEweight_id==%d)*LHEweight_wgt/LHE_originalWeight*weight*%0.f", scale.c_str(), v+1+firstId, lumi);
	bgYields[index]->get(thisRegion)->tree->Project( h_bg_v[index][v+1]->GetName(), var.c_str(), scale.c_str() );

	//normalization=h_bg_v[index][0]->Integral();

        error=0.;
	thisIntegral = h_bg_v[index][v+1]->IntegralAndError(1, -1, error);
	thisScale = normalization/thisIntegral;
	
	h_bg_v[index][v+1]->Scale(thisScale);

	double err0=0.;
	double integral0 = h_bg_v[index][0]->IntegralAndError(1, -1, err0);
	
	bool doAdd = true;
	//if( (error/thisIntegral)/(err0/integral0) > 2. ) doAdd = false;

	for( int b=0; b<nBins; ++b ){
	  
	  binContent[v][b] = h_bg_v[index][v]->GetBinContent(b+1);
	  binContent[v+1][b] = h_bg_v[index][v+1]->GetBinContent(b+1);
	  //sigma2[b] = (binContent[v][b]-binContent[0][b])*(binContent[v][b]-binContent[0][b]);
	  
	  float localSigma = 0;
	  float localSigmap = 0;
	  float localSigmam = 0;
	  
	  if(binContent[v][b]-binContent[0][b] > 0. || binContent[v+1][b]-binContent[0][b] > 0.){
	    
	    localSigmap = (binContent[v][b]-binContent[0][b] > binContent[v+1][b]-binContent[0][b]) ? binContent[v][b]-binContent[0][b] : binContent[v+1][b]-binContent[0][b];
	   	    
	    ++Nplus[b];
	    
	    sigma2p[b] += localSigmap*localSigmap;
	    
	    
	  } else if(binContent[v][b]-binContent[0][b] < 0. || binContent[v+1][b]-binContent[0][b] < 0.){
	    
	    localSigmam = (binContent[0][b]-binContent[v][b] > binContent[0][b]-binContent[v+1][b]) ? binContent[0][b]-binContent[v][b] : binContent[0][b]-binContent[v+1][b];
	    
	    ++Nminus[b];
	    
	    sigma2m[b] += localSigmam*localSigmam;
	    
	  }  
	 
//	  std::cout << sigma2p[b] << std::endl;
//	  std::cout << sigma2m[b] << std::endl;
	  
	}
	
	h_bg_v[index][v]->SetLineColor(colors_v[v]);
	h_bg_v[index][v]->SetLineWidth(2);
	h_bg_v[index][v]->SetMarkerColor(colors_v[v]);
	
	h_bg_r[index][v]=(TH1D*) h_bg_v[index][v]->Clone(Form("h_%s_%s_r%d", var.c_str(), bgYields[i]->getName().c_str(), v+firstId));
	h_bg_r[index][v]->Sumw2();
	
	//	h_bg_r[index][v]->Divide(h_bg_v[index][0]);
	h_bg_r[index][v]->Divide(h_bg[index]);

	h_bg_v[index][v+1]->SetLineColor(colors_v[v+1]);
	h_bg_v[index][v+1]->SetLineWidth(2);
	h_bg_v[index][v+1]->SetMarkerColor(colors_v[v+1]);
	
	h_bg_r[index][v+1]=(TH1D*) h_bg_v[index][v]->Clone(Form("h_%s_%s_r%d", var.c_str(), bgYields[i]->getName().c_str(), v+1+firstId));
	h_bg_r[index][v+1]->Sumw2();
	
	//	h_bg_r[index][v+1]->Divide(h_bg_v[index][0]);
	h_bg_r[index][v+1]->Divide(h_bg[index]);
	
	//	std::cout << "variations: " << v <<","<< v+1 << std::endl;

	v+=1;
		

      }
      
      h_bg[index]->SetFillColor( colors[index] );
      h_bg[index]->SetLineColor( kBlack );
      bgStack.Add(h_bg[index]);
      
    }
    
    g_v[0][0]=new TGraphAsymmErrors(h_bg_v[0][0]);
    g_r[0][0]=new TGraphAsymmErrors(h_bg_r[0][0]);
    
    TH1F* horiginal = (TH1F*) h_bg[0]->Clone(Form("horiginal_%s", thisRegion.getName().c_str()));
    TH1F* hcentral = (TH1F*) h_bg_v[0][0]->Clone(Form("hcentral_%s", thisRegion.getName().c_str()));
    TH1F* hup=new TH1F(Form("hup_%s", thisRegion.getName().c_str()), "", nBins, bins);
    TH1F* hdn=new TH1F(Form("hdn_%s", thisRegion.getName().c_str()), "", nBins, bins);
    
    for( int b=0; b<nBins; ++b){

      std::cout << Nplus[b] << "\t" << Nminus[b] << std::endl;

//      //NNPDF
//      //      sigma[b]=sqrt(1./(nVariations*1.-2.)*sigma2[b]);
//      sigmap[b]=sqrt(1./(Nplus[b]*1.-1.)*sigma2p[b]);
//      sigmam[b]=sqrt(1./(Nminus[b]*1.-1.)*sigma2m[b]);

      //MMHT2014
      //      sigma[b]=sqrt(sigma2[b]);
      sigmap[b]=sqrt(sigma2p[b]);
      sigmam[b]=sqrt(sigma2m[b]);
      
//      //CT10
//      //      sigma[b]=sqrt(sigma2[b]);
//      sigmap[b]=sqrt(sigma2p[b])/1.645;
//      sigmam[b]=sqrt(sigma2m[b])/1.645;
      
      //h_bg_v[0][0]->SetBinError(b+1, sigma[b]);
      //h_bg_r[0][0]->SetBinError(b+1, sigma[b]/binContent[0][b]);
      
//      g_v[0][0]->SetPointEYhigh(b, sigma[b]);
//      g_v[0][0]->SetPointEYlow(b, sigma[b]);
      g_v[0][0]->SetPointEYhigh(b, sigmap[b]);
      g_v[0][0]->SetPointEYlow(b, sigmam[b]);

      h_bg_r[0][0]->SetBinError(b+1, 0.);
//      g_r[0][0]->SetPointEYhigh(b, sigma[b]/binContent[0][b]);
//      g_r[0][0]->SetPointEYlow(b, sigma[b]/binContent[0][b]);
      g_r[0][0]->SetPointEYhigh(b, sigmap[b]/binContent[0][b]);
      g_r[0][0]->SetPointEYlow(b, sigmam[b]/binContent[0][b]);

//      hup->SetBinContent(b+1, binContent[0][b]+sigma[b]);
//      hdn->SetBinContent(b+1, binContent[0][b]-sigma[b]);
      hup->SetBinContent(b+1, binContent[0][b]+sigmap[b]);
      hdn->SetBinContent(b+1, binContent[0][b]-sigmam[b]);
      
      //  std::cout << sigmap[b] << "\t" << sigmam[b] << std::endl;
      
    }

    fout->cd();
    hup->Write();
    hdn->Write();
    hcentral->Write();
    horiginal->Write();
    
    TCanvas* c1 = new TCanvas( "c1", "", 600, 700 );
    c1->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3-0.1,1,1);
    pad1->SetBottomMargin(0.15);
    pad1->Draw();
    pad1->cd();

    float yMin = 0.;
    float yMax1 = h_bg_v[0][0]->GetMaximum()*1.5;
    float yMax2 = 1.5*(h_bg_v[0][0]->GetMaximum() + sqrt(h_bg_v[0][0]->GetMaximum()));
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( h_bg_v[0][0]->GetNbinsX()<2 ) yMax *=3.;
    
    if(logY) {
      gPad->SetLogy();
      yMin=1e-1;
      yMax*=50.;
    }

    std::string labelY = "Events";
    TH2D* h_axes = new TH2D("axes", "", 10, 200, 1500, 10, yMin, yMax );
    h_axes->GetXaxis()->SetTitle(label.c_str());
    h_axes->GetYaxis()->SetTitle(labelY.c_str());
    h_axes->GetYaxis()->SetTitleOffset(1.5);
    h_axes->GetYaxis()->SetLabelSize(0.04);
    h_axes->Draw();

    std::vector<std::string> niceNames = thisRegion.getNiceNames();
    //  niceNames.push_back("M_{T2} > 200 GeV");
    for( unsigned i=0; i<niceNames.size(); ++i ) {

      float yMax = 0.9-(float)i*0.04;
      float yMin = yMax - 0.04;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      //TPaveText* regionText = new TPaveText( 0.18, yMin, 0.5, yMax, "brNDC" );
      //regionText->SetTextSize(0.035);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
      regionText->AddText( niceNames[i].c_str() );
      regionText->Draw("same");

    }
    
    TLegend* legend = new TLegend( 0.6, 0.9-(bgYields.size()+1)*0.04, 0.93, 0.9 );
    //legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);

    h_bg[0]->SetLineColor(1);
    h_bg[0]->SetFillColor(0);
      
    legend->AddEntry( h_bg[0], "Top", "L");
    legend->AddEntry( h_bg_v[0][0], variations[0].c_str(), "L");

    legend->Draw("same");
    
    h_bg[0]->Draw("histo,same");
    //    h_bg_v[0][0]->Draw("same");

    g_v[0][0]->Draw("P,same");

    fout2->cd();
    for(int v=0; v<nVariations*2-1; ++v){
      h_bg_v[0][v]->Draw("same");
      h_bg_v[0][v]->Write();
      
    }
    

    TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation(lumi); //Not drawing data yet
    labelTop->Draw("same");

    gPad->RedrawAxis();

	
    c1->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.21);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.1);
    pad2->Draw();
    pad2->cd();
    
    pad2->SetGridy();
    
    for( unsigned i=0; i<bgYields.size(); ++i )
      for(int v=0; v<nVariations*2-1; ++v){

	h_bg_r[i][v]->SetStats(0);
	h_bg_r[i][v]->GetXaxis()->SetLabelSize(0.00);
	h_bg_r[i][v]->GetXaxis()->SetTickLength(0.09);
	h_bg_r[i][v]->GetYaxis()->SetNdivisions(5,5,0);
	h_bg_r[i][v]->GetYaxis()->SetRangeUser(0.8,1.2);
	//h_bg_r[i][v]->GetYaxis()->SetRangeUser(0.,2.);
	h_bg_r[i][v]->GetYaxis()->SetTitleSize(0.17);
	h_bg_r[i][v]->GetYaxis()->SetTitleOffset(0.4);
	h_bg_r[i][v]->GetYaxis()->SetLabelSize(0.17);
	h_bg_r[i][v]->GetYaxis()->SetTitle("ratio");
       	
      }
    
    h_bg_r[0][0]->Draw("P");
    g_r[0][0]->Draw("P,same");
    
    gPad->RedrawAxis();
      
    c1->cd();
    
    
    c1->SaveAs( Form("%s/%s_%s.eps", fullPath.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.png", fullPath.c_str(), var.c_str(), thisRegion.getName().c_str()) );
    c1->SaveAs( Form("%s/%s_%s.pdf", fullPath.c_str(), var.c_str(), thisRegion.getName().c_str()) );

    delete c1;
    delete h_axes;
    
    delete h_data;
    
    delete bins;
    
    for(unsigned i=0; i<bgYields.size(); ++i){
      //for(unsigned i=1; i<bgYields.size(); ++i){ //Not drawing QCD
    
      delete h_bg[i];
      
      for(int v=0; v<nVariations*2-1; ++v){
	delete h_bg_v[i][v];
	delete h_bg_r[i][v];
      }
    
    }
    
    TH1F::AddDirectory(kFALSE);

  }// for MT2 regions
  
  fout->Close();
  fout2->Close();
  
}
