#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"

#include "interface/MT2Config.h"
#include "../interface/MT2DrawTools.h"


bool drawInclusivePlot = false;


void drawMT2vsVar( const std::string& outdir, TTree* tree, const std::string& varName, std::vector< std::pair<int,int> >& bins, std::pair<float,float>& htBin );
void drawYield( const std::string& outdir, TTree* tree, const std::string& suffix, const std::string& cuts, float binWidth=200 );

int main(  int argc, char* argv[]) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./drawMT2shape [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  MT2Config cfg(configFileName);


  std::string dir = cfg.getEventYieldDir();

  bool doZinv =0;

  MT2DrawTools::setStyle();

  TFile* file;
  TTree* tree;
  if( !doZinv ){
    file = TFile::Open(Form("%s/zllControlRegion/mc.root", dir.c_str()) );
    // file = TFile::Open("EventYields_data_Run2016_12p9ifb_40ifbScaled_noZinvExtrapol/zllControlRegion/mc.root");
    tree = (TTree*)file->Get("zllCR/HT250toInf_j1toInf_b0toInf/tree_zllCR_HT250toInf_j1toInf_b0toInf");
  }else{
    file = TFile::Open(Form("%s/ZJetsIncl.root", dir.c_str()) );
    //file = TFile::Open("EventYields_data_2016_SnTMC_362ifb_noBtagSF_ll1jpt/ZJetsIncl.root");
    tree = (TTree*)file->Get("ZJets_inclusive/HT250toInf_j1toInf_b0toInf/tree_ZJets_inclusive_HT250toInf_j1toInf_b0toInf");
  }
  // TFile* file = TFile::Open("EventYields_data_2016_SnTMC_362ifb_18ifbUB/zllControlRegion/mc.root");
  //TFile* file = TFile::Open("EventYields_data_Run2016_12p9ifb_40ifbScaled_noZinvExtrapol/zllControlRegion/mc.root");


  //  TTree* tree = (TTree*)file->Get("zllCR/HT250toInf_j1toInf_b0toInf/tree_zllCR_HT250toInf_j1toInf_b0toInf");
  //TTree* tree = (TTree*)file->Get("zllCR/HT200toInf_j1toInf_b0toInf/tree_zllCR_HT200toInf_j1toInf_b0toInf");

  tree->GetEntries();


  std::vector< std::pair<float,float> > htBins;
  if( drawInclusivePlot ) {
    htBins.push_back( std::pair<float,float>(250.,13000.) );
  } else {
    htBins.push_back( std::pair<float,float>(250.,450.) );
    htBins.push_back( std::pair<float,float>(450.,575.) );
    htBins.push_back( std::pair<float,float>(575.,1000.) );
    htBins.push_back( std::pair<float,float>(1000.,1500.) );
    htBins.push_back( std::pair<float,float>(1500.,13000.) );
  }

  std::string outdir = dir + "/mt2ShapeComparisons";

  if( doZinv )
    outdir += "_zinv";
  
  //  std::string outdir = "mt2ShapeComparisons";
  system( Form("mkdir -p %s", outdir.c_str()) );
  
  for( unsigned i=0; i<htBins.size(); ++i ) {

    std::vector< std::pair<int,int> > nJetsBins;
    nJetsBins.push_back( std::pair<int,int>(2,3) );
    nJetsBins.push_back( std::pair<int,int>(4,6) );
    nJetsBins.push_back( std::pair<int,int>(7,13000) );
    if( !drawInclusivePlot )
      drawMT2vsVar( outdir, tree, "nJets", nJetsBins, htBins[i] );

    std::vector< std::pair<int,int> > nBJetsBins;
    nBJetsBins.push_back( std::pair<int,int>(0,0) );
    nBJetsBins.push_back( std::pair<int,int>(1,1) );
    nBJetsBins.push_back( std::pair<int,int>(2,2) );
    if( drawInclusivePlot )
      nBJetsBins.push_back( std::pair<int,int>(3,13000) );
    drawMT2vsVar( outdir, tree, "nBJets", nBJetsBins, htBins[i] );

  }

  drawYield( outdir, tree, "ht1500", "ht>1500.", 200 );
  drawYield( outdir, tree, "ht1000_nJ23", "ht>1000. && ht<1500. && nJets>=2 && nJets<=3", 200 );
  drawYield( outdir, tree, "ht1000_nJ46", "ht>1000. && ht<1500. && nJets>=4 && nJets<=6", 200 );
  drawYield( outdir, tree, "ht1000_nJ7", "ht>1000. && ht<1500. && nJets>=7", 200 );
  drawYield( outdir, tree, "ht575_nJ23", "ht>575. && ht<1000. && nJets>=2 && nJets<=3", 200 );
  drawYield( outdir, tree, "ht575_nJ46", "ht>575. && ht<1000. && nJets>=4 && nJets<=6", 200 );
  drawYield( outdir, tree, "ht575_nJ7", "ht>575. && ht<1000. && nJets>=7", 200 );
  
  return 0;

}



void drawMT2vsVar( const std::string& outdir, TTree* tree, const std::string& varName, std::vector< std::pair<int,int> >& bins, std::pair<float,float>& htBin ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 800);
  c1->cd();
  c1->SetLogy();

  std::vector<int> colors;
  //colors.push_back(kBlack);
  colors.push_back(46);
  colors.push_back(38);
  colors.push_back(29);
  colors.push_back(42);


  TPad* pad1 = new TPad("pad1", "", 0, 0.3-0.1, 1, 0.99);
  pad1->SetBottomMargin(0.15);
  pad1->SetLogy();
  pad1->Draw();
  TPad* pad2 = new TPad("pad2", "", 0, 0, 1, 0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);
  pad2->Draw();


  pad2->cd();

  TH2D* h2_axes_ratio = new TH2D("axes_ratio", "", 10, 200., 1700., 10, 0.4, 1.6);
  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.4);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->SetYTitle("Ratio");
  h2_axes_ratio->Draw();

  TLine* line_one = new TLine( 200., 1., 1700., 1. );
  line_one->Draw("same");

  pad2->Draw();

  pad1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 200., 1700., 10, 0.0001, 1);
  h2_axes->SetXTitle("M_{T2} [GeV]");
  h2_axes->SetYTitle("Normalized to Unity");
  h2_axes->Draw();


  TLegend* legend = new TLegend( 0.5, 0.62, 0.9, 0.92 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);

  if( htBin.second<10000. )
    legend->SetHeader( Form("%.0f < H_{T} < %.0f GeV", htBin.first, htBin.second) );
  else {
    if( drawInclusivePlot ) {
      legend->SetHeader( Form("H_{T} > %.0f GeV, NJets #geq 3", htBin.first) );
    } else {
      legend->SetHeader( Form("H_{T} > %.0f GeV", htBin.first) );
    }
  }

  int numberOfBins = 10;

  TH1D* h1_ref = new TH1D( Form("ref_%s_ht%.0f", varName.c_str(), htBin.first), "", numberOfBins, 200., 1700. );
  h1_ref->Sumw2();

  for( unsigned i=0; i<bins.size(); ++i ) {

    TH1D* h1_mt2 = new TH1D( Form("mt2_%s_ht%.0f_%d", varName.c_str(), htBin.first, i), "", numberOfBins, 200., 1700. );
    h1_mt2->Sumw2();

    if( drawInclusivePlot ) {
      //  tree->Project( h1_mt2->GetName(), "mt2", Form("weight*( nJets>=3 && ht>%f && ht<=%f && %s>=%d && %s<=%d)", htBin.first, htBin.second, varName.c_str(), bins[i].first, varName.c_str(), bins[i].second) );
      tree->Project( h1_mt2->GetName(), "mt2", Form("weight*( nJets>=3 && ht>%f && ht<=%f && %s>=%d && %s<=%d)", htBin.first, htBin.second, varName.c_str(), bins[i].first, varName.c_str(), bins[i].second) );
      //tree->Project( h1_mt2->GetName(), "mt2", Form("weight*( Z_pt>200 && fabs(Z_mass-91.19)<20. && nJets>=3 && ht>%f && ht<=%f && %s>=%d && %s<=%d)", htBin.first, htBin.second, varName.c_str(), bins[i].first, varName.c_str(), bins[i].second) );
    } else {
      tree->Project( h1_mt2->GetName(), "mt2", Form("weight*( nJets>=2 && ht>%f && ht<=%f && %s>=%d && %s<=%d)", htBin.first, htBin.second, varName.c_str(), bins[i].first, varName.c_str(), bins[i].second) );
      //      tree->Project( h1_mt2->GetName(), "mt2", Form("weight*(  Z_pt>200 && fabs(Z_mass-91.19)<20. && nJets>=2 && ht>%f && ht<=%f && %s>=%d && %s<=%d)", htBin.first, htBin.second, varName.c_str(), bins[i].first, varName.c_str(), bins[i].second) );
    }

    h1_mt2->SetLineColor(colors[i]);
    h1_mt2->SetLineWidth(2);

    h1_mt2->Scale(1./h1_mt2->Integral());
    h1_mt2->Draw("same");
    //h1_mt2->DrawNormalized("same");

    if( i==0 ) {
      for( int ibin=0; ibin<h1_mt2->GetXaxis()->GetNbins()+1; ++ibin ) {
        h1_ref->SetBinContent(ibin, h1_mt2->GetBinContent(ibin) );
        h1_ref->SetBinError(ibin, h1_mt2->GetBinError(ibin) );
      }
    } else {
      pad2->cd();
      TH1D* h1_ratio = new TH1D( Form("ratio_%s_ht%.0f_%d", varName.c_str(), htBin.first, i), "", numberOfBins, 200., 1700. );
      h1_ratio->Sumw2();
      for( int ibin=0; ibin<h1_ratio->GetXaxis()->GetNbins()+1; ++ibin ) {
        float num = h1_mt2->GetBinContent(ibin);
        float den = h1_ref->GetBinContent(ibin);
        float err_num = h1_mt2->GetBinError(ibin);
        float err_den = h1_ref->GetBinError(ibin);
        if( den>0. ) {
          h1_ratio->SetBinContent(ibin, num/den );
          h1_ratio->SetBinError(ibin, sqrt( err_num*err_num/(den*den) + err_den*err_den*num*num/(den*den*den*den) ));
        }
      }

      h1_ratio->SetLineColor(colors[i]);
      h1_ratio->SetLineWidth(2);
      h1_ratio->Draw("same");
      pad1->cd();
    }

    if( bins[i].first==bins[i].second )
      legend->AddEntry( h1_mt2, Form("%s = %d", varName.c_str(), bins[i].second), "L" );
    else if( bins[i].second<10000 )
      legend->AddEntry( h1_mt2, Form("%d #leq %s #leq %d", bins[i].first, varName.c_str(), bins[i].second), "L" );
    else
      legend->AddEntry( h1_mt2, Form("%s #geq %d", varName.c_str(), bins[i].first), "L" );

  }

  legend->Draw("same");

  c1->cd();
  pad1->Draw();
  pad2->Draw();

  MT2DrawTools::addLabels( c1, 40., "CMS Simulation" );
  gPad->RedrawAxis();
  if( drawInclusivePlot ) {
    c1->SaveAs( Form("%s/mt2_vs_%s_ht%.0f_all.eps", outdir.c_str(), varName.c_str(), htBin.first) );
    c1->SaveAs( Form("%s/mt2_vs_%s_ht%.0f_all.png", outdir.c_str(), varName.c_str(), htBin.first) );
    c1->SaveAs( Form("%s/mt2_vs_%s_ht%.0f_all.pdf", outdir.c_str(), varName.c_str(), htBin.first) );
  } else {
    c1->SaveAs( Form("%s/mt2_vs_%s_ht%.0f.eps", outdir.c_str(), varName.c_str(), htBin.first) );
    c1->SaveAs( Form("%s/mt2_vs_%s_ht%.0f.png", outdir.c_str(), varName.c_str(), htBin.first) );
    c1->SaveAs( Form("%s/mt2_vs_%s_ht%.0f.pdf", outdir.c_str(), varName.c_str(), htBin.first) );
  }

  delete c1;
  delete h2_axes;
  delete h2_axes_ratio;

}



void drawYield( const std::string& outdir, TTree* tree, const std::string& suffix, const std::string& cuts, float binWidth ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float xMin=200.;
  float xMax=1400.;
  int nBins = (int)((xMax-xMin)/binWidth);

  TH1D* h1_mt2 = new TH1D("mt2", "", nBins, xMin, xMax);
  h1_mt2->Sumw2();
  h1_mt2->SetXTitle( "M_{T2} [GeV]" );
  h1_mt2->SetYTitle( "Events" );

  //TH1D* h1_mt2_rebin = new TH1D("mt2_rebin", "", 6, xMin, xMax);
  //h1_mt2_rebin->Sumw2();

  tree->Project( h1_mt2->GetName(), "mt2", Form("40.*weight*(%s)", cuts.c_str()) );
  //tree->Project( h1_mt2_rebin->GetName(), "mt2", Form("40.*weight*(%s)", cuts.c_str()) );

  h1_mt2->SetLineColor(46);
  h1_mt2->SetLineWidth(2);

  //h1_mt2_rebin->SetLineColor(46);
  //h1_mt2_rebin->SetLineWidth(2);

  h1_mt2->Draw("");
  //h1_mt2_rebin->Draw("same");

  TLine* line20 = new TLine( xMin, 20., xMax, 20.);
  line20->SetLineStyle(2);
  line20->Draw("same");

  TLine* line10 = new TLine( xMin, 10., xMax, 10.);
  line10->SetLineStyle(2);
  line10->SetLineColor(kGray);
  line10->Draw("same");

  //TLegend* legend = new TLegend( 0.6, 0.75, 0.9, 0.9);
  //legend->SetFillColor(0);
  //legend->SetTextSize(0.035);
  ////legend->SetHeader(cuts.c_str());
  //legend->AddEntry( h1_mt2, "100 GeV bins", "L" );
  //legend->AddEntry( h1_mt2_rebin, "200 GeV bins", "L" );
  //legend->Draw("same");

  MT2DrawTools::addLabels( c1, 40., "CMS Simulation" );
  gPad->RedrawAxis();
  c1->SaveAs(Form("%s/yield_%s.eps", outdir.c_str(), suffix.c_str()));
  c1->SaveAs(Form("%s/yield_%s.pdf", outdir.c_str(), suffix.c_str()));

  c1->SetLogy();

  h1_mt2->Draw("");

  line20->Draw("same");
  line10->Draw("same");

  MT2DrawTools::addLabels( c1, 40., "CMS Simulation" );
  gPad->RedrawAxis();
  c1->SaveAs(Form("%s/yield_%s_log.eps", outdir.c_str(), suffix.c_str()));
  c1->SaveAs(Form("%s/yield_%s_log.pdf", outdir.c_str(), suffix.c_str()));

  
  delete c1;
  delete h1_mt2;
  //delete h1_mt2_rebin;
  
}
