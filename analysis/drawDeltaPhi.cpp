#include <iostream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "THStack.h"
#include "TGraphErrors.h"

#include "../interface/MT2DrawTools.h"






int main( int argc, char* argv[] ) {


  MT2DrawTools::setStyle();


  TFile* file_zinv = TFile::Open("EventYields_data_Run2015D_25nsGolden_miniAODv2/analyses.root");
  TTree* tree_zinv = (TTree*)file_zinv->Get("ZJets/HT200toInf_j1toInf_b0toInf/tree_ZJets_HT200toInf_j1toInf_b0toInf");

  TFile* file_gjet = TFile::Open("EventYields_data_Run2015D_25nsGolden_miniAODv2/gammaControlRegion/mc.root");
  TTree* tree_gjet = (TTree*)file_gjet->Get("gammaCRtree/HT200toInf_j1toInf_b0toInf/tree_gammaCRtree_HT200toInf_j1toInf_b0toInf");


  std::vector<float> htCuts;
  htCuts.push_back( 200. );
  htCuts.push_back( 450. );
  htCuts.push_back( 1000. );

  for( unsigned i=0; i<htCuts.size(); ++i ) {


    TH1D* h1_zinv   = new TH1D( "zinv"  , "", 32, 0., 3.2 );
    TH1D* h1_prompt = new TH1D( "prompt", "", 32, 0., 3.2 );
    TH1D* h1_frag   = new TH1D( "frag"  , "", 32, 0., 3.2 );
    TH1D* h1_fake   = new TH1D( "fake"  , "", 32, 0., 3.2 );

    h1_zinv   ->Sumw2();
    h1_prompt ->Sumw2();
    h1_frag   ->Sumw2();
    h1_fake   ->Sumw2();

    tree_zinv->Project( "zinv"  , "deltaPhiMin", Form("weight*(ht>%f && mt2>200. && met>200. && nJets>1)", htCuts[i]) );
    tree_gjet->Project( "prompt", "deltaPhiMin", Form("weight*(ht>%f && mt2>200. && met>200. && nJets>1 && prompt==2)", htCuts[i] ));
    tree_gjet->Project( "frag"  , "deltaPhiMin", Form("weight*(ht>%f && mt2>200. && met>200. && nJets>1 && prompt==1)", htCuts[i] ));
    tree_gjet->Project( "fake"  , "deltaPhiMin", Form("weight*(ht>%f && mt2>200. && met>200. && nJets>1 && prompt==0)", htCuts[i] ));

    float zinv_int = h1_zinv->Integral();
    float prompt_int = h1_prompt->Integral();
    float frag_int = h1_frag->Integral();
    float fake_int = h1_fake->Integral();
    float gjet_int = prompt_int + frag_int + fake_int;

    h1_zinv->Scale( 1./zinv_int );  
    h1_prompt->Scale( 1./gjet_int );  
    h1_frag->Scale( 1./gjet_int );  
    h1_fake->Scale( 1./gjet_int );  

    h1_prompt->SetFillColor(18);
    h1_frag  ->SetFillColor(38);
    h1_fake  ->SetFillColor(46);


    THStack* stack = new THStack(0);
    stack->Add( h1_fake );
    stack->Add( h1_frag );
    stack->Add( h1_prompt );

    


    TCanvas* c1 = new TCanvas("c1", "", 600., 600.);
    c1->cd();

    float yMax = h1_zinv->GetMaximum()/h1_zinv->Integral();

    TH2D* h2_axes = new TH2D("axes", "", 10, 0., 3.2, 10, 0., 1.2*yMax);
    h2_axes->SetXTitle( "Minimum #Delta#Phi (ME_{T}, jets)" );
    h2_axes->SetYTitle( "Normalized to Unity" );

    h2_axes->Draw("");

    h1_zinv->SetMarkerStyle(20);
    h1_zinv->SetMarkerColor(kBlack);
    h1_zinv->SetMarkerSize(1.2);

    float xMin_legend = (htCuts[i]<700.) ? 0.2 : 0.65;
    float xMax_legend = (htCuts[i]<700.) ? 0.55 : 0.9;

    TLegend* legend = new TLegend( xMin_legend, 0.6, xMax_legend, 0.9, Form("H_{T} > %.0f GeV", htCuts[i]) );
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry( h1_zinv, "Z#rightarrow#nu#nu MC", "P" );
    legend->AddEntry( h1_prompt, "Prompt #gamma", "F" );
    legend->AddEntry( h1_frag, "Fragm. #gamma", "F" );
    legend->AddEntry( h1_fake, "Fakes", "F" );
    legend->Draw("same");

    stack->Draw("histo same");
    h1_zinv->Draw("p same");

    TPaveText* labelTop = MT2DrawTools::getLabelTopSimulation();
    labelTop->Draw("same");

    gPad->RedrawAxis();

    c1->SaveAs(Form("deltaPhiMin_ht%.0f.eps", htCuts[i]) );
    c1->SaveAs(Form("deltaPhiMin_ht%.0f.pdf", htCuts[i]) );

    delete c1;
    delete h2_axes;
    delete stack;
    delete h1_prompt;
    delete h1_frag;
    delete h1_fake;
    delete h1_zinv;

  } // for ht cuts

  return 0;

}
