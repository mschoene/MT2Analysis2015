#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

#include "interface/MT2Sample.h"
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateTree.h"
#include "interface/MT2DrawTools.h"

#define mt2_cxx
#include "../interface/mt2.h"

#include "TEventList.h"



//power = TF1("power",[[0]*TMath::Power(x,[1])) 
//lowHT: p0 = 73792.1; p1 = -2.53584 
//medHT: p0 = 57785.8; p1 = -2.49123 
//hiHT : p0 = 1958.03; p1 = -1.70589 
//extHT: p0 = 385.588; p1 = -1.30116


Bool_t doInclusive=false;

MT2Analysis<MT2Estimate>* estimateYield( const std::string& regionsSetSR, MT2Analysis<MT2EstimateTree>* analysisCR);

TF1* getPower(float ht);
void getRhat(int nj, float &r0, float &r1, float &r2, float &r3);
void getFnj (float ht, float &f23, float &f46, float &f7);

int main( int argc, char* argv[] ) {

  std::string regionsSetCR = "zurich_onlyHT";
  std::string regionsSetSR = "zurich";
  MT2Analysis<MT2EstimateTree>* qcdCR = MT2Analysis<MT2EstimateTree>::readFromFile("qcdCRtree.root", "QCD");
  qcdCR->setName("qcdCR");
  MT2Analysis<MT2EstimateTree>* qcdSR = MT2Analysis<MT2EstimateTree>::readFromFile("qcdSRtree.root", "QCD");
  qcdSR->setName("qcdSR");

  MT2Analysis<MT2Estimate>* qcdSRpred;
  if (doInclusive)
    qcdSRpred = estimateYield("zurich_onlyHT", qcdCR);
  else
    qcdSRpred = estimateYield("zurich", qcdCR);

  qcdSRpred->writeToFile(doInclusive ? "qcdINCpred.root" : "qcdSRpred.root");


}

MT2Analysis<MT2Estimate>* estimateYield( const std::string& regionsSetSR, MT2Analysis<MT2EstimateTree>* analysisCR ) {

  MT2Analysis<MT2Estimate>* analysisSR =  new MT2Analysis<MT2Estimate>( "estimateSR", regionsSetSR );

  std::set<MT2Region> MT2Regions = analysisCR->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree = analysisCR->get(thisRegion)->tree;
    int nentries = tree->GetEntries();

    float mt2;
    float weight;

    tree->SetBranchAddress("mt2"    , &mt2);
    tree->SetBranchAddress("weight" , &weight);
  
    float ht = thisRegion.htMin()+1.0;

    TF1* power = getPower(ht);
    float f23, f46, f7;
    getFnj (ht, f23, f46, f7);

    
    for( int iEntry=0; iEntry<nentries; ++iEntry ) {
      if( iEntry % 25000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
      tree->GetEntry(iEntry);

      if(mt2<200) continue;


      float rFit = power->Eval(mt2);

      if (doInclusive)
	analysisSR->get( ht, 4, 0,  -1, mt2 )->yield->Fill(mt2,rFit*weight);  // doesnt matter where to fill, it's inclusive
      else{
	float r0, r1, r2, r3;
	getRhat (2, r0,r1,r2,r3);
	analysisSR->get( ht, 2, 0,  -1, mt2 )->yield->Fill(mt2,  f23     *r0*rFit*weight);
	analysisSR->get( ht, 2, 1,  -1, mt2 )->yield->Fill(mt2,  f23     *r1*rFit*weight);
	analysisSR->get( ht, 2, 2,  -1, mt2 )->yield->Fill(mt2,  f23     *r2*rFit*weight);
	
	getRhat (4, r0,r1,r2,r3);
	analysisSR->get( ht, 4, 0,  -1, mt2 )->yield->Fill(mt2,  f46     *r0*rFit*weight);
	analysisSR->get( ht, 4, 1,  -1, mt2 )->yield->Fill(mt2,  f46     *r1*rFit*weight);
	analysisSR->get( ht, 4, 2,  -1, mt2 )->yield->Fill(mt2,  f46     *r2*rFit*weight);
	analysisSR->get( ht, 4, 3,  -1, mt2 )->yield->Fill(mt2, (f23+f46)*r3*rFit*weight);
	
	getRhat (7, r0,r1,r2,r3);
	analysisSR->get( ht, 7, 0,  -1, mt2 )->yield->Fill(mt2,  f7      *r0*rFit*weight);
	analysisSR->get( ht, 7, 1,  -1, mt2 )->yield->Fill(mt2,  f7      *r1*rFit*weight);
	analysisSR->get( ht, 7, 2,  -1, mt2 )->yield->Fill(mt2,  f7      *r2*rFit*weight);
	analysisSR->get( ht, 7, 3,  -1, mt2 )->yield->Fill(mt2,  f7      *r3*rFit*weight);
      }	
    }
    delete tree;
    delete power;
  }
  analysisSR->finalize();

  return analysisSR;
}


TF1* getPower(float ht){
  TF1* power = new TF1("power","[0]*TMath::Power(x,[1])");
  if (ht<575){
    power->SetParameter(0, 73792.1);
    power->SetParameter(1,-2.53584);
  }
  else if (ht<1000){
    power->SetParameter(0, 57785.8);
    power->SetParameter(1,-2.49123);
  }
  else if (ht<1500){
    power->SetParameter(0, 1958.03);
    power->SetParameter(1,-1.70589);
  }
  else {
    power->SetParameter(0, 385.588);
    power->SetParameter(1,-1.30116);
  }

  return power;
}

void getRhat(int nj, float &r0, float &r1, float &r2, float &r3){
// fit from 80 GeV
//r_0b = {'h450toInf_j2to3': 0.81  , 'h450toInf_j4to6': 0.74  , 'h450toInf_j7toInf': 0.64 }
//r_1b = {'h450toInf_j2to3': 0.16  , 'h450toInf_j4to6': 0.20  , 'h450toInf_j7toInf': 0.24 }
//r_2b = {'h450toInf_j2to3': 0.027 , 'h450toInf_j4to6': 0.047 , 'h450toInf_j7toInf': 0.073}
//r_3b = {'h450toInf_j2to3': 0.0016, 'h450toInf_j4to6': 0.0062, 'h450toInf_j7toInf': 0.014}
//  if (nj<4){
//    r0 = 0.81;
//    r1 = 0.16;
//    r2 = 0.027;
//    r3 = 0.0016;
//  }
//  else if(nj<7){
//    r0 = 0.74;
//    r1 = 0.20;
//    r2 = 0.047;
//    r3 = 0.0062;
//  }
//  else{
//    r0 = 0.64;
//    r1 = 0.24;
//    r2 = 0.073;
//    r3 = 0.0014;
//  }

// fit from 100 GeV
//r_0b = {'h450toInf_j2to3': 0.81  , 'h450toInf_j4to6': 0.73  , 'h450toInf_j7toInf': 0.62 }
//r_1b = {'h450toInf_j2to3': 0.15  , 'h450toInf_j4to6': 0.21  , 'h450toInf_j7toInf': 0.27 }
//r_2b = {'h450toInf_j2to3': 0.027 , 'h450toInf_j4to6': 0.051 , 'h450toInf_j7toInf': 0.085}
//r_3b = {'h450toInf_j2to3': 0.0014, 'h450toInf_j4to6': 0.0060, 'h450toInf_j7toInf': 0.013}
//  if (nj<4){
//    r0 = 0.81;
//    r1 = 0.15;
//    r2 = 0.027;
//    r3 = 0.0014;
//  }
//  else if(nj<7){
//    r0 = 0.73;
//    r1 = 0.21;
//    r2 = 0.051;
//    r3 = 0.0060;
//  }
//  else{
//    r0 = 0.62;
//    r1 = 0.27;
//    r2 = 0.085;
//    r3 = 0.0013;
//  }


// with cut on deltaPhi and on diffmetmht and fit from mt2>200
//r_0b = {'h450toInf_j2to3': 0.781 , 'h450toInf_j4to6': 0.578, 'h450toInf_j7toInf': 0.322}
//r_1b = {'h450toInf_j2to3': 0.162 , 'h450toInf_j4to6': 0.282, 'h450toInf_j7toInf': 0.365}
//r_2b = {'h450toInf_j2to3': 0.037 , 'h450toInf_j4to6': 0.118, 'h450toInf_j7toInf': 0.250}
//r_3b = {'h450toInf_j2to3': 0.0032, 'h450toInf_j4to6': 0.021, 'h450toInf_j7toInf': 0.060}

//high ht
  if (nj<4){
    r0 = 0.772;
    r1 = 0.202;
    r2 = 0.026;
    r3 = 0.002;
  }
  else if(nj<7){
    r0 = 0.534;
    r1 = 0.289;
    r2 = 0.129;
    r3 = 0.026;
  }
  else{
    r0 = 0.359;
    r1 = 0.335;
    r2 = 0.236;
    r3 = 0.057;
  }
//inclusive ht
//  if (nj<4){
//    r0 = 0.781;
//    r1 = 0.162;
//    r2 = 0.037;
//    r3 = 0.003;
//  }
//  else if(nj<7){
//    r0 = 0.578;
//    r1 = 0.282;
//    r2 = 0.118;
//    r3 = 0.021;
//  }
//  else{
//    r0 = 0.322;
//    r1 = 0.365;
//    r2 = 0.250;
//    r3 = 0.060;
//  }

}

void  getFnj (float ht, float &f23, float &f46, float &f7){
  // from MT2>50
//f(23) = 7.391954e-01 +- 2.767147e-01
//f(46) = 2.608047e-01 +- 4.231398e-02
//f(7)   = 0.0 +- 0.0
//Medium HT:
//f(23)  = 3.163028e-01 +- 7.205447e-02
//f(46)  = 6.551170e-01 +- 2.373133e-01
//f(7)    = 2.858040e-02 +- 3.405140e-04
//High HT:
//f(23)  = 1.187736e-01 +- 2.216757e-04
//f(46)  = 6.424398e-01 +- 3.592727e-03
//f(7)    = 2.387865e-01 +- 7.733932e-04
//Extreme HT:
//f(23)  = 1.552706e-01 +- 6.294285e-03 
//f(46)  = 5.661131e-01 +- 1.307790e-02
//f(7)    = 2.786162e-01 +- 3.244173e-03

//  // from MT2>100
//Low HT
//f(23)=2.142857e-01
//f(46)=7.857146e-01
//f(7)=0
//
//Medium HT
//f(23)=1.596447e-01
//f(46)=6.134987e-01
//f(7)=2.268566e-01
//
//High HT
//f(23)=7.583390e-02
//f(46)=6.413972e-01
//f(7)=2.855961e-01
//
//Extreme HT
//f(23)=8.151883e-02
//f(46)=5.805375e-01
//f(7)=3.451356e-01

//  // from MT2>200

  if (ht<575){
    f23 = 0.30;
    f46 = 0.70;
    f7  = 0.01;
  }
  else if (ht<1000){
    f23 = 0.14;
    f46 = 0.50;
    f7  = 0.36;
  }
  else if (ht<1500){
    f23 = 0.18;
    f46 = 0.56;
    f7  = 0.26; 
  }
  else {
    f23 = 0.12;
    f46 = 0.59;
    f7  = 0.29;
  }

}
