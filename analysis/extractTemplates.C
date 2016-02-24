#include "TPaveStats.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFitResultPtr.h" 
#include "TFitResult.h" 

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooAddPdf.h"
#include "rooDoubleCB.h"

using namespace RooFit;
using namespace std;


#include <sstream>


int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}


void extractTemplates()
{ 

  TString fileName;
  fileName="/shome/casal/templatesRandS/QCD_13TeV_MGMLM_Spring15_bestMatching_angles_withNeutrinos.root";


   



  //=========  settings ====================
  //gROOT->LoadMacro("~/tdrstyle.C");
  //setTDRStyle();


 //gROOT->SetStyle("Plain");
 gStyle->SetPadGridX(kTRUE);
 gStyle->SetPadGridY(kTRUE);
 gStyle->SetPadRightMargin(0.07);
 gStyle->SetPadLeftMargin(0.13);
 //gStyle->SetOptStat("eour");
 gStyle->SetOptStat("");
 //gStyle->SetTitleXSize(0.07); 
 //gStyle->SetTitleXOffset(0.6); 
 //tyle->SetTitleYSize(0.3);
 //gStyle->SetLabelSize(0.6) 
 //gStyle->SetTextSize(0.5);
  gStyle->SetPalette(1);
  


 //=============================================


 TCanvas *canvas;
 canvas = new TCanvas("Prova","Prova",500,500);
 canvas->cd();


 //==============

 string output;


 
 TFile * file = new TFile(fileName); 
 TIter next(file->GetListOfKeys());
 TKey* key;
 string compString = "TH1F";
 int counter=0;
 while ( (key=(TKey*)next()) )
   {
     counter ++;
     if(key->GetClassName() == compString){ 
       //cout << "class name: " << key->GetClassName() << endl;
       TH1F* h = (TH1F*) key->ReadObj();
       //string histoName = h->GetName() ;
       //if(histoName.find("ResponsePt")==string::npos) continue;
       //if(histoName.find("h_tot")==string::npos) continue;
       TString histoName = h->GetName() ;
       if(!histoName.Contains("ResponsePt")) continue;
       if(!histoName.Contains("h_tot")) continue;
       cout << "name: " << histoName << endl;

       bool doFit=true;
       if(doFit){
	 double nentries = h->GetEntries();
	 double rms = h->GetRMS();
	 double mean = h->GetMean();
	 double leftBound,rightBound;
	 
	 leftBound = mean-1.0*rms;
	 rightBound = mean+1.0*rms;


	 TF1* f1 = new TF1("f1","gaus",leftBound,rightBound);
	 if(!(nentries>3)) {
	   h->Draw(); gPad->Update(); 
	   stringstream stream; stream << h->GetName() << ".pdf";
	   string outputfileName = "output/"+stream.str();
	   canvas->Print(outputfileName.c_str());
	   continue;
	 }


	 TFitResultPtr r = h->Fit(f1,"SMRL");   
	 h->Draw(); gPad->Update(); 
	 //Wait();
	 
	 double tmpMean = r->Parameter(1);
	 double tmpSigma = r->Parameter(2);    
	 
	 double xMin,xMax;
	 //xMin = 0;
	 //xMax = 3.;
	 xMin = TMath::Min(0.0,tmpMean - tmpSigma*10.);
	 xMax = TMath::Max(3.0,tmpMean + tmpSigma*10.);
	 RooRealVar x("x","x",xMin,xMax) ;

	 double meanRangeMin;
	 double meanRangeMax;
	 if(tmpMean<0){
	   meanRangeMin = 3*tmpMean;
	   meanRangeMax = 0.1*tmpMean;
	 }else{
	   meanRangeMin = 0.1*tmpMean;
	   meanRangeMax = 3*tmpMean;
	 }

	 RooRealVar meanRoo("mean","mean of gaussian",tmpMean,meanRangeMin,meanRangeMax) ;
	 RooRealVar sigmaRoo("sigma","width of gaussian",tmpSigma,tmpSigma*0.5,tmpSigma*1.5); 
	 
	 RooRealVar a("a","a",3.,1.,20.);
	 RooRealVar aDx("aDx","aDx",3.,0.5,20.);
	 RooRealVar n("n","n",5.,0.,50.);   
	 RooRealVar nDx("nDx","nDx",5.,0.,50.);   
	 
	 RooDoubleCB func1("cb","cb PDF",x,meanRoo,sigmaRoo,a,n,aDx,nDx) ;
	 RooPlot* xframe = x.frame(Title("CB p.d.f.")) ;
    

	 RooDataHist  dh("dh","dh",x,Import(*h));

	 func1.fitTo(dh,NumCPU(4));

	 RooWorkspace *w = new RooWorkspace("wsp_"+histoName,"workspace "+histoName) ;
	 w->import(func1);
	 w->importClassCode("RooDoubleCB",kTRUE);
	 w->writeToFile("responseTemplates.root",kFALSE);

	 dh.plotOn(xframe);
	 func1.plotOn(xframe) ;
	 //
	 //xframe->GetXaxis()->SetRangeUser(xMin*10,xMax*10);
	 xframe->Draw(); gPad->Update();
	 //Wait();

	 string outputfileName;

	 gPad->SetLogy(0);
	 stringstream stream1; stream1 << h->GetName() << "_lin.pdf";
	 stringstream stream2; stream2 << h->GetName() << "_lin.png";
	 outputfileName= "output/"+stream1.str();  canvas->Print(outputfileName.c_str());
	 outputfileName= "output/"+stream2.str();  canvas->Print(outputfileName.c_str());

	 
	 gPad->SetLogy(1);
	 stringstream stream1log; stream1log << h->GetName() << "_log.pdf";
	 stringstream stream2log; stream2log << h->GetName() << "_log.png";
	 outputfileName= "output/"+stream1log.str();  canvas->Print(outputfileName.c_str());
	 outputfileName= "output/"+stream2log.str();  canvas->Print(outputfileName.c_str());


       }
     }
     //if(counter==100) break;
   }




 //delete canvas;
}



// ----------------------------------




