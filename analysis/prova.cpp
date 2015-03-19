#include <iostream>

#include "interface/MT2Analysis.h"
#include "interface/MT2EstimateSyst.h"



int main() {

  system("rm prova.root prova2.root");

  // create a signal region structure with the 8TeV definitions of signal regions:
  // "pippo" is a name that will figure in the histogram names
  // "8TeV" is the name of a given signal region set (for now i've implemented the one used at 8 tev and the recent proposal by mario - called "13TeV")
  // 702 is an optional argument: it's the id of the sample (so 702 in this case is DYJetsToLL_M-50_HT-100to200)
  MT2Analysis<MT2EstimateSyst> *analysis = new MT2Analysis<MT2EstimateSyst>("pippo", "13TeV_CSA14", 702);
  // ^ the above command creates the whole HT x SR structure, and instantiates one MT2EstimateSyst per region

  float ht = 500.;
  float met = 300.;
  int njets = 2;
  int nbjets = 0;
  float mt2 = 305.;

  float weight = 1.;
  float weight_btagUp   = 1.05;
  float weight_btagDown = 0.95;

  MT2EstimateSyst* est = analysis->get(ht,njets,nbjets,met);
  // ^ the above command looks for the region identified by ht,njets,nbjets,met, and gets the corresponding pointer to the MT2EstimateSyst
  // it will return a null pointer if the event doesnt fall in any region (eg if it has too low HT)
  if( est!=0 ) {
    // if you found a valid region, you can use the pointer, for example filling histograms
    std::cout << "filling region: " << est->region->getName() << std::endl;
    est->yield->Fill( mt2, weight );
    est->yield_systUp->Fill( mt2, weight_btagUp );
    est->yield_systDown->Fill( mt2, weight_btagDown );
  }

  // the command below calls the method addOverflow() for all elements in the MT2Analysis
  // so make sure that if you implement a new class, it has that method implemented
  analysis->finalize();

  // write the histograms to a file:
  analysis->writeToFile("prova.root");
  // ^ the above will call the method write() for all objects, so again make sure it's implemented

  MT2Analysis<MT2EstimateSyst> *analysis2 = new MT2Analysis<MT2EstimateSyst>("pippo2", "13TeV_CSA14", 702);
  *analysis2 = *analysis;

  analysis2->writeToFile("prova2.root");

  analysis2->setName("paperino");
  analysis2->writeToFile("paperino.root");


  return 0;

}

