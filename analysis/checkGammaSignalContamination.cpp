#include <iostream>
#include "interface/MT2Region.h"
#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"



int main() {

  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this


  std::string gammaCR_file = "GammaControlRegion_PHYS14_v5_skimprune_zurich_4fb/mc.root";

  MT2Analysis<MT2Estimate>* prompt = MT2Analysis<MT2Estimate>::readFromFile(gammaCR_file, "prompt_pass");

  std::set<MT2Region> regions = prompt->getRegions();

  std::vector< MT2Analysis<MT2Estimate>* > signals = MT2Analysis<MT2Estimate>::readAllFromFile(gammaCR_file, "SMS");

  for( unsigned i=0; i<signals.size(); ++i ) {

    std::cout << std::endl;
    std::cout << signals[i]->getName() << std::endl; 

    for( std::set<MT2Region>::iterator iR = regions.begin(); iR != regions.end(); ++iR ) {

      if( iR->nBJetsMin()>1 ) continue;

      MT2Estimate* this_prompt = prompt->get(*iR);
      MT2Estimate* this_signal = signals[i]->get(*iR);

      float yield_prompt = this_prompt->yield->Integral();
      float yield_signal = this_signal->yield->Integral();
      float contamination100 =  yield_signal/yield_prompt*100.;

      if( contamination100<0.1 ) continue;

      std::cout << "    " << iR->getName() << "  \t" << contamination100 << std::endl;

    }


  }


  return 0;

}

