#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>




int main( int argc, char* argv[] ) {



  if( argc!=2 ) {

    std::cout << "USAGE: ./printYield [EventYieldsDirectory]" << std::endl;
    exit(11);

  }


  std::string dir(argv[1]);

  std::string firstInputFile  = dir + "/analyses.root";

  std::vector< MT2Analysis<MT2Estimate>* > analyses_bg;
  analyses_bg.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "Top") );
  analyses_bg.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "WJets") );
  analyses_bg.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "ZJets") );
  //analyses_bg.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "QCD") );


  for(unsigned a =0;  a < analyses_bg.size(); ++a){
    std::string ofs = dir + "/yieldtable_" + analyses_bg[a]->getName() + ".log";
    analyses_bg[a]->print(ofs);
  }


  MT2Analysis<MT2Estimate>* analysisSum = new MT2Analysis<MT2Estimate>( *(analyses_bg[0]) );
  for(unsigned a =1;  a < analyses_bg.size(); ++a)
    (*analysisSum) += (*(analyses_bg[a]));

  analysisSum->setName("FullBkg");

  analysisSum->print(dir + "/yieldtable_FullBkg.log");


  std::vector< MT2Analysis<MT2Estimate>* > estimates;
  estimates.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "zinv") );
  estimates.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "llep") );
  for(unsigned a =0;  a < estimates.size(); ++a){
    if( estimates[a]==0 ) {
      std::cout << "-> Looks like you didn't create datacards yet for dir '" << dir << "'. Will not print the estimates." << std::endl;
      break;
    }
    std::string ofs = dir + "/yieldtable_" + estimates[a]->getName() + ".log";
    estimates[a]->print(ofs);
  }
  
  //MT2Analysis<MT2Estimate>* analysisTop = MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "LostLepton");


  std::vector < MT2Analysis<MT2Estimate>* > analysesSignal = MT2Analysis<MT2Estimate>::readAllFromFile(firstInputFile.c_str(), "SMS");
  
  for(unsigned a =0;  a < analysesSignal.size(); ++a) {
    std::string ofs = dir + "/yieldtable_" + analysesSignal[a]->getName() + ".log";
    analysesSignal[a]->print(ofs);
  }


  return 0;

}
