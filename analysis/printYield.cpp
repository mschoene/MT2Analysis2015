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
  analyses_bg[0]->setName("Top");

  analyses_bg.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "ZJets") );
  analyses_bg[1]->setName("Z+jets");

  analyses_bg.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "WJets") );
  analyses_bg[2]->setName("W+jets");

  analyses_bg.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "QCD") );
  analyses_bg[3]->setName("QCD");
  
  MT2Analysis<MT2Estimate>* analysisSum = new MT2Analysis<MT2Estimate>( *(analyses_bg[0]) );
  for(unsigned a =1;  a < analyses_bg.size(); ++a)
    (*analysisSum) += (*(analyses_bg[a]));
  analysisSum->setName("Total SM");

  std::vector< MT2Analysis<MT2Estimate>* > estimates;
  estimates.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "zinv") );
  estimates[0]->setName("Z$\rightarrow\nu\nu$");
  
  estimates.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "llep") );
  estimates[1]->setName("Lost lepton");

  //std::vector < MT2Analysis<MT2Estimate>* > analysesSignal = MT2Analysis<MT2Estimate>::readAllFromFile(firstInputFile.c_str(), "SMS");
  std::vector < MT2Analysis<MT2Estimate>* > analysesSignal;
  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T1tttt_mGluino1500_mLSP100") );
  analysesSignal[0]->setName("T1tttt 1500,100");
  
  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T1tttt_mGluino1200_mLSP800") );
  analysesSignal[1]->setName("T1tttt 1200,800");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T1bbbb_mGluino1500_mLSP100") );
  analysesSignal[2]->setName("T1bbbb 1500,100");
  
  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T1bbbb_mGluino1000_mLSP900") );
  analysesSignal[3]->setName("T1bbbb 1000,900");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T1qqqq_mGluino1400_mLSP100") );
  analysesSignal[4]->setName("T1qqqq 1400,100");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T1qqqq_mGluino1000_mLSP800") );
  analysesSignal[5]->setName("T1qqqq 1000,800");

//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T2tt_2J_mStop850_mLSP100") );
//  analysesSignal[6]->setName("T2tt 850,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T2tt_2J_mStop650_mLSP325") );
//  analysesSignal[7]->setName("T2tt 650,325");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T2tt_2J_mStop500_mLSP325") );
//  analysesSignal[8]->setName("T2tt 500,325");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T2tt_2J_mStop425_mLSP325") );
//  analysesSignal[9]->setName("T2tt 425,325");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T2bb_2J_mStop900_mLSP100") );
//  analysesSignal[10]->setName("T2bb 900,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T2bb_2J_mStop600_mLSP580") );
//  analysesSignal[11]->setName("T2bb 600,580");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T2qq_2J_mStop1200_mLSP100") );
//  analysesSignal[12]->setName("T2qq 1200,100");
//
//  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile(firstInputFile.c_str(), "SMS_T2qq_2J_mStop600_mLSP550") );
//  analysesSignal[13]->setName("T2qq 600,550");
  
  std::set<MT2HTRegion> htRegions = analyses_bg[0]->getHTRegions();
  std::set<MT2SignalRegion> sigRegions = analyses_bg[0]->getSignalRegions();
  
  std::set<MT2Region> mt2Regions = analyses_bg[0]->getRegions();
  
  MT2HTRegion* oldHTRegion = 0;
  for( std::set<MT2HTRegion>::iterator iHT=htRegions.begin(); iHT!=htRegions.end(); ++iHT )  {
    
    if( oldHTRegion && *(iHT) == *(oldHTRegion) ) continue;
    oldHTRegion = new MT2HTRegion( *(iHT) );
    
    std::cout << "Printing table for region: " << iHT->getName() << std::endl;

    std::string ofs = dir + "/yieldtable_full_" + iHT->getName() + ".tex";
    
    std::ofstream ofs_file;
    if ( std::ifstream(ofs) )
      system(Form("rm -f %s", ofs.c_str()));
    
    ofs_file.open( ofs, std::ofstream::app );
    
    ofs_file << "\\begin{table}[htbp]" << std::endl;
    ofs_file << "\\caption{Background estimate yields for " << iHT->getNiceName().c_str() << ".}" << std::endl;
    ofs_file << "\\centering" << std::endl;
    ofs_file << "\\begin{tabular}{r";
    for( unsigned int b=0; b<sigRegions.size(); ++b ){
      ofs_file << "|c";
    }
    ofs_file << "}" << std::endl;
    
    ofs_file << "\\hline" << std::endl;
    ofs_file << "\\multicolumn{" << sigRegions.size()+1 << "}{c}{" << iHT->getNiceName().c_str() << "}\\\\" << std::endl;
    ofs_file <<"\\hline \\hline" << std::endl;

    ofs_file << "Process";    
    for( std::set<MT2SignalRegion>::iterator iSig=sigRegions.begin(); iSig!=sigRegions.end(); ++iSig )
      ofs_file << " & " << iSig->getNiceName().c_str();
    ofs_file << "\\\\" << std::endl;    
    ofs_file <<"\\hline \\hline" << std::endl;

    MT2HTRegion* thisHTRegion = new MT2HTRegion( *(iHT) );
    
    for(unsigned a =0;  a < analyses_bg.size(); ++a){
      ofs_file << analyses_bg[a]->getName().c_str();
      analyses_bg[a]->print(ofs_file, thisHTRegion );
    }
    
    ofs_file <<"\\hline" << std::endl;
    ofs_file << analysisSum->getName().c_str();
    analysisSum->print(ofs_file, thisHTRegion );
    ofs_file <<"\\hline" << std::endl;


//    for(unsigned a =0;  a < estimates.size(); ++a){
//      if( estimates[a]==0 ) {
//	std::cout << "-> Looks like you didn't create datacards yet for dir '" << dir << "'. Will not print the estimates." << std::endl;
//	break;
//      }
//      ofs_file << estimates[a]->getName.c_str();
//      estimates[a]->print(ofs_file, thisHTRegion );    
//    }
  
    
    for(unsigned a =0;  a < analysesSignal.size(); ++a) {
      ofs_file << analysesSignal[a]->getName().c_str();
      analysesSignal[a]->print(ofs_file, thisHTRegion );
    }
    
    ofs_file << "\\hline" << std::endl;
    ofs_file << "\\end{tabular}" << std::endl;
    ofs_file << "\\end{table}" << std::endl;
    ofs_file << std::endl;
  
  } //loop over HT regions

  
  for( std::set<MT2Region>::iterator iMT2=mt2Regions.begin(); iMT2!=mt2Regions.end(); ++iMT2 )  {

    std::cout << "Printing table for region: " << iMT2->getName() << std::endl;
    
    std::string ofs_tr = dir + "/yieldtable_full_" + iMT2->getName() + ".tex";
    std::ofstream ofs_file_tr;
    if ( std::ifstream(ofs_tr) )
      system(Form("rm -f %s", ofs_tr.c_str()));
    
    ofs_file_tr.open( ofs_tr, std::ofstream::app );
    
    std::vector< std::string > names = iMT2->getNiceNamesLatex();
    
    MT2Region* thisRegion = new MT2Region( *(iMT2) );
    int nBins;
    double* bins;
    thisRegion->getBins(nBins, bins);

    ofs_file_tr << "\\begin{table}[htbp]" << std::endl;
    ofs_file_tr << "\\caption{Background estimate yields for " << names[0].c_str() << ", " << names[1].c_str() << ".}" << std::endl;
    ofs_file_tr << "\\centering" << std::endl;
    ofs_file_tr << "\\begin{tabular}{r";
    for( int b=0; b<nBins; ++b ){
      ofs_file_tr << "|c";
    }
    ofs_file_tr << "}" << std::endl;

    ofs_file_tr << "\\hline" << std::endl;
    ofs_file_tr << "\\multicolumn{" << nBins+1 << "}{c}{" << names[0].c_str() << ", " << names[1].c_str() << "}\\\\" << std::endl;
    ofs_file_tr <<"\\hline \\hline" << std::endl;
    
    ofs_file_tr << "Process";

    std::vector< std::string > binNames = thisRegion->getBinNames();
    for( unsigned b=0; b < binNames.size(); ++b )
      ofs_file_tr << " & " << binNames[b].c_str();
    ofs_file_tr << "\\\\" << std::endl;
    ofs_file_tr <<"\\hline \\hline" << std::endl;

    for(unsigned a =0;  a < analyses_bg.size(); ++a){
      ofs_file_tr << analyses_bg[a]->getName().c_str();
      analyses_bg[a]->print(ofs_file_tr, thisRegion );
    }
   
    ofs_file_tr <<"\\hline" << std::endl;
    ofs_file_tr << analysisSum->getName().c_str();
    analysisSum->print(ofs_file_tr, thisRegion );
    ofs_file_tr <<"\\hline" << std::endl;

//    for(unsigned a =0;  a < estimates.size(); ++a){
//      if( estimates[a]==0 ) {
//	std::cout << "-> Looks like you didn't create datacards yet for dir '" << dir << "'. Will not print the estimates." << std::endl;
//	break;
//      }
//      ofs_file_tr << estimates[a]->getName.c_str(); 
//      estimates[a]->print(ofs_file_tr, thisRegion );
//    }
  
    
    ofs_file_tr <<"\\hline" << std::endl;
    for(unsigned a =0;  a < analysesSignal.size(); ++a) {
      ofs_file_tr << analysesSignal[a]->getName().c_str();
      analysesSignal[a]->print(ofs_file_tr, thisRegion );
    }
    ofs_file_tr <<"\\hline" << std::endl;
    ofs_file_tr << "\\end{tabular}" << std::endl;
    ofs_file_tr << "\\end{table}" << std::endl;
    ofs_file_tr << std::endl;
    
  } //loop over TR's


  return 0;
  
}
