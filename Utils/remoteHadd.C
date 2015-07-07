/*
How to run:
root -l
.L remoteHadd.C+
hadd("/mmasciov/CSA14_production/XXX/", "_PU_S14_POSTLS170", "/scratch/mmasciov/ToPostProcess/", "_PU_S14_POSTLS170_babytree.root")
*/

#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"

using namespace std;

int hadd(string inputPath = "/mmasciov/CSA14v2_production/06Feb2015_csa14v2/",
	 string outputPath = "/scratch/mmasciov/CSA14v2/ToPostProcess/",
	 string fileExtension = "_PU_S14_POSTLS170_babytree.root"){
  
  int start_s=clock();

  system(Form("mkdir -p %s", outputPath.c_str()));

  system(Form("ls /pnfs/psi.ch/cms/trivcat/store/user/%s > listOfSamples_hadd.txt", inputPath.c_str()));

  std::cout << "Starting..." << std:: endl;
  
  ifstream haddSampleIN("listOfSamples_hadd.txt");
  char haddBuffer[300];
  
  int sampleCounter(0);
  while( haddSampleIN.getline(haddBuffer, 300, '\n') ) {

    if (haddBuffer[0] == '\n' || haddBuffer[0] == ' ' || haddBuffer[0] == '\0' || strlen(haddBuffer)==1 ) /*|| strlen(haddBuffer) < 18)*/
      continue;

    ++sampleCounter;

    char sampleToHadd[300];
    sscanf(haddBuffer, "%s\n", sampleToHadd);
    std::cout << sampleToHadd << std::endl;

    std::string outputFile = outputPath + sampleToHadd + fileExtension;
    FILE* rootfile;
    rootfile = fopen(outputFile.c_str(), "r");
    if (rootfile == NULL);
    else {

      std::cout << "ATTENTION: Target file already exists!!! I will skip this sample..." << std::endl;
      continue;
    
    }

    system(Form("ls /pnfs/psi.ch/cms/trivcat/store/user/%s%s > listOfFiles_preHadd.txt", inputPath.c_str(), sampleToHadd));
    
    std::string fullInput_hadd = "dcap://t3se01.psi.ch:22125/" + std::string(sampleToHadd);
    //std::string fullInput_hadd = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/" + inputPath + sampleToHadd;

    ifstream haddFileIN("listOfFiles_preHadd.txt");
    char haddFile[300];

    int fileCounter(0);
    std::ofstream listToHadd;
    listToHadd.open("listOfFiles_hadd.txt");

    while( haddFileIN.getline(haddFile, 300, '\n') ){

      if (haddFile[0] == '\n' || haddFile[0] == ' ' || haddFile[0] == '\0')
	continue;
      
      ++fileCounter;
      
      char fileToHadd[300];
      sscanf(haddFile, "%s\n", fileToHadd);
      
      listToHadd << fullInput_hadd.c_str() << "/" << fileToHadd << std::endl;

    }

    system(Form("hadd %s @listOfFiles_hadd.txt", outputFile.c_str()));
   
    std::cout << "# files hadd-ed for sample " << sampleToHadd << ": " << fileCounter << std::endl;

  }
  
  std::cout << "# of samples hadd-ed: " << sampleCounter;

  int stop_s = clock();
  std::cout << std::endl << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;
  
  system("rm -f listOfSamples_hadd.txt");
  system("rm -f listOfFiles_preHadd.txt");
  system("rm -f listOfFiles_hadd.txt");

  return 0;
}
