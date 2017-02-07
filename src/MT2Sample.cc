#include "../interface/MT2Sample.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"



MT2Sample::MT2Sample() {};



MT2Sample::~MT2Sample() {};




std::vector<MT2Sample> MT2Sample::loadSamples(const std::string& filename, int idMin, int idMax) {

  return loadSamples(filename, "", idMin, idMax);

}



std::vector<MT2Sample> MT2Sample::loadSamples(const std::string& filename, const std::string& filter, int idMin, int idMax) {


  std::vector<MT2Sample> fSamples;

  char buffer[200];
  std::ifstream IN(filename.c_str());


  std::cout << "------------------------------------" << std::endl;
  std::cout << "Sample File  " << filename << std::endl;

  int counter(0);
  TString fPath;


  while( IN.getline(buffer, 500, '\n') ) {

    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'                                                                                                                                                                                 
    }

    std::string sampleFilePath(buffer);
    TString sampleFilePath_tstr(sampleFilePath);
    if( !(sampleFilePath_tstr.EndsWith(".root")) ) continue;

    MT2Sample s;
    s.file = sampleFilePath;

    sampleFilePath_tstr.ReplaceAll( "/", " " );

    std::istringstream iss(sampleFilePath_tstr.Data());
    std::vector<std::string> pathDirs;
    do {
      std::string sub;
      iss >> sub;
      pathDirs.push_back(sub);
    } while (iss);

    TString rootFileName(pathDirs[pathDirs.size()-2]); // last element of vector always empty: take -2

    if( filter!="" && !(rootFileName.Contains(filter.c_str())) ) continue;

    rootFileName.ReplaceAll( "_", " " );
    rootFileName.ReplaceAll( ".", " " );
 
    std::vector<std::string> fileNameParts;
    std::istringstream iss2(rootFileName.Data());
    do {
      std::string sub;
      iss2 >> sub;
      fileNameParts.push_back(sub);
    } while (iss2);

    if( fileNameParts.size()==0 ) {
      std::cout << "[MT2Sample::loadSamples] ERROR! Illegal file name: " << rootFileName << std::endl;
      exit(1101);
    }
     

    TFile* file = TFile::Open(s.file.c_str());
    TTree* tree = (TTree*)file->Get("mt2");

    
    ULong64_t evt_nEvts;
    int evt_id;
    float evt_filter, evt_kfactor, evt_xsec, evt_scale1fb;

    tree->SetBranchAddress("evt_id",&evt_id);
    tree->SetBranchAddress("evt_nEvts",&evt_nEvts);
    tree->SetBranchAddress("evt_filter",&evt_filter);
    tree->SetBranchAddress("evt_kfactor",&evt_kfactor);
    tree->SetBranchAddress("evt_xsec",&evt_xsec);
    tree->SetBranchAddress("evt_scale1fb",&evt_scale1fb);
    
    tree->GetEntry(0);

    s.id       = evt_id;

    if(s.id<=0){

      if( rootFileName.Contains("JetHT") )               s.id = 1;
      else if( rootFileName.Contains("HTMHT") )          s.id = 2;
      else if( rootFileName.Contains("MET") )            s.id = 3;
      else if( rootFileName.Contains("DoubleEG") )       s.id = 4;
      else if( rootFileName.Contains("DoubleMuon") )     s.id = 5;
      else if( rootFileName.Contains("MuonEG") )         s.id = 6;
      else if( rootFileName.Contains("SinglePhoton") )   s.id = 7;
      else if( rootFileName.Contains("SingleMuon") )     s.id = 8;
      else if( rootFileName.Contains("SingleElectron") ) s.id = 9;
      else if( rootFileName.Contains("singletop") )	 s.id = 403;
      else if( rootFileName.Contains("mg_lo") )	         s.id = 330;
      else if( rootFileName.Contains("wjets") )	         s.id = 504;
      else if( rootFileName.Contains("dyjetsll") )	 s.id = 705;

    }

    //   std::cout << "s.id  = " << s.id << std::endl;
    //   std::cout << "evtid = " << evt_id << std::endl;
    
    s.nevents  = evt_nEvts;
    s.xsection = evt_xsec;
    s.filter   = evt_filter;
    s.kfact    = evt_kfactor;
    s.scale1fb = evt_scale1fb;

    if( (idMin>=0 && s.id<idMin) || (idMax>=0 && s.id>idMax) ) {
      file->Close();
      continue;
    }

    bool isData = s.id>0 && s.id<100;

    s.sname = fileNameParts[0];
    s.name  = fileNameParts[0];

    if( !isData ) {

      bool foundPU=false;
      int lastOnes = 3;
      for( unsigned i=1; i<fileNameParts.size() && lastOnes>0; ++i ) {
        if( fileNameParts[i] == "root" ) break;
        if( fileNameParts[i] == "babytree" ) continue;
        if( fileNameParts[i] == "prune" ) continue;
        if( fileNameParts[i] == "skim" ) continue;
	if( fileNameParts[i] == "post" ) continue;
        s.name += "_" + fileNameParts[i];
        if( fileNameParts[i]=="PU" ) foundPU = true;
        if( !foundPU )
          s.sname += "_" + fileNameParts[i];
        else
          lastOnes--;
      }

    } else {

      for( unsigned i=1; i<fileNameParts.size(); ++i ) {
        if( fileNameParts[i] == "root" ) break;
        if( fileNameParts[i] == "babytree" ) continue;
        if( fileNameParts[i] == "prune" ) continue;
        if( fileNameParts[i] == "skim" ) continue;
	if( fileNameParts[i] == "post" ) continue;
        s.name += "_" + fileNameParts[i];
        s.sname += "_" + fileNameParts[i];
      }

    }


    // Manually set variables
    s.PU_avg_weight = 1;
    s.lumi = 1000; // 1/fb
    
    std::cout << " ---- " << std::endl;
    std::cout << "  New sample added: " << s.name << std::endl;
    std::cout << "   Sample no.      " << counter << std::endl;
    std::cout << "   Short name:     " << s.sname << std::endl;
    std::cout << "   File:           " << s.file << std::endl;
    std::cout << "   ID:             " << s.id << std::endl;
    std::cout << "   Events:         " << s.nevents  << std::endl;
    std::cout << "   Xsection:       " << s.xsection << std::endl;
    //std::cout << "   Lumi:           " << s.lumi << std::endl;
    std::cout << "   Kfactor:        " << s.kfact << std::endl;
    std::cout << "   Filter:         " << s.filter << std::endl;
    std::cout << "   Scale to 1/fb:  " << s.scale1fb << std::endl;
    //std::cout << "   Avg PU weight:  " << s.PU_avg_weight << std::endl;
    //std::cout << "   type:           " << s.type << std::endl;
    //std::cout << "   Color:          " << s.color << std::endl;
    fSamples.push_back(s);
    file->Close();
    tree = 0;
    counter++;
    
  }

  std::cout << "------------------------------------" << std::endl;

  return fSamples;

}


