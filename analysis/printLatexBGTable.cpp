#include <iostream>
#include <fstream>
#include <sstream>

#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"



struct BGTable {

  float zinv;
  float zinv_statUp;
  float zinv_statDn;
  float zinv_systUp;
  float zinv_systDn;

  float llep;
  float llep_statUp;
  float llep_statDn;
  float llep_systUp;
  float llep_systDn;

  float qcd;
  float qcd_statUp;
  float qcd_statDn;
  float qcd_systUp;
  float qcd_systDn;

};



BGTable getTable( const std::string& tableFileName );
std::string makeSingleLine( float yield, float statUp, float statDn, float systUp, float systDn );
std::string getSingleErrPart( float up, float dn );


int main( int argc, char* argv[] ) {

  if( argc==1 ) {

    std::cout << "Usage: ./printLatexBGTable.cpp [dir]" << std::endl;
    exit(191);

  }

  std::string dir(argv[1]);

  MT2Analysis<MT2Estimate>* analysis = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "ZJets" ); // any one is good, just need to know the regions

  std::set<MT2Region> regions = analysis->getRegions();
  std::set<MT2HTRegion> htregions = analysis->getHTRegions();

  ofstream ofs(Form("latexBGTable_%s.tex", dir.c_str()));

  for( std::set<MT2HTRegion>::iterator iHT = htregions.begin(); iHT != htregions.end(); ++iHT ) {

    // find corresponding regions:
    std::set<MT2Region> selectedRegions;
    for( std::set<MT2Region>::iterator iR = regions.begin(); iR != regions.end(); ++iR ) {

      if( iR->htRegion()->isIncluded( &(*iHT)) ) continue;
      selectedRegions.insert( *iR );

    }

    for( std::set<MT2Region>::iterator iR = selectedRegions.begin(); iR != selectedRegions.end(); ++iR )
      std::cout << iR->getName() << std::endl;


    ofs << "\\begin{table}[htbp]" << std::endl; 
    ofs << "\\caption{Background estimate yields for $" << iHT->getNiceName() << "$}" << std::endl;
    ofs << "\\caption{Background estimate yields for $" << iHT->getNiceName() << "$}" << std::endl;
    ofs << "\\centering" << std::endl;
    ofs << "\\begin{tabular}{r";
    for( std::set<MT2Region>::iterator iR = selectedRegions.begin(); iR != selectedRegions.end(); ++iR )
      ofs << "|c";
    ofs << "}" << std::endl;
    
    std::string zinvLine = "Invisible Z ";
    std::string llepLine = "Lost Lepton ";
    std::string qcdLine = "QCD ";

    for( std::set<MT2Region>::iterator iR = selectedRegions.begin(); iR != selectedRegions.end(); ++iR ) {

      int nBins;
      double* bins;
      iR->getBins(nBins, bins );

      for( int iBin=0; iBin<nBins-1; ++iBin ) {

        std::string tableName;
        if( bins[iBin+1]>0. )
          tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lfto%.0lf.txt", dir.c_str(), iR->getName().c_str(), bins[iBin], bins[iBin+1]) );
        else 
          tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lftoInf.txt", dir.c_str(), iR->getName().c_str(), bins[iBin] ));

        BGTable thisTable = getTable(tableName);

        std::string thisZinvLine = makeSingleLine( thisTable.zinv, thisTable.zinv_statUp, thisTable.zinv_statDn, thisTable.zinv_systUp, thisTable.zinv_systDn );
        std::string thisllepLine = makeSingleLine( thisTable.llep, thisTable.llep_statUp, thisTable.llep_statDn, thisTable.llep_systUp, thisTable.llep_systDn );
        std::string thisqcdLine = makeSingleLine( thisTable.qcd, thisTable.qcd_statUp, thisTable.qcd_statDn, thisTable.qcd_systUp, thisTable.qcd_systDn );
        
        zinvLine = zinvLine + " & " + thisZinvLine;
        llepLine = llepLine + " & " + thisllepLine;
        qcdLine = qcdLine + " & " + thisqcdLine;

      } // for mt2 bins

      ofs << zinvLine << std::endl;
      ofs << llepLine << std::endl;
      ofs << qcdLine << std::endl;

    } // for selected regions

  } // for ht regions

  ofs << "\\hline" << std::endl;
  ofs << "\\end{tabular}" << std::endl;
  ofs << "\\end{table}" << std::endl;

  return 0;

}




BGTable getTable( const std::string& tableFileName ) {

  ifstream ifs( tableFileName.c_str() );

  BGTable table;

  while( ifs.good() ) {


    char thisLine[256];
    ifs.getline( thisLine, 256 );
    if( thisLine[0]=='#' ) continue;

    std::istringstream thisLine_iss(thisLine);
    
    std::string name;
    float yield, statUp, statDn, systUp, systDn;
    thisLine_iss >> name >> yield >> statUp >> statDn >> systUp >> systDn;

    if( name=="zinv" ) {
      table.zinv = yield;
      table.zinv_statUp = statUp;
      table.zinv_statDn = statDn;
      table.zinv_systUp = systUp;
      table.zinv_systDn = systDn;
    } else if( name=="llep" ) {
      table.llep = yield;
      table.llep_statUp = statUp;
      table.llep_statDn = statDn;
      table.llep_systUp = systUp;
      table.llep_systDn = systDn;
    } else if( name=="qcd" ) {
      table.qcd = yield;
      table.qcd_statUp = statUp;
      table.qcd_statDn = statDn;
      table.qcd_systUp = systUp;
      table.qcd_systDn = systDn;
    } else {
      continue;
    }

  }

  return table;

}


std::string makeSingleLine( float yield, float statUp, float statDn, float systUp, float systDn ) {


  std::string statPart = getSingleErrPart(statUp, statDn);
  std::string systPart = getSingleErrPart(systUp, systDn);

  std::string returnLine(Form("%.3f $%s$(stat.) $%s$(syst.)", yield, statPart.c_str(), systPart.c_str()) );

  return returnLine;

}



std::string getSingleErrPart( float up, float dn ) {

  std::string thisPart;

  if( up==dn ) 
    thisPart = std::string(Form("#pm %.3f", up) );
  else
    thisPart = std::string(Form("^{+%.3f}_{-%.3f}", up, dn) );

  return thisPart;

}

