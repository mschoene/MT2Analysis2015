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


int round(float d) {
  return (int)(floor(d + 0.5));
}

int main( int argc, char* argv[] ) {

  if( argc==1 ) {

    std::cout << "Usage: ./printLatexBGTable.cpp [dir]" << std::endl;
    exit(191);

  }

  std::string dir(argv[1]);

  MT2Analysis<MT2Estimate>* analysis = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "ZJets" ); // any one is good, just need to know the regions

  std::set<MT2Region> regions = analysis->getRegions();

  std::string ofs_name(Form("latexBGTable_%s.tex", dir.c_str()));
  ofstream ofs;
  ofs.open( ofs_name, std::ofstream::app );

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR != regions.end(); ++iR ){
    std::cout << iR->getName() << std::endl;
    std::vector< std::string > names = iR->getNiceNamesLatex();

    ofs << "\\begin{table}[htbp]" << std::endl; 
    ofs << "\\caption{Background estimate yields for " << names[0].c_str() << ", " << names[1].c_str() << ".}" << std::endl;
    ofs << "\\centering" << std::endl;
    ofs << "\\begin{tabular}{r";
    int nBins;
    double *bins;
    iR->getBins(nBins, bins);
    for( int b=0; b<nBins; ++b ){
      ofs << "|c";
      //std::cout<< b <<"out of"<< nBins <<": " << bins[b] << std::endl;
    }
    ofs << "}" << std::endl;
    
    ofs << "\\hline" << std::endl;
    ofs << "\\multicolumn{" << nBins+1 << "}{c}{" << names[0].c_str() << ", " << names[1].c_str() << "}\\\\" << std::endl;
    ofs <<"\\hline \\hline" << std::endl;

    ofs << "Process";

    std::vector< std::string > binNames = iR->getBinNamesLatex();
    for( unsigned b=0; b < binNames.size(); ++b )
      ofs << " & " << binNames[b].c_str();
    ofs << "\\\\" << std::endl;
    ofs <<"\\hline \\hline" << std::endl;

    std::string zinvLine = "Invisible Z ";
    std::string llepLine = "Lost Lepton ";
    std::string qcdLine = "QCD ";

    for( int iBin=0; iBin<nBins; ++iBin ) {
      
      std::string tableName;
      if( iBin < nBins-1 )
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
    
    zinvLine = zinvLine + " \\\\";
    llepLine = llepLine + " \\\\";
    qcdLine  = qcdLine  + " \\\\";
    
    ofs << zinvLine << std::endl;
    ofs << llepLine << std::endl;
    ofs << qcdLine << std::endl;

    ofs << "\\hline" << std::endl;
    ofs << "\\end{tabular}" << std::endl;
    ofs << "\\end{table}" << std::endl;
  
  } // for selected regions
  
  
  ofs.close();
  
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

  std::string returnLine;
  if(yield<0.05 && statUp<0.05 && statDn<0.05 && systUp<0.05 && systDn<0.05)
    returnLine="---";
  else if(round(statUp) > 9 && round(statDn) > 9 && round(systUp) > 9 && round(systDn) > 9)
    returnLine=Form("%.0f $%s$(stat.) $%s$(syst.)", yield, statPart.c_str(), systPart.c_str());
  else
    returnLine=Form("%.1f $%s$(stat.) $%s$(syst.)", yield, statPart.c_str(), systPart.c_str());

  return returnLine;

}



std::string getSingleErrPart( float up, float dn ) {

  std::string thisPart;
  
  int UP=round(up*10.);
  int DN=round(dn*10.);
  
  if(round(up) > 9 && round(dn) > 9){
    
    UP=round(up);
    DN=round(dn);
    if( UP==DN )
      thisPart = std::string(Form("\\pm %.0f", up) );
    else
      thisPart = std::string(Form("^{+%.0f}_{-%.0f}", up, dn) );
  
  }
  else{
    
    if( UP==DN )
      thisPart = std::string(Form("\\pm %.1f", up) );
    else
      thisPart = std::string(Form("^{+%.1f}_{-%.1f}", up, dn) );

  }
 
  return thisPart;
 
}

