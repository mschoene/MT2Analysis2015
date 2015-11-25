#include <iostream>
#include <fstream>
#include <sstream>

#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "TMath.h"


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

  int data;

};



BGTable getTable( const std::string& tableFileName );
std::string makeSingleLine( float yield, float statUp, float statDn, float systUp, float systDn );
std::string makeSingleLineData( int yield );
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

  //  MT2Analysis<MT2Estimate>* analysisData = MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "data" ); // any one is good, just need to know the regions
  
  std::vector < MT2Analysis<MT2Estimate>* > analysesSignal;
  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1tttt_mGluino1500_mLSP100") );
  analysesSignal[0]->setName("T1tttt 1500,100");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1tttt_mGluino1200_mLSP800") );
  analysesSignal[1]->setName("T1tttt 1200,800");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1bbbb_mGluino1500_mLSP100") );
  analysesSignal[2]->setName("T1bbbb 1500,100");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1bbbb_mGluino1000_mLSP900") );
  analysesSignal[3]->setName("T1bbbb 1000,900");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1qqqq_mGluino1400_mLSP100") );
  analysesSignal[4]->setName("T1qqqq 1400,100");

  analysesSignal.push_back( MT2Analysis<MT2Estimate>::readFromFile( dir + "/analyses.root", "SMS_T1qqqq_mGluino1000_mLSP800") );
  analysesSignal[5]->setName("T1qqqq 1000,800");

  std::set<MT2Region> regions = analysis->getRegions();

  std::string ofs_name(Form("latexBGTable_%s.tex", dir.c_str()));
  std::cout << ofs_name << std::endl;
  std::ofstream ofs;
  ofs.open( ofs_name , std::ofstream::app );

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR != regions.end(); ++iR ){
    std::cout << iR->getName() << std::endl;
    std::vector< std::string > names = iR->getNiceNamesLatex();

    MT2Region* thisRegion = new MT2Region( *(iR) );

    ofs << "\\begin{table}[htbp]" << std::endl; 
    ofs << "\\caption{Background estimate, observation, and signal yields in bins of \\mttwo for " << names[0].c_str() << ", " << names[1].c_str() << ". The yields are normalized to $1.26~\\mathrm{fb}^{-1}$.}" << std::endl;
    ofs << "\\scriptsize" << std::endl;
    ofs << "\\centering" << std::endl;
    ofs << "\\makebox[\\textwidth][c]{" << std::endl;
    ofs << "\\begin{tabular}{r";
    int nBins;
    double *bins;
    iR->getBins(nBins, bins);
    for( int b=0; b<nBins; ++b ){
      ofs << "|c";
      //      std::cout<< b <<"out of"<< nBins <<": " << bins[b] << std::endl;
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
    std::string totLine = "Total ";
    std::string dataLine = "Observation ";

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
      std::string thistotLine =  makeSingleLine( thisTable.qcd+thisTable.llep+thisTable.zinv, TMath::Sqrt(thisTable.qcd_statUp*thisTable.qcd_statUp + thisTable.llep_statUp*thisTable.llep_statUp + thisTable.zinv_statUp*thisTable.zinv_statUp), TMath::Sqrt(thisTable.qcd_statDn*thisTable.qcd_statDn + thisTable.llep_statDn*thisTable.llep_statDn + thisTable.zinv_statDn*thisTable.zinv_statDn), TMath::Sqrt(thisTable.qcd_systUp*thisTable.qcd_systUp + thisTable.llep_systUp*thisTable.llep_systUp + thisTable.zinv_systUp*thisTable.zinv_systUp), TMath::Sqrt(thisTable.qcd_systDn*thisTable.qcd_systDn + thisTable.llep_systDn*thisTable.llep_systDn + thisTable.zinv_systDn*thisTable.zinv_systDn) );
      
      std::string thisDataLine =  makeSingleLineData( thisTable.data );
      
      zinvLine = zinvLine + " & " + thisZinvLine;
      llepLine = llepLine + " & " + thisllepLine;
      qcdLine = qcdLine + " & " + thisqcdLine;
      totLine = totLine + " & " + thistotLine;

      dataLine = dataLine + " & " + thisDataLine;
      
      
    } // for mt2 bins
    
    zinvLine = zinvLine + " \\\\";
    llepLine = llepLine + " \\\\";
    qcdLine  = qcdLine  + " \\\\";
    totLine  = totLine  + " \\\\";
    dataLine  = dataLine  + " \\\\";
    
    ofs << zinvLine << std::endl;
    ofs << llepLine << std::endl;
    ofs << qcdLine << std::endl;

    ofs << "\\hline" << std::endl;

    ofs << totLine << std::endl;

    ofs << "\\hline" << std::endl;

    ofs << dataLine << std::endl;
    
//    ofs << analysisData->getName().c_str();
//    analysisData->print( ofs, thisRegion );

    ofs << "\\hline" << std::endl;
    ofs << "\\hline" << std::endl;

    for(unsigned a =0;  a < analysesSignal.size(); ++a) {
      ofs << analysesSignal[a]->getName().c_str();
      analysesSignal[a]->print(ofs, thisRegion );
    }
    
    ofs << "\\hline" << std::endl;

    ofs << "\\end{tabular}}" << std::endl;
    ofs << "\\end{table}" << std::endl;
  
  } // for selected regions
  
  
  ofs.close();
  
  return 0;
  
}


BGTable getTable( const std::string& tableFileName ) {

  std::ifstream ifs( tableFileName.c_str() );

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
    } else if( name=="data" ){
      table.data = yield;
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
    returnLine=Form("%.2f $%s$(stat.) $%s$(syst.)", yield, statPart.c_str(), systPart.c_str());
  //    returnLine="---";
  else if(round(statUp) > 9 && round(statDn) > 9 && round(systUp) > 9 && round(systDn) > 9)
    returnLine=Form("%.0f $%s$(stat.) $%s$(syst.)", yield, statPart.c_str(), systPart.c_str());
  else
    returnLine=Form("%.1f $%s$(stat.) $%s$(syst.)", yield, statPart.c_str(), systPart.c_str());

  return returnLine;

}


std::string makeSingleLineData( int yield ) {

  std::string returnLine;
  returnLine=Form("%d", yield);

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
  else if( up<0.05 || dn < 0.05){
    
    if( UP==DN )
      thisPart = std::string(Form("\\pm %.2f", up) );
    else
      thisPart = std::string(Form("^{+%.2f}_{-%.2f}", up, dn) );

  }
  else{
    
    if( UP==DN )
      thisPart = std::string(Form("\\pm %.1f", up) );
    else
      thisPart = std::string(Form("^{+%.1f}_{-%.1f}", up, dn) );

  }
 
  return thisPart;
 
}

