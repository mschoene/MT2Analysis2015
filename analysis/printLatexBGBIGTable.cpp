#include <iostream>
#include <fstream>
#include <sstream>

#include "interface/MT2Analysis.h"
#include "interface/MT2Estimate.h"
#include "interface/MT2EstimateSyst.h"
#include "interface/MT2EstimateSigContSyst.h"
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

  int zinv_nCR;
  int llep_nCR;
  int qcd_nCR;

};



BGTable getTable( const std::string& tableFileName );
std::string makeSingleLine( float yield, float statUp, float statDn, float systUp, float systDn );
std::string makeSingleLineData( int yield );
std::string makeSingleLinePost( float yield, float stat );
std::string getSingleErrPart( float yield, float up, float dn, float otherUp, float otherDn );

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

  std::string ofs_name(Form("latexBGBIGTable_%s.tex", dir.c_str()));
  std::cout << ofs_name << std::endl;
  std::ofstream ofs;
  ofs.open( ofs_name , std::ofstream::app );

  int regionIndex=0;
  for( std::set<MT2Region>::iterator iR = regions.begin(); iR != regions.end(); ++iR ){
    std::cout << iR->getName() << std::endl;
    std::vector< std::string > names = iR->getNiceNamesLatex();

    MT2Region* thisRegion = new MT2Region( *(iR) );

    std::string oldHTname="";
    std::string oldJname ="";

    if( ( iR->nJetsMin() > 1 && (regionIndex==0 || regionIndex>10)) || ( iR->nJetsMax()==1 && iR->htMin()<200+1 && iR->nBJetsMax()==0 ) ){
      
      oldHTname=names[0];
      oldJname=names[1];
      regionIndex=0;
      
      ofs << std::endl;
      ofs << "\\clearpage" << std::endl;
      ofs << "\\begin{table}[htbp]" << std::endl; 
      if( iR->nJetsMin() > 1 )
	ofs << "\\caption{Background estimate and observation in bins of \\mttwo for " << names[0].c_str() <<". The yields are normalized to $2.3~\\mathrm{fb}^{-1}$.}" << std::endl;
      else
	ofs << "\\caption{Background estimate and observation in bins of jet \\Pt for the monojet regions. The yields are normalized to $2.3~\\mathrm{fb}^{-1}$.}" << std::endl;
      
      ofs << "\\scriptsize" << std::endl;
      ofs << "\\centering" << std::endl;
      ofs << "\\makebox[\\textwidth][c]{" << std::endl;
      ofs << "\\begin{tabular}{c|c|c|c|c|c|c}" << std::endl;
      
      ofs << "\\hline" << std::endl;
      if( iR->nJetsMin() > 1 )
	ofs << "\\multicolumn{7}{c}{" << names[0].c_str() << "}\\\\" << std::endl;
      else
	ofs << "\\multicolumn{7}{c}{\\njets $= 1$}\\\\" << std::endl;

      ofs <<"\\hline" << std::endl;
      
      if(iR->nJetsMin() > 1)
	ofs << "\\njets, \\nbtags & \\mttwo [GeV] & Z $\\rightarrow \\nu\\bar{\\nu}$ & Lost lepton & Multijet & Total background & Data \\\\" << std::endl;
      else
	ofs << "\\njets, \\nbtags & Jet \\Pt [GeV] & Z $\\rightarrow \\nu\\bar{\\nu}$ & Lost lepton & Multijet & Total background & Data \\\\"<< std::endl;
      
      ofs <<"\\hline" << std::endl;
      
    }
    
    int nBins;
    double *bins;
    iR->getBins(nBins, bins);
    
    regionIndex+=1;

    if( iR->nJetsMin() > 1 ){

      ofs << "\\multirow{" << nBins << "}{*}{" << names[1].c_str() << "}" << std::endl;
      
      for( int iBin=0; iBin<nBins; ++iBin ) {
	
	std::string tableName;
	if( iR->nJetsMax()==1 )
	  tableName = std::string(Form("%s/datacard_templates/table_%s_m0toInf.txt", dir.c_str(), iR->getName().c_str()) );
	else{
	  if( iBin < nBins-1 )
	    tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lfto%.0lf.txt", dir.c_str(), iR->getName().c_str(), bins[iBin], bins[iBin+1]) );
	  else 
	    tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lftoInf.txt", dir.c_str(), iR->getName().c_str(), bins[iBin] ));
	}
	
	std::cout << "Reading table " << tableName << std::endl;
	
	BGTable thisTable = getTable(tableName);
	
	std::string thisZinvLine = makeSingleLine( thisTable.zinv, thisTable.zinv_statUp, thisTable.zinv_statDn, thisTable.zinv_systUp, thisTable.zinv_systDn );
	std::string thisllepLine = makeSingleLine( thisTable.llep, thisTable.llep_statUp, thisTable.llep_statDn, thisTable.llep_systUp, thisTable.llep_systDn );
	std::string thisqcdLine = makeSingleLine( thisTable.qcd, thisTable.qcd_statUp, thisTable.qcd_statDn, thisTable.qcd_systUp, thisTable.qcd_systDn );
	std::string thistotLine =  makeSingleLine( thisTable.qcd+thisTable.llep+thisTable.zinv, TMath::Sqrt(thisTable.qcd_statUp*thisTable.qcd_statUp + thisTable.llep_statUp*thisTable.llep_statUp + thisTable.zinv_statUp*thisTable.zinv_statUp), TMath::Sqrt(thisTable.qcd_statDn*thisTable.qcd_statDn + thisTable.llep_statDn*thisTable.llep_statDn + thisTable.zinv_statDn*thisTable.zinv_statDn), TMath::Sqrt(thisTable.qcd_systUp*thisTable.qcd_systUp + thisTable.llep_systUp*thisTable.llep_systUp + thisTable.zinv_systUp*thisTable.zinv_systUp), TMath::Sqrt(thisTable.qcd_systDn*thisTable.qcd_systDn + thisTable.llep_systDn*thisTable.llep_systDn + thisTable.zinv_systDn*thisTable.zinv_systDn) );
	
	std::string thisDataLine =  makeSingleLineData( thisTable.data );
	
	std::string thisLine;
	if( iBin < nBins-1 )
	  thisLine = Form(" & $%.0lf-%.0lf$ &", bins[iBin], bins[iBin+1]);
	else
	  thisLine = Form(" & $\\gt %.0lf$ &", bins[iBin]);

	thisLine = thisLine  + thisZinvLine + " & " + thisllepLine + " & " + thisqcdLine + " & " + thistotLine + " & " + thisDataLine + " \\\\";
	ofs << thisLine << std::endl;
	
      } // for mt2 bins
      

      ofs << "\\hline" << std::endl;
      
      if(regionIndex > 10){

	ofs << "\\end{tabular}}" << std::endl;
	ofs << "\\end{table}" << std::endl;
	ofs << std::endl;

      }
      
    }
 
    else if( iR->nJetsMax() == 1 ){

      if( regionIndex==1 )
	ofs << "\\multirow{" << 7 << "}{*}{" << names[1].c_str() << "}" << std::endl;
      else if( regionIndex==8 ) 
	ofs << "\\multirow{" << 5 << "}{*}{" << names[1].c_str() << "}" << std::endl;

      for( int iBin=0; iBin<nBins; ++iBin ) {
	
	std::string tableName;
	if( iR->nJetsMax()==1 )
	  tableName = std::string(Form("%s/datacard_templates/table_%s_m0toInf.txt", dir.c_str(), iR->getName().c_str()) );
	else{
	  if( iBin < nBins-1 )
	    tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lfto%.0lf.txt", dir.c_str(), iR->getName().c_str(), bins[iBin], bins[iBin+1]) );
	  else 
	    tableName = std::string(Form("%s/datacard_templates/table_%s_m%.0lftoInf.txt", dir.c_str(), iR->getName().c_str(), bins[iBin] ));
	}
	
	std::cout << "Reading table " << tableName << std::endl;
	
	BGTable thisTable = getTable(tableName);
	
	std::string thisZinvLine = makeSingleLine( thisTable.zinv, thisTable.zinv_statUp, thisTable.zinv_statDn, thisTable.zinv_systUp, thisTable.zinv_systDn );
	std::string thisllepLine = makeSingleLine( thisTable.llep, thisTable.llep_statUp, thisTable.llep_statDn, thisTable.llep_systUp, thisTable.llep_systDn );
	std::string thisqcdLine = makeSingleLine( thisTable.qcd, thisTable.qcd_statUp, thisTable.qcd_statDn, thisTable.qcd_systUp, thisTable.qcd_systDn );
	std::string thistotLine =  makeSingleLine( thisTable.qcd+thisTable.llep+thisTable.zinv, TMath::Sqrt(thisTable.qcd_statUp*thisTable.qcd_statUp + thisTable.llep_statUp*thisTable.llep_statUp + thisTable.zinv_statUp*thisTable.zinv_statUp), TMath::Sqrt(thisTable.qcd_statDn*thisTable.qcd_statDn + thisTable.llep_statDn*thisTable.llep_statDn + thisTable.zinv_statDn*thisTable.zinv_statDn), TMath::Sqrt(thisTable.qcd_systUp*thisTable.qcd_systUp + thisTable.llep_systUp*thisTable.llep_systUp + thisTable.zinv_systUp*thisTable.zinv_systUp), TMath::Sqrt(thisTable.qcd_systDn*thisTable.qcd_systDn + thisTable.llep_systDn*thisTable.llep_systDn + thisTable.zinv_systDn*thisTable.zinv_systDn) );
	
	std::string thisDataLine =  makeSingleLineData( thisTable.data );
	
	std::string thisLine;
	if( iR->htMax() > 0 )
	  thisLine = Form(" & $%.0lf-%.0lf$ &", iR->htMin(), iR->htMax());
	else
	  thisLine = Form(" & $\\gt %.0lf$ &", iR->htMin());

	thisLine = thisLine  + thisZinvLine + " & " + thisllepLine + " & " + thisqcdLine + " & " + thistotLine + " & " + thisDataLine + " \\\\";
	ofs << thisLine << std::endl;
	
      } // for mt2 bins
      

      if(regionIndex==7 || regionIndex==12)
	ofs << "\\hline" << std::endl;
      
      if(regionIndex==12){

	ofs << "\\end{tabular}}" << std::endl;
	ofs << "\\end{table}" << std::endl;
	ofs << std::endl;

      }
      
    }
   
  

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
//    float yield, statUp, statDn, systUp, systDn;
    thisLine_iss >> name;// >> yield >> statUp >> statDn >> systUp >> systDn;

    if( name=="zinv" ) {
      float yield, statUp, statDn, systUp, systDn;
      thisLine_iss >> yield >> statUp >> statDn >> systUp >> systDn;
      table.zinv = yield;
      table.zinv_statUp = statUp;
      table.zinv_statDn = statDn;
      table.zinv_systUp = systUp;
      table.zinv_systDn = systDn;
    } else if( name=="llep" ) {
      float yield, statUp, statDn, systUp, systDn;
      thisLine_iss >> yield >> statUp >> statDn >> systUp >> systDn;
      table.llep = yield;
      table.llep_statUp = statUp;
      table.llep_statDn = statDn;
      table.llep_systUp = systUp;
      table.llep_systDn = systDn;
    } else if( name=="qcd" ) {
      float yield, statUp, statDn, systUp, systDn;
      thisLine_iss >> yield >> statUp >> statDn >> systUp >> systDn;
      table.qcd = yield;
      table.qcd_statUp = statUp;
      table.qcd_statDn = statDn;
      table.qcd_systUp = systUp;
      table.qcd_systDn = systDn;
    } else if( name=="data" ){
      int yield;
      thisLine_iss >> yield;
      table.data = yield;
    } else if( name=="zinv_nCR" ){
      int yield;
      thisLine_iss >> yield;
      table.zinv_nCR = yield;
    } else if( name=="llep_nCR" ){
      int yield;
      thisLine_iss  >> yield;
      table.llep_nCR = yield;
    } else if( name=="qcd_nCR" ){
      int yield;
      thisLine_iss  >> yield;
      table.qcd_nCR = yield;
    } else {
      continue;
    }

  }

  return table;

}


std::string makeSingleLine( float yield, float statUp, float statDn, float systUp, float systDn ) {


  std::string statPart = getSingleErrPart(yield, statUp, statDn, systUp, systDn);
  std::string systPart = getSingleErrPart(yield, systUp, systDn, statUp, statDn);

  std::string returnLine;
  if(yield<0.05 && statUp<0.05 && statDn<0.05 && systUp<0.05 && systDn<0.05)
    returnLine=Form("%.2f $%s$(stat.) $%s$(syst.)", yield, statPart.c_str(), systPart.c_str());
  //    returnLine="---";
  else if(round(yield)>9 || round(statUp) > 9 || round(statDn) > 9 || round(systUp) > 9 || round(systDn) > 9)
    returnLine=Form("%.0f $%s$(stat.) $%s$(syst.)", yield, statPart.c_str(), systPart.c_str());
  else
    returnLine=Form("%.1f $%s$(stat.) $%s$(syst.)", yield, statPart.c_str(), systPart.c_str());

  return returnLine;

}

std::string makeSingleLinePost( float yield, float stat ) {


  std::string statPart = getSingleErrPart(yield, stat, stat, stat, stat);

  std::string returnLine;
  if(yield<0.05 && stat<0.05)
    returnLine=Form("%.2f $%s$", yield, statPart.c_str());
  else if( round(stat) > 9 )
    returnLine=Form("%.0f $%s$", yield, statPart.c_str());
  else
    returnLine=Form("%.1f $%s$", yield, statPart.c_str());

  return returnLine;

}


std::string makeSingleLineData( int yield ) {

  std::string returnLine;
  returnLine=Form("%d", yield);

  return returnLine;

}



std::string getSingleErrPart( float yield, float up, float dn, float otherUp, float otherDn ) {

  std::string thisPart;
  
  int Y=round(yield*10.);
  int UP=round(up*10.);
  int DN=round(dn*10.);
  
  int otherUP=round(otherUp*10.);
  int otherDN=round(otherDn*10.);
  
  if(round(yield)>9 || round(up) > 9 || round(dn) > 9 || round(otherUp) > 9 || round(otherDn) > 9){
    
    UP=round(up);
    DN=round(dn);
    Y=round(yield);

    if( UP==DN && Y-DN>=0. )
      thisPart = std::string(Form("\\pm %.0f", up) );
    else{
      if(Y-DN>=0.)
	thisPart = std::string(Form("^{+%.0f}_{-%.0f}", up, dn) );
      else
	thisPart = std::string(Form("^{+%.0f}_{-%.0f}", up, yield) );
    }
  }
  else if( yield<0.05 && up<0.05 && dn < 0.05 && otherUp < 0.05 && otherDn < 0.05 ){
    
    Y=round(yield*100.);
    UP=round(up*100.);
    DN=round(dn*100.);

    if( UP==DN && Y-DN>=0. )
      thisPart = std::string(Form("\\pm %.2f", up) );
    else{
      if(Y-DN>=0.)
	thisPart = std::string(Form("^{+%.2f}_{-%.2f}", up, dn) );
      else
	thisPart = std::string(Form("^{+%.2f}_{-%.2f}", up, yield) );
    }
  }
  else{
    
    if( UP==DN && Y-DN>=0. )
      thisPart = std::string(Form("\\pm %.1f", up) );
    else{
      if(Y-DN>=0.)
	thisPart = std::string(Form("^{+%.1f}_{-%.1f}", up, dn) );
      else
	thisPart = std::string(Form("^{+%.1f}_{-%.1f}", up, yield) );
    }
  }
 
  return thisPart;
 
}

