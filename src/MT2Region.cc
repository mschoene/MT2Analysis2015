#include "../interface/MT2Region.h"

#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <cstdlib>

#include "TH1F.h"


////////////////////////////////////////////////////////
//
//                    MT2HTRegion
//
////////////////////////////////////////////////////////



MT2HTRegion::MT2HTRegion( const std::string& name ) {

  // this constructor parses the name
  // the name has to be passed in the format:
  // HT[htMin]to[htMax]
  // where :
  //   - htMin and htMax have to be integers
  //   - htMax is allowed to be "Inf"


  std::stringstream ss(name);
  std::vector<std::string> parts;
  std::string item;
  while(std::getline(ss, item, '_')) {
    parts.push_back(item);
  }


  if( parts.size()!=1 ) {
    std::cout << "[MT2HTRegion]::MT2HTRegion ERROR! Unrecognized MT2HTRegion name: " << name << std::endl;
    exit(455);
  }


  int htMin(-1), htMax(-1);

  char htMax_char[100];
  sscanf( parts[0].c_str(), "HT%dto%s", &htMin, htMax_char);
  std::string htMax_str(htMax_char);
  if( htMax_str=="Inf" ) 
    htMax = -1;
  else
    sscanf( parts[0].c_str(), "HT%dto%d", &htMin, &htMax);


  this->htMin  = htMin;
  this->htMax  = htMax;

}





MT2HTRegion::MT2HTRegion( const MT2HTRegion& rhs ) {

  htMin = rhs.htMin;
  htMax = rhs.htMax;
  
}



MT2HTRegion::MT2HTRegion( float ahtMin, float ahtMax ) {

  htMin = ahtMin;
  htMax = ahtMax;

}



//float MT2HTRegion::metMin( float ht ) const {
//
//  float metMin = (ht<1000.) ? 200. : 30.;
//  return metMin;
//
//}






//float MT2HTRegion::metMin() const {
// 
//  float metMin = (htMin<999.) ? 200. : 30.;
//  return metMin;
// 
//}
//bool MT2HTRegion::isInclusiveHT() const {
//
//  if( htMin<999. && ( htMax < 0 || htMax >= 1000. )) return true;
//  else return false;
//  
//}
//
//
//float MT2HTRegion::metMinInclusiveHT( float ht ) const {
//  
//  float metMin = (ht<1000.) ? 200. : 30.;
//  return metMin;
//
//}

std::string MT2HTRegion::getName() const {

  std::string htMax_str(Form("%.0f", htMax));
  if( htMax==-1 ) htMax_str = "Inf";

  char n[512];
  sprintf( n, "HT%.0fto%s", htMin, htMax_str.c_str() );
  std::string n_str(n);

  return n_str;

}


std::string MT2HTRegion::getNiceName() const {

  std::string htMax_str(Form("%.0f", htMax));
  char htPart[500];
  if( htMax==-1 ) 
    sprintf( htPart, "H_{T} > %.0f GeV", htMin );
  else
    sprintf( htPart, "%.0f < H_{T} < %.0f GeV", htMin, htMax );
  std::string htPart_str(htPart);

  return htPart_str;

}

std::string MT2HTRegion::getNiceNameLatex() const {

  std::string htMax_str(Form("%.0f", htMax));
  char htPart[500];
  if( htMax==-1 ) 
    sprintf( htPart, "H$_{\\mathrm{T}} > %.0f$~GeV", htMin );
  else
    sprintf( htPart, "$%.0f < $H$_{\\mathrm{T}} < %.0f$~GeV", htMin, htMax );
  std::string htPart_str(htPart);

  return htPart_str;

}


std::string MT2HTRegion::getCuts() const {

  std::string cuts;
  if( htMax<0. ) 
    cuts = std::string( Form( "ht >= %f", htMin) );
  else
    cuts = std::string( Form( "ht >= %f && ht<%f", htMin, htMax) );

  return cuts;

}


bool MT2HTRegion::operator==( const MT2HTRegion& rhs ) const {

  return (this->getName()==rhs.getName());

}


bool MT2HTRegion::operator!=( const MT2HTRegion& rhs ) const {

  return (this->getName()!=rhs.getName());

}


bool MT2HTRegion::operator<( const MT2HTRegion& rhs ) const {

  if( *this==rhs ) return true;

  float thisHtMax = (htMax>=0.) ? htMax : 99999.;
  float rhsHtMax = (rhs.htMax>=0.) ? rhs.htMax : 99999.;

  bool returnBool;

  if( htMin==rhs.htMin ) {
    
    returnBool = thisHtMax<rhsHtMax;

  } else {

    returnBool = ( htMin<rhs.htMin );
  }

  return returnBool;

}



bool MT2HTRegion::isIncluded( const MT2HTRegion* htRegion ) const {

  bool returnBool = true;

  if( htMin < htRegion->htMin ) returnBool = false;
  if( htMax > htRegion->htMax && htRegion->htMax >=0. ) returnBool = false;
  if( htMax < htRegion->htMax && htMax < 0 ) returnBool = false; //To solve case where htMax is Inf and htRegion->htMax is finite


  return returnBool;

}


////////////////////////////////////////////////////////
//
//                  MT2SignalRegion
//
////////////////////////////////////////////////////////


MT2SignalRegion::MT2SignalRegion( const std::string& name ) {

  std::stringstream ss(name);
  std::vector<std::string> parts;
  std::string item;
  while(std::getline(ss, item, '_')) {
    parts.push_back(item);
  }


  if( parts.size()!=2 && parts.size()!=3 && parts.size()!=4 ) {
    std::cout << "[MT2SignalRegion]::MT2SignalRegion ERROR! Unrecognized MT2SignalRegion name: " << name << std::endl;
    exit(457);
  }


  int jMin(-1), jMax(-1);

  if( parts[0].size()<= 3 ) {
    sscanf( parts[0].c_str(), "j%d", &jMin );
    jMax = jMin;
  } else {
    char jMax_char[100];
    sscanf( parts[0].c_str(), "j%dto%s", &jMin, jMax_char);
    std::string jMax_str(jMax_char);
    if( jMax_str=="Inf" ) 
      jMax = -1;
    else
      sscanf( parts[0].c_str(), "j%dto%d", &jMin, &jMax);
  }


  int bMin(-1), bMax(-1);

  if( parts[1].size()<= 3 ) {
    sscanf( parts[1].c_str(), "b%d", &bMin );
    bMax = bMin;
  } else {
    char bMax_char[100];
    sscanf( parts[1].c_str(), "b%dto%s", &bMin, bMax_char);
    std::string bMax_str(bMax_char);
    if( bMax_str=="Inf" ) 
      bMax = -1;
    else
      sscanf( parts[1].c_str(), "b%dto%d", &bMin, &bMax);
  }


  nJetsMin  = jMin;
  nJetsMax  = jMax;
  nBJetsMin = bMin;
  nBJetsMax = bMax;
  
  
  mtCut = "";
  if( parts.size()>2 ) mtCut = parts[2];

  if( mtCut!="" && mtCut!="loMT" && mtCut!="hiMT" ) {
    std::cout << "[MT2SignalRegion::MT2SignalRegion] ERROR! Unkown mtCut '" << mtCut << "'!" << std::endl;
    exit(3535);
  }


}




MT2SignalRegion::MT2SignalRegion(int njmin, int njmax, int nbmin, int nbmax, const std::string& mtcut ) {

  nJetsMin = njmin;
  nJetsMax = njmax;
  nBJetsMin = nbmin;
  nBJetsMax = nbmax;
  //if( nJetsMin<nBJetsMin ) nJetsMin = nBJetsMin;

  mtCut = mtcut;

}




MT2SignalRegion::MT2SignalRegion( const MT2SignalRegion& rhs ) {

  nJetsMin = rhs.nJetsMin;
  nJetsMax = rhs.nJetsMax;
  nBJetsMin = rhs.nBJetsMin;
  nBJetsMax = rhs.nBJetsMax;

  mtCut = rhs.mtCut;

}



std::string MT2SignalRegion::getName() const {
 
  std::string jString = getSingleJetString( "j", nJetsMin,  nJetsMax  );
  std::string bString = getSingleJetString( "b", nBJetsMin, nBJetsMax );

  std::string signal_region = jString + "_" + bString;
  if( mtCut!="" ) signal_region = signal_region + "_" + mtCut;

  return signal_region;

}



std::string MT2SignalRegion::getSingleJetString( const std::string& prefix, int n_min , int n_max ) const {

  std::string n_min_str(Form("%d", n_min));
  std::string n_max_str(Form("%d", n_max));

  if( n_max<0 ) n_max_str="Inf";
  if( n_min<0 ) n_min_str=(prefix=="j") ? "2" : "0";
  
  std::string signal_region;
  if( n_min_str!=n_max_str )
    signal_region = prefix + n_min_str + "to" +  n_max_str;
  else
    signal_region = prefix + n_min_str;

  return signal_region;

}




std::string MT2SignalRegion::getNiceName() const {

  std::string niceName_j = getNiceJetName( "j", nJetsMin,  nJetsMax  );
  std::string niceName_b = getNiceJetName( "b", nBJetsMin,  nBJetsMax  );

  std::string niceName = niceName_j;
  if( niceName!="" && niceName_b!="" ) niceName += ", " + niceName_b;

  if( mtCut=="loMT"  ) niceName += " (low M_{T})";
  else if( mtCut=="hiMT" ) niceName += " (high M_{T})";

  return niceName;

}


std::string MT2SignalRegion::getNiceJetName( const std::string& pedix, int nmin, int nmax ) const {

  if( nmin==-1 && nmax==-1 ) return std::string("");

  char n[500];
//  if( nmax==nmin )
//    sprintf( n, "N(%s) = %d", pedix.c_str(), nmin );
//  else {
//    if( nmax==-1 )
//      sprintf( n, "N(%s) #geq %d", pedix.c_str(), nmin );
//    else
//      sprintf( n, "%d #leq N(%s) #leq %d", nmin, pedix.c_str(), nmax );
//  }
  if( nmax==nmin )
    sprintf( n, "%d%s", nmin, pedix.c_str() );
  else {
    if( nmax==-1 )
      sprintf( n, "#geq%d%s", nmin, pedix.c_str() );
    else
      sprintf( n, "%d-%d%s", nmin, nmax, pedix.c_str() );
  }

  std::string nicename(n);

  return nicename;

}

std::string MT2SignalRegion::getNiceNameLatex() const {

  std::string niceName_j = getNiceJetNameLatex( "j", nJetsMin,  nJetsMax  );
  std::string niceName_b = getNiceJetNameLatex( "b", nBJetsMin,  nBJetsMax  );

  std::string niceName = niceName_j;
  if( niceName!="" && niceName_b!="" ) niceName += ", " + niceName_b;

  if( mtCut=="loMT"  ) niceName += " (low M_{T})";
  else if( mtCut=="hiMT" ) niceName += " (high M_{T})";

  return niceName;

}


std::string MT2SignalRegion::getNiceJetNameLatex( const std::string& pedix, int nmin, int nmax ) const {

  if( nmin==-1 && nmax==-1 ) return std::string("");

  char n[500];
//  if( nmax==nmin )
//    sprintf( n, "N(%s) = %d", pedix.c_str(), nmin );
//  else {
//    if( nmax==-1 )
//      sprintf( n, "N(%s) #geq %d", pedix.c_str(), nmin );
//    else
//      sprintf( n, "%d #leq N(%s) #leq %d", nmin, pedix.c_str(), nmax );
//  }
  if( nmax==nmin )
    sprintf( n, "$%d$%s", nmin, pedix.c_str() );
  else {
    if( nmax==-1 )
      sprintf( n, "$\\geq%d$%s", nmin, pedix.c_str() );
    else
      sprintf( n, "$%d-%d$%s", nmin, nmax, pedix.c_str() );
  }

  std::string nicename(n);

  return nicename;

}



std::string MT2SignalRegion::getCuts() const {

  std::string  jetCuts = this->getJetCuts();
  std::string bjetCuts = this->getBJetCuts();
  std::string fullCuts = jetCuts + " && " + bjetCuts;

  return fullCuts;

}



std::string MT2SignalRegion::getJetCuts() const {

  std::string cuts( Form("nJets>=%d", nJetsMin) );
  if( nJetsMax>=0 ) {
    std::string addition( Form(" && nJets<=%d", nJetsMax) );
    cuts += addition;
  }

  return cuts;

}


std::string MT2SignalRegion::getBJetCuts() const {

  std::string cuts( Form("nBJets>=%d", nBJetsMin) );
  if( nBJetsMax>=0 ) {
    std::string addition( Form(" && nBJets<=%d", nBJetsMax) );
    cuts += addition;
  }

  return cuts;

}





bool MT2SignalRegion::operator==( const MT2SignalRegion& rhs ) const {

  return (this->getName()==rhs.getName());
 
}


bool MT2SignalRegion::operator!=( const MT2SignalRegion& rhs ) const {

  return (this->getName()!=rhs.getName());
 
}
    

bool MT2SignalRegion::operator<( const MT2SignalRegion& rhs ) const {

  if( *this == rhs ) return false;

  int  thisNJmax = (nJetsMax>=0) ? nJetsMax : 99999;
  int  rhsNJmax = (rhs.nJetsMax>=0) ? rhs.nJetsMax : 99999;

  bool returnBool;
  
  if( nJetsMax == 1 && (rhs.nJetsMax > 1 || rhs.nJetsMax < 0) ){
    
    returnBool = true;

  }
  else if ( rhs.nJetsMax == 1 && (nJetsMax > 1 || nJetsMax < 0) ){

    returnBool = false;
    
  } 
  else if ( nJetsMax == 1 && rhs.nJetsMax == 1 ){
    
    returnBool = ( nBJetsMin < rhs.nBJetsMin );

  }
  else{
    if( (nBJetsMax >= 0 && rhs.nBJetsMax >= 0) || (rhs.nBJetsMax < 0 && nBJetsMax < 0) ) {

      if( thisNJmax == rhsNJmax ) {
	
	if( nBJetsMin!=rhs.nBJetsMin ) {
	  
	  returnBool = ( nBJetsMin<rhs.nBJetsMin );
	  
	} else {
	  
	  if( mtCut!=rhs.mtCut ) {
	    
	    if( mtCut=="loMT" ) {
	      returnBool = true;
	    } else {
	      returnBool = false;
	    }
	    
	  } else { // everything is the same
	    
	    returnBool = false;
	    
	  }
	  
	} //if nbjetsmin
	
      } else {
	
	returnBool = thisNJmax<rhsNJmax;
	
      } // if njetsmax
      
    } // if nbjetsmax
    else {
      if( nBJetsMax < 0 )
	returnBool = false;
      else returnBool = true;
      
    }
  }
  
  return returnBool;
  
}




bool MT2SignalRegion::operator>( const MT2SignalRegion& rhs ) const {
  
  return !( (*this) <= rhs );

}


bool MT2SignalRegion::operator>=( const MT2SignalRegion& rhs ) const {

  return ( (*this) > rhs || (*this) == rhs );

}


bool MT2SignalRegion::operator<=( const MT2SignalRegion& rhs ) const {

  return ( (*this) < rhs || (*this) == rhs );

}


bool MT2SignalRegion::isIncluded( const MT2SignalRegion* sigRegion ) const {

  bool returnBool = true;

  if( nJetsMin < sigRegion->nJetsMin ) returnBool = false;
  if( nJetsMax > sigRegion->nJetsMax && sigRegion->nJetsMax>=0 ) returnBool = false;
  if( nJetsMax < sigRegion->nJetsMax && nJetsMax < 0 ) returnBool = false; //To solve case where nJetsMax is inf and sigRegion->nJetsMax is finite
  if( nBJetsMin < sigRegion->nBJetsMin ) returnBool = false;
  if( nBJetsMax > sigRegion->nBJetsMax && sigRegion->nBJetsMax>=0 ) returnBool = false;
  if( nBJetsMax < sigRegion->nBJetsMax && nBJetsMax < 0 ) returnBool = false; //To solve case where nBJetsMax is inf and sigRegion->nBJetsMax is finite 
  if( sigRegion->mtCut != "" && mtCut != sigRegion->mtCut ) return false;


  return returnBool;

}








////////////////////////////////////////////////////////
//
//                  MT2Region
//
////////////////////////////////////////////////////////



MT2Region::MT2Region( const std::string& regionName ) {

  // this constructor parses the name
  // the name has to be passed in the format:
  // [htRegionName]_[signalRegionName]

  std::stringstream ss(regionName);
  std::vector<std::string> parts;
  std::string item;
  while(std::getline(ss, item, '_')) {
    parts.push_back(item);
  }


  if( parts.size()==0 ) {
    std::cout << "[MT2Region]::MT2Region ERROR! Unrecognized MT2Region name: " << regionName << std::endl;
    exit(459);
  }


  std::string htRegionName = parts[0];
  std::string signalRegionName = "";
  for( unsigned i=1; i<parts.size(); ++i ) {
    if( i==1 ) signalRegionName = parts[i] + "_";
    else if( i==parts.size()-1 ) signalRegionName += parts[i];
    else signalRegionName = signalRegionName + parts[i] + "_";
  }

  htRegion_ = new MT2HTRegion( htRegionName );
  sigRegion_ = new MT2SignalRegion( signalRegionName );

}

  




void MT2Region::getBins( int &nBins, double*& bins) const {

  std::string regionName = this->getName();



  if( regionName == "HT450toInf_j2toInf_b0toInf" ) {  // this is the inclusive region

    //    const int nBins_tmp                        = 7;
    //    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 600., 800., 1000. };
    //    const int nBins_tmp                        = 5;
    //    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 1000., 1500.};
    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 1000.};
    //    const int nBins_tmp                        = 3;
    //    bins = new double[nBins_tmp+1]{200., 300., 400., 600.};
    nBins = nBins_tmp;

  } 
  // Plot
  else if( regionName == "HT200to1000_j2to3_b0" ){
    
    const int nBins_tmp                        = 36;
    bins = new double[nBins_tmp+1]{200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT200to1000_j2to3_b1toInf" ){
    
    const int nBins_tmp                        = 36;
    bins = new double[nBins_tmp+1]{200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT200to1000_j4toInf_b0" ){
    
    const int nBins_tmp                        = 36;
    bins = new double[nBins_tmp+1]{200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT200to1000_j4toInf_b1toInf" ){
    
    const int nBins_tmp                        = 36;
    bins = new double[nBins_tmp+1]{200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.};
    nBins = nBins_tmp;
    
  }
  //
  else if( regionName == "HT1000toInf_j2to3_b0" ){
    
    const int nBins_tmp                        = 36;
    bins = new double[nBins_tmp+1]{200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT1000toInf_j2to3_b1toInf" ){
    
    const int nBins_tmp                        = 36;
    bins = new double[nBins_tmp+1]{200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT1000toInf_j4toInf_b0" ){
    
    const int nBins_tmp                        = 36;
    bins = new double[nBins_tmp+1]{200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT1000toInf_j4toInf_b1toInf" ){
    
    const int nBins_tmp                        = 36;
    bins = new double[nBins_tmp+1]{200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.};
    nBins = nBins_tmp;
    
  }  
  // Monojet
  else if( regionName == "HT200toInf_j1_b0toInf" ){ // monojet inclusive
    
    const int nBins_tmp                        = 7;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT200toInf_j1_b0" ){ // monojet 0b
    
    const int nBins_tmp                        = 7;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT200toInf_j1_b1toInf" ){ // monojet 1b
    
    const int nBins_tmp                        = 6;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 600., 800., 1500.};
    nBins = nBins_tmp;
    
  }

  else if( regionName == "HT200to250_j1_b0" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT250to350_j1_b0" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT350to450_j1_b0" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT450to575_j1_b0" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT575to700_j1_b0" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT700to1000_j1_b0" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT1000to1200_j1_b0" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }  else if( regionName == "HT1200toInf_j1_b0" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};   //NEW EOY 2016 binning
    nBins = nBins_tmp;
    
  }

  else if( regionName == "HT200to250_j1_b1toInf" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT250to350_j1_b1toInf" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT350to450_j1_b1toInf" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT450to575_j1_b1toInf" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT575to700_j1_b1toInf" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }
  else if( regionName == "HT700toInf_j1_b1toInf" ){ // monojet inclusive
    
    const int nBins_tmp                        = 1;
    bins = new double[nBins_tmp+1]{0., 1500.};
    nBins = nBins_tmp;
    
  }

  
  else if( regionName == "HT200to450_j2to3_b0" ){ // new MT2 binning
    
    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j2to3_b1" ){ // new MT2 binning
    
    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j2to3_b2" ){ // new MT2 binning
    
    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j4to6_b0" ){ // new MT2 binning
    
    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j4to6_b1" ){ // new MT2 binning
    
    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j4to6_b2" ){ // new MT2 binning
    
    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j7toInf_b0" ){ // new MT2 binning
    
    const int nBins_tmp                        = 2;
    bins = new double[nBins_tmp+1]{200., 300., 1500.};   //NEW EOY 2016 binning
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j7toInf_b1" ){ // new MT2 binning
    
    const int nBins_tmp                        = 2;
    bins = new double[nBins_tmp+1]{200., 300., 1500.};   //NEW EOY 2016 binning
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j7toInf_b2" ){ // new MT2 binning
    
    const int nBins_tmp                        = 2;
    bins = new double[nBins_tmp+1]{200., 300., 1500.};   //NEW EOY 2016 binning
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j2to6_b3toInf" ){ // new MT2 binning
    
    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};   //NEW EOY 2016 binning
    nBins = nBins_tmp;
    
  } 
  else if( regionName == "HT200to450_j7toInf_b3toInf" ){ // new MT2 binning
    
    const int nBins_tmp                        = 2;
    bins = new double[nBins_tmp+1]{200., 300., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;
    
  } 
  //////End of mono-jet
  
//else if( regionName == "HT450to575_j2toInf_b0toInf" ) {  // inclusive low HT trigger region

//  const int nBins_tmp                        = 50;
//  bins = new double[nBins_tmp+1];
//  for( unsigned i=0; i<nBins_tmp+1; ++i ) 
//    bins[i] = (double)i*10.;
//  nBins = nBins_tmp;

//} 

//else if( regionName == "HT575to1000_j2toInf_b0toInf" ) {  // inclusive medium HT trigger region

//  const int nBins_tmp                        = 50;
//  bins = new double[nBins_tmp+1];
//  for( unsigned i=0; i<nBins_tmp+1; ++i ) 
//    bins[i] = (double)i*10.;
//  nBins = nBins_tmp;

//} 

//else if( regionName == "HT575to900_j2toInf_b0toInf" ) {  // modified inclusive medium HT trigger region

//  const int nBins_tmp                        = 50;
//  bins = new double[nBins_tmp+1];
//  for( unsigned i=0; i<nBins_tmp+1; ++i ) 
//    bins[i] = (double)i*10.;
//  nBins = nBins_tmp;

//} 

//else if( regionName == "HT1000toInf_j2toInf_b0toInf" ) {  // inclusive high HT trigger region

//  const int nBins_tmp                        = 50;
//  bins = new double[nBins_tmp+1];
//  for( unsigned i=0; i<nBins_tmp+1; ++i ) 
//    bins[i] = (double)i*10.;
//  nBins = nBins_tmp;

//} 

//else if( regionName == "HT900toInf_j2toInf_b0toInf" ) {  // modified inclusive high HT trigger region

//  const int nBins_tmp                        = 50;
//  bins = new double[nBins_tmp+1];
//  for( unsigned i=0; i<nBins_tmp+1; ++i ) 
//    bins[i] = (double)i*10.;
//  nBins = nBins_tmp;

//} 

////HERE (is the important stuff)

  else if( regionName == "HT450to575_j2to3_b0" ){ // new MT2 binning
    
    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;
    
  } else if( regionName == "HT450to575_j2to3_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j4to6_b0" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j4to6_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};   //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j7toInf_b0" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};   //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j7toInf_b1" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};   //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j2to3_b2" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};  //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j4to6_b2" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};  //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j7toInf_b2" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 300., 400., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j2to6_b3toInf" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j7toInf_b3toInf" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200.,300., 400., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j2to3_b0" ){ 

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j2to3_b1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};  //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j4to6_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};  //NEW EOY 2016 binning (same)
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j4to6_b1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j7toInf_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j7toInf_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j2to3_b2" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j4to6_b2" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j7toInf_b2" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j2to6_b3toInf" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 1500.};   //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j7toInf_b3toInf" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j2to3_b0" ){

    const int nBins_tmp                        = 6;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1200., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j2to3_b1" ){

    const int nBins_tmp                        = 6;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1200., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j4to6_b0" ){

    const int nBins_tmp                        = 6;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1200., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j4to6_b1" ){

    const int nBins_tmp                        = 6;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1200., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j7toInf_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j7toInf_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j2to3_b2" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j4to6_b2" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j7toInf_b2" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1500.}; //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j2to6_b3toInf" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 400., 600., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j7toInf_b3toInf" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 400., 600., 1500.}; //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j2to3_b0" ){

    const int nBins_tmp                        = 6;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1400., 1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j2to3_b1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.}; //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j4to6_b0" ){

    const int nBins_tmp                        = 6;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1400., 1800.}; //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j4to6_b1" ){

    const int nBins_tmp                        = 6;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1400., 1800.}; //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j7toInf_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};//NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j7toInf_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1500.};//NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j2to3_b2" ){

    const int nBins_tmp                        = 2;
    bins = new double[nBins_tmp+1]{200., 400.,  1500.};  //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j4to6_b2" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800.,1500.};//NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j7toInf_b2" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1500.};//NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j2to6_b3toInf" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 400., 600., 1500.}; //NEW EOY 2016 binning
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j7toInf_b3toInf" ){

    const int nBins_tmp                        = 2;
    bins = new double[nBins_tmp+1]{200., 400., 1500.};
    nBins = nBins_tmp;

  } 












  ////// DARK MATTER
  else if( regionName == "HT450to575_j2_b0" ){
    
    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;
    
  } else if( regionName == "HT450to575_j2_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j2_b0toInf" ){
    
    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;
    
  } else if( regionName == "HT450to575_j2_b0to1" ){
    
    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;
    
  } else if( regionName == "HT450to575_j3_b0" ){ 
    
    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;
    
  } else if( regionName == "HT450to575_j3_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j3_b0to1" ){ 
    
    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;
    
  } else if( regionName == "HT450to575_j2_b0toInf" ){
    
    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;
    
  } 

  else if( regionName == "HT575to1000_j2_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j2_b1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j2_b0to1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j2_b0toInf" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j3_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j3_b1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j3_b0to1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j3_b0toInf" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } 

  else if( regionName == "HT1000to1500_j2_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j2_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j2_b0to1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j2_b0toInf" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j3_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j3_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j3_b0to1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j3_b0toInf" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  }

  else if( regionName == "HT1500toInf_j2_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j2_b1" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 400., 600., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j2_b0to1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j2_b0toInf" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j3_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j3_b1" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 400., 600., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j3_b0to1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j3_b0toInf" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } 

  else if( regionName == "HT450to575_j4toInf_b0" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j4toInf_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j4toInf_b0to1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT450to575_j4toInf_b0toInf" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 1500.};
    nBins = nBins_tmp;

  } 

  else if( regionName == "HT575to1000_j4toInf_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j4toInf_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j4toInf_b0to1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT575to1000_j4toInf_b0toInf" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 300., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } 
  
  else if( regionName == "HT1000to1500_j4toInf_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j4toInf_b1" ){

    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j4toInf_b0to1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1000to1500_j4toInf_b0toInf" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } 
  
  else if( regionName == "HT1500toInf_j4toInf_b0" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j4toInf_b1" ){

    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{200., 400., 600., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j4toInf_b0to1" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

  } else if( regionName == "HT1500toInf_j4toInf_b0toInf" ){

    const int nBins_tmp                        = 5;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1000., 1500.};
    nBins = nBins_tmp;

   
  } else if( regionName == "HT450toInf_j1_b0toInf" ){ // monojet old inclusive

//    const int nBins_tmp                        = 3;
//    bins = new double[nBins_tmp+1]{450., 575., 1000., 1500.};
    const int nBins_tmp                        = 7;
    bins = new double[nBins_tmp+1]{200., 250., 300., 400., 450., 575., 1000., 1500.};
    nBins = nBins_tmp;

  } 
  //////END DARKMATTER

  else { // default binning
 
    //    const int nBins_tmp                        = 7;
    //    bins = new double[nBins_tmp+1]{200., 300., 400., 500., 600., 800., 1000., 1500. };
    const int nBins_tmp                        = 4;
    bins = new double[nBins_tmp+1]{200., 400., 600., 800., 1500. };
//    const int nBins_tmp                        = 9;
//    bins = new double[nBins_tmp+1]{0., 100., 200., 300., 400., 500., 600., 800., 1000., 1500. };
    //const int nBins_tmp                        = 4;
    //bins = new double[nBins_tmp+1]{200., 300., 400., 600., 1000.};
    //const int nBins_tmp                        = 5;
    //bins = new double[nBins_tmp+1]{200., 300., 400., 600., 1000., 1500.};
    nBins = nBins_tmp;

  }


}

void MT2Region::getBins_qcdCR( int &nBins, double*& bins) const {
    const int nBins_tmp = 18;
    bins = new double[nBins_tmp+1]{40,45,50,55,60,65,70,75,80,85,90,95,100,125,200,300,450,800, 1500}; //nBins_tmp = 18;
    //bins = new double[nBins_tmp+1]{40,45,50,55,60,65,70,75,80,88,100,125,200,300,450,800, 1500}; //nBins_tmp = 16;
    //bins = new double[nBins_tmp+1]{40,45,50,55,60,65,70,75,80,88,100,125,180,250,450,800};
    //bins = new double[nBins_tmp+1]{30,35,40,45,50,55,60,65,70,75,80,88,100,125,180,250,450,800};
    //bins = new double[nBins_tmp+1]{30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200};
    nBins = nBins_tmp;
}

std::vector< std::string > MT2Region::getNiceNames() const {

  std::vector< std::string > names;
  names.push_back(htRegion_->getNiceName());
  names.push_back(sigRegion_->getNiceName());

  return names;

}

std::vector< std::string > MT2Region::getNiceNamesLatex() const {

  std::vector< std::string > names;
  names.push_back(htRegion_->getNiceNameLatex());
  names.push_back(sigRegion_->getNiceNameLatex());

  return names;

}


std::string MT2Region::getBinName(double& min, double& max) const {

  char binName[500];
  if( max < 0 )
    sprintf( binName, "M_{T2} > %.0f GeV", min );
  else
    sprintf( binName, "%.0f < M_{T2} < %.0f GeV", min, max );
  
  std::string binName_str(binName);
  
  return binName_str;
  
}

std::string MT2Region::getBinNameLatex(double& min, double& max) const {

  char binName[500];
  if( max < 0 )
    sprintf( binName, "M$_{\\mathrm{T2}} > %.0f$~GeV", min );
  else
    sprintf( binName, "$%.0f < $M$_{\\mathrm{T2}} < %.0f$~GeV", min, max );
  
  std::string binName_str(binName);
  
  return binName_str;
  
}

std::vector< std::string> MT2Region::getBinNames() const {
  
  int nBins;
  double* bins;
  this->getBins(nBins, bins);
  
  std::vector< std::string > names;
  for( int i=1; i<nBins; ++i )
    names.push_back( getBinName(bins[i-1], bins[i]) );
  double lastbin=-1.;
  names.push_back( getBinName(bins[nBins-1], lastbin) );

  return names;

}

std::vector< std::string> MT2Region::getBinNamesLatex() const {
  
  int nBins;
  double* bins;
  this->getBins(nBins, bins);
  
  std::vector< std::string > names;
  for( int i=1; i<nBins; ++i )
    names.push_back( getBinNameLatex(bins[i-1], bins[i]) );
  double lastbin=-1.;
  names.push_back( getBinNameLatex(bins[nBins-1], lastbin) );

  return names;

}


std::string MT2Region::getRegionCuts() const {

  std::string sigCuts = sigRegion_->getCuts();
  std::string htCuts  = htRegion_->getCuts();

  std::string cuts = sigCuts + " && " + htCuts;
  return cuts;

}


bool MT2Region::isIncluded( const MT2Region* region ) const {

  return ( ( sigRegion_->isIncluded( region->sigRegion() ) ) && ( htRegion_->isIncluded( region->htRegion() ) ) );

}


bool MT2Region::operator==( const MT2Region& rhs ) const {

  return ( (*htRegion_)==(*(rhs.htRegion())) && (*sigRegion_)==(*(rhs.sigRegion())) );

}




bool MT2Region::operator!=( const MT2Region& rhs ) const {

  return !( *this == rhs );

}




bool MT2Region::operator<( const MT2Region& rhs ) const {

  if( (sigRegion_->nJetsMax > 1 || sigRegion_->nJetsMax < 0) && (rhs.sigRegion()->nJetsMax > 1 || rhs.sigRegion()->nJetsMax < 0)){ 
    if( (*htRegion_)!=(*(rhs.htRegion())) ) {
      return (*htRegion_)<(*(rhs.htRegion()));
    } else {
      return (*sigRegion_)<(*(rhs.sigRegion()));
    }
  }
  else{
    
    if( (*sigRegion_)!=(*(rhs.sigRegion())) ) {
      return (*sigRegion_)<(*(rhs.sigRegion()));
    } else {
      return (*htRegion_)<(*(rhs.htRegion()));
    }

  }

  return true;

}




bool MT2Region::operator>( const MT2Region& rhs ) const {

    return !(*this <= rhs);

}




bool MT2Region::operator>=( const MT2Region& rhs ) const {

  return (*this > rhs || *this == rhs);

}




bool MT2Region::operator<=( const MT2Region& rhs ) const {

  return (*this < rhs || *this == rhs);

}



//  LocalWords:  tmp





