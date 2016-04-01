#include "../interface/MT2EstimateSyst.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TEfficiency.h"




void MT2EstimateSyst::rebinYields( MT2Analysis<MT2EstimateSyst>* analysis, int nBins, double* bins) {

  std::set<MT2Region> regions = analysis->getRegions();

  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateSyst* estimate = analysis->get(*iR);
    TH1D* thisYield = estimate->yield;

    std::string oldName(thisYield->GetName());
    delete thisYield;
    thisYield = new TH1D( oldName.c_str(), "", nBins, bins );
  

    TH1D* this_systUp = estimate->yield_systUp;
    TH1D* this_systDown = estimate->yield_systDown;


    std::string oldName_systUp(this_systUp->GetName());
    std::string oldName_systDown(this_systDown->GetName());
 
    delete this_systUp;
    this_systUp = new TH1D( oldName_systUp.c_str(), "", nBins, bins );
    this_systUp->Sumw2();

    delete this_systDown;
    this_systDown = new TH1D( oldName_systDown.c_str(), "", nBins, bins );
    this_systDown->Sumw2();
    
  }
  
}





MT2EstimateSyst::MT2EstimateSyst( const std::string& aname, const MT2Region& aregion ) : MT2Estimate( aname, aregion ) {

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

  yield_systUp = new TH1D( this->getHistoName("yield_systUp").c_str(), "", nBins, bins);
  yield_systUp->Sumw2();
  yield_systDown = new TH1D( this->getHistoName("yield_systDown").c_str(), "", nBins, bins);
  yield_systDown->Sumw2();

}




MT2EstimateSyst::MT2EstimateSyst( const std::string& aname, const MT2Region& aregion, const MT2Estimate& pass, const MT2Estimate& tot ) : MT2Estimate( aname, aregion ) {

  //get the right bins from the input files
  int nBins;
  double* bins;
  pass.getYieldBins(nBins, bins);

  //Rebinning the yield
  TH1D* thisYield = this->yield;
  std::string oldName(thisYield->GetName());
  if( thisYield!=0 )
    delete thisYield;
  thisYield = new TH1D( oldName.c_str(), "", nBins, bins );
  
 
  yield_systUp = new TH1D( this->getHistoName("yield_systUp").c_str(), "", nBins, bins);
  yield_systUp->Sumw2();
  yield_systDown = new TH1D( this->getHistoName("yield_systDown").c_str(), "", nBins, bins);
  yield_systDown->Sumw2();


  TEfficiency eff( *(pass.yield), *(tot.yield) );

  for( int i=1; i<this->yield->GetNbinsX()+1; ++i ) {

    yield         ->SetBinContent( i, eff.GetEfficiency(i) );
    yield_systDown->SetBinContent( i, eff.GetEfficiency(i) - eff.GetEfficiencyErrorLow(i) );
    yield_systUp  ->SetBinContent( i, eff.GetEfficiency(i) + eff.GetEfficiencyErrorUp(i) );

  }


}




MT2EstimateSyst::MT2EstimateSyst( const MT2Estimate& rhs ) : MT2Estimate(rhs) {

  this->yield_systUp = new TH1D(*(rhs.yield));
  this->yield_systUp   ->SetName(this->getHistoName("yield_systUp").c_str()); 

  this->yield_systDown = new TH1D(*(rhs.yield));
  this->yield_systDown ->SetName(this->getHistoName("yield_systDown").c_str()); 

}




MT2EstimateSyst::MT2EstimateSyst( const MT2EstimateSyst& rhs ) : MT2Estimate(rhs) {

  this->yield_systUp = new TH1D(*(rhs.yield_systUp));
  this->yield_systDown = new TH1D(*(rhs.yield_systDown));

}



MT2EstimateSyst::~MT2EstimateSyst() {

  delete yield_systUp;
  delete yield_systDown;

}




MT2Analysis<MT2EstimateSyst>* MT2EstimateSyst::makeEfficiencyAnalysis( const std::string& aname, MT2Analysis<MT2Estimate>* pass, MT2Analysis<MT2Estimate>* all ) {

  std::set<MT2Region> regions = pass->getRegions();

  std::set<MT2EstimateSyst*> data;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate* thisPass = pass->get( *iR );
    MT2Estimate* thisAll  = all ->get( *iR );

    MT2EstimateSyst* thisEff  = new MT2EstimateSyst( aname, *iR, *thisPass, *thisAll );
    data.insert( thisEff );

  } // for regions


  MT2Analysis<MT2EstimateSyst>* analysis = new MT2Analysis<MT2EstimateSyst>( aname, data );

  return analysis;

}



MT2Analysis<MT2EstimateSyst>* MT2EstimateSyst::makeAnalysisFromEstimate( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2Estimate>* estimate ) {

  std::set<MT2Region> regions = estimate->getRegions();

  std::set<MT2EstimateSyst*> data;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2Estimate*  thisEstimate = estimate->get( *iR );
    MT2EstimateSyst* thisEstimateSyst = new MT2EstimateSyst( *thisEstimate );
    data.insert( thisEstimateSyst );

  } // for regions


  MT2Analysis<MT2EstimateSyst>* analysis = new MT2Analysis<MT2EstimateSyst>( aname, data );

  return analysis;

}


MT2Analysis<MT2EstimateSyst>* MT2EstimateSyst::makeIntegralAnalysisFromEstimate( const std::string& aname, const std::string& regionsSet, MT2Analysis<MT2EstimateSyst>* estimate ) {

  std::set<MT2Region> regions = estimate->getRegions();

  std::set<MT2EstimateSyst*> data;

  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateSyst*  thisEstimate = new MT2EstimateSyst( estimate->get( *iR )->getName(), *iR );
    MT2EstimateSyst*  tempEstimate = estimate->get( *iR ); 

    double error;
    double integral = tempEstimate->yield->IntegralAndError(1, tempEstimate->yield->GetNbinsX()+1, error);
    double integralUp = tempEstimate->yield_systUp->Integral(1, tempEstimate->yield->GetNbinsX()+1);
    double integralDown = tempEstimate->yield_systDown->Integral(1, tempEstimate->yield->GetNbinsX()+1);
    for( int iBin = 1; iBin < thisEstimate->yield->GetNbinsX()+1; ++iBin ) {
      thisEstimate->yield->SetBinContent(iBin, integral);
      thisEstimate->yield_systUp->SetBinContent(iBin, integralUp);
      thisEstimate->yield_systDown->SetBinContent(iBin, integralDown);

//  for( std::set<MT2Region>::iterator iR=regions.begin(); iR!=regions.end(); ++iR ) {
//
//    MT2EstimateSyst*  thisEstimate = estimate->get( *iR );
//    
//    int nBins = thisEstimate->yield_systUp->GetNbinsX();
//    
//    double error;
//    double integral = thisEstimate->yield->IntegralAndError(1, nBins+1, error);
//    double integralUp = thisEstimate->yield_systUp->Integral(1, nBins+1);
//    double integralDown = thisEstimate->yield_systDown->Integral(1, nBins+1);
//    
//    for( int iBin = 1; iBin < nBins+1; ++iBin ){
//      thisEstimate->yield->SetBinContent(iBin, integral);
//      thisEstimate->yield_systUp->SetBinContent(iBin, integralUp);
//      thisEstimate->yield_systDown->SetBinContent(iBin, integralDown);

      for( int jBin = 1; jBin < thisEstimate->yield3d->GetNbinsY()+1; ++jBin )
        for( int kBin = 1; kBin < thisEstimate->yield3d->GetNbinsZ()+1; ++kBin ){
          thisEstimate->yield3d->SetBinContent(iBin, jBin, kBin, integral);
          thisEstimate->yield3d->SetBinError(iBin, jBin, kBin, error);
        }
    }

    data.insert( thisEstimate );

  } // for regions                                                                                                                                                                                                               

  MT2Analysis<MT2EstimateSyst>* analysis = new MT2Analysis<MT2EstimateSyst>( aname, data );

  return analysis;

}



//void MT2EstimateSyst::rebinYields( MT2Analysis<MT2EstimateSyst>* analysis, int nBins, double* bins) {
//
//  std::set<MT2Region> regions = analysis->getRegions();
//
//  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {
//
//    MT2EstimateSyst* estimate = analysis->get(*iR);
//    TH1D* thisYield = estimate->yield;
//
//    std::string oldName(thisYield->GetName());
//    if( thisYield!=0 )
//      delete thisYield;
//    thisYield = new TH1D( oldName.c_str(), "", nBins, bins );
//  
//
//    TH1D* this_systUp = estimate->yield_systUp;
//    TH1D* this_systDown = estimate->yield_systDown;
//
//
//    std::string oldName_systUp(this_systUp->GetName());
//    std::string oldName_systDown(this_systDown->GetName());
// 
//    if( this_systUp!=0 )
//      delete this_systUp;
//    this_systUp = new TH1D( oldName_systUp.c_str(), "", nBins, bins );
//    this_systUp->Sumw2();
//
//    if( this_systDown!=0 )
//      delete this_systDown;
//    this_systDown = new TH1D( oldName_systDown.c_str(), "", nBins, bins );
//    this_systDown->Sumw2();
//    
//  }
//  
//}


TGraphAsymmErrors* MT2EstimateSyst::getGraph() const {

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(0);
  graph->SetName( this->getHistoName("grSyst").c_str() );


  for( int iBin=1; iBin<yield->GetNbinsX()+1; ++iBin ) {

    float x = yield->GetBinCenter(iBin);
    float x_minus = yield->GetBinLowEdge(iBin);
    float x_plus  = yield->GetBinLowEdge(iBin+1);
    float y = yield->GetBinContent(iBin);
    float y_plus  = yield_systUp  ->GetBinContent(iBin);
    float y_minus = yield_systDown->GetBinContent(iBin);
    
    int iPoint = iBin-1;
    graph->SetPoint( iPoint, x, y );

    graph->SetPointError( iPoint, x-x_minus, x_plus-x, y-y_minus, y_plus-y );

  }


  return graph;

}




void MT2EstimateSyst::setName( const std::string& newName ) {

  MT2Estimate::setName(newName);

  yield_systUp->SetName( this->getHistoName("yield_systUp").c_str() );
  yield_systDown->SetName( this->getHistoName("yield_systDown").c_str() );

}



void MT2EstimateSyst::finalize( ) {

  MT2Estimate::finalize();

}



void MT2EstimateSyst::getShit( TFile* file, const std::string& path ) {

  MT2Estimate::getShit(file, path);
  yield_systUp = (TH1D*)file->Get(Form("%s/%s", path.c_str(), yield_systUp->GetName()));
  yield_systDown = (TH1D*)file->Get(Form("%s/%s", path.c_str(), yield_systDown->GetName()));


}



void MT2EstimateSyst::write() const {

  MT2Estimate::write();
  yield_systUp->SetLineColor(kRed);
  yield_systDown->SetLineColor(kGreen);
  yield_systUp->Write();
  yield_systDown->Write();

}



void MT2EstimateSyst::print(const std::string& ofs){

  Int_t binXmin=1;
  Int_t binXmax=-1;

  Double_t integral = yield->Integral(binXmin, binXmax);
  Double_t integral_up = yield_systUp->Integral(binXmin, binXmax);
  Double_t integral_down = yield_systDown->Integral(binXmin, binXmax);
  

  std::ofstream ofs_file;
  ofs_file.open( ofs, std::ofstream::app );
  if(integral >= 10)
    ofs_file << std::fixed << std::setprecision(1);
  else if(integral < 10)
    ofs_file << std::fixed << std::setprecision(2);
  ofs_file << " & " << integral << " $^{+" << fabs(integral_up-integral) << "}_{-" << fabs(integral_down-integral) << "}$";

}


const MT2EstimateSyst& MT2EstimateSyst::operator=( const MT2EstimateSyst& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield_systUp = new TH1D(*(rhs.yield_systUp));
  this->yield_systDown = new TH1D(*(rhs.yield_systDown));

  this->setName( this->getName() );


  return *this;

}




const MT2EstimateSyst& MT2EstimateSyst::operator=( const MT2Estimate& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));
  this->yield3d = new TH3D(*(rhs.yield3d));
  this->yield_systUp = new TH1D(*(rhs.yield));
  this->yield_systDown = new TH1D(*(rhs.yield));

  this->setName( this->getName() );


  return *this;

}






MT2EstimateSyst MT2EstimateSyst::operator+( const MT2EstimateSyst& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSyst::operator*] ERROR! Can't multiply MT2EstimateSyst with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin+otherBin;
    float newErrUp = sqrt( thisErrUp*thisErrUp + otherErrUp*otherErrUp );
    float newErrDown = sqrt( thisErrDown*thisErrDown + otherErrDown*otherErrDown );

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );
    
    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}



MT2EstimateSyst MT2EstimateSyst::operator-( const MT2EstimateSyst& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSyst::operator*] ERROR! Can't multiply MT2EstimateSyst with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin-otherBin;
    float newErrUp = sqrt( thisErrUp*thisErrUp + otherErrUp*otherErrUp );
    float newErrDown = sqrt( thisErrDown*thisErrDown + otherErrDown*otherErrDown );

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}




MT2EstimateSyst MT2EstimateSyst::operator/( const MT2EstimateSyst& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSyst::operator*] ERROR! Can't multiply MT2EstimateSyst with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin/otherBin;
    //Uncertainty down influences the new uncertainty up and vice versa
    float newErrUp = sqrt( thisErrUp*thisErrUp/(otherBin*otherBin) + otherErrDown*otherErrDown*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );
    float newErrDown = sqrt( thisErrDown*thisErrDown/(otherBin*otherBin) + otherErrUp*otherErrUp*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );
    //float newErrUp = sqrt( thisErrUp*thisErrUp/(otherBin*otherBin) + otherErrUp*otherErrUp*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );
    //float newErrDown = sqrt( thisErrDown*thisErrDown/(otherBin*otherBin) + otherErrDown*otherErrDown*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}




MT2EstimateSyst MT2EstimateSyst::operator*( const MT2EstimateSyst& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSyst::operator*] ERROR! Can't multiply MT2EstimateSyst with different MT2Regions!" << std::endl;
    exit(113);
  }


  MT2EstimateSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = (thisBinUp - thisBin)/thisBin;
    float thisErrDown = (thisBin - thisBinDown)/thisBin;
    float otherErrUp = (otherBinUp - otherBin)/otherBin;
    float otherErrDown = (otherBin - otherBinDown)/otherBin;

    float newBin = thisBin*otherBin;
    float newErrUp = sqrt( otherErrUp*otherErrUp + thisErrUp*thisErrUp );
    float newErrDown = sqrt( otherErrDown*otherErrDown + thisErrDown*thisErrDown );

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp*newBin );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown*newBin );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );
    
  }

  return result;

}



MT2EstimateSyst MT2EstimateSyst::operator+( const MT2Estimate& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSyst::operator*] ERROR! Can't multiply MT2EstimateSyst with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin     = thisBin+otherBin;
    float newErrUp   = thisErrUp;
    float newErrDown = thisErrDown;

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}


MT2EstimateSyst MT2EstimateSyst::operator-( const MT2Estimate& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSyst::operator*] ERROR! Can't multiply MT2EstimateSyst with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin     = thisBin-otherBin;
    float newErrUp   = thisErrUp;
    float newErrDown = thisErrDown;

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for ( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}


MT2EstimateSyst MT2EstimateSyst::operator*( const MT2Estimate& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSyst::operator*] ERROR! Can't multiply MT2EstimateSyst with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin     = thisBin*otherBin;
    float newErrUp   = thisErrUp*otherBin;
    float newErrDown = thisErrDown*otherBin;

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}


MT2EstimateSyst MT2EstimateSyst::operator/( const MT2Estimate& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateSyst::operator*] ERROR! Can't multiply MT2EstimateSyst with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateSyst result(*this);

  for( int iBin=1; iBin<result.yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = result.yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = result.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = result.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin     = thisBin/otherBin;
    float newErrUp   = thisErrUp/otherBin;
    float newErrDown = thisErrDown/otherBin;

    result.yield         ->SetBinContent( iBin, newBin );
    result.yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    result.yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<result.yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<result.yield3d->GetNbinsZ()+1; ++iBinZ)
	result.yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return result;

}





MT2EstimateSyst MT2EstimateSyst::operator*( float k ) const{

  MT2EstimateSyst result(this->getName(), *(this->region) );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(k);
  
  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(k);

  result.yield_systUp = new TH1D(*(this->yield_systUp));
  result.yield_systUp->Scale(k);

  result.yield_systDown = new TH1D(*(this->yield_systDown));
  result.yield_systDown->Scale(k);

  return result;

}



MT2EstimateSyst MT2EstimateSyst::operator/( float k ) const{

  MT2EstimateSyst result(this->getName(), *(this->region) );
  result.yield = new TH1D(*(this->yield));
  result.yield->Scale(1./k);

  result.yield3d = new TH3D(*(this->yield3d));
  result.yield3d->Scale(1./k);

  result.yield_systUp = new TH1D(*(this->yield_systUp));
  result.yield_systUp->Scale(1./k);

  result.yield_systDown = new TH1D(*(this->yield_systDown));
  result.yield_systDown->Scale(1./k);

  return result;


}



const MT2EstimateSyst& MT2EstimateSyst::operator+=( const MT2EstimateSyst& rhs ) {

  this->yield->Add(rhs.yield);
  this->yield3d->Add(rhs.yield3d);
  this->yield_systUp->Add(rhs.yield_systUp);
  this->yield_systDown->Add(rhs.yield_systDown);
  return (*this);

}

const MT2EstimateSyst& MT2EstimateSyst::operator/=( const MT2EstimateSyst& rhs ) {


  for( int iBin=1; iBin<this->yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = this->yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = this->yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = this->yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin/otherBin;    
    //Uncertainty down influences the new uncertainty up and vice versa
    float newErrUp = sqrt( thisErrUp*thisErrUp/(otherBin*otherBin) + otherErrDown*otherErrDown*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );
    float newErrDown = sqrt( thisErrDown*thisErrDown/(otherBin*otherBin) + otherErrUp*otherErrUp*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );

    // float newErrUp = sqrt( thisErrUp*thisErrUp/(otherBin*otherBin) + otherErrUp*otherErrUp*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );
    // float newErrDown = sqrt( thisErrDown*thisErrDown/(otherBin*otherBin) + otherErrDown*otherErrDown*thisBin*thisBin/(otherBin*otherBin*otherBin*otherBin) );

    this->yield         ->SetBinContent( iBin, newBin );
    this->yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    this->yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<this->yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<this->yield3d->GetNbinsZ()+1; ++iBinZ)
        this->yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }

  return (*this);

}



const MT2EstimateSyst& MT2EstimateSyst::operator*=( const MT2EstimateSyst& rhs ) {


  for( int iBin=1; iBin<this->yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = this->yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = this->yield_systUp->GetBinContent(iBin);
    float otherBinUp = rhs.yield_systUp->GetBinContent(iBin);
    float thisBinDown  = this->yield_systDown->GetBinContent(iBin);
    float otherBinDown = rhs.yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;
    float otherErrUp = otherBinUp - otherBin;
    float otherErrDown = otherBin - otherBinDown;

    float newBin = thisBin*otherBin;
    float newErrUp = sqrt( thisBin*thisBin*otherErrUp*otherErrUp + otherBin*otherBin*thisErrUp*thisErrUp );
    float newErrDown = sqrt( thisBin*thisBin*otherErrDown*otherErrDown + otherBin*otherBin*thisErrDown*thisErrDown );

    this->yield         ->SetBinContent( iBin, newBin );
    this->yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    this->yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<this->yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<this->yield3d->GetNbinsZ()+1; ++iBinZ)
        this->yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }


  return (*this);

}



const MT2EstimateSyst& MT2EstimateSyst::operator*=( const MT2Estimate& rhs ) {

  for( int iBin=1; iBin<this->yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = this->yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = this->yield_systUp->GetBinContent(iBin);
    float thisBinDown  = this->yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin = thisBin*otherBin;
    float newErrUp   = thisErrUp*otherBin;
    float newErrDown = thisErrDown*otherBin;

    this->yield         ->SetBinContent( iBin, newBin );
    this->yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    this->yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<this->yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<this->yield3d->GetNbinsZ()+1; ++iBinZ)
        this->yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );

  }


  return (*this);

}



const MT2EstimateSyst& MT2EstimateSyst::operator/=( const MT2Estimate& rhs ) {

  for( int iBin=1; iBin<this->yield->GetNbinsX()+1; ++iBin ) {

    float thisBin  = this->yield->GetBinContent(iBin);
    float otherBin = rhs.yield->GetBinContent(iBin);

    float thisBinUp  = this->yield_systUp->GetBinContent(iBin);
    float thisBinDown  = this->yield_systDown->GetBinContent(iBin);

    float thisErrUp = thisBinUp - thisBin;
    float thisErrDown = thisBin - thisBinDown;

    float newBin = thisBin/otherBin;
    float newErrUp   = thisErrUp/otherBin;
    float newErrDown = thisErrDown/otherBin;

    this->yield         ->SetBinContent( iBin, newBin );
    this->yield_systUp  ->SetBinContent( iBin, newBin + newErrUp );
    this->yield_systDown->SetBinContent( iBin, newBin - newErrDown );

    for( int iBinY=1; iBinY<this->yield3d->GetNbinsY()+1; ++iBinY)
      for( int iBinZ=1; iBinZ<this->yield3d->GetNbinsZ()+1; ++iBinZ)
        this->yield3d   ->SetBinContent( iBin, iBinY, iBinZ, newBin );
    
  }


  return (*this);

}



const MT2EstimateSyst& MT2EstimateSyst::operator*=( float k ) {

  this->yield->Scale(k);
  this->yield3d->Scale(k);
  this->yield_systUp->Scale(k);
  this->yield_systDown->Scale(k);
  return (*this);

}

const MT2EstimateSyst& MT2EstimateSyst::operator/=( float k ) {

  this->yield->Scale(1./k);
  this->yield3d->Scale(1./k);
  this->yield_systUp->Scale(1./k);
  this->yield_systDown->Scale(1./k);
  return (*this);

}




// friend functions


MT2EstimateSyst operator*( float k, const MT2EstimateSyst& rhs ) {

  return rhs*k;

}


MT2EstimateSyst operator/( float k, const MT2EstimateSyst& rhs ) {

  return rhs/k;

}


