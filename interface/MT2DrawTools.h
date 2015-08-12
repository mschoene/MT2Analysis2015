#ifndef MT2DrawTools_h
#define MT2DrawTools_h

#include "TStyle.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TColor.h"



class MT2DrawTools {

 public:

  static TStyle* setStyle();

  static TPaveText* getLabelTop( float lumi );
  static TPaveText* getLabelTopSimulation( float lumi );
  static TPaveText* getLabelTop( const std::string& text="CMS Preliminary, #sqrt{s} = 13 TeV" );
  static TPaveText* getLabelTopSimulation( const std::string& text="CMS Simulation, #sqrt{s} = 13 TeV" );

  static std::string getLumiText( float lumi );

  static TGraphAsymmErrors* getPoissonGraph( TH1D* h1, bool drawZeros=true, const std::string& xerrType="0", float nSigma=1. );
  static TGraphAsymmErrors* getRatioGraph( TH1D* h1, TH1D* h2 );
  
  static TPad* getCanvasMainPad( bool logY=false );
  static TPad* getCanvasRatioPad( bool logY=false );
  static TH2D* getRatioAxes( float xMin, float xMax, float yMin=0., float yMax=2. );
  
  static TPaveText*  getRatioText( double integral_data, double integral_mc, double error_datamc );
  
  static double getSFError(double integral_data, double error_data, double integral_mc, double error_mc);
  static TLine* getSFLine(double integral_data, double integral_mc, float xMin, float xMax);
  static TGraphErrors* getSFBand(double integral_data, double error_data, double integral_mc, double error_mc, float xMin, float xMax);

  static TGraphErrors* getSystBand(float xMin, float xMax, double SystErr=0.0);
  static TH1D* getMCBandHisto( TH1D* histo_mc, double SystErr=0.0 );

  static void addOverflowSingleHisto( TH1D* yield );

 private:

};

#endif
