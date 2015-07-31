#ifndef MT2DrawTools_h
#define MT2DrawTools_h

#include "TStyle.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
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

 private:

};

#endif
