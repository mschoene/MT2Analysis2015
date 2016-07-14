#ifndef MT2DrawTools_h
#define MT2DrawTools_h

#include "TStyle.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TColor.h"
#include "TMatrixD.h"

#include "../interface/MT2Analysis.h"



#define kQCD 401
#define kWJets 417
#define kZJets 419
#define kTop 855

#define kQCDest 402
#define kZinv 430
#define kLostLepton 418




class MT2EstimateTree;


class MT2DrawTools {

 public:

  MT2DrawTools( const std::string& outputdir="plots_tmp", float lumi=0. );
  //MT2DrawTools( const MT2Config& cfg );

  void set_outDir( const std::string& outdir );
  void set_data( MT2Analysis<MT2EstimateTree>* data );
  void set_mc( std::vector< MT2Analysis<MT2EstimateTree>* >* mc );
  void set_lumi( float lumi );
  void set_lumiErr( float lumiErr );
  void set_shapeNorm( bool shapeNorm );
  void set_mcSF( float mcsf );
  void set_addOverflow( bool addOver );
  void set_displaySF( bool displaySF );
  void set_doPaperPlots( bool doPaperPlots );

  bool twoPads() const;

  static TStyle* setStyle();

  static void addLabels( TCanvas* c1, float lumi, const std::string& text="CMS Preliminary");

  static TPaveText* getLabelTop( float lumi );
  static TPaveText* getLabelTop( float lumi2015, float lumi2016 );
  static TPaveText* getLabelTopSimulation( float lumi );
  static TPaveText* getLabelCMS( const std::string& text="CMS" );
  static TPaveText* getLabelTop( const std::string& text="CMS Preliminary, #sqrt{s} = 13 TeV" );
  static TPaveText* getLabelTopSimulation( const std::string& text="CMS Simulation, #sqrt{s} = 13 TeV" );

  static std::string getLumiText( float lumi );

  static TGraphAsymmErrors* getPoissonGraph( TH1D* h1, bool drawZeros=true, const std::string& xerrType="0", float nSigma=1. );
  static TGraphAsymmErrors* getRatioGraph( TH1D* h1, TH1D* h2, const std::string& xerrType="0" );

  static float getDataMCSF( TCanvas* c1 );
  static float graphIntegral( TGraphAsymmErrors* graph, float xMin = -99999., float xMax=999999. );
  static TList* getCorrectList( TCanvas* c1 );
  
  static TPad* getCanvasMainPad( bool logY=false );
  static TPad* getCanvasRatioPad( bool logY=false );
  static TH2D* getRatioAxes( float xMin, float xMax, float yMin=0., float yMax=2. );
  
  static TPaveText*  getRatioText( double integral_data, double integral_mc, double error_datamc );
  static TPaveText*  getFitText( TF1* f );
  
  static double getSFError(double integral_data, double error_data, double integral_mc, double error_mc);
  static TLine* getSFLine(double integral_data, double integral_mc, float xMin, float xMax);
  static TGraphErrors* getSFBand(double integral_data, double error_data, double integral_mc, double error_mc, float xMin, float xMax);

  static TF1* getSFFit(TGraphAsymmErrors* g_ratio, float xMin, float xMax);
  static void getSFFitParameters(TF1* f, double &sf, double &sfErr, double &chi2, int &ndof);
  static TGraphErrors* getSFFitBand(TF1* f, float xMin, float xMax);

  static TGraphErrors* getSystBand(float xMin, float xMax, double SystErr=0.0);
  static TH1D* getMCBandHisto( TH1D* histo_mc, double SystErr=0.0 );
  static TH1D* getBandAtOne( TH1D* h );

  static void addOverflowSingleHisto( TH1D* yield );
  static void addOverflowSingleHisto( TH3D* yield3d );

  static TH1D* getBand(TF1* f, const std::string& name ); // the TH1D then needs to be draw with the option "C E3"
  static TH1D* getBand(TF1 *f, TMatrixD const& m, std::string name, bool getRelativeBand=false, int npx=100); 


  std::vector<TCanvas*> drawRegionYields_fromTree( const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="", const std::string& kinCuts="", const std::string& topoCuts="" );
  //void drawRegionYields_fromTree( MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName="", const std::string& units="", const std::string& cutsLabel="" );

 private:

  std::string outdir_;
  float lumi_;
  float lumiErr_;
  bool shapeNorm_;
  bool addOverflow_;
  bool displaySF_;
  bool doPaperPlots_;

  MT2Analysis<MT2EstimateTree>* data_;
  std::vector< MT2Analysis<MT2EstimateTree>* >* mc_;
  float mcSF_;
  

};

#endif
