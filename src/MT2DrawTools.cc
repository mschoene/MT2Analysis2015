#include "../interface/MT2DrawTools.h"

#include "RooHistError.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMinuit.h"

#include "../interface/MT2EstimateTree.h"



MT2DrawTools::MT2DrawTools( const std::string& outDir, float lumi ) {

  lumi_    = lumi;
  lumiErr_ = 0.026;
  shapeNorm_ = false;
  outdir_ = outDir;

  data_ = 0;
  data2_ = 0;
  mc_ = 0;

  mcSF_ = 1.;
  mcSF2_ = 1.;

  addOverflow_ = true;

  displaySF_ = true;

  doPaperPlots_ = false;

  std::cout << "[MT2DrawTools] Initiating: " << std::endl;
  std::cout << "     lumi: " << lumi_ << std::endl;
  std::cout << "     lumiErr: " << lumiErr_ << std::endl;
  std::cout << "     shapeNorm: " << shapeNorm_ << std::endl;
  std::cout << "     mcSF: " << mcSF_ << std::endl;
  std::cout << "     outDir: " << outdir_ << std::endl;
  std::cout << "     doPaperPlots: " << doPaperPlots_ << std::endl;

}


void MT2DrawTools::set_data( MT2Analysis<MT2EstimateTree>* data ) {

  data_ = data;

}

void MT2DrawTools::set_data2( MT2Analysis<MT2EstimateTree>* data2 ) {

  data2_ = data2;

}


void MT2DrawTools::set_mc( std::vector< MT2Analysis<MT2EstimateTree>* >* mc ) {

  mc_ = mc;

}


void MT2DrawTools::set_outDir( const std::string& outdir ) {

  std::cout << "[MT2DrawTools] Setting outdir to: " << outdir << std::endl;
  system( Form("mkdir -p %s", outdir.c_str()) );
  outdir_ = outdir;

}


void MT2DrawTools::set_lumi( float lumi) { 

  std::cout << "[MT2DrawTools] Setting lumi to: " << lumi<< std::endl;
  lumi_ = lumi; 

}


void MT2DrawTools::set_lumiErr( float lumiErr ) { 

  std::cout << "[MT2DrawTools] Setting lumi error to: " << lumiErr << std::endl;
  lumiErr_ = lumiErr; 

}


void MT2DrawTools::set_shapeNorm( bool shapeNorm ) { 

  if( shapeNorm )
    std::cout << "[MT2DrawTools] Using shape normalization." << std::endl;
  else
    std::cout << "[MT2DrawTools] Using lumi normalization." << std::endl;
  shapeNorm_ = shapeNorm;

}


void MT2DrawTools::set_addOverflow( bool addOver ) { 

  if( addOver )
    std::cout << "[MT2DrawTools] Adding overflow bins." << std::endl;
  else
    std::cout << "[MT2DrawTools] Disabled adding overflow bins." << std::endl;
  addOverflow_ = addOver;

}



void MT2DrawTools::set_displaySF( bool displaySF ) { 

  if( displaySF )
    std::cout << "[MT2DrawTools] Setting display SF: ON" << std::endl;
  else
    std::cout << "[MT2DrawTools] Setting display SF: OFF" << std::endl;
  displaySF_ = displaySF;

}



void MT2DrawTools::set_mcSF( float mcsf ) {

  std::cout << "[MT2DrawTools] Setting MC SF to: " << mcsf << std::endl;
  mcSF_ = mcsf;

}

void MT2DrawTools::set_mcSF2( float mcsf2 ) {

  std::cout << "[MT2DrawTools] Setting MC SF to: " << mcsf2 << std::endl;
  mcSF2_ = mcsf2;

}


TStyle* MT2DrawTools::setStyle() {

  // set the TStyle
  TStyle* style = new TStyle("DrawBaseStyle", "");
  style->SetCanvasColor(0);
  style->SetPadColor(0);
  style->SetFrameFillColor(0);
  style->SetStatColor(0);
  style->SetOptStat(0);
  style->SetTitleFillColor(0);
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);
  style->SetFrameBorderMode(0);
  style->SetPadBottomMargin(0.12);
  style->SetPadLeftMargin(0.12);
  style->cd();
  // For the canvas:
  style->SetCanvasBorderMode(0);
  style->SetCanvasColor(kWhite);
  style->SetCanvasDefH(600); //Height of canvas
  style->SetCanvasDefW(600); //Width of canvas
  style->SetCanvasDefX(0); //POsition on screen
  style->SetCanvasDefY(0);
  // For the Pad:
  style->SetPadBorderMode(0);
  style->SetPadColor(kWhite);
  style->SetPadGridX(false);
  style->SetPadGridY(false);
  style->SetGridColor(0);
  style->SetGridStyle(3);
  style->SetGridWidth(1);
  // For the frame:
  style->SetFrameBorderMode(0);
  style->SetFrameBorderSize(1);
  style->SetFrameFillColor(0);
  style->SetFrameFillStyle(0);
  style->SetFrameLineColor(1);
  style->SetFrameLineStyle(1);
  style->SetFrameLineWidth(1);
  // Margins:
  style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.15);//0.13);
  style->SetPadLeftMargin(0.15);//0.16);
  style->SetPadRightMargin(0.05);//0.02);
  // For the Global title:
  style->SetOptTitle(0);
  style->SetTitleFont(42);
  style->SetTitleColor(1);
  style->SetTitleTextColor(1);
  style->SetTitleFillColor(10);
  style->SetTitleFontSize(0.05);
  // For the axis titles:
  style->SetTitleColor(1, "XYZ");
  style->SetTitleFont(42, "XYZ");
  style->SetTitleSize(0.05, "XYZ");
  style->SetTitleXOffset(1.15);//0.9);
  style->SetTitleYOffset(1.5); // => 1.15 if exponents
  // For the axis labels:
  style->SetLabelColor(1, "XYZ");
  style->SetLabelFont(42, "XYZ");
  style->SetLabelOffset(0.007, "XYZ");
  style->SetLabelSize(0.045, "XYZ");
  // For the axis:
  style->SetAxisColor(1, "XYZ");
  style->SetStripDecimals(kTRUE);
  style->SetTickLength(0.03, "XYZ");
  style->SetNdivisions(510, "XYZ");
  //  style->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  style->SetPadTickY(1);
  // for histograms:
  style->SetHistLineColor(1);
  // for the pallete
  Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red  [5] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[5] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue [5] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 100);
  style->SetNumberContours(100);

  style->cd();
  return style;

}


std::string MT2DrawTools::getLumiText( float lumi ) {

  std::string returnText;
  if( lumi>=1.0 )
    returnText = (std::string)Form(" %.1f fb^{-1}", lumi);
  else if( lumi>0.01 )
    returnText = (std::string)Form(" %.0f pb^{-1}", 1000.*lumi);
  else 
    returnText = (std::string)Form(" %.1f pb^{-1}", 1000.*lumi);

  return returnText;

}


TPaveText* MT2DrawTools::getLabelTop( float lumi ) {

  char text[300];
  sprintf( text, "%s (13 TeV)", getLumiText(lumi).c_str() );
  //  sprintf( text, "CMS Preliminary, %s at #sqrt{s} = 13 TeV", getLumiText(lumi).c_str() );
  std::string text_str(text);
  return getLabelTop(text_str);

}

TPaveText* MT2DrawTools::getLabelTop( float lumi2015, float lumi2016 ) {

  char text[300];
  sprintf( text, "%s + %s (13 TeV)", getLumiText(lumi2015).c_str(), getLumiText(lumi2016).c_str() );
  //  sprintf( text, "CMS Preliminary, %s at #sqrt{s} = 13 TeV", getLumiText(lumi).c_str() );
  std::string text_str(text);
  return getLabelTop(text_str);

}


TPaveText* MT2DrawTools::getLabelTopSimulation( float lumi ) {

  char text[300];
  sprintf( text, "CMS Simulation, %.1f fb^{-1} at #sqrt{s} = 13 TeV", lumi );
  std::string text_str(text);
  return getLabelTopSimulation(text_str);

}

TPaveText* MT2DrawTools::getLabelTop( const std::string& text ) {

  TPaveText* label_top = new TPaveText(0.4,0.959,0.975,0.963, "brNDC");
  //  TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top->SetBorderSize(0);
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(42);  // label_top->SetTextFont(62);
  label_top->AddText(text.c_str());

  return label_top;

}


TPaveText* MT2DrawTools::getLabelTopSimulation( const std::string& text ) {

  TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top->SetBorderSize(0);
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right                                                                                                                                        
  label_top->SetTextFont(62);
  label_top->AddText(text.c_str());

  return label_top;

}


TPaveText* MT2DrawTools::getLabelCMS( const std::string& text ) {

  TPaveText* label_cms = new TPaveText(0.143,0.96,0.27,0.965, "brNDC");
  label_cms->SetBorderSize(0);
  label_cms->SetFillColor(kWhite);
  label_cms->SetTextSize(0.042);
  label_cms->SetTextAlign(11); // align left
  label_cms->SetTextFont(61);
  label_cms->AddText( text.c_str() );

  return label_cms;

}



void MT2DrawTools::addLabels( TCanvas* c1, float lumi, const std::string& text  ) {

  c1->cd();
  TPaveText* labelTop = MT2DrawTools::getLabelTop( lumi );
  labelTop->Draw("same");
  TPaveText* labelCMS = MT2DrawTools::getLabelCMS( text.c_str() );
  labelCMS->Draw("same");

}



TGraphAsymmErrors* MT2DrawTools::getPoissonGraph( TH1D* histo, bool drawZeros, const std::string& xerrType, float nSigma ) {

  //  histo->SetBinErrorOption(TH1::kPoisson);
  const double alpha = 1 - 0.6827;

  unsigned int nBins = histo->GetNbinsX();
  int emptyBins=0;
  for( unsigned i=1; i < nBins; ++i ) {
    if( histo->GetBinContent(i)==0 ) emptyBins += 1;
  }
  if( (float)emptyBins/(float)nBins > 0.4 ) drawZeros=false;

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(0);

  for( int iBin=1; iBin<(histo->GetXaxis()->GetNbins()+1); ++iBin ) {

    int y; // these are data histograms, so y has to be integer
    double x, xerr, yerrplus, yerrminus;
    x = histo->GetBinCenter(iBin);
    if( xerrType=="0" )
      xerr = 0.;
    else if( xerrType=="binWidth" )
      xerr = histo->GetBinWidth(iBin)/2.;
    else if( xerrType=="sqrt12" )
      xerr = histo->GetBinWidth(iBin)/sqrt(12.);
    else {
      std::cout << "[MT2DrawTools::getPoissonGraph] Unkown xerrType '" << xerrType << "'. Setting to bin width." << std::endl;
      xerr = histo->GetBinWidth(iBin);
    }

    y = (int)histo->GetBinContent(iBin);

    if( y==0 && !drawZeros ) continue;     
    
    double ym =  (y==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,y,1.));
    double yp =  ROOT::Math::gamma_quantile_c(alpha/2,y+1,1) ;
    
    //    yerrminus = histo->GetBinErrorLow(iBin);
    //    yerrplus = histo->GetBinErrorUp(iBin);
    
    //double ym, yp;
    //RooHistError::instance().getPoissonInterval(y,ym,yp,nSigma);
    
    yerrplus = yp - y;
    yerrminus = y - ym;

    //    if(y==0)
    //      std::cout << yerrplus << "\t" << yerrminus << std::endl;

    int thisPoint = graph->GetN();
    graph->SetPoint( thisPoint, x, y );
    graph->SetPointError( thisPoint, xerr, xerr, yerrminus, yerrplus );

  }

  return graph;

}



TGraphAsymmErrors* MT2DrawTools::getRatioGraph( TH1D* histo_data, TH1D* histo_mc, const std::string& xerrType){

  if( !histo_data || !histo_mc ) return 0;

  TGraphAsymmErrors* graph  = new TGraphAsymmErrors();
  
  //  TGraphAsymmErrors* graph_data = MT2DrawTools::getPoissonGraph(histo_data, false);
  TGraphAsymmErrors* graph_data = MT2DrawTools::getPoissonGraph(histo_data, true);
  
  for( int i=0; i < graph_data->GetN(); ++i){
    
    Double_t x_tmp, data;
    graph_data->GetPoint( i, x_tmp, data );

    Double_t data_errUp = graph_data->GetErrorYhigh(i);
    Double_t data_errDn = graph_data->GetErrorYlow(i);
    
    int iBin = histo_mc->FindBin(x_tmp);
    float mc = histo_mc->GetBinContent(iBin);
    float mc_err = histo_mc->GetBinError(iBin);


    float ratio = data/mc;
    float ratio_errUp = sqrt( data_errUp*data_errUp/(mc*mc) + mc_err*mc_err*data*data/(mc*mc*mc*mc) );
    float ratio_errDn = sqrt( data_errDn*data_errDn/(mc*mc) + mc_err*mc_err*data*data/(mc*mc*mc*mc) );

    double xerr;
    
    if( xerrType=="0" )
      xerr = 0.;
    else if( xerrType=="binWidth" )
      xerr = histo_mc->GetBinWidth(iBin)/2.;
    else if( xerrType=="sqrt12" )
      xerr = histo_mc->GetBinWidth(iBin)/sqrt(12.);
    else {
      std::cout << "[MT2DrawTools::getPoissonGraph] Unkown xerrType '" << xerrType << "'. Setting to bin width." << std::endl;
      xerr = histo_mc->GetBinWidth(iBin);
    }

    graph->SetPoint(i, x_tmp, ratio );
    graph->SetPointEYhigh(i, ratio_errUp );
    graph->SetPointEYlow(i, ratio_errDn );
    graph->SetPointEXhigh(i, xerr );
    graph->SetPointEXlow(i, xerr );

  }

  graph->SetLineColor(1);
  graph->SetMarkerColor(1);
  graph->SetMarkerStyle(20);

  return graph;

}


// TGraphAsymmErrors* MT2DrawTools::getRatioGraph( TH1D* histo_data, TH1D* histo_mc ){

//   if( !histo_data || !histo_mc ) return 0;

//   TGraphAsymmErrors* graph  = new TGraphAsymmErrors();
  
//   TGraphAsymmErrors* graph_data = MT2DrawTools::getPoissonGraph(histo_data, false, "binWidth");
  
//   for( int i=0; i < graph_data->GetN(); ++i){
    
//     Double_t x_tmp, data;
//     graph_data->GetPoint( i, x_tmp, data );

//     Double_t data_errUp = graph_data->GetErrorYhigh(i);
//     Double_t data_errDn = graph_data->GetErrorYlow(i);
    
//     int iBin = histo_mc->FindBin(x_tmp);
//     float mc = histo_mc->GetBinContent(iBin);
//     float mc_err = histo_mc->GetBinError(iBin);


//     float ratio = data/mc;
//     float ratio_errUp = sqrt( data_errUp*data_errUp/(mc*mc) + mc_err*mc_err*data*data/(mc*mc*mc*mc) );
//     float ratio_errDn = sqrt( data_errDn*data_errDn/(mc*mc) + mc_err*mc_err*data*data/(mc*mc*mc*mc) );

//     graph->SetPoint(i, x_tmp, ratio );
//     graph->SetPointEYhigh(i, ratio_errUp );
//     graph->SetPointEYlow(i, ratio_errDn );

//     float  xerr = histo_data->GetBinWidth(iBin)/2.;
//     graph->SetPointEXhigh(i, xerr);    
//     graph->SetPointEXlow(i, xerr);    

//   }

//   graph->SetLineColor(1);
//   graph->SetMarkerColor(1);
//   graph->SetMarkerStyle(20);

//   return graph;

// }

TH1D* MT2DrawTools::getBandAtOne( TH1D* h ){

  TH1D* h_band = (TH1D*)h->Clone( Form("%s_band", h->GetName()) );
  h_band->SetMarkerSize(0);
  h_band->SetFillColor ( h->GetLineColor()-4 );
  h_band->SetFillStyle (3001);
  for ( int iBin=1; iBin <= h->GetNbinsX(); iBin++){
    h_band->SetBinContent(iBin,1);
    double error = h->GetBinContent(iBin) ? h->GetBinError(iBin)/h->GetBinContent(iBin) : 0.0;
    h_band->SetBinError(iBin, error);
  }
  
  return h_band;

}


TPad* MT2DrawTools::getCanvasMainPad( bool logY ){
  
  std::string padApp = "";
  if( logY )
    padApp = "_log";
  TPad* pad1 = new TPad(Form("pad1%s", padApp.c_str()), Form("pad1%s", padApp.c_str()), 0, 0.3-0.1, 1, 1);
  pad1->SetBottomMargin(0.15);
  if( logY )
    pad1->SetLogy();

  return pad1;

}

TPad* MT2DrawTools::getCanvasRatioPad( bool logY ){

  std::string padApp = "";
  if( logY )
    padApp = "_log";
  TPad* pad2 = new TPad(Form("pad2%s", padApp.c_str()), Form("pad2%s", padApp.c_str()), 0, 0, 1, 0.21);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.1);

  return pad2;

}


TH2D* MT2DrawTools::getRatioAxes( float xMin, float xMax, float yMin, float yMax ){

  TH2D* h2_axes_ratio = new TH2D("axes_ratio", "", 10, xMin, xMax, 10, yMin, yMax );
  h2_axes_ratio->SetStats(0);
  h2_axes_ratio->GetXaxis()->SetLabelSize(0.00);
  h2_axes_ratio->GetXaxis()->SetTickLength(0.09);
  h2_axes_ratio->GetYaxis()->SetNdivisions(5,5,0);
  h2_axes_ratio->GetYaxis()->SetTitleSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitleOffset(0.4);
  h2_axes_ratio->GetYaxis()->SetLabelSize(0.17);
  h2_axes_ratio->GetYaxis()->SetTitle("2016/2017");
  //  h2_axes_ratio->GetYaxis()->SetTitle("sel 0/i");
  //  h2_axes_ratio->GetYaxis()->SetTitle("Data / MC");

  return h2_axes_ratio;

}


double MT2DrawTools::getSFError(double integral_data, double error_data, double integral_mc, double error_mc){

  double error_datamc = integral_data/integral_mc*(sqrt( (error_data/integral_mc)*(error_data/integral_mc) + (integral_data*error_mc/(integral_data*integral_data))*(integral_data*error_mc/(integral_data*integral_data)) ));
  
  return error_datamc;

}

TPaveText* MT2DrawTools::getRatioText( double integral_data, double integral_mc, double error_datamc ){
 
  TPaveText* ratioText = new TPaveText( 0.133, -0.051, 0.4, 0.1 , "brNDC" );
  ratioText->SetTextSize(0.035);
  ratioText->SetTextFont(62);
  ratioText->SetTextColor(2);
  ratioText->SetFillColor(0);
  ratioText->SetTextAlign(11);
  ratioText->AddText( Form("2016/2017 = %.2f +/- %.2f", integral_data/integral_mc, error_datamc) );
  //  ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", integral_data/integral_mc, error_datamc) );

  return ratioText;

}
  
TLine* MT2DrawTools::getSFLine(double integral_data, double integral_mc, float xMin, float xMax){

  double scaleFactor = integral_data/integral_mc;
  TLine* lineSF = new TLine(xMin, scaleFactor, xMax, scaleFactor);
  lineSF->SetLineColor(kRed);

  return lineSF;

}

TGraphErrors* MT2DrawTools::getSFBand(double integral_data, double error_data, double integral_mc, double error_mc, float xMin, float xMax){
  
  double error_datamc = MT2DrawTools::getSFError(integral_data, error_data, integral_mc, error_mc);

  double x[2]={(double)xMin, (double)xMax};
  double xerr[2]={0., 0.};
  double yerr[2]={error_datamc, error_datamc};
  double y[2]={integral_data/integral_mc, integral_data/integral_mc};

  TGraphErrors* SFband = new TGraphErrors(2, x, y, xerr, yerr);
  SFband->SetLineColor(0);
  SFband->SetFillColor(kRed);
  SFband->SetFillStyle(3244);
  
  return SFband;
  
}


TF1* MT2DrawTools::getSFFit(TGraphAsymmErrors* g_ratio, float xMin, float xMax){

  if( !g_ratio ) return 0;

  TF1* f=new TF1("f", "[0]", xMin, xMax);
  f->SetLineColor(kRed);
  f->SetParameter(0, 1.0);
  g_ratio->Fit(f, "0");
  
  return f;

}

void MT2DrawTools::getSFFitParameters(TF1* f, double &sf, double &sfErr, double &chi2, int &ndof){

  chi2  = f->GetChisquare();
  ndof     = f->GetNDF();
  sf    = f->GetParameter(0);
  sfErr = f->GetParError(0);

}

TGraphErrors* MT2DrawTools::getSFFitBand(TF1* f, float xMin, float xMax){
  
  double chi2, sf, sfErr;
  int ndof;
  MT2DrawTools::getSFFitParameters(f, sf, sfErr, chi2, ndof);

  double x[2]    ={(double)xMin, (double)xMax};
  double y[2]    ={sf, sf};

  double xerr[2] ={0., 0.};
  double yerr[2] ={sfErr, sfErr};

  TGraphErrors* SFband = new TGraphErrors(2, x, y, xerr, yerr);
  SFband->SetLineColor(0);
  SFband->SetFillColor(kRed);
  SFband->SetFillStyle(3244);
  
  return SFband;
  
}

TPaveText* MT2DrawTools::getFitText( TF1* f ){
 

  TPaveText* ratioText = new TPaveText( 0.135, -0.051, 0.4, 0.1 , "brNDC" );
  //ratioText->SetTextSize(0.025);
  ratioText->SetTextSize(0.031);
  ratioText->SetTextFont(62);
  ratioText->SetTextColor(2);
  ratioText->SetFillColor(0);
  ratioText->SetTextAlign(11);
  
  double chi2, sf, sfErr;
  int ndof;
  MT2DrawTools::getSFFitParameters(f, sf, sfErr, chi2, ndof);
  //ratioText->AddText( Form("Data/MC = %.2f #pm %.2f (#chi^{2}/ndof = %.2f / %d)", sf, sfErr, chi2, ndof) );
  ratioText->AddText( Form("Data/MC = %.2f #pm %.2f", sf, sfErr) );

  return ratioText;

}


TGraphErrors* MT2DrawTools::getSystBand(float xMin, float xMax, double SystErr){
  
  double x[2]={(double)xMin, (double)xMax};
  double xerr[2]={0., 0.};
  double yerr[2]={SystErr, SystErr};
  double y[2]={1.0, 1.0};

  TGraphErrors* SystBand = new TGraphErrors(2, x, y, xerr, yerr);
  SystBand->SetLineColor(0);
  SystBand->SetFillColor(kGray+2);
  SystBand->SetFillStyle(3244);
  
  return SystBand;
  
}


TH1D* MT2DrawTools::getMCBandHisto( TH1D* histo_mc, double SystErr ){

  TH1D* histoBand = (TH1D*) histo_mc->Clone("histo_band");
  for( int b=1; b <= histoBand->GetNbinsX(); ++b ){

    float thisStatErr = histoBand->GetBinError(b);
    float thisStats = histoBand->GetBinContent(b);
    float thisSystErr = thisStats*SystErr;
    float thisErr = sqrt(thisStatErr*thisStatErr+thisSystErr*thisSystErr);
    histoBand->SetBinError(b, thisErr);

  }

  histoBand->SetLineColor(0);
  histoBand->SetFillColor(kGray+2);
  histoBand->SetFillStyle(3244);

  return histoBand;

}


void MT2DrawTools::addOverflowSingleHisto( TH1D* yield ) {

  yield->SetBinContent(yield->GetNbinsX(),
		       yield->GetBinContent(yield->GetNbinsX()  )+
		       yield->GetBinContent(yield->GetNbinsX()+1)  );
  yield->SetBinError(  yield->GetNbinsX(),
		       sqrt(yield->GetBinError(yield->GetNbinsX() )*
			    yield->GetBinError(yield->GetNbinsX() )+
			    yield->GetBinError(yield->GetNbinsX()+1)*
			    yield->GetBinError(yield->GetNbinsX()+1)  ));

  yield->SetBinContent(yield->GetNbinsX()+1, 0.);
  yield->SetBinError  (yield->GetNbinsX()+1, 0.);

}



void MT2DrawTools::addOverflowSingleHisto( TH3D* yield3d ) {
  
  for (int y=1; y<=yield3d->GetNbinsY()+1; ++y)
    for (int z=1; z<=yield3d->GetNbinsZ()+1; ++z){

      yield3d->SetBinContent(yield3d->GetNbinsX(), y, z,
			   yield3d->GetBinContent(yield3d->GetNbinsX(), y, z  )+
			   yield3d->GetBinContent(yield3d->GetNbinsX()+1, y, z)  );
      yield3d->SetBinError(  yield3d->GetNbinsX(), y, z,
			   sqrt(yield3d->GetBinError(yield3d->GetNbinsX(), y, z  )*
				yield3d->GetBinError(yield3d->GetNbinsX(), y, z  )+
				yield3d->GetBinError(yield3d->GetNbinsX()+1, y, z)*
				yield3d->GetBinError(yield3d->GetNbinsX()+1, y, z)  ));
      
      yield3d->SetBinContent(yield3d->GetNbinsX()+1, y, z, 0.);
      yield3d->SetBinError  (yield3d->GetNbinsX()+1, y, z, 0.);
    }
  
}




TH1D* MT2DrawTools::getBand( TF1* f, const std::string& name ) {

 const int ndim_resp_q = f->GetNpar();
 TMatrixD emat_resp_q(ndim_resp_q, ndim_resp_q);
 gMinuit->mnemat(&emat_resp_q[0][0], ndim_resp_q);

 return getBand(f, emat_resp_q, name);

}



// Create uncertainty band (histogram) for a given function and error matrix
// in the range of the function.
TH1D* MT2DrawTools::getBand(TF1 *f, TMatrixD const& m, std::string name, bool getRelativeBand, int npx) {

 Bool_t islog = true;
 //double xmin = f->GetXmin()*0.9;
 //double xmax = f->GetXmax()*1.1; //fixes problem in drawing with c option
 double xmin = f->GetXmin();
 double xmax = f->GetXmax()*1.1; //fixes problem in drawing with c option
 int npar = f->GetNpar();
 //TString formula = f->GetExpFormula();

 // Create binning (linear or log)
 Double_t xvec[npx];
 xvec[0] = xmin;
 double dx = (islog ? pow(xmax/xmin, 1./npx) : (xmax-xmin)/npx);
 for (int i = 0; i != npx; ++i) {
   xvec[i+1] = (islog ? xvec[i]*dx : xvec[i]+dx);
 }


 //
 // Compute partial derivatives numerically
 // can be used with any fit function
 //
 Double_t sigmaf[npx];
 TH1D* h1_band = new TH1D(name.c_str(), "", npx, xvec);

 for( int ipx=0; ipx<npx; ++ipx ) {

   sigmaf[ipx] = 0.;
   Double_t partDeriv[npar];

   //compute partial derivatives of f wrt its parameters:
   for( int ipar=0; ipar<npar; ++ipar ) {

     Float_t pi = f->GetParameter(ipar);
     Float_t dpi = sqrt(m[ipar][ipar])*0.01; //small compared to the par sigma
     f->SetParameter(ipar, pi+dpi);
     Float_t fplus = f->Eval(xvec[ipx]); 
     f->SetParameter(ipar, pi-dpi);
     Float_t fminus = f->Eval(xvec[ipx]); 
     f->SetParameter(ipar, pi); //put it back as it was

     partDeriv[ipar] = (fplus-fminus)/(2.*dpi);

   } //for params

   //compute sigma(f) at x:
   for( int ipar=0; ipar<npar; ++ipar ) {
     for( int jpar=0; jpar<npar; ++jpar ) {
       sigmaf[ipx] += partDeriv[ipar]*partDeriv[jpar]*m[ipar][jpar];
     }
   }
   sigmaf[ipx] = sqrt(sigmaf[ipx]); //absolute band

   h1_band->SetBinContent( ipx, f->Eval(xvec[ipx]) );
   if( getRelativeBand )
     h1_band->SetBinError( ipx, sigmaf[ipx]/f->Eval(xvec[ipx]) );
   else
     h1_band->SetBinError( ipx, sigmaf[ipx] );

 } //for points

 h1_band->SetMarkerStyle(20);
 h1_band->SetMarkerSize(0);
 h1_band->SetFillColor(18);
 h1_band->SetFillStyle(3001);


 //TGraph* h1_statError = new TGraph(npx, xvec, sigmaf);
//TH2D* h2_axesStat = new TH2D("axesStat", "", 10, 20., 1400., 10, 0., 10.);
//h2_axesStat->GetXaxis()->SetNoExponent();
//h2_axesStat->GetXaxis()->SetMoreLogLabels();
//TCanvas* cStat = new TCanvas("cStat", "cStat", 600, 600);
//cStat->cd();
//cStat->SetLogx();
//h2_axesStat->Draw();
//h1_band->Draw("psame");
//std::string canvasName = "stat/" + name + ".eps";
//cStat->SaveAs(canvasName.c_str());

//delete h2_axesStat;
//delete cStat;

 return h1_band;

} //getband




std::vector<TCanvas*> MT2DrawTools::drawRegionYields_fromTree( const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units, const std::string& kinCuts, const std::string& topoCuts ) {
//void MT2DrawTools::drawRegionYields_fromTree( MT2Analysis<MT2EstimateTree>* data, std::vector<MT2Analysis<MT2EstimateTree>* >  bgYields, const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units, const std::string& kinCuts ) {


  TString sel_tstr(selection);
  if( sel_tstr.Contains("weight") ) {
    std::cout << "[MT2DrawTools::drawRegionYields_fromTree] WARNING!! Selection contains 'weight'!! Are you sure you know what you're doing??" << std::endl;
  }

  std::vector<TCanvas*> returnVector;

  system( Form("mkdir -p %s", outdir_.c_str()) );


  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;




  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = mc_->at(0)->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    TCanvas *old_c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(Form("c1_%s", iMT2->getName().c_str()));
    if( old_c1 != 0 ) delete old_c1;
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_data = (data_) ? data_->get(thisRegion)->tree : 0;
    TH1D* h1_data = 0;
    TGraphAsymmErrors* gr_data = 0;
    if( tree_data ) {
      h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
      //tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );
        tree_data->Project( "h1_data", varName.c_str(), Form("%f*weight*(%s)", lumi_, selection.c_str()) );
	//tree_data->Project( "h1_data", varName.c_str(), Form("%f*(%s)", data_->getWeight(), selection.c_str()) );
      if( addOverflow_ )
        MT2DrawTools::addOverflowSingleHisto(h1_data);
      gr_data = MT2DrawTools::getPoissonGraph(h1_data, false, "binWidth");
      gr_data->SetMarkerStyle(20);
      gr_data->SetMarkerSize(1.2);

      h1_data->SetFillColor( 855 );
      h1_data->SetLineColor( kBlack );
    }

    // TTree* tree_data2 = (data2_) ? data2_->get(thisRegion)->tree : 0;
    // TH1D* h1_data2 = 0;
    // TGraphAsymmErrors* gr_data2 = 0;
    // if( tree_data2 ) {
    //   h1_data2 = new TH1D("h1_data2", "", nBins, xMin, xMax );
    //   tree_data2->Project( "h1_data2", varName.c_str(), selection.c_str() );
    //   //tree_data2->Project( "h1_data2", varName.c_str(), Form("%f*(%s)", data2_->getWeight(), selection.c_str()) );
    //   if( addOverflow_ )
    //     MT2DrawTools::addOverflowSingleHisto(h1_data2);
    //   gr_data2 = MT2DrawTools::getPoissonGraph(h1_data2, false, "binWidth");
    //   gr_data2->SetMarkerStyle(24);
    //   gr_data2->SetMarkerSize(1.2);
    // }
    

    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<4; ++i ) { 
      //org    for( unsigned i=0; i<mc_->size(); ++i ) { 
      TTree* tree_mc = (mc_->at(i)->get(thisRegion)->tree);
      std::string thisName = "h1_" + mc_->at(i)->getName() + "_" + thisRegion.getName();
      TObject* obj = gROOT->FindObject(thisName.c_str());
      if( obj ) delete obj;
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_mc->Sumw2();
      if( selection!="" )
        tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi_, selection.c_str()) );
      //tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s)", lumi_*mc_->at(i)->getWeight(), selection.c_str()) );
      else
        tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight", lumi_) );
      //tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f", lumi_*mc_->at(i)->getWeight()) );

      if( addOverflow_ )
        MT2DrawTools::addOverflowSingleHisto(h1_mc);

      histos_mc.push_back(h1_mc);
    }



   std::vector< TH1D* > histos_mc_sig;
   for( unsigned i=4; i<mc_->size(); ++i ) { 
      TTree* tree_mc = (mc_->at(i)->get(thisRegion)->tree);
      std::string thisName = "h1_sig_" + mc_->at(i)->getName() + "_" + thisRegion.getName();
      TObject* obj = gROOT->FindObject(thisName.c_str());
      if( obj ) delete obj;
      TH1D* h1_sig_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_sig_mc->Sumw2();
      if( selection!="" )
        tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi_, selection.c_str()) );
      else
        tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight", lumi_) );
     
      if( addOverflow_ )
        MT2DrawTools::addOverflowSingleHisto(h1_sig_mc);

      h1_sig_mc->SetLineColor(mc_->at(i)->getColor() ); 
      h1_sig_mc->SetLineStyle(2 ); 
      h1_sig_mc->SetLineWidth(2 ); 
      histos_mc_sig.push_back(h1_sig_mc);
   }


    TH1D* mc_sum;
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      if( i==0 ) {
        mc_sum = new TH1D( *histos_mc[i] );
        mc_sum->SetName("mc_sum");
      } else {
        mc_sum->Add( histos_mc[i] );
      }
    }

    float scaleFactor = mcSF_;
    if( data_ ) {
      float sf;
      if( addOverflow_ ){
	//adding the overflow bin to the SF calculation
	std::cout << "Integrals data: " << h1_data->Integral(1, nBins+1) << "\t mc: " << mc_sum->Integral(1, nBins+1) << std::endl;
	sf  = h1_data->Integral(1, nBins+1)/mc_sum->Integral(1, nBins+1);

      }else{ 
	//not adding the overflow bin
	std::cout << "Integrals data: " << h1_data->Integral(1, nBins) << "\t mc: " << mc_sum->Integral(1, nBins) << std::endl;
	sf  = h1_data->Integral(1, nBins)/mc_sum->Integral(1, nBins);
      }
      std::cout << "SF: " << sf << std::endl;
      if( shapeNorm_ ) scaleFactor *= sf;
    }


    // float scaleFactor2 = mcSF2_;
    // if( data2_ ) {
    //   float sf;
    //   if( addOverflow_ ){
    // 	//adding the overflow bin to the SF calculation
    // 	std::cout << "Integrals: " << h1_data2->Integral(1, nBins+1) << "\t" << mc_sum->Integral(1, nBins+1) << std::endl;
    // 	sf  = h1_data2->Integral(1, nBins+1)/mc_sum->Integral(1, nBins+1);

    //   }else{ 
    // 	//not adding the overflow bin
    // 	std::cout << "Integrals: " << h1_data2->Integral(1, nBins) << "\t" << mc_sum->Integral(1, nBins) << std::endl;
    // 	sf  = h1_data2->Integral(1, nBins)/mc_sum->Integral(1, nBins);
    //   }
    //   std::cout << "SF: " << sf << std::endl;
    //   if( shapeNorm_ ) scaleFactor2 *= sf;
    // }

    TObject* oldStack = gROOT->FindObject(Form("bgStack_%s", iMT2->getName().c_str()));
    //TObject* oldStack = gROOT->FindObject("bgStack");
    if( oldStack ) delete oldStack;
    TH1D* histo_mc;
    //THStack* bgStack = new THStack( "bgStack", Form("bgStack_%s", iMT2->getName().c_str()) );
   
    THStack* bgStack = new THStack( Form("bgStack_%s", iMT2->getName().c_str()),"" );
    // THStack* bgStack = new THStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = 4 - i - 1;
      //org      int index = mc_->size() - i - 1;
      histos_mc[index]->SetFillColor( mc_->at(index)->getColor() );
      histos_mc[index]->SetLineColor( kBlack );
      //if( shapeNorm_ && data_ )
      histos_mc[index]->Scale( scaleFactor );

      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack->Add(histos_mc[index]);
    }


    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr_ );
    
    TGraphAsymmErrors* g_ratio = 0;
    TGraphAsymmErrors* g_ratio2  = 0;
    if( data_ )  g_ratio = MT2DrawTools::getRatioGraph(h1_data, mcBand);
    // if( data2_ ) g_ratio2 = MT2DrawTools::getRatioGraph(h1_data2, mcBand);
    // if( data2_ ) g_ratio2->SetMarkerStyle(24);    

    // TH1D* h_ratio = (TH1D*) h1_data->Clone("h_ratio");
    // h_ratio->Divide( mcBand );

    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    //TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr_);

    TF1* fSF = (data_) ? MT2DrawTools::getSFFit(g_ratio, xMin, xMax) : 0;
    //    TF1* fSF2 = (data2_) ? MT2DrawTools::getSFFit(g_ratio2, xMin, xMax) : 0;
    TGraphErrors* SFFitBand = (fSF) ? MT2DrawTools::getSFFitBand(fSF, xMin, xMax) : 0;
    
    // double error_data;
    // double integral_data = h1_data->IntegralAndError(0, nBins+1, error_data);

    double error_mc;
    double integral_mc = mc_sum->IntegralAndError(0, nBins+1, error_mc);

    //    double error_datamc = MT2DrawTools::getSFError(integral_data, error_data, integral_mc, error_mc );



    TCanvas* c1 = new TCanvas(Form("c1_%s", iMT2->getName().c_str()), "", 600, 600);
    c1->cd();
    
    TCanvas* c1_log = new TCanvas(Form("c1_log_%s", iMT2->getName().c_str()), "", 600, 600);

    float yMaxScale = 1.4;
    float yMax1 = (data_) ? h1_data->GetMaximum()*yMaxScale : 0.;
    float yMax2 = (data_) ? yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum())) : 0.;
    float yMax3 = yMaxScale*(bgStack->GetMaximum());
    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    if( histo_mc->GetNbinsX()<2 ) yMax *=3.;
    yMax*=1.25;
    
    std::string xAxisTitle;
    if( units!="" ) 
      xAxisTitle = (std::string)(Form("%s [%s]", axisName.c_str(), units.c_str()) );
    else
      xAxisTitle = (std::string)(Form("%s", axisName.c_str()) );

    std::string yAxisTitle;
    if(binWidth>0.99){
      if( units!="" ) 
	yAxisTitle = (std::string)(Form("Events / (%.0f %s)", binWidth, units.c_str()));
      else
	yAxisTitle = (std::string)(Form("Events / (%.0f)", binWidth));
    }
    else if(binWidth>0.099){
      if( units!="" ) 
	yAxisTitle = (std::string)(Form("Events / (%.2f %s)", binWidth, units.c_str()));
      else
	yAxisTitle = (std::string)(Form("Events / (%.2f)", binWidth));
    }
    else{
      if( units!="" ) 
	yAxisTitle = (std::string)(Form("Events / (%.4f %s)", binWidth, units.c_str()));
      else
	yAxisTitle = (std::string)(Form("Events / (%.4f)", binWidth));
    }

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
  
    TPad* pad1 = 0;
    if( this->twoPads() ) {
      pad1 = MT2DrawTools::getCanvasMainPad();
      pad1->Draw();
      pad1->cd();
    }


    h2_axes->Draw();

    float yMin_log = (data_ && h1_data->GetMinimum()>2.) ? 1. : 0.1;

    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, yMin_log, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();

    TPad* pad1_log = 0;
    if( this->twoPads() ) {
      pad1_log = MT2DrawTools::getCanvasMainPad( true );
      pad1_log->Draw();
      pad1_log->cd();
    } else {
      c1_log->SetLogy();
    }

    h2_axes_log->Draw();
   


    std::vector<std::string> niceNames = thisRegion.getNiceNames();
 //   for( unsigned i=0; i<niceNames.size(); ++i ) {
 //
 //     float yMax = 0.9-(float)i*0.05;
 //     float yMin = yMax - 0.05;
 //     TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
 //     regionText->SetTextSize(0.035);
 //     regionText->SetTextFont(42);
 //     regionText->SetFillColor(0);
 //     regionText->SetTextAlign(11);
 //     regionText->AddText( niceNames[i].c_str() );
 //
 //     //c1->cd();
 //     //regionText->Draw("same");
 // 
 //     //c1_log->cd();
 //     //regionText->Draw("same");
 // 
 //   }


    
    for( unsigned i=0; i<niceNames.size(); ++i ) { 
      
      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.030);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);

    
      if( i==0 ) {

        if(kinCuts!="") {
          regionText->AddText( kinCuts.c_str() );
        } else {
          regionText->AddText( niceNames[i].c_str() );
        }

      } else if( i==1 ) {
    
        if(topoCuts!="") {
          regionText->AddText( topoCuts.c_str() );
        } else {
          regionText->AddText( niceNames[i].c_str() );
        }

      }
    
      if( this->twoPads() )
        pad1->cd();
      else
        c1->cd();
      regionText->Draw("same");
      
      if( this->twoPads() )
        pad1_log->cd();
      else
        c1_log->cd();
      regionText->Draw("same");
      
    }


    TPaveText* normText = new TPaveText( 0.47, 0.7, 0.62, 0.82, "brNDC" );
    normText->SetFillColor(0);
    normText->SetTextSize(0.035);
    if( scaleFactor!=1. ) {
      if( displaySF_ )
        normText->AddText( Form("#splitline{2017 scaled}{by %.2f}", scaleFactor) );
      //        normText->AddText( Form("#splitline{MC scaled}{by %.2f}", scaleFactor) );
      else
	normText->AddText( "#splitline{Shape}{Norm.}" );
      if( this->twoPads() ) 
        pad1->cd();
      else
        c1->cd();
      normText->Draw("same");
      if( this->twoPads() ) 
        pad1_log->cd();
      else
        c1_log->cd();
      normText->Draw("same");
    }



    int addLines = (data_) ? 2 : 0;

    TLegend* legend2 = new TLegend( 0.5, 0.94-(mc_->size()+addLines-5)*0.05, 0.74, 0.94 );
    for( unsigned i=0; i<histos_mc_sig.size(); ++i ) { 
      legend2->AddEntry( histos_mc_sig[i], mc_->at(i+4)->getFullName().c_str(), "L" );
    }

    //    addLines = (data2_&&data_) ? 4 : 2;
    TLegend* legend = new TLegend( 0.75, 0.94-(mc_->size()+addLines-5)*0.05, 0.95, 0.94 );
    //    TLegend* legend = new TLegend( 0.67, 0.9-(mc_->size()+addLines)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    if( data_ ) {
      //      legend->AddEntry( gr_data, "Data", "P" );
      legend->AddEntry( h1_data, "Higgs 2016", "P" );
    }
    // if( data2_ ) {
    //   legend->AddEntry( gr_data2, "Data 2017", "P" );
    //   //std::cout<< h1_data2->GetName() << std::endl;
    //   //std::cout << "integral: " << h1_data2->Integral(3, -1) << std::endl;
    //   //std::cout << "data2(30-60): "<< h1_data2->GetBinContent(3) + h1_data2->GetBinContent(4) << std::endl;
    //   //data2_int = h1_data2->GetBinContent(3) + h1_data2->GetBinContent(4);
    // }
    //float mc_int = 0.;
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], mc_->at(i)->getFullName().c_str(), "F" );
    }


    //float qcdFrac = (histos_mc[0]->GetBinContent(3) + histos_mc[0]->GetBinContent(4) )/ mc_int;
    //std::cout << "qcd frac: " << qcdFrac << std::endl;
    //std::cout << "qcd estimate: " << data_int*qcdFrac << std::endl;
    //std::cout << "qcd true (0-30): " <<  1.27*(histos_mc[0]->GetBinContent(1) + histos_mc[0]->GetBinContent(2)) << std::endl;
    //std::cout << "qcd true (30-60): " << 1.27*(histos_mc[0]->GetBinContent(3) + histos_mc[0]->GetBinContent(3)) << std::endl;
    //    if( data_ ) 
    legend->AddEntry( mcBand, "Uncert.", "F" );

   
    
    
    TPaveText* fitText = (fSF) ? MT2DrawTools::getFitText( fSF ) : 0;
    
    //    TPaveText* ratioText = MT2DrawTools::getRatioText( integral_data, integral_mc, error_datamc );
    //    TLine* lineSF = MT2DrawTools::getSFLine(integral_data, integral_mc, xMin, xMax);
    //    TGraphErrors* SFband = MT2DrawTools::getSFBand(integral_data, error_data, integral_mc, error_mc, xMin, xMax);
    
    float yMinR=0.0;
    float yMaxR=2.0;
    
    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );

    std::string CMStext = doPaperPlots_ ? "CMS" : "CMS Preliminary";

    c1->cd();
    if( this->twoPads() )
      pad1->cd();
    legend->Draw("same");
    legend2->Draw("same");
    bgStack->Draw("histo same");

    for( unsigned i=0; i<histos_mc_sig.size(); ++i ) { 
      histos_mc_sig[i]->Draw("L same");
    }

    mcBand->Draw("E2 same");

    if( data_ ) {
      mcBand->Draw("E2 same");
      h1_data->Draw("E2 same");
      //     gr_data->Draw("p same");
    }
    // if( data2_ ) {
    //   //      mcBand->Draw("E2 same");
    //   gr_data2->Draw("p same");
    // }
    if( !shapeNorm_ && fitText )
      fitText->Draw("same");
    //    ratioText->Draw("same");
 

    (data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)c1, lumi_, "CMS Simulation"); 
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, "CMS Simulation"); 

    gPad->RedrawAxis();

    c1_log->cd();
    if( this->twoPads() )
      pad1_log->cd();
    legend->Draw("same");
    legend2->Draw("same");
    bgStack->Draw("histo same");

    for( unsigned i=0; i<histos_mc_sig.size(); ++i ) { 
      histos_mc_sig[i]->Draw("L same");
    }

    mcBand->Draw("E2 same");

    if( data_ ) {
      mcBand->Draw("E2 same");
      gr_data->Draw("p same");
    }
    // if( data2_ ) {
    //   //      mcBand->Draw("E2 same");
    //   gr_data2->Draw("p same");
    // }
    if( !shapeNorm_ && fitText )
      fitText->Draw("same");
    //    ratioText->Draw("same");
    (data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)c1_log, lumi_, "CMS Simulation"); 
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, "CMS Simulation"); 

    gPad->RedrawAxis();
    
    c1->cd();

    if( twoPads() ) {
      TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
      pad2->Draw();
      pad2->cd();

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      if( !shapeNorm_ ){

        //systBand->Draw("3,same");
        lineCentral->Draw("same");

        if( data_ ) {
          SFFitBand->Draw("3,same");
          fSF->Draw("same");
        }

//        SFband->Draw("3,same");
//        lineSF->Draw("same");

      }

      //      h_ratio->Draw("P,same" );
      // if( g_ratio ) g_ratio->Draw("PE,same");    
      // if( g_ratio2 ) g_ratio2->Draw("PE,same");    
      gPad->RedrawAxis();


      c1_log->cd();
      TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
      pad2_log->Draw();
      pad2_log->cd();

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      if( !shapeNorm_ ){

        //systBand->Draw("3,same");
        lineCentral->Draw("same");

        if( data_ ) {
          SFFitBand->Draw("3,same");
          fSF->Draw("same");
        }
        
//        SFband->Draw("3,same");
//        lineSF->Draw("same");

      }
      //      h_ratio->Draw("P,same" );
      if( g_ratio ) g_ratio->Draw("PE,same");
      if( g_ratio2 ) g_ratio2->Draw("PE,same");
      gPad->RedrawAxis();

    } // if twoPads


    std::string regionSaveName = (MT2Regions.size()!=1) ? "_" + thisRegion.getName() : "";

    c1->SaveAs( Form("%s/%s%s.eps", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1->SaveAs( Form("%s/%s%s.png", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1->SaveAs( Form("%s/%s%s.pdf", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );

    c1_log->SaveAs( Form("%s/%s%s_log.eps", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1_log->SaveAs( Form("%s/%s%s_log.png", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1_log->SaveAs( Form("%s/%s%s_log.pdf", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );

    returnVector.push_back(c1);

    delete h2_axes;
    delete h2_axes_log;
    delete h2_axes_ratio;
    
    delete h1_data;
    //    delete h1_data2;
  
    //for( unsigned i=0; i<histos_mc.size(); ++i )
    //  delete histos_mc[i];

    delete c1_log;

  }// for MT2 regions

  return returnVector;

}




float MT2DrawTools::getDataMCSF( TCanvas* c1 ) {


  TList* list = getCorrectList( c1 );
  

  TGraphAsymmErrors* gr_data = (TGraphAsymmErrors*)list->FindObject("Graph");

  TH1D* mcBand = (TH1D*)list->FindObject("histo_band");
  float mcInt = mcBand->Integral();

  float dataInt = MT2DrawTools::graphIntegral( gr_data );

  return dataInt/mcInt;

}




float MT2DrawTools::graphIntegral( TGraphAsymmErrors* graph, float xMin, float xMax ) {

  float integral = 0.;

  for( int iPoint=0; iPoint<graph->GetN(); ++iPoint ) {

    Double_t x, y;
    graph->GetPoint( iPoint, x, y );

    if( x >= xMin && x <= xMax ) 
      integral += y;

  }

  return integral;

}

    




TList* MT2DrawTools::getCorrectList( TCanvas* c1 ) {

  TList* list = c1->GetListOfPrimitives();

  if( list->GetSize()==2 ) { // means that it's a two-pad canvas

    std::cout << "Get list for ratio pad version" << std::endl;
    TPad* pad1 = (TPad*)list->FindObject("pad1");
    list = pad1->GetListOfPrimitives();

  }

  return list;

}

  


bool MT2DrawTools::twoPads() const {

  return data_ && mc_;

}



































































//Data will take over the place of the ratio one, so signal or w/e you want to put in, now WITH weight!
//Fuck that for a bit actually, can't be bothered, just draw the fractions, no ratio for now
std::vector<TCanvas*> MT2DrawTools::drawRegionFractions_fromTree( const std::string& saveName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, std::string axisName, const std::string& units, const std::string& kinCuts, const std::string& topoCuts ) {


  TString sel_tstr(selection);
  if( sel_tstr.Contains("weight") ) {
    std::cout << "[MT2DrawTools::drawRegionYields_fromTree] WARNING!! Selection contains 'weight'!! Are you sure you know what you're doing??" << std::endl;
  }

  std::vector<TCanvas*> returnVector;

  system( Form("mkdir -p %s", outdir_.c_str()) );


  //  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;


  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = mc_->at(0)->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    TCanvas *old_c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(Form("c1_%s", iMT2->getName().c_str()));
    if( old_c1 != 0 ) delete old_c1;
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_data = (data_) ? data_->get(thisRegion)->tree : 0;
    TH1D* h1_data = 0;
    TGraphAsymmErrors* gr_data = 0;
    if( tree_data ) {
      h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
      //      tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );
      //tree_data->Project( "h1_data", varName.c_str(), Form("%f*(%s)", data_->getWeight(), selection.c_str()) );
      tree_data->Project( "h1_data", varName.c_str(), Form("%f*weight*(%s)", lumi_, selection.c_str()) );
      if( addOverflow_ )
        MT2DrawTools::addOverflowSingleHisto(h1_data);

      //      h1_data->Scale(1./ h1_data->Integral());
      h1_data->SetMarkerStyle(20);
      h1_data->SetMarkerSize(1.2);


      h1_data->SetFillColor( 855 );
      h1_data->SetLineColor( kBlack );

      gr_data = MT2DrawTools::getPoissonGraph(h1_data, false, "binWidth");
      gr_data->SetMarkerStyle(20);
      gr_data->SetMarkerSize(1.2);
    }
    

    TGraphAsymmErrors* gr_mc = 0;

    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<mc_->size(); ++i ) { 
      TTree* tree_mc = (mc_->at(i)->get(thisRegion)->tree);
      std::string thisName = "h1_" + mc_->at(i)->getName() + "_" + thisRegion.getName();
      TObject* obj = gROOT->FindObject(thisName.c_str());
      if( obj ) delete obj;
      TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
      h1_mc->Sumw2();
      if( selection!="" )
        tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi_, selection.c_str()) );
      else
        tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight", lumi_) );

      if( addOverflow_ )
        MT2DrawTools::addOverflowSingleHisto(h1_mc);

      //h1_mc->SetFillColor( mc_->at(index)->getColor() );
      h1_mc->SetLineColor(  mc_->at(i)->getColor() );

      h1_mc->SetLineWidth( 2 );
      if(i>3)
	h1_mc->SetLineStyle( 2 );

      histos_mc.push_back(h1_mc);

      gr_mc = MT2DrawTools::getPoissonGraph(h1_mc, false, "binWidth");
      //if( varName == "nBJets"){
      // for(int iBin=1; iBin< h1_mc->GetNbinsX(); iBin++){
      //   std::cout << "iBin " << i << " = " << h1_mc->GetBinContent(iBin) << std::endl;
      //}}
    }





    TObject* oldStack = gROOT->FindObject(Form("bgStack_%s", iMT2->getName().c_str()));
    //TObject* oldStack = gROOT->FindObject("bgStack");
    if( oldStack ) delete oldStack;
    TH1D* histo_mc;
    //THStack* bgStack = new THStack( "bgStack", Form("bgStack_%s", iMT2->getName().c_str()) );
   
    THStack* bgStack = new THStack( Form("bgStack_%s", iMT2->getName().c_str()),"" );
    // THStack* bgStack = new THStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 
      int index = mc_->size() - i - 1;
      //   histos_mc[index]->SetFillColor( mc_->at(index)->getColor() );
      //histos_mc[index]->SetLineColor( kBlack );
      //if( shapeNorm_ && data_ )
      //      histos_mc[index]->Scale( scaleFactor );

      if(i==0) histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else histo_mc->Add(histos_mc[index]);
      bgStack->Add(histos_mc[index]);
    }



    // TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

    float scaleFactor = mcSF_;
    if( data_ ) {
      float sf;
      if( addOverflow_ ){
    	//adding the overflow bin to the SF calculation
    	std::cout << "Integrals: " << h1_data->Integral(1, nBins+1) << "\t" << histo_mc->Integral(1, nBins+1) << std::endl;
    	sf  = h1_data->Integral(1, nBins+1)/histo_mc->Integral(1, nBins+1);

      }else{ 
    	//not adding the overflow bin
    	std::cout << "Integrals: " << h1_data->Integral(1, nBins) << "\t" << histo_mc->Integral(1, nBins) << std::endl;
    	sf  = h1_data->Integral(1, nBins)/histo_mc->Integral(1, nBins);
      }
      std::cout << "SF: " << sf << std::endl;
      if( shapeNorm_ ) scaleFactor *= sf;
    }

    histo_mc->Scale(1./histo_mc->Integral(1, nBins+1));
    if( data_)    h1_data->Scale(1./h1_data->Integral(1, nBins+1));

    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr_ );
    
    TH1D* g_ratio = 0;
    if( data_ ) g_ratio = (TH1D*) h1_data->Clone("h1_data");
    if( data_)  g_ratio->Divide( histo_mc );
    //    TGraphAsymmErrors* g_ratio = 0;
    // if( data_ ) g_ratio = MT2DrawTools::getRatioGraph(h1_data, mcBand);
    
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    //TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr_);
    
    // TF1* fSF = (data_) ? MT2DrawTools::getSFFit(g_ratio, xMin, xMax) : 0;
    // TGraphErrors* SFFitBand = (fSF) ? MT2DrawTools::getSFFitBand(fSF, xMin, xMax) : 0;

    TCanvas* c1 = new TCanvas(Form("c1_%s", iMT2->getName().c_str()), "", 600, 600);
    c1->cd();
    
    TCanvas* c1_log = new TCanvas(Form("c1_log_%s", iMT2->getName().c_str()), "", 600, 600);

    float yMaxScale = 1.1;
    //  float yMax1 = (data_) ? h1_data->GetMaximum()*yMaxScale : 0.;
    //float yMax2 = (data_) ? yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum())) : 0.;
    float yMax1 = yMaxScale*(histos_mc[0]->GetMaximum()/ histos_mc[0]->Integral());
    float yMax2 = yMaxScale*(histos_mc[0]->GetMaximum()/ histos_mc[0]->Integral());
    float yMax3 = yMaxScale*(histos_mc[0]->GetMaximum()/ histos_mc[0]->Integral());
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      if ( yMaxScale*(histos_mc[i]->GetMaximum()/ histos_mc[i]->Integral() ) > yMax3 )
	yMax3 = yMaxScale*(histos_mc[i]->GetMaximum()/ histos_mc[i]->Integral());
    }

    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    if( histos_mc[0]->GetNbinsX()<2 ) yMax *=3.;
    yMax*=1.25;
    
    std::string xAxisTitle;
    if( units!="" ) 
      xAxisTitle = (std::string)(Form("%s [%s]", axisName.c_str(), units.c_str()) );
    else
      xAxisTitle = (std::string)(Form("%s", axisName.c_str()) );

    // std::string yAxisTitle;
    // if(binWidth>0.99){
    //   if( units!="" ) 
    // 	yAxisTitle = (std::string)(Form("Events / (%.0f %s)", binWidth, units.c_str()));
    //   else
    // 	yAxisTitle = (std::string)(Form("Events / (%.0f)", binWidth));
    // }
    // else if(binWidth>0.099){
    //   if( units!="" ) 
    // 	yAxisTitle = (std::string)(Form("Events / (%.2f %s)", binWidth, units.c_str()));
    //   else
    // 	yAxisTitle = (std::string)(Form("Events / (%.2f)", binWidth));
    // }
    // else{
    //   if( units!="" ) 
    // 	yAxisTitle = (std::string)(Form("Events / (%.4f %s)", binWidth, units.c_str()));
    //   else
    // 	yAxisTitle = (std::string)(Form("Events / (%.4f)", binWidth));
    // }

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    //    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
  
    TPad* pad1 = 0;
    if( this->twoPads() ) {
      std::cout << "two pads" << std::endl;
      pad1 = MT2DrawTools::getCanvasMainPad();
      pad1->Draw();
      pad1->cd();
    }

    h2_axes->Draw();

    float yMin_log = (data_ && h1_data->GetMinimum()>2.) ? 1. : 0.0001;
    //   float yMin_log = (data_ && h1_data->GetMinimum()>2.) ? 1. : 0.001;

    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, yMin_log, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    //    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();

    TPad* pad1_log = 0;
    if( this->twoPads() ) {
      pad1_log = MT2DrawTools::getCanvasMainPad( true );
      pad1_log->Draw();
      pad1_log->cd();
    } else {
      c1_log->SetLogy();
    }

    h2_axes_log->Draw();
   


    std::vector<std::string> niceNames = thisRegion.getNiceNames();
 //   for( unsigned i=0; i<niceNames.size(); ++i ) {
 //
 //     float yMax = 0.9-(float)i*0.05;
 //     float yMin = yMax - 0.05;
 //     TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
 //     regionText->SetTextSize(0.035);
 //     regionText->SetTextFont(42);
 //     regionText->SetFillColor(0);
 //     regionText->SetTextAlign(11);
 //     regionText->AddText( niceNames[i].c_str() );
 //
 //     //c1->cd();
 //     //regionText->Draw("same");
 // 
 //     //c1_log->cd();
 //     //regionText->Draw("same");
 // 
 //   }
 

    
    for( unsigned i=0; i<niceNames.size(); ++i ) { 
      
      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.030);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);

    
      if( i==0 ) {

        if(kinCuts!="") {
          regionText->AddText( kinCuts.c_str() );
        } else {
          regionText->AddText( niceNames[i].c_str() );
        }

      } else if( i==1 ) {
    
        if(topoCuts!="") {
          regionText->AddText( topoCuts.c_str() );
        } else {
          regionText->AddText( niceNames[i].c_str() );
        }

      }
    
      if( this->twoPads() )
        pad1->cd();
      else
        c1->cd();
      regionText->Draw("same");
      
      if( this->twoPads() )
        pad1_log->cd();
      else
        c1_log->cd();
      regionText->Draw("same");
      
    }


    TPaveText* normText = new TPaveText( 0.47, 0.7, 0.62, 0.82, "brNDC" );
    normText->SetFillColor(0);
    normText->SetTextSize(0.035);
    if( scaleFactor!=1. ) {
      if( displaySF_ )
        normText->AddText( Form("#splitline{Norm diff}{by %.2f}", scaleFactor) );
      else
    	normText->AddText( "#splitline{Shape}{Norm.}" );
      if( this->twoPads() ) 
        pad1->cd();
      else
        c1->cd();
        normText->Draw("same");
      if( this->twoPads() ) 
        pad1_log->cd();
      else
        c1_log->cd();
      normText->Draw("same");
    }


    int addLines = (data_) ? 2 : 0;
    TLegend* legend = new TLegend( 0.67, 0.9-(mc_->size()+addLines)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    //float data_int = 0.;
    if( data_ ) {
      //      legend->AddEntry( gr_data, "Data 2017", "P" );
      legend->AddEntry( h1_data, "H #gamma#gamma 2016", "P" );
      //      legend->AddEntry( gr_data, "Data", "P" );
    }
    //float mc_int = 0.;
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], mc_->at(i)->getFullName().c_str(), "F" );
    }
    //    if( data_ ) 
      //      legend->AddEntry( mcBand, "MC Uncert.", "F" );

   
    
      //  TPaveText* fitText = (fSF) ? MT2DrawTools::getFitText( fSF ) : 0;
    
    float yMinR=0.0;
    float yMaxR=2.0;
    
    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );

    std::string CMStext = doPaperPlots_ ? "CMS" : "CMS Preliminary";

    c1->cd();
    if( this->twoPads() )
      pad1->cd();
    legend->Draw("same");
    for( unsigned i=0; i<histos_mc.size(); ++i ) 
      histos_mc[i]->DrawNormalized("histo same");

    mcBand->Draw("E2 same");

    //    gr_mc->Draw("L samm
    if( data_ ) {
      mcBand->Draw("E2 same");
      h1_data->DrawNormalized("p same");
      //      gr_data->Draw("p same");
    }
    // if( !shapeNorm_ && fitText )
      //      fitText->Draw("same");
    //    ratioText->Draw("same");


    //    (data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)c1, lumi_, "CMS Simulation"); 
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, "CMS Simulation"); 

    gPad->RedrawAxis();

    c1_log->cd();
    if( this->twoPads() )
      pad1_log->cd();
    legend->Draw("same");
    for( unsigned i=0; i<histos_mc.size(); ++i ) 
      histos_mc[i]->DrawNormalized("histo same");
    mcBand->Draw("E2 same");
    if( data_ ) {
      //      mcBand->Draw("E2 same");
      h1_data->DrawNormalized("p same");
      //      gr_data->Draw("p same");
    }
    //    if( !shapeNorm_ && fitText )
    //  fitText->Draw("same");
    //    ratioText->Draw("same");
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)c1_log, lumi_, "CMS Simulation"); 
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, "CMS Simulation"); 

    gPad->RedrawAxis();
    
    c1->cd();

    if( twoPads() ) {
      TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
      pad2->Draw();
      pad2->cd();

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      if( !shapeNorm_ ){

        //systBand->Draw("3,same");
        lineCentral->Draw("same");

        if( data_ ) {
	  //          SFFitBand->Draw("3,same");
	  //          fSF->Draw("same");
        }

//        SFband->Draw("3,same");
//        lineSF->Draw("same");

      }

      // h1_data->DrawNormalized("p same");
      if( g_ratio ) g_ratio->Draw("PE,same");    
      gPad->RedrawAxis();


      c1_log->cd();
      TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
      pad2_log->Draw();
      pad2_log->cd();

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      if( !shapeNorm_ ){

        //systBand->Draw("3,same");
        lineCentral->Draw("same");

        if( data_ ) {
	  //          SFFitBand->Draw("3,same");
          //fSF->Draw("same");
        }
        
//        SFband->Draw("3,same");
//        lineSF->Draw("same");

      }
      if( g_ratio ) g_ratio->Draw("PE,same");
      gPad->RedrawAxis();

    } // if twoPads


    std::string regionSaveName = (MT2Regions.size()!=1) ? "_" + thisRegion.getName() : "";

    c1->SaveAs( Form("%s/norm_%s%s.eps", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1->SaveAs( Form("%s/norm_%s%s.png", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1->SaveAs( Form("%s/norm_%s%s.pdf", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );

    c1_log->SaveAs( Form("%s/norm_%s%s_log.eps", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1_log->SaveAs( Form("%s/norm_%s%s_log.png", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1_log->SaveAs( Form("%s/norm_%s%s_log.pdf", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );

    returnVector.push_back(c1);

    if( g_ratio )         delete h2_axes;
    if( g_ratio )         delete h2_axes_log;
    if( g_ratio )     delete h2_axes_ratio;
    
    if(data_)    delete h1_data;
    delete c1_log;

    if( g_ratio)    delete g_ratio;

  }// for MT2 regions

  return returnVector;

}













































//Draw 2D histograms from a tree
std::vector<TCanvas*> MT2DrawTools::drawRegion2D_fromTree( const std::string& saveName, const std::string& varName, const std::string& varName2, const std::string& selection, int nBins, float xMin, float xMax, int nBins2, float yMin, float yMax, std::string axisName, std::string axisName2, const std::string& units, const std::string& units2, const std::string& kinCuts, const std::string& topoCuts ) {


  TString sel_tstr(selection);
  if( sel_tstr.Contains("weight") ) {
    std::cout << "[MT2DrawTools::drawRegion2D_fromTree] WARNING!! Selection contains 'weight'!! Are you sure you know what you're doing??" << std::endl;
  }

  std::vector<TCanvas*> returnVector;

  system( Form("mkdir -p %s", outdir_.c_str()) );

  if( axisName=="" ) axisName = varName;
  if( axisName2=="" ) axisName2 = varName2;


  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = mc_->at(0)->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    TCanvas *old_c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(Form("c1_%s", iMT2->getName().c_str()));
    if( old_c1 != 0 ) delete old_c1;
  
    MT2Region thisRegion( (*iMT2) );

    // TTree* tree_data = (data_) ? data_->get(thisRegion)->tree : 0;
    // TH1D* h1_data = 0;
    // TGraphAsymmErrors* gr_data = 0;
    // if( tree_data ) {
    //   h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
    //   tree_data->Project( "h1_data", varName.c_str(), selection.c_str() );
    //   //tree_data->Project( "h1_data", varName.c_str(), Form("%f*(%s)", data_->getWeight(), selection.c_str()) );
    //   if( addOverflow_ )
    //     MT2DrawTools::addOverflowSingleHisto(h1_data);
    //   gr_data = MT2DrawTools::getPoissonGraph(h1_data, false, "binWidth");
    //   gr_data->SetMarkerStyle(20);
    //   gr_data->SetMarkerSize(1.2);
    // }
    


    std::vector< TH2D* > histos_mc;
    for( unsigned i=0; i<mc_->size(); ++i ) { 
      TTree* tree_mc = (mc_->at(i)->get(thisRegion)->tree);
      std::string thisName = "h1_" + mc_->at(i)->getName() + "_" + thisRegion.getName();
      TObject* obj = gROOT->FindObject(thisName.c_str());
      if( obj ) delete obj;
      TH2D* h1_mc = new TH2D( thisName.c_str(), "", nBins, xMin, xMax, nBins2, yMin, yMax );
      h1_mc->Sumw2();
      if( selection!="" )
        tree_mc->Project( thisName.c_str(), Form("%s:%s", varName.c_str(), varName2.c_str()) , Form("%f*weight*(%s)", lumi_, selection.c_str()) );
      //tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*(%s)", lumi_*mc_->at(i)->getWeight(), selection.c_str()) );
      else
        tree_mc->Project( thisName.c_str(), Form("%s:%s", varName.c_str(), varName2.c_str()), Form("%f*weight", lumi_)  );
      //tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f", lumi_*mc_->at(i)->getWeight()) );

      //      if( addOverflow_ )
      // MT2DrawTools::addOverflowSingleHisto(h1_mc);

      h1_mc->SetMarkerStyle( 20 );//      h1_mc->SetMarkerStyle( 6 );  
      h1_mc->SetMarkerSize( 0.4);
      h1_mc->SetMarkerColor( mc_->at(i)->getColor() );
      h1_mc->SetFillColor(   mc_->at(i)->getColor() );

      // h1_mc->SetLineWidth( 2 );
      // if(i>3)
      // h1_mc->SetLineStyle( 2 );

      histos_mc.push_back(h1_mc);
    }

    
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    //TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr_);
    
    // TF1* fSF = (data_) ? MT2DrawTools::getSFFit(g_ratio, xMin, xMax) : 0;
    // TGraphErrors* SFFitBand = (fSF) ? MT2DrawTools::getSFFitBand(fSF, xMin, xMax) : 0;

    TCanvas* c1 = new TCanvas(Form("c1_%s", iMT2->getName().c_str()), "", 600, 600);
    c1->cd();
    
    TCanvas* c1_log = new TCanvas(Form("c1_log_%s", iMT2->getName().c_str()), "", 600, 600);

    // float yMaxScale = 1.1;
    // float yMax1 = (data_) ? h1_data->GetMaximum()*yMaxScale : 0.;
    // float yMax2 = (data_) ? yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum())) : 0.;
    // float yMax3 = yMaxScale*(histos_mc[0]->GetMaximum()/ histos_mc[0]->Integral());
    // for( unsigned i=0; i<histos_mc.size(); ++i ) {  
    //   if ( yMaxScale*(histos_mc[i]->GetMaximum()/ histos_mc[i]->Integral() ) > yMax3 )
    // 	yMax3 = yMaxScale*(histos_mc[i]->GetMaximum()/ histos_mc[i]->Integral());
    // }


    // float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    // if( yMax3 > yMax ) yMax = yMax3;
    // if( histos_mc[0]->GetNbinsX()<2 ) yMax *=3.;
    // yMax*=1.25;
    
    std::string xAxisTitle;
    std::string yAxisTitle;
    if( units!="" ) 
      xAxisTitle = (std::string)(Form("%s [%s]", axisName.c_str(), units.c_str()) );
    else
      xAxisTitle = (std::string)(Form("%s", axisName.c_str()) );

    if( units2!="" ) 
      yAxisTitle = (std::string)(Form("%s [%s]", axisName2.c_str(), units2.c_str()) );
    else
      yAxisTitle = (std::string)(Form("%s", axisName2.c_str()) );

    TH2D* h2_axes = new TH2D("axes", "", nBins, xMin, xMax, nBins2, yMin, yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
  
    TPad* pad1 = 0;
    if( this->twoPads() ) {
      pad1 = MT2DrawTools::getCanvasMainPad();
      pad1->Draw();
      pad1->cd();
    }


    h2_axes->Draw();

    // float yMin_log = (data_ && h1_data->GetMinimum()>2.) ? 1. : 0.001;

    TH2D* h2_axes_log = new TH2D("axes_log", "", nBins, xMin, xMax, nBins2, yMin, yMax );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    //    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();

    TPad* pad1_log = 0;
    if( this->twoPads() ) {
      pad1_log = MT2DrawTools::getCanvasMainPad( true );
      pad1_log->Draw();
      pad1_log->cd();
    } else {
      c1_log->SetLogz();//      c1_log->SetLogy();
    }

    //    h2_axes_log->Draw();
   


    std::vector<std::string> niceNames = thisRegion.getNiceNames();
    for( unsigned i=0; i<niceNames.size(); ++i ) { 
      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.030);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);
    
      if( i==0 ) {
        if(kinCuts!="") 
          regionText->AddText( kinCuts.c_str() );
	else 
          regionText->AddText( niceNames[i].c_str() );
      } else if( i==1 ) {
	if(topoCuts!="") 
          regionText->AddText( topoCuts.c_str() );
        else 
          regionText->AddText( niceNames[i].c_str() );
      }
    
      if( this->twoPads() )
        pad1->cd();
      else
        c1->cd();
      regionText->Draw("same");
      
      if( this->twoPads() )
        pad1_log->cd();
      else
        c1_log->cd();
      regionText->Draw("same");
      
    }

    //    int addLines = (data_) ? 2 : 0;
    TLegend* legend = new TLegend( 0.67, 0.9-(mc_->size())*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    //float data_int = 0.;
    // if( data_ ) {
    //   legend->AddEntry( gr_data, "Data", "P" );
    // }
    //float mc_int = 0.;
    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      legend->AddEntry( histos_mc[i], mc_->at(i)->getFullName().c_str(), "F" );
    }
    //    if( data_ ) 
      //      legend->AddEntry( mcBand, "MC Uncert.", "F" );

   
    
      //  TPaveText* fitText = (fSF) ? MT2DrawTools::getFitText( fSF ) : 0;
    
    float yMinR=0.0;
    float yMaxR=2.0;
    
    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );

    std::string CMStext = doPaperPlots_ ? "CMS" : "CMS Preliminary";

    c1->cd();
    if( this->twoPads() )
      pad1->cd();

    // //for( unsigned i=0; i<histos_mc.size(); ++i ) {
    // for( int i=histos_mc.size()-1; i>=0; --i ) {
    //   if(i==(histos_mc.size()-1) )
    //   	histos_mc[i]->Draw("scat");
    //   else
    // 	histos_mc[i]->Draw("scat same");
    // }

    for( unsigned i=0; i<histos_mc.size(); ++i ) {
      if(i==0)
      	histos_mc[i]->Draw("scat");
      else
	histos_mc[i]->Draw("scat same");
    }
    legend->Draw("same");
    //  if( data_ ) {
    //   //      mcBand->Draw("E2 same");
    //   gr_data->Draw("p same");
    // }
    // if( !shapeNorm_ && fitText )
    //      fitText->Draw("same");
    //    ratioText->Draw("same");


    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)c1, lumi_, "CMS Simulation"); 
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, "CMS Simulation"); 

    gPad->RedrawAxis();

    c1_log->cd();
    if( this->twoPads() )
      pad1_log->cd();

    for( unsigned i=0; i<histos_mc.size(); ++i ) 
      histos_mc[i]->Draw("scat same");
    legend->Draw("same");

    //if( data_ ) {
      //      mcBand->Draw("E2 same");
    //  gr_data->Draw("p same");
    // }
    //    if( !shapeNorm_ && fitText )
    //  fitText->Draw("same");
    //    ratioText->Draw("same");
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)c1_log, lumi_, "CMS Simulation"); 
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, "CMS Simulation"); 

    gPad->RedrawAxis();
    
    c1->cd();

    if( twoPads() ) {
      TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
      pad2->Draw();
      pad2->cd();

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      if( !shapeNorm_ ){
        //systBand->Draw("3,same");
        lineCentral->Draw("same");
      }
      //      if( g_ratio ) g_ratio->Draw("PE,same");    
      gPad->RedrawAxis();

      c1_log->cd();
      TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
      pad2_log->Draw();
      pad2_log->cd();

      h2_axes_ratio->Draw("");
      lineCentral->Draw("same");
      if( !shapeNorm_ ){
        //systBand->Draw("3,same");
        lineCentral->Draw("same");
      }
      //      if( g_ratio ) g_ratio->Draw("PE,same");
      gPad->RedrawAxis();
    } // if twoPads

    std::string regionSaveName = (MT2Regions.size()!=1) ? "_" + thisRegion.getName() : "";

    c1->SaveAs( Form("%s/2d_%s%s.eps", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1->SaveAs( Form("%s/2d_%s%s.png", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1->SaveAs( Form("%s/2d_%s%s.pdf", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );

    // c1_log->SaveAs( Form("%s/2d_%s%s_log.eps", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    // c1_log->SaveAs( Form("%s/2d_%s%s_log.png", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    // c1_log->SaveAs( Form("%s/2d_%s%s_log.pdf", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );

    returnVector.push_back(c1);

    delete h2_axes;
    delete h2_axes_log;
    delete h2_axes_ratio;
    
    // delete h1_data;
    delete c1_log;

  }// for MT2 regions

  return returnVector;

}
































//Data will take over the place of the ratio one, so signal or w/e you want to put in, now WITH weight!
std::vector<TCanvas*> MT2DrawTools::drawSelections_fromTree( const std::string& saveName, const std::string& varName, std::vector<std::string>& selection, std::vector<std::string>& sel_name, int nBins, float xMin, float xMax, std::string axisName, const std::string& units, const std::string& kinCuts, const std::string& topoCuts ) {


  TString sel_tstr(selection[0]);
  if( sel_tstr.Contains("weight") ) {
    std::cout << "[MT2DrawTools::drawRegionYields_fromTree] WARNING!! Selection contains 'weight'!! Are you sure you know what you're doing??" << std::endl;
  }

  std::vector<TCanvas*> returnVector;

  system( Form("mkdir -p %s", outdir_.c_str()) );


  //  float binWidth = (xMax-xMin)/nBins;
  if( axisName=="" ) axisName = varName;


  TH1::AddDirectory(kTRUE); // stupid ROOT memory allocation needs this

  std::set<MT2Region> MT2Regions = mc_->at(0)->getRegions();
  
  for( std::set<MT2Region>::iterator iMT2 = MT2Regions.begin(); iMT2!=MT2Regions.end(); ++iMT2 ) {

    TCanvas *old_c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(Form("c1_%s", iMT2->getName().c_str()));
    if( old_c1 != 0 ) delete old_c1;
  
    MT2Region thisRegion( (*iMT2) );

    TTree* tree_data = (data_) ? data_->get(thisRegion)->tree : 0;
    TH1D* h1_data = 0;
    TGraphAsymmErrors* gr_data = 0;
    if( tree_data ) {
      h1_data = new TH1D("h1_data", "", nBins, xMin, xMax );
      tree_data->Project( "h1_data", varName.c_str(), selection[0].c_str() );
      //tree_data->Project( "h1_data", varName.c_str(), Form("%f*(%s)", data_->getWeight(), selection[0].c_str()) );
      if( addOverflow_ )
        MT2DrawTools::addOverflowSingleHisto(h1_data);

      //      h1_data->Scale(1./ h1_data->Integral());
      h1_data->SetMarkerStyle(20);
      h1_data->SetMarkerSize(1.2);

      gr_data = MT2DrawTools::getPoissonGraph(h1_data, false, "binWidth");
      gr_data->SetMarkerStyle(20);
      gr_data->SetMarkerSize(1.2);
    }
    



    std::vector< TH1D* > histos_mc_ratios;

    std::vector< TH1D* > histos_mc;
    for( unsigned i=0; i<mc_->size(); ++i ) { 
      for( unsigned count_sel=0; count_sel<selection.size(); ++count_sel ) { 

	TTree* tree_mc = (mc_->at(i)->get(thisRegion)->tree);
	std::string thisName = "h1_" + mc_->at(i)->getName() + "_" + thisRegion.getName()+"_"+ Form("%d",count_sel);

	TObject* obj = gROOT->FindObject(thisName.c_str());
	if( obj ) delete obj;
	TH1D* h1_mc = new TH1D( thisName.c_str(), "", nBins, xMin, xMax );
	h1_mc->Sumw2();
	if( selection[0]!="" )
	  tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight*(%s)", lumi_, selection[count_sel].c_str()) );
	else
	  tree_mc->Project( thisName.c_str(), varName.c_str(), Form("%f*weight", lumi_) );

	if( addOverflow_ )
	  MT2DrawTools::addOverflowSingleHisto(h1_mc);

	//h1_mc->SetFillColor( mc_->at(index)->getColor() );
	
	h1_mc->SetLineColor(  mc_->at(i)->getColor()+count_sel );

	h1_mc->SetLineWidth( 2 );
	if(i>3)
	  h1_mc->SetLineStyle( 2 );

	h1_mc->SetLineStyle( count_sel+1 );

	histos_mc.push_back(h1_mc);
	//if( varName == "nBJets"){
	// for(int iBin=1; iBin< h1_mc->GetNbinsX(); iBin++){
	//   std::cout << "iBin " << i << " = " << h1_mc->GetBinContent(iBin) << std::endl;

	TH1D* h1_mc_ratio = (TH1D*) h1_mc ->Clone(Form("%s_ratio",thisName.c_str() ));
	if((i+count_sel)>0)
	  histos_mc_ratios.push_back(h1_mc_ratio);

	//}}
      }
    }

  
    for( unsigned i=0; i<histos_mc_ratios.size(); ++i ) { 

      float integral = histos_mc[0]->Integral();
      float integral_r = histos_mc_ratios[i]->Integral();
      histos_mc_ratios[i]->Divide( histos_mc[0] );
      //      histos_mc_ratios[i]->Scale( integral / integral_r );

    }


    TObject* oldStack = gROOT->FindObject(Form("bgStack_%s", iMT2->getName().c_str()));
    //TObject* oldStack = gROOT->FindObject("bgStack");
    if( oldStack ) delete oldStack;
    TH1D* histo_mc;
    //THStack* bgStack = new THStack( "bgStack", Form("bgStack_%s", iMT2->getName().c_str()) );
   
    THStack* bgStack = new THStack( Form("bgStack_%s", iMT2->getName().c_str()),"" );
    // THStack* bgStack = new THStack("bgStack", "");
    for( unsigned i=0; i<histos_mc.size(); ++i ) { 

      int index = ( (mc_->size()) * (selection.size()) ) - i -1 ;

      //   histos_mc[index]->SetFillColor( mc_->at(index)->getColor() );
      //histos_mc[index]->SetLineColor( kBlack );
      //if( shapeNorm_ && data_ )
      //      histos_mc[index]->Scale( scaleFactor );

      if(i==0) 
	histo_mc = (TH1D*) histos_mc[index]->Clone("histo_mc");
      else 
	histo_mc->Add(histos_mc[index]);

      bgStack->Add(histos_mc[index]);
    }

    // TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

    float scaleFactor = mcSF_;
    if( data_ ) {
      float sf;
      if( addOverflow_ ){
    	//adding the overflow bin to the SF calculation
    	std::cout << "Integrals: " << h1_data->Integral(1, nBins+1) << "\t" << histo_mc->Integral(1, nBins+1) << std::endl;
    	sf  = h1_data->Integral(1, nBins+1)/histo_mc->Integral(1, nBins+1);

      }else{ 
    	//not adding the overflow bin
    	std::cout << "Integrals: " << h1_data->Integral(1, nBins) << "\t" << histo_mc->Integral(1, nBins) << std::endl;
    	sf  = h1_data->Integral(1, nBins)/histo_mc->Integral(1, nBins);
      }
      std::cout << "SF: " << sf << std::endl;
      if( shapeNorm_ ) scaleFactor *= sf;
    }

    //    histo_mc->Scale(1./histo_mc->Integral(1, nBins+1));
    if(data_)    h1_data->Scale(1./h1_data->Integral(1, nBins+1));


    TH1D* mcBand = MT2DrawTools::getMCBandHisto( histo_mc, lumiErr_ );
    
    TH1D* g_ratio = 0;
    if( data_ ) g_ratio = (TH1D*) h1_data->Clone("h1_data");
    if (data_) g_ratio->Divide( histo_mc );
    //    TGraphAsymmErrors* g_ratio = 0;
    // if( data_ ) g_ratio = MT2DrawTools::getRatioGraph(h1_data, mcBand);
    
    TLine* lineCentral = new TLine(xMin, 1.0, xMax, 1.0);
    lineCentral->SetLineColor(1);
    //TGraphErrors* systBand = MT2DrawTools::getSystBand(xMin, xMax, lumiErr_);
    
    // TF1* fSF = (data_) ? MT2DrawTools::getSFFit(g_ratio, xMin, xMax) : 0;
    // TGraphErrors* SFFitBand = (fSF) ? MT2DrawTools::getSFFitBand(fSF, xMin, xMax) : 0;

    TCanvas* c1 = new TCanvas(Form("c1_%s", iMT2->getName().c_str()), "", 600, 600);
    c1->cd();
    
    TCanvas* c1_log = new TCanvas(Form("c1_log_%s", iMT2->getName().c_str()), "", 600, 600);

    float yMaxScale = 1.1;
    //  float yMax1 = (data_) ? h1_data->GetMaximum()*yMaxScale : 0.;
    //float yMax2 = (data_) ? yMaxScale*(h1_data->GetMaximum() + sqrt(h1_data->GetMaximum())) : 0.;
    float yMax1 = yMaxScale*(histos_mc[0]->GetMaximum());
    float yMax2 = yMaxScale*(histos_mc[0]->GetMaximum());
    float yMax3 = yMaxScale*(histos_mc[0]->GetMaximum());
			     
    //  float yMax1 = yMaxScale*(histos_mc[0]->GetMaximum()/ histos_mc[0]->Integral());
    // float yMax2 = yMaxScale*(histos_mc[0]->GetMaximum()/ histos_mc[0]->Integral());
    // float yMax3 = yMaxScale*(histos_mc[0]->GetMaximum()/ histos_mc[0]->Integral());
 
   for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      if ( yMaxScale*(histos_mc[i]->GetMaximum()/ histos_mc[i]->Integral() ) > yMax3 )
	yMax3 = yMaxScale*(histos_mc[i]->GetMaximum()/ histos_mc[i]->Integral());
    }

    float yMax = (yMax1>yMax2) ? yMax1 : yMax2;
    if( yMax3 > yMax ) yMax = yMax3;
    if( histos_mc[0]->GetNbinsX()<2 ) yMax *=3.;
    yMax*=1.25;
    
    std::string xAxisTitle;
    if( units!="" ) 
      xAxisTitle = (std::string)(Form("%s [%s]", axisName.c_str(), units.c_str()) );
    else
      xAxisTitle = (std::string)(Form("%s", axisName.c_str()) );

    // std::string yAxisTitle;
    // if(binWidth>0.99){
    //   if( units!="" ) 
    // 	yAxisTitle = (std::string)(Form("Events / (%.0f %s)", binWidth, units.c_str()));
    //   else
    // 	yAxisTitle = (std::string)(Form("Events / (%.0f)", binWidth));
    // }
    // else if(binWidth>0.099){
    //   if( units!="" ) 
    // 	yAxisTitle = (std::string)(Form("Events / (%.2f %s)", binWidth, units.c_str()));
    //   else
    // 	yAxisTitle = (std::string)(Form("Events / (%.2f)", binWidth));
    // }
    // else{
    //   if( units!="" ) 
    // 	yAxisTitle = (std::string)(Form("Events / (%.4f %s)", binWidth, units.c_str()));
    //   else
    // 	yAxisTitle = (std::string)(Form("Events / (%.4f)", binWidth));
    // }

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
    h2_axes->SetXTitle(xAxisTitle.c_str());
    //    h2_axes->SetYTitle(yAxisTitle.c_str());

    c1->cd();
  


    TPad* pad1 = 0;
    if( true ) {
      //    if( this->twoPads() ) {
      pad1 = MT2DrawTools::getCanvasMainPad();
      pad1->Draw();
      pad1->cd();
    }


    h2_axes->Draw();

    float yMin_log = (data_ && h1_data->GetMinimum()>2.) ? 1. : 0.0001;
    //   float yMin_log = (data_ && h1_data->GetMinimum()>2.) ? 1. : 0.001;

    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, yMin_log, yMax*2.0 );
    h2_axes_log->SetXTitle(xAxisTitle.c_str());
    //    h2_axes_log->SetYTitle(yAxisTitle.c_str());

    c1_log->cd();

    TPad* pad1_log = 0;
    //   if( this->twoPads() ) {
      pad1_log = MT2DrawTools::getCanvasMainPad( true );
      pad1_log->Draw();
      pad1_log->cd();
    // } else {
    //   c1_log->SetLogy();
    // }

    h2_axes_log->Draw();
   


    std::vector<std::string> niceNames = thisRegion.getNiceNames();
 //   for( unsigned i=0; i<niceNames.size(); ++i ) {
 //
 //     float yMax = 0.9-(float)i*0.05;
 //     float yMin = yMax - 0.05;
 //     TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
 //     regionText->SetTextSize(0.035);
 //     regionText->SetTextFont(42);
 //     regionText->SetFillColor(0);
 //     regionText->SetTextAlign(11);
 //     regionText->AddText( niceNames[i].c_str() );
 //
 //     //c1->cd();
 //     //regionText->Draw("same");
 // 
 //     //c1_log->cd();
 //     //regionText->Draw("same");
 // 
 //   }
 

    
    for( unsigned i=0; i<niceNames.size(); ++i ) { 
      
      float yMax = 0.9-(float)i*0.05;
      float yMin = yMax - 0.05;
      TPaveText* regionText = new TPaveText( 0.18, yMin, 0.55, yMax, "brNDC" );
      regionText->SetTextSize(0.030);
      regionText->SetTextFont(42);
      regionText->SetFillColor(0);
      regionText->SetTextAlign(11);

    
      if( i==0 ) {

        if(kinCuts!="") {
          regionText->AddText( kinCuts.c_str() );
        } else {
          regionText->AddText( niceNames[i].c_str() );
        }

      } else if( i==1 ) {
    
        if(topoCuts!="") {
          regionText->AddText( topoCuts.c_str() );
        } else {
          regionText->AddText( niceNames[i].c_str() );
        }

      }
    
      //      if( this->twoPads() )
        pad1->cd();
	//else
        //c1->cd();
      //regionText->Draw("same");
      
	//      if( this->twoPads() )
        pad1_log->cd();
	//      else
	//        c1_log->cd();
      //regionText->Draw("same");
      
    }


    TPaveText* normText = new TPaveText( 0.47, 0.7, 0.62, 0.82, "brNDC" );
    normText->SetFillColor(0);
    normText->SetTextSize(0.035);
    if( scaleFactor!=1. ) {
      if( displaySF_ )
        normText->AddText( Form("#splitline{Norm diff}{by %.2f}", scaleFactor) );
      else
    	normText->AddText( "#splitline{Shape}{Norm.}" );
      //      if( this->twoPads() ) 
        pad1->cd();
	//else
	//        c1->cd();

        normText->Draw("same");
	//      if( this->twoPads() ) 
        pad1_log->cd();
	//      else
	//        c1_log->cd();

      normText->Draw("same");
    }


    int addLines = (data_) ? 2 : 0;
    TLegend* legend = new TLegend( 0.67, 0.9-(mc_->size()*selection.size() +addLines)*0.06, 0.93, 0.9 );
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    //float data_int = 0.;
    if( data_ ) {
      legend->AddEntry( gr_data, "Data 2017", "P" );
      //      legend->AddEntry( gr_data, "Data", "P" );
    }
    //float mc_int = 0.;


 

    for( unsigned i=0; i<histos_mc.size(); ++i ) {  
      if(i==0)
	legend->SetHeader( Form("%s", mc_->at(i%mc_->size())->getFullName().c_str()), "C");
      legend->AddEntry( histos_mc[i], Form("%s", sel_name[i].c_str()), "F" );
      //    legend->AddEntry( histos_mc[i], mc_->at(i%mc_->size())->getFullName().c_str(), "F" );
    }

    //    if( data_ ) 
    //      legend->AddEntry( mcBand, "MC Uncert.", "F" );
    
    //  TPaveText* fitText = (fSF) ? MT2DrawTools::getFitText( fSF ) : 0;
    
    float yMinR=0.50;
    float yMaxR=1.50;

    // float yMinR=0.0;
    // float yMaxR=2.0;
    
    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );
    //    TH2D* h2_axes_ratio = MT2DrawTools::getRatioAxes( xMin, xMax, yMinR, yMaxR );

    std::string CMStext = doPaperPlots_ ? "CMS" : "CMS Preliminary";

    c1->cd();
    //    if( this->twoPads() )
      pad1->cd();
    legend->Draw("same");
    for( unsigned i=0; i<histos_mc.size(); ++i ) 
      histos_mc[i]->Draw("histo same");
    //      histos_mc[i]->DrawNormalized("histo same");
    if( data_ ) {
      //      mcBand->Draw("E2 same");
      h1_data->DrawNormalized("p same");
      //      gr_data->Draw("p same");
    }
    // if( !shapeNorm_ && fitText )
      //      fitText->Draw("same");
    //    ratioText->Draw("same");


    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)c1, lumi_, "CMS Simulation"); 
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)pad1, lumi_, "CMS Simulation"); 

    gPad->RedrawAxis();

    c1_log->cd();
    //    if( this->twoPads() )
    pad1_log->cd();
    legend->Draw("same");
    for( unsigned i=0; i<histos_mc.size(); ++i ) 
      histos_mc[i]->Draw("histo same");
    //      histos_mc[i]->DrawNormalized("histo same");
    if( data_ ) {
      //      mcBand->Draw("E2 same");
      h1_data->DrawNormalized("p same");
      //      gr_data->Draw("p same");
    }
    //    if( !shapeNorm_ && fitText )
    //  fitText->Draw("same");
    //    ratioText->Draw("same");
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)c1_log, lumi_, "CMS Simulation"); 
    //(data_) ? MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, CMStext.c_str() ) : MT2DrawTools::addLabels( (TCanvas*)pad1_log, lumi_, "CMS Simulation"); 

    gPad->RedrawAxis();
    
    c1->cd();

    //    if( twoPads() ) {
    TPad* pad2 = MT2DrawTools::getCanvasRatioPad();
    pad2->Draw();
    pad2->cd();

    h2_axes_ratio->Draw("");
    lineCentral->Draw("same");
    if( !shapeNorm_ ){

      //systBand->Draw("3,same");
      lineCentral->Draw("same");

      if( data_ ) {
	//          SFFitBand->Draw("3,same");
	//          fSF->Draw("same");
      }

      //        SFband->Draw("3,same");
      //        lineSF->Draw("same");

    }

    // h1_data->DrawNormalized("p same");
    if( g_ratio ) g_ratio->Draw("PE,same");  

    std::cout << "Added a ratio histo" << histos_mc_ratios.size() <<  std::endl;

    for( unsigned i=0; i<histos_mc_ratios.size(); ++i ) { 
      histos_mc_ratios[i]->Draw("L same");
      std::cout << "Added a ratio histo" << std::endl;
    }
  
    gPad->RedrawAxis();


    c1_log->cd();
    TPad* pad2_log = MT2DrawTools::getCanvasRatioPad( true );
    pad2_log->Draw();
    pad2_log->cd();

    h2_axes_ratio->Draw("");
    lineCentral->Draw("same");
      if( !shapeNorm_ ){

        //systBand->Draw("3,same");
        lineCentral->Draw("same");

        if( data_ ) {
	  //          SFFitBand->Draw("3,same");
          //fSF->Draw("same");
        }
        
//        SFband->Draw("3,same");
//        lineSF->Draw("same");

      }
      if( g_ratio ) g_ratio->Draw("PE,same");
      for( unsigned i=0; i<histos_mc_ratios.size(); ++i ) { 
	histos_mc_ratios[i]->Draw("L same");
      }

      gPad->RedrawAxis();

      //    } // if twoPads


    std::string regionSaveName = (MT2Regions.size()!=1) ? "_" + thisRegion.getName() : "";

    c1->SaveAs( Form("%s/norm_%s%s.eps", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1->SaveAs( Form("%s/norm_%s%s.png", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1->SaveAs( Form("%s/norm_%s%s.pdf", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );

    c1_log->SaveAs( Form("%s/norm_%s%s_log.eps", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1_log->SaveAs( Form("%s/norm_%s%s_log.png", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );
    c1_log->SaveAs( Form("%s/norm_%s%s_log.pdf", outdir_.c_str(), saveName.c_str(), regionSaveName.c_str()) );

    returnVector.push_back(c1);

    delete h2_axes;
    delete h2_axes_log;
    delete h2_axes_ratio;
    
    delete h1_data;
    delete c1_log;

    delete g_ratio;

  }// for MT2 regions

  return returnVector;

}
