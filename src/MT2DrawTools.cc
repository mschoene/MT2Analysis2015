#include "../interface/MT2DrawTools.h"

#include "RooHistError.h"




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
  style->SetTitleYOffset(1.3); // => 1.15 if exponents
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
  style->SetPadTickX(1); // To get tick marks on the opposite side of the frame
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
  sprintf( text, "CMS Preliminary, %s at #sqrt{s} = 13 TeV", getLumiText(lumi).c_str() );
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

  TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top->SetBorderSize(0);
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(62);
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


TGraphAsymmErrors* MT2DrawTools::getPoissonGraph( TH1D* histo, bool drawZeros, const std::string& xerrType, float nSigma ) {


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

    double ym, yp;
    RooHistError::instance().getPoissonInterval(y,ym,yp,nSigma);
    yerrplus = yp - y;
    yerrminus = y - ym;

    int thisPoint = graph->GetN();
    graph->SetPoint( thisPoint, x, y );
    graph->SetPointError( thisPoint, xerr, xerr, yerrminus, yerrplus );

  }

  return graph;

}


TGraphAsymmErrors* MT2DrawTools::getRatioGraph( TH1D* histo_data, TH1D* histo_mc ){
  
  TGraphAsymmErrors* graph_ = MT2DrawTools::getPoissonGraph(histo_data, false);
  TGraphAsymmErrors* graph  = new TGraphAsymmErrors();
  
  for( int i=0; i < graph_->GetN(); ++i){
    
    Double_t x_tmp, y_tmp, errUp, errDown;       
    graph_->GetPoint( i, x_tmp, y_tmp );

    errUp   = graph_->GetErrorYhigh(i);
    errDown = graph_->GetErrorYlow(i);
    
    int iBin = histo_mc->FindBin(x_tmp);
    float mc = histo_mc->GetBinContent(iBin);
    graph->SetPoint(i, x_tmp, y_tmp/mc);
    graph->SetPointEYhigh(i, errUp/mc);
    graph->SetPointEYlow(i, errDown/mc);
    

  }

  graph->SetLineColor(1);
  graph->SetMarkerColor(1);
  graph->SetMarkerStyle(20);

  return graph;

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
  h2_axes_ratio->GetYaxis()->SetTitle("Data / MC");

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
  ratioText->AddText( Form("Data/MC = %.2f +/- %.2f", integral_data/integral_mc, error_datamc) );

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
 
  double chi2, sf, sfErr;
  int ndof;
  MT2DrawTools::getSFFitParameters(f, sf, sfErr, chi2, ndof);

  TPaveText* ratioText = new TPaveText( 0.135, -0.051, 0.4, 0.1 , "brNDC" );
  ratioText->SetTextSize(0.025);
  ratioText->SetTextFont(62);
  ratioText->SetTextColor(2);
  ratioText->SetFillColor(0);
  ratioText->SetTextAlign(11);
  
  ratioText->AddText( Form("Data/MC = %.2f #pm %.2f (#chi^{2}/ndof = %.2f / %d)", sf, sfErr, chi2, ndof) );

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
