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



TPaveText* MT2DrawTools::getLabelTop( float lumi ) {

  char text[300];
  sprintf( text, "CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13 TeV", lumi );
  std::string text_str(text);
  return getLabelTop(text_str);

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




TGraphAsymmErrors* MT2DrawTools::getPoissonGraph( TH1D* histo, bool drawZeros, const std::string& xerrType, float nSigma ) {

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
