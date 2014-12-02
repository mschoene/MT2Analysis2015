#ifndef MT2DrawTools_h
#define MT2DrawTools_h

#include "TStyle.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"






class MT2DrawTools {

 public:

  static TStyle* setStyle();

  static TPaveText* getLabelTop( const std::string& text="CMS Preliminary, #sqrt{s} = 13 TeV" );

  static TGraphAsymmErrors* getPoissonGraph( TH1D* h1, bool drawZeros=true, const std::string& xerrType="0", float nSigma=1. );

 private:

};

#endif