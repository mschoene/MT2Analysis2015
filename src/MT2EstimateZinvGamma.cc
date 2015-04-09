#include "../interface/MT2EstimateZinvGamma.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TRandom3.h"
#include "TTree.h"





MT2EstimateZinvGamma::MT2EstimateZinvGamma( const std::string& aname, const MT2Region& aregion ) : MT2Estimate( aname, aregion ) {


  sietaieta = new TH1D( this->getHistoName("sietaieta").c_str(), "", 200, 0., 0.035 );
  sietaieta->Sumw2();


  int nbins = 8;
  int xmax = 20.;

  // this histo will be used to create histogram templates:
  //iso = new TH1D( this->getHistoName("iso").c_str(), "", nbins-1, bins );
  iso = new TH1D( this->getHistoName("iso").c_str(), "", nbins, 0., xmax );
  iso->Sumw2();

  x_ = new RooRealVar( "x", "", 0., xmax );
  w_ = new RooRealVar( "w", "", 0., 1000. );

  int nbins_mt2;
  double *bins_mt2;
  aregion.getBins( nbins_mt2, bins_mt2 );


  for( unsigned i=0; i<nbins_mt2; ++i ) {

    RooDataSet* isoDataset = new RooDataSet( this->getHistoName(Form("iso_bin%d", i)).c_str(), "", RooArgSet(*x_,*w_), w_->GetName() );
    iso_bins.push_back(isoDataset);

    //TH1D* this_iso_hist = new TH1D( this->getHistoName(Form("iso_bin%d_hist", i)).c_str() , "", nbins-1, bins );
    TH1D* this_iso_hist = new TH1D( this->getHistoName(Form("iso_bin%d_hist", i)).c_str() , "", nbins, 0., xmax );
    this_iso_hist->Sumw2();
    iso_bins_hist.push_back(this_iso_hist);

  }

}




MT2EstimateZinvGamma::MT2EstimateZinvGamma( const MT2EstimateZinvGamma& rhs ) : MT2Estimate(rhs) {

  this->iso = new TH1D(*(rhs.iso));
  this->sietaieta = new TH1D(*(rhs.sietaieta));

  int xmax = this->iso->GetXaxis()->GetXmax();

  this->x_ = new RooRealVar( "x", "", 0., xmax );
  this->w_ = new RooRealVar( "w", "", 0., 1000. );

  for( unsigned i=0; i<rhs.iso_bins.size(); ++i ) {
    RooDataSet* newDataSet = new RooDataSet( *(rhs.iso_bins[i]) );
    this->iso_bins.push_back(newDataSet);
    TH1D* this_iso_hist = new TH1D( *(rhs.iso_bins_hist[i]) );
    iso_bins_hist.push_back(this_iso_hist);
  }

}



// destructor

MT2EstimateZinvGamma::~MT2EstimateZinvGamma() {

  delete iso;
  delete sietaieta;

  for( unsigned i=0; i<iso_bins.size(); ++i ) { 
    delete iso_bins[i];
    delete iso_bins_hist[i];
  }


}



void MT2EstimateZinvGamma::setName( const std::string& newName ) {

  MT2Estimate::setName(newName);

  iso->SetName( this->getHistoName("iso").c_str() );
  sietaieta->SetName( this->getHistoName("sietaieta").c_str() );

  for( unsigned i=0; i<iso_bins.size(); ++i ) {
    iso_bins[i]->SetName( this->getHistoName(Form("iso_bin%d", i)).c_str() );
    iso_bins_hist[i]->SetName( this->getHistoName(Form("iso_bin%d_hist", i)).c_str() );
  }

}



void MT2EstimateZinvGamma::fillIso( float iso, float weight, float mt2 ) {

  this->iso->Fill( iso, weight );


  if( mt2>0. ) {
    int foundBin = this->yield->FindBin(mt2);
    if( foundBin > this->yield->GetNbinsX() ) foundBin=this->yield->GetNbinsX(); // overflow will go in last bin
    foundBin-=1; // want first bin to be 0 (fuck you root)
    if( foundBin>=0 ) {
      x_->setVal(iso);
      w_->setVal(weight);
      iso_bins[foundBin]->add( RooArgList(*x_, *w_), weight );
      iso_bins_hist[foundBin]->Fill( iso, weight );
    }
  }
      
}
 


RooDataSet* MT2EstimateZinvGamma::isoData() const {

  RooDataSet* dataset = new RooDataSet( this->getHistoName("isoData").c_str(), "", RooArgSet(*x_,*w_), w_->GetName() );
  
  for( unsigned i=0; i<iso_bins.size(); ++i ) 
    dataset->append( *(iso_bins[i]) );
  
  return dataset;

}
 

void MT2EstimateZinvGamma::fakeDatasetsFromHistos(int seed) {

  TRandom3 rand(seed);

  for( unsigned i=0; i<iso_bins_hist.size(); ++i ) {

    RooDataSet* newdataset = new RooDataSet( this->getHistoName(Form("iso_bin%d", i)).c_str(), "", RooArgSet(*x_,*w_), w_->GetName() );

    for( unsigned ibin=1; ibin<iso_bins_hist[i]->GetNbinsX()+1; ++ibin ) {

      float isoValue = iso_bins_hist[i]->GetBinLowEdge(ibin);
      //float isoValue = iso_bins_hist[i]->GetBinCenter(ibin);

      //int entries = (int)iso_bins_hist[i]->GetBinContent(ibin);
      int entries = rand.Poisson( iso_bins_hist[i]->GetBinContent(ibin) );

      for( unsigned ientry=0; ientry<entries; ++ientry ) {

        x_->setVal(isoValue);
        w_->setVal(1.);

        newdataset->add( RooArgList(*x_, *w_), 1. );

      } // for entries

    } // for bins

    iso_bins[i] = new RooDataSet( *newdataset );

  } // for histos

}
    



void MT2EstimateZinvGamma::finalize() {

  MT2Estimate::finalize();

  MT2Estimate::addOverflowSingleHisto( iso );
  MT2Estimate::addOverflowSingleHisto( sietaieta );

  for( unsigned i=0; i<iso_bins.size(); ++i ) {
    MT2Estimate::addOverflowSingleHisto(iso_bins_hist[i]);
  }

}



void MT2EstimateZinvGamma::getShit( TFile* file, const std::string& path ) {

  MT2Estimate::getShit(file, path);
  iso = (TH1D*)file->Get(Form("%s/%s", path.c_str(), iso->GetName()));
  sietaieta = (TH1D*)file->Get(Form("%s/%s", path.c_str(), sietaieta->GetName()));

  for( unsigned i=0; i<iso_bins.size(); ++i ) {
    iso_bins[i] = (RooDataSet*)file->Get(Form("%s/%s", path.c_str(), iso_bins[i]->GetName()));
    iso_bins_hist[i] = (TH1D*)file->Get(Form("%s/%s", path.c_str(), iso_bins_hist[i]->GetName()));
  }


}


void MT2EstimateZinvGamma::print(const std::string& ofs){

  MT2Estimate::print( ofs );

}



void MT2EstimateZinvGamma::write() const {

  MT2Estimate::write();

  iso->Write();
  sietaieta->Write();

  for( unsigned i=0; i<iso_bins.size(); ++i ) {
    iso_bins[i]->Write();
    iso_bins_hist[i]->Write();
  }

}



const MT2EstimateZinvGamma& MT2EstimateZinvGamma::operator=( const MT2EstimateZinvGamma& rhs ) {


  this->region = new MT2Region(*(rhs.region));

  this->yield = new TH1D(*(rhs.yield));

  this->iso = new TH1D(*(rhs.iso));

  this->sietaieta = new TH1D(*(rhs.sietaieta));

  for( unsigned i=0; i<iso_bins.size(); ++i ) {
    this->iso_bins[i] = new RooDataSet( *(rhs.iso_bins[i]) );
    this->iso_bins_hist[i] = new TH1D( *(rhs.iso_bins_hist[i]) );
  }

  this->setName(this->getName());

  return *this;

}




MT2EstimateZinvGamma MT2EstimateZinvGamma::operator+( const MT2EstimateZinvGamma& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateZinvGamma::operator+] ERROR! Can't add MT2EstimateZinvGamma with different MT2Regions!" << std::endl;
    exit(113);
  }

  MT2EstimateZinvGamma result(*this);
  result.yield->Add(rhs.yield);

  result.iso->Add(rhs.iso);
  result.sietaieta->Add(rhs.sietaieta);

  for( unsigned i=0; i<iso_bins.size(); ++i ) {

    result.iso_bins[i]->append( *(rhs.iso_bins[i]) );
    result.iso_bins_hist[i]->Add( rhs.iso_bins_hist[i] );

  }
  
  return result;

}



MT2EstimateZinvGamma MT2EstimateZinvGamma::operator-( const MT2EstimateZinvGamma& rhs ) const{


  if( *(this->region) != *(rhs.region) ) {
    std::cout << "[MT2EstimateZinvGamma::operator-] ERROR! Can't add MT2EstimateZinvGamma with different MT2Regions!" << std::endl;
    exit(113);
  }

  std::cout << "[MT2EstimateZinvGamma::operator-] CAREFUL!! RooDataSets will not be subtracted but appended!!" << std::endl;

  MT2EstimateZinvGamma result(*this);

  result.yield->Add(rhs.yield, -1.);
  result.iso->Add(rhs.iso, -1.);
  result.sietaieta->Add(rhs.sietaieta, -1.);

  for( unsigned i=0; i<iso_bins.size(); ++i ) {

    result.iso_bins[i]->append( *(rhs.iso_bins[i]) ); // should put negative weights!
    result.iso_bins_hist[i]->Add( rhs.iso_bins_hist[i], -1. );

  }
  
  return result;

}




const MT2EstimateZinvGamma& MT2EstimateZinvGamma::operator+=( const MT2EstimateZinvGamma& rhs ) {

  this->yield->Add(rhs.yield);

  this->iso->Add(rhs.iso);
  this->sietaieta->Add(rhs.sietaieta);

  for( unsigned i=0; i<iso_bins.size(); ++i ) {

    this->iso_bins[i]->append( *(rhs.iso_bins[i]) );
    this->iso_bins_hist[i]->Add( rhs.iso_bins_hist[i] );

  }

  return (*this);

}




const MT2EstimateZinvGamma& MT2EstimateZinvGamma::operator-=( const MT2EstimateZinvGamma& rhs ) {

  this->yield->Add(rhs.yield, -1.);

  this->iso->Add(rhs.iso, -1.);
  this->sietaieta->Add(rhs.sietaieta, -1.);

  for( unsigned i=0; i<iso_bins.size(); ++i ) {

    this->iso_bins[i]->append( *(rhs.iso_bins[i]) );
    this->iso_bins_hist[i]->Add( rhs.iso_bins_hist[i], -1. );

  }

  return (*this);

}



MT2EstimateZinvGamma MT2EstimateZinvGamma::operator*( float k ) const{


  MT2EstimateZinvGamma result(*this);
  result.yield->Scale(k);
  result.iso->Scale(k);
  result.sietaieta->Scale(k);

  std::vector<RooDataSet*> new_iso_bins;

  for( unsigned i=0; i<iso_bins.size(); ++i ) {

    result.iso_bins_hist[i]->Scale(k);

    std::string oldName(result.iso_bins[i]->GetName());
    RooDataSet* newDataset = new RooDataSet( "newDataset", "", RooArgSet(*x_,*w_), w_->GetName() );

    for (int iEntry = 0; iEntry < result.iso_bins[i]->numEntries(); iEntry++) {
      x_->setVal(result.iso_bins[i]->get(iEntry)->getRealValue("x"));
      float newWeight = k* result.iso_bins[i]->weight();
      w_->setVal(newWeight);
      newDataset->add( RooArgList(*x_, *w_), newWeight );
    }

    delete result.iso_bins[i];
    result.iso_bins[i] = 0;

    newDataset->SetName(oldName.c_str());
    new_iso_bins.push_back(newDataset);

  }

  
  result.iso_bins.clear();
  
  for( unsigned i=0; i<new_iso_bins.size(); ++i )
    result.iso_bins.push_back(new_iso_bins[i]);

  return result;

}


  



const MT2EstimateZinvGamma& MT2EstimateZinvGamma::operator*=( float k ) {

  this->yield->Scale(k);
  this->iso->Scale(k);
  this->sietaieta->Scale(k);

  std::vector<RooDataSet*> new_iso_bins;

  for( unsigned i=0; i<iso_bins.size(); ++i ) {

    this->iso_bins_hist[i]->Scale(k);

    std::string oldName(this->iso_bins[i]->GetName());
    RooDataSet* newDataset = new RooDataSet( "newDataset", "", RooArgSet(*x_,*w_), w_->GetName() );

    for (int iEntry = 0; iEntry < this->iso_bins[i]->numEntries(); iEntry++) {
      x_->setVal(this->iso_bins[i]->get(iEntry)->getRealValue("x"));
      float newWeight = k* this->iso_bins[i]->weight();
      w_->setVal(newWeight);
      newDataset->add( RooArgList(*x_, *w_), newWeight );
    }

    delete this->iso_bins[i];
    this->iso_bins[i] = 0;

    newDataset->SetName(oldName.c_str());
    new_iso_bins.push_back(newDataset);

  }

  this->iso_bins.clear();
  
  for( unsigned i=0; i<new_iso_bins.size(); ++i )
    this->iso_bins.push_back(new_iso_bins[i]);

  return (*this);

}




MT2EstimateZinvGamma MT2EstimateZinvGamma::operator/( float k ) const{


  MT2EstimateZinvGamma result(*this);
  result.yield->Scale(1./k);
  result.iso->Scale(1./k);
  result.sietaieta->Scale(1./k);

  std::vector<RooDataSet*> new_iso_bins;

  for( unsigned i=0; i<iso_bins.size(); ++i ) {

    result.iso_bins_hist[i]->Scale(1./k);

    std::string oldName(result.iso_bins[i]->GetName());
    RooDataSet* newDataset = new RooDataSet( "newDataset", "", RooArgSet(*x_,*w_), w_->GetName() );

    for (int iEntry = 0; iEntry < result.iso_bins[i]->numEntries(); iEntry++) {
      x_->setVal(result.iso_bins[i]->get(iEntry)->getRealValue("x"));
      float newWeight = result.iso_bins[i]->weight()/k;
      w_->setVal(newWeight);
      newDataset->add( RooArgList(*x_, *w_), newWeight );
    }

    delete result.iso_bins[i];
    result.iso_bins[i] = 0;

    newDataset->SetName(oldName.c_str());
    new_iso_bins.push_back(newDataset);

  }

  
  result.iso_bins.clear();
  
  for( unsigned i=0; i<new_iso_bins.size(); ++i )
    result.iso_bins.push_back(new_iso_bins[i]);

  return result;

}




const MT2EstimateZinvGamma& MT2EstimateZinvGamma::operator/=( float k ) {

  this->yield->Scale(1./k);
  this->iso->Scale(1./k);
  this->sietaieta->Scale(1./k);

  std::vector<RooDataSet*> new_iso_bins;

  for( unsigned i=0; i<iso_bins.size(); ++i ) {

    this->iso_bins_hist[i]->Scale(1./k);

    std::string oldName(this->iso_bins[i]->GetName());
    RooDataSet* newDataset = new RooDataSet( "newDataset", "", RooArgSet(*x_,*w_), w_->GetName() );

    for (int iEntry = 0; iEntry < this->iso_bins[i]->numEntries(); iEntry++) {
      x_->setVal(this->iso_bins[i]->get(iEntry)->getRealValue("x"));
      float newWeight = this->iso_bins[i]->weight()/k;
      w_->setVal(newWeight);
      newDataset->add( RooArgList(*x_, *w_), newWeight );
    }

    delete this->iso_bins[i];
    this->iso_bins[i] = 0;

    newDataset->SetName(oldName.c_str());
    new_iso_bins.push_back(newDataset);

  }

  this->iso_bins.clear();
  
  for( unsigned i=0; i<new_iso_bins.size(); ++i )
    this->iso_bins.push_back(new_iso_bins[i]);

  return (*this);

}

  
  
  
  
// friend functions:

MT2EstimateZinvGamma operator*( float k, const MT2EstimateZinvGamma& rhs ) {

  return rhs*k;

}



MT2EstimateZinvGamma operator/( float k, const MT2EstimateZinvGamma& rhs ) {

  return rhs/k;

}


