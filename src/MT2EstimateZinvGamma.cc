#include "../interface/MT2EstimateZinvGamma.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TRandom3.h"
#include "TTree.h"









//some function to create a newly binned MT2EstimateZinGamma
void MT2EstimateZinvGamma::rebinYields( MT2Analysis<MT2EstimateZinvGamma>* thisEstimate, int nBins, Double_t* bins ){
  std::set<MT2Region> regions = thisEstimate->getRegions();
  for( std::set<MT2Region>::iterator iR = regions.begin(); iR!=regions.end(); ++iR ) {

    MT2EstimateZinvGamma* estimate = thisEstimate->get(*iR);
    TH1D* thisYield = estimate->yield;
    std::string oldName(thisYield->GetName());
    delete thisYield;
    thisYield = new TH1D( oldName.c_str(), "", nBins, bins );

    std::vector<TH1D*> this_iso_bins_hist = estimate->iso_bins_hist;
    unsigned int old_size =    this_iso_bins_hist.size();
    std::vector<RooDataSet*> this_iso_bins = estimate->iso_bins;
 

    std::string oldNames_hist[old_size];
    std::string oldNames[old_size];

    for( unsigned int i=0; i<old_size; ++i ){
      std::string temp(this_iso_bins[i]->GetName());
        oldNames[i] = temp;
      std::string temp_hist(this_iso_bins_hist[i]->GetName());
        oldNames_hist[i] = temp_hist;
    }   

    for( unsigned i=0; i<old_size; ++i ) { 
      delete this_iso_bins_hist[i];
      delete this_iso_bins[i];
    }
    this_iso_bins_hist.clear();
    this_iso_bins.clear();
  

    int nbins_mt2 = nBins;
    Double_t bins_mt2[nBins];
    for(int i = 0; i<nBins; i++)
      bins_mt2[i] = bins[i];
 

    int nbins = 8;
    int xmax = 10.;


    RooRealVar* this_x_ = estimate->x_;
    delete this_x_;
    RooRealVar* this_w_ = estimate->w_;
    delete this_w_;
    this_x_ = new RooRealVar( "x", "", 0., xmax );
    this_w_ = new RooRealVar( "w", "", 0., 1000. );

    for( int i=0; i<nbins_mt2; ++i ) {

      // RooDataSet* isoDataset = new RooDataSet( oldNames[i].c_str(), "", RooArgSet(*x_,*w_), w_->GetName() );
      //  iso_bins.push_back(isoDataset);
      //TH1D* this_iso_hist = new TH1D( this->getHistoName(Form("iso_bin%d_hist", i)).c_str() , "", nbins-1, bins );
      RooDataSet* this_isoDataset = new RooDataSet( oldNames[i].c_str() , "", RooArgSet(*this_x_,*this_w_), this_w_->GetName() );
      this_iso_bins.push_back(this_isoDataset);
    std::cout << "FUCKKKKKKKKKKKK" << std::endl;

      TH1D* this_iso_hist = new TH1D( oldNames_hist[i].c_str() , "", nbins, 0., xmax );
      this_iso_hist->Sumw2();
      this_iso_bins_hist.push_back(this_iso_hist);
  }
  
   
  }//end regions

}





MT2Analysis<MT2EstimateZinvGamma>*  MT2EstimateZinvGamma::makeInclusiveEstimateFromInclusiveTree( const std::string& aname, MT2Analysis<MT2EstimateTree>* analysis, const std::string& selectionTree, const std::string& var, int nBins, Double_t* bins ){

  std::set<MT2Region> regions = analysis->getRegions();
 if( regions.size()!=1 ) {
    std::cout << "[MT2EstimateTree::makeAnalysisFromEstimateTreeInclusive] ERROR!! You need to pass an inclusive MT2EstimateTree Analysis to use this function!" << std::endl;
    exit(19191);
  }

  MT2EstimateTree* treeInclusive = analysis->get( *(regions.begin()) );


  // will create a new analysis with same region as original one (i.e. inclusive)
  //BUT different binning for the yield and iso histogram
  std::set<MT2Region> newRegions = analysis->getRegions();
 
  //MT2Analysis<MT2EstimateZinvGamma>* analysis = new MT2Analysis<MT2EstimateTree>( aname, regionsSet );
  MT2Analysis<MT2EstimateZinvGamma>* thisEstimate  = new MT2Analysis<MT2EstimateZinvGamma>( aname, newRegions );//  = analysis->get( *(regions.begin())  );

  if ( nBins!=0 ){
    MT2EstimateZinvGamma::rebinYields( (MT2Analysis<MT2EstimateZinvGamma>*)thisEstimate, nBins, bins );
  }

  MT2EstimateZinvGamma* theEst = thisEstimate->get( *(regions.begin()));

  //fill the new iso with the treeee
  int nentries = treeInclusive->tree->GetEntries();

 
  Float_t ht;
  treeInclusive->tree->SetBranchAddress("ht", &ht);
  Float_t weight;
  treeInclusive->tree->SetBranchAddress("weight", &weight);
  Float_t iso;
  treeInclusive->tree->SetBranchAddress("iso", &iso);


  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    treeInclusive->tree->GetEntry(iEntry);

    if( iEntry % 5000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    std::cout << "FUUUUUUUUUCK1 " << std::endl;

    theEst->fillIso( iso, weight, ht );
    std::cout << "FUUUUUUUUUCK2 " << std::endl;

 
  }

  thisEstimate->finalize();

  //done and over
  return thisEstimate;
}









MT2EstimateZinvGamma::MT2EstimateZinvGamma( const std::string& aname, const MT2Region& aregion ) : MT2Estimate( aname, aregion ) {


  sietaieta = new TH1D( this->getHistoName("sietaieta").c_str(), "", 200, 0., 0.035 );
  sietaieta->Sumw2();


  int nbins = 8;
  int xmax = 10.;

  // this histo will be used to create histogram templates:
  //iso = new TH1D( this->getHistoName("iso").c_str(), "", nbins-1, bins );
  iso = new TH1D( this->getHistoName("iso").c_str(), "", nbins, 0., xmax );
  iso->Sumw2();

  x_ = new RooRealVar( "x", "", 0., xmax );
  w_ = new RooRealVar( "w", "", 0., 1000. );

  int nbins_mt2;
  double *bins_mt2;
  aregion.getBins( nbins_mt2, bins_mt2 );


  for( int i=0; i<nbins_mt2; ++i ) {

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
std::cout << "FUU in fill " << std::endl;
 
      w_->setVal(weight);
      iso_bins[foundBin]->add( RooArgList(*x_, *w_), weight );
 
      iso_bins_hist[foundBin]->Fill( iso, weight );
 std::cout << "FUU in fill " << std::endl;
 
        
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

    for( int ibin=1; ibin<iso_bins_hist[i]->GetNbinsX()+1; ++ibin ) {

      float isoValue = iso_bins_hist[i]->GetBinLowEdge(ibin);
      //float isoValue = iso_bins_hist[i]->GetBinCenter(ibin);

      //int entries = (int)iso_bins_hist[i]->GetBinContent(ibin);
      int entries = rand.Poisson( iso_bins_hist[i]->GetBinContent(ibin) );

      for( int ientry=0; ientry<entries; ++ientry ) {

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

//MT2Estimate::addOverflowSingleHisto( iso );
//MT2Estimate::addOverflowSingleHisto( sietaieta );

//for( unsigned i=0; i<iso_bins.size(); ++i ) {
//  MT2Estimate::addOverflowSingleHisto(iso_bins_hist[i]);
//}

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


