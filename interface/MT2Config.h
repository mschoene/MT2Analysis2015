#ifndef MT2Config_h
#define MT2Config_h


#include "MT2Region.h"
#include "TH1D.h"
#include "TFile.h"

#include <string>



// this is the basic Estimate class: 
// it refers to a region, and it has a yield histogram
// other more complex classes (like LostLepton estimates)
// should inherit from this one, 
// and add further specialized data members




class MT2Config {
 public:
  MT2Config( const std::string& configFileName );
  std::string regionsSet()      const { return regionsSet_; };
  std::string mcSamples()       const { return mcSamples_; };
  std::string sigSamples()      const { return sigSamples_; };
  std::string dataSamples()     const { return dataSamples_; };
  std::string lostLeptonTag()   const { return lostLeptonTag_; };
  std::string qcdTag()          const { return qcdTag_; };
  std::string zinvTag()         const { return zinvTag_; };
  std::string additionalStuff() const { return additionalStuff_; };

  bool useMC() {
    bool useEstimates = lostLeptonTag_!="" && qcdTag_!="" && zinvTag_!="";
    return !useEstimates; }

 private:

  std::string regionsSet_;
  std::string mcSamples_;
  std::string sigSamples_;
  std::string dataSamples_;
  std::string lostLeptonTag_;
  std::string qcdTag_;
  std::string zinvTag_;
  std::string additionalStuff_;
};

#endif
