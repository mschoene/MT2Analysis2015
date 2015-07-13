#ifndef MT2Config_h
#define MT2Config_h

#include <string>

class MT2Config {

 public:

<<<<<<< HEAD
  MT2Config( const std::string& configFileName );

  std::string name() const { return name_; };

  float lumi()      const { return lumi_; };
=======
// this is the basic Estimate class: 
// it refers to a region, and it has a yield histogram
// other more complex classes (like LostLepton estimates)
// should inherit from this one, 
// and add further specialized data members




class MT2Config {
 public:
  MT2Config( const std::string& configFileName );
>>>>>>> e9c4691... CherryPick2, shall conflicts rise?
  std::string regionsSet()      const { return regionsSet_; };
  std::string mcSamples()       const { return mcSamples_; };
  std::string sigSamples()      const { return sigSamples_; };
  std::string dataSamples()     const { return dataSamples_; };
  std::string additionalStuff() const { return additionalStuff_; };

  std::string gammaTemplateType() const { return gammaTemplateType_; };
  std::string gammaTemplateRegions() const { return gammaTemplateRegions_; };

  bool useMC() const;

  bool dummyAnalysis() const;

  std::string getEventYieldDir() const;

  void saveAs( const std::string& filename ) const;

 private:

  std::string name_;

  float lumi_;
  std::string regionsSet_;
  std::string mcSamples_;
  std::string sigSamples_;
  std::string dataSamples_;
  std::string additionalStuff_;

  std::string gammaTemplateType_;
  std::string gammaTemplateRegions_;

};



#endif
