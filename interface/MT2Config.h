#ifndef MT2Config_h
#define MT2Config_h

#include <string>

class MT2Config {

 public:

  MT2Config( const std::string& configFileName );

  std::string name() const { return name_; };

  float lumi()      const { return lumi_; };
  float lumi_JetHT()         const;
  float lumi_HTMHT()         const;
  float lumi_SinglePhoton()  const;
  float lumi_DoubleEG()      const;
  float lumi_DoubleMu()      const;

 
  std::string regionsSet()             const { return regionsSet_; };
  std::string mcSamples()              const { return mcSamples_; };
  std::string qcdMCSamples()           const { return qcdMCSamples_; };
  std::string sigSamples()             const { return sigSamples_; };
  std::string dataSamples()            const { return dataSamples_; };
  std::string qcdDataSamples()         const { return qcdDataSamples_; };
  std::string qcdMonoJetMCSamples()    const { return qcdMonoJetMCSamples_; };
  std::string qcdMonoJetDataSamples()  const { return qcdMonoJetDataSamples_; };
  std::string additionalStuff()        const { return additionalStuff_; };
  std::string analysisType()           const { return analysisType_; };
  std::string crRegionsSet()           const { return crRegionsSet_; };
  std::string qcdRegionsSet()          const { return qcdRegionsSet_; };

  std::string gammaTemplateType() const { return gammaTemplateType_; };
  std::string gammaTemplateRegions() const { return gammaTemplateRegions_; };
  float gammaIsoCut() const { return gammaIsoCut_; };
  std::string gamma2bMethod() const { return gamma2bMethod_; };

  std::string zllRegions() const { return zllRegions_; };

  void set_gammaTemplateType( const std::string& newType ) { gammaTemplateType_ = newType; };

  bool useMC() const;

  bool dummyAnalysis() const;

  std::string getEventYieldDir() const;
  std::string getGammaCRdir() const;
  std::string getDiPhotonCRdir() const;

  void saveAs( const std::string& filename ) const;

 private:

  float defaultLumi( float lumi) const;

  std::string name_;

  float lumi_;
  float lumi_JetHT_;
  float lumi_HTMHT_;
  float lumi_SinglePhoton_;
  float lumi_DoubleEG_;
  float lumi_DoubleMu_;

  std::string regionsSet_;
  std::string mcSamples_;
  std::string qcdMCSamples_;
  std::string sigSamples_;
  std::string dataSamples_;
  std::string qcdDataSamples_;
  std::string qcdMonoJetDataSamples_;
  std::string qcdMonoJetMCSamples_;
  std::string additionalStuff_;
  std::string analysisType_;
  std::string crRegionsSet_;
  std::string qcdRegionsSet_;

  std::string gammaTemplateType_;
  std::string gammaTemplateRegions_;
  float gammaIsoCut_;
  std::string gamma2bMethod_;

  std::string zllRegions_;

};



#endif
