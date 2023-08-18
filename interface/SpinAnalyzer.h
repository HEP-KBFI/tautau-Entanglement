#ifndef TauAnalysis_Entanglement_SpinAnalyzer_h
#define TauAnalysis_Entanglement_SpinAnalyzer_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"               // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/EntanglementDataset.h"   // EntanglementDataset
#include "TauAnalysis/Entanglement/interface/Measurement.h"           // spin::Measurement
#include "TauAnalysis/Entanglement/interface/SpinAlgoBase.h"          // SpinAlgoBase

#include <TRandom3.h>                                                 // TRandom3

#include <vector>                                                     // std::vector<>

class SpinAnalyzer
{
 public:
  SpinAnalyzer(const edm::ParameterSet& cfg);
  ~SpinAnalyzer();

  spin::Measurement
  operator()(const EntanglementDataset& dataset);

 protected:
  std::string spinAnalyzer_;
  SpinAlgoBase* algo_;
  unsigned numBootstrapSamples_;
  int maxEvents_afterCuts_;

  TRandom3 rnd_;

  int verbosity_;
};

#endif // TauAnalysis_Entanglement_SpinAnalyzer_h
