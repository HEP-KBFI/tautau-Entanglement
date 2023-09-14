#ifndef TauAnalysis_Entanglement_SpinAnalyzer_h
#define TauAnalysis_Entanglement_SpinAnalyzer_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"           // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/BinnedDataset.h"     // spin::BinnedDataset
#include "TauAnalysis/Entanglement/interface/BinnedMeasurement.h" // spin::BinnedMeasurement1d, spin::BinnedMeasurement2d
#include "TauAnalysis/Entanglement/interface/Dataset.h"           // spin::Dataset
#include "TauAnalysis/Entanglement/interface/Measurement.h"       // spin::Measurement
#include "TauAnalysis/Entanglement/interface/SpinAlgoBase.h"      // SpinAlgoBase

#include <TRandom3.h>                                             // TRandom3

#include <vector>                                                 // std::vector<>

namespace spin
{

class SpinAnalyzer
{
 public:
  SpinAnalyzer(const edm::ParameterSet& cfg);
  ~SpinAnalyzer();

  void
  set_verbosity(int verbosity);

  spin::Measurement
  operator()(const spin::Dataset& dataset) const;

  spin::BinnedMeasurement1d
  operator()(const spin::BinnedDataset1d& dataset) const;
  spin::BinnedMeasurement2d
  operator()(const spin::BinnedDataset2d& dataset) const;

 protected:
  spin::Measurement
  build_measurement(const spin::Dataset& dataset, int maxEvents_afterCuts, int verbosity) const;

  std::string spinAnalyzer_;
  SpinAlgoBase* algo_;
  unsigned numBootstrapSamples_;
  int maxEvents_afterCuts_;

  mutable TRandom3 rnd_;

  int verbosity_;
};

}

#endif // TauAnalysis_Entanglement_SpinAnalyzer_h
