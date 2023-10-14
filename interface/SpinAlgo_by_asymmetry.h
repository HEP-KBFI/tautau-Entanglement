#ifndef TauAnalysis_Entanglement_SpinAlgo_by_asymmetry_h
#define TauAnalysis_Entanglement_SpinAlgo_by_asymmetry_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"      // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/Dataset.h"      // spin::Dataset
#include "TauAnalysis/Entanglement/interface/Measurement.h"  // spin::Measurement
#include "TauAnalysis/Entanglement/interface/SpinAlgoBase.h" // SpinAlgoBase

namespace spin
{

class SpinAlgo_by_asymmetry : public SpinAlgoBase
{
 public:
  SpinAlgo_by_asymmetry(const edm::ParameterSet& cfg);
  ~SpinAlgo_by_asymmetry();

  spin::Measurement
  operator()(const spin::DatasetPtrs& dataset);
};

}

#endif // TauAnalysis_Entanglement_SpinAlgo_by_asymmetry_h
