#ifndef TauAnalysis_Entanglement_SpinAlgoBase_h
#define TauAnalysis_Entanglement_SpinAlgoBase_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"     // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/Dataset.h"     // spin::Dataset
#include "TauAnalysis/Entanglement/interface/Measurement.h" // spin::Measurement

namespace spin
{

class SpinAlgoBase
{
 public:
  SpinAlgoBase(const edm::ParameterSet& cfg);
  virtual ~SpinAlgoBase();

  virtual
  void
  set_verbosity(int verbosity);

  virtual
  spin::Measurement
  operator()(const spin::Dataset& dataset) = 0;

 protected:
  void
  addEntanglementVariables(spin::Measurement& measurement);

  int verbosity_;
};

}

#endif // TauAnalysis_Entanglement_SpinAlgoBase_h
