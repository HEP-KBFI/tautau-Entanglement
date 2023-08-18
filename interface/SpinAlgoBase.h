#ifndef TauAnalysis_Entanglement_SpinAlgoBase_h
#define TauAnalysis_Entanglement_SpinAlgoBase_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"             // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/EntanglementDataset.h" // EntanglementDataset
#include "TauAnalysis/Entanglement/interface/Measurement.h"         // spin::Measurement

class SpinAlgoBase
{
 public:
  SpinAlgoBase(const edm::ParameterSet& cfg);
  virtual ~SpinAlgoBase();

  virtual
  spin::Measurement
  operator()(const EntanglementDataset& dataset) = 0;

 protected:
  int verbosity_;
};

#endif // TauAnalysis_Entanglement_SpinAlgoBase_h
