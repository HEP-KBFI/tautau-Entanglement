#ifndef TauAnalysis_Entanglement_SpinAnalyzerBase_h
#define TauAnalysis_Entanglement_SpinAnalyzerBase_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent

class SpinAnalyzerBase
{
 public:
  SpinAnalyzerBase(const edm::ParameterSet& cfg);
  virtual ~SpinAnalyzerBase();

  enum { kTauPlus, kTauMinus };

  virtual
  reco::Candidate::Vector
  operator()(const KinematicEvent& kineEvt, int tau) = 0;

 protected:
  int hAxis_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_SpinAnalyzerBase_h
