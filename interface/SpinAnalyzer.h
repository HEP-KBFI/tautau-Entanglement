#ifndef TauAnalysis_Entanglement_SpinAnalyzer_h
#define TauAnalysis_Entanglement_SpinAnalyzer_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent

class SpinAnalyzer
{
 public:
  SpinAnalyzer(const edm::ParameterSet& cfg);
  virtual
  ~SpinAnalyzer();

  enum { kTauPlus, kTauMinus };

  virtual
  reco::Candidate::Vector
  operator()(const KinematicEvent& evt, int tau) = 0;

 protected:
  int hAxis_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_SpinAnalyzer_h
