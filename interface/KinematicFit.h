#ifndef TauAnalysis_Entanglement_KinematicFit_h
#define TauAnalysis_Entanglement_KinematicFit_h

#include "DataFormats/Candidate/interface/Candidate.h"         // reco::Candidate::LorentzVector

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent
#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"   // SpinAnalyzer

class KinematicFit
{
 public:
  KinematicFit(const edm::ParameterSet& cfg);
  ~KinematicFit();

  KinematicEvent
  operator()(const KinematicEvent& evt);

 private:
  SpinAnalyzer spinAnalyzer_;

  int applyTauMassConstraint_;
  bool applyLifetimeConstraint_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_KinematicFit_h
