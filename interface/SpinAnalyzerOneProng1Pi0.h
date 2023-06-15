#ifndef TauAnalysis_Entanglement_SpinAnalyzerOneProng1Pi0_h
#define TauAnalysis_Entanglement_SpinAnalyzerOneProng1Pi0_h

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent
#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"   // SpinAnalyzer

class SpinAnalyzerOneProng1Pi0 : public SpinAnalyzer
{
 public:
  SpinAnalyzerOneProng1Pi0(const edm::ParameterSet& cfg);
  ~SpinAnalyzerOneProng1Pi0();

  reco::Candidate::Vector
  operator()(const KinematicEvent& evt, int tau);
};

#endif // TauAnalysis_Entanglement_SpinAnalyzerOneProng1Pi0_h
