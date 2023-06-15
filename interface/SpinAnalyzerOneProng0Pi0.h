#ifndef TauAnalysis_Entanglement_SpinAnalyzerOneProng0Pi0_h
#define TauAnalysis_Entanglement_SpinAnalyzerOneProng0Pi0_h

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent
#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"   // SpinAnalyzer

class SpinAnalyzerOneProng0Pi0 : public SpinAnalyzer
{
 public:
  SpinAnalyzerOneProng0Pi0(const edm::ParameterSet& cfg);
  ~SpinAnalyzerOneProng0Pi0();

  reco::Candidate::Vector
  operator()(const KinematicEvent& evt, int tau);
};

#endif // TauAnalysis_Entanglement_SpinAnalyzerOneProng0Pi0_h
