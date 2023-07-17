#ifndef TauAnalysis_Entanglement_SpinAnalyzerThreeProng0Pi0_h
#define TauAnalysis_Entanglement_SpinAnalyzerThreeProng0Pi0_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"          // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"   // KinematicEvent
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerBase.h" // SpinAnalyzerBase

class SpinAnalyzerThreeProng0Pi0 : public SpinAnalyzerBase
{
 public:
  SpinAnalyzerThreeProng0Pi0(const edm::ParameterSet& cfg);
  ~SpinAnalyzerThreeProng0Pi0();

  reco::Candidate::Vector
  operator()(const KinematicEvent& evt, int tau);
};

#endif // TauAnalysis_Entanglement_SpinAnalyzerThreeProng0Pi0_h
