#ifndef TauAnalysis_Entanglement_SpinAnalyzer_h
#define TauAnalysis_Entanglement_SpinAnalyzer_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"                    // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"             // KinematicEvent
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerOneProng0Pi0.h"   // SpinAnalyzerOneProng0Pi0
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerOneProng1Pi0.h"   // SpinAnalyzerOneProng1Pi0
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerThreeProng0Pi0.h" // SpinAnalyzerThreeProng0Pi0

class SpinAnalyzer
{
 public:
  SpinAnalyzer(const edm::ParameterSet& cfg);
  ~SpinAnalyzer();

  reco::Candidate::Vector
  operator()(const KinematicEvent& evt, int tau);

 protected:
  SpinAnalyzerOneProng0Pi0   spinAnalyzerOneProng0Pi0_;
  SpinAnalyzerOneProng1Pi0   spinAnalyzerOneProng1Pi0_;
  SpinAnalyzerThreeProng0Pi0 spinAnalyzerThreeProng0Pi0_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_SpinAnalyzer_h
