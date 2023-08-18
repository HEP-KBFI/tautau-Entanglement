#ifndef TauAnalysis_Entanglement_PolarimetricVectorAlgoThreeProng0Pi0_h
#define TauAnalysis_Entanglement_PolarimetricVectorAlgoThreeProng0Pi0_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"                    // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"             // KinematicEvent
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoBase.h" // PolarimetricVectorAlgoBase

class PolarimetricVectorAlgoThreeProng0Pi0 : public PolarimetricVectorAlgoBase
{
 public:
  PolarimetricVectorAlgoThreeProng0Pi0(const edm::ParameterSet& cfg);
  ~PolarimetricVectorAlgoThreeProng0Pi0();

  reco::Candidate::Vector
  operator()(const KinematicEvent& evt, int tau);
};

#endif // TauAnalysis_Entanglement_PolarimetricVectorAlgoThreeProng0Pi0_h
