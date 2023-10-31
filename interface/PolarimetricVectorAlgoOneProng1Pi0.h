#ifndef TauAnalysis_Entanglement_PolarimetricVectorAlgoOneProng1Pi0_h
#define TauAnalysis_Entanglement_PolarimetricVectorAlgoOneProng1Pi0_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"                    // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"             // KinematicEvent
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoBase.h" // PolarimetricVectorAlgoBase

class PolarimetricVectorAlgoOneProng1Pi0 : public PolarimetricVectorAlgoBase
{
 public:
  PolarimetricVectorAlgoOneProng1Pi0(const edm::ParameterSet& cfg);
  ~PolarimetricVectorAlgoOneProng1Pi0();

  reco::Candidate::Vector
  operator()(const KinematicEvent& evt, int tau) const;
};

#endif // TauAnalysis_Entanglement_PolarimetricVectorAlgoOneProng1Pi0_h
