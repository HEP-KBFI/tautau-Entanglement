#ifndef TauAnalysis_Entanglement_PolarimetricVectorAlgoBase_h
#define TauAnalysis_Entanglement_PolarimetricVectorAlgoBase_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent

namespace pol
{
  enum { kTauPlus, kTauMinus };
}

class PolarimetricVectorAlgoBase
{
 public:
  PolarimetricVectorAlgoBase(const edm::ParameterSet& cfg);
  virtual ~PolarimetricVectorAlgoBase();

  virtual
  reco::Candidate::Vector
  operator()(const KinematicEvent& kineEvt, int tau) const = 0;

 protected:
  int hAxis_;

  int collider_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_PolarimetricVectorAlgoBase_h
