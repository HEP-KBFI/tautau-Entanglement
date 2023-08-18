#ifndef TauAnalysis_Entanglement_PolarimetricVector_h
#define TauAnalysis_Entanglement_PolarimetricVector_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"                              // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"                       // KinematicEvent
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoOneProng0Pi0.h"   // PolarimetricVectorAlgoOneProng0Pi0
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoOneProng1Pi0.h"   // PolarimetricVectorAlgoOneProng1Pi0
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoThreeProng0Pi0.h" // PolarimetricVectorAlgoThreeProng0Pi0

class PolarimetricVector
{
 public:
  PolarimetricVector(const edm::ParameterSet& cfg);
  ~PolarimetricVector();

  reco::Candidate::Vector
  operator()(const KinematicEvent& evt, int tau);

 protected:
  PolarimetricVectorAlgoOneProng0Pi0   algoOneProng0Pi0_;
  PolarimetricVectorAlgoOneProng1Pi0   algoOneProng1Pi0_;
  PolarimetricVectorAlgoThreeProng0Pi0 algoThreeProng0Pi0_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_PolarimetricVector_h
