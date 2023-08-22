#ifndef TauAnalysis_Entanglement_StartPosFinder_h
#define TauAnalysis_Entanglement_StartPosFinder_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"            // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"     // KinematicEvent
#include "TauAnalysis/Entanglement/interface/PolarimetricVector.h" // PolarimetricVector
#include "TauAnalysis/Entanglement/interface/StartPosAlgoBase.h"   // StartPosAlgoBase

class StartPosFinder
{
 public:
  explicit StartPosFinder(const edm::ParameterSet&);
  ~StartPosFinder();
   
  KinematicEvent
  operator()(const KinematicEvent& kineEvt);
 
 private:
  StartPosAlgoBase* algo_;

  PolarimetricVector polarimetricVector_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_StartPosFinder_h
