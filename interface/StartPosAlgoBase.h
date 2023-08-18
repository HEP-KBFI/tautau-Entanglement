#ifndef TauAnalysis_Entanglement_StartPosAlgoBase_h
#define TauAnalysis_Entanglement_StartPosAlgoBase_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent

class StartPosAlgoBase
{
 public:
  explicit StartPosAlgoBase(const edm::ParameterSet& cfg);
  virtual ~StartPosAlgoBase();
   
  virtual
  KinematicEvent
  operator()(const KinematicEvent& kineEvt) = 0;
 
 protected:
  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_StartPosAlgoBase_h
