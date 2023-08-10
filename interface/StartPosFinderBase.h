#ifndef TauAnalysis_Entanglement_StartPosFinderBase_h
#define TauAnalysis_Entanglement_StartPosFinderBase_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent

class StartPosFinderBase
{
 public:
  explicit StartPosFinderBase(const edm::ParameterSet& cfg);
  virtual ~StartPosFinderBase();
   
  virtual
  KinematicEvent
  operator()(const KinematicEvent& kineEvt) = 0;
 
 protected:
  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_StartPosFinderBase_h
