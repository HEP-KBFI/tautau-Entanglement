#ifndef TauAnalysis_Entanglement_StartPosAlgo1_h
#define TauAnalysis_Entanglement_StartPosAlgo1_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"          // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/Resolutions.h"      // Resolutions
#include "TauAnalysis/Entanglement/interface/StartPosAlgoBase.h" // StartPosAlgoBase

class StartPosAlgo1 : public StartPosAlgoBase
{
 public:
  explicit StartPosAlgo1(const edm::ParameterSet&);
  ~StartPosAlgo1();
   
  KinematicEvent
  operator()(const KinematicEvent& kineEvt);
 
 private:
  Resolutions* resolutions_;

  bool applyHiggsMassConstraint_;
};

#endif // TauAnalysis_Entanglement_StartPosAlgo1_h
