#ifndef TauAnalysis_Entanglement_StartPosAlgo2_h
#define TauAnalysis_Entanglement_StartPosAlgo2_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"          // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/Resolutions.h"      // Resolutions
#include "TauAnalysis/Entanglement/interface/StartPosAlgoBase.h" // StartPosAlgoBase

class StartPosAlgo2 : public StartPosAlgoBase
{
 public:
  explicit StartPosAlgo2(const edm::ParameterSet&);
  ~StartPosAlgo2();
   
  KinematicEvent
  operator()(const KinematicEvent& kineEvt);

 private:
  Resolutions* resolutions_;

  int collider_;

  bool applyRecoilEnergy_and_PzConstraint_;
};

#endif // TauAnalysis_Entanglement_StartPosAlgo2_h
