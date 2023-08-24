#ifndef TauAnalysis_Entanglement_StartPosAlgo2_h
#define TauAnalysis_Entanglement_StartPosAlgo2_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"          // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/StartPosAlgoBase.h" // StartPosAlgoBase

#include <vector>                                                // std::vector<>

class StartPosAlgo2 : public StartPosAlgoBase
{
 public:
  explicit StartPosAlgo2(const edm::ParameterSet&);
  ~StartPosAlgo2();
   
  std::vector<KinematicEvent>
  operator()(const KinematicEvent& kineEvt);

 private:
  bool applyRecoilEnergy_and_PzConstraint_;
};

#endif // TauAnalysis_Entanglement_StartPosAlgo2_h
