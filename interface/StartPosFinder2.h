#ifndef TauAnalysis_Entanglement_StartPosFinder2_h
#define TauAnalysis_Entanglement_StartPosFinder2_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"            // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/Resolutions.h"        // Resolutions
#include "TauAnalysis/Entanglement/interface/StartPosFinderBase.h" // StartPosFinderBase

class StartPosFinder2 : public StartPosFinderBase
{
 public:
  explicit StartPosFinder2(const edm::ParameterSet&);
  ~StartPosFinder2();
   
  KinematicEvent
  operator()(const KinematicEvent& kineEvt);

 private:
  Resolutions* resolutions_;

  bool applyHiggsMassConstraint_;
};

#endif // TauAnalysis_Entanglement_StartPosFinder2_h
