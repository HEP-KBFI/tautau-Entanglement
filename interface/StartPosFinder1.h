#ifndef TauAnalysis_Entanglement_StartPosFinder1_h
#define TauAnalysis_Entanglement_StartPosFinder1_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"            // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/Resolutions.h"        // Resolutions
#include "TauAnalysis/Entanglement/interface/StartPosFinderBase.h" // StartPosFinderBase

class StartPosFinder1 : public StartPosFinderBase
{
 public:
  explicit StartPosFinder1(const edm::ParameterSet&);
  ~StartPosFinder1();
   
  KinematicEvent
  operator()(const KinematicEvent& kineEvt);
 
 private:
  Resolutions* resolutions_;

  bool applyHiggsMassConstraint_;
};

#endif // TauAnalysis_Entanglement_StartPosFinder1_h
