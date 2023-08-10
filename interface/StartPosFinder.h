#ifndef TauAnalysis_Entanglement_StartPosFinder_h
#define TauAnalysis_Entanglement_StartPosFinder_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"            // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"     // KinematicEvent
#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"       // SpinAnalyzer
#include "TauAnalysis/Entanglement/interface/StartPosFinderBase.h" // StartPosFinderBase

class StartPosFinder
{
 public:
  explicit StartPosFinder(const edm::ParameterSet&);
  ~StartPosFinder();
   
  KinematicEvent
  operator()(const KinematicEvent& kineEvt);
 
 private:
  StartPosFinderBase* algo_;

  int mode_;

  SpinAnalyzer spinAnalyzer_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_StartPosFinder_h
