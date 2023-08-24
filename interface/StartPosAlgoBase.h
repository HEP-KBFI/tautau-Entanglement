#ifndef TauAnalysis_Entanglement_StartPosAlgoBase_h
#define TauAnalysis_Entanglement_StartPosAlgoBase_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent
#include "TauAnalysis/Entanglement/interface/Resolutions.h"    // Resolutions

#include <vector>                                              // std::vector<>

class StartPosAlgoBase
{
 public:
  explicit StartPosAlgoBase(const edm::ParameterSet& cfg);
  virtual ~StartPosAlgoBase();
   
  int
  get_algo() const;

  virtual
  std::vector<KinematicEvent>
  operator()(const KinematicEvent& kineEvt) = 0;
 
 protected:
  int algo_;

  Resolutions* resolutions_;

  int collider_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_StartPosAlgoBase_h
