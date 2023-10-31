#ifndef TauAnalysis_Entanglement_StartPosTIPCompatibility_h
#define TauAnalysis_Entanglement_StartPosTIPCompatibility_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent

class StartPosTIPCompatibility
{
 public:
  explicit StartPosTIPCompatibility(const edm::ParameterSet& cfg);
  ~StartPosTIPCompatibility();
   
  // CV: compute compatibility of StartPosFinder solutions with transverse impact parameters,
  //     using the procedure described in the paper arXiv:hep-ph/9307269
  double
  operator()(const KinematicEvent& kineEvt);
 
 private:
  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_StartPosTIPCompatibility_h
