#ifndef TauAnalysis_Entanglement_KinematicFitOneProng0Pi0_h
#define TauAnalysis_Entanglement_KinematicFitOneProng0Pi0_h

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"    // KinematicEvent
#include "TauAnalysis/Entanglement/interface/KinematicFit.h"      // KinematicFit
#include "TauAnalysis/Entanglement/interface/KinematicParticle.h" // KinematicParticle

class KinematicFitOneProng0Pi0 : public KinematicFit
{
 public:
  KinematicFitOneProng0Pi0();
  ~KinematicFitOneProng0Pi0();

  KinematicEvent
  operator()(const KinematicEvent& evt);

 protected:
  void
  findStartPosition();

  const KinematicParticle* piMinus_;
  const KinematicParticle* piPlus_;
};

#endif // TauAnalysis_Entanglement_KinematicFitOneProng0Pi0_h
