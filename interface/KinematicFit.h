#ifndef TauAnalysis_Entanglement_KinematicFit_h
#define TauAnalysis_Entanglement_KinematicFit_h

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent

class KinematicFit
{
 public:
  KinematicFit();
  virtual
  ~KinematicFit();

  virtual
  KinematicEvent
  operator()(const KinematicEvent& evt) = 0;

 protected:
  virtual
  void
  findStartPosition() = 0;
};

#endif // TauAnalysis_Entanglement_KinematicFit_h
