#ifndef TauAnalysis_Entanglement_StartPosFinder_h
#define TauAnalysis_Entanglement_StartPosFinder_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"    // reco::GenParticle
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h" // reco::GenParticleCollection

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"   // KinematicEvent
#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"     // SpinAnalyzer

class KinematicFitStartPosFinder
{
 public:
  explicit KinematicFitStartPosFinder(const edm::ParameterSet&);
  ~KinematicFitStartPosFinder();
   
  KinematicEvent
  operator()(const KinematicEvent& kineEvt);
 
 private:
  SpinAnalyzer spinAnalyzer_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_KinematicFitStartPosFinder_h
