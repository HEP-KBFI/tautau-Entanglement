#ifndef TauAnalysis_Entanglement_GenKinematicEventBuilder_h
#define TauAnalysis_Entanglement_GenKinematicEventBuilder_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"      // reco::GenParticle
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"   // reco::GenParticleCollection

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"     // KinematicEvent
#include "TauAnalysis/Entanglement/interface/Resolutions.h"        // Resolutions
#include "TauAnalysis/Entanglement/interface/Smearing.h"           // Smearing
#include "TauAnalysis/Entanglement/interface/PolarimetricVector.h" // PolarimetricVector

class GenKinematicEventBuilder
{
 public:
  explicit GenKinematicEventBuilder(const edm::ParameterSet& cfg);
  ~GenKinematicEventBuilder();
   
  KinematicEvent
  operator()(const reco::GenParticleCollection& genParticles);
 
 private:
  Resolutions* resolutions_;

  Smearing smearing_;
  bool applySmearing_;

  PolarimetricVector polarimetricVector_;

  int collider_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_GenKinematicEventBuilder_h
