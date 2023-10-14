#ifndef TauAnalysis_Entanglement_get_leadTrack_h
#define TauAnalysis_Entanglement_get_leadTrack_h

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"     // reco::GenParticle
#pragma GCC diagnostic pop

#include "TauAnalysis/Entanglement/interface/KinematicParticle.h" // KinematicParticle

#include <vector>                                                 // std::vector<>

const reco::GenParticle*
get_leadTrack(const std::vector<const reco::GenParticle*>& decayProducts);

const KinematicParticle*
get_leadTrack(const std::vector<KinematicParticle>& daughters);

#endif // TauAnalysis_Entanglement_get_leadTrack_h
