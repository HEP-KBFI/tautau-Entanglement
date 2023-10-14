#ifndef TauAnalysis_Entanglement_printGenParticles_h
#define TauAnalysis_Entanglement_printGenParticles_h

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle
#pragma GCC diagnostic pop

#include <vector>                                             // std::vector<>

void
printGenParticles(const std::vector<const reco::GenParticle*>& decayProducts,
                  bool cartesian = true);

#endif // TauAnalysis_Entanglement_printGenParticles_h

