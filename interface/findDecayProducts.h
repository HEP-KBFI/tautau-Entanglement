#ifndef TauAnalysis_Entanglement_findDecayProducts_h
#define TauAnalysis_Entanglement_findDecayProducts_h

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle
#pragma GCC diagnostic pop

#include <vector>                                             // std::vector<>

void
findDecayProducts(const reco::GenParticle* mother, 
                  std::vector<const reco::GenParticle*>& daughters,
                  bool treatNeutralPionsAsStable = true);

#endif // TauAnalysis_Entanglement_findDecayProducts_h
