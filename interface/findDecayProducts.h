#ifndef TauAnalysis_Entanglement_findDecayProducts_h
#define TauAnalysis_Entanglement_findDecayProducts_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle

#include <vector>                                             // std::vector<>

void
findDecayProducts(const reco::GenParticle* mother, 
                  std::vector<const reco::GenParticle*>& daughters);

#endif // TauAnalysis_Entanglement_findDecayProducts_h
