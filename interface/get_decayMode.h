#ifndef TauAnalysis_Entanglement_get_decayMode_h
#define TauAnalysis_Entanglement_get_decayMode_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle

#include <vector>                                             // std::vector<>

int
get_decayMode(const std::vector<const reco::GenParticle*>& tau_ch,
              const std::vector<const reco::GenParticle*>& tau_pi0,
              const std::vector<const reco::GenParticle*>& tau_nu);

bool
is1Prong(int decayMode);

bool
is3Prong(int decayMode);

#endif // TauAnalysis_Entanglement_get_decayMode_h
