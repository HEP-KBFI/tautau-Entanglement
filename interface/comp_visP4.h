#ifndef TauAnalysis_Entanglement_comp_visP4_h
#define TauAnalysis_Entanglement_comp_visP4_h

#include "DataFormats/Candidate/interface/Candidate.h"        // reco::Candidate::LorentzVector
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle

#include <vector>                                             // std::vector<>

reco::Candidate::LorentzVector
comp_visP4(const std::vector<const reco::GenParticle*>& decayProducts);

#endif // TauAnalysis_Entanglement_comp_visP4_h
