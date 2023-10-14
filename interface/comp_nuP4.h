#ifndef TauAnalysis_Entanglement_comp_nuP4_h
#define TauAnalysis_Entanglement_comp_nuP4_h

#include "DataFormats/Candidate/interface/Candidate.h"        // reco::Candidate::LorentzVector
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle
#pragma GCC diagnostic pop

#include <vector>                                             // std::vector<>

reco::Candidate::LorentzVector
comp_nuP4(const std::vector<const reco::GenParticle*>& decayProducts);

#endif // TauAnalysis_Entanglement_comp_nuP4_h
