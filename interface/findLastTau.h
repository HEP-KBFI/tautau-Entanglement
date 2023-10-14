#ifndef TauAnalysis_Entanglement_findLastTau_h
#define TauAnalysis_Entanglement_findLastTau_h

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle
#pragma GCC diagnostic pop

const reco::GenParticle*
findLastTau(const reco::GenParticle* mother);

#endif // TauAnalysis_Entanglement_findLastTau_h
