#ifndef TauAnalysis_Entanglement_findLastTau_h
#define TauAnalysis_Entanglement_findLastTau_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle

const reco::GenParticle*
findLastTau(const reco::GenParticle* mother);

#endif // TauAnalysis_Entanglement_findLastTau_h
