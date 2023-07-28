#ifndef TauAnalysis_Entanglement_getP4_hf_h
#define TauAnalysis_Entanglement_getP4_hf_h

#include "DataFormats/Candidate/interface/Candidate.h" // Candidate::LorentzVector, Candidate::Vector

reco::Candidate::LorentzVector
getP4_hf(const reco::Candidate::LorentzVector& p4,
         const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k);

#endif // TauAnalysis_Entanglement_getP4_hf_h
