#ifndef TauAnalysis_Entanglement_get_tipPerp_h
#define TauAnalysis_Entanglement_get_tipPerp_h

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::LorentzVector, reco::Candidate::Point

double
get_tipPerp(const reco::Candidate::LorentzVector& tauP4, 
            const reco::Candidate::LorentzVector& visTauP4, 
            const reco::Candidate::Point& pv, const reco::Candidate::Point& tipPCA);

#endif // TauAnalysis_Entanglement_get_tipPerp_h
