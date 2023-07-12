#ifndef TauAnalysis_Entanglement_getP4_ttrf_hf_trf_h
#define TauAnalysis_Entanglement_getP4_ttrf_hf_trf_h

#include "DataFormats/Candidate/interface/Candidate.h" // Candidate::LorentzVector, Candidate::Vector

#include <Math/Boost.h>                                // ROOT::Math::Boost

reco::Candidate::LorentzVector
getP4_ttrf_hf_trf(const reco::Candidate::LorentzVector& p4,
                  const ROOT::Math::Boost& boost_ttrf,
                  const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                  const ROOT::Math::Boost& boost_trf);

#endif // TauAnalysis_Entanglement_getP4_ttrf_hf_trf_h
