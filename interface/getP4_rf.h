#ifndef TauAnalysis_Entanglement_getP4_rf_h
#define TauAnalysis_Entanglement_getP4_rf_h

#include "DataFormats/Candidate/interface/Candidate.h" // Candidate::LorentzVector

#include <Math/Boost.h>                                // ROOT::Math::Boost

reco::Candidate::LorentzVector
getP4_rf(const reco::Candidate::LorentzVector& p4,
         const ROOT::Math::Boost& boost);

#endif // TauAnalysis_Entanglement_getP4_rf_h
