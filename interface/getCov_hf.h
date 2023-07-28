#ifndef TauAnalysis_Entanglement_getCov_hf_h
#define TauAnalysis_Entanglement_getCov_hf_h

#include "DataFormats/Candidate/interface/Candidate.h"            // Candidate::Vector

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3

math::Matrix3x3
getCov_hf(double dk, double dr, double dn,
          const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k);

#endif // TauAnalysis_Entanglement_getCov_hf_h
