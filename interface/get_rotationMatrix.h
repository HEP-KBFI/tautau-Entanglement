#ifndef TauAnalysis_Entanglement_get_rotationMatrix_h
#define TauAnalysis_Entanglement_get_rotationMatrix_h

#include "DataFormats/Candidate/interface/Candidate.h"            // reco::Candidate::Vector

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3

math::Matrix3x3
get_rotationMatrix(const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k);

math::Matrix3x3
get_rotationMatrixInv(const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k);

#endif // TauAnalysis_Entanglement_get_rotationMatrix_h
