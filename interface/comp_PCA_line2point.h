#ifndef TauAnalysis_Entanglement_comp_PCA_line2point_h
#define TauAnalysis_Entanglement_comp_PCA_line2point_h

#include "DataFormats/Candidate/interface/Candidate.h"            // reco::Candidate::Vector, reco::Candidate::Point

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3

reco::Candidate::Point
comp_PCA_line2point(const reco::Candidate::Point& P1, const reco::Candidate::Vector& V1,
                    const reco::Candidate::Point& P2,
                    const math::Matrix3x3* cov = nullptr,
                    double lambdaMin = -1.e+6, double lambdaMax = +1.e+6, 
                    int verbosity = -1);

#endif // TauAnalysis_Entanglement_comp_PCA_line2point_h

