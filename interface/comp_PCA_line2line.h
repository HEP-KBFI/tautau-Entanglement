#ifndef TauAnalysis_Entanglement_comp_PCA_line2line_h
#define TauAnalysis_Entanglement_comp_PCA_line2line_h

#include "DataFormats/Candidate/interface/Candidate.h"            // reco::Candidate::Vector, reco::Candidate::Point

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3

std::pair<reco::Candidate::Point, reco::Candidate::Point>
comp_PCA_line2line(const reco::Candidate::Point& P1, const reco::Candidate::Vector& V1,
                   const reco::Candidate::Point& P2, const reco::Candidate::Vector& V2,
                   const math::Matrix3x3* cov = nullptr,
                   double lambda1Min = -1.e+6, double lambda1Max = +1.e+6, double lambda2Min = -1.e+6, double lambda2Max = +1.e+6, 
                   int verbosity = -1);

#endif // TauAnalysis_Entanglement_comp_PCA_line2line_h

