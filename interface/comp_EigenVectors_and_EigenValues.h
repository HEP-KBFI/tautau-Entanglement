#ifndef TauAnalysis_Entanglement_comp_EigenVectors_and_EigenValues_h
#define TauAnalysis_Entanglement_comp_EigenVectors_and_EigenValues_h

#include "FWCore/Utilities/interface/EDMException.h" // edm::Exception

#include <TMatrixD.h>                                // TMatrixD
#include <TVectorD.h>                                // TVectorD

#include <utility>                                   // std::make_pair(), std::pair<>

std::vector<std::pair<TVectorD, double>>
comp_EigenVectors_and_EigenValues(const TMatrixD& symmMatrix, int verbosity = -1);

#endif // TauAnalysis_Entanglement_comp_EigenVectors_and_EigenValues_h
