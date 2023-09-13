#ifndef TauAnalysis_Entanglement_printEigenVectors_and_EigenValues_h
#define TauAnalysis_Entanglement_printEigenVectors_and_EigenValues_h

#include "TauAnalysis/Entanglement/interface/comp_EigenVectors_and_EigenValues.h" // comp_EigenVectors_and_EigenValues()

#include <TVectorD.h>                                                             // TVectorD

#include <utility>                                                                // std::pair<>

void
printEigenVectors_and_EigenValues(const std::vector<std::pair<TVectorD, double>>& EigenVectors_and_EigenValues);

template <typename T>
void
printEigenVectors_and_EigenValues(const T& symmMatrix, int verbosity = -1)
{
  // CV: matrix passed as function argument needs to be symmetric,
  //     because ROOT can only handle the case that all EigenValues are real
  std::vector<std::pair<TVectorD, double>> EigenVectors_and_EigenValues = comp_EigenVectors_and_EigenValues(symmMatrix, verbosity);
  printEigenVectors_and_EigenValues(EigenVectors_and_EigenValues);
}

#endif // TauAnalysis_Entanglement_printEigenVectors_and_EigenValues_h
