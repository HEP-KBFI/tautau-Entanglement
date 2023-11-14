#include "TauAnalysis/Entanglement/interface/printEigenVectors_and_EigenValues.h"

#include <iostream> // std::cout

void
printEigenVectors_and_EigenValues(const std::vector<std::pair<TVectorD, double>>& EigenVectors_and_EigenValues)
{
  size_t dim = EigenVectors_and_EigenValues.size();
  for ( size_t idxEigenVector = 0; idxEigenVector < dim; ++idxEigenVector )
  {
    const std::pair<TVectorD, double>& EigenVector_and_EigenValue = EigenVectors_and_EigenValues.at(idxEigenVector);
    const TVectorD& EigenVector = EigenVector_and_EigenValue.first;
    double EigenValue = EigenVector_and_EigenValue.second;
    std::cout << "EigenVector #" << idxEigenVector << " (EigenValue = " << EigenValue << "):\n";
    EigenVector.Print();
  }
}
