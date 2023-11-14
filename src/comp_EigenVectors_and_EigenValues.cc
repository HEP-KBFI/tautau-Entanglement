#include "TauAnalysis/Entanglement/interface/comp_EigenVectors_and_EigenValues.h"

#include "FWCore/Utilities/interface/EDMException.h" // edm::Exception

#include <assert.h>                                  // assert()
#include <iostream>                                  // std::cout

std::vector<std::pair<TVectorD, double>>
comp_EigenVectors_and_EigenValues(const TMatrixD& symmMatrix, int verbosity)
{
  if ( verbosity >= 3 )
  {
    std::cout << "<comp_EigenVectors_and_EigenValues>:\n";
    std::cout << "symmMatrix:\n";
    symmMatrix.Print();
  }

  // CV: matrix passed as function argument needs to be symmetric,
  //     because ROOT can only handle the case that all EigenValues are real
  assert(symmMatrix.GetNcols() == symmMatrix.GetNrows());
  int dim = symmMatrix.GetNrows();

  std::vector<std::pair<TVectorD, double>> EigenVectors_and_EigenValues;

  TVectorD EigenValues(dim);
  TMatrixD EigenVectors(dim,dim);
  try
  { 
    EigenVectors = symmMatrix.EigenVectors(EigenValues);
  }
  catch ( const edm::Exception& )
  {
    std::cerr << "Error in <comp_EigenVectors_and_EigenValues>: Failed to compute Eigenvectors and Eigenvalues !!\n";
    return std::vector<std::pair<TVectorD, double>>();
  }
  for ( int idxEigenVector = 0; idxEigenVector < dim; ++idxEigenVector )
  {
    TVectorD EigenVector(dim);
    for ( int idxComponent = 0; idxComponent < dim; ++idxComponent )
    {
      EigenVector(idxComponent) = EigenVectors(idxComponent,idxEigenVector);
    }
    double EigenValue = EigenValues(idxEigenVector);
    EigenVectors_and_EigenValues.push_back(std::make_pair(EigenVector, EigenValue));
  }

  return EigenVectors_and_EigenValues;
}
