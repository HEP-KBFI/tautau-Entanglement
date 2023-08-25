#ifndef TauAnalysis_Entanglement_comp_EigenVectors_and_EigenValues_h
#define TauAnalysis_Entanglement_comp_EigenVectors_and_EigenValues_h

#include "DataFormats/Math/interface/Matrix.h"       // math::Matrix
#include "DataFormats/Math/interface/Vector.h"       // math::Vector
#include "FWCore/Utilities/interface/EDMException.h" // edm::Exception

#include <TMatrixD.h>                                // TMatrixD
#include <TVectorD.h>                                // TVectorD

#include <utility>                                   // std::make_pair(), std::pair<>

template <typename T>
std::vector<std::pair<TVectorD, double>>
comp_EigenVectors_and_EigenValues(const T& symmMatrix, int verbosity = -1)
{
  if ( verbosity >= 3 )
  {
    std::cout << "<comp_EigenVectors_and_EigenValues>:\n";
  }

  // CV: matrix passed as function argument needs to be symmetric,
  //     because ROOT can only handle the case that all EigenValues are real
  assert(T::kRows == T::kCols);
  int dim = T::kRows;

  std::vector<std::pair<TVectorD, double>> EigenVectors_and_EigenValues;

  TMatrixD tmp(dim,dim);
  for ( int idxRow = 0; idxRow < dim; ++idxRow )
  {
    for ( int idxColumn = 0; idxColumn < dim; ++idxColumn )
    {
      tmp(idxRow,idxColumn) = symmMatrix(idxRow,idxColumn);
    }
  }
  if ( verbosity >= 3 )
  {
    std::cout << "tmp:\n";
    tmp.Print();
  }

  TVectorD EigenValues(dim);
  TMatrixD EigenVectors(dim,dim);
  try
  { 
    EigenVectors = tmp.EigenVectors(EigenValues);
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

#endif // TauAnalysis_Entanglement_comp_EigenVectors_and_EigenValues_h
