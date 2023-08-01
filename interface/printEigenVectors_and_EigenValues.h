#ifndef TauAnalysis_Entanglement_printEigenVectors_and_EigenValues_h
#define TauAnalysis_Entanglement_printEigenVectors_and_EigenValues_h

#include "FWCore/Utilities/interface/Exception.h" // edm::Exception

#include <TMatrixD.h>                             // TMatrixD
#include <TVectorD.h>                             // TVectorD

template <typename T>
void
printEigenVectors_and_EigenValues(const T& cov)
{
  assert(T::kRows == T::kCols);
  int dim =  T::kRows;

  TMatrixD tmp(dim,dim);
  for ( int idxRow = 0; idxRow < dim; ++idxRow )
  {
    for ( int idxColumn = 0; idxColumn < dim; ++idxColumn )
    {
      tmp(idxRow,idxColumn) = cov(idxRow,idxColumn);
    }
  }

  TVectorD eigenValues(dim);
  TMatrixD eigenVectors(dim,dim);
  try
  { 
    eigenVectors = tmp.EigenVectors(eigenValues);
  }
  catch ( edm::Exception )
  {
    std::cerr << "Error in <printEigenVectors_and_EigenValues>: Failed to compute Eigenvectors and Eigenvalues !!\n";
    return;
  }
  for ( int idxEigenVector = 0; idxEigenVector < dim; ++idxEigenVector )
  {
    TVectorD eigenVector(dim);
    for ( int idxComponent = 0; idxComponent < dim; ++idxComponent )
    {
      eigenVector(idxComponent) = eigenVectors(idxComponent,idxEigenVector);
    }
    std::cout << "EigenVector #" << idxEigenVector << " (EigenValue = " << eigenValues(idxEigenVector) << "):\n";
    eigenVector.Print();
  }
}

#endif // TauAnalysis_Entanglement_printEigenVectors_and_EigenValues_h
