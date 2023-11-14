#ifndef TauAnalysis_Entanglement_convert_to_SMatrix_h
#define TauAnalysis_Entanglement_convert_to_SMatrix_h

#include "DataFormats/Math/interface/Matrix.h"               // math::Matrix

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException

#include <TMatrixD.h>                                        // TMatrixD

template <unsigned int numRows, unsigned int numColumns>
typename math::Matrix<numRows,numColumns>::type
convert_to_SMatrix(const TMatrixD& matrix, int verbosity = -1)
{
  if ( !(matrix.GetNrows() == numRows && matrix.GetNcols() == numColumns) )
    throw cmsException("KinematicFit", __LINE__)
      << "Dimension of matrix given as function argument (" << matrix.GetNrows() << "," << matrix.GetNcols() << ")" 
      << " does not match template parameters (" << numRows << "," << numColumns << ") !!\n";
  typename math::Matrix<numRows,numColumns>::type retVal;
  for ( unsigned int idxRow = 0; idxRow < numRows; ++idxRow )
  {
    for ( unsigned int idxColumn = 0; idxColumn < numColumns; ++idxColumn )
    {
      retVal(idxRow,idxColumn) = matrix(idxRow,idxColumn);
    }
  }
  return retVal;
}

#endif // TauAnalysis_Entanglement_convert_to_SMatrix_h
