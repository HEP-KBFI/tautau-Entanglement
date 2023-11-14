#ifndef TauAnalysis_Entanglement_convert_to_TMatrixD_h
#define TauAnalysis_Entanglement_convert_to_TMatrixD_h

#include <TMatrixD.h> // TMatrixD

template <typename T>
TMatrixD
convert_to_TMatrixD(const T& matrix, int verbosity = -1)
{
  int numRows = T::kRows;
  int numColumns = T::kCols;
  TMatrixD retVal(numRows,numColumns);
  for ( int idxRow = 0; idxRow < numRows; ++idxRow )
  {
    for ( int idxColumn = 0; idxColumn < numColumns; ++idxColumn )
    {
      retVal(idxRow,idxColumn) = matrix(idxRow,idxColumn);
    }
  }
  return retVal;
}

#endif // TauAnalysis_Entanglement_convert_to_TMatrixD_h
