#ifndef TauAnalysis_Entanglement_printMatrix_h
#define TauAnalysis_Entanglement_printMatrix_h

#include <TMatrixD.h> // TMatrixD

#include <iostream>   // std::cout

template <typename T>
void
printMatrix(const std::string& label, const T& m, bool printDet = false)
{
  std::cout << label << ":\n";
  std::cout << m << "\n";
  if ( printDet )
  {
    double det = -1.;
    m.Det2(det);
    std::cout << " det = " << det << "\n";
  }
}

void
printMatrix(const std::string& label, const TMatrixD& m);

#endif // TauAnalysis_Entanglement_printMatrix_h
