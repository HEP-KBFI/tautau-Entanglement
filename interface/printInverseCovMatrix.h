#ifndef TauAnalysis_Entanglement_printInverseCovMatrix_h
#define TauAnalysis_Entanglement_printInverseCovMatrix_h

#include "TauAnalysis/Entanglement/interface/printCovMatrix.h" // printCovMatrix()

#include <TString.h> // Form()

template <typename T>
void
printInverseCovMatrix(const std::string& label, const T& cov)
{
  int errorFlag = 0;
  T covInv = cov.Inverse(errorFlag);
  if ( errorFlag != 0 )
  {
    std::cerr << "Erorr in <printCovMatrix>: Failed to invert matrix cov !!\n";
    return;
  }
  std::cout << label << "^-1:\n";
  std::cout << covInv << "\n";
  double det = -1.;
  covInv.Det2(det);
  std::cout << " det = " << det << "\n";
}

#endif // TauAnalysis_Entanglement_printInverseCovMatrix_h
