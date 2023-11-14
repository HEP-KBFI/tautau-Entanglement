#ifndef TauAnalysis_Entanglement_printInverseMatrix_h
#define TauAnalysis_Entanglement_printInverseMatrix_h

#include "TauAnalysis/Entanglement/interface/printMatrix.h" // printMatrix()

#include <TMatrixD.h>                                       // TMatrixD
#include <TString.h>                                        // Form()

template <typename T>
void
printInverseMatrix(const std::string& label, const T& m, bool printDet = false)
{
  int errorFlag = 0;
  T mInv = cov.Inverse(errorFlag);
  if ( errorFlag != 0 )
  {
    std::cerr << "Erorr in <printCovMatrix>: Failed to invert matrix m !!\n";
    return;
  }
  printMatrix(Form("%s^-1", label.c_str()), mInv, printDet);
}

void
printInverseMatrix(const std::string& label, const TMatrixD& m)
{
  if ( m.Determinant() == 0. )
  {
    std::cerr << "Erorr in <printCovMatrix>: Failed to invert matrix m !!\n";
    return;
  }
  TMatrixD mInv(TMatrixD::kInverted, m);
  printMatrix(Form("%s^-1", label.c_str()), mInv);
}

#endif // TauAnalysis_Entanglement_printInverseMatrix_h
