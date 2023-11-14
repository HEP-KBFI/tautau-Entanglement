#include "TauAnalysis/Entanglement/interface/printMatrix.h"

void
printMatrix(const std::string& label, const TMatrixD& m)
{
  std::cout << label << ":\n";
  m.Print();
}
