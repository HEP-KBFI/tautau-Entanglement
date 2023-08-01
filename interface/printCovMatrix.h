#ifndef TauAnalysis_Entanglement_printCovMatrix_h
#define TauAnalysis_Entanglement_printCovMatrix_h

template <typename T>
void
printCovMatrix(const std::string& label, const T& cov)
{
  std::cout << label << ":\n";
  std::cout << cov << "\n";
  double det = -1.;
  cov.Det2(det);
  std::cout << " det = " << det << "\n";
}

#endif // TauAnalysis_Entanglement_printCovMatrix_h
