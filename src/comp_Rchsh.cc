#include "TauAnalysis/Entanglement/interface/comp_Rchsh.h"

#include "FWCore/Utilities/interface/Exception.h"                                 // cms::Exception

#include "TauAnalysis/Entanglement/interface/comp_EigenVectors_and_EigenValues.h" // comp_EigenVectors_and_EigenValues()
#include "TauAnalysis/Entanglement/interface/printEigenVectors_and_EigenValues.h" // printEigenVectors_and_EigenValues()

#include "Math/Functions.h"                                                       // ROOT::Math::Transpose() 

double
comp_Rchsh(const math::Matrix3x3& C, int verbosity)
{
  // CV: compute observable Rchsh according to Eq. (10)
  //     in the paper arXiv:2211.10513
  // Karl (18/10/23): we quote its square in the paper, hence no square root operation in [*]
  if ( verbosity >= 3 )
  {
    std::cout << "<comp_Rchsh>:\n";
    std::cout << "C:\n";
    std::cout << C << "\n";
  }
  math::Matrix3x3 CT = ROOT::Math::Transpose(C);
  math::Matrix3x3 CT_times_C = CT*C;
  if ( verbosity >= 3 )
  {
    std::cout << "C^T*C:\n";
    std::cout << CT_times_C << "\n";
  }
  std::vector<std::pair<TVectorD, double>> EigenVectors_and_EigenValues;
  try
  { 
    EigenVectors_and_EigenValues = comp_EigenVectors_and_EigenValues(CT_times_C, verbosity);
  }
  catch ( const cms::Exception& )
  {
    std::cerr << "Error in <comp_Rchsh>: Caught exception from comp_EigenVectors_and_EigenValues !!\n";
    return -1.;
  }
  if ( verbosity >= 3 )
  {
    printEigenVectors_and_EigenValues(EigenVectors_and_EigenValues);
  }
  assert(EigenVectors_and_EigenValues.size() == 3);
  double Rchsh = EigenVectors_and_EigenValues[0].second + EigenVectors_and_EigenValues[1].second; // [*]
  if ( verbosity >= 3 )
  {
    std::cout << "Rchsh = " << Rchsh << "\n";
  }
  return Rchsh;
}
