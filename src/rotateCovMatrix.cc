#include "TauAnalysis/Entanglement/interface/rotateCovMatrix.h"

math::Matrix3x3
rotateCovMatrix(const math::Matrix3x3& cov, const math::Matrix3x3& rotMatrix)
{
  // compute elements of covariance matrix cov in rotated coordinate system, given by the matrix rotMatrix,
  // as described in Example 4.25 in Section 4.9.2 of
  //   V. Blobel and E. Lohrmann "Statistische und numerische Methoden der Datenanalyse".
  math::Matrix3x3 rotMatrixT = ROOT::Math::Transpose(rotMatrix);
  return rotMatrix*cov*rotMatrixT;
}
