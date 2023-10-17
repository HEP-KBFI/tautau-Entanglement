#ifndef TauAnalysis_Entanglement_rotateVector_h
#define TauAnalysis_Entanglement_rotateVector_h

#include "DataFormats/Candidate/interface/Candidate.h"            // reco::Candidate::Vector

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3, math::Vector3

template <typename T>
T
rotateVector(const T& v, const math::Matrix3x3& rotMatrix)
{
  // compute elements of vector v in rotated coordinate system, given by the matrix rotMatrix,
  // as described in Example 4.25 in Section 4.9.2 of
  //   V. Blobel and E. Lohrmann "Statistische und numerische Methoden der Datenanalyse".
  math::Vector3 tmp;
  tmp(0) = v.x();
  tmp(1) = v.y();
  tmp(2) = v.z();
  tmp = rotMatrix*tmp;  
  return T(tmp(0), tmp(1), tmp(2));
}

#endif // TauAnalysis_Entanglement_rotateVector_h
