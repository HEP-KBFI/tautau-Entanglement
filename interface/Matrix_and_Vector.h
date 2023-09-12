#ifndef TauAnalysis_Entanglement_Matrix_and_Vector_h
#define TauAnalysis_Entanglement_Matrix_and_Vector_h

#include "DataFormats/Math/interface/Matrix.h" // math::Matrix
#include "DataFormats/Math/interface/Vector.h" // math::Vector

namespace math
{
  typedef Matrix<2,2>::type Matrix2x2;
  typedef Vector<2>::type   Vector2;

  typedef Matrix<3,3>::type Matrix3x3;
  typedef Vector<3>::type   Vector3;

  typedef Matrix<4,4>::type Matrix4x4;
  typedef Vector<4>::type   Vector4;

  typedef Matrix<5,5>::type Matrix5x5;
  typedef Vector<5>::type   Vector5;

  typedef Matrix<5,7>::type Matrix5x7;
  typedef Matrix<7,5>::type Matrix7x5;

  typedef Matrix<7,7>::type Matrix7x7;
  typedef Vector<7>::type   Vector7;
}

#endif // TauAnalysis_Entanglement_Matrix_and_Vector_h
