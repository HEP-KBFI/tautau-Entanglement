#ifndef TauAnalysis_Entanglement_Matrix_and_Vector_h
#define TauAnalysis_Entanglement_Matrix_and_Vector_h

#include "DataFormats/Math/interface/Matrix.h" // math::Matrix
#include "DataFormats/Math/interface/Vector.h" // math::Vector

namespace kinFit
{
  const int numParameters = 19;
  // CV: the measured parameters are defined in the following order:
  //       primary vertex                        (3)
  //       momentum of neutrino from tau+        (3)
  //       decay vertex of tau+                  (3)
  //       momentum of neutrino from tau-        (3)
  //       decay vertex of tau-                  (3)
  //       recoil four-vector                    (4)
  // where:
  //     the energy and momentum components of four-vectors are given in the order:
  //      (px, py, pz, E)
  //     the position of vertices are given in the order:
  //      (x, y, z)
  //     cf. Section II of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
  //
  //     The four-vectors of tau+ and tau- are not really "measured";
  //     we use the "huge error method" described in Section 6 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
  //     and set their covariance matrix to diagonal matrix with large values on the diagonal

  const int numConstraints  = 11;
  // CV: constrains are defined in the following order:
  //       Higgs mass constraint                 (1)
  //       neutrino mass constraint for tau+     (1)
  //       "parallelism" constraint for tau+ [1] (2)
  //       neutrino mass constraint for tau-     (1)
  //       "parallelism" constraint for tau- [1] (2) 
  //       constraint that recoil = tau+ + tau-  (4)
  //  [1] cf. Section 4.1.3.3 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
}

namespace math
{
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

  const int P = kinFit::numParameters;
  const int C = kinFit::numConstraints;

  typedef Matrix<P,P>::type MatrixPxP;
  typedef Matrix<C,P>::type MatrixCxP;
  typedef Matrix<P,C>::type MatrixPxC;
  typedef Matrix<C,C>::type MatrixCxC;
  typedef Vector<P>::type   VectorP;
  typedef Vector<C>::type   VectorC;
}

#endif // TauAnalysis_Entanglement_Matrix_and_Vector_h
