#ifndef TauAnalysis_Entanglement_KinFitParameters_h
#define TauAnalysis_Entanglement_KinFitParameters_h

#include "DataFormats/Math/interface/Matrix.h" // math::Matrix
#include "DataFormats/Math/interface/Vector.h" // math::Vector

namespace kinFit
{
  const int numParameters = 17;
  // CV: the measured parameters are defined in the following order:
  //       primary vertex position (x,y,z)       (3)
  //       Px, Py of neutrino from tau+          (2)
  //       decay vertex position (x,y,z) of tau+ (3)
  //       Px, Py of neutrino from tau-          (2)
  //       decay vertex position (x,y,z) of tau- (3)
  //       recoil four-vector (Px,Py,Pz,E)       (4)
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
}

namespace math
{
  const int P = kinFit::numParameters;
  typedef Matrix<P,P>::type MatrixPxP;
  typedef Vector<P>::type   VectorP;
}

#endif // TauAnalysis_Entanglement_KinFitParameters_h
