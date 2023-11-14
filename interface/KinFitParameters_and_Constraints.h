#ifndef TauAnalysis_Entanglement_KinFitParameters_and_Constraints_h
#define TauAnalysis_Entanglement_KinFitParameters_and_Constraints_h

namespace kinFit
{
  // CV: the measured parameters are defined in the following order:
  //       primary vertex position (x,y,z)                                     (3)
  //       Px, Py of neutrino from tau+                                        (2)
  //       decay vertex position (x,y,z) of tau+                               (3)
  //       Px, Py of neutrino from tau-                                        (2)
  //       decay vertex position (x,y,z) of tau-                               (3)
  //       recoil four-vector (Px,Py,Pz,E)                                     (4)
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
  const int numParameters = 17;
  
  //     The constraints are defined in the following order:
  //       "parallelism" constraint for tau+ [1] (2)
  //       "parallelism" constraint for tau- [1] (2) 
  //       constraint that recoil = tau+ + tau-  (4)
  //       Higgs mass constraint                 (1)
  //     Note that the Higgs mass constraintis applied at the LHC,
  //     but not at the SuperKEKB collider (Belle)
  //  [1] cf. Section 4.1.3.3 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
  //
  //     All these constraints are equality constraints.
  //
  //     In addition to these equality constraints, two inequality constraints are used:
  //       mT(visible decay products of tau+, neutrino from tau+ decay) < mTau (1)
  //       mT(visible decay products of tau-, neutrino from tau- decay) < mTau (1)
  // where the symbol mT denotes the transverse mass
  const int numConstraints_LHC = 9;
  const int numConstraints_SuperKEKB = 8;
}

#endif // TauAnalysis_Entanglement_KinFitParameters_and_Constraints_h
