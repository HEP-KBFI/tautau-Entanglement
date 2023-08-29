#ifndef TauAnalysis_Entanglement_comp_steerability_h
#define TauAnalysis_Entanglement_comp_steerability_h

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3

/**
 * @brief Computes steerability according to Eq. 21 of the paper arXiv:2211.10513
 * @param C Spin correlation matrix
 * @param n_theta Number of samples in polar angle, default 200
 * @param n_phi Number of samples in azimuthal angle, default 200
 * @return Steerability S
 *
 * @todo Is this still valid in case the polarization vectors are nonzero?
 * @todo Return MC integration error as well?
 * @todo Kahan summation?
 */
double
comp_steerability(const math::Matrix3x3& C,
                  unsigned n_theta = 200,
                  unsigned n_phi = 200);

#endif // TauAnalysis_Entanglement_comp_steerability_h
