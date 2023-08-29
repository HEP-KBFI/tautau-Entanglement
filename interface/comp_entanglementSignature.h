#ifndef TauAnalysis_Entanglement_comp_entanglementSignature_h
#define TauAnalysis_Entanglement_comp_entanglementSignature_h

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3

/**
 * @brief Computes entanglement signature according to Eq. 4 of the paper arXiv:2211.10513
 * @param C Spin correlation matrix
 * @return Entanglement signature E_k
 */
double
comp_entanglementSignature(const math::Matrix3x3& C);

#endif // TauAnalysis_Entanglement_comp_entanglementSignature_h
