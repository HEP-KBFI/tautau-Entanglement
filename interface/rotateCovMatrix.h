#ifndef TauAnalysis_Entanglement_rotateCovMatrix_h
#define TauAnalysis_Entanglement_rotateCovMatrix_h

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3

math::Matrix3x3
rotateCovMatrix(const math::Matrix3x3& cov, const math::Matrix3x3& rotMatrix);

#endif // TauAnalysis_Entanglement_rotateCovMatrix_h
