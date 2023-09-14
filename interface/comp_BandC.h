#ifndef TauAnalysis_Entanglement_comp_BandC_h
#define TauAnalysis_Entanglement_comp_BandC_h

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3, math::Vector3

math::Vector3
comp_Bp(double hPlus_n, double hPlus_r, double hPlus_k);

math::Vector3
comp_Bm(double hMinus_n, double hMinus_r, double hMinus_k);

math::Matrix3x3
comp_C(double hPlus_n, double hPlus_r, double hPlus_k,
       double hMinus_n, double hMinus_r, double hMinus_k);

#endif // TauAnalysis_Entanglement_comp_BandC_h
