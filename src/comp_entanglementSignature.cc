#include "TauAnalysis/Entanglement/interface/comp_entanglementSignature.h"

double
comp_entanglementSignature(const math::Matrix3x3& C)
{
  const double trace = C.Trace();
  double E = 0.;
  for(int iRow = 0; iRow < 3; ++iRow)
  {
    const double E_trial = std::fabs(trace - C(iRow, iRow)) - C(iRow, iRow);
    if(E_trial > E)
    {
      E = E_trial;
    }
  }
  return E;
}
