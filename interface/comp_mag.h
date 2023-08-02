#ifndef TauAnalysis_Entanglement_comp_mag_h
#define TauAnalysis_Entanglement_comp_mag_h

#include "TauAnalysis/Entanglement/interface/square.h" // square()

template <typename T>
double
comp_mag(const T& v)
{
  double mag2 = 0.;
  for ( int idx = 0; idx < T::kSize; ++idx )
  {
    mag2 += square(v(idx));
  }
  return std::sqrt(mag2);
}

#endif // TauAnalysis_Entanglement_comp_mag_h
