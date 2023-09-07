#include "TauAnalysis/Entanglement/interface/comp_BandC.h"

namespace
{
  math::Vector3
  comp_B(double h_n, double h_r, double h_k,
         double b, double evtWeight)
  {
    // CV: compute polarization vectors B+ and B- for tau+ and tau- according to text following Eq. (4.18)
    //     in the paper arXiv:1508.05271
    math::Vector3 B;
    B(0) = b*evtWeight*h_n;
    B(1) = b*evtWeight*h_r;
    B(2) = b*evtWeight*h_k;
    return B;
  }
}

math::Vector3
comp_Bp(double hPlus_n, double hPlus_r, double hPlus_k,
        double evtWeight)
{
  // CV: not sure if for tau+ the constant factor b should be +3 or -3 ?!
  double b = 3.;
  return comp_B(hPlus_n, hPlus_r, hPlus_k, b, evtWeight);
}

math::Vector3
comp_Bm(double hMinus_n, double hMinus_r, double hMinus_k,
        double evtWeight)
{
  // CV: not sure if for tau- the constant factor b should be +3 or -3 ?!
  double b = 3.;
  return comp_B(hMinus_n, hMinus_r, hMinus_k, b, evtWeight);
}

math::Matrix3x3
comp_C(double hPlus_n, double hPlus_r, double hPlus_k,
       double hMinus_n, double hMinus_r, double hMinus_k,
       double evtWeight)
{
  // CV: compute spin correlation matrix C according to Eq. (25)
  //     in the paper arXiv:2211.10513.
  //     The ordering of rows vs columns for tau+ and tau- has been agreed with Luca on 06/09/2023.
  math::Matrix3x3 C;
  double c = -9.;
  C(0,0) = c*evtWeight*hPlus_n*hMinus_n;
  C(0,1) = c*evtWeight*hPlus_r*hMinus_n;
  C(0,2) = c*evtWeight*hPlus_k*hMinus_n;
  C(1,0) = c*evtWeight*hPlus_n*hMinus_r;
  C(1,1) = c*evtWeight*hPlus_r*hMinus_r;
  C(1,2) = c*evtWeight*hPlus_k*hMinus_r;
  C(2,0) = c*evtWeight*hPlus_n*hMinus_k;
  C(2,1) = c*evtWeight*hPlus_r*hMinus_k;
  C(2,2) = c*evtWeight*hPlus_k*hMinus_k;
  return C;
}
