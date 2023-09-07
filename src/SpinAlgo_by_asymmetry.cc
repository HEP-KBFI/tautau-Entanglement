#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_asymmetry.h"

using namespace spin;

SpinAlgo_by_asymmetry::SpinAlgo_by_asymmetry(const edm::ParameterSet& cfg)
  : SpinAlgoBase(cfg)
{}

SpinAlgo_by_asymmetry::~SpinAlgo_by_asymmetry()
{}

namespace
{
  // CV: The class Asymmetry_C computes the asymmetry observable A[B_i], 
  //     given by Eq. (16) in the paper arXiv:2109.09345,
  //     which is used to determine the polarization vector B
  //    (In Eq. (16), the polarization vector B is denoted by the symbol p_i^A)
  class Asymmetry_B
  {
    public:
    Asymmetry_B()
      : N_positive_(0.)
      , N_negative_(0.)
    {}
    ~Asymmetry_B()
    {}

    void
    fill(double cosTheta, double evtWeight)
    {
      if ( cosTheta > 0. ) N_positive_ += evtWeight;
      else                 N_negative_ += evtWeight;
    }

    double
    get_A() const
    {
      double numerator   = N_positive_ - N_negative_;
      double denominator = N_positive_ + N_negative_;
      if ( denominator > 0. )
      {
        return numerator/denominator;
      }
      else
      {
        return 0.;
      }
    }
   private:
    double N_positive_;
    double N_negative_;
  };

  // CV: The class Asymmetry_C computes the asymmetry observable A_ij, 
  //     given by Eq. (25) in the paper arXiv:2205.00542,
  //     which is used to determine the spin correlation matrix C
  class Asymmetry_C
  { 
   public:
    Asymmetry_C()
      : N_positive_(0.)
      , N_negative_(0.)
    {}
    ~Asymmetry_C()
    {}

    void
    fill(double cosTheta_p, double cosTheta_m, double evtWeight)
    {
      double cosTheta_p_times_cosTheta_m = cosTheta_p*cosTheta_m;
      if ( cosTheta_p_times_cosTheta_m > 0. ) N_positive_ += evtWeight;
      else                                    N_negative_ += evtWeight;
    }

    double
    get_A() const
    {
      double numerator   = N_positive_ - N_negative_;
      double denominator = N_positive_ + N_negative_;
      if ( denominator > 0. )
      {
        return numerator/denominator;
      }
      else
      {
        return 0.;
      }
    }
   private:
    double N_positive_;
    double N_negative_;
  };
}

spin::Measurement
SpinAlgo_by_asymmetry::operator()(const spin::Dataset& dataset)
{
  if ( verbosity_ >= 3 )
  {
    std::cout << "<SpinAlgo_by_asymmetry::operator()>:\n";
    std::cout << " #entries = " << dataset.size() << "\n";
  }

  Asymmetry_B Bp_n, Bp_r, Bp_k, Bm_n, Bm_r, Bm_k;
  Asymmetry_C C_nn, C_rn, C_kn, C_nr, C_rr, C_kr, C_nk, C_rk, C_kk;
  size_t numEntries = dataset.size();
  for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
  {
    const spin::Data& entry = dataset.at(idxEntry);

    double hPlus_n = entry.get_hPlus_n();
    double hPlus_r = entry.get_hPlus_r();
    double hPlus_k = entry.get_hPlus_k();
    
    double hMinus_n = entry.get_hMinus_n();
    double hMinus_r = entry.get_hMinus_r();
    double hMinus_k = entry.get_hMinus_k();

    double evtWeight = entry.get_evtWeight();

    if ( verbosity_ >= 3 )
    {
      std::cout << "entry #" << idxEntry << ":\n";
      std::cout << " hPlus: n = " << hPlus_n  << ", r = " << hPlus_r  << ", k = " << hPlus_k  << "\n";
      std::cout << " hMinus: n = " << hMinus_n << ", r = " << hMinus_r << ", k = " << hMinus_k << "\n";
      std::cout << " evtWeight = " << evtWeight << "\n";
    }

    // CV: compute polarization vectors B+ and B- for tau+ and tau- according to Eq. (16)
    //     in the paper arXiv:2109.09345
    Bp_n.fill(hPlus_n, evtWeight);
    Bp_r.fill(hPlus_r, evtWeight);
    Bp_k.fill(hPlus_k, evtWeight);

    Bm_n.fill(hMinus_n, evtWeight);
    Bm_r.fill(hMinus_r, evtWeight);
    Bm_k.fill(hMinus_k, evtWeight);

    // CV: compute spin correlation matrix C according to Eq. (25)
    //     in the paper arXiv:2205.00542
    C_nn.fill(hPlus_n, hMinus_n, evtWeight);
    C_rn.fill(hPlus_r, hMinus_n, evtWeight);
    C_kn.fill(hPlus_k, hMinus_n, evtWeight);
    C_nr.fill(hPlus_n, hMinus_r, evtWeight);
    C_rr.fill(hPlus_r, hMinus_r, evtWeight);
    C_kr.fill(hPlus_k, hMinus_r, evtWeight);
    C_nk.fill(hPlus_n, hMinus_k, evtWeight);
    C_rk.fill(hPlus_r, hMinus_k, evtWeight);
    C_kk.fill(hPlus_k, hMinus_k, evtWeight);
  }

  // CV: factors +2 for tau+ and -2 for tau- taken from Eq. (16) in the paper arXiv:2109.09345
  //    (the spin analyzing power alpha is +1 for tau+ and -1 for tau-,
  //     cf. text following Eq. (24) in the paper arXiv:2205.00542)
  const double bp = +2.;
  math::Vector3 Bp;
  Bp(0) = bp*Bp_n.get_A();
  Bp(1) = bp*Bp_r.get_A();
  Bp(2) = bp*Bp_k.get_A();

  const double bm = -2.;
  math::Vector3 Bm;
  Bm(0) = bm*Bm_n.get_A();
  Bm(1) = bm*Bm_r.get_A();
  Bm(2) = bm*Bm_k.get_A();

  // CV: factor -4 taken from Eq. (25) in the paper arXiv:2205.00542
  //    (the spin analyzing power alpha is +1 for tau+ and -1 for tau-,
  //     cf. text following Eq. (24) in the same paper)
  const double c = -4.;
  math::Matrix3x3 C;
  C(0,0) = c*C_nn.get_A();
  C(0,1) = c*C_rn.get_A();
  C(0,2) = c*C_kn.get_A();
  C(1,0) = c*C_nr.get_A();
  C(1,1) = c*C_rr.get_A();
  C(1,2) = c*C_kr.get_A();
  C(2,0) = c*C_nk.get_A();
  C(2,1) = c*C_rk.get_A();
  C(2,2) = c*C_kk.get_A();

  spin::Measurement measurement(Bp, Bm, C);
  addEntanglementVariables(measurement);
  return measurement;
}
