#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_summation.h"

#include "TauAnalysis/Entanglement/interface/comp_BandC.h" // comp_Bm(), comp_Bp(), comp_C()

using namespace spin;

SpinAlgo_by_summation::SpinAlgo_by_summation(const edm::ParameterSet& cfg)
  : SpinAlgoBase(cfg)
{}

SpinAlgo_by_summation::~SpinAlgo_by_summation()
{}

spin::Measurement
SpinAlgo_by_summation::operator()(const spin::Dataset& dataset)
{
  if ( verbosity_ >= 3 )
  {
    std::cout << "<SpinAlgo_by_summation::operator()>:\n";
    std::cout << " #entries = " << dataset.size() << "\n";
  }

  math::Vector3 Bp_sum, Bm_sum;
  math::Matrix3x3 C_sum;
  double evtWeight_sum = 0.;
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

    // CV: compute polarization vectors B+ and B- for tau+ and tau- according to text following Eq. (4.18)
    //     in the paper arXiv:1508.05271
    math::Vector3 Bp = comp_Bp(hPlus_n, hPlus_r, hPlus_k);
    Bp_sum += evtWeight*Bp;

    math::Vector3 Bm = comp_Bm(hMinus_n, hMinus_r, hMinus_k);
    Bm_sum += evtWeight*Bm;
    
    // CV: compute spin correlation matrix C according to Eq. (25)
    //     in the paper arXiv:2211.10513
    math::Matrix3x3 C = comp_C(hPlus_n, hPlus_r, hPlus_k, hMinus_n, hMinus_r, hMinus_k);
    C_sum += evtWeight*C;

    evtWeight_sum += evtWeight;
  }

  Bp_sum *= (1./evtWeight_sum);
  Bm_sum *= (1./evtWeight_sum);

  C_sum *= (1./evtWeight_sum);

  spin::Measurement measurement(Bp_sum, Bm_sum, C_sum);
  addEntanglementVariables(measurement);
  return measurement;
}
