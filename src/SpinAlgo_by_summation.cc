#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_summation.h"

#include "TauAnalysis/Entanglement/interface/comp_Rchsh.h" // comp_Rchsh()

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

  math::Vector3 Bp, Bm;
  math::Matrix3x3 C;
  double evtWeight_sum = 0.;
  size_t numEntries = dataset.size();
  for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
  {
    const spin::Data& entry = dataset.at(idxEntry);

    double hPlus_r = entry.get_hPlus_r();
    double hPlus_n = entry.get_hPlus_n();
    double hPlus_k = entry.get_hPlus_k();
    
    double hMinus_r = entry.get_hMinus_r();
    double hMinus_n = entry.get_hMinus_n();
    double hMinus_k = entry.get_hMinus_k();

    double evtWeight = entry.get_evtWeight();

    if ( verbosity_ >= 3 )
    {
      std::cout << "entry #" << idxEntry << ":\n";
      std::cout << " hPlus: r = " << hPlus_r  << ", n = " << hPlus_n  << ", k = " << hPlus_k  << "\n";
      std::cout << " hMinus: r = " << hMinus_r << ", n = " << hMinus_n << ", k = " << hMinus_k << "\n";
      std::cout << " evtWeight = " << evtWeight << "\n";
    }

    Bp(0) += evtWeight*hPlus_r;
    Bp(1) += evtWeight*hPlus_n;
    Bp(2) += evtWeight*hPlus_k;

    Bm(0) += evtWeight*hMinus_r;
    Bm(1) += evtWeight*hMinus_n;
    Bm(2) += evtWeight*hMinus_k;
    
    // CV: compute matrix C according to Eq. (25)
    //     in the paper arXiv:2211.10513
    double c = -9.*evtWeight;
    C(0,0) += c*hPlus_r*hMinus_r;
    C(0,1) += c*hPlus_r*hMinus_n;
    C(0,2) += c*hPlus_r*hMinus_k;
    C(1,0) += c*hPlus_n*hMinus_r;
    C(1,1) += c*hPlus_n*hMinus_n;
    C(1,2) += c*hPlus_n*hMinus_k;
    C(2,0) += c*hPlus_k*hMinus_r;
    C(2,1) += c*hPlus_k*hMinus_n;
    C(2,2) += c*hPlus_k*hMinus_k;

    evtWeight_sum += evtWeight;
  }

  Bp *= (1./evtWeight_sum);
  Bm *= (1./evtWeight_sum);

  C *= (1./evtWeight_sum);

  double Rchsh = comp_Rchsh(C, verbosity_);

  spin::Measurement measurement(Bp, Bm, C, Rchsh);
  return measurement;
}
