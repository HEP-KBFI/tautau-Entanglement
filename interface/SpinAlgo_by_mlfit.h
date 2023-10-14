#ifndef TauAnalysis_Entanglement_SpinAlgo_by_mlfit_h
#define TauAnalysis_Entanglement_SpinAlgo_by_mlfit_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"               // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/Dataset.h"               // spin::Dataset
#include "TauAnalysis/Entanglement/interface/Measurement.h"           // spin::Measurement
#include "TauAnalysis/Entanglement/interface/SpinAlgoBase.h"          // SpinAlgoBase
#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_summation.h" // SpinAlgo_by_summation

#include <Minuit2/Minuit2Minimizer.h>                                 // ROOT::Minuit2::Minuit2Minimizer

namespace spin
{

class SpinAlgo_by_mlfit : public SpinAlgoBase
{
 public:
  SpinAlgo_by_mlfit(const edm::ParameterSet& cfg);
  ~SpinAlgo_by_mlfit();

  void
  set_verbosity(int verbosity);

  static
  void
  set_dataset_norm_passed(const spin::Dataset* dataset);

  static
  void
  set_dataset_norm_failed(const spin::Dataset* dataset);

  spin::Measurement
  operator()(const spin::DatasetPtrs& dataset);

 private:
  ROOT::Math::Minimizer* mlfit_;
  std::vector<double> par_gen_;
  bool scanLikelihood_;
  std::string outputFileName_;

  SpinAlgo_by_summation* algo_by_summation_;
};

}

#endif // TauAnalysis_Entanglement_SpinAlgo_by_mlfit_h
