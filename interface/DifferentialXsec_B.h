#ifndef TauAnalysis_Entanglement_DifferentialXsec_B_h
#define TauAnalysis_Entanglement_DifferentialXsec_B_h

#include <TH1.h>  // TH1

#include <string> // std::string

namespace spin
{

// CV: The class DifferentialXsec_B represents the differential cross section
//       dsigma/dcosTheta
//     given by Eq. (4.18) in the paper arXiv:1508.05271
class DifferentialXsec_B
{
 public:
  DifferentialXsec_B(const std::string& label, int numBins = 40);
  ~DifferentialXsec_B();

  void
  fill(double cosTheta, double evtWeight);

  const std::string&
  get_label() const;

  const TH1*
  get_histogram() const;

 private:
  std::string label_;
  TH1* histogram_;
};

double
fit_DifferentialXsec_B(const DifferentialXsec_B& B_i, int verbosity);

}

#endif // TauAnalysis_Entanglement_DifferentialXsec_B_h
