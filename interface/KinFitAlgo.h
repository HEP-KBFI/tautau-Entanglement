#ifndef TauAnalysis_Entanglement_KinFitAlgo_h
#define TauAnalysis_Entanglement_KinFitAlgo_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"              // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/cmsException.h"         // cmsException
#include "TauAnalysis/Entanglement/interface/KinFitConstraintBase.h" // KinFitConstraintBase

#include <TMatrixD.h>                                                // TMatrixD
#include <TVectorD.h>                                                // TVectorD

class KinFitSummary
{
 public:
  KinFitSummary();
  KinFitSummary(int iteration, const TVectorD& alpha, const TMatrixD& V_alpha, double chi2, int status);
  ~KinFitSummary();

  KinFitSummary&
  operator=(const KinFitSummary& fit);

  int
  get_iteration() const;

  const TVectorD&
  get_alpha() const;

  const TMatrixD&
  get_V_alpha() const;

  double
  get_chi2() const;

  int
  get_status() const;

 protected:
  int iteration_;
  TVectorD alpha_;
  TMatrixD V_alpha_;
  double chi2_;
  int status_;
};

class KinFitAlgo
{
 public:
  KinFitAlgo(const edm::ParameterSet& cfg);
  ~KinFitAlgo();

  KinFitSummary
  operator()(const TVectorD& alpha0, const TMatrixD& V_alpha0, KinFitConstraintBase& constraint, const TVectorD& startPos,
             std::vector<KinFitSummary>* fitHistory = nullptr);

 private:
  int verbosity_;
};

#endif // TauAnalysis_Entanglement_KinFitAlgo_h
