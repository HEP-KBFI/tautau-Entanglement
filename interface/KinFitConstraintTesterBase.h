#ifndef TauAnalysis_Entanglement_KinFitConstraintTesterBase_h
#define TauAnalysis_Entanglement_KinFitConstraintTesterBase_h

#include "TauAnalysis/Entanglement/interface/KinFitConstraintBase.h" // KinFitConstraintBase()

#include <TMatrixD.h>                                                // TMatrixD
#include <TRandom3.h>                                                // TRandom3
#include <TVectorD.h>                                                // TVectorD

#include <string>                                                    // std::string

class KinFitConstraintTesterBase
{
 public:
  KinFitConstraintTesterBase(const TVectorD& alpha0, int verbosity = -1);
  ~KinFitConstraintTesterBase();

  void
  operator()(KinFitConstraintBase& constraint, const std::string& outputFileName);

 protected:
  TVectorD alpha0_;
  
  TRandom3 rnd_;
  TVectorD rndMean_;
  TVectorD rndWidth_;

  TMatrixD plotRange_;

  int verbosity_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraintTesterBase_h
