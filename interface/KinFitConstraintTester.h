#ifndef TauAnalysis_Entanglement_KinFitConstraintTester_h
#define TauAnalysis_Entanglement_KinFitConstraintTester_h

#include "TauAnalysis/Entanglement/interface/KinFitConstraintTesterBase.h" // KinFitConstraintTesterBase

class KinFitConstraintTester : public KinFitConstraintTesterBase
{
 public:
  KinFitConstraintTester(const TVectorD& alpha0, int verbosity = -1);

  void
  operator()(KinFitConstraintBase& constraint, const std::string& outputFileName);
};

#endif // TauAnalysis_Entanglement_KinFitConstraintTester_h
