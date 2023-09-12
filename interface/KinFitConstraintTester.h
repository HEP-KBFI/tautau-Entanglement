#ifndef TauAnalysis_Entanglement_KinFitConstraintTester_h
#define TauAnalysis_Entanglement_KinFitConstraintTester_h

#include "TauAnalysis/Entanglement/interface/KinFitConstraintTesterBase.h" // KinFitConstraintTesterBase

template <unsigned int P, unsigned int C>
class KinFitConstraintTester : public KinFitConstraintTesterBase<P,C>
{
 public:
  KinFitConstraintTester(int verbosity = -1)
    : KinFitConstraintTesterBase<P,C>(verbosity)
  {
    rndMean_( 0) = 0.; rndWidth_( 0) = 0.0100; // pvX
    rndMean_( 1) = 0.; rndWidth_( 1) = 0.0100; // pvY
    rndMean_( 2) = 0.; rndWidth_( 2) = 0.0100; // pvX
    rndMean_( 3) = 0.; rndWidth_( 3) = 1.;     // nuTauPlusPx
    rndMean_( 4) = 0.; rndWidth_( 4) = 1.;     // nuTauPlusPy
    rndMean_( 5) = 0.; rndWidth_( 5) = 0.0100; // svTauPlusX
    rndMean_( 6) = 0.; rndWidth_( 6) = 0.0100; // svTauPlusY
    rndMean_( 7) = 0.; rndWidth_( 7) = 0.0100; // svTauPlusZ
    rndMean_( 8) = 0.; rndWidth_( 8) = 1.;     // nuTauMinusPx
    rndMean_( 9) = 0.; rndWidth_( 9) = 1.;     // nuTauMinusPy
    rndMean_(10) = 0.; rndWidth_(10) = 0.0100; // svTauPlusX
    rndMean_(11) = 0.; rndWidth_(11) = 0.0100; // svTauPlusY
    rndMean_(12) = 0.; rndWidth_(12) = 0.0100; // svTauPlusZ
    rndMean_(13) = 0.; rndWidth_(13) = 1.;     // recoilPx
    rndMean_(14) = 0.; rndWidth_(14) = 1.;     // recoilPy
    rndMean_(15) = 0.; rndWidth_(15) = 1.;     // recoilPz
    rndMean_(16) = 0.; rndWidth_(16) = 1.;     // recoilE
 
    for ( unsigned idxParameter = 0; idxParameter < P; ++idxParameter )
    {
      for ( unsigned idxConstraint = 0; idxConstraint < C; ++idxConstraint )
      {
        plotRange_(idxConstraint, idxParameter) = 3.*rndWidth_(idxParameter);
      }
    }
  }

 private:
  using KinFitConstraintTesterBase<P,C>::rndMean_;
  using KinFitConstraintTesterBase<P,C>::rndWidth_;
  using KinFitConstraintTesterBase<P,C>::plotRange_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraintTester_h
