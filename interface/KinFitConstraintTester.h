#ifndef TauAnalysis_Entanglement_KinFitConstraintTester_h
#define TauAnalysis_Entanglement_KinFitConstraintTester_h

#include "TauAnalysis/Entanglement/interface/KinFitConstraintTesterBase.h" // KinFitConstraintTesterBase

#include <vector>                                                          // std::vector<>

bool
isIn(unsigned int givenIndex, const std::vector<unsigned int>& selIndices)
{
  for ( unsigned int selIndex : selIndices )
  {
    if ( givenIndex == selIndex ) return true;
  }
  return false;
}

template <unsigned int P, unsigned int C>
class KinFitConstraintTester : public KinFitConstraintTesterBase<P,C>
{
 public:
  KinFitConstraintTester(const typename KinFitConstraintBase<P,C>::VectorP& alpha0, int verbosity = -1)
    : KinFitConstraintTesterBase<P,C>(alpha0, verbosity)
  {
    rndMean_ = alpha0_;
    
    rndWidth_( 0) = 0.0100; // pvX
    rndWidth_( 1) = 0.0100; // pvY
    rndWidth_( 2) = 0.0100; // pvX
    rndWidth_( 3) = 0.1;    // nuTauPlusPx
    rndWidth_( 4) = 0.1;    // nuTauPlusPy
    rndWidth_( 5) = 0.0100; // svTauPlusX
    rndWidth_( 6) = 0.0100; // svTauPlusY
    rndWidth_( 7) = 0.0100; // svTauPlusZ
    rndWidth_( 8) = 0.1;    // nuTauMinusPx
    rndWidth_( 9) = 0.1;    // nuTauMinusPy
    rndWidth_(10) = 0.0100; // svTauPlusX
    rndWidth_(11) = 0.0100; // svTauPlusY
    rndWidth_(12) = 0.0100; // svTauPlusZ
    rndWidth_(13) = 1.;     // recoilPx
    rndWidth_(14) = 1.;     // recoilPy
    rndWidth_(15) = 1.;     // recoilPz
    rndWidth_(16) = 1.;     // recoilE
 
    for ( unsigned int idxParameter = 0; idxParameter < P; ++idxParameter )
    {
      for ( unsigned int idxConstraint = 0; idxConstraint < C; ++idxConstraint )
      {
        if (  isIn(idxParameter, std::vector<unsigned int>{ 0, 1, 2, 5, 6, 7, 10, 11, 12 })                        ||
             (isIn(idxParameter, std::vector<unsigned int>{ 3, 4, 8, 9                   }) && idxConstraint >= 4) )
        {
          plotRange_(idxConstraint, idxParameter) = 1.e+2*rndWidth_(idxParameter);
        }
        else
        {
          plotRange_(idxConstraint, idxParameter) = 1.e+1*rndWidth_(idxParameter);
        }
      }
    }
  }

 private:
  using KinFitConstraintTesterBase<P,C>::alpha0_;
  using KinFitConstraintTesterBase<P,C>::rndMean_;
  using KinFitConstraintTesterBase<P,C>::rndWidth_;
  using KinFitConstraintTesterBase<P,C>::plotRange_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraintTester_h
