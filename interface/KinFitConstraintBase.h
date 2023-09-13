#ifndef TauAnalysis_Entanglement_KinFitConstraintBase_h
#define TauAnalysis_Entanglement_KinFitConstraintBase_h

#include "DataFormats/Math/interface/Matrix.h" // math::Matrix
#include "DataFormats/Math/interface/Vector.h" // math::Vector

// CV: The template parameters P and C represent 
//     the number of parameters and the number of constraints, respectively.
//     It is neccessary to make the number of parameters and constraints template parameters,
//     because the ROOT::Math::SVector and ROOT::Math::SMatrix classes
//     on which the kinematic fitting code is based require these template parameters.

template <unsigned int P, unsigned int C>
class KinFitConstraintBase
{
 public:
  typedef typename math::Matrix<P,P>::type MatrixPxP;
  typedef typename math::Matrix<P,C>::type MatrixPxC;
  typedef typename math::Matrix<C,P>::type MatrixCxP;
  typedef typename math::Matrix<C,C>::type MatrixCxC;
  typedef typename math::Vector<P>::type VectorP;
  typedef typename math::Vector<C>::type VectorC;

  KinFitConstraintBase(int verbosity = -1)
    : errorFlag_(false)
    , verbosity_(verbosity)
  {}
  virtual ~KinFitConstraintBase()
  {}

  virtual void
  set_alphaA(const KinFitConstraintBase::VectorP& alphaA) = 0;

  const MatrixCxP&
  get_D() const
  {
    return D_;
  }

  const VectorC&
  get_d() const
  {
    return d_;
  }

  const MatrixCxC&
  get_d_metric() const
  {
    return d_metric_;
  }

  bool
  get_errorFlag() const
  {
    return errorFlag_;
  }
  
 protected:
  MatrixCxP D_;
  VectorC d_;
  MatrixCxC d_metric_;

  bool errorFlag_;
  
  int verbosity_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraintBase_h
