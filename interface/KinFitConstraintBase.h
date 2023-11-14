#ifndef TauAnalysis_Entanglement_KinFitConstraintBase_h
#define TauAnalysis_Entanglement_KinFitConstraintBase_h

#include <TMatrixD.h> // TMatrixD
#include <TVectorD.h> // TVectorD

class KinFitConstraintBase
{
 public:
  KinFitConstraintBase(unsigned int Np, unsigned int Nc_eq, unsigned int Nc_ineq = 0, int verbosity = -1);
  virtual ~KinFitConstraintBase();

  unsigned int
  get_Np() const;

  unsigned int
  get_Nc_eq() const;

  unsigned int
  get_Nc_ineq() const;

  virtual void
  set_alphaA(const TVectorD& alphaA) = 0;

  const TMatrixD&
  get_D_eq() const;

  const TVectorD&
  get_d_eq() const;

  const TMatrixD&
  get_d_eq_metric() const;

  const TMatrixD&
  get_D_ineq() const;

  const TVectorD&
  get_d_ineq() const;

  const TMatrixD&
  get_d_ineq_metric() const;

  bool
  get_errorFlag() const;
  
 protected:
  unsigned int Np_;
  unsigned int Nc_eq_;
  unsigned int Nc_ineq_;

  TMatrixD D_eq_;
  TVectorD d_eq_;
  TMatrixD d_eq_metric_;

  TMatrixD D_ineq_;
  TVectorD d_ineq_;
  TMatrixD d_ineq_metric_;

  bool errorFlag_;
  
  int verbosity_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraintBase_h
