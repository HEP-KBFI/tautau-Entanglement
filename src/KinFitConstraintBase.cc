#include "TauAnalysis/Entanglement/interface/KinFitConstraintBase.h"

KinFitConstraintBase::KinFitConstraintBase(unsigned int Np, unsigned int Nc_eq, unsigned int Nc_ineq, int verbosity)
  : Np_(Np)
  , Nc_eq_(Nc_eq)
  , Nc_ineq_(Nc_ineq)
  , D_eq_(Nc_eq,Np)
  , d_eq_(Nc_eq)
  , d_eq_metric_(Nc_eq,Nc_eq)
  , D_ineq_(Nc_ineq,Np)
  , d_ineq_(Nc_ineq)
  , d_ineq_metric_(Nc_ineq,Nc_ineq)
  , errorFlag_(false)
  , verbosity_(verbosity)
{}

KinFitConstraintBase::~KinFitConstraintBase()
{}

unsigned int
KinFitConstraintBase::get_Np() const
{
  return Np_;
}

unsigned int
KinFitConstraintBase::get_Nc_eq() const
{
  return Nc_eq_;
}

unsigned int
KinFitConstraintBase::get_Nc_ineq() const
{
  return Nc_ineq_;
}

const TMatrixD&
KinFitConstraintBase::get_D_eq() const
{
  return D_eq_;
}

const TVectorD&
KinFitConstraintBase::get_d_eq() const
{
  return d_eq_;
}

const TMatrixD&
KinFitConstraintBase::get_d_eq_metric() const
{
  return d_eq_metric_;
}

const TMatrixD&
KinFitConstraintBase::get_D_ineq() const
{
  return D_ineq_;
}

const TVectorD&
KinFitConstraintBase::get_d_ineq() const
{
  return d_ineq_;
}

const TMatrixD&
KinFitConstraintBase::get_d_ineq_metric() const
{
  return d_ineq_metric_;
}

bool
KinFitConstraintBase::get_errorFlag() const
{
  return errorFlag_;
}
  
