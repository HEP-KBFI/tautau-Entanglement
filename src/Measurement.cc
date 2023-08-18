#include "TauAnalysis/Entanglement/interface/Measurement.h"

using namespace spin;

Measurement::Measurement(const math::Vector3& Bp, const math::Vector3& Bm, const math::Matrix3x3& C,
                         double Rchsh)
  : Bp_(Bp)
  , Bm_(Bm)
  , C_(C)
  , Rchsh_(Rchsh)
{}

Measurement::~Measurement()
{}

void
Measurement::set_BpErr(const math::Vector3& BpErr)
{
  BpErr_ = BpErr;
}

void
Measurement::set_BmErr(const math::Vector3& BmErr)
{
  BmErr_ = BmErr;
}

void
Measurement::set_CErr(const math::Matrix3x3& CErr)
{
  CErr_ = CErr;
}

void
Measurement::set_RchshErr(double Rchsh)
{
  Rchsh_ = Rchsh;
}

const math::Vector3&
Measurement::get_Bp() const
{
  return Bp_;
}

const math::Vector3&
Measurement::get_BpErr() const
{
  return BpErr_;
}

const math::Vector3&
Measurement::get_Bm() const
{
  return Bm_;
}

const math::Vector3&
Measurement::get_BmErr() const
{
  return BmErr_;
}

const math::Matrix3x3&
Measurement::get_C() const
{
  return C_;
}
 
const math::Matrix3x3&
Measurement::get_CErr() const
{
  return CErr_;
}

double
Measurement::get_Rchsh() const
{
  return Rchsh_;
}

double
Measurement::get_RchshErr() const
{
  return RchshErr_;
}
