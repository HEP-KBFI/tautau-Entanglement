#include "TauAnalysis/Entanglement/interface/Measurement.h"

using namespace spin;

Measurement::Measurement()
  : Rchsh_(0.)
  , RchshErr_(0.)
{}

Measurement::Measurement(const math::Vector3& Bp, const math::Vector3& Bm, const math::Matrix3x3& C)
  : Bp_(Bp)
  , Bm_(Bm)
  , C_(C)
  , concurrence_(0.)
  , concurrenceErr_(0.)
  , Ek_(0.)
  , EkErr_(0.)
  , Rchsh_(0.)
  , RchshErr_(0.)
  , steerability_(0.)
  , steerabilityErr_(0.)
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
Measurement::set_concurrenceErr(double concurrenceErr)
{
  concurrenceErr_ = concurrenceErr;
}

void
Measurement::set_EkErr(double EkErr)
{
  EkErr_ = EkErr;
}

void
Measurement::set_RchshErr(double RchshErr)
{
  RchshErr_ = RchshErr;
}

void
Measurement::set_steerabilityErr(double steerabilityErr)
{
  steerabilityErr_ = steerabilityErr;
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
Measurement::get_concurrence() const
{
  return concurrence_;
}

double
Measurement::get_concurrenceErr() const
{
  return concurrenceErr_;
}

double
Measurement::get_Ek() const
{
  return Ek_;
}

double
Measurement::get_EkErr() const
{
  return EkErr_;
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

double
Measurement::get_steerability() const
{
  return steerability_;
}

double
Measurement::get_steerabilityErr() const
{
  return steerabilityErr_;
}
