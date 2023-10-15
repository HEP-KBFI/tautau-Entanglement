#include "TauAnalysis/Entanglement/interface/Measurement.h"

using namespace spin;

Measurement::Measurement()
  : concurrence_(0.)
  , concurrenceMedian_(0.)
  , concurrenceErr_(0.)
  , Ek_(0.)
  , EkMedian_(0.)
  , EkErr_(0.)
  , Rchsh_(0.)
  , RchshMedian_(0.)
  , RchshErr_(0.)
  , steerability_(0.)
  , steerabilityMedian_(0.)
  , steerabilityErr_(0.)
  , sample_size_(0)
  , num_bootstrap_samples_(0)
  , boostrap_size_(0)
  , num_concurrence_failures_(0)
{}

Measurement::Measurement(const math::Vector3& Bp, const math::Vector3& Bm, const math::Matrix3x3& C)
  : Bp_(Bp)
  , Bm_(Bm)
  , C_(C)
  , concurrence_(0.)
  , concurrenceMedian_(0.)
  , concurrenceErr_(0.)
  , Ek_(0.)
  , EkMedian_(0.)
  , EkErr_(0.)
  , Rchsh_(0.)
  , RchshMedian_(0.)
  , RchshErr_(0.)
  , steerability_(0.)
  , steerabilityMedian_(0.)
  , steerabilityErr_(0.)
  , sample_size_(0)
  , num_bootstrap_samples_(0)
  , boostrap_size_(0)
  , num_concurrence_failures_(0)
{}

Measurement::~Measurement()
{}

void
Measurement::set_Bp_medianErr(const math::Vector3& BpMedian,
                              const math::Vector3& BpErr)
{
  BpMedian_ = BpMedian;
  BpErr_ = BpErr;
}

void
Measurement::set_Bm_medianErr(const math::Vector3& BmMedian,
                              const math::Vector3& BmErr)
{
  BmMedian_ = BmMedian;
  BmErr_ = BmErr;
}

void
Measurement::set_C_medianErr(const math::Matrix3x3& CMedian,
                             const math::Matrix3x3& CErr)
{
  CMedian_ = CMedian;
  CErr_ = CErr;
}

void
Measurement::set_concurrence_medianErr(double concurrenceMedian,
                                       double concurrenceErr)
{
  concurrenceMedian_ = concurrenceMedian;
  concurrenceErr_ = concurrenceErr;
}

void
Measurement::set_Ek_medianErr(double EkMedian,
                              double EkErr)
{
  EkMedian_ = EkMedian;
  EkErr_ = EkErr;
}

void
Measurement::set_Rchsh_medianErr(double RchshMedian,
                                 double RchshErr)
{
  RchshMedian_ = RchshMedian;
  RchshErr_ = RchshErr;
}

void
Measurement::set_steerability_medianErr(double steerabilityMedian,
                                        double steerabilityErr)
{
  steerabilityMedian_ = steerabilityMedian;
  steerabilityErr_ = steerabilityErr;
}

void
Measurement::set_metadata(std::size_t sample_size,
             std::size_t num_bootstrap_samples,
             std::size_t boostrap_size,
             std::size_t num_concurrence_failures)
{
  sample_size_ = sample_size;
  num_bootstrap_samples_ = num_bootstrap_samples;
  boostrap_size_ = boostrap_size;
  num_concurrence_failures_ = num_concurrence_failures;
}

const math::Vector3&
Measurement::get_Bp(bool isMedian) const
{
  return isMedian ? BpMedian_ : Bp_;
}

const math::Vector3&
Measurement::get_BpErr() const
{
  return BpErr_;
}

const math::Vector3&
Measurement::get_Bm(bool isMedian) const
{
  return isMedian ? BmMedian_ : Bm_;
}

const math::Vector3&
Measurement::get_BmErr() const
{
  return BmErr_;
}

const math::Matrix3x3&
Measurement::get_C(bool isMedian) const
{
  return isMedian ? CMedian_ : C_;
}
 
const math::Matrix3x3&
Measurement::get_CErr() const
{
  return CErr_;
}

double
Measurement::get_concurrence(bool isMedian) const
{
  return isMedian ? concurrenceMedian_ : concurrence_;
}

double
Measurement::get_concurrenceErr() const
{
  return concurrenceErr_;
}

double
Measurement::get_Ek(bool isMedian) const
{
  return isMedian ? EkMedian_ : Ek_;
}

double
Measurement::get_EkErr() const
{
  return EkErr_;
}

double
Measurement::get_Rchsh(bool isMedian) const
{
  return isMedian ? RchshMedian_ : Rchsh_;
}

double
Measurement::get_RchshErr() const
{
  return RchshErr_;
}

double
Measurement::get_steerability(bool isMedian) const
{
  return isMedian ? steerabilityMedian_ : steerability_;
}

double
Measurement::get_steerabilityErr() const
{
  return steerabilityErr_;
}

std::size_t
Measurement::get_sample_size() const
{
  return sample_size_;
}

std::size_t
Measurement::get_num_bootstrap_samples() const
{
  return num_bootstrap_samples_;
}

std::size_t
Measurement::get_boostrap_size() const
{
  return boostrap_size_;
}

std::size_t
Measurement::get_num_concurrence_failures() const
{
  return num_concurrence_failures_;
}
