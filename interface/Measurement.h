#ifndef TauAnalysis_Entanglement_Measurement_h
#define TauAnalysis_Entanglement_Measurement_h

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3, math::Vector3

namespace spin
{

class Measurement
{
 public:
  Measurement();
  Measurement(const math::Vector3& Bp, const math::Vector3& Bm, const math::Matrix3x3& C);
  ~Measurement();

  void
  set_Bp_medianErr(const math::Vector3& BpMedian,
                   const math::Vector3& BpErr);

  void
  set_Bm_medianErr(const math::Vector3& BmMedian,
                   const math::Vector3& BmErr);

  void
  set_C_medianErr(const math::Matrix3x3& CMedian,
                  const math::Matrix3x3& CErr);
  
  void
  set_concurrence_medianErr(double concurrenceMedian,
                            double concurrenceErr);

  void
  set_Ek_medianErr(double EkMedian,
                   double EkErr);

  void
  set_Rchsh_medianErr(double RchshMedian,
                      double RchshErr);

  void
  set_steerability_medianErr(double steerabilityMedian,
                             double steerabilityErr);

  void
  set_metadata(std::size_t sample_size,
               std::size_t num_bootstrap_samples,
               std::size_t boostrap_size,
               std::size_t num_concurrence_failures);

  const math::Vector3&
  get_Bp(bool isMedian = false) const;

  const math::Vector3&
  get_BpErr() const;

  const math::Vector3&
  get_Bm(bool isMedian = false) const;

  const math::Vector3&
  get_BmErr() const;

  const math::Matrix3x3&
  get_C(bool isMedian = false) const;

  const math::Matrix3x3&
  get_CErr() const;

  double
  get_concurrence(bool isMedian = false) const;

  double
  get_concurrenceErr() const;

  double
  get_Ek(bool isMedian = false) const;

  double
  get_EkErr() const;

  double
  get_Rchsh(bool isMedian = false) const;

  double
  get_RchshErr() const;

  double
  get_steerability(bool isMedian = false) const;

  double
  get_steerabilityErr() const;

  std::size_t
  get_sample_size() const;

  std::size_t
  get_num_bootstrap_samples() const;

  std::size_t
  get_boostrap_size() const;

  std::size_t
  get_num_concurrence_failures() const;

  friend class SpinAlgoBase;

 private:
  math::Vector3 Bp_;
  math::Vector3 BpMedian_;
  math::Vector3 BpErr_;
  math::Vector3 Bm_;
  math::Vector3 BmMedian_;
  math::Vector3 BmErr_;
  math::Matrix3x3 C_;
  math::Matrix3x3 CMedian_;
  math::Matrix3x3 CErr_;

  double concurrence_;
  double concurrenceMedian_;
  double concurrenceErr_;
  double Ek_;
  double EkMedian_;
  double EkErr_;
  double Rchsh_;
  double RchshMedian_;
  double RchshErr_;
  double steerability_;
  double steerabilityMedian_;
  double steerabilityErr_;

  std::size_t sample_size_;
  std::size_t num_bootstrap_samples_;
  std::size_t boostrap_size_;
  std::size_t num_concurrence_failures_;
};

}

#endif // TauAnalysis_Entanglement_Measurement_h
