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
  set_BpErr(const math::Vector3& BpErr);

  void
  set_BmErr(const math::Vector3& BmErr);

  void
  set_CErr(const math::Matrix3x3& CErr);
  
  void
  set_concurrenceErr(double concurrenceErr);

  void
  set_EkErr(double EkErr);

  void
  set_RchshErr(double RchshErr);

  void
  set_steerabilityErr(double steerabilityErr);

  const math::Vector3&
  get_Bp() const;

  const math::Vector3&
  get_BpErr() const;

  const math::Vector3&
  get_Bm() const;

  const math::Vector3&
  get_BmErr() const;

  const math::Matrix3x3&
  get_C() const;

  const math::Matrix3x3&
  get_CErr() const;

  double
  get_concurrence() const;

  double
  get_concurrenceErr() const;

  double
  get_Ek() const;

  double
  get_EkErr() const;

  double
  get_Rchsh() const;

  double
  get_RchshErr() const;

  double
  get_steerability() const;

  double
  get_steerabilityErr() const;

  friend class SpinAlgoBase;

 private:
  math::Vector3 Bp_;
  math::Vector3 BpErr_;
  math::Vector3 Bm_;
  math::Vector3 BmErr_;
  math::Matrix3x3 C_;
  math::Matrix3x3 CErr_;
  double concurrence_;
  double concurrenceErr_;
  double Ek_;
  double EkErr_;
  double Rchsh_;
  double RchshErr_;
  double steerability_;
  double steerabilityErr_;
};

}

#endif // TauAnalysis_Entanglement_Measurement_h
