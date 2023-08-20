#ifndef TauAnalysis_Entanglement_Measurement_h
#define TauAnalysis_Entanglement_Measurement_h

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3, math::Vector3

namespace spin
{

class Measurement
{
 public:
  Measurement();
  Measurement(const math::Vector3& Bp, const math::Vector3& Bm, const math::Matrix3x3& C, 
              double Rchsh);
  ~Measurement();

  void
  set_BpErr(const math::Vector3& BpErr);

  void
  set_BmErr(const math::Vector3& BmErr);

  void
  set_CErr(const math::Matrix3x3& CErr);

  void
  set_RchshErr(double RchshErr);

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
  get_Rchsh() const;

  double
  get_RchshErr() const;

 private:
  math::Vector3 Bp_;
  math::Vector3 BpErr_;
  math::Vector3 Bm_;
  math::Vector3 BmErr_;
  math::Matrix3x3 C_;
  math::Matrix3x3 CErr_;
  double Rchsh_;
  double RchshErr_;
};

}

#endif // TauAnalysis_Entanglement_Measurement_h
