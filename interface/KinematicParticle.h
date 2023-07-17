#ifndef TauAnalysis_Entanglement_KinematicParticle_h
#define TauAnalysis_Entanglement_KinematicParticle_h

/** KinematicParticle
 *
 * Representation of a single particle, used as input and output data format in kinematic fits.
 *
 * The origin (vertex) and direction (momentum) of the particle 
 * can be given in either of the following two parametrizations:
 *  - In the 5-parameter ("C") format
 *      alpha_C = (c, phi0, d0, lamba, z0)
 *    defined in Section 1.1 of [1], where
 *      c = 1/(2R), with R being the radius of curvature
 *      phi0 is the azimuthal angle of the track momentum vector at the point of closest approach to the reference point,
 *      d0 is the signed distance of closest approach to the reference point in the x~y plane,
 *      lambda = cot(theta), with theta being the polar angle measured from the +z axis,
 *      z0 is the z position of the track at the point of closest approach, in the x~y plane, to the reference point
 *  - In the 7-parameter ("W") format
 *      alpha_W = (px, py, pz, E, x, y, z)
 *    defined in Section 1.1 of [1], where
 *      px, py, pz, and E are the three momentum components and the energy of the particle
 *      x, y, z are the coordinates of the particle's production vertex
 *
 * [1] http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
 *
 * \authors Christian Veelken, Tallinn
 *
 */

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::LorentzVector, reco::Candidate::Point
#include "DataFormats/Math/interface/Matrix.h"         // math::Matrix
#include "DataFormats/Math/interface/Vector.h"         // math::Vector

#include <TDatabasePDG.h>                              // TDatabasePDG

#include <cmath>                                       // cos, sin, sqrt

namespace math
{
  typedef Matrix<5,5>::type Matrix5x5;
  typedef Vector<5>::type   Vector5;

  typedef Matrix<7,7>::type Matrix7x7;
  typedef Vector<7>::type   Vector7;
}

class KinematicParticle
{
 public:
  KinematicParticle(int pdgId);
  ~KinematicParticle();

  void 
  set_params5(const math::Vector5& params5, const math::Matrix5x5& cov5x5);

  void 
  set_params7(const math::Vector7& params7, const math::Matrix7x7& cov7x7);

  const reco::Candidate::LorentzVector&
  get_p4() const;

  const reco::Candidate::Point&
  get_vertex() const;

  int
  get_pdgId() const;

  double
  get_mass() const;

  int
  get_charge() const;

  const math::Vector5&
  get_params5() const;

  const math::Matrix5x5&
  get_cov5x5() const;

  const math::Vector7&
  get_params7() const;

  const math::Matrix7x7&
  get_cov7x7() const;

  friend class Smearing;

 private:
  static TDatabasePDG* pdg_;

  reco::Candidate::LorentzVector p4_;

  reco::Candidate::Point vertex_;

  int pdgId_;
  double mass_;
  int charge_;

  math::Vector5 params5_;
  math::Matrix5x5 cov5x5_;
  bool params5_isValid_ = false;

  math::Vector7 params7_;
  math::Matrix7x7 cov7x7_;
  bool params7_isValid_ = false;
};

void
printKinematicParticle(const std::string& label,
                       const KinematicParticle& particle,
                       bool cartesian = true);

#endif // TauAnalysis_Entanglement_KinematicParticle_h
