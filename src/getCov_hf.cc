#include "TauAnalysis/Entanglement/interface/getCov_hf.h"

#include "TauAnalysis/Entanglement/interface/square.h" // square()

math::Matrix3x3
getCov_hf(double dk, double dr, double dn,
          const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k)
{
  double dk2 = square(dk);
  double dr2 = square(dr);
  double dn2 = square(dn);
  reco::Candidate::Vector x(1., 0., 0.);
  reco::Candidate::Vector y(0., 1., 0.);
  reco::Candidate::Vector z(0., 0., 1.);
  double k_x = k.Dot(x);
  double k_y = k.Dot(y);
  double k_z = k.Dot(z);
  double r_x = r.Dot(x);
  double r_y = r.Dot(y);
  double r_z = r.Dot(z);
  double n_x = n.Dot(x);
  double n_y = n.Dot(y);
  double n_z = n.Dot(z);
  // CV: computation of covariance matrix in rotated coordinates taken from the paper
  //       "On transformation of covariance matrices between local Cartesian cooridinate systems and commutative diagrams",
  //       T. Soler and M. Chin; 
  //     also see Appendix I of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
  math::Matrix3x3 cov; 
  cov(0,0) = dk2*k_x*k_x + dr2*r_x*r_x + dn2*n_x*n_x;
  cov(0,1) = dk2*k_x*k_y + dr2*r_x*r_y + dn2*n_x*n_y;
  cov(0,2) = dk2*k_x*k_z + dr2*r_x*r_z + dn2*n_x*n_z;
  cov(1,0) = dk2*k_y*k_x + dr2*r_y*r_x + dn2*n_y*n_x;
  cov(1,1) = dk2*k_y*k_y + dr2*r_y*r_y + dn2*n_y*n_y;
  cov(1,2) = dk2*k_y*k_z + dr2*r_y*r_z + dn2*n_y*n_z;
  cov(2,0) = dk2*k_z*k_x + dr2*r_z*r_x + dn2*n_z*n_x;
  cov(2,1) = dk2*k_z*k_y + dr2*r_z*r_y + dn2*n_z*n_y;
  cov(2,2) = dk2*k_z*k_z + dr2*r_z*r_z + dn2*n_z*n_z;
  return cov;
}
