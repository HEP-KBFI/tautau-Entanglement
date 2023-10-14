#include "TauAnalysis/Entanglement/interface/get_rotationMatrix.h"

math::Matrix3x3
get_rotationMatrix(const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k)
{
  // compute rotation matrix for transformation from laboratory frame { x, y, z } to helicity frame { r, n, k }
  math::Matrix3x3 rotMatrix;
  std::vector<reco::Candidate::Vector> rnk = { r, n, k };
  reco::Candidate::Vector x(1.,0.,0.);
  reco::Candidate::Vector y(0.,1.,0.);
  reco::Candidate::Vector z(0.,0.,1.);
  std::vector<reco::Candidate::Vector> xyz = { x, y, z };
  for ( unsigned int i = 0; i < 3; ++i )
  {
    for ( unsigned int j = 0; j < 3; ++j )
    {
      const reco::Candidate::Vector& e_i = rnk[i];
      const reco::Candidate::Vector& e_j = xyz[j];
      rotMatrix(i,j) = e_i.Dot(e_j);
    }
  }
  return rotMatrix;
}

math::Matrix3x3
get_rotationMatrixInv(const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k)
{
  // compute rotation matrix for transformation from helicity frame { r, n, k } back to laboratory frame { x, y, z }
  math::Matrix3x3 rotMatrixInv;
  std::vector<reco::Candidate::Vector> rnk = { r, n, k };
  reco::Candidate::Vector x(1.,0.,0.);
  reco::Candidate::Vector y(0.,1.,0.);
  reco::Candidate::Vector z(0.,0.,1.);
  std::vector<reco::Candidate::Vector> xyz = { x, y, z };
  for ( unsigned int i = 0; i < 3; ++i )
  {
    for ( unsigned int j = 0; j < 3; ++j )
    {
      const reco::Candidate::Vector& e_i = xyz[i];
      const reco::Candidate::Vector& e_j = rnk[j];
      rotMatrixInv(i,j) = e_i.Dot(e_j);
    }
  }
  return rotMatrixInv;
}
