#include "TauAnalysis/Entanglement/interface/comp_PCA_line2line.h"

#include "DataFormats/Math/interface/Matrix.h"                // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                // math::Vector

#include "TauAnalysis/Entanglement/interface/constants.h"     // ct, mTau
#include "TauAnalysis/Entanglement/interface/cmsException.h"  // cmsException
#include "TauAnalysis/Entanglement/interface/printDistance.h" // printDistance()
#include "TauAnalysis/Entanglement/interface/printPoint.h"    // printPoint()

#include <cmath>                                              // std::sqrt()

namespace math
{
  typedef Matrix<3,3>::type Matrix3x3;
  typedef Vector<3>::type   Vector3;
}

std::pair<reco::Candidate::Point, reco::Candidate::Point>
comp_PCA_line2line(const reco::Candidate::Point& P1,
                   const reco::Candidate::Vector& V1,
                   const reco::Candidate::Point& P2,
                   const reco::Candidate::Vector& V2,
                   int verbosity)
{
  // CV: compute point of closest approach (PCA) between two straight lines in three dimensions;
  //     code based on https://math.stackexchange.com/questions/1993953/closest-points-between-two-lines
  if ( verbosity >= 4 )
  {
    std::cout << "<comp_PCA_line2line>:" << std::endl;
  }
  auto e1 = V1.unit();
  auto e2 = V2.unit();
  reco::Candidate::Vector d = e1.Cross(e2).unit();
  math::Matrix3x3 v;
  v(0,0) =  e1.x();
  v(0,1) = -e2.x();
  v(0,2) =  d.x();
  v(1,0) =  e1.y();
  v(1,1) = -e2.y();
  v(1,2) =  d.y();
  v(2,0) =  e1.z();
  v(2,1) = -e2.z();
  v(2,2) =  d.z();
  if ( verbosity >= 4 )
  {
    std::cout << "v:\n";
    std::cout << v << "\n";

  }
  math::Vector3 r;
  r(0) = -P1.x() + P2.x();
  r(1) = -P1.y() + P2.y();
  r(2) = -P1.z() + P2.z();
  if ( verbosity >= 2 )
  {
    std::cout << "r:\n";
    std::cout << r << "\n";

  }
  // CV: invert matrix x;
  //     see Section "Linear algebra functions" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html for the syntax
  int errorFlag = 0;
  math::Matrix3x3 vinv = v.Inverse(errorFlag);
  if ( errorFlag != 0 )
    throw cmsException("comp_PCA_line2line", __LINE__)
       << "Failed to invert matrix v !!\n";
  math::Vector3 t = vinv*r;
  if ( verbosity >= 4 )
  {
    std::cout << "t:\n";
    std::cout << t << "\n";
  }

  reco::Candidate::Point Q1 = P1 + t(0)*e1;
  reco::Candidate::Point Q2 = P2 + t(1)*e2;
  if ( verbosity >= 4 )
  {
    printPoint("Q1", Q1);
    printDistance("Q1 - P1", Q1 - P1, true);
    printDistance("Q1 - P1", Q1 - P1, false);
    printPoint("Q2", Q2);
    printDistance("Q2 - P2", Q2 - P2, true);
    printDistance("Q2 - P2", Q2 - P2, false);
  }

  return std::pair<reco::Candidate::Point, reco::Candidate::Point>(Q1, Q2);
}

reco::Candidate::Point
comp_PCA_line2line(const reco::Candidate::Point& pv,
                   const reco::Candidate::LorentzVector& tauP4,
                   const reco::Candidate::Point& sv,
                   const reco::Candidate::LorentzVector& visTauP4,
                   int verbosity)
{
  std::pair<reco::Candidate::Point, reco::Candidate::Point> Qs = comp_PCA_line2line(pv, tauP4.Vect(), sv, visTauP4.Vect(), verbosity);
  const reco::Candidate::Point& Q1 = Qs.first;
  double t = std::sqrt((Q1 - pv).mag2());

  auto e_tau = tauP4.Vect().unit();
  double min_t = 1.e-2*(tauP4.energy()/mTau)*ct;
  if ( t < min_t ) t = min_t;
  if ( verbosity >= 4 )
  {
    std::cout << "min_t = " << min_t << "\n";
    std::cout << "t = " << t << "\n";
  }
  reco::Candidate::Point pca = pv + t*e_tau;
  if ( verbosity >= 4 )
  {
    printPoint("pca", pca);
    printDistance("pca - pv", pca - pv, true);
    printDistance("pca - pv", pca - pv, false);
  }
  return pca;
}
