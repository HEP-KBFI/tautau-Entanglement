#include "TauAnalysis/Entanglement/interface/comp_PCA_line2line.h"

#include "DataFormats/Math/interface/Matrix.h"                // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                // math::Vector

#include "TauAnalysis/Entanglement/interface/constants.h"     // ct, mTau
#include "TauAnalysis/Entanglement/interface/cmsException.h"  // cmsException
#include "TauAnalysis/Entanglement/interface/printDistance.h" // printDistance()
#include "TauAnalysis/Entanglement/interface/printPoint.h"    // printPoint()

namespace math
{
  typedef Matrix<3,3>::type Matrix3x3;
  typedef Vector<3>::type   Vector3;
}

reco::Candidate::Point
comp_PCA_line2line(const reco::Candidate::Point& pv,
                   const reco::Candidate::LorentzVector& tauP4,
                   const reco::Candidate::Point& sv,
                   const reco::Candidate::LorentzVector& visTauP4,
                   int verbosity)
{
  // CV: compute point of closest approach (PCA) between two straight lines in three dimensions;
  //     code based on https://math.stackexchange.com/questions/1993953/closest-points-between-two-lines
  if ( verbosity >= 2 )
  {
    std::cout << "<comp_PCA_line2line>:" << std::endl;
  }
  auto e_tau = tauP4.Vect().unit();
  auto e_vis = visTauP4.Vect().unit();
  reco::Candidate::Vector d = e_tau.Cross(e_vis).unit();
  math::Matrix3x3 v;
  v(0,0) =  e_tau.x();
  v(0,1) = -e_vis.x();
  v(0,2) =  d.x();
  v(1,0) =  e_tau.y();
  v(1,1) = -e_vis.y();
  v(1,2) =  d.y();
  v(2,0) =  e_tau.z();
  v(2,1) = -e_vis.z();
  v(2,2) =  d.z();
  if ( verbosity >= 2 )
  {
    std::cout << "v:\n";
    std::cout << v << "\n";

  }
  math::Vector3 r;
  r(0) = -pv.x() + sv.x();
  r(1) = -pv.y() + sv.y();
  r(2) = -pv.z() + sv.z();
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
  math::Vector3 lambda = vinv*r;
  if ( verbosity >= 2 )
  {
    std::cout << "lambda:\n";
    std::cout << lambda << "\n";
    reco::Candidate::Point pca1 = pv + lambda(0)*e_tau;
    printPoint("pca1", pca1);
    printDistance("pca1 - pv", pca1 - pv, true);
    printDistance("pca1 - pv", pca1 - pv, false);
    reco::Candidate::Point pca2 = sv + lambda(1)*e_vis;
    printPoint("pca2", pca2);
    printDistance("pca2 - sv", pca2 - sv, true);
    printDistance("pca2 - sv", pca2 - sv, false);
    std::cout << "|d| = " << std::sqrt((lambda(2)*d).mag2()) << "\n";      
  }
  double min_lambda0 = 1.e-2*(tauP4.energy()/mTau)*ct;
  double lambda0 = ( lambda(0) >= min_lambda0 ) ? lambda(0) : min_lambda0;
  if ( verbosity >= 2 )
  {
    std::cout << "min_lambda0 = " << min_lambda0 << "\n";
    std::cout << "lambda0 = " << lambda0 << "\n";
  }
  reco::Candidate::Point pca = pv + lambda0*e_tau;
  if ( verbosity >= 2 )
  {
    printPoint("pca", pca);
    printDistance("pca - pv", pca - pv, true);
    printDistance("pca - pv", pca - pv, false);
  }
  return pca;
}
