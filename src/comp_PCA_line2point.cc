#include "TauAnalysis/Entanglement/interface/comp_PCA_line2point.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/printMatrix.h"  // printMatrix()
#include "TauAnalysis/Entanglement/interface/printPoint.h"   // printPoint()

namespace
{
  template <typename T>
  math::Vector3
  convert_to_mathVector(const T& V)
  {
    math::Vector3 retVal;
    retVal(0) = V.x();
    retVal(1) = V.y();
    retVal(2) = V.z();
    return retVal;
  }

  reco::Candidate::Vector
  convert_to_recoVector(const math::Vector3& V)
  {
    return reco::Candidate::Vector(V(0), V(1), V(2));
  }
}

reco::Candidate::Point
comp_PCA_line2point(const reco::Candidate::Point& P1, const reco::Candidate::Vector& V1,
                    const reco::Candidate::Point& P2,
                    const math::Matrix3x3* cov,
                    double lambdaMin, double lambdaMax, 
                    int verbosity)
{
  // CV: compute point of closest approach (PCA) between line (P1 + lambda*V1) and point (P2)
  if ( verbosity >= 4 )
  {
    std::cout << "<comp_PCA_line2point>:" << std::endl;
  }

  math::Matrix3x3 covInv;
  if ( cov )
  {
    // CV: compute inverse of covariance matrix
    int errorFlag = 0;
    covInv = cov->Inverse(errorFlag);
    if ( errorFlag != 0 )
    {
      if ( verbosity >= 0 )
      {
        printMatrix("cov", *cov, true);
      }
      throw cmsException("comp_PCA_line2point", __LINE__) 
        << "Failed to invert matrix cov !!\n";
    }
  }
  else
  {
    // CV: set inverse of covariance matrix to unit matrix
    for ( int idxRow = 0; idxRow < 3; ++idxRow )
    {
      for ( int idxColumn = 0; idxColumn < 3; ++idxColumn )
      {
        if ( idxRow == idxColumn ) covInv(idxRow,idxColumn) = 1.;
        else                       covInv(idxRow,idxColumn) = 0.;
      }
    }
  }

  math::Vector3 e1 = convert_to_mathVector(V1.unit());

  double lambda = ROOT::Math::Dot(convert_to_mathVector(P2 - P1), covInv*e1)/ROOT::Math::Dot(e1, covInv*e1);
  if ( lambda < lambdaMin ) lambda = lambdaMin;
  if ( lambda > lambdaMax ) lambda = lambdaMax;
  if ( verbosity >= 4 )
  {
    std::cout << "lambda = " << lambda << "\n";
  }

  reco::Candidate::Point pca = P1 + lambda*convert_to_recoVector(e1);
  if ( verbosity >= 4 )
  {
    printPoint("pca", pca);
  }
  return pca;
}
