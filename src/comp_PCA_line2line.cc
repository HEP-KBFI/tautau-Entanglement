#include "TauAnalysis/Entanglement/interface/comp_PCA_line2line.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"   // cmsException
#include "TauAnalysis/Entanglement/interface/printCovMatrix.h" // printCovMatrix()
#include "TauAnalysis/Entanglement/interface/printPoint.h"     // printPoint()

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

std::pair<reco::Candidate::Point, reco::Candidate::Point>
comp_PCA_line2line(const reco::Candidate::Point& P1, const reco::Candidate::Vector& V1,
                   const reco::Candidate::Point& P2, const reco::Candidate::Vector& V2,
                   const math::Matrix3x3* cov,
                   double lambda1Min, double lambda1Max, double lambda2Min, double lambda2Max, 
                   int verbosity)
{
  // CV: compute point of closest approach (PCA) between two straight lines in three dimensions;
  //     code based on https://math.stackexchange.com/questions/1993953/closest-points-between-two-lines
  if ( verbosity >= 4 )
  {
    std::cout << "<comp_PCA_line2line>:" << std::endl;
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
        printCovMatrix("cov", *cov);
      }
      throw cmsException("comp_PCA_line2line", __LINE__) 
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
  math::Vector3 e2 = convert_to_mathVector(V2.unit());

  math::Matrix2x2 A;
  A(0,0) =  ROOT::Math::Dot(e1, covInv*e1);
  A(0,1) = -ROOT::Math::Dot(e1, covInv*e2);
  A(1,0) =  ROOT::Math::Dot(e2, covInv*e1);
  A(1,1) = -ROOT::Math::Dot(e2, covInv*e2);
  int errorFlag = 0;
  math::Matrix2x2 Ainv = A.Inverse(errorFlag);
  if ( errorFlag != 0 )
  {
    if ( verbosity >= 0 )
    {
      printCovMatrix("A", A);
    }
    throw cmsException("comp_PCA_line2line", __LINE__) 
      << "Failed to invert matrix A !!\n";
  }

  math::Vector2 c;
  c(0) = ROOT::Math::Dot(e1, covInv*convert_to_mathVector(P2 - P1));
  c(1) = ROOT::Math::Dot(e2, covInv*convert_to_mathVector(P2 - P1));

  math::Vector2 lambda = Ainv*c;
  double lambda1 = lambda(0);
  if ( lambda1 < lambda1Min ) lambda1 = lambda1Min;
  if ( lambda1 > lambda1Max ) lambda1 = lambda1Max;
  double lambda2 = lambda(1);
  if ( lambda2 < lambda2Min ) lambda2 = lambda2Min;
  if ( lambda2 > lambda2Max ) lambda2 = lambda2Max;
  if ( verbosity >= 4 )
  {
    std::cout << "lambda1 = " << lambda1 << "\n";
    std::cout << "lambda2 = " << lambda2 << "\n";
  }

  reco::Candidate::Point pca1 = P1 + lambda1*convert_to_recoVector(e1);
  reco::Candidate::Point pca2 = P2 + lambda2*convert_to_recoVector(e2);
  if ( verbosity >= 4 )
  {
    printPoint("pca1", pca1);
    printPoint("pca2", pca2);
  }
  return std::pair<reco::Candidate::Point, reco::Candidate::Point>(pca1, pca2);
}
