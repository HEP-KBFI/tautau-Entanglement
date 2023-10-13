#ifndef TauAnalysis_Entanglement_KinFitAlgo_h
#define TauAnalysis_Entanglement_KinFitAlgo_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"              // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/cmsException.h"         // cmsException
#include "TauAnalysis/Entanglement/interface/KinFitConstraintBase.h" // KinFitConstraintBase
#include "TauAnalysis/Entanglement/interface/printCovMatrix.h"       // printCovMatrix()

#include "DataFormats/Math/interface/Matrix.h"                       // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                       // math::Vector

#include <cmath>                                                     // std::isnan(), std::sqrt()

// CV: The template parameters P and C represent 
//     the number of parameters and the number of constraints, respectively.
//     It is neccessary to make the number of parameters and constraints template parameters,
//     because the ROOT::Math::SVector and ROOT::Math::SMatrix classes
//     on which the kinematic fitting code is based require these template parameters.

template <unsigned int P>
class KinFitSummary
{
 public:
  typedef typename math::Matrix<P,P>::type MatrixPxP;
  typedef typename math::Vector<P>::type VectorP;

  KinFitSummary()
    : iteration_(-1)
    , chi2_(-1.)
    , status_(-1)
  {}
  
  KinFitSummary(int iteration, const VectorP& alpha, const MatrixPxP& V_alpha, double chi2, int status)
    : iteration_(iteration)
    , alpha_(alpha)
    , V_alpha_(V_alpha)
    , chi2_(chi2)
    , status_(status)
  {}
  ~KinFitSummary()
  {}

  int
  get_iteration() const
  {
    return iteration_;
  }

  const VectorP&
  get_alpha() const
  {
    return alpha_;
  }

  const MatrixPxP&
  get_V_alpha() const
  {
    return V_alpha_;
  }

  double
  get_chi2() const
  {
    return chi2_;
  }

  int
  get_status() const
  {
    return status_;
  }

 protected:
  int iteration_;
  VectorP alpha_;
  MatrixPxP V_alpha_;
  double chi2_;
  int status_;
};

namespace
{
  template <typename T>
  bool
  isNaN(const T& alpha)
  {
    size_t numElements = alpha.Dim();
    for ( size_t idxElement = 0; idxElement < numElements; ++idxElement )
    {
      if ( std::isnan(alpha(idxElement)) ) return true;
    }
    return false;
  }
}

template <unsigned int P, unsigned int C>
class KinFitAlgo
{
 public:
  typedef typename math::Matrix<P,P>::type MatrixPxP;
  typedef typename math::Matrix<P,C>::type MatrixPxC;
  typedef typename math::Matrix<C,P>::type MatrixCxP;
  typedef typename math::Matrix<C,C>::type MatrixCxC;
  typedef typename math::Vector<P>::type VectorP;
  typedef typename math::Vector<C>::type VectorC;

  KinFitAlgo(const edm::ParameterSet& cfg)
    : verbosity_(cfg.getParameter<int>("verbosity"))
  {}
  ~KinFitAlgo()
  {}

  KinFitSummary<P>
  operator()(const VectorP& alpha0, const MatrixPxP& V_alpha0, KinFitConstraintBase<P,C>& constraint, const VectorP& startPos,
             std::vector<KinFitSummary<P>>* fitHistory = nullptr)
  {
    int Vinv_alpha0_errorFlag = 0;
    MatrixPxP Vinv_alpha0 = V_alpha0.Inverse(Vinv_alpha0_errorFlag);
    if ( Vinv_alpha0_errorFlag != 0 ) 
    {
      printCovMatrix("V_alpha0", V_alpha0);
      throw cmsException("KinFitAlgo::operator()", __LINE__) 
        << "Failed to invert matrix V_alpha0 !!\n";
    }

    if ( fitHistory )
    {
      fitHistory->push_back(KinFitSummary<P>(-1, startPos, V_alpha0, -1., -1));
    }

    const int numParameters = P;
    const int numConstraints = C;

    KinFitSummary<P> bestfit;
    double min_chi2 = -1.;
    int status = -1;
    const int max_iterations = 25;
    bool hasConverged = false;
    int iteration = 0;
    VectorP alphaA = startPos;
    while ( !hasConverged && iteration < max_iterations )
    {
      if ( verbosity_ >= 1 )
      {
        std::cout << "iteration #" << iteration << ":\n";
        std::cout << "alphaA:\n";
        std::cout << alphaA << "\n";
      }

      constraint.set_alphaA(alphaA);
      MatrixCxP D = constraint.get_D();
      VectorC d = constraint.get_d();
      if ( verbosity_ >= 1 )
      {
        std::cout << "D (@alphaA):\n";
        std::cout << D << "\n";
        std::cout << "d (@alphaA):\n";
        std::cout << d << "\n";
        double d_mag = std::sqrt(ROOT::Math::Dot(d, constraint.get_d_metric()*d));
        std::cout << "|d (@alphaA)| = " << d_mag << "\n";
      }

      MatrixPxC DT = ROOT::Math::Transpose(D);

      MatrixCxC Vinv_D = D*V_alpha0*DT;
      int V_D_errorFlag = 0;
      MatrixCxC V_D = Vinv_D.Inverse(V_D_errorFlag);
      if ( V_D_errorFlag != 0 )
      {
        printCovMatrix("Vinv_D", Vinv_D);
        throw cmsException("KinFitAlgo::operator()", __LINE__) 
          << "Failed to invert matrix Vinv_D !!\n";
      }

      VectorP dalpha0 = alpha0 - alphaA;
      VectorC lambda = V_D*(D*dalpha0 + d);
      if ( verbosity_ >= 1 )
      {
        // CV: Compute "distance from satisfaction" according to formula
        //     given at the end of Section III in https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
        //     The "distance from satisfaction" corresponds to how much the fit paramaters need to get "pulled"
        //     in order to satisfy the constraint equations
        VectorC D_times_dalpha0 = D*dalpha0;
        VectorC DfS;
        for ( int idx = 0; idx < numConstraints; ++idx )
        {
          DfS(idx) = (D_times_dalpha0(idx) + d(idx))/std::sqrt(Vinv_D(idx,idx));
        }
        std::cout << "DfS:\n";
        std::cout << DfS << "\n";
      }

      //const double sfStepSize = +1.;
      double sfStepSize = std::min(1., 2.*iteration/max_iterations);
      VectorP alpha = alpha0 - sfStepSize*V_alpha0*DT*lambda;

      MatrixPxP V_alpha = V_alpha0 - V_alpha0*DT*V_D*D*V_alpha0;
      if ( verbosity_ >= 1 )
      {
        printCovMatrix("V_alpha", V_alpha);
      }

      VectorP dalpha = alpha - alphaA;
      double dalpha_mag = std::sqrt(ROOT::Math::Dot(dalpha, dalpha));
      if ( verbosity_ >= 1 )
      {
        std::cout << "dalpha:\n";
        std::cout << dalpha << "\n";
        std::cout << "|alpha - alphaA| = " << dalpha_mag << "\n";
      }

      VectorP alpha_minus_alpha0 = alpha - alpha0;
      if ( verbosity_ >= 1 )
      {
        VectorP pulls;
        for ( int idx = 0; idx < numParameters; ++idx )
        {      
          pulls(idx) = alpha_minus_alpha0(idx)/std::sqrt(V_alpha0(idx,idx));
        }
        std::cout << "pulls:\n";
        std::cout << pulls << "\n";
      }

      double chi2 = ROOT::Math::Dot(alpha_minus_alpha0, Vinv_alpha0*alpha_minus_alpha0) + ROOT::Math::Dot(lambda, D*dalpha + d);
      chi2 /= numParameters;
      if ( verbosity_ >= 1 )
      {
        std::cout << "chi^2/DoF = " << chi2 << "\n";
      }

      VectorC residuals = D*dalpha + d;
      double residuals_sum = 0.;
      for ( int idx = 0; idx < numConstraints; ++idx )
      {
        residuals_sum += std::fabs(residuals(idx));
      }
      if ( verbosity_ >= 1 )
      {
        std::cout << "residuals of constraint equations:\n";
        std::cout << residuals << "\n";
        std::cout << "sum_i |residual[i]| = " << residuals_sum << "\n";
      }

      bool alpha_isNaN = isNaN(alpha);
      constraint.set_alphaA(alpha);
      d = constraint.get_d();
      double d_mag = std::sqrt(ROOT::Math::Dot(d, constraint.get_d_metric()*d));
      if ( verbosity_ >= 1 )
      {
        std::cout << "isNaN(alpha) = " << alpha_isNaN << "\n";
        std::cout << "d (@alpha):\n";
        std::cout << d << "\n";
        std::cout << "|d (@alpha)| = " << d_mag << "\n";
      }
      if ( !alpha_isNaN && !constraint.get_errorFlag() && d_mag < 1.e-1 )
      {
        if ( status == -1 || chi2 < min_chi2 )
        {
          min_chi2 = chi2;
          status = 0;
          bestfit = KinFitSummary<P>(iteration, alpha, V_alpha, chi2, status);
        }
        if ( status == 0 && dalpha_mag < 1.e-3 )
        {
          status = 1;
          bestfit = KinFitSummary<P>(iteration, alpha, V_alpha, chi2, status);
          hasConverged = true;
        }
      }

      if ( fitHistory )
      {
        fitHistory->push_back(KinFitSummary<P>(iteration, alpha, V_alpha, chi2, status));
      }

      alphaA = alpha;

      ++iteration;
    }

    if ( bestfit.get_status() == 0 && bestfit.get_chi2() < 1.e+2 )
    {
      hasConverged = true;
    }
    //if ( !hasConverged )
    //{
    //  std::cerr << "WARNING: KinematicFit failed to converge !!" << std::endl;
    //}

    if ( verbosity_ >= 1 )
    {
      std::cout << "best fit (@KinFitAlgo):\n";
      std::cout << "iteration = " << bestfit.get_iteration() << " (max_iterations = " << max_iterations << ")\n";
      std::cout << "alpha:\n";
      std::cout << bestfit.get_alpha() << "\n";
      std::cout << "status = " << bestfit.get_status() << "\n";
      std::cout << "min(chi^2/DoF) = " << bestfit.get_chi2() << "\n";
    }

    return bestfit;
  }

 private:
  int verbosity_;
};

#endif // TauAnalysis_Entanglement_KinFitAlgo_h
