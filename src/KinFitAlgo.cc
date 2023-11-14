#include "TauAnalysis/Entanglement/interface/KinFitAlgo.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/format_vT.h"    // format_vint()
#include "TauAnalysis/Entanglement/interface/printMatrix.h"  // printMatrix()
#include "TauAnalysis/Entanglement/interface/printVector.h"  // printVector()

#include <cmath>                                             // std::isnan(), std::sqrt()

KinFitSummary::KinFitSummary()
  : iteration_(-1)
  , chi2_(-1.)
  , status_(-1)
{}
  
KinFitSummary::KinFitSummary(int iteration, const TVectorD& alpha, const TMatrixD& V_alpha, double chi2, int status)
  : iteration_(iteration)
  , alpha_(alpha)
  , V_alpha_(V_alpha)
  , chi2_(chi2)
  , status_(status)
{}

KinFitSummary::~KinFitSummary()
{}

KinFitSummary&
KinFitSummary::operator=(const KinFitSummary& fit)
{
  iteration_ = fit.iteration_;
  alpha_.ResizeTo(fit.alpha_.GetNrows());
  alpha_ = fit.alpha_;
  V_alpha_.ResizeTo(fit.V_alpha_.GetNrows(),fit.V_alpha_.GetNcols());
  V_alpha_ = fit.V_alpha_;
  chi2_ = fit.chi2_;
  status_ =  fit.status_;
  return *this;
}

int
KinFitSummary::get_iteration() const
{
  return iteration_;
}

const TVectorD&
KinFitSummary::get_alpha() const
{
  return alpha_;
}

const TMatrixD&
KinFitSummary::get_V_alpha() const
{
  return V_alpha_;
}

double
KinFitSummary::get_chi2() const
{
  return chi2_;
}

int
KinFitSummary::get_status() const
{
  return status_;
}

namespace
{
  bool
  isNaN(const TVectorD& alpha)
  {
    size_t numElements = alpha.GetNrows();
    for ( size_t idxElement = 0; idxElement < numElements; ++idxElement )
    {
      if ( std::isnan(alpha(idxElement)) ) return true;
    }
    return false;
  }
}

KinFitAlgo::KinFitAlgo(const edm::ParameterSet& cfg)
  : verbosity_(cfg.getParameter<int>("verbosity"))
{}

KinFitAlgo::~KinFitAlgo()
{}

KinFitSummary
KinFitAlgo::operator()(const TVectorD& alpha0, const TMatrixD& V_alpha0, KinFitConstraintBase& constraint, const TVectorD& startPos,
                       std::vector<KinFitSummary>* fitHistory)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<KinFitAlgo::operator()>:\n";
  }

  if ( V_alpha0.Determinant() == 0. )
  {
    printMatrix("V_alpha0", V_alpha0);
    throw cmsException("KinFitAlgo::operator()", __LINE__) 
      << "Failed to invert matrix V_alpha0 !!\n";
  }
  TMatrixD Vinv_alpha0(TMatrixD::kInverted, V_alpha0);

  if ( fitHistory )
  {
    fitHistory->push_back(KinFitSummary(-1, startPos, V_alpha0, -1., -1));
  }

  const unsigned int Np = constraint.get_Np();
  const unsigned int Nc_eq = constraint.get_Nc_eq();

  KinFitSummary bestfit;
  double min_chi2 = -1.;
  int status = -1;
  const int max_iterations = 25;
  bool hasConverged = false;
  int iteration = 0;
  TVectorD alphaA = startPos;
  while ( !hasConverged && iteration < max_iterations )
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "iteration #" << iteration << ":\n";
      printVector("alphaA", alphaA);
    }

    constraint.set_alphaA(alphaA);

    const TMatrixD D_eq = constraint.get_D_eq();
    const TVectorD d_eq = constraint.get_d_eq();
    if ( verbosity_ >= 3 )
    {
      printMatrix("D_eq (@alphaA)", D_eq);
      printVector("d_eq (@alphaA)", d_eq);
      double d_eq_mag = std::sqrt(d_eq*(constraint.get_d_eq_metric()*d_eq));
      std::cout << "|d_eq (@alphaA)| = " << d_eq_mag << "\n";
    }

    TMatrixD D_ineq;
    TVectorD d_ineq;
    if ( constraint.get_Nc_ineq() > 0 )
    {
      D_ineq.ResizeTo(constraint.get_Nc_ineq(),Np);
      D_ineq = constraint.get_D_ineq();
      d_ineq.ResizeTo(constraint.get_Nc_ineq());
      d_ineq = constraint.get_d_ineq();
      if ( verbosity_ >= 3 )
      {
        printMatrix("D_ineq (@alphaA)", D_ineq);
        printVector("d_ineq (@alphaA)", d_ineq);
      }
    }

    // CV: different cases to handle inequality constraints,
    //     cf. https://machinelearningmastery.com/lagrange-multiplier-approach-with-inequality-constraints/
    double bestcase_chi2 = -1.;
    int bestcase_Nc_ineq = -1;
    TVectorD bestcase_alpha(Np);
    TMatrixD bestcase_V_alpha(Np,Np);
    double bestcase_d_eq_mag = -1.;
    double bestcase_dalpha_mag = -1.;
    bool isFirst = true;
    unsigned int numCases = (1 << constraint.get_Nc_ineq());
    for ( unsigned int idxCase = 0; idxCase < numCases; ++idxCase )
    {
      std::vector<int> rowIndices;
      for ( unsigned int idxConstraint_ineq = 0; idxConstraint_ineq < constraint.get_Nc_ineq(); ++idxConstraint_ineq )
      {
        if ( idxCase & (1 << idxConstraint_ineq) ) rowIndices.push_back(idxConstraint_ineq);
      }
      const unsigned int Nc_ineq = rowIndices.size();

      const unsigned int Nc = Nc_eq + Nc_ineq;
      if ( verbosity_ >= 1 )
      {
        std::cout << "case #" << idxCase << ": Nc = " << Nc << "\n";
        std::cout << "(inequality constraints = " << format_vint(rowIndices) << ")\n";
      }

      TMatrixD D(Nc, constraint.get_Np());
      D.SetSub(0, 0, D_eq);
      TVectorD d(Nc);
      d.SetSub(0, d_eq);
      for ( unsigned int idxRow = 0; idxRow < rowIndices.size(); ++idxRow )
      {
        for ( unsigned int idxColumn = 0; idxColumn < constraint.get_Np(); ++idxColumn )
        {
          D(constraint.get_Nc_eq()+idxRow,idxColumn) = D_ineq(rowIndices[idxRow],idxColumn);
        }
        d(constraint.get_Nc_eq()+idxRow) = d_ineq(rowIndices[idxRow]);
      }
      if ( verbosity_ >= 2 )
      {
        printMatrix("D (@alphaA)", D);
        printVector("d (@alphaA)", d);
      }

      TMatrixD DT(TMatrixD::kTransposed, D);

      TMatrixD Vinv_D = D*V_alpha0*DT;
      if ( Vinv_D.Determinant() == 0. )
      {
        if ( verbosity_ >= 0 )
        {
          printMatrix("Vinv_D", Vinv_D);
        }
        //throw cmsException("KinFitAlgo::operator()", __LINE__) 
        //  << "Failed to invert matrix Vinv_D !!\n";
        std::cerr << "WARNING: Failed to invert matrix Vinv_D !!\n";
        continue;
      }
      TMatrixD V_D(TMatrixD::kInverted, Vinv_D);

      TVectorD dalpha0 = alpha0 - alphaA;
      TVectorD lambda = V_D*(D*dalpha0 + d);
      if ( verbosity_ >= 2 )
      {
        // CV: Compute "distance from satisfaction" according to formula
        //     given at the end of Section III in https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
        //     The "distance from satisfaction" corresponds to how much the fit paramaters need to get "pulled"
        //     in order to satisfy the constraint equations
        TVectorD D_times_dalpha0 = D*dalpha0;
        TVectorD DfS(Nc);
        for ( unsigned int idx = 0; idx < Nc; ++idx )
        {
          DfS(idx) = (D_times_dalpha0(idx) + d(idx))/std::sqrt(Vinv_D(idx,idx));
        }
        printVector("DfS", DfS);
      }

      //const double sfStepSize = +1.;
      double sfStepSize = std::min(1., 2.*iteration/max_iterations);
      TVectorD alpha = alpha0 - sfStepSize*V_alpha0*DT*lambda;

      TMatrixD V_alpha = V_alpha0 - V_alpha0*DT*V_D*D*V_alpha0;
      if ( verbosity_ >= 2 )
      {
        printMatrix("V_alpha", V_alpha);
      }

      TVectorD dalpha = alpha - alphaA;
      double dalpha_mag = std::sqrt(dalpha*dalpha);
      if ( verbosity_ >= 1 )
      {
        printVector("dalpha", dalpha);
        std::cout << "|alpha - alphaA| = " << dalpha_mag << "\n";
      }

      TVectorD alpha_minus_alpha0 = alpha - alpha0;
      if ( verbosity_ >= 2 )
      {
        TVectorD pulls(Np);
        for ( unsigned int idx = 0; idx < Np; ++idx )
        {
          pulls(idx) = alpha_minus_alpha0(idx)/std::sqrt(V_alpha0(idx,idx));
        }
        printVector("pulls", pulls);
      }

      double chi2 = alpha_minus_alpha0*(Vinv_alpha0*alpha_minus_alpha0) + lambda*(D*dalpha + d);
      chi2 /= Np;
      if ( verbosity_ >= 1 )
      {
        std::cout << "chi^2/DoF = " << chi2 << "\n";
      }

      TVectorD residuals = D*dalpha + d;
      double residuals_sum = 0.;
      for ( unsigned int idx = 0; idx < Nc_eq; ++idx )
      {
        residuals_sum += std::fabs(residuals(idx));
      }
      if ( verbosity_ >= 2 )
      {
        printVector("residuals of constraint equations", residuals);
        std::cout << "sum_i |residual[i]| = " << residuals_sum << "\n";
      }

      bool alpha_isNaN = isNaN(alpha);
      constraint.set_alphaA(alpha);
      const TVectorD dA_eq = constraint.get_d_eq();
      double dA_eq_mag = std::sqrt(dA_eq*(constraint.get_d_eq_metric()*dA_eq));
      if ( verbosity_ >= 1 )
      {
        std::cout << "isNaN(alpha) = " << alpha_isNaN << "\n";
        printVector("d_eq (@alpha)", dA_eq);
        std::cout << "|d_eq (@alpha)| = " << dA_eq_mag << "\n";
      }
      bool skip = false;
      if ( constraint.get_Nc_ineq() > 0 )
      {
        const TVectorD& dA_ineq = constraint.get_d_ineq();
        for ( unsigned int idxRow = 0; idxRow < rowIndices.size(); ++idxRow )
        {
          if ( !(dA_ineq(rowIndices[idxRow]) <= 0.) )
          {
            if ( verbosity_ >= 1 )
            {
              std::cout << "--> skipping this solution, because inequality constraint #" << rowIndices[idxRow] << " is violated !!\n";
              std::cout << "   (H = " << dA_ineq(rowIndices[idxRow]) << ") !!\n";
            }
            skip = true;
          }
        }
      }
      // CV: do not require inequality constraints to be satisfied in case all inequality constraints are "active",
      //     cf. https://machinelearningmastery.com/lagrange-multiplier-approach-with-inequality-constraints/
      //     in order to ensure that there is at least one case (of handling the inequality constraints)
      //     that yields a solution for each iteration of the kinematic fit,
      //     so that the kinematic fit does not get "stuck"
      if ( Nc_ineq == constraint.get_Nc_ineq() )
      {
        skip = false;
      }
      if ( !alpha_isNaN && !skip && ((int)Nc_ineq < bestcase_Nc_ineq || ((int)Nc_ineq == bestcase_Nc_ineq && chi2 < bestcase_chi2) || isFirst) )
      {
std::cout << "iteration #" << iteration << ": setting bestcase...\n";
std::cout << "chi2 = " << chi2 << ", Nc_ineq = " << Nc_ineq << ", dA_eq_mag = " << dA_eq_mag << ", dalpha_mag = " << dalpha_mag << "\n";
        bestcase_chi2 = chi2;
        bestcase_Nc_ineq = Nc_ineq;
        bestcase_alpha = alpha;
        bestcase_V_alpha = V_alpha;
        bestcase_d_eq_mag = dA_eq_mag;
        bestcase_dalpha_mag = dalpha_mag;
        isFirst = false;
      }
    }
std::cout << "iteration #" << iteration << ": bestcase_d_eq_mag = " << bestcase_d_eq_mag << "\n";
    if ( bestcase_d_eq_mag >= 0. && bestcase_d_eq_mag < 1.e-2 )
    {
      if ( status == -1 || bestcase_chi2 < min_chi2 )
      {
std::cout << "iteration #" << iteration << ": setting bestfit...\n";
std::cout << "bestcase_chi2 = " << bestcase_chi2 << "\n";
        min_chi2 = bestcase_chi2;
        status = 0;
        bestfit = KinFitSummary(iteration, bestcase_alpha, bestcase_V_alpha, bestcase_chi2, status);
      }
      if ( status == 0 && bestcase_dalpha_mag >= 0. && bestcase_dalpha_mag < 1.e-3 )
      {
std::cout << "iteration #" << iteration << ": updating bestfit...\n";
std::cout << "bestcase_chi2 = " << bestcase_chi2 << "\n";
        status = 1;
        bestfit = KinFitSummary(iteration, bestcase_alpha, bestcase_V_alpha, bestcase_chi2, status);
        hasConverged = true;
      }
    }

    if ( fitHistory )
    {
      fitHistory->push_back(KinFitSummary(iteration, bestcase_alpha, bestcase_V_alpha, bestcase_chi2, status));
    }

    alphaA = bestcase_alpha;

    ++iteration;
  }

  if ( bestfit.get_status() == 0 && bestfit.get_chi2() < 1.e+2 )
  {
    hasConverged = true;
  }
  if ( !hasConverged && verbosity_ >= 0 )
  {
    std::cerr << "WARNING: KinematicFit failed to converge !!\n";
  }

  if ( verbosity_ >= 1 )
  {
    std::cout << "best fit (@KinFitAlgo):\n";
    std::cout << "iteration = " << bestfit.get_iteration() << " (max_iterations = " << max_iterations << ")\n";
    if ( bestfit.get_status() != -1 )
    {
      printVector("alpha", bestfit.get_alpha());
    }
    std::cout << "status = " << bestfit.get_status() << "\n";
    std::cout << "min(chi^2/DoF) = " << bestfit.get_chi2() << "\n";
  }

  return bestfit;
}
