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
  int bestfit_status = -1;
  double bestfit_chi2 = -1.;
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
    int bestcase_Nc_ineq = -1;
    double bestcase_residuals_eq_mag = -1.;
    double bestcase_dalpha_mag = -1.;
    double bestcase_chi2 = -1.;
    TVectorD bestcase_alpha(Np);
    TMatrixD bestcase_V_alpha(Np,Np);
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

      TMatrixD V_D = D*V_alpha0*DT;
      if ( V_D.Determinant() == 0. )
      {
        if ( verbosity_ >= 0 )
        {
          printMatrix("V_D", V_D);
        }
        //throw cmsException("KinFitAlgo::operator()", __LINE__) 
        //  << "Failed to invert matrix V_D !!\n";
        std::cerr << "WARNING: Failed to invert matrix V_D !!\n";
        continue;
      }
      TMatrixD Vinv_D(TMatrixD::kInverted, V_D);

      TVectorD dalpha0 = alpha0 - alphaA;
      TVectorD lambda = Vinv_D*(D*dalpha0 + d);
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
          DfS(idx) = (D_times_dalpha0(idx) + d(idx))/std::sqrt(V_D(idx,idx));
        }
        printVector("DfS", DfS);
      }

      // CV: take direction for changing parameter vector alpha from system of linear equations,
      //     but optimize length of change vector separately by scanning stepsize between 0 and 2 times 
      //     the stepsize obtained from the system of linear equations.
      //     Take as optimal stepsize the stepsize which best satisfies all constraint equations, regardless of chi^2.
      double best_sfStepSize = -1.;
      double best_residuals_eq_mag = -1.;
      for ( double sfStepSize = 0.; sfStepSize <= 2.5; sfStepSize += 0.1 )
      {
        TVectorD alpha = alpha0 - sfStepSize*V_alpha0*DT*lambda;
        constraint.set_alphaA(alpha);
        const TVectorD residuals_eq = constraint.get_d_eq();
        double residuals_eq_mag = std::sqrt(residuals_eq*(constraint.get_d_eq_metric()*residuals_eq));
        if ( verbosity_ >= 1 )
        {
          std::cout << "sfStepSize = " << sfStepSize << ": |residuals| = " << residuals_eq_mag << "\n";
        }
        if ( best_sfStepSize < 0 || residuals_eq_mag < best_residuals_eq_mag )
        {
          best_sfStepSize = sfStepSize;
          best_residuals_eq_mag = residuals_eq_mag;
        }
      }
//std::cout << "iteration #" << iteration << ": best_sfStepSize = " << best_sfStepSize << ", best_residuals_eq_mag = " << best_residuals_eq_mag << "\n";
      TVectorD alpha = alpha0 - best_sfStepSize*V_alpha0*DT*lambda;

      TMatrixD V_alpha = V_alpha0 - V_alpha0*DT*Vinv_D*D*V_alpha0;
      if ( verbosity_ >= 2 )
      {
        printMatrix("V_alpha", V_alpha);
      }

      TVectorD dalpha = alpha - alphaA;
      double dalpha_mag = std::sqrt(dalpha*(Vinv_alpha0*dalpha));
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
          pulls(idx) = alpha_minus_alpha0(idx)/std::sqrt(V_alpha0(idx,idx) - V_alpha(idx,idx));
        }
        printVector("pulls", pulls);
      }

      double chi2 = alpha_minus_alpha0*(Vinv_alpha0*alpha_minus_alpha0) + 2.*lambda*(D*dalpha + d);
      chi2 /= Np;
      if ( verbosity_ >= 1 )
      {
        std::cout << "chi^2/DoF = " << chi2 << "\n";
      }

      bool alpha_isNaN = isNaN(alpha);
      if ( verbosity_ >= 1 )
      {
        std::cout << "isNaN(alpha) = " << alpha_isNaN << "\n";
      }

      constraint.set_alphaA(alpha);
      TVectorD residuals_eq = constraint.get_d_eq();
      double residuals_eq_mag = std::sqrt(residuals_eq*(constraint.get_d_eq_metric()*residuals_eq));
      if ( verbosity_ >= 1 )
      {
        printVector("residuals of constraint equations", residuals_eq);
        std::cout << "|residuals| = " << residuals_eq_mag << "\n";
      }

      bool skip = false;
      if ( constraint.get_Nc_ineq() > 0 )
      {
        const TVectorD& residuals_ineq = constraint.get_d_ineq();
        for ( unsigned int idxRow = 0; idxRow < rowIndices.size(); ++idxRow )
        {
          if ( !(residuals_ineq(rowIndices[idxRow]) <= 0.) )
          {
            if ( verbosity_ >= 1 )
            {
              std::cout << "--> skipping this solution, because inequality constraint #" << rowIndices[idxRow] << " is violated !!\n";
              std::cout << "   (H = " << residuals_ineq(rowIndices[idxRow]) << ") !!\n";
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
      if ( !alpha_isNaN && !skip && ((int)Nc_ineq < bestcase_Nc_ineq || ((int)Nc_ineq == bestcase_Nc_ineq && residuals_eq_mag < bestcase_residuals_eq_mag) || isFirst) )
      {
//std::cout << "iteration #" << iteration << ": setting bestcase...\n";
//std::cout << "Nc_ineq = " << Nc_ineq << ", residuals_eq_mag = " << residuals_eq_mag << ", dalpha_mag = " << dalpha_mag << ", chi2 = " << chi2 << "\n";
        bestcase_Nc_ineq = Nc_ineq;
        bestcase_residuals_eq_mag = residuals_eq_mag;
        bestcase_dalpha_mag = dalpha_mag;
        bestcase_chi2 = chi2;
        bestcase_alpha = alpha;
        bestcase_V_alpha = V_alpha;
        isFirst = false;
      }
    }
//std::cout << "iteration #" << iteration << ":\n";
//std::cout << "bestcase_Nc_ineq = " << bestcase_Nc_ineq << ", bestcase_residuals_eq_mag = " << bestcase_residuals_eq_mag << ", bestcase_dalpha_mag = " << bestcase_dalpha_mag << ", bestcase_chi2 = " << bestcase_chi2 << "\n";
//printVector("bestcase_alpha", bestcase_alpha);

    int bestcase_status = -1;
    if ( bestcase_residuals_eq_mag >= 0. && bestcase_residuals_eq_mag < 1. )
    {
      bestcase_status = 0;
    }
    if ( bestcase_status == 0 && bestcase_dalpha_mag > 0. && bestcase_dalpha_mag < 1. )
    {
      bestcase_status = 1;
    }
    if ( bestcase_status == 1 && bestcase_residuals_eq_mag < 1.e-1 && bestcase_dalpha_mag < 1.e-1 )
    {
      bestcase_status = 2;
    }
//std::cout << "iteration #" << iteration << ": bestcase_status = " << bestcase_status << "\n";

    if ( bestcase_status > bestfit_status || (bestcase_status == bestfit_status && bestcase_chi2 < bestfit_chi2) )
    {
//if ( bestfit_status == -1 )
//{
//  std::cout << "iteration #" << iteration << ": setting bestfit...\n";
//}
//else
//{
//  std::cout << "iteration #" << iteration << ": updating bestfit...\n";
//}
      bestfit = KinFitSummary(iteration, bestcase_alpha, bestcase_V_alpha, bestcase_chi2, bestcase_status);
      bestfit_status = bestcase_status;
      bestfit_chi2 = bestcase_chi2;
      if ( bestfit_status == 2 )
      {
        hasConverged = true;
      }
    }

    if ( fitHistory )
    {
      fitHistory->push_back(KinFitSummary(iteration, bestcase_alpha, bestcase_V_alpha, bestcase_chi2, bestcase_status));
    }

    alphaA = bestcase_alpha;
//std::cout << "iteration #" << iteration << ":\n";
//printVector("alphaA", alphaA);

    ++iteration;
  }

  if ( !(bestfit.get_status() >= 0) && verbosity_ >= 0 )
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
