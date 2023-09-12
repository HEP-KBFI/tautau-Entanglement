#include "TauAnalysis/Entanglement/interface/KinematicFit.h"

#include "DataFormats/Math/interface/deltaPhi.h"                          // deltaPhi()
#include "DataFormats/Math/interface/Matrix.h"                            // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                            // math::Vector
#include "DataFormats/TauReco/interface/PFTau.h"                          // reco::PFTau::kOneProng0PiZero

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/comp_nuPz.h"                 // build_nuP4(), comp_nuPz()
#include "TauAnalysis/Entanglement/interface/constants.h"                 // kLHC, kSuperKEKB
#include "TauAnalysis/Entanglement/interface/KinFitAlgo.h"                // KinFitAlgo<>, KinFitSummary<>
#include "TauAnalysis/Entanglement/interface/KinFitConstraint.h"          // KinFitConstraint<>
#include "TauAnalysis/Entanglement/interface/KinFitParameters.h"          // kinFit::numParameters, math::MatrixPxP
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h"         // math::Matrix*, math::Vector*
#include "TauAnalysis/Entanglement/interface/printCovMatrix.h"            // printCovMatrix()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"        // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/square.h"                    // square()

#include "Math/Functions.h"                                               // ROOT::Math::Dot(), ROOT::Math::Similarity(), ROOT::Math::Transpose() 
#include "TMath.h"                                                        // TMath::Pi() 

#include <cmath>                                                          // std::abs(), std::exp(), std::fabs(), std::isnan(), std::sin(), std::sqrt()
#include <iostream>                                                       // std::cout
#include <string>                                                         // std::string

using namespace kinFit;
using namespace math;

const int numConstraints_LHC = 9;
const int numConstraints_SuperKEKB = 8;

KinematicFit::KinematicFit(const edm::ParameterSet& cfg)
  : polarimetricVector_(cfg)
  , skip_(cfg.getParameter<bool>("skip"))
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  std::string collider = cfg.getParameter<std::string>("collider");
  if      ( collider == "LHC"       ) collider_ = kLHC;
  else if ( collider == "SuperKEKB" ) collider_ = kSuperKEKB;
  else throw cmsException("KinematicFit", __LINE__)
    << "Invalid Configuration parameter 'collider' = " << collider << " !!\n";
}

KinematicFit::~KinematicFit()
{}

namespace kinFit
{
  math::Matrix3x3
  build_nuCov(const math::Matrix2x2& cov2x2, double nu_dPzdPx, double nu_dPzdPy)
  {
    math::Matrix3x3 cov3x3;
    cov3x3.Place_at(cov2x2, 0, 0);
    cov3x3(0,2) = cov2x2(0,0)*nu_dPzdPx + cov2x2(1,0)*nu_dPzdPy;
    cov3x3(1,2) = cov2x2(0,1)*nu_dPzdPx + cov2x2(1,1)*nu_dPzdPy;
    cov3x3(2,0) = cov3x3(0,2);
    cov3x3(2,1) = cov3x3(1,2);
    cov3x3(2,2) = cov2x2(0,0)*square(nu_dPzdPx) + cov2x2(0,1)*nu_dPzdPx*nu_dPzdPy + cov2x2(1,0)*nu_dPzdPy*nu_dPzdPx + cov2x2(1,1)*square(nu_dPzdPy);
    return cov3x3;
  }
}

KinematicEvent
KinematicFit::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<KinematicFit::operator()>:\n"; 
  }

  if ( skip_ )
  {
    return kineEvt;
  }

  KinematicEvent kineEvt_kinFit = kineEvt;

  double bestfit_chi2 = -1.;
  int bestfit_status = -1;
  for ( int idxSignTauPlus = 0; idxSignTauPlus <= 1; ++idxSignTauPlus )
  {
    double signTauPlus = 0.;
    if      ( idxSignTauPlus == 0 ) signTauPlus = -1.;
    else if ( idxSignTauPlus == 1 ) signTauPlus = +1.;
    else assert(0);
    for ( int idxSignTauMinus = 0; idxSignTauMinus <= 1; ++idxSignTauMinus )
    {
      double signTauMinus = 0.;
      if      ( idxSignTauMinus == 0 ) signTauMinus = -1.;
      else if ( idxSignTauMinus == 1 ) signTauMinus = +1.;
      else assert(0);

      if ( verbosity_ >= 1 )
      {
        std::cout << "sign(tau+) = " << signTauPlus << ", sign(tau-) = " << signTauMinus << "\n";
      }

      const reco::Candidate::Point& startPos_pv = kineEvt.pv();
      const Matrix3x3& pvCov = kineEvt.pvCov();
      const reco::Candidate::LorentzVector& startPos_visTauPlusP4 = kineEvt.visTauPlusP4();
      const reco::Candidate::LorentzVector& startPos_nuTauPlusP4 = kineEvt.nuTauPlusP4();
      Matrix2x2 nuTauPlusCov = kineEvt.nuTauPlusCov().Sub<Matrix2x2>(0,0);
      const reco::Candidate::Point& startPos_svTauPlus = kineEvt.svTauPlus();
      const Matrix3x3& svTauPlusCov = kineEvt.svTauPlusCov();
      const reco::Candidate::LorentzVector& startPos_visTauMinusP4 = kineEvt.visTauMinusP4();
      const reco::Candidate::LorentzVector& startPos_nuTauMinusP4 = kineEvt.nuTauMinusP4();
      Matrix2x2 nuTauMinusCov = kineEvt.nuTauMinusCov().Sub<Matrix2x2>(0,0);
      const reco::Candidate::Point& startPos_svTauMinus = kineEvt.svTauMinus();
      const Matrix3x3& svTauMinusCov = kineEvt.svTauMinusCov();
      const reco::Candidate::LorentzVector& startPos_recoilP4 = kineEvt.recoilP4();
      const Matrix4x4& recoilCov = kineEvt.recoilCov();

      // CV: Check that chosen signs for tau+ and tau- reproduce neutrino Pz of start position; skip sign combination if not.
      //     Skipping these sign combinations saves computing time and avoids running into unphysical solutions.
      double nuTauPlusPx = startPos_nuTauPlusP4.px();
      double nuTauPlusPy = startPos_nuTauPlusP4.py();
      double nuTauPlus_dPzdPx, nuTauPlus_dPzdPy;
      double nuTauPlusPz = comp_nuPz(startPos_visTauPlusP4, nuTauPlusPx, nuTauPlusPy, signTauPlus, nuTauPlus_dPzdPx, nuTauPlus_dPzdPy, verbosity_);
      double nuTauMinusPx = startPos_nuTauMinusP4.px();
      double nuTauMinusPy = startPos_nuTauMinusP4.py();
      double nuTauMinus_dPzdPx, nuTauMinus_dPzdPy;
      double nuTauMinusPz = comp_nuPz(startPos_visTauMinusP4, nuTauMinusPx, nuTauMinusPy, signTauMinus, nuTauMinus_dPzdPx, nuTauMinus_dPzdPy, verbosity_);
      if ( std::fabs(nuTauPlusPz  - startPos_nuTauPlusP4.pz())  > std::max(1., 0.10*startPos_nuTauPlusP4.pz())  ||
           std::fabs(nuTauMinusPz - startPos_nuTauMinusP4.pz()) > std::max(1., 0.10*startPos_nuTauMinusP4.pz()) )
      {
        if ( verbosity_ >= 1 )
        {
          std::cout << "--> skipping this sign combination, because it does not reproduce the neutrino Pz of the start position !!\n";
        }
        continue;
      }

      VectorP alpha0;
      alpha0( 0) = startPos_pv.x();
      alpha0( 1) = startPos_pv.y();
      alpha0( 2) = startPos_pv.z();
      alpha0( 3) = startPos_nuTauPlusP4.px();
      alpha0( 4) = startPos_nuTauPlusP4.py();
      alpha0( 5) = startPos_svTauPlus.x();
      alpha0( 6) = startPos_svTauPlus.y();
      alpha0( 7) = startPos_svTauPlus.z();
      alpha0( 8) = startPos_nuTauMinusP4.px();
      alpha0( 9) = startPos_nuTauMinusP4.py();
      alpha0(10) = startPos_svTauMinus.x();
      alpha0(11) = startPos_svTauMinus.y();
      alpha0(12) = startPos_svTauMinus.z();
      alpha0(13) = startPos_recoilP4.px();
      alpha0(14) = startPos_recoilP4.py();
      alpha0(15) = startPos_recoilP4.pz();
      alpha0(16) = startPos_recoilP4.energy();
      if ( verbosity_ >= 1 )
      {
        std::cout << "alpha0:\n";
        std::cout << alpha0 << "\n";
      }

      // CV: build covariance matrix of all measured parameters;
      //     for the syntax of "embedding" a small covariance matrix into a larger one, 
      //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
      MatrixPxP V_alpha0;
      V_alpha0.Place_at(pvCov        ,  0,  0);
      V_alpha0.Place_at(nuTauPlusCov ,  3,  3);
      V_alpha0.Place_at(svTauPlusCov ,  5,  5);
      V_alpha0.Place_at(nuTauMinusCov,  8,  8);
      V_alpha0.Place_at(svTauMinusCov, 10, 10);
      V_alpha0.Place_at(recoilCov    , 13, 13);
      if ( verbosity_ >= 1 )
      {
        printCovMatrix("V_alpha0", V_alpha0);
      }

      VectorP alphaA = alpha0;
      if ( verbosity_ >= 1 )
      {
        std::cout << "alphaA:\n";
        std::cout << alphaA << "\n";
      }

      // CV: compute solution to minimization problem
      //     using formulas given in Section 3 of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
      KinFitSummary<P> fitResult;
      if ( collider_ == kLHC ) 
      {
        const int C = numConstraints_LHC;
        KinFitConstraint<P,C> constraint(collider_, kineEvt, signTauPlus, signTauMinus, verbosity_);
        edm::ParameterSet cfg_kinFit;
        cfg_kinFit.addParameter<int>("verbosity", -1);
        KinFitAlgo<P,C> kinFit(cfg_kinFit);
        fitResult = kinFit(alpha0, V_alpha0, constraint, alphaA);
      }
      else if ( collider_ == kSuperKEKB )
      {
        const int C = numConstraints_SuperKEKB;
        KinFitConstraint<P,C> constraint(collider_, kineEvt, signTauPlus, signTauMinus, verbosity_);
        edm::ParameterSet cfg_kinFit;
        cfg_kinFit.addParameter<int>("verbosity", -1);
        KinFitAlgo<P,C> kinFit(cfg_kinFit);
        fitResult = kinFit(alpha0, V_alpha0, constraint, alphaA);
      }
      else assert(0);

      if (  fitResult.get_status() >  bestfit_status                                         ||
           (fitResult.get_status() == bestfit_status && fitResult.get_chi2() < bestfit_chi2) )
       {
        VectorP bestfit_alpha = fitResult.get_alpha();
        MatrixPxP bestfit_V_alpha = fitResult.get_V_alpha();

        // CV: store results of kinematic fit in KinematicEvent class;
        //     for the syntax of retrieving a small covariance matrix from a larger one, 
        //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
        kineEvt_kinFit.pv_ = reco::Candidate::Point(bestfit_alpha(0), bestfit_alpha(1), bestfit_alpha(2));
        kineEvt_kinFit.pvCov_ = bestfit_V_alpha.Sub<Matrix3x3>(0,0);
        kineEvt_kinFit.recoilP4_ = reco::Candidate::LorentzVector(bestfit_alpha(13), bestfit_alpha(14), bestfit_alpha(15), bestfit_alpha(16));
        kineEvt_kinFit.recoilCov_ = bestfit_V_alpha.Sub<Matrix4x4>(13,13);
        double bestfit_nuTauPlusPx = bestfit_alpha(3);
        double bestfit_nuTauPlusPy = bestfit_alpha(4);
        double bestfit_nuTauPlus_dPzdPx, bestfit_nuTauPlus_dPzdPy;
        double bestfit_nuTauPlusPz = comp_nuPz(
                 kineEvt.visTauPlusP4_, bestfit_nuTauPlusPx, bestfit_nuTauPlusPy, signTauPlus,
                 bestfit_nuTauPlus_dPzdPx, bestfit_nuTauPlus_dPzdPy,
                 verbosity_);
        kineEvt_kinFit.nuTauPlusP4_ = build_nuP4(bestfit_nuTauPlusPx, bestfit_nuTauPlusPy, bestfit_nuTauPlusPz);
        kineEvt_kinFit.nuTauPlusP4_isValid_ = true;
        kineEvt_kinFit.nuTauPlusCov_ = build_nuCov(bestfit_V_alpha.Sub<Matrix2x2>(3,3), bestfit_nuTauPlus_dPzdPx, bestfit_nuTauPlus_dPzdPy);
        kineEvt_kinFit.tauPlusP4_ = kineEvt.visTauPlusP4_ + kineEvt.nuTauPlusP4_;
        kineEvt_kinFit.tauPlusP4_isValid_ = true;
        kineEvt_kinFit.svTauPlus_ = reco::Candidate::Point(bestfit_alpha(5), bestfit_alpha(6), bestfit_alpha(7));
        kineEvt_kinFit.svTauPlusCov_ = bestfit_V_alpha.Sub<Matrix3x3>(5,5);
        double bestfit_nuTauMinusPx = bestfit_alpha(8);
        double bestfit_nuTauMinusPy = bestfit_alpha(9);
        double bestfit_nuTauMinus_dPzdPx, bestfit_nuTauMinus_dPzdPy;
        double bestfit_nuTauMinusPz = comp_nuPz(
                 kineEvt.visTauMinusP4_, bestfit_nuTauMinusPx, bestfit_nuTauMinusPy, signTauMinus,
                 bestfit_nuTauMinus_dPzdPx, bestfit_nuTauMinus_dPzdPy,
                 verbosity_);
        kineEvt_kinFit.nuTauMinusP4_ = build_nuP4(bestfit_nuTauMinusPx, bestfit_nuTauMinusPy, bestfit_nuTauMinusPz);
        kineEvt_kinFit.nuTauMinusP4_isValid_ = true;       
        kineEvt_kinFit.nuTauMinusCov_ = build_nuCov(bestfit_V_alpha.Sub<Matrix2x2>(8,8), bestfit_nuTauMinus_dPzdPx, bestfit_nuTauMinus_dPzdPy);
        kineEvt_kinFit.tauMinusP4_ = kineEvt.visTauMinusP4_ + kineEvt.nuTauMinusP4_;
        kineEvt_kinFit.tauMinusP4_isValid_ = true;
        kineEvt_kinFit.svTauMinus_ = reco::Candidate::Point(bestfit_alpha(10), bestfit_alpha(11), bestfit_alpha(12));
        kineEvt_kinFit.svTauMinusCov_ = bestfit_V_alpha.Sub<Matrix3x3>(10,10);
        
        kineEvt_kinFit.kinFitCov_ = bestfit_V_alpha;
        kineEvt_kinFit.kinFitStatus_ = fitResult.get_status();
        kineEvt_kinFit.kinFitChi2_ = fitResult.get_chi2();
        kineEvt_kinFit.kinFit_isValid_ = true;
 
        bestfit_status = fitResult.get_status();
        bestfit_chi2 = fitResult.get_chi2();
      }
    }
  }
  
  if ( verbosity_ >= 1 )
  {
    std::cout << "best fit (among all tau+, tau- sign combinations):\n";
    std::cout << "status = " << bestfit_status << "\n";
    std::cout << "min(chi^2/DoF) = " << bestfit_chi2 << "\n";
  }

  kineEvt_kinFit.recoilP4_ = kineEvt_kinFit.tauPlusP4_ + kineEvt_kinFit.tauMinusP4_;
  if ( kineEvt_kinFit.kinFit_isValid_ )
  {
    reco::Candidate::Vector hPlus = polarimetricVector_(kineEvt_kinFit, pol::kTauPlus);
    kineEvt_kinFit.hPlus_ = hPlus;
    kineEvt_kinFit.hPlus_isValid_ = true;
    reco::Candidate::Vector hMinus = polarimetricVector_(kineEvt_kinFit, pol::kTauMinus);
    kineEvt_kinFit.hMinus_ = hMinus;
    kineEvt_kinFit.hMinus_isValid_ = true;
  }

  return kineEvt_kinFit;
}

