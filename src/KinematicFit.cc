#include "TauAnalysis/Entanglement/interface/KinematicFit.h"

#include "DataFormats/Math/interface/deltaPhi.h"                          // deltaPhi()
#include "DataFormats/Math/interface/Matrix.h"                            // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                            // math::Vector
#include "DataFormats/TauReco/interface/PFTau.h"                          // reco::PFTau::kOneProng0PiZero

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/comp_cosThetaGJ.h"           // comp_cosThetaGJ(), comp_cosThetaGJ_solution()
#include "TauAnalysis/Entanglement/interface/comp_mag.h"                  // comp_mag()
#include "TauAnalysis/Entanglement/interface/constants.h"                 // ct, kLHC, kSuperKEKB, mHiggs, mTau
#include "TauAnalysis/Entanglement/interface/cube.h"                      // cube()
#include "TauAnalysis/Entanglement/interface/fixMass.h"                   // fixHiggsMass(), fixTauMass()
#include "TauAnalysis/Entanglement/interface/getCov_hf.h"                 // getCov_hf()
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h"         // math::Matrix*, math::Vector*
#include "TauAnalysis/Entanglement/interface/printCovMatrix.h"            // printCovMatrix()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"        // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printPoint.h"                // printPoint()
#include "TauAnalysis/Entanglement/interface/square.h"                    // square()

#include "Math/Functions.h"                                               // ROOT::Math::Dot(), ROOT::Math::Similarity(), ROOT::Math::Transpose() 
#include "TMath.h"                                                        // TMath::Pi() 

#include <cmath>                                                          // std::abs(), std::exp(), std::fabs(), std::isnan(), std::sin(), std::sqrt()
#include <iostream>                                                       // std::cout
#include <string>                                                         // std::string

using namespace kinFit;
using namespace math;

KinematicFit::KinematicFit(const edm::ParameterSet& cfg)
  : polarimetricVector_(cfg)
  , applyLifetimeConstraint_(cfg.getParameter<bool>("applyLifetimeConstraint"))
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
  double
  comp_nuPz(const reco::Candidate::LorentzVector& visP4, double nuPx, double nuPy, double sign, 
            double& nu_dPzdPx, double& nu_dPzdPy,
            int verbosity = -1)
  {
    if ( verbosity >= 3 )
    {
      std::cout << "<comp_nuPz>:\n";
    }

    double visPx    = visP4.px();
    double visPy    = visP4.py();
    double visPz    = visP4.pz();
    double visE     = visP4.energy();
    double visMass  = ( visMass < mTau ) ? visP4.mass() : mTau;
    double visMass2 = square(visMass);
    double mTau2    = square(mTau); 
    double term1    = visMass2 + square(visPx) + square(visPy);
    double term2    = (mTau2 - visMass2)*visPz + 2.*(nuPx*visPx + nuPy*visPy)*visPz;
    double term3    = square(mTau2) + square(visMass2) - 4.*square(visPx*nuPy - nuPx*visPy) 
                     + mTau2*(-2.*visMass2 + 4.*nuPx*visPx + 4.*nuPy*visPy) 
                     - 4.*visMass2*(nuPx*(nuPx + visPx) + nuPy*(nuPy + visPy));
    if ( verbosity >= 3 )
    {
      std::cout << "term1 = " << term1 << "\n";
      std::cout << "term2 = " << term2 << "\n";
      std::cout << "term3 = " << term3 << "\n";
    }

    double nuPz = (1./(2.*term1))*(term2 + sign*visE*sqrt(std::max(0., term3)));
    if ( verbosity >= 3 )
    {
      std::cout << "nuPz = " << nuPz << "\n";
    }

    nu_dPzdPx = (1./term1)*(visPx*visPz
               + sign*(visE/std::max(1.e-2, std::sqrt(term3)))*(mTau2*visPx - visMass2*(2.*nuPx + visPx) + 2.*visPx*nuPy*visPy - 2.*nuPx*square(visPy)));
    nu_dPzdPy = (1./term1)*(visPy*visPz
               + sign*(visE/std::max(1.e-2, std::sqrt(term3)))*(mTau2*visPy - visMass2*(2.*nuPy + visPy) + 2.*visPy*nuPx*visPx - 2.*nuPy*square(visPx)));
    if ( verbosity >= 3 )
    {
      std::cout << "nu_dPzdPx = " << nu_dPzdPx << "\n";
      std::cout << "nu_dPzdPy = " << nu_dPzdPy << "\n";
    }

    return nuPz;
  }

  reco::Candidate::LorentzVector
  build_nuP4(double nuPx, double nuPy, double nuPz)
  {
    double nuE = std::sqrt(square(nuPx) + square(nuPy) + square(nuPz));
    reco::Candidate::LorentzVector nuP4(nuPx, nuPy, nuPz, nuE);
    return nuP4;
  }

  math::Matrix3x3
  build_nuCov(const math::Matrix2x2& cov2x2, double nu_dPzdPx, double nu_dPzdPy)
  {
    math::Matrix3x3 cov3x3;
    cov3x3.Place_at(cov2x2, 0, 0);
    return cov3x3;
  }

  bool
  isNaN(const KinematicEvent& kineEvt, int verbosity = -1)
  {
    // CV: check if any parameters of kinematic fit are "not-a-number" (NaN)
    if ( verbosity >= 1 )
    {
      std::cout << "<isNaN>:\n";
    }
    const reco::Candidate::Point& pv = kineEvt.pv();
    if ( std::isnan(pv.X()) || std::isnan(pv.Y()) || std::isnan(pv.Z()) )
    {
      if ( verbosity >= 1 )
      {
        std::cout << "PV contains NaNs.\n";
      }
      return true;
    }
    const reco::Candidate::LorentzVector& nuTauPlusP4 = kineEvt.nuTauPlusP4();
    if ( std::isnan(nuTauPlusP4.px()) || std::isnan(nuTauPlusP4.py()) || std::isnan(nuTauPlusP4.pz()) || std::isnan(nuTauPlusP4.energy()) )
    {
      if ( verbosity >= 1 )
      {
        std::cout << "neutrino from tau+ decay contains NaNs.\n";
      }
      return true;
    }
    const reco::Candidate::Point& svTauPlus = kineEvt.svTauPlus();
    if ( std::isnan(svTauPlus.X()) || std::isnan(svTauPlus.Y()) || std::isnan(svTauPlus.Z()) )
    { 
      if ( verbosity >= 1 )
      {
        std::cout << "SV(tau+) contains NaNs.\n";
      }
      return true;
    }
    const reco::Candidate::LorentzVector& nuTauMinusP4 = kineEvt.nuTauMinusP4();
    if ( std::isnan(nuTauMinusP4.px()) || std::isnan(nuTauMinusP4.py()) || std::isnan(nuTauMinusP4.pz()) || std::isnan(nuTauMinusP4.energy()) )
    { 
      if ( verbosity >= 1 )
      {
        std::cout << "neutrino from tau- decay contains NaNs.\n";
      }
      return true;
    }
    const reco::Candidate::Point& svTauMinus = kineEvt.svTauMinus();
    if ( std::isnan(svTauMinus.X()) || std::isnan(svTauMinus.Y()) || std::isnan(svTauMinus.Z()) )
    {
      if ( verbosity >= 1 )
      {
        std::cout << "SV(tau-) contains NaNs.\n";
      }
      return true;
    }
    const reco::Candidate::LorentzVector& recoilP4 = kineEvt.recoilP4();
    if ( std::isnan(recoilP4.px()) || std::isnan(recoilP4.py()) || std::isnan(recoilP4.pz()) || std::isnan(recoilP4.energy()) )
    {
      if ( verbosity >= 1 )
      {
        std::cout << "recoil contains NaNs.\n";
      }
      return true;
    }
    return false;
  }

  bool
  isPhysicalSolution(const KinematicEvent& kineEvt, int collider, int verbosity = -1)
  {
    // CV: check that constraint equations are satisfied
    //    (not just their linear approximation)
    if ( verbosity >= 1 )
    {
      std::cout << "<isPhysicalSolution>:\n";
    }
    const double max_dphi = (30./180.)*TMath::Pi();
    const double max_dtheta = (30./180.)*TMath::Pi();
    if ( collider != kSuperKEKB )
    {
      reco::Candidate::LorentzVector higgsP4 = kineEvt.tauPlusP4() + kineEvt.tauMinusP4();
      if ( std::fabs(higgsP4.mass() - mHiggs) > 10. )
      {
        if ( verbosity >= 1 )
        {
          std::cout << "fails Higgs mass constraint (mass = " << higgsP4.mass() << ").\n";
        }
        return false;
      }
    }
    if ( !(kineEvt.tauPlusP4().mass() >= 0. && kineEvt.tauPlusP4().mass() <= 5.) )
    { 
      if ( verbosity >= 1 )
      {
        std::cout << "fails tau+ mass constraint (mass = " << kineEvt.tauPlusP4().mass() << ").\n";
      }
      return false;
    }
    auto tauPlusD3 = kineEvt.svTauPlus() - kineEvt.pv();
    if ( std::fabs(reco::deltaPhi(kineEvt.tauPlusP4().phi(), tauPlusD3.phi())) > max_dphi   || 
         std::fabs(kineEvt.tauPlusP4().theta() - tauPlusD3.theta())            > max_dtheta )
    {
      if ( verbosity >= 1 )
      {
        std::cout << "fails parallelism constraint for tau+.\n";
        std::cout << " phi of tau+: four-vector = " << kineEvt.tauPlusP4().phi() << ", decay vertex = " << tauPlusD3.phi() << "\n";
        std::cout << " theta of tau+: four-vector = " << kineEvt.tauPlusP4().theta() << ", decay vertex = " << tauPlusD3.theta() << "\n";
      }
      return false;
    }
    if ( !(kineEvt.tauMinusP4().mass() >= 0. && kineEvt.tauMinusP4().mass() <= 5.) )
    { 
      if ( verbosity >= 1 )
      {
        std::cout << "fails tau- mass constraint (mass = " << kineEvt.tauMinusP4().mass() << ").\n";
      }
      return false;
    }
    auto tauMinusD3 = kineEvt.svTauMinus() - kineEvt.pv();
    if ( std::fabs(reco::deltaPhi(kineEvt.tauMinusP4().phi(), tauMinusD3.phi())) > max_dphi   || 
         std::fabs(kineEvt.tauMinusP4().theta() - tauMinusD3.theta())            > max_dtheta )
    {
      if ( verbosity >= 1 )
      {
        std::cout << "fails parallelism constraint for tau-.\n";
        std::cout << " phi of tau-: four-vector = " << kineEvt.tauMinusP4().phi() << ", decay vertex = " << tauMinusD3.phi() << "\n";
        std::cout << " theta of tau-: four-vector = " << kineEvt.tauMinusP4().theta() << ", decay vertex = " << tauMinusD3.theta() << "\n";
      }
      return false;
    }
    return true;
  }

  template <unsigned int numParameters, unsigned int numConstraints>
  void
  fit(const KinematicEvent& kineEvt, double signTauPlus, double signTauMinus, bool applyLifetimeConstraint, int collider,
      KinematicEvent& kineEvt_kinFit, int& status, double& min_chi2, int& iteration_bestfit, int max_iterations, bool& hasConverged,
      int verbosity)
  {
    const int C = numConstraints;  

    typedef typename Matrix<C,P>::type MatrixCxP;
    typedef typename Matrix<P,C>::type MatrixPxC;
    typedef typename Matrix<C,C>::type MatrixCxC;
    typedef typename Vector<C>::type   VectorC;

    int iteration = 0;
    VectorP alpha0;
    KinematicEvent kineEvtA = kineEvt;
    while ( !hasConverged && iteration < max_iterations )
    {
      if ( verbosity >= 1 )
      {
        std::cout << "sign(tau+) = " << signTauPlus << ", sign(tau-) = " << signTauMinus << ": iteration #" << iteration << ":\n";
      }

      //-----------------------------------------------------------------------------------------------
      // CV: take all covariance matrix from kineEvt and NOT from kineEvtA,
      //     as otherwise uncertainties will continuously shrink as iterations progress !!
      //-----------------------------------------------------------------------------------------------

      const reco::Candidate::Point& pv = kineEvtA.pv();
      double pvX = pv.X();
      double pvY = pv.Y();
      double pvZ = pv.Z();
      const Matrix3x3& pvCov = kineEvt.pvCov();

      const reco::Candidate::LorentzVector& recoilP4 = kineEvtA.recoilP4();
      double recoilPx = recoilP4.px();
      double recoilPy = recoilP4.py();
      double recoilPz = recoilP4.pz();
      double recoilE  = recoilP4.energy();
      const Matrix4x4& recoilCov = kineEvt.recoilCov();

      const reco::Candidate::LorentzVector& visTauPlusP4 = kineEvtA.visTauPlusP4();
      double visTauPlusPx = visTauPlusP4.px();
      double visTauPlusPy = visTauPlusP4.py();
      double visTauPlusPz = visTauPlusP4.pz();

      assert(kineEvtA.nuTauPlusP4_isValid());
      double nuTauPlusPx = kineEvtA.nuTauPlusP4().px();
      double nuTauPlusPy = kineEvtA.nuTauPlusP4().py();
      double nuTauPlus_dPzdPx, nuTauPlus_dPzdPy;
      double nuTauPlusPz = comp_nuPz(visTauPlusP4, nuTauPlusPx, nuTauPlusPy, signTauPlus, nuTauPlus_dPzdPx, nuTauPlus_dPzdPy, verbosity);
      reco::Candidate::LorentzVector nuTauPlusP4 = build_nuP4(nuTauPlusPx, nuTauPlusPy, nuTauPlusPz);
      Matrix2x2 nuTauPlusCov = kineEvt.nuTauPlusCov().Sub<Matrix2x2>(0,0);
      if ( verbosity >= 1 )
      {
        printLorentzVector("nuTauPlusP4", nuTauPlusP4);
      }

      double tauPlusPx = visTauPlusPx + nuTauPlusPx;
      double tauPlusPy = visTauPlusPy + nuTauPlusPy;
      double tauPlusPz = visTauPlusPz + nuTauPlusPz;
      double tauPlusE  = std::sqrt(square(tauPlusPx) + square(tauPlusPy) + square(tauPlusPz) + square(mTau));
      reco::Candidate::LorentzVector tauPlusP4(tauPlusPx, tauPlusPy, tauPlusPz, tauPlusE);
      double tauPlusPt = tauPlusP4.pt();
      double tauPlusP  = tauPlusP4.P();
      assert(kineEvtA.svTauPlus_isValid());
      const reco::Candidate::Point& svTauPlus = kineEvtA.svTauPlus();
      double svTauPlusX = svTauPlus.X();
      double svTauPlusY = svTauPlus.Y();
      double svTauPlusZ = svTauPlus.Z();
      const Matrix3x3& svTauPlusCov = kineEvt.svTauPlusCov();
      auto tauPlusD3 = svTauPlus - pv;
      double tauPlusDx = tauPlusD3.X();
      double tauPlusDy = tauPlusD3.Y();
      double tauPlusDz = tauPlusD3.Z();
      double tauPlusDt = std::sqrt(tauPlusD3.Perp2());
      double tauPlusD  = std::sqrt(tauPlusD3.Mag2());
      if ( verbosity >= 1 )
      {
        printLorentzVector("tauPlusP4", tauPlusP4);
      }

      const reco::Candidate::LorentzVector& visTauMinusP4 = kineEvtA.visTauMinusP4();
      double visTauMinusPx = visTauMinusP4.px();
      double visTauMinusPy = visTauMinusP4.py();
      double visTauMinusPz = visTauMinusP4.pz();
    
      assert(kineEvtA.nuTauMinusP4_isValid());
      double nuTauMinusPx = kineEvtA.nuTauMinusP4().px();
      double nuTauMinusPy = kineEvtA.nuTauMinusP4().py();
      double nuTauMinus_dPzdPx, nuTauMinus_dPzdPy;
      double nuTauMinusPz = comp_nuPz(visTauMinusP4, nuTauMinusPx, nuTauMinusPy, signTauMinus, nuTauMinus_dPzdPx, nuTauMinus_dPzdPy, verbosity);
      reco::Candidate::LorentzVector nuTauMinusP4 = build_nuP4(nuTauMinusPx, nuTauMinusPy, nuTauMinusPz);
      Matrix2x2 nuTauMinusCov = kineEvt.nuTauMinusCov().Sub<Matrix2x2>(0,0);
      if ( verbosity >= 1 )
      {
        printLorentzVector("nuTauMinusP4", nuTauMinusP4);
      }

      double tauMinusPx = visTauMinusPx + nuTauMinusPx;
      double tauMinusPy = visTauMinusPy + nuTauMinusPy;
      double tauMinusPz = visTauMinusPz + nuTauMinusPz;
      double tauMinusE  = std::sqrt(square(tauMinusPx) + square(tauMinusPy) + square(tauMinusPz) + square(mTau));
      reco::Candidate::LorentzVector tauMinusP4(tauMinusPx, tauMinusPy, tauMinusPz, tauMinusE);
      double tauMinusPt = tauMinusP4.pt();
      double tauMinusP  = tauMinusP4.P();
      assert(kineEvtA.svTauMinus_isValid());
      const reco::Candidate::Point& svTauMinus = kineEvtA.svTauMinus();
      double svTauMinusX = svTauMinus.X();
      double svTauMinusY = svTauMinus.Y();
      double svTauMinusZ = svTauMinus.Z();
      const Matrix3x3& svTauMinusCov = kineEvt.svTauMinusCov();
      auto tauMinusD3 = svTauMinus - pv;
      double tauMinusDx = tauMinusD3.X();
      double tauMinusDy = tauMinusD3.Y();
      double tauMinusDz = tauMinusD3.Z();
      double tauMinusDt = std::sqrt(tauMinusD3.Perp2());
      double tauMinusD  = std::sqrt(tauMinusD3.Mag2());
      if ( verbosity >= 1 )
      {
        printLorentzVector("tauMinusP4", tauMinusP4);
      }

      // CV: Check that chosen signs for tau+ and tau- reproduce neutrino Pz of start position; skip sign combination if not.
      //     Skipping these sign combinations saves computing time and avoids running into unphysical solutions.
      if ( iteration == 0 )
      {
        if ( std::fabs(nuTauPlusPz  - kineEvt.nuTauPlusP4().pz())  > std::max(1., 0.10*kineEvt.nuTauPlusP4().pz())  ||
             std::fabs(nuTauMinusPz - kineEvt.nuTauMinusP4().pz()) > std::max(1., 0.10*kineEvt.nuTauMinusP4().pz()) )
        {
          if ( verbosity >= 1 )
          {
            std::cout << "--> skipping this sign combination, because it does not reproduce the neutrino Pz of the start position !!\n";
          }
          return;
        }
      }

      //-----------------------------------------------------------------------------------------------
      // CV: The nomenclature of the different matrices and vectors used for the kinematic fit
      //     follows the one introduced in Section 3 of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
      //     For the first iteration, the vector alpha0 refers to the starting point of the kinematic fit,
      //     computed by the function KinematicFitStartPosFinder::operator()
      //     For subsequent iterations, the vector alpha0 refers to the result alpha obtained by the previous iteration
      //     The symbol H in the comments refers to the constraint equations,
      //     cf. Eq. (3) in https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
      //-----------------------------------------------------------------------------------------------

      // CV: Build vector of all measured parameters at point A 
      //     where constrained equations are linearized.
      //
      //     The measured parameters are defined in the following order:
      //       primary vertex position (x,y,z)       (3)
      //       Px, Py of neutrino from tau+          (2)
      //       decay vertex position (x,y,z) of tau+ (3)
      //       Px, Py of neutrino from tau-          (2)
      //       decay vertex position (x,y,z) of tau- (3)
      //       recoil four-vector (Px,Py,Pz,E)       (4)
      // where:
      //     the energy and momentum components of four-vectors are given in the order:
      //      (px, py, pz, E)
      //     the position of vertices are given in the order:
      //      (x, y, z)
      //     cf. Section II of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
      //
      //     The four-vectors of tau+ and tau- are not really "measured";
      //     we use the "huge error method" described in Section 6 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
      //     and set their covariance matrix to diagonal matrix with large values on the diagonal
      VectorP alphaA;
      alphaA( 0) = pvX;
      alphaA( 1) = pvY;
      alphaA( 2) = pvZ;
      alphaA( 3) = nuTauPlusPx;
      alphaA( 4) = nuTauPlusPy;
      alphaA( 5) = svTauPlusX;
      alphaA( 6) = svTauPlusY;
      alphaA( 7) = svTauPlusZ;
      alphaA( 8) = nuTauMinusPx;
      alphaA( 9) = nuTauMinusPy;
      alphaA(10) = svTauMinusX;
      alphaA(11) = svTauMinusY;
      alphaA(12) = svTauMinusZ;
      alphaA(13) = recoilPx;
      alphaA(14) = recoilPy;
      alphaA(15) = recoilPz;
      alphaA(16) = recoilE;
      if ( verbosity >= 1 )
      {
        std::cout << "alphaA:\n";
        std::cout << alphaA << "\n";
      }

      if ( iteration == 0 )
      {
        alpha0 = alphaA;
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
      if ( verbosity >= 1 )
      {
        printCovMatrix("V_alpha0", V_alpha0);
      }

      int Vinv_alpha0_errorFlag = 0;
      MatrixPxP Vinv_alpha0 = V_alpha0.Inverse(Vinv_alpha0_errorFlag);
      if ( Vinv_alpha0_errorFlag != 0 )
      {
        printCovMatrix("V_alpha0", V_alpha0);
        throw cmsException("KinematicFit::operator()", __LINE__)
          << "Failed to invert matrix V_alpha0 !!\n";
      }
      if ( verbosity >= 1 )
      {
        std::cout << "Vinv_alpha0:\n";
        std::cout << Vinv_alpha0 << "\n";
      }

      // CV: Define constraints.
      // 
      //     The constraints are defined in the following order:
      //       "parallelism" constraint for tau+ [1] (2)
      //       "parallelism" constraint for tau- [1] (2) 
      //       constraint that recoil = tau+ + tau-  (4)
      //       Higgs mass constraint                 (1)
      //     Note that the Higgs mass constraintis applied at the LHC,
      //     but not at the SuperKEKB collider (Belle)
      //  [1] cf. Section 4.1.3.3 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
      MatrixCxP D;
      VectorC   d;        
      // CV: add "parallelism" constraint for tau+
      //     The constraint equations H_Pphi and H_Ptheta are given by Eqs. (4.10) and (4.11)
      //     in Section 4.1.3.3 of [1].
      //     The small correction for the helix bend between tau lepton production and decay in Eq. (4.13) is not used, 
      //     as this correction is expected to be on the level of 50 microrad, which we consider negligible.
      //     Following [1], we perform a transformation of variables from phi and theta to 1/2*pi - phi and 1/2*pi - theta,
      //     in order to reduce the magnitude of derivatives (by avoiding the division by small numbers)
      //
      //   [1] https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
      VectorP D0_1;
      D0_1(0) =  (1. - tauPlusDx/tauPlusDt)/tauPlusDy;                                                                                     // dH_Pphi/dpvX
      D0_1(1) = -(tauPlusDx/tauPlusDy)*D0_1(0);                                                                                            // dH_Pphi/dpvY
      D0_1(3) =  (1. - tauPlusPx/tauPlusPt)/tauPlusPy;                                                                                     // dH_Pphi/dnuTauPlusPx
      D0_1(4) = -(tauPlusPx/tauPlusPy)*D0_1(3);                                                                                            // dH_Pphi/dnuTauPlusPy     
      D0_1(5) = -D0_1(0);                                                                                                                  // dH_Pphi/dsvTauPlusX
      D0_1(6) = -D0_1(1);                                                                                                                  // dH_Pphi/dsvTauPlusY
      double D0_1_mag = comp_mag(D0_1);
      VectorP D0_2;
      D0_2(0) = -(tauPlusDy/tauPlusDx)*(1. - tauPlusDy/tauPlusDt)/tauPlusDx;                                                               // dH_Pphi/dpvX
      D0_2(1) =  (1. - tauPlusDy/tauPlusDt)/tauPlusDx;                                                                                     // dH_Pphi/dpvY
      D0_2(3) = -(tauPlusPy/tauPlusPx)*(1. - tauPlusPy/tauPlusPt)/tauPlusPx;                                                               // dH_Pphi/dnuTauPlusPx
      D0_2(4) =  (1. - tauPlusPy/tauPlusPt)/tauPlusPx;                                                                                     // dH_Pphi/dnuTauPlusPy
      D0_2(5) = -D0_2(0);                                                                                                                  // dH_Pphi/dsvTauPlusX
      D0_2(6) = -D0_2(1);                                                                                                                  // dH_Pphi/dsvTauPlusY
      double D0_2_mag = comp_mag(D0_2);
      if ( verbosity >= 1 )
      {
        std::cout << "D0 (1st):\n";
        std::cout << D0_1 << "\n";
        std::cout << " mag = " << D0_1_mag << "\n";
        std::cout << "D0 (2nd):\n";
        std::cout << D0_2 << "\n";
        std::cout << " mag = " << D0_2_mag << "\n";
      }
      //if ( D0_1_mag < D0_2_mag )
      //{
        for ( int idx = 0; idx < kinFit::numParameters; ++idx )
        {
          D(0,idx) = D0_1(idx);
        }
        d(0) = (tauPlusDt - tauPlusDx)/tauPlusDy + (tauPlusPx - tauPlusPt)/tauPlusPy;
      //}
      //else
      //{
      //  for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      //  {
      //    D(0,idx) = D0_2(idx);
      //  }
      //  d(0) = (tauPlusDt - tauPlusDy)/tauPlusDx + (tauPlusPy - tauPlusPt)/tauPlusPx;
      //}
      VectorP D1_1;
      D1_1(0) =  (tauPlusDx/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);                                                                       // dH_Ptheta/dpvX
      D1_1(1) =  (tauPlusDy/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);                                                                       // dH_Ptheta/dpvY
      D1_1(2) = -1./tauPlusD + (tauPlusD - tauPlusDt)/square(tauPlusDz);                                                                   // dH_Ptheta/dpvZ
      D1_1(3) =  (1./square(tauPlusPz))*nuTauPlus_dPzdPx*(tauPlusP - tauPlusPt)                                                            // dH_Ptheta/dnuTauPlusPx
                 + (1./tauPlusPz)*(tauPlusPx/tauPlusPt - (tauPlusPx + nuTauPlus_dPzdPx*tauPlusPz)/tauPlusP);
      D1_1(4) =  (1./square(tauPlusPz))*nuTauPlus_dPzdPy*(tauPlusP - tauPlusPt)                                                            // dH_Ptheta/dnuTauPlusPy
                 + (1./tauPlusPz)*(tauPlusPy/tauPlusPt - (tauPlusPy + nuTauPlus_dPzdPy*tauPlusPz)/tauPlusP);
      D1_1(5) = -D1_1(0);                                                                                                                  // dH_Ptheta/dsvTauPlusX
      D1_1(6) = -D1_1(1);                                                                                                                  // dH_Ptheta/dsvTauPlusY
      D1_1(7) = -D1_1(2);                                                                                                                  // dH_Ptheta/dsvTauPlusZ
      double D1_1_mag = comp_mag(D1_1);
      VectorP D1_2;
      D1_2(0) = -(tauPlusDx*tauPlusDz/cube(tauPlusDt))*(1. - tauPlusDz/tauPlusD);                                                          // dH_Ptheta/dpvX
      D1_2(1) = -(tauPlusDy*tauPlusDz/cube(tauPlusDt))*(1. - tauPlusDz/tauPlusD);                                                          // dH_Ptheta/dpvY
      D1_2(2) =  (1. - tauPlusDz/tauPlusD)/tauPlusDt;                                                                                      // dH_Ptheta/dpvZ
      D1_2(3) =  (tauPlusPx/cube(tauPlusPt))*(tauPlusP - tauPlusPz)                                                                        // dH_Ptheta/dnuTauPlusPx
                 + (1./tauPlusPt)*(nuTauPlus_dPzdPx - (tauPlusPx + nuTauPlus_dPzdPx*tauPlusPz)/tauPlusP);
      D1_2(4) =  (tauPlusPy/cube(tauPlusPt))*(tauPlusP - tauPlusPz)                                                                        // dH_Ptheta/dnuTauPlusPy
                 + (1./tauPlusPt)*(nuTauPlus_dPzdPy - (tauPlusPy + nuTauPlus_dPzdPy*tauPlusPz)/tauPlusP);
      D1_2(5) = -D1_2(0);                                                                                                                  // dH_Ptheta/dsvTauPlusX
      D1_2(6) = -D1_2(1);                                                                                                                  // dH_Ptheta/dsvTauPlusY
      D1_2(7) = -D1_2(2);                                                                                                                  // dH_Ptheta/dsvTauPlusZ
      double D1_2_mag = comp_mag(D1_2);
      if ( verbosity >= 1 )
      {
        std::cout << "D1 (1st):\n";
        std::cout << D1_1 << "\n";
        std::cout << " mag = " << D1_1_mag << "\n";
        std::cout << "D1 (2nd):\n";
        std::cout << D1_2 << "\n";
        std::cout << " mag = " << D1_2_mag << "\n";
      }
      //if ( D1_1_mag < D1_2_mag )
      //{
        for ( int idx = 0; idx < kinFit::numParameters; ++idx )
        {
          D(1,idx) = D1_1(idx);
        }
        d(1) = (tauPlusD - tauPlusDt)/tauPlusDz + (tauPlusPt - tauPlusP)/tauPlusPz;
      //}
      //else
      //{
      //  for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      //  {
      //     D(1,idx) = D1_2(idx);
      //  }
      //  d(1) = (tauPlusD - tauPlusDz)/tauPlusDt + (tauPlusPz - tauPlusP)/tauPlusPt;
      //}
      // CV: add "parallelism" constraint for tau-
      //     The constraint equations H_Pphi and H_Ptheta are given by Eqs. (4.10) and (4.11)
      //     in Section 4.1.3.3 of [1].
      //     The small correction for the helix bend between tau lepton production and decay in Eq. (4.13) is not used, 
      //     as this correction is expected to be on the level of 50 microrad, which we consider negligible.
      //     Following [1], we perform a transformation of variables from phi and theta to 1/2*pi - phi and 1/2*pi - theta,
      //     in order to reduce the magnitude of derivatives (by avoiding the division by small numbers)
      //
      //   [1] https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
      VectorP D2_1;
      D2_1( 0) =  (1. - tauMinusDx/tauMinusDt)/tauMinusDy;                                                                                 // dH_Pphi/dpvX
      D2_1( 1) = -(tauMinusDx/tauMinusDy)*D2_1(0);                                                                                         // dH_Pphi/dpvY
      D2_1( 8) =  (1. - tauMinusPx/tauMinusPt)/tauMinusPy;                                                                                 // dH_Pphi/dnuTauMinusPx
      D2_1( 9) = -(tauMinusPx/tauMinusPy)*D2_1(8);                                                                                         // dH_Pphi/dnuTauMinusPy
      D2_1(10) = -D2_1(0);                                                                                                                 // dH_Pphi/dsvTauMinusX
      D2_1(11) = -D2_1(1);                                                                                                                 // dH_Pphi/dsvTauMinusY
      double D2_1_mag = comp_mag(D2_1);
      VectorP D2_2;
      D2_2( 0) = -(tauMinusDy/tauMinusDx)*(1. - tauMinusDy/tauMinusDt)/tauMinusDx;                                                         // dH_Pphi/dpvX
      D2_2( 1) =  (1. - tauMinusDy/tauMinusDt)/tauMinusDx;                                                                                 // dH_Pphi/dpvY
      D2_2( 8) = -(tauMinusPy/tauMinusPx)*(1. - tauMinusPy/tauMinusPt)/tauMinusPx;                                                         // dH_Pphi/dnuTauMinusPx
      D2_2( 9) =  (1. - tauMinusPy/tauMinusPt)/tauMinusPx;                                                                                 // dH_Pphi/dnuTauMinusPy
      D2_2(10) = -D2_2(0);                                                                                                                 // dH_Pphi/dsvTauMinusX
      D2_2(11) = -D2_2(1);                                                                                                                 // dH_Pphi/dsvTauMinusY
      double D2_2_mag = comp_mag(D2_2);
      if ( verbosity >= 1 )
      {
        std::cout << "D2 (1st):\n";
        std::cout << D2_1 << "\n";
        std::cout << " mag = " << D2_2_mag << "\n";
        std::cout << "D2 (2nd):\n";
        std::cout << D2_2 << "\n";
        std::cout << " mag = " << D2_2_mag << "\n";
      }
      //if ( D2_1_mag < D2_2_mag )
      //{
        for ( int idx = 0; idx < kinFit::numParameters; ++idx )
        {
          D(2,idx) = D2_1(idx);
        }
        d(2) =  (tauMinusDt - tauMinusDx)/tauMinusDy + (tauMinusPx - tauMinusPt)/tauMinusPy;
      //}
      //else
      //{
      //  for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      //  {
      //    D(2,idx) = D2_2(idx);
      //  }
      //  d(2) =  (tauMinusDt - tauMinusDy)/tauMinusDx + (tauMinusPy - tauMinusPt)/tauMinusPx;
      //}
      VectorP D3_1;
      D3_1( 0) =  (tauMinusDx/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);                                                                  // dH_Ptheta/dpvX
      D3_1( 1) =  (tauMinusDy/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);                                                                  // dH_Ptheta/dpvY
      D3_1( 2) = -1./tauMinusD + (tauMinusD - tauMinusDt)/square(tauMinusDz);                                                              // dH_Ptheta/dpvZ
      D3_1( 8) =  (1./square(tauMinusPz))*nuTauMinus_dPzdPx*(tauMinusP - tauMinusPt)                                                       // dH_Ptheta/dnuTauMinusPx
                 + (1./tauMinusPz)*(tauMinusPx/tauMinusPt - (tauMinusPx + nuTauMinus_dPzdPx*tauMinusPz)/tauMinusP);
      D3_1( 9) =  (1./square(tauMinusPz))*nuTauMinus_dPzdPy*(tauMinusP - tauMinusPt)                                                       // dH_Ptheta/dnuTauMinusPy
                 + (1./tauMinusPz)*(tauMinusPy/tauMinusPt - (tauMinusPy + nuTauMinus_dPzdPy*tauMinusPz)/tauMinusP);
      D3_1(10) = -D3_1(0);                                                                                                                 // dH_Ptheta/dsvTauMinusX
      D3_1(11) = -D3_1(1);                                                                                                                 // dH_Ptheta/dsvTauMinusY
      D3_1(12) = -D3_1(2);                                                                                                                 // dH_Ptheta/dsvTauMinusZ
      double D3_1_mag = comp_mag(D3_1);
      VectorP D3_2;
      D3_2( 0) = -(tauMinusDx*tauMinusDz/cube(tauMinusDt))*(1. - tauMinusDz/tauMinusD);                                                    // dH_Ptheta/dpvX
      D3_2( 1) = -(tauMinusDy*tauMinusDz/cube(tauMinusDt))*(1. - tauMinusDz/tauMinusD);                                                    // dH_Ptheta/dpvY
      D3_2( 2) =  (1. - tauMinusDz/tauMinusD)/tauMinusDt;                                                                                  // dH_Ptheta/dpvZ
      D3_2( 8) = (tauMinusPx/cube(tauMinusPt))*(tauMinusP - tauMinusPz)                                                                    // dH_Ptheta/dnuTauMinusPx
                 + (1./tauMinusPt)*(nuTauMinus_dPzdPx - (tauMinusPx + nuTauMinus_dPzdPx*tauMinusPz)/tauMinusP);
      D3_2( 9) = (tauMinusPy/cube(tauMinusPt))*(tauMinusP - tauMinusPz)                                                                    // dH_Ptheta/dnuTauMinusPy
                 + (1./tauMinusPt)*(nuTauMinus_dPzdPy - (tauMinusPy + nuTauMinus_dPzdPy*tauMinusPz)/tauMinusP);
      D3_2(10) = -D3_2(0);                                                                                                                 // dH_Ptheta/dsvTauMinusX
      D3_2(11) = -D3_2(1);                                                                                                                 // dH_Ptheta/dsvTauMinusY
      D3_2(12) = -D3_2(2);                                                                                                                 // dH_Ptheta/dsvTauMinusZ
      double D3_2_mag = comp_mag(D3_2);
      if ( verbosity >= 1 )
      {
        std::cout << "D3 (1st):\n";
        std::cout << D3_1 << "\n";
        std::cout << " mag = " << D3_1_mag << "\n";
        std::cout << "D3 (2nd):\n";
        std::cout << D3_2 << "\n";
        std::cout << " mag = " << D3_2_mag << "\n";
      }
      //if ( D3_1_mag < D3_2_mag )
      //{
        for ( int idx = 0; idx < kinFit::numParameters; ++idx )
        {
          D(3,idx) = D3_1(idx);
        }
        d(3) =  (tauMinusD - tauMinusDt)/tauMinusDz + (tauMinusPt - tauMinusP)/tauMinusPz;
      //}
      //else
      //{
      //  for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      //  {
      //    D(3,idx) = D3_2(idx);
      //  }
      //  d(3) =  (tauMinusD - tauMinusDz)/tauMinusDt + (tauMinusPz - tauMinusP)/tauMinusPt;
      //}
      // CV: add constraint that recoil = tau+ + tau-
      //       H_px     = nuTauPlusPx + visTauPlusPx + nuTauMinusPx + visTauMinusPx - recoilPx = 0
      //       H_py     = nuTauPlusPy + visTauPlusPy + nuTauMinusPy + visTauMinusPy - recoilPy = 0
      //       H_pz     = nuTauPlusPz + visTauPlusPz + nuTauMinusPz + visTauMinusPz - recoilPz = 0
      //       H_energy = nuTauPlusE  + visTauPlusE  + nuTauMinusE  + visTauMinusE  - recoilE  = 0
      //     where:
      //       nuTauPlusE  = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2 )
      //       nuTauMinusE = sqrt(nuTauMisusPx^2 + nuTauMinusPy^2 + nuTauMinusPz^2)
      D(4, 3) = +1.;                                                                                                                       // dH_px/dnuTauPlusPx
      D(4, 8) = +1.;                                                                                                                       // dH_px/dnuTauMinusPx
      D(4,13) = -1.;                                                                                                                       // dH_px/drecoilPx
      d(4) = tauPlusPx + tauMinusPx - recoilPx;
      D(5, 4) = +1.;                                                                                                                       // dH_py/dnuTauPlusPy
      D(5, 9) = +1.;                                                                                                                       // dH_py/dnuTauMinusPy
      D(5,14) = -1.;                                                                                                                       // dH_py/drecoilPy
      d(5) = tauPlusPy + tauMinusPy - recoilPy;
      D(6, 3) = nuTauPlus_dPzdPx;                                                                                                          // dH_pz/dnuTauPlusPx
      D(6, 4) = nuTauPlus_dPzdPy;                                                                                                          // dH_pz/dnuTauPlusPy
      D(6, 8) = nuTauMinus_dPzdPx;                                                                                                         // dH_pz/dnuTauMinusPx
      D(6, 9) = nuTauMinus_dPzdPy;                                                                                                         // dH_pz/dnuTauMinusPy
      D(6,15) = -1.;                                                                                                                       // dH_pz/drecoilPz
      d(6) = visTauPlusPz + nuTauPlusPz + visTauMinusPz + nuTauMinusPz - recoilPz;
      D(7, 3) = (tauPlusPx  + nuTauPlus_dPzdPx*tauPlusPz)/tauPlusE;                                                                        // dH_energy/dnuTauPlusPx
      D(7, 4) = (tauPlusPy  + nuTauPlus_dPzdPy*tauPlusPz)/tauPlusE;                                                                        // dH_energy/dnuTauPlusPx
      D(7, 8) = (tauMinusPx + nuTauMinus_dPzdPx*tauMinusPz)/tauMinusE;                                                                     // dH_energy/dnuTauPlusPx
      D(7, 9) = (tauMinusPy + nuTauMinus_dPzdPy*tauMinusPz)/tauMinusE;                                                                     // dH_energy/dnuTauPlusPx
      D(7,16) = -1.;                                                                                                                       // dH_energy/drecoilE
      d(7) = tauPlusE + tauMinusE - recoilE;
      if ( collider != kSuperKEKB )
      {
        // CV: add Higgs mass constraint
        //       H = sqrt(2*mTau^2 
        //              + 2*(visTauPlusE  + nuTauPlusE )*(visTauMinusE  + nuTauMinusE ) 
        //              - 2*(visTauPlusPx + nuTauPlusPx)*(visTauMinusPx + nuTauMinusPx)
        //              - 2*(visTauPlusPy + nuTauPlusPy)*(visTauMinusPy + nuTauMinusPy)
        //              - 2*(visTauPlusPz + nuTauPlusPz)*(visTauMinusPz + nuTauMinusPz)) - mHiggs = 0
        //     where
        //       nuTauPlusE  = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2 )
        //       nuTauMinusE = sqrt(nuTauMisusPx^2 + nuTauMinusPy^2 + nuTauMinusPz^2)
        double denominator0 = (tauPlusP4 + tauMinusP4).mass();
        if ( verbosity >= 1 )
        {
          std::cout << "mTauTau = " << (tauPlusP4 + tauMinusP4).mass() << "\n";
        }
        // CV: Setting idx to 8 generates compile-time errors that array indices are out-of-bounds
        //     for ROOT::Math::SMatrix and ROOT::Math::SVector objects.
        //     The error disappears when setting idx to numConstraints - 1 instead.
        //     Note that the index will be wrong when collider = SuperKEKB,
        //     but this does not matter as the code will actually be executed for collider = LHC only.
        const int idx = numConstraints - 1;
        D(idx,3) = ((tauMinusE/tauPlusE)*(tauPlusPx  + nuTauPlus_dPzdPx*tauPlusPz)   - tauMinusPx - nuTauPlus_dPzdPx*tauMinusPz)/denominator0; // dH/dnuTauPlusPx
        D(idx,4) = ((tauMinusE/tauPlusE)*(tauPlusPy  + nuTauPlus_dPzdPy*tauPlusPz)   - tauMinusPy - nuTauPlus_dPzdPy*tauMinusPz)/denominator0; // dH/dnuTauPlusPy
        D(idx,8) = ((tauPlusE/tauMinusE)*(tauMinusPx + nuTauMinus_dPzdPx*tauMinusPz) - tauPlusPx  - nuTauMinus_dPzdPx*tauPlusPz)/denominator0; // dH/dnuTauMinusPx
        D(idx,9) = ((tauPlusE/tauMinusE)*(tauMinusPy + nuTauMinus_dPzdPy*tauMinusPz) - tauPlusPy  - nuTauMinus_dPzdPy*tauPlusPz)/denominator0; // dH/dnuTauMinusPy
        d(idx) = (tauPlusP4 + tauMinusP4).mass() - mHiggs;
      }
      if ( verbosity >= 1 )
      {
        std::cout << "D:\n";
        std::cout << D << "\n";
        std::cout << "d:\n";
        std::cout << d << "\n";
      }

      VectorP p;
      if ( applyLifetimeConstraint )
      {
        // CV: add "soft" constraint on tau decay probability
        //             P  = exp(-d*mTau/(ct*tauPlusP))*exp(-d*mTau/(ct*tauMinusP))
        //    => -2*ln(P) = 2*d*mtau/(ct*(visTauPlusP + nuTauPlusP)) + 2*d*mtau/(ct*(visTauMinusP + nuTauMinusP))
        //
        //     where
        //       nuTauPlusP  = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2 )
        //       nuTauMinusP = sqrt(nuTauMisusPx^2 + nuTauMinusPy^2 + nuTauMinusPz^2)        
        double term1TauPlus  = 2.*mTau/(ct*tauPlusD*tauPlusP);
        double term2TauPlus  = 2.*mTau*tauPlusD/(ct*cube(tauPlusP));
        double term1TauMinus = 2.*mTau/(ct*tauMinusD*tauMinusP);
        double term2TauMinus = 2.*mTau*tauMinusD/(ct*cube(tauMinusP));
        p( 0) = -(term1TauPlus*tauPlusDx + term1TauMinus*tauMinusDx);                      // d(-2*ln(P))/pvX
        p( 1) = -(term1TauPlus*tauPlusDy + term1TauMinus*tauMinusDy);                      // d(-2*ln(P))/pvY
        p( 2) = -(term1TauPlus*tauPlusDz + term1TauMinus*tauMinusDz);                      // d(-2*ln(P))/pvZ
        p( 3) = -term2TauPlus*(tauPlusPx + nuTauPlus_dPzdPx*tauPlusPz);                    // d(-2*ln(P))/dnuTauPlusPx
        p( 4) = -term2TauPlus*(tauPlusPy + nuTauPlus_dPzdPy*tauPlusPz);                    // d(-2*ln(P))/dnuTauPlusPy
        p( 5) =  term1TauPlus*tauPlusDx;                                                   // d(-2*ln(P))/dsvTauPlusX
        p( 6) =  term1TauPlus*tauPlusDy;                                                   // d(-2*ln(P))/dsvTauPlusY
        p( 7) =  term1TauPlus*tauPlusDz;                                                   // d(-2*ln(P))/dsvTauPlusZ
        p( 8) = -term2TauMinus*(tauMinusPx + nuTauMinus_dPzdPx*tauMinusPz);                // d(-2*ln(P))/dnuTauMinusPx
        p( 9) = -term2TauMinus*(tauMinusPy + nuTauMinus_dPzdPy*tauMinusPz);                // d(-2*ln(P))/dnuTauMinusPy
        p(10) =  term1TauMinus*tauMinusDx;                                                 // d(-2*ln(P))/dsvTauMinusX
        p(11) =  term1TauMinus*tauMinusDy;                                                 // d(-2*ln(P))/dsvTauMinusY
        p(12) =  term1TauMinus*tauMinusDz;                                                 // d(-2*ln(P))/dsvTauMinusZ
        if ( verbosity >= 1 )
        {
          std::cout << "p:\n";
          std::cout << p << "\n";
        }
      }

      //-----------------------------------------------------------------------------------------------
      // CV: compute solution to minimization problem
      //     using formulas given in Section 3 of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
      //     For the syntax of matrix operations,
      //     see Sections "Matrix and vector functions" and "Linear algebra functions"
      //     of the ROOT documentation at
      //       https://root.cern.ch/doc/v608/MatVecFunctions.html 
      //     and
      //       https://root.cern.ch/doc/v608/SMatrixDoc.html 
      //     respectively.
      //-----------------------------------------------------------------------------------------------

      MatrixPxC DT = ROOT::Math::Transpose(D);    
      if ( verbosity >= 1 )
      {
        std::cout << "DT:\n";
        std::cout << DT << "\n";
      }

      MatrixCxC Vinv_D = D*V_alpha0*DT;
      int V_D_errorFlag = 0;
      MatrixCxC V_D = Vinv_D.Inverse(V_D_errorFlag);
      if ( V_D_errorFlag != 0 )
      {
        printCovMatrix("Vinv_D", Vinv_D);
        throw cmsException("KinematicFit::operator()", __LINE__)
          << "Failed to invert matrix Vinv_D !!\n";
      }
      if ( verbosity >= 1 )
      {
        printCovMatrix("V_D", V_D);
      }

      VectorP dalpha0 = alpha0 - alphaA;
      VectorC lambda = V_D*(D*dalpha0 + d);
      if ( applyLifetimeConstraint )
      {
        lambda -= V_D*D*V_alpha0*p;
      }
      if ( verbosity >= 1 )
      {
        std::cout << "dalpha0:\n";
        std::cout << dalpha0 << "\n";
        std::cout << "D*dalpha0:\n";
        std::cout << D*dalpha0 << "\n";
        std::cout << "lambda:\n";
        std::cout << lambda << "\n";
        std::cout << "contributions to lambda:\n";
        VectorC D_times_dalpha0 = D*dalpha0;
        for ( int idxRow = 0; idxRow < (int)numConstraints; ++idxRow )
        {
          std::cout << "lambda[" << idxRow << "]:";
          for ( int idxColumn = 0; idxColumn < (int)numConstraints; ++idxColumn )
          {
            if ( idxColumn != 0 ) std::cout << ",";
            std::cout << " " << V_D(idxRow,idxColumn)*D_times_dalpha0(idxColumn) << " + " << V_D(idxRow,idxColumn)*d(idxColumn);
          }
          std::cout << "\n";
        }
        std::cout << "DT*lambda:\n";
        std::cout << DT*lambda << "\n";
        std::cout << "distance from satisfaction:\n";
        for ( int idx = 0; idx < (int)numConstraints; ++idx )
        {
          double pull = (D_times_dalpha0(idx) + d(idx))/std::sqrt(Vinv_D(idx,idx));
          std::cout << "pull[" << idx << "] = " << pull << "\n";
        }
      }

      VectorP alpha = alpha0 - V_alpha0*DT*lambda;
      if ( applyLifetimeConstraint )
      {
        alpha -= V_alpha0*p;
      }
      if ( verbosity >= 1 )
      {
        std::cout << "alpha0:\n";
        std::cout << alpha0 << "\n";
        std::cout << "alpha:\n";
        std::cout << alpha << "\n";
      }
 
      MatrixPxP V_alpha = V_alpha0 - V_alpha0*DT*V_D*D*V_alpha0;
      if ( verbosity >= 1 )
      {
        printCovMatrix("V_alpha", V_alpha);
      }

      VectorP dalpha = alpha - alphaA;
      double dalpha_mag = std::sqrt(ROOT::Math::Dot(dalpha, dalpha));
      if ( verbosity >= 1 )
      {
        std::cout << "dalpha:\n";
        std::cout << dalpha << "\n";
        std::cout << "|alpha - alphaA| = " << dalpha_mag << "\n";
      }

      // CV: compute chi^2
      VectorP alpha_minus_alpha0 = alpha - alpha0;
      double chi2 = ROOT::Math::Dot(alpha_minus_alpha0, Vinv_alpha0*alpha_minus_alpha0) + ROOT::Math::Dot(lambda, Vinv_D*lambda);
      if ( applyLifetimeConstraint )
      {
        chi2 += 2.*tauPlusD*mTau/(ct*tauPlusP) + 2.*tauMinusD*mTau/(ct*tauMinusP);
      }
      chi2 /= kinFit::numParameters;
      if ( verbosity >= 1 )
      {
        std::cout << "chi^2/DoF = " << chi2 << "\n";
      }

      // CV: compute residuals
      VectorC residuals = D*dalpha + d;
      double residuals_sum = 0.;
      for ( int idx = 0; idx < (int)numConstraints; ++idx )
      {
        residuals_sum += std::fabs(residuals(idx));
      }
      if ( verbosity >= 1 )
      {
        std::cout << "residuals of constraint equations:\n";
        std::cout << residuals << "\n";
        std::cout << "sum_i |residual[i]| = " << residuals_sum << "\n";
      }

      if ( verbosity >= 1 )
      {
        VectorP pulls;
        for ( int idx = 0; idx < (int)numParameters; ++idx )
        {      
          pulls(idx) = alpha_minus_alpha0(idx)/std::sqrt(V_alpha0(idx,idx));
        }
        std::cout << "pulls:\n";
        std::cout << pulls << "\n";
      }

      // CV: store results of kinematic fit in KinematicEvent class;
      //     for the syntax of retrieving a small covariance matrix from a larger one, 
      //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
      kineEvtA.pv_ = reco::Candidate::Point(alpha(0), alpha(1), alpha(2));
      kineEvtA.pvCov_ = V_alpha.Sub<Matrix3x3>(0,0);
      //kineEvtA.recoilP4_ = fixHiggsMass(reco::Candidate::LorentzVector(alpha(13), alpha(14), alpha(15), alpha(16)));
      kineEvtA.recoilP4_ = reco::Candidate::LorentzVector(alpha(13), alpha(14), alpha(15), alpha(16));
      kineEvtA.recoilCov_ = V_alpha.Sub<Matrix4x4>(13,13);
      double nuTauPlusPxA = alpha(3);
      double nuTauPlusPyA = alpha(4);
      double nuTauPlus_dPzdPxA, nuTauPlus_dPzdPyA;
      double nuTauPlusPzA = comp_nuPz(visTauPlusP4, nuTauPlusPxA, nuTauPlusPyA, signTauPlus, nuTauPlus_dPzdPxA, nuTauPlus_dPzdPyA, verbosity);
      kineEvtA.nuTauPlusP4_ = build_nuP4(nuTauPlusPxA, nuTauPlusPyA, nuTauPlusPzA);
      kineEvtA.nuTauPlusP4_isValid_ = true;
      kineEvtA.nuTauPlusCov_ = build_nuCov(V_alpha.Sub<Matrix2x2>(3,3), nuTauPlus_dPzdPxA, nuTauPlus_dPzdPyA);
      //kineEvtA.tauPlusP4_ = fixTauMass(kineEvtA.visTauPlusP4_ + kineEvtA.nuTauPlusP4_);    
      kineEvtA.tauPlusP4_ = kineEvtA.visTauPlusP4_ + kineEvtA.nuTauPlusP4_;
      kineEvtA.tauPlusP4_isValid_ = true;
      kineEvtA.svTauPlus_ = reco::Candidate::Point(alpha(5), alpha(6), alpha(7));
      kineEvtA.svTauPlusCov_ = V_alpha.Sub<Matrix3x3>(5,5);
      double nuTauMinusPxA = alpha(8);
      double nuTauMinusPyA = alpha(9);
      double nuTauMinus_dPzdPxA, nuTauMinus_dPzdPyA;
      double nuTauMinusPzA = comp_nuPz(visTauMinusP4, nuTauMinusPxA, nuTauMinusPyA, signTauMinus, nuTauMinus_dPzdPxA, nuTauMinus_dPzdPyA, verbosity);
      kineEvtA.nuTauMinusP4_ = build_nuP4(nuTauMinusPxA, nuTauMinusPyA, nuTauMinusPzA);
      kineEvtA.nuTauMinusP4_isValid_ = true;       
      kineEvtA.nuTauMinusCov_ = build_nuCov(V_alpha.Sub<Matrix2x2>(8,8), nuTauMinus_dPzdPxA, nuTauMinus_dPzdPyA);
      //kineEvtA.tauMinusP4_ = fixTauMass(kineEvtA.visTauMinusP4_ + kineEvtA.nuTauMinusP4_);
      kineEvtA.tauMinusP4_ = kineEvtA.visTauMinusP4_ + kineEvtA.nuTauMinusP4_;
      kineEvtA.tauMinusP4_isValid_ = true;
      kineEvtA.svTauMinus_ = reco::Candidate::Point(alpha(10), alpha(11), alpha(12));
      kineEvtA.svTauMinusCov_ = V_alpha.Sub<Matrix3x3>(10,10);
      if ( verbosity >= 1 )
      {
        std::cout << "constraint equations:\n";
        printPoint("PV", kineEvtA.pv());
        reco::Candidate::LorentzVector higgsP4 = kineEvtA.tauPlusP4() + kineEvtA.tauMinusP4();
        printLorentzVector("Higgs", higgsP4);
        std::cout << " mass = " << higgsP4.mass() << "\n";
        printLorentzVector("tau+", kineEvtA.tauPlusP4());
        std::cout << " mass = " << kineEvtA.tauPlusP4().mass() << "\n";
        printLorentzVector("neutrino from tau+ decay", kineEvtA.nuTauPlusP4_);
        std::cout << " mass = " << kineEvtA.nuTauPlusP4().mass() << "\n";
        printPoint("SV(tau+)", kineEvtA.svTauPlus());
        auto tauPlusD3 = kineEvtA.svTauPlus() - kineEvtA.pv();
        std::cout << "phi of tau+: four-vector = " << kineEvtA.tauPlusP4().phi() << ", decay vertex = " << tauPlusD3.phi() << "\n";
        std::cout << "theta of tau+: four-vector = " << kineEvtA.tauPlusP4().theta() << ", decay vertex = " << tauPlusD3.theta() << "\n";
        printLorentzVector("tau-", kineEvtA.tauMinusP4());
        std::cout << " mass = " << kineEvtA.tauMinusP4().mass() << "\n";
        printLorentzVector("neutrino from tau- decay", kineEvtA.nuTauMinusP4());
        std::cout << " mass = " << kineEvtA.nuTauMinusP4().mass() << "\n";
        printPoint("SV(tau-)", kineEvtA.svTauMinus());
        auto tauMinusD3 = kineEvtA.svTauMinus() - kineEvtA.pv();
        std::cout << "phi of tau-: four-vector = " << kineEvtA.tauMinusP4().phi() << ", decay vertex = " << tauMinusD3.phi() << "\n";
        std::cout << "theta of tau-: four-vector = " << kineEvtA.tauMinusP4().theta() << ", decay vertex = " << tauMinusD3.theta() << "\n";
        printLorentzVector("recoil", kineEvtA.recoilP4());
        std::cout << " mass = " << kineEvtA.recoilP4().mass() << "\n";
        std::cout << "Higgs - recoil:"
                  << " dE = "  << higgsP4.energy() - kineEvtA.recoilP4().energy() << ","
                  << " dPx = " << higgsP4.px()     - kineEvtA.recoilP4().px()     << ","
                  << " dPy = " << higgsP4.py()     - kineEvtA.recoilP4().py()     << ","
                  << " dPz = " << higgsP4.pz()     - kineEvtA.recoilP4().pz()     << "\n";
      }

      if ( verbosity >= 1 )
      {
        std::cout << "isNaN = " << isNaN(kineEvtA, verbosity) << "\n";
        std::cout << "isPhysicalSolution = " << isPhysicalSolution(kineEvtA, collider, verbosity) << "\n";
      }

      if ( !isNaN(kineEvtA) && isPhysicalSolution(kineEvtA, collider) )
      {
        if ( status == -1 || chi2 < min_chi2 )
        {
          kineEvt_kinFit = kineEvtA;
          kineEvt_kinFit.kinFitCov_ = V_alpha;
          kineEvt_kinFit.kinFitChi2_ = chi2;
          kineEvt_kinFit.kinFit_isValid_ = true;
          min_chi2 = chi2;
          status = 0;
          iteration_bestfit = iteration;

          if ( dalpha_mag < 1.e-1 )
          {
            status = 1;
            hasConverged = true;
          }
        }
      }

      if ( isNaN(kineEvtA) )
      {
        // CV: there is no hope that the KinematicFit may still converge,
        //     once one of the paramaters becomes "not-a-number" (NaN)
        std::cerr << "WARNING: Parameters contain NaNs -> aborting KinematicFit at iteration #" << iteration << " !!" << std::endl;
        break;
      }

      ++iteration;
    }
  }

  const int numConstraints_LHC = 9;
  const int numConstraints_SuperKEKB = 8;
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

  double min_chi2 = -1.;
  int status = -1;
  int iteration_bestfit = -1;
  const int max_iterations = 10;
  bool hasConverged = false;
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

      if ( collider_ == kLHC ) 
      {
        fit<numParameters,numConstraints_LHC>(
          kineEvt, signTauPlus, signTauMinus, applyLifetimeConstraint_, collider_,
          kineEvt_kinFit, status, min_chi2, iteration_bestfit, max_iterations, hasConverged,
          verbosity_);
      }
      else if ( collider_ == kSuperKEKB )
      {
        fit<numParameters,numConstraints_SuperKEKB>(
          kineEvt, signTauPlus, signTauMinus, applyLifetimeConstraint_, collider_,
          kineEvt_kinFit, status, min_chi2, iteration_bestfit, max_iterations, hasConverged,
          verbosity_);
      }
      else assert(0);
    }
  }
  if ( status == 0 && min_chi2 < 1.e+2 )
  {
    hasConverged = true;
  }
  if ( !hasConverged )
  {
    std::cerr << "WARNING: KinematicFit failed to converge !!" << std::endl;
  }
  if ( verbosity_ >= 1 )
  {
    std::cout << "best fit:\n";
    std::cout << " iteration = " << iteration_bestfit << " (max_iterations = " << max_iterations << ")\n";
    std::cout << " status = " << status << "\n";
    std::cout << " min(chi^2) = " << min_chi2 << "\n";
  }

  kineEvt_kinFit.kinFitStatus_ = status;
  kineEvt_kinFit.recoilP4_ = kineEvt_kinFit.tauPlusP4_ + kineEvt_kinFit.tauMinusP4_;
  if ( hasConverged )
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

