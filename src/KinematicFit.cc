#include "TauAnalysis/Entanglement/interface/KinematicFit.h"

#include "DataFormats/Math/interface/deltaPhi.h"                          // deltaPhi()
#include "DataFormats/Math/interface/Matrix.h"                            // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                            // math::Vector
#include "DataFormats/TauReco/interface/PFTau.h"                          // reco::PFTau::kOneProng0PiZero

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/comp_cosThetaGJ.h"           // comp_cosThetaGJ(), comp_cosThetaGJ_solution(), kPlus, kMinus
#include "TauAnalysis/Entanglement/interface/comp_mag.h"                  // comp_mag()
#include "TauAnalysis/Entanglement/interface/constants.h"                 // ct, mHiggs, mTau
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

KinematicFit::KinematicFit(const edm::ParameterSet& cfg)
  : spinAnalyzer_(cfg)
  , applyTauMassConstraint_(cfg.getParameter<int>("applyTauMassConstraint"))
  , applyLifetimeConstraint_(cfg.getParameter<bool>("applyLifetimeConstraint"))
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{}

KinematicFit::~KinematicFit()
{}

namespace
{
  reco::Candidate::LorentzVector
  build_nuP4(double nuPx, double nuPy, double nuPz)
  {
    double nuE = std::sqrt(square(nuPx) + square(nuPy) + square(nuPz));
    reco::Candidate::LorentzVector nuP4(nuPx, nuPy, nuPz, nuE);
    return nuP4;
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
  isPhysicalSolution(const KinematicEvent& kineEvt, int verbosity = -1)
  {
    // CV: check that constraint equations are satisfied
    //    (not just their linear approximation)
    if ( verbosity >= 1 )
    {
      std::cout << "<isPhysicalSolution>:\n";
    }
    const double max_dphi = (10./180.)*TMath::Pi();
    const double max_dtheta = (10./180.)*TMath::Pi();
    reco::Candidate::LorentzVector higgsP4 = kineEvt.tauPlusP4() + kineEvt.tauMinusP4();
    if ( std::fabs(higgsP4.mass() - mHiggs) > 10. )
    {
      if ( verbosity >= 1 )
      {
        std::cout << "fails Higgs mass constraint (mass = " << higgsP4.mass() << ").\n";
      }
      return false;
    }
    if ( std::fabs(kineEvt.tauPlusP4().mass() - mTau) > mTau )
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
    if ( std::fabs(kineEvt.tauMinusP4().mass() - mTau) > mTau )
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
}

KinematicEvent
KinematicFit::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<KinematicFit::operator()>:\n"; 
  }

  KinematicEvent kineEvt_kinfit = kineEvt;

  int iteration = 0;
  const int max_iterations = 10;
  math::VectorP alpha0;
  KinematicEvent kineEvtA = kineEvt;
  double min_chi2 = -1.;
  int status = -1;
  bool hasConverged = false;
  while ( !hasConverged && iteration < max_iterations )
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "iteration #" << iteration << ":\n";
    }

    //-----------------------------------------------------------------------------------------------
    // CV: take all covariance matrix from kineEvt and NOT from kineEvtA,
    //     as otherwise uncertainties will continuously shrink as iterations progress !!
    //-----------------------------------------------------------------------------------------------

    const reco::Candidate::Point& pv = kineEvtA.pv();
    double pvX = pv.X();
    double pvY = pv.Y();
    double pvZ = pv.Z();
    const math::Matrix3x3& pvCov = kineEvt.pvCov();

    const reco::Candidate::LorentzVector& recoilP4 = kineEvtA.recoilP4();
    double recoilPx = recoilP4.px();
    double recoilPy = recoilP4.py();
    double recoilPz = recoilP4.pz();
    double recoilE  = recoilP4.energy();
    const math::Matrix4x4& recoilCov = kineEvt.recoilCov();

    assert(kineEvtA.tauPlusP4_isValid());
    const reco::Candidate::LorentzVector& tauPlusP4 = kineEvtA.tauPlusP4();
    double tauPlusPx = tauPlusP4.px();
    double tauPlusPy = tauPlusP4.py();
    double tauPlusPz = tauPlusP4.pz();
    double tauPlusPt = tauPlusP4.pt();
    double tauPlusP  = tauPlusP4.P();
    double tauPlusE  = tauPlusP4.energy();
    assert(kineEvtA.svTauPlus_isValid());
    const reco::Candidate::Point& svTauPlus = kineEvtA.svTauPlus();
    double svTauPlusX = svTauPlus.X();
    double svTauPlusY = svTauPlus.Y();
    double svTauPlusZ = svTauPlus.Z();
    const math::Matrix3x3& svTauPlusCov = kineEvt.svTauPlusCov();
    auto tauPlusD3 = svTauPlus - pv;
    double tauPlusDx = tauPlusD3.X();
    double tauPlusDy = tauPlusD3.Y();
    double tauPlusDz = tauPlusD3.Z();
    double tauPlusDt = std::sqrt(tauPlusD3.Perp2());
    double tauPlusD  = std::sqrt(tauPlusD3.Mag2());

    const reco::Candidate::LorentzVector& visTauPlusP4 = kineEvtA.visTauPlusP4();
    double visTauPlusPx = visTauPlusP4.px();
    double visTauPlusPy = visTauPlusP4.py();
    double visTauPlusPz = visTauPlusP4.pz();
    double visTauPlusP  = visTauPlusP4.P();
    double visTauPlusE  = visTauPlusP4.energy();

    assert(kineEvtA.nuTauPlusP4_isValid());
    const reco::Candidate::LorentzVector& nuTauPlusP4 = kineEvtA.nuTauPlusP4();
    double nuTauPlusPx = nuTauPlusP4.px();
    double nuTauPlusPy = nuTauPlusP4.py();
    double nuTauPlusPz = nuTauPlusP4.pz();
    double nuTauPlusE  = nuTauPlusP4.energy();
    const math::Matrix3x3& nuTauPlusCov = kineEvt.nuTauPlusCov();

    assert(kineEvtA.tauMinusP4_isValid());
    const reco::Candidate::LorentzVector& tauMinusP4 = kineEvtA.tauMinusP4();
    double tauMinusPx = tauMinusP4.px();
    double tauMinusPy = tauMinusP4.py();
    double tauMinusPz = tauMinusP4.pz();
    double tauMinusPt = tauMinusP4.pt();
    double tauMinusP  = tauMinusP4.P();
    double tauMinusE  = tauMinusP4.energy();
    assert(kineEvtA.svTauMinus_isValid());
    const reco::Candidate::Point& svTauMinus = kineEvtA.svTauMinus();
    double svTauMinusX = svTauMinus.X();
    double svTauMinusY = svTauMinus.Y();
    double svTauMinusZ = svTauMinus.Z();
    const math::Matrix3x3& svTauMinusCov = kineEvt.svTauMinusCov();
    auto tauMinusD3 = svTauMinus - pv;
    double tauMinusDx = tauMinusD3.X();
    double tauMinusDy = tauMinusD3.Y();
    double tauMinusDz = tauMinusD3.Z();
    double tauMinusDt = std::sqrt(tauMinusD3.Perp2());
    double tauMinusD  = std::sqrt(tauMinusD3.Mag2());

    const reco::Candidate::LorentzVector& visTauMinusP4 = kineEvtA.visTauMinusP4();
    double visTauMinusPx = visTauMinusP4.px();
    double visTauMinusPy = visTauMinusP4.py();
    double visTauMinusPz = visTauMinusP4.pz();
    double visTauMinusP  = visTauMinusP4.P();
    double visTauMinusE  = visTauMinusP4.energy();

    assert(kineEvtA.nuTauMinusP4_isValid());
    const reco::Candidate::LorentzVector& nuTauMinusP4 = kineEvtA.nuTauMinusP4();
    double nuTauMinusPx = nuTauMinusP4.px();
    double nuTauMinusPy = nuTauMinusP4.py();
    double nuTauMinusPz = nuTauMinusP4.pz();
    double nuTauMinusE  = nuTauMinusP4.energy();
    const math::Matrix3x3& nuTauMinusCov = kineEvt.nuTauMinusCov();

    //-----------------------------------------------------------------------------------------------
    // CV: the nomenclature of the different matrices and vectors used for the kinematic fit
    //     follows the one introduced in Section 3 of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
    //     For the first iteration, the vector alpha0 refers to the starting-point of the kinematic fit,
    //     computed by the function KinematicFitStartPosFinder::operator()
    //     For subsequent iterations, the vector alpha0 refers to the result alpha obtained by the previous iteration
    //     The symbol H in the comments refers to the constraint equations,
    //     cf. Eq. (3) in https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
    //-----------------------------------------------------------------------------------------------

    // CV: build vector of all measured parameters at point A 
    //     where constrained equations are linearized
    math::VectorP alphaA;
    alphaA( 0) = pvX;
    alphaA( 1) = pvY;
    alphaA( 2) = pvZ;
    alphaA( 3) = nuTauPlusPx;
    alphaA( 4) = nuTauPlusPy;
    alphaA( 5) = nuTauPlusPz;
    alphaA( 6) = svTauPlusX;
    alphaA( 7) = svTauPlusY;
    alphaA( 8) = svTauPlusZ;
    alphaA( 9) = nuTauMinusPx;
    alphaA(10) = nuTauMinusPy;
    alphaA(11) = nuTauMinusPz;
    alphaA(12) = svTauMinusX;
    alphaA(13) = svTauMinusY;
    alphaA(14) = svTauMinusZ;
    alphaA(15) = recoilPx;
    alphaA(16) = recoilPy;
    alphaA(17) = recoilPz;
    alphaA(18) = recoilE;
    if ( verbosity_ >= 1 )
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
    math::MatrixPxP V_alpha0;
    V_alpha0.Place_at(pvCov        ,  0,  0);
    V_alpha0.Place_at(nuTauPlusCov ,  3,  3);
    V_alpha0.Place_at(svTauPlusCov ,  6,  6);
    V_alpha0.Place_at(nuTauMinusCov,  9,  9);
    V_alpha0.Place_at(svTauMinusCov, 12, 12);
    V_alpha0.Place_at(recoilCov    , 15, 15);
    if ( verbosity_ >= 1 )
    {
      printCovMatrix("V_alpha0", V_alpha0);
    }

    int Vinv_alpha0_errorFlag = 0;
    math::MatrixPxP Vinv_alpha0 = V_alpha0.Inverse(Vinv_alpha0_errorFlag);
    if ( Vinv_alpha0_errorFlag != 0 )
    {
      printCovMatrix("V_alpha0", V_alpha0);
      throw cmsException("KinematicFit::operator()", __LINE__)
        << "Failed to invert matrix V_alpha0 !!\n";
    }
    if ( verbosity_ >= 1 )
    {
      std::cout << "Vinv_alpha0:\n";
      std::cout << Vinv_alpha0 << "\n";
    }
  
    math::MatrixCxP D;
    math::VectorC   d;
    // CV: add Higgs mass constraint
    //       H = sqrt(2*mTau^2 
    //              + 2*(visTauPlusE  + nuTauPlusE )*(visTauMinusE  + nuTauMinusE ) 
    //              - 2*(visTauPlusPx + nuTauPlusPx)*(visTauMinusPx + nuTauminusPx)
    //              - 2*(visTauPlusPy + nuTauPlusPy)*(visTauMinusPy + nuTauminusPy)
    //              - 2*(visTauPlusPz + nuTauPlusPz)*(visTauMinusPz + nuTauminusPz)) - mHiggs = 0
    //     where
    //       nuTauPlusE  = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2 )
    //       nuTauMinusE = sqrt(nuTauMisusPx^2 + nuTauMinusPy^2 + nuTauMinusPz^2)
    double denominator0 = (tauPlusP4 + tauMinusP4).mass();
    D( 0, 3) = ((tauPlusPx/tauPlusE)*tauMinusE  - tauMinusPx)/denominator0;              // dH/dnuTauPlusPx
    D( 0, 4) = ((tauPlusPy/tauPlusE)*tauMinusE  - tauMinusPy)/denominator0;              // dH/dnuTauPlusPy
    D( 0, 5) = ((tauPlusPz/tauPlusE)*tauMinusE  - tauMinusPz)/denominator0;              // dH/dnuTauPlusPz
    D( 0, 9) = ((tauMinusPx/tauMinusE)*tauPlusE - tauPlusPx )/denominator0;              // dH/dnuTauMinusPx
    D( 0,10) = ((tauMinusPy/tauMinusE)*tauPlusE - tauPlusPy )/denominator0;              // dH/dnuTauMinusPy
    D( 0,11) = ((tauMinusPz/tauMinusE)*tauPlusE - tauPlusPz )/denominator0;              // dH/dnuTauMinusPz
    d( 0) = (visTauPlusP4 + nuTauPlusP4 + visTauMinusP4 + nuTauMinusP4).mass() - mHiggs;
    // CV: add mass constraint for tau+
    if ( applyTauMassConstraint_ == 0 )
    {
      // CV: "regular" tau mass constraint
      //       H = sqrt(mVis^2 + 2*visTauPlusE*nuTauPlusE 
      //              - 2*visTauPlusPx*nuTauPlusPx - 2*visTauPlusPy*nuTauPlusPy - 2*visTauPlusPz*nuTauPlusPz) - mTau = 0
      //     where
      //       nuTauPlusE = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2)
      double denominator1 = (visTauPlusP4 + nuTauPlusP4).mass();
      D( 1, 3) = ((nuTauPlusPx/nuTauPlusE)*visTauPlusE - visTauPlusPx)/denominator1;     // dH/dnuTauPlusPx
      D( 1, 4) = ((nuTauPlusPy/nuTauPlusE)*visTauPlusE - visTauPlusPy)/denominator1;     // dH/dnuTauPlusPy
      D( 1, 5) = ((nuTauPlusPz/nuTauPlusE)*visTauPlusE - visTauPlusPz)/denominator1;     // dH/dnuTauPlusPz
      d( 1) = (visTauPlusP4 + nuTauPlusP4).mass() - mTau;
    }
    else if ( applyTauMassConstraint_ == 1 )
    {
      // CV: constraint on Gottfried-Jackson angle
      //       H = (((visTauPlusPx + nuTauPlusPx)*visTauPlusPx + (visTauPlusPy + nuTauPlusPy)*visTauPlusPy + (visTauPlusPz + nuTauPlusPz)*visTauPlusPz)/
      //            (sqrt((visTauPlusPx + nuTauPlusPx)^2 + (visTauPlusPy + nuTauPlusPy)^2 + (visTauPlusPz + nuTauPlusPz)^2)*visTauPlusP))
      //          - cos(thetaGJ) = 0
      //     where
      //       cos(thetaGJ) = (sqrt((visTauPlusPx + nuTauPlusPx)^2 + (visTauPlusPy + nuTauPlusPy)^2 + (visTauPlusPz + nuTauPlusPz)^2)
      //                      + visTauPlusP - ((visTauPlusE + nuTauPlusE) +/- visTauPlusE)^2)/
      //                      (2*sqrt((visTauPlusPx + nuTauPlusPx)^2 + (visTauPlusPy + nuTauPlusPy)^2)*visTauPlusP)
      //     is derived in analogy to Eq. (4.1) in Section 4.1.2 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
      //     and
      //       nuTauPlusE = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2)
      double denominator1 = 2.*visTauPlusP*cube(tauPlusP)*tauPlusE;
      D( 1, 3) = (+2.*square(mTau)*visTauPlusE*tauPlusPx - square(visTauPlusE)*tauPlusE*tauPlusPx
                 + tauPlusE*(-square(mTau)*tauPlusPx
                            + nuTauPlusPx*(square(visTauPlusPx) - square(visTauPlusPy) - 2.*nuTauPlusPy*visTauPlusPy
                                                                - square(visTauPlusPz) - 2.*nuTauPlusPz*visTauPlusPz)
                            + visTauPlusPx*(square(visTauPlusP) + 2.*square(nuTauPlusPy) + 2.*nuTauPlusPy*visTauPlusPy
                                                                + 2.*square(nuTauPlusPz) + 2.*nuTauPlusPz*visTauPlusPz)))/denominator1;
      D( 1, 4) = (+2.*square(mTau)*visTauPlusE*tauPlusPy - square(visTauPlusE)*tauPlusE*tauPlusPy
                 + tauPlusE*(-square(mTau)*tauPlusPy
                            + nuTauPlusPy*(square(visTauPlusPy) - square(visTauPlusPz) - 2.*nuTauPlusPx*visTauPlusPx
                                                                - square(visTauPlusPz) - 2.*nuTauPlusPz*visTauPlusPz)
                            + visTauPlusPy*(square(visTauPlusP) + 2.*square(nuTauPlusPx) + 2.*nuTauPlusPx*visTauPlusPx
                                                                + 2.*square(nuTauPlusPz) + 2.*nuTauPlusPz*visTauPlusPz)))/denominator1;
      D( 1, 5) = (+2.*square(mTau)*visTauPlusE*tauPlusPz - square(visTauPlusE)*tauPlusE*tauPlusPz
                 + tauPlusE*(-square(mTau)*tauPlusPz
                            + nuTauPlusPz*(square(visTauPlusPz) - square(visTauPlusPx) - 2.*nuTauPlusPx*visTauPlusPx
                                                                - square(visTauPlusPy) - 2.*nuTauPlusPy*visTauPlusPy)
                            + visTauPlusPz*(square(visTauPlusP) + 2.*square(nuTauPlusPx) + 2.*nuTauPlusPx*visTauPlusPx 
                                                                + 2.*square(nuTauPlusPy) + 2.*nuTauPlusPy*visTauPlusPy)))/denominator1;
      d( 1) = comp_cosThetaGJ(tauPlusP4, visTauPlusP4) - comp_cosThetaGJ_solution(tauPlusP4, visTauPlusP4);
    }
    else assert(0);
    // CV: add "parallelism" constraint for tau+
    //     The constraint equations H_Pphi and H_Ptheta are given by Eqs. (4.10) and (4.11)
    //     in Section 4.1.3.3 of [1].
    //     The small correction for the helix bend between tau lepton production and decay in Eq. (4.13) is not used, 
    //     as this correction is expected to be on the level of 50 microrad, which we consider negligible.
    //     Following [1], we perform a transformation of variables from phi and theta to 1/2*pi - phi and 1/2*pi - theta,
    //     in order to reduce the magnitude of derivatives (by avoiding the division by small numbers)
    //
    //   [1] https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
    math::VectorP D2_1;
    D2_1( 0) =  (1. - tauPlusDx/tauPlusDt)/tauPlusDy;                                      // dH_Pphi/dpvX
    D2_1( 1) = -(tauPlusDx/tauPlusDy)*D2_1(0);                                             // dH_Pphi/dpvY
    D2_1( 3) =  (1. - tauPlusPx/tauPlusPt)/tauPlusPy;                                      // dH_Pphi/dnuTauPlusPx
    D2_1( 4) = -(tauPlusPx/tauPlusPy)*D2_1(3);                                             // dH_Pphi/dnuTauPlusPy     
    D2_1( 6) = -D2_1(0);                                                                   // dH_Pphi/dsvTauPlusX
    D2_1( 7) = -D2_1(1);                                                                   // dH_Pphi/dsvTauPlusY
    double D2_1_mag = comp_mag(D2_1);
    math::VectorP D2_2;
    D2_2( 0) = -(tauPlusDy/tauPlusDx)*(1. - tauPlusDy/tauPlusDt)/tauPlusDx;                // dH_Pphi/dpvX
    D2_2( 1) =  (1. - tauPlusDy/tauPlusDt)/tauPlusDx;                                      // dH_Pphi/dpvY
    D2_2( 3) = -(tauPlusPy/tauPlusPx)*(1. - tauPlusPy/tauPlusPt)/tauPlusPx;                // dH_Pphi/dnuTauPlusPx
    D2_2( 4) =  (1. - tauPlusPy/tauPlusPt)/tauPlusPx;                                      // dH_Pphi/dnuTauPlusPy
    D2_2( 6) = -D2_2(0);                                                                   // dH_Pphi/dsvTauPlusX
    D2_2( 7) = -D2_2(1);                                                                   // dH_Pphi/dsvTauPlusY
    double D2_2_mag = comp_mag(D2_2);
    if ( verbosity_ >= 1 )
    {
      std::cout << "D2 (1st):\n";
      std::cout << D2_1 << "\n";
      std::cout << " mag = " << D2_1_mag << "\n";
      std::cout << "D2 (2nd):\n";
      std::cout << D2_2 << "\n";
      std::cout << " mag = " << D2_2_mag << "\n";
    }
    if ( D2_1_mag < D2_2_mag )
    {
      for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      {
        D( 2,idx) = D2_1(idx);
      }
      d( 2) = (tauPlusDt - tauPlusDx)/tauPlusDy + (tauPlusPx - tauPlusPt)/tauPlusPy;
    }
    else
    {
      for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      {
        D( 2,idx) = D2_2(idx);
      }
      d( 2) = (tauPlusDt - tauPlusDy)/tauPlusDx + (tauPlusPy - tauPlusPt)/tauPlusPx;
    }
    math::VectorP D3_1;
    D3_1( 0) =  (tauPlusDx/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);                        // dH_Ptheta/dpvX
    D3_1( 1) =  (tauPlusDy/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);                        // dH_Ptheta/dpvY
    D3_1( 2) = -1./tauPlusD + (tauPlusD - tauPlusDt)/square(tauPlusDz);                    // dH_Ptheta/dpvZ
    D3_1( 3) =  (tauPlusPx/tauPlusPz)*(1./tauPlusPt - 1./tauPlusP);                        // dH_Ptheta/dnuTauPlusPx
    D3_1( 4) =  (tauPlusPy/tauPlusPz)*(1./tauPlusPt - 1./tauPlusP);                        // dH_Ptheta/dnuTauPlusPy
    D3_1( 5) = -1./tauPlusP + (tauPlusP - tauPlusPt)/square(tauPlusPz);                    // dH_Ptheta/dnuTauPlusPz      
    D3_1( 6) = -D3_1(0);                                                                   // dH_Ptheta/dsvTauPlusX
    D3_1( 7) = -D3_1(1);                                                                   // dH_Ptheta/dsvTauPlusY
    D3_1( 8) = -D3_1(2);                                                                   // dH_Ptheta/dsvTauPlusZ
    double D3_1_mag = comp_mag(D3_1);
    math::VectorP D3_2;
    D3_2( 0) = -(tauPlusDx*tauPlusDz/cube(tauPlusDt))*(1. - tauPlusDz/tauPlusD);           // dH_Ptheta/dpvX
    D3_2( 1) = -(tauPlusDy*tauPlusDz/cube(tauPlusDt))*(1. - tauPlusDz/tauPlusD);           // dH_Ptheta/dpvY
    D3_2( 2) =  (1. - tauPlusDz/tauPlusD)/tauPlusDt;                                       // dH_Ptheta/dpvZ
    D3_2( 3) = -(tauPlusPx*tauPlusPz/cube(tauPlusPt))*(1. - tauPlusPz/tauPlusP);           // dH_Ptheta/dnuTauPlusPx
    D3_2( 4) = -(tauPlusPy*tauPlusPz/cube(tauPlusPt))*(1. - tauPlusPz/tauPlusP);           // dH_Ptheta/dnuTauPlusPy
    D3_2( 5) =  (1. - tauPlusPz/tauPlusP)/tauPlusPt;                                       // dH_Ptheta/dnuTauPlusPz
    D3_2( 6) = -D3_2(0);                                                                   // dH_Ptheta/dsvTauPlusX
    D3_2( 7) = -D3_2(1);                                                                   // dH_Ptheta/dsvTauPlusY
    D3_2( 8) = -D3_2(2);                                                                   // dH_Ptheta/dsvTauPlusZ
    double D3_2_mag = comp_mag(D3_2);
    if ( verbosity_ >= 1 )
    {
      std::cout << "D3 (1st):\n";
      std::cout << D3_1 << "\n";
      std::cout << " mag = " << D3_1_mag << "\n";
      std::cout << "D3 (2nd):\n";
      std::cout << D3_2 << "\n";
      std::cout << " mag = " << D3_2_mag << "\n";
    }
    if ( D3_1_mag < D3_2_mag )
    {
      for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      {
        D( 3,idx) = D3_1(idx);
      }
      d( 3) = (tauPlusD - tauPlusDt)/tauPlusDz + (tauPlusPt - tauPlusP)/tauPlusPz;
    }
    else
    {
      for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      {
        D( 3,idx) = D3_2(idx);
      }
      d( 3) = (tauPlusD - tauPlusDz)/tauPlusDt + (tauPlusPz - tauPlusP)/tauPlusPt;
    }
    // CV: add mass constraint for tau-
    if ( applyTauMassConstraint_ == 0 )
    {
      // CV: "regular" tau mass constraint
      //       H = sqrt(mVis^2 + 2*visTauMinusE*nuTauMinusE 
      //              - 2*visTauMinusPx*nuTauMinusPx - 2*visTauMinusPy*nuTauMinusPy - 2*visTauMinusPz*nuTauMinusPz) - mTau = 0
      //     where
      //       nuTauMinusE = sqrt(nuTauMinusPx^2  + nuTauMinusPy^2  + nuTauMinusPz^2)
      double denominator1 = (visTauMinusP4 + nuTauMinusP4).mass();
      D( 4, 9) = ((nuTauMinusPx/nuTauMinusE)*visTauMinusE - visTauMinusPx)/denominator1;   // dH/dnuTauMinusPx
      D( 4,10) = ((nuTauMinusPy/nuTauMinusE)*visTauMinusE - visTauMinusPy)/denominator1;   // dH/dnuTauMinusPy
      D( 4,11) = ((nuTauMinusPz/nuTauMinusE)*visTauMinusE - visTauMinusPz)/denominator1;   // dH/dnuTauMinusPz
      d( 4) = (visTauMinusP4 + nuTauMinusP4).mass() - mTau;
    }
    else if ( applyTauMassConstraint_ == 1 )
    {
      // CV: constraint on Gottfried-Jackson angle
      //       H = (((visTauMinusPx + nuTauMinusPx)*visTauMinusPx + (visTauMinusPy + nuTauMinusPy)*visTauMinusPy + (visTauMinusPz + nuTauMinusPz)*visTauMinusPz)/
      //            (sqrt((visTauMinusPx + nuTauMinusPx)^2 + (visTauMinusPy + nuTauMinusPy)^2 + (visTauMinusPz + nuTauMinusPz)^2)*visTauMinusP))
      //          - cos(thetaGJ) = 0
      //     where
      //       cos(thetaGJ) = (sqrt((visTauMinusPx + nuTauMinusPx)^2 + (visTauMinusPy + nuTauMinusPy)^2 + (visTauMinusPz + nuTauMinusPz)^2)
      //                      + visTauMinusP - ((visTauMinusE + nuTauMinusE) +/- visTauMinusE)^2)/
      //                      (2*sqrt((visTauMinusPx + nuTauMinusPx)^2 + (visTauMinusPy + nuTauMinusPy)^2)*visTauMinusP)
      //     is derived in analogy to Eq. (4.1) in Section 4.1.2 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
      //     and
      //       nuTauMinusE = sqrt(nuTauMinusPx^2  + nuTauMinusPy^2  + nuTauMinusPz^2)
      double denominator1 = 2.*visTauMinusP*cube(tauMinusP)*tauMinusE;
      D( 4, 9) = (+2.*square(mTau)*visTauMinusE*tauMinusPx - square(visTauMinusE)*tauMinusE*tauMinusPx
                 + tauMinusE*(-square(mTau)*tauMinusPx
                             + nuTauMinusPx*(square(visTauMinusPx) - square(visTauMinusPy) - 2.*nuTauMinusPy*visTauMinusPy
                                                                   - square(visTauMinusPz) - 2.*nuTauMinusPz*visTauMinusPz)
                             + visTauMinusPx*(square(visTauMinusP) + 2.*square(nuTauMinusPy) + 2.*nuTauMinusPy*visTauMinusPy
                                                                   + 2.*square(nuTauMinusPz) + 2.*nuTauMinusPz*visTauMinusPz)))/denominator1;
      D( 4,10) = (+2.*square(mTau)*visTauMinusE*tauMinusPy - square(visTauMinusE)*tauMinusE*tauMinusPy
                 + tauMinusE*(-square(mTau)*tauMinusPy
                             + nuTauMinusPy*(square(visTauMinusPy) - square(visTauMinusPz) - 2.*nuTauMinusPx*visTauMinusPx
                                                                   - square(visTauMinusPz) - 2.*nuTauMinusPz*visTauMinusPz)
                             + visTauMinusPy*(square(visTauMinusP) + 2.*square(nuTauMinusPx) + 2.*nuTauMinusPx*visTauMinusPx
                                                                   + 2.*square(nuTauMinusPz) + 2.*nuTauMinusPz*visTauMinusPz)))/denominator1;
      D( 4,11) = (+2.*square(mTau)*visTauMinusE*tauMinusPz - square(visTauMinusE)*tauMinusE*tauMinusPz
                 + tauMinusE*(-square(mTau)*tauMinusPz
                             + nuTauMinusPz*(square(visTauMinusPz) - square(visTauMinusPx) - 2.*nuTauMinusPx*visTauMinusPx
                                                                   - square(visTauMinusPy) - 2.*nuTauMinusPy*visTauMinusPy)
                             + visTauMinusPz*(square(visTauMinusP) + 2.*square(nuTauMinusPx) + 2.*nuTauMinusPx*visTauMinusPx 
                                                                   + 2.*square(nuTauMinusPy) + 2.*nuTauMinusPy*visTauMinusPy)))/denominator1;
      d( 4) = comp_cosThetaGJ(tauMinusP4, visTauMinusP4) - comp_cosThetaGJ_solution(tauMinusP4, visTauMinusP4);
    }
    else assert(0);
    // CV: add "parallelism" constraint for tau-
    //     The constraint equations H_Pphi and H_Ptheta are given by Eqs. (4.10) and (4.11)
    //     in Section 4.1.3.3 of [1].
    //     The small correction for the helix bend between tau lepton production and decay in Eq. (4.13) is not used, 
    //     as this correction is expected to be on the level of 50 microrad, which we consider negligible.
    //     Following [1], we perform a transformation of variables from phi and theta to 1/2*pi - phi and 1/2*pi - theta,
    //     in order to reduce the magnitude of derivatives (by avoiding the division by small numbers)
    //
    //   [1] https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
    math::VectorP D5_1;
    D5_1( 0) =  (1. - tauMinusDx/tauMinusDt)/tauMinusDy;                                   // dH_Pphi/dpvX
    D5_1( 1) = -(tauMinusDx/tauMinusDy)*D5_1(0);                                           // dH_Pphi/dpvY
    D5_1( 9) =  (1. - tauMinusPx/tauMinusPt)/tauMinusPy;                                   // dH_Pphi/dnuTauMinusPx
    D5_1(10) = -(tauMinusPx/tauMinusPy)*D5_1(9);                                           // dH_Pphi/dnuTauMinusPy      
    D5_1(12) = -D5_1(0);                                                                   // dH_Pphi/dsvTauMinusX
    D5_1(13) = -D5_1(1);                                                                   // dH_Pphi/dsvTauMinusY
    double D5_1_mag = comp_mag(D5_1);
    math::VectorP D5_2;
    D5_2( 0) = -(tauMinusDy/tauMinusDx)*(1. - tauMinusDy/tauMinusDt)/tauMinusDx;           // dH_Pphi/dpvX
    D5_2( 1) =  (1. - tauMinusDy/tauMinusDt)/tauMinusDx;                                   // dH_Pphi/dpvY
    D5_2( 9) = -(tauMinusPy/tauMinusPx)*(1. - tauMinusPy/tauMinusPt)/tauMinusPx;           // dH_Pphi/dnuTauMinusPx
    D5_2(10) =  (1. - tauMinusPy/tauMinusPt)/tauMinusPx;                                   // dH_Pphi/dnuTauMinusPy
    D5_2(12) = -D5_2(0);                                                                   // dH_Pphi/dsvTauMinusX
    D5_2(13) = -D5_2(1);                                                                   // dH_Pphi/dsvTauMinusY
    double D5_2_mag = comp_mag(D5_2);
    if ( verbosity_ >= 1 )
    {
      std::cout << "D5 (1st):\n";
      std::cout << D5_1 << "\n";
      std::cout << " mag = " << D5_1_mag << "\n";
      std::cout << "D5 (2nd):\n";
      std::cout << D5_2 << "\n";
      std::cout << " mag = " << D5_2_mag << "\n";
    }
    if ( D5_1_mag < D5_2_mag )
    {
      for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      {
        D( 5,idx) = D5_1(idx);
      }
      d( 5) =  (tauMinusDt - tauMinusDx)/tauMinusDy + (tauMinusPx - tauMinusPt)/tauMinusPy;
    }
    else
    {
      for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      {
        D( 5,idx) = D5_2(idx);
      }
      d( 5) =  (tauMinusDt - tauMinusDy)/tauMinusDx + (tauMinusPy - tauMinusPt)/tauMinusPx;
    }
    math::VectorP D6_1;
    D6_1( 0) =  (tauMinusDx/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);                    // dH_Ptheta/dpvX
    D6_1( 1) =  (tauMinusDy/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);                    // dH_Ptheta/dpvY
    D6_1( 2) = -1./tauMinusD + (tauMinusD - tauMinusDt)/square(tauMinusDz);                // dH_Ptheta/dpvZ
    D6_1( 9) =  (tauMinusPx/tauMinusPz)*(1./tauMinusPt - 1./tauMinusP);                    // dH_Ptheta/dnuTauMinusPx
    D6_1(10) =  (tauMinusPy/tauMinusPz)*(1./tauMinusPt - 1./tauMinusP);                    // dH_Ptheta/dnuTauMinusPy
    D6_1(11) = -1./tauMinusP + (tauMinusP - tauMinusPt)/square(tauMinusPz);                // dH_Ptheta/dnuTauMinusPz      
    D6_1(12) = -D6_1(0);                                                                   // dH_Ptheta/dsvTauMinusX
    D6_1(13) = -D6_1(1);                                                                   // dH_Ptheta/dsvTauMinusY
    D6_1(14) = -D6_1(2);                                                                   // dH_Ptheta/dsvTauMinusZ
    double D6_1_mag = comp_mag(D6_1);
    math::VectorP D6_2;
    D6_2( 0) = -(tauMinusDx*tauMinusDz/cube(tauMinusDt))*(1. - tauMinusDz/tauMinusD);      // dH_Ptheta/dpvX
    D6_2( 1) = -(tauMinusDy*tauMinusDz/cube(tauMinusDt))*(1. - tauMinusDz/tauMinusD);      // dH_Ptheta/dpvY
    D6_2( 2) =  (1. - tauMinusDz/tauMinusD)/tauMinusDt;                                    // dH_Ptheta/dpvZ
    D6_2( 9) = -(tauMinusPx*tauMinusPz/cube(tauMinusPt))*(1. - tauMinusPz/tauMinusP);      // dH_Ptheta/dnuTauMinusPx
    D6_2(10) = -(tauMinusPy*tauMinusPz/cube(tauMinusPt))*(1. - tauMinusPz/tauMinusP);      // dH_Ptheta/dnuTauMinusPy
    D6_2(11) =  (1. - tauMinusPz/tauMinusP)/tauMinusPt;                                    // dH_Ptheta/dnuTauMinusPz
    D6_2(12) = -D6_2(0);                                                                   // dH_Ptheta/dsvTauMinusX
    D6_2(13) = -D6_2(1);                                                                   // dH_Ptheta/dsvTauMinusY
    D6_2(14) = -D6_2(2);                                                                   // dH_Ptheta/dsvTauMinusZ
    double D6_2_mag = comp_mag(D6_2);
    if ( verbosity_ >= 1 )
    {
      std::cout << "D6 (1st):\n";
      std::cout << D6_1 << "\n";
      std::cout << " mag = " << D6_1_mag << "\n";
      std::cout << "D6 (2nd):\n";
      std::cout << D6_2 << "\n";
      std::cout << " mag = " << D6_2_mag << "\n";
    }
    if ( D6_1_mag < D6_2_mag )
    {
      for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      {
        D( 6,idx) = D6_1(idx);
      }
      d( 6) =  (tauMinusD - tauMinusDt)/tauMinusDz + (tauMinusPt - tauMinusP)/tauMinusPz;
    }
    else
    {
      for ( int idx = 0; idx < kinFit::numParameters; ++idx )
      {
        D( 6,idx) = D6_2(idx);
      }
      d( 6) =  (tauMinusD - tauMinusDz)/tauMinusDt + (tauMinusPz - tauMinusP)/tauMinusPt;
    }
    // CV: add constraint that recoil = tau+ + tau-
    //       H_px     = nuTauPlusPx + visTauPlusPx + nuTauMinusPx + visTauMinusPx - recoilPx = 0
    //       H_py     = nuTauPlusPy + visTauPlusPy + nuTauMinusPy + visTauMinusPy - recoilPy = 0
    //       H_pz     = nuTauPlusPz + visTauPlusPz + nuTauMinusPz + visTauMinusPz - recoilPz = 0
    //       H_energy = nuTauPlusE  + visTauPlusE  + nuTauMinusE  + visTauMinusE  - recoilE  = 0
    //     where:
    //       nuTauPlusE  = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2 )
    //       nuTauMinusE = sqrt(nuTauMisusPx^2 + nuTauMinusPy^2 + nuTauMinusPz^2)
    D( 7, 3) = +1.;                                                                      // dH_px/dnuTauPlusPx
    D( 7, 9) = +1.;                                                                      // dH_px/dnuTauMinusPx
    D( 7,15) = -1.;                                                                      // dH_px/drecoilPx
    d( 7) = visTauPlusPx + nuTauPlusPx + visTauMinusPx + nuTauMinusPx - recoilPx;
    D( 8, 4) = +1.;                                                                      // dH_py/dnuTauPlusPy
    D( 8,10) = +1.;                                                                      // dH_py/dnuTauMinusPy
    D( 8,16) = -1.;                                                                      // dH_py/drecoilPy
    d( 8) = visTauPlusPy + nuTauPlusPy + visTauMinusPy + nuTauMinusPy - recoilPy;
    D( 9, 5) = +1.;                                                                      // dH_pz/dnuTauPlusPz
    D( 9,11) = +1.;                                                                      // dH_pz/dnuTauMinusPz
    D( 9,17) = -1.;                                                                      // dH_pz/drecoilPz
    d( 9) = visTauPlusPz + nuTauPlusPz + visTauMinusPz + nuTauMinusPz - recoilPz;
    D(10, 3) = tauPlusPx/tauPlusE;                                                       // dH_energy/dnuTauPlusPx
    D(10, 4) = tauPlusPy/tauPlusE;                                                       // dH_energy/dnuTauPlusPy
    D(10, 5) = tauPlusPz/tauPlusE;                                                       // dH_energy/dnuTauPlusPz
    D(10, 9) = tauMinusPx/tauMinusE;                                                     // dH_energy/dnuTauPlusPx
    D(10,10) = tauMinusPy/tauMinusE;                                                     // dH_energy/dnuTauPlusPy
    D(10,11) = tauMinusPz/tauMinusE;                                                     // dH_energy/dnuTauPlusPz
    D(10,18) = -1.;                                                                      // dH_energy/drecoilE
    d(10) = visTauPlusE  + nuTauPlusE  + visTauMinusE  + nuTauMinusE  - recoilE;
    if ( verbosity_ >= 1 )
    {
      std::cout << "D:\n";
      std::cout << D << "\n";
      std::cout << "d:\n";
      std::cout << d << "\n";
    }

    math::VectorP p;
    if ( applyLifetimeConstraint_ )
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
      p( 3) = -term2TauPlus*tauPlusPx;                                                   // d(-2*ln(P))/dnuTauPlusPx
      p( 4) = -term2TauPlus*tauPlusPy;                                                   // d(-2*ln(P))/dnuTauPlusPy
      p( 5) = -term2TauPlus*tauPlusPz;                                                   // d(-2*ln(P))/dnuTauPlusPz
      p( 6) =  term1TauPlus*tauPlusDx;                                                   // d(-2*ln(P))/dsvTauPlusX
      p( 7) =  term1TauPlus*tauPlusDy;                                                   // d(-2*ln(P))/dsvTauPlusY
      p( 8) =  term1TauPlus*tauPlusDz;                                                   // d(-2*ln(P))/dsvTauPlusZ
      p( 9) = -term2TauMinus*tauMinusPx;                                                 // d(-2*ln(P))/dnuTauMinusPx
      p(10) = -term2TauMinus*tauMinusPy;                                                 // d(-2*ln(P))/dnuTauMinusPy
      p(11) = -term2TauMinus*tauMinusPz;                                                 // d(-2*ln(P))/dnuTauMinusPz
      p(12) =  term1TauMinus*tauMinusDx;                                                 // d(-2*ln(P))/dsvTauMinusX
      p(13) =  term1TauMinus*tauMinusDy;                                                 // d(-2*ln(P))/dsvTauMinusY
      p(14) =  term1TauMinus*tauMinusDz;                                                 // d(-2*ln(P))/dsvTauMinusZ
      if ( verbosity_ >= 1 )
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

    math::MatrixPxC DT = ROOT::Math::Transpose(D);    
    if ( verbosity_ >= 1 )
    {
      std::cout << "DT:\n";
      std::cout << DT << "\n";
    }

    math::MatrixCxC Vinv_D = D*V_alpha0*DT;
    int V_D_errorFlag = 0;
    math::MatrixCxC V_D = Vinv_D.Inverse(V_D_errorFlag);
    if ( V_D_errorFlag != 0 )
    {
      printCovMatrix("Vinv_D", Vinv_D);
      throw cmsException("KinematicFit::operator()", __LINE__)
        << "Failed to invert matrix Vinv_D !!\n";
    }
    if ( verbosity_ >= 1 )
    {
      printCovMatrix("V_D", V_D);
    }

    math::VectorP dalpha0 = alpha0 - alphaA;
    math::VectorC lambda = V_D*(D*dalpha0 + d);
    if ( applyLifetimeConstraint_ )
    {
      lambda -= V_D*D*V_alpha0*p;
    }
    if ( verbosity_ >= 1 )
    {
      std::cout << "dalpha0:\n";
      std::cout << dalpha0 << "\n";
      std::cout << "D*dalpha0:\n";
      std::cout << D*dalpha0 << "\n";
      std::cout << "lambda:\n";
      std::cout << lambda << "\n";
      std::cout << "contributions to lambda:\n";
      math::VectorC D_times_dalpha0 = D*dalpha0;
      for ( int idxRow = 0; idxRow < kinFit::numConstraints; ++idxRow )
      {
        std::cout << "lambda[" << idxRow << "]:";
        for ( int idxColumn = 0; idxColumn < kinFit::numConstraints; ++idxColumn )
        {
          if ( idxColumn != 0 ) std::cout << ",";
          std::cout << " " << V_D(idxRow,idxColumn)*D_times_dalpha0(idxColumn) << " + " << V_D(idxRow,idxColumn)*d(idxColumn);
        }
        std::cout << "\n";
      }
      std::cout << "DT*lambda:\n";
      std::cout << DT*lambda << "\n";
      std::cout << "distance from satisfaction:\n";
      for ( int idx = 0; idx < kinFit::numConstraints; ++idx )
      {
        double pull = (D_times_dalpha0(idx) + d(idx))/std::sqrt(Vinv_D(idx,idx));
        std::cout << "pull[" << idx << "] = " << pull << "\n";
      }
    }

    math::VectorP alpha = alpha0 - V_alpha0*DT*lambda;
    if ( applyLifetimeConstraint_ )
    {
      alpha -= V_alpha0*p;
    }
    if ( verbosity_ >= 1 )
    {
      std::cout << "alpha0:\n";
      std::cout << alpha0 << "\n";
      std::cout << "alpha:\n";
      std::cout << alpha << "\n";
    }
 
    math::MatrixPxP V_alpha = V_alpha0 - V_alpha0*DT*V_D*D*V_alpha0;
    if ( verbosity_ >= 1 )
    {
      printCovMatrix("V_alpha", V_alpha);
    }

    math::VectorP dalpha = alpha - alphaA;
    double dalpha_mag = std::sqrt(ROOT::Math::Dot(dalpha, dalpha));
    if ( verbosity_ >= 1 )
    {
      std::cout << "dalpha:\n";
      std::cout << dalpha << "\n";
      std::cout << "|alpha - alphaA| = " << dalpha_mag << "\n";
    }

    // CV: compute chi^2
    math::VectorP alpha_minus_alpha0 = alpha - alpha0;
    double chi2 = ROOT::Math::Dot(alpha_minus_alpha0, Vinv_alpha0*alpha_minus_alpha0) + ROOT::Math::Dot(lambda, Vinv_D*lambda);
    if ( applyLifetimeConstraint_ )
    {
      chi2 += 2.*tauPlusD*mTau/(ct*tauPlusP) + 2.*tauMinusD*mTau/(ct*tauMinusP);
    }
    chi2 /= kinFit::numParameters;
    if ( verbosity_ >= 1 )
    {
      std::cout << "chi^2 = " << chi2 << "\n";
    }

    // CV: compute residuals
    math::VectorC residuals = D*dalpha + d;
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

    if ( verbosity_ >= 1 )
    {
      math::VectorP pulls;
      for ( int idx = 0; idx < numParameters; ++idx )
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
    kineEvtA.pvCov_ = V_alpha.Sub<math::Matrix3x3>(0,0);
    //kineEvtA.recoilP4_ = fixHiggsMass(reco::Candidate::LorentzVector(alpha(15), alpha(16), alpha(17), alpha(18)));
    kineEvtA.recoilP4_ = reco::Candidate::LorentzVector(alpha(15), alpha(16), alpha(17), alpha(18));
    kineEvtA.recoilCov_ = V_alpha.Sub<math::Matrix4x4>(15,15);
    kineEvtA.nuTauPlusP4_ = build_nuP4(alpha(3), alpha(4), alpha(5));
    kineEvtA.nuTauPlusP4_isValid_ = true;
    kineEvtA.nuTauPlusCov_ = V_alpha.Sub<math::Matrix3x3>(3,3);
    //kineEvtA.tauPlusP4_ = fixTauMass(kineEvtA.visTauPlusP4_ + kineEvtA.nuTauPlusP4_);    
    kineEvtA.tauPlusP4_ = kineEvtA.visTauPlusP4_ + kineEvtA.nuTauPlusP4_;
    kineEvtA.tauPlusP4_isValid_ = true;
    kineEvtA.svTauPlus_ = reco::Candidate::Point(alpha(6), alpha(7), alpha(8));
    kineEvtA.svTauPlusCov_ = V_alpha.Sub<math::Matrix3x3>(6,6);
    kineEvtA.nuTauMinusP4_ = build_nuP4(alpha(9), alpha(10), alpha(11));
    kineEvtA.nuTauMinusP4_isValid_ = true;
    kineEvtA.nuTauMinusCov_ = V_alpha.Sub<math::Matrix3x3>(9,9);
    //kineEvtA.tauMinusP4_ = fixTauMass(kineEvtA.visTauMinusP4_ + kineEvtA.nuTauMinusP4_);
    kineEvtA.tauMinusP4_ = kineEvtA.visTauMinusP4_ + kineEvtA.nuTauMinusP4_;
    kineEvtA.tauMinusP4_isValid_ = true;
    kineEvtA.svTauMinus_ = reco::Candidate::Point(alpha(12), alpha(13), alpha(14));
    kineEvtA.svTauMinusCov_ = V_alpha.Sub<math::Matrix3x3>(12,12);
    if ( verbosity_ >= 1 )
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
      reco::Candidate::LorentzVector recoilP4(alpha(15), alpha(16), alpha(17), alpha(18));
      printLorentzVector("recoil", recoilP4);
      std::cout << " mass = " << recoilP4.mass() << "\n";
      std::cout << "Higgs - recoil:"
                << " dE = "  << higgsP4.energy() - recoilP4.energy() << ","
                << " dPx = " << higgsP4.px()     - recoilP4.px()     << ","
                << " dPy = " << higgsP4.py()     - recoilP4.py()     << ","
                << " dPz = " << higgsP4.pz()     - recoilP4.pz()     << "\n";
    }

    if ( verbosity_ >= 1 )
    {
      std::cout << "isNaN = " << isNaN(kineEvtA, verbosity_) << "\n";
      std::cout << "isPhysicalSolution = " << isPhysicalSolution(kineEvtA, verbosity_) << "\n";
    }

    if ( !isNaN(kineEvtA) && isPhysicalSolution(kineEvtA) )
    {
      if ( status == -1 || chi2 < min_chi2 )
      {
        kineEvt_kinfit = kineEvtA;
        kineEvt_kinfit.kinFitCov_ = V_alpha;
        kineEvt_kinfit.kinFitChi2_ = chi2;
        kineEvt_kinfit.kinFit_isValid_ = true;
        min_chi2 = chi2;
        status = 0;
      }
      if ( status == 0 && dalpha_mag < 1.e-1 )
      {
        status = 1;
        hasConverged = true;
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
    std::cout << "#iterations = " << iteration << " (max_iterations = " << max_iterations << ")\n";
    std::cout << " status = " << status << "\n";
    std::cout << " min(chi^2) = " << min_chi2 << "\n";
  }

  kineEvt_kinfit.kinFitStatus_ = status;
  kineEvt_kinfit.recoilP4_ = kineEvt_kinfit.tauPlusP4_ + kineEvt_kinfit.tauMinusP4_;
  if ( hasConverged )
  {
    reco::Candidate::Vector hPlus = spinAnalyzer_(kineEvt_kinfit, SpinAnalyzerBase::kTauPlus);
    kineEvt_kinfit.hPlus_ = hPlus;
    kineEvt_kinfit.hPlus_isValid_ = true;
    reco::Candidate::Vector hMinus = spinAnalyzer_(kineEvt_kinfit, SpinAnalyzerBase::kTauMinus);
    kineEvt_kinfit.hMinus_ = hMinus;
    kineEvt_kinfit.hMinus_isValid_ = true;
  }

  return kineEvt_kinfit;
}

