#include "TauAnalysis/Entanglement/interface/KinematicFit.h"

#include "FWCore/Utilities/interface/Exception.h"                         // cms::Exception

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
#include "TauAnalysis/Entanglement/interface/rotateCovMatrix.h"           // rotateCovMatrix()
#include "TauAnalysis/Entanglement/interface/rotateVector.h"              // rotateVector()
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
  double
  get_sign(const reco::Candidate::LorentzVector& visP4, const reco::Candidate::LorentzVector& nuP4, bool& errorFlag, int verbosity = -1)
  {
    const double signPlus  = +1.;
    const double signMinus = -1.;
    double dummy = 0.;
    errorFlag = false;
    double nuPz_signPlus  = comp_nuPz(visP4, nuP4.px(), nuP4.py(), signPlus,  dummy, dummy, errorFlag, verbosity);
    double nuPz_signMinus = comp_nuPz(visP4, nuP4.px(), nuP4.py(), signMinus, dummy, dummy, errorFlag, verbosity);
    if ( std::fabs(nuPz_signPlus - nuP4.pz()) < std::fabs(nuPz_signMinus - nuP4.pz()) )
    {
      if ( std::fabs(nuPz_signPlus  - nuP4.pz()) > 1.e-1*std::max(1., nuP4.pz()) ) errorFlag = true; 
      return signPlus;
    }
    else
    {
      if ( std::fabs(nuPz_signMinus - nuP4.pz()) > 1.e-1*std::max(1., nuP4.pz()) ) errorFlag = true;
      return signMinus;
    }
  }

  reco::Candidate::Point
  get_svStartPos(const reco::Candidate::Point& pv, const reco::Candidate::Point& sv, const math::Matrix3x3& rotMatrix_xyz2rnk, int verbosity)
  {
    auto flightlength_xyz = sv - pv;
    auto flightlength_rnk = rotateVector(flightlength_xyz, rotMatrix_xyz2rnk);
    // CV: set flightlength to at least 10 micrometer
    //     in direction of visible decay products
    const double epsilon = 1.e-3;
    if ( flightlength_rnk.z() < epsilon )
    {
      double sf = epsilon/flightlength_rnk.z();
      if ( verbosity >= 3 )
      {
        std::cout << "WARNING: Component k of flightlength vector small or negative "
                  << " -> scaling flightlength vector by factor = " << sf << " !!\n";
        std::cout << "flightlength (before scaling): r = " << flightlength_rnk.x() << "," 
                  << " n = " << flightlength_rnk.y() << ", k = " << flightlength_rnk.z() << "\n";
        std::cout << "flightlength (after  scaling): r = " << sf*flightlength_rnk.x() << "," 
                  << " n = " << sf*flightlength_rnk.y() << ", k = " << sf*flightlength_rnk.z() << "\n";
      }
      flightlength_rnk *= sf;
    }
    return reco::Candidate::Point(flightlength_rnk.x(), flightlength_rnk.y(), flightlength_rnk.z());
  }

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

  bool signTauPlus_errorFlag = false;
  double signTauPlus = get_sign(kineEvt.visTauPlusP4(), kineEvt.nuTauPlusP4(), signTauPlus_errorFlag, verbosity_);
  bool signTauMinus_errorFlag = false;
  double signTauMinus = get_sign(kineEvt.visTauMinusP4(), kineEvt.nuTauMinusP4(), signTauMinus_errorFlag, verbosity_);
  if ( signTauPlus_errorFlag || signTauMinus_errorFlag )
  {
    std::cerr << "WARNING: Failed to reproduce the neutrino Pz of the start position !!\n";
    return kineEvt_kinFit;
  }

  if ( verbosity_ >= 1 )
  {
    std::cout << "sign(tau+) = " << signTauPlus << ", sign(tau-) = " << signTauMinus << "\n";
  }

  const reco::Candidate::Point& startPos_pv = kineEvt.pv();
  const Matrix3x3& pvCov = kineEvt.pvCov();
  const reco::Candidate::LorentzVector& startPos_visTauPlusP4 = kineEvt.visTauPlusP4();
  const reco::Candidate::LorentzVector& startPos_nuTauPlusP4 = kineEvt.nuTauPlusP4();
  Matrix2x2 nuTauPlusCov = kineEvt.nuTauPlusCov().Sub<Matrix2x2>(0,0);
  // CV: tau+ decay vertex is represented by flightlength = sv - pv in the kinematic fit
  //     and given in helicity-frame coordinates { r, n, k }
  reco::Candidate::Point startPos_svTauPlus = get_svStartPos(kineEvt.pv(), kineEvt.svTauPlus(), kineEvt.tauPlus_rotMatrix_xyz2rnk(), verbosity_);
  const Matrix3x3& svTauPlusCov = rotateCovMatrix(kineEvt.svTauPlusCov(), kineEvt.tauPlus_rotMatrix_xyz2rnk());
  const reco::Candidate::LorentzVector& startPos_visTauMinusP4 = kineEvt.visTauMinusP4();
  const reco::Candidate::LorentzVector& startPos_nuTauMinusP4 = kineEvt.nuTauMinusP4();
  Matrix2x2 nuTauMinusCov = kineEvt.nuTauMinusCov().Sub<Matrix2x2>(0,0);
  // CV: tau- decay vertex is represented by flightlength = sv - pv in the kinematic fit
  //     and given in helicity-frame coordinates { r, n, k }
  reco::Candidate::Point startPos_svTauMinus = get_svStartPos(kineEvt.pv(), kineEvt.svTauMinus(), kineEvt.tauMinus_rotMatrix_xyz2rnk(), verbosity_);
  const Matrix3x3& svTauMinusCov = rotateCovMatrix(kineEvt.svTauMinusCov(), kineEvt.tauMinus_rotMatrix_xyz2rnk());
  const reco::Candidate::LorentzVector& startPos_recoilP4 = kineEvt.recoilP4();
  const Matrix4x4& recoilCov = kineEvt.recoilCov();

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
  if ( verbosity_ >= 2 )
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
  if ( verbosity_ >= 2 )
  {
    printCovMatrix("V_alpha0", V_alpha0);
  }

  VectorP alphaA = alpha0;
  if ( verbosity_ >= 2 )
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
    cfg_kinFit.addParameter<int>("verbosity", verbosity_);
    KinFitAlgo<P,C> kinFit(cfg_kinFit);        
    try
    { 
      fitResult = kinFit(alpha0, V_alpha0, constraint, alphaA);
    }
    catch ( const cms::Exception& )
    {
      std::cerr << "Error in <KinematicFit::operator()>: Caught exception from KinFitAlgo !!\n";
      return kineEvt_kinFit;
    }
  }
  else if ( collider_ == kSuperKEKB )
  {
    const int C = numConstraints_SuperKEKB;
    KinFitConstraint<P,C> constraint(collider_, kineEvt, signTauPlus, signTauMinus, verbosity_);
    edm::ParameterSet cfg_kinFit;
    cfg_kinFit.addParameter<int>("verbosity", verbosity_);
    KinFitAlgo<P,C> kinFit(cfg_kinFit);
    try
    { 
      fitResult = kinFit(alpha0, V_alpha0, constraint, alphaA);
    }
    catch ( const cms::Exception& )
    {
      std::cerr << "Error in <KinematicFit::operator()>: Caught exception from KinFitAlgo !!\n";
      return kineEvt_kinFit;
    }
  }
  else assert(0);

  VectorP fitResult_alpha = fitResult.get_alpha();
  MatrixPxP fitResult_V_alpha = fitResult.get_V_alpha();

  double fitResult_nuTauPlusPx = fitResult_alpha(3);
  double fitResult_nuTauPlusPy = fitResult_alpha(4);
  double fitResult_nuTauPlus_dPzdPx, fitResult_nuTauPlus_dPzdPy;
  bool fitResult_nuTauPlus_errorFlag = false;
  double fitResult_nuTauPlusPz = comp_nuPz(
           startPos_visTauPlusP4, fitResult_nuTauPlusPx, fitResult_nuTauPlusPy, signTauPlus,
           fitResult_nuTauPlus_dPzdPx, fitResult_nuTauPlus_dPzdPy,
	   fitResult_nuTauPlus_errorFlag,
           verbosity_);
  double fitResult_nuTauMinusPx = fitResult_alpha(8);
  double fitResult_nuTauMinusPy = fitResult_alpha(9);
  double fitResult_nuTauMinus_dPzdPx, fitResult_nuTauMinus_dPzdPy;
  bool fitResult_nuTauMinus_errorFlag = false;
  double fitResult_nuTauMinusPz = comp_nuPz(
           startPos_visTauMinusP4, fitResult_nuTauMinusPx, fitResult_nuTauMinusPy, signTauMinus,
           fitResult_nuTauMinus_dPzdPx, fitResult_nuTauMinus_dPzdPy,
	   fitResult_nuTauMinus_errorFlag,
           verbosity_);

  if ( !(fitResult_nuTauPlus_errorFlag || fitResult_nuTauMinus_errorFlag) )
  {
    // CV: store results of kinematic fit in KinematicEvent class;
    //     for the syntax of retrieving a small covariance matrix from a larger one, 
    //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
    kineEvt_kinFit.pv_ = reco::Candidate::Point(fitResult_alpha(0), fitResult_alpha(1), fitResult_alpha(2));
    kineEvt_kinFit.pvCov_ = fitResult_V_alpha.Sub<Matrix3x3>(0,0);
    kineEvt_kinFit.recoilP4_ = reco::Candidate::LorentzVector(fitResult_alpha(13), fitResult_alpha(14), fitResult_alpha(15), fitResult_alpha(16));
    kineEvt_kinFit.recoilCov_ = fitResult_V_alpha.Sub<Matrix4x4>(13,13);
    kineEvt_kinFit.nuTauPlusP4_ = build_nuP4(fitResult_nuTauPlusPx, fitResult_nuTauPlusPy, fitResult_nuTauPlusPz);
    kineEvt_kinFit.nuTauPlusP4_isValid_ = true;
    kineEvt_kinFit.nuTauPlusCov_ = build_nuCov(fitResult_V_alpha.Sub<Matrix2x2>(3,3), fitResult_nuTauPlus_dPzdPx, fitResult_nuTauPlus_dPzdPy);
    kineEvt_kinFit.tauPlusP4_ = startPos_visTauPlusP4 + kineEvt_kinFit.nuTauPlusP4_;
    kineEvt_kinFit.tauPlusP4_isValid_ = true;
    reco::Candidate::Point fitResult_svTauPlus = reco::Candidate::Point(fitResult_alpha(5), fitResult_alpha(6), fitResult_alpha(7));
    kineEvt_kinFit.svTauPlus_ = rotateVector(fitResult_svTauPlus, kineEvt.tauPlus_rotMatrix_rnk2xyz()) 
      + reco::Candidate::Vector(startPos_pv.x(), startPos_pv.y(), startPos_pv.z());
    kineEvt_kinFit.svTauPlusCov_ = rotateCovMatrix(fitResult_V_alpha.Sub<Matrix3x3>(5,5), kineEvt.tauPlus_rotMatrix_rnk2xyz());        
    kineEvt_kinFit.nuTauMinusP4_ = build_nuP4(fitResult_nuTauMinusPx, fitResult_nuTauMinusPy, fitResult_nuTauMinusPz);
    kineEvt_kinFit.nuTauMinusP4_isValid_ = true;       
    kineEvt_kinFit.nuTauMinusCov_ = build_nuCov(fitResult_V_alpha.Sub<Matrix2x2>(8,8), fitResult_nuTauMinus_dPzdPx, fitResult_nuTauMinus_dPzdPy);
    kineEvt_kinFit.tauMinusP4_ = startPos_visTauMinusP4 + kineEvt_kinFit.nuTauMinusP4_;
    kineEvt_kinFit.tauMinusP4_isValid_ = true;
    reco::Candidate::Point fitResult_svTauMinus = reco::Candidate::Point(fitResult_alpha(10), fitResult_alpha(11), fitResult_alpha(12));
    kineEvt_kinFit.svTauMinus_ = rotateVector(fitResult_svTauMinus, kineEvt.tauMinus_rotMatrix_rnk2xyz())
      + reco::Candidate::Vector(startPos_pv.x(), startPos_pv.y(), startPos_pv.z());
    kineEvt_kinFit.svTauMinusCov_ = rotateCovMatrix(fitResult_V_alpha.Sub<Matrix3x3>(10,10), kineEvt.tauMinus_rotMatrix_rnk2xyz());
    kineEvt_kinFit.kinFitCov_ = fitResult_V_alpha;
    kineEvt_kinFit.kinFitStatus_ = fitResult.get_status();
    kineEvt_kinFit.kinFitChi2_ = fitResult.get_chi2();
    kineEvt_kinFit.kinFit_isValid_ = true;
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

