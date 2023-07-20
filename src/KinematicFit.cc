#include "TauAnalysis/Entanglement/interface/KinematicFit.h"

#include "DataFormats/Math/interface/deltaR.h"                    // deltaR2()
#include "DataFormats/Math/interface/Matrix.h"                    // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                    // math::Vector
#include "DataFormats/TauReco/interface/PFTau.h"                  // reco::PFTau::kOneProng0PiZero

#include "TauAnalysis/Entanglement/interface/cmsException.h"      // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"         // mHiggs, mTau
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix*, math::Vector*, printCovMatrix()
#include "TauAnalysis/Entanglement/interface/square.h"            // square()

#include "Math/Functions.h"                                       // ROOT::Math::Dot(), ROOT::Math::Similarity(), ROOT::Math::Transpose() 

#include <cmath>                                                  // std::abs(), std::fabs(), std::sqrt()
#include <iostream>                                               // std::cout
#include <string>                                                 // std::string

using namespace kinFit;

KinematicFit::KinematicFit(const edm::ParameterSet& cfg)
  : spinAnalyzer_(cfg)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{}

KinematicFit::~KinematicFit()
{}

namespace
{
  math::Matrix4x4
  get_tauCov()
  {
    // CV: the four-vectors of tau+ and tau- are not really "measured";
    //     we acccount for this by setting their covariance matrix to diagonal matrix with large values on the diagonal,
    //     following the "huge error method" described in Section 6 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
    math::Matrix4x4 tauCov;
    tauCov(0,0) = square(1.e+3);
    tauCov(1,1) = square(1.e+3);
    tauCov(2,2) = square(1.e+3);
    tauCov(3,3) = square(1.e+3);
    return tauCov;
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

  bool hasConverged = false;
  bool hasFailed = false;
  int iteration = 1;
  const int max_iterations = 10;
  double min_chi2 = -1.;
  while ( !hasConverged && !hasFailed && iteration <= max_iterations )
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "iteration #" << iteration << ":\n";
    }

    const reco::Candidate::Point& pv = kineEvt_kinfit.pv();
    double pvX = pv.X();
    double pvY = pv.Y();
    double pvZ = pv.Z();
    const math::Matrix3x3& pvCov = kineEvt_kinfit.pvCov();

    const reco::Candidate::LorentzVector& recoilP4 = kineEvt_kinfit.recoilP4();
    double recoilPx = recoilP4.px();
    double recoilPy = recoilP4.py();
    double recoilPz = recoilP4.pz();
    double recoilE  = recoilP4.energy();
    const math::Matrix4x4& recoilCov = kineEvt_kinfit.recoilCov();

    reco::Candidate::LorentzVector tauPlusP4 = kineEvt_kinfit.tauPlusP4();
    double tauPlusPx = tauPlusP4.px();
    double tauPlusPy = tauPlusP4.py();
    double tauPlusPz = tauPlusP4.pz();
    double tauPlusPt = tauPlusP4.pt();
    double tauPlusP  = tauPlusP4.P();
    double tauPlusE  = tauPlusP4.energy();
    math::Matrix4x4 tauPlusCov = get_tauCov();
    assert(kineEvt_kinfit.svTauPlus_isValid());
    const reco::Candidate::Point& svTauPlus = kineEvt_kinfit.svTauPlus();
    double svTauPlusX = svTauPlus.X();
    double svTauPlusY = svTauPlus.Y();
    double svTauPlusZ = svTauPlus.Z();
    const math::Matrix3x3& svTauPlusCov = kineEvt_kinfit.svTauPlusCov();
    auto tauPlusD3 = svTauPlus - pv;
    double tauPlusDx = tauPlusD3.X();
    double tauPlusDy = tauPlusD3.Y();
    double tauPlusDz = tauPlusD3.Z();
    double tauPlusDt = std::sqrt(tauPlusD3.Perp2());
    double tauPlusD  = std::sqrt(tauPlusD3.Mag2());

    reco::Candidate::LorentzVector visTauPlusP4 = kineEvt_kinfit.visTauPlusP4();
    double visTauPlusPx = visTauPlusP4.px();
    double visTauPlusPy = visTauPlusP4.py();
    double visTauPlusPz = visTauPlusP4.pz();
    double visTauPlusE  = visTauPlusP4.energy();

    reco::Candidate::LorentzVector tauMinusP4 = kineEvt_kinfit.tauMinusP4();
    double tauMinusPx = tauMinusP4.px();
    double tauMinusPy = tauMinusP4.py();
    double tauMinusPz = tauMinusP4.pz();
    double tauMinusPt = tauMinusP4.pt();
    double tauMinusP  = tauMinusP4.P();
    double tauMinusE  = tauMinusP4.energy();
    math::Matrix4x4 tauMinusCov = get_tauCov();
    assert(kineEvt_kinfit.svTauMinus_isValid());
    const reco::Candidate::Point& svTauMinus = kineEvt_kinfit.svTauMinus();
    double svTauMinusX = svTauMinus.X();
    double svTauMinusY = svTauMinus.Y();
    double svTauMinusZ = svTauMinus.Z();
    const math::Matrix3x3& svTauMinusCov = kineEvt_kinfit.svTauMinusCov();
    auto tauMinusD3 = svTauMinus - pv;
    double tauMinusDx = tauMinusD3.X();
    double tauMinusDy = tauMinusD3.Y();
    double tauMinusDz = tauMinusD3.Z();
    double tauMinusDt = std::sqrt(tauMinusD3.Perp2());
    double tauMinusD  = std::sqrt(tauMinusD3.Mag2());

    reco::Candidate::LorentzVector visTauMinusP4 = kineEvt_kinfit.visTauMinusP4();
    double visTauMinusPx = visTauMinusP4.px();
    double visTauMinusPy = visTauMinusP4.py();
    double visTauMinusPz = visTauMinusP4.pz();
    double visTauMinusE  = visTauMinusP4.energy();

    //-----------------------------------------------------------------------------------------------
    // CV: the nomenclature of the different matrices and vectors used for the kinematic fit
    //     follows the one introduced in Section 3 of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
    //     For the first iteration, the vector alpha0 refers to the starting-point of the kinematic fit,
    //     computed by the function KinematicFitStartPosFinder::operator()
    //     For subsequent iterations, the vector alpha0 refers to the result alpha obtained by the previous iteration
    //     The symbol H in the comments refers to the constraint equations,
    //     cf. Eq. (3) in https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
    //-----------------------------------------------------------------------------------------------

    // CV: build covariance matrix of all measured parameters;
    //     for the syntax of "embedding" a small covariance matrix into a larger one, 
    //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
    math::MatrixPxP V_alpha0;
    V_alpha0.Place_at(pvCov        ,  0,  0);
    V_alpha0.Place_at(tauPlusCov   ,  3,  3);
    V_alpha0.Place_at(svTauPlusCov ,  7,  7);
    V_alpha0.Place_at(tauMinusCov  , 10, 10);
    V_alpha0.Place_at(svTauMinusCov, 14, 14);
    V_alpha0.Place_at(recoilCov    , 17, 17);

    math::MatrixCxP D;
    math::VectorC   d;
    // CV: add Higgs mass constraint
    //       H = (tauPlusE + tauMinusE)^2 - (tauPlusPx + tauMinusPx)^2 - (tauPlusPy + tauMinusPy)^2 - (tauPlusPz + tauMinusPz)^2 - mHiggs^2 = 0
    D( 0, 3) = -2.*(tauPlusPx + tauMinusPx);                                // dH/dtauPlusPx
    D( 0, 4) = -2.*(tauPlusPy + tauMinusPy);                                // dH/dtauPlusPy
    D( 0, 5) = -2.*(tauPlusPz + tauMinusPz);                                // dH/dtauPlusPz
    D( 0, 6) = +2.*(tauPlusE  + tauMinusE );                                // dH/dtauPlusE
    D( 0,10) = -2.*(tauPlusPx + tauMinusPx);                                // dH/dtauMinusPx
    D( 0,11) = -2.*(tauPlusPy + tauMinusPy);                                // dH/dtauMinusPy
    D( 0,12) = -2.*(tauPlusPz + tauMinusPz);                                // dH/dtauMinusPz
    D( 0,13) = +2.*(tauPlusE  + tauMinusE );                                // dH/dtauMinusE
    d( 0) = square((tauPlusP4 + tauMinusP4).mass()) - square(mHiggs);
    // CV: add tau+ mass constraint
    //       H = tauPlusE^2 - tauPlusPx^2 - tauPlusPy^2 - tauPlusPz^2 - mTau^2 = 0
    D( 1, 3) = -2.*tauPlusPx;                                               // dH/dtauPlusPx
    D( 1, 4) = -2.*tauPlusPy;                                               // dH/dtauPlusPy
    D( 1, 5) = -2.*tauPlusPz;                                               // dH/dtauPlusPz
    D( 1, 6) = +2.*tauPlusE;                                                // dH/dtauPlusE
    d( 1) = square(tauPlusP4.mass()) - square(mTau);
    // CV: add neutrino mass constraint for tau+
    //       H = (tauPlusE - visTauPlusE)^2 - (tauPlusPx - visTauPlusPx)^2 - (tauPlusPy - visTauPlusPy)^2 - (tauPlusPz - visTauPlusPz)^2 = 0
    D( 2, 3) = -2.*(tauPlusPx - visTauPlusPx );                             // dH/dtauPlusPx
    D( 2, 4) = -2.*(tauPlusPy - visTauPlusPy );                             // dH/dtauPlusPy
    D( 2, 5) = -2.*(tauPlusPz - visTauPlusPz );                             // dH/dtauPlusPz
    D( 2, 6) = +2.*(tauPlusE  - visTauPlusE  );                             // dH/dtauPlusE
    d( 2) = square((tauPlusP4 - visTauPlusP4).mass());
    // CV: add "parallelism" constraint for tau+
    //     The constraint equations H_Pphi and H_Ptheta are given by Eqs. (4.10) and (4.11)
    //     in Section 4.1.3.3 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
    //     The small correction for the helix bend between tau lepton production and decay in Eq. (4.13) is not used, 
    //     as this correction is expected to be on the level of 50 microrad, which we consider negligible
    D( 3, 3) =  (1. - tauPlusPx/tauPlusPt)/tauPlusPy;                       // dH_Pphi/dtauPlusPx
    D( 3, 4) = -(tauPlusPx/tauPlusPy)*D(3,3);                               // dH_Pphi/dtauPlusPy
    D( 3, 0) =  (1. - tauPlusDx/tauPlusDt)/tauPlusDy;                       // dH_Pphi/dpvX
    D( 3, 1) = -(tauPlusDx/tauPlusDy)*D(3,0);                               // dH_Pphi/dpvY
    D( 3, 7) = -D(3,0);                                                     // dH_Pphi/dsvTauPlusX
    D( 3, 8) = -D(3,1);                                                     // dH_Pphi/dsvTauPlusY
    d( 3) =  (tauPlusDt - tauPlusDx)/tauPlusDy + (tauPlusPx - tauPlusPt)/tauPlusPy;
    D( 4, 3) =  (tauPlusPx/tauPlusPz)*(1./tauPlusPt - 1./tauPlusP);         // dH_Ptheta/dtauPlusPx
    D( 4, 4) =  (tauPlusPy/tauPlusPz)*(1./tauPlusPt - 1./tauPlusP);         // dH_Ptheta/dtauPlusPy
    D( 4, 5) = -1./tauPlusP + (tauPlusP - tauPlusPt)/square(tauPlusPz);     // dH_Ptheta/dtauPlusPz
    D( 4, 0) =  (tauPlusDx/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);         // dH_Ptheta/dpvX
    D( 4, 1) =  (tauPlusDy/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);         // dH_Ptheta/dpvY
    D( 4, 2) = -1./tauPlusD + (tauPlusD - tauPlusDt)/square(tauPlusDz);     // dH_Ptheta/dpvZ
    D( 4, 7) = -D(4,0);                                                     // dH_Ptheta/dsvTauPlusX
    D( 4, 8) = -D(4,1);                                                     // dH_Ptheta/dsvTauPlusY
    D( 4, 9) = -D(4,2);                                                     // dH_Ptheta/dsvTauPlusZ
    d( 4) =  (tauPlusD - tauPlusDt)/tauPlusDz + (tauPlusPt - tauPlusP)/tauPlusPz;
    // CV: add tau- mass constraint
    //       H = tauMinusE^2 - tauMinusPx^2 - tauMinusPy^2 - tauMinusPz^2 - mTau^2 = 0
    D( 5,10) = -2.*tauMinusPx;                                              // dH/dtauMinusPx
    D( 5,11) = -2.*tauMinusPy;                                              // dH/dtauMinusPy
    D( 5,12) = -2.*tauMinusPz;                                              // dH/dtauMinusPz
    D( 5,13) = +2.*tauMinusE;                                               // dH/dtauMinusE
    d( 5) = square(tauMinusP4.mass()) - square(mTau);
    // CV: add neutrino mass constraint for tau-
    //       H = (tauMinusE - visTauMinusE)^2 - (tauMinusPx - visTauMinusPx)^2 - (tauMinusPy - visTauMinusPy)^2 - (tauMinusPz - visTauMinusPz)^2 = 0
    D( 6,10) = -2.*(tauMinusPx - visTauMinusPx);                            // dH/dtauMinusPx
    D( 6,11) = -2.*(tauMinusPy - visTauMinusPy);                            // dH/dtauMinusPy
    D( 6,12) = -2.*(tauMinusPz - visTauMinusPz);                            // dH/dtauMinusPz
    D( 6,13) = +2.*(tauMinusE  - visTauMinusE );                            // dH/dtauMinusE
    d( 6) = square((tauMinusP4 - visTauMinusP4).mass());
    // CV: add "parallelism" constraint for tau+
    //     The constraint equations H_Pphi and H_Ptheta are given by Eqs. (4.10) and (4.11)
    //     in Section 4.1.3.3 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
    //     The small correction for the helix bend between tau lepton production and decay in Eq. (4.13) is not used, 
    //     as this correction is expected to be on the level of 50 microrad, which we consider negligible
    D( 7,10) =  (1. - tauMinusPx/tauMinusPt)/tauMinusPy;                    // dH_Pphi/dtauMinusPx
    D( 7,11) = -(tauMinusPx/tauMinusPy)*D(7,10);                            // dH_Pphi/dtauMinusPy
    D( 7, 0) =  (1. - tauMinusDx/tauMinusDt)/tauMinusDy;                    // dH_Pphi/dpvX
    D( 7, 1) = -(tauMinusDx/tauMinusDy)*D(7,0);                             // dH_Pphi/dpvY
    D( 7,14) = -D(7,0);                                                     // dH_Pphi/dsvTauMinusX
    D( 7,15) = -D(7,1);                                                     // dH_Pphi/dsvTauMinusY
    d( 7) =  (tauMinusDt - tauMinusDx)/tauMinusDy + (tauMinusPx - tauMinusPt)/tauMinusPy;
    D( 8,10) =  (tauMinusPx/tauMinusPz)*(1./tauMinusPt - 1./tauMinusP);     // dH_Ptheta/dtauMinusPx
    D( 8,11) =  (tauMinusPy/tauMinusPz)*(1./tauMinusPt - 1./tauMinusP);     // dH_Ptheta/dtauMinusPy
    D( 8,12) = -1./tauMinusP + (tauMinusP - tauMinusPt)/square(tauMinusPz); // dH_Ptheta/dtauMinusPz
    D( 8, 0) =  (tauMinusDx/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);     // dH_Ptheta/dpvX
    D( 8, 1) =  (tauMinusDy/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);     // dH_Ptheta/dpvY
    D( 8, 2) = -1./tauMinusD + (tauMinusD - tauMinusDt)/square(tauMinusDz); // dH_Ptheta/dpvZ
    D( 8,14) = -D(8,0);                                                     // dH_Ptheta/dsvTauMinusX
    D( 8,15) = -D(8,1);                                                     // dH_Ptheta/dsvTauMinusY
    D( 8,16) = -D(8,2);                                                     // dH_Ptheta/dsvTauMinusZ
    d(8) =  (tauMinusD - tauMinusDt)/tauMinusDz + (tauMinusPt - tauMinusP)/tauMinusPz;
    // CV: add constraint that recoil = tau+ + tau-
    //       H_px     = tauPlusPx + tauMinusPx - recoilPx = 0
    //       H_py     = tauPlusPy + tauMinusPy - recoilPy = 0
    //       H_pz     = tauPlusPz + tauMinusPz - recoilPz = 0
    //       H_energy = tauPlusE  + tauMinusE  - recoilE  = 0
    D( 9, 3) = +1.;                                                       // dH_px/dtauPlusPx
    D( 9,10) = +1.;                                                       // dH_px/dtauMinusPx
    D( 9,17) = -1.;                                                       // dH_px/drecoilPx
    d( 9) = tauPlusPx + tauMinusPx - recoilPx;
    D(10, 4) = +1.;                                                       // dH_py/dtauPlusPy
    D(10,11) = +1.;                                                       // dH_py/dtauMinusPy
    D(10,18) = -1.;                                                       // dH_py/drecoilPy
    d(10) = tauPlusPy + tauMinusPy - recoilPy;
    D(11, 5) = +1.;                                                       // dH_pz/dtauPlusPz
    D(11,12) = +1.;                                                       // dH_pz/dtauMinusPz
    D(11,19) = -1.;                                                       // dH_pz/drecoilPz
    d(11) = tauPlusPz + tauMinusPz - recoilPz;
    D(12, 6) = +1.;                                                       // dH_energy/dtauPlusE
    D(12,13) = +1.;                                                       // dH_energy/dtauMinusE
    D(12,20) = -1.;                                                       // dH_energy/drecoilE
    d(12) = tauPlusE  + tauMinusE  - recoilE;
    if ( verbosity_ >= 1 )
    {
      printCovMatrix("V_alpha0", V_alpha0);
      std::cout << "D:\n";
      std::cout << D << "\n";
      std::cout << "d:\n";
      std::cout << d << "\n";
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

    math::MatrixCxC Vinv_D = D*V_alpha0*DT;
    int V_D_errorFlag = 0;
    math::MatrixCxC V_D = Vinv_D.Inverse(V_D_errorFlag);
    if ( V_D_errorFlag != 0 )
    {
      printCovMatrix("Vinv_D", Vinv_D);
      throw cmsException("KinematicFit::operator()", __LINE__)
        << "Failed to invert matrix Vinv_D !!\n";
    }

    math::VectorC lambda = V_D*d;
    if ( verbosity_ >= 1 )
    {
      std::cout << "lambda:\n";
      std::cout << lambda << "\n";
    }

    math::VectorP dalpha = -V_alpha0*DT*lambda;
    if ( verbosity_ >= 1 )
    {
      std::cout << "dalpha:\n";
      std::cout << dalpha << "\n";
    }
 
    math::MatrixPxP V_alpha = V_alpha0 - V_alpha0*DT*V_D*D*V_alpha0;
    if ( verbosity_ >= 1 )
    {
      printCovMatrix("V_alpha", V_alpha);
    }

    // CV: compute chi^2
    double chi2 = ROOT::Math::Dot(lambda, d);
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
      std::cout << "sum = " << residuals_sum << "\n";
    }

    if ( verbosity_ >= 1 )
    {
      math::VectorP pulls;
      for ( int idx = 0; idx < numParameters; ++idx )
      {      
        pulls(idx) = dalpha(idx)/std::sqrt(V_alpha0(idx,idx));
      }
      std::cout << "pulls:\n";
      std::cout << pulls << "\n";
    }

    // CV: store results of kinematic fit in KinematicEvent class;
    //     for the syntax of retrieving a small covariance matrix from a larger one, 
    //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
    const double epsilon = 1.e-2;
    if ( residuals_sum < epsilon && (iteration == 1 || chi2 < min_chi2) )
    {
      reco::Candidate::Point pv_kinfit(pvX + dalpha(0), pvY + dalpha(1), pvZ + dalpha(2));
      kineEvt_kinfit.pv_ = pv;
      kineEvt_kinfit.pvCov_ = V_alpha.Sub<math::Matrix3x3>(0,0);
      reco::Candidate::LorentzVector recoilP4_kinfit(recoilPx + dalpha(17), recoilPy + dalpha(18), recoilPz + dalpha(19), recoilE + dalpha(20));
      kineEvt_kinfit.recoilP4_ = recoilP4_kinfit;
      kineEvt_kinfit.recoilCov_ = V_alpha.Sub<math::Matrix4x4>(17,17);
      reco::Candidate::LorentzVector tauPlusP4_kinfit(tauPlusPx + dalpha(3), tauPlusPy + dalpha(4), tauPlusPz + dalpha(5), tauPlusE + dalpha(6));
      kineEvt_kinfit.tauPlusP4_ = tauPlusP4_kinfit;
      kineEvt_kinfit.tauPlusP4_isValid_ = true;
      kineEvt_kinfit.tauPlusCov_ = V_alpha.Sub<math::Matrix4x4>(3,3);
      reco::Candidate::Point svTauPlus_kinfit(svTauPlusX + dalpha(7), svTauPlusY + dalpha(8), svTauPlusZ + dalpha(9));
      kineEvt_kinfit.svTauPlus_ = svTauPlus_kinfit;
      kineEvt_kinfit.svTauPlusCov_ = V_alpha.Sub<math::Matrix3x3>(7,7);
      reco::Candidate::LorentzVector tauMinusP4_kinfit(tauMinusPx + dalpha(10), tauMinusPy + dalpha(11), tauMinusPz + dalpha(12), tauMinusE + dalpha(13));
      kineEvt_kinfit.tauMinusP4_ = tauMinusP4_kinfit;
      kineEvt_kinfit.tauMinusP4_isValid_ = true;
      kineEvt_kinfit.tauMinusCov_ = V_alpha.Sub<math::Matrix4x4>(10,10);
      reco::Candidate::Point svTauMinus_kinfit(svTauMinusX + dalpha(14), svTauMinusY + dalpha(15), svTauMinusZ + dalpha(16));
      kineEvt_kinfit.svTauMinus_ = svTauMinus_kinfit;
      kineEvt_kinfit.svTauMinusCov_ = V_alpha.Sub<math::Matrix3x3>(14,14);
      kineEvt_kinfit.kinFitCov_ = V_alpha;
      kineEvt_kinfit.kinFitChi2_ = chi2;
      kineEvt_kinfit.kinFit_isValid_ = true;
      if ( verbosity_ >= 1 )
      {
        std::cout << "constraint equations:\n";
        reco::Candidate::LorentzVector higgsP4_kinfit = tauPlusP4_kinfit + tauMinusP4_kinfit;
        std::cout << "Higgs mass = " << higgsP4_kinfit.mass() << "\n";
        std::cout << "tau+ mass = " << tauPlusP4_kinfit.mass() << "\n";
        double nuTauPlusPx   = tauPlusP4_kinfit.px()     - visTauPlusP4.px();
        double nuTauPlusPy   = tauPlusP4_kinfit.py()     - visTauPlusP4.py();
        double nuTauPlusPz   = tauPlusP4_kinfit.pz()     - visTauPlusP4.pz();
        double nuTauPlusE    = tauPlusP4_kinfit.energy() - visTauPlusP4.energy();
        std::cout << "neutrino from tau+ decay:" 
                  << " E = " << nuTauPlusE << ", Px = " << nuTauPlusPx << ", Py = " << nuTauPlusPy << ", Pz = " << nuTauPlusPz << ","
                  << " mass = " << (tauPlusP4 - visTauPlusP4).mass() << "\n";
        auto tauPlusD3_kinfit = svTauPlus_kinfit - pv_kinfit;
        std::cout << "phi of tau+: four-vector = " << tauPlusP4_kinfit.phi() << ", decay vertex = " << tauPlusD3_kinfit.phi() << "\n";
        std::cout << "theta of tau+: four-vector = " << tauPlusP4_kinfit.theta() << ", decay vertex = " << tauPlusD3_kinfit.theta() << "\n";
        std::cout << "tau- mass = " << tauMinusP4_kinfit.mass() << "\n";
        double nuTauMinusPx   = tauMinusP4_kinfit.px()     - visTauMinusP4.px();
        double nuTauMinusPy   = tauMinusP4_kinfit.py()     - visTauMinusP4.py();
        double nuTauMinusPz   = tauMinusP4_kinfit.pz()     - visTauMinusP4.pz();
        double nuTauMinusE    = tauMinusP4_kinfit.energy() - visTauMinusP4.energy();
        std::cout << "neutrino from tau- decay:" 
                  << " E = " << nuTauMinusE << ", Px = " << nuTauMinusPx << ", Py = " << nuTauMinusPy << ", Pz = " << nuTauMinusPz << ","
                  << " mass = " << (tauMinusP4 - visTauMinusP4).mass() << "\n";
        auto tauMinusD3_kinfit = svTauMinus_kinfit - pv_kinfit;
        std::cout << "phi of tau-: four-vector = " << tauMinusP4_kinfit.phi() << ", decay vertex = " << tauMinusD3_kinfit.phi() << "\n";
        std::cout << "theta of tau-: four-vector = " << tauMinusP4_kinfit.theta() << ", decay vertex = " << tauMinusD3_kinfit.theta() << "\n";
        std::cout << "Higgs - recoil:"
                  " dE = "  << higgsP4_kinfit.energy() - recoilP4_kinfit.energy() << ","
                  " dPx = " << higgsP4_kinfit.px()     - recoilP4_kinfit.px()     << ","
                  " dPy = " << higgsP4_kinfit.py()     - recoilP4_kinfit.py()     << ","
                  " dPz = " << higgsP4_kinfit.pz()     - recoilP4_kinfit.pz()     << "\n";
      }
    }

    double diff = std::sqrt(ROOT::Math::Dot(dalpha, dalpha));
    if ( residuals_sum > epsilon )
    {
      hasFailed = true;
    }
    else if ( diff < 1.e-2 )
    {
      hasConverged = true;
    }
  }
  if ( hasFailed || iteration > max_iterations )
  {
    std::cerr << "WARNING: KinematicFit failed to converge !!" << std::endl;
  }
  if ( verbosity_ >= 1 )
  {
    std::cout << "#iterations = " << iteration << " (max_iterations = " << max_iterations << ")\n";
  }

  reco::Candidate::Vector hPlus = spinAnalyzer_(kineEvt_kinfit, SpinAnalyzerBase::kTauPlus);
  kineEvt_kinfit.hPlus_ = hPlus;
  kineEvt_kinfit.hPlus_isValid_ = true;
  reco::Candidate::Vector hMinus = spinAnalyzer_(kineEvt_kinfit, SpinAnalyzerBase::kTauMinus);
  kineEvt_kinfit.hMinus_ = hMinus;
  kineEvt_kinfit.hMinus_isValid_ = true;

  return kineEvt_kinfit;
}

