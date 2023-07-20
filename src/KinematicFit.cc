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
    //     follows the one introduced in Section 6.1 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
    //     For the first iteration, the vector eta0 refers to the starting-point of the kinematic fit,
    //     computed by the function KinematicFitStartPosFinder::operator()
    //     For subsequent iterations, the vector eta0 refers to the result eta obtained by the previous iteration
    //     The symbol H in the comments refers to the constraint equation,
    //     cf. Section 5 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
    //-----------------------------------------------------------------------------------------------

    // CV: build covariance matrix of all measured parameters;
    //     for the syntax of "embedding" a small covariance matrix into a larger one, 
    //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
    math::MatrixMxM V_eta0;
    V_eta0.Place_at(pvCov        , 0, 0);
    V_eta0.Place_at(svTauPlusCov , 3, 3);
    V_eta0.Place_at(svTauMinusCov, 6, 6);
    V_eta0.Place_at(recoilCov    , 9, 9);

    math::MatrixCxM D;
    math::MatrixCxP E;
    math::VectorC   d;
    // CV: add Higgs mass constraint
    E(0,0) = -2.*(tauPlusPx + tauMinusPx);                                // dH/dtauPlusPx
    E(0,1) = -2.*(tauPlusPy + tauMinusPy);                                // dH/dtauPlusPy
    E(0,2) = -2.*(tauPlusPz + tauMinusPz);                                // dH/dtauPlusPz
    E(0,3) = +2.*(tauPlusE  + tauMinusE );                                // dH/dtauPlusE
    E(0,4) = -2.*(tauPlusPx + tauMinusPx);                                // dH/dtauMinusPx
    E(0,5) = -2.*(tauPlusPy + tauMinusPy);                                // dH/dtauMinusPy
    E(0,6) = -2.*(tauPlusPz + tauMinusPz);                                // dH/dtauMinusPz
    E(0,7) = +2.*(tauPlusE  + tauMinusE );                                // dH/dtauMinusE
    d(0) = (tauPlusP4 + tauMinusP4).mass() - mHiggs;
    // CV: add tau+ mass constraint
    E(1,0) = -2.*tauPlusPx;                                               // dH/dtauPlusPx
    E(1,1) = -2.*tauPlusPy;                                               // dH/dtauPlusPy
    E(1,2) = -2.*tauPlusPz;                                               // dH/dtauPlusPz
    E(1,3) = +2.*tauPlusE;                                                // dH/dtauPlusE
    d(1) = tauPlusP4.mass() - mTau;
    // CV: add tau- mass constraint
    E(2,4) = -2.*tauMinusPx;                                              // dH/dtauMinusPx
    E(2,5) = -2.*tauMinusPy;                                              // dH/dtauMinusPy
    E(2,6) = -2.*tauMinusPz;                                              // dH/dtauMinusPz
    E(2,7) = +2.*tauMinusE;                                               // dH/dtauMinusE
    d(2) = tauMinusP4.mass() - mTau;
    // CV: add neutrino mass constraint for tau+
    E(3,0) = -2.*(tauPlusPx  - visTauPlusPx );                            // dH/dtauPlusPx
    E(3,1) = -2.*(tauPlusPy  - visTauPlusPy );                            // dH/dtauPlusPy
    E(3,2) = -2.*(tauPlusPz  - visTauPlusPz );                            // dH/dtauPlusPz
    E(3,3) = +2.*(tauPlusE   - visTauPlusE  );                            // dH/dtauPlusE
    d(3) = (tauPlusP4 - visTauPlusP4).mass();
    // CV: add neutrino mass constraint for tau-
    E(4,4) = -2.*(tauMinusPx - visTauMinusPx);                            // dH/dtauMinusPx
    E(4,5) = -2.*(tauMinusPy - visTauMinusPy);                            // dH/dtauMinusPy
    E(4,6) = -2.*(tauMinusPz - visTauMinusPz);                            // dH/dtauMinusPz
    E(4,7) = +2.*(tauMinusE  - visTauMinusE );                            // dH/dtauMinusE
    d(4) = (tauMinusP4 - visTauMinusP4).mass();
    // CV: add "parallelism" constraint for tau+
    //     cf. Section 4.1.3.3 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
    E(5,0) =  (1. - tauPlusPx/tauPlusPt)/tauPlusPy;                       // dH_phi/dtauPlusPx
    E(5,1) = -(tauPlusPx/tauPlusPy)*E(5,0);                               // dH_phi/dtauPlusPy
    D(5,0) =  (1. - tauPlusDx/tauPlusDt)/tauPlusDy;                       // dH_phi/dpvX
    D(5,1) = -(tauPlusDx/tauPlusDy)*D(5,0);                               // dH_phi/dpvY
    D(5,3) = -D(5,0);                                                     // dH_phi/dsvTauPlusX
    D(5,4) = -D(5,1);                                                     // dH_phi/dsvTauPlusY
    d(5) =  (tauPlusDt - tauPlusDx)/tauPlusDy + (tauPlusPx - tauPlusPt)/tauPlusPy;
    E(6,0) =  (tauPlusPx/tauPlusPz)*(1./tauPlusPt - 1./tauPlusP);         // dH_theta/dtauPlusPx
    E(6,1) =  (tauPlusPy/tauPlusPz)*(1./tauPlusPt - 1./tauPlusP);         // dH_theta/dtauPlusPy
    E(6,2) = -1./tauPlusP + (tauPlusP - tauPlusPt)/square(tauPlusPz);     // dH_theta/dtauPlusPz
    D(6,0) =  (tauPlusDx/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);         // dH_theta/dpvX
    D(6,1) =  (tauPlusDy/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);         // dH_theta/dpvY
    D(6,2) = -1./tauPlusD + (tauPlusD - tauPlusDt)/square(tauPlusDz);     // dH_theta/dpvZ
    D(6,3) = -D(6,0);                                                     // dH_theta/dsvTauPlusX
    D(6,4) = -D(6,1);                                                     // dH_theta/dsvTauPlusY
    D(6,5) = -D(6,2);                                                     // dH_theta/dsvTauPlusZ
    d(6) =  (tauPlusD - tauPlusDt)/tauPlusDz + (tauPlusPt - tauPlusP)/tauPlusPz;
    // CV: add "parallelism" constraint for tau-
    //     cf. Section 4.1.3.3 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
    E(7,4) =  (1. - tauMinusPx/tauMinusPt)/tauMinusPy;                    // dH_phi/dtauMinusPx
    E(7,5) = -(tauMinusPx/tauMinusPy)*E(7,4);                             // dH_phi/dtauMinusPy
    D(7,0) =  (1. - tauMinusDx/tauMinusDt)/tauMinusDy;                    // dH_phi/dpvX
    D(7,1) = -(tauMinusDx/tauMinusDy)*D(7,0);                             // dH_phi/dpvY
    D(7,6) = -D(7,0);                                                     // dH_phi/dsvTauMinusX
    D(7,7) = -D(7,1);                                                     // dH_phi/dsvTauMinusY
    d(7) =  (tauMinusDt - tauMinusDx)/tauMinusDy + (tauMinusPx - tauMinusPt)/tauMinusPy;
    E(8,4) =  (tauMinusPx/tauMinusPz)*(1./tauMinusPt - 1./tauMinusP);     // dH_theta/dtauMinusPx
    E(8,5) =  (tauMinusPy/tauMinusPz)*(1./tauMinusPt - 1./tauMinusP);     // dH_theta/dtauMinusPy
    E(8,6) = -1./tauMinusP + (tauMinusP - tauMinusPt)/square(tauMinusPz); // dH_theta/dtauMinusPz
    D(8,0) =  (tauMinusDx/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);     // dH_theta/dpvX
    D(8,1) =  (tauMinusDy/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);     // dH_theta/dpvY
    D(8,2) = -1./tauMinusD + (tauMinusD - tauMinusDt)/square(tauMinusDz); // dH_theta/dpvZ
    D(8,6) = -D(8,0);                                                     // dH_theta/dsvTauMinusX
    D(8,7) = -D(8,1);                                                     // dH_theta/dsvTauMinusY
    D(8,8) = -D(8,2);                                                     // dH_theta/dsvTauMinusZ
    d(8) =  (tauMinusD - tauMinusDt)/tauMinusDz + (tauMinusPt - tauMinusP)/tauMinusPz;
    // CV: add constraint that recoil = tau+ + tau-
    E( 9, 0) = +1.;                                                       // dH_px/dtauPlusPx
    E( 9, 4) = +1.;                                                       // dH_px/dtauMinusPx
    D( 9, 9) = -1.;                                                       // dH_px/drecoilPx
    d( 9) = tauPlusPx + tauMinusPx - recoilPx;
    E(10, 1) = +1.;                                                       // dH_py/dtauPlusPy
    E(10, 5) = +1.;                                                       // dH_py/dtauMinusPy
    D(10,10) = -1.;                                                       // dH_py/drecoilPy
    d(10) = tauPlusPy + tauMinusPy - recoilPy;
    E(11, 2) = +1.;                                                       // dH_pz/dtauPlusPz
    E(11, 6) = +1.;                                                       // dH_pz/dtauMinusPz
    D(11,11) = -1.;                                                       // dH_pz/drecoilPz
    d(11) = tauPlusPz + tauMinusPz - recoilPz;
    E(12, 3) = +1.;                                                       // dH_energy/dtauPlusE
    E(12, 7) = +1.;                                                       // dH_energy/dtauMinusE
    D(12,12) = -1.;                                                       // dH_energy/drecoilE
    d(12) = tauPlusE  + tauMinusE  - recoilE;
    if ( verbosity_ >= 1 )
    {
      printCovMatrix("V_eta0", V_eta0);
      std::cout << "D:\n";
      std::cout << D << "\n";
      std::cout << "E:\n";
      std::cout << E << "\n";
      std::cout << "d:\n";
      std::cout << d << "\n";
    }

    //-----------------------------------------------------------------------------------------------
    // CV: compute solution to minimization problem
    //     using formulas given in Section 6.1 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
    //     For the syntax of matrix operations,
    //     see Sections "Matrix and vector functions" and "Linear algebra functions"
    //     of the ROOT documentation at
    //       https://root.cern.ch/doc/v608/MatVecFunctions.html 
    //     and
    //       https://root.cern.ch/doc/v608/SMatrixDoc.html 
    //     respectively.
    //-----------------------------------------------------------------------------------------------

    math::MatrixMxC DT = ROOT::Math::Transpose(D);
    
    //math::MatrixCxC Vinv_D = D*V_eta0*DT;
    //int V_D_errorFlag = 0;
    //math::MatrixCxC V_D = Vinv_D.Inverse(V_D_errorFlag);
    //if ( V_D_errorFlag != 0 )
    //{
    //  printCovMatrix("Vinv_D", Vinv_D);
    //  throw cmsException("KinematicFit::operator()", __LINE__)
    //    << "Failed to invert matrix Vinv_D !!\n";
    //}

    // CV: The matrix Vinv_D is rank-deficient and cannot be inverted,
    //     because the first five rows of matrix D (and hence of matrix Vinv_D) are zero.
    //     We fix this by inverting the (C - 5) x (C - 5) submatrix.
    //     The approach is motivated by Section 7.7.3 of the book 
    //       V. Blobel and E. Lohrmann,
    //       "Statistische und numerische Methoden der Datenanalyse",
    //       Teubner, 1998
    math::MatrixCxC Vinv_D = D*V_eta0*DT;
    math::MatrixCm5xCm5 Vinv_D_sub = Vinv_D.Sub<math::MatrixCm5xCm5>(5,5);
    int V_D_sub_errorFlag = 0;
    math::MatrixCm5xCm5 V_D_sub = Vinv_D_sub.Inverse(V_D_sub_errorFlag);
    if ( V_D_sub_errorFlag != 0 )
    {
      printCovMatrix("Vinv_D_sub", Vinv_D_sub);
      throw cmsException("KinematicFit::operator()", __LINE__)
        << "Failed to invert matrix Vinv_D_sub !!\n";
    }
    math::MatrixCxC V_D;
    V_D.Place_at(V_D_sub, 5, 5);
    if ( verbosity_ >= 1 )
    {
      printCovMatrix("V_D", V_D);
    }

    math::MatrixPxC ET = ROOT::Math::Transpose(E);

    math::MatrixPxP Vinv_E = ET*V_D*E;
    int V_E_errorFlag = 0;
    math::MatrixPxP V_E = Vinv_E.Inverse(V_E_errorFlag);
    if ( V_E_errorFlag != 0 )
    {
      printCovMatrix("Vinv_E", Vinv_E);
      throw cmsException("KinematicFit::operator()", __LINE__)
        << "Failed to invert matrix Vinv_E !!\n";
    }

    math::MatrixCxC V_lambda = V_D - V_D*E*V_E*ET*V_D;
    math::VectorC lambda = V_lambda*d;
    math::VectorM deta = -V_eta0*DT*lambda;
    math::VectorP dz = -V_E*ET*V_D*d;
    math::MatrixMxM V_eta = V_eta0 - V_eta0*DT*V_lambda*D*V_eta0;
    math::MatrixPxP V_z = V_E;
    math::MatrixPxM cov_z_eta = -V_E*ET*V_D*D*V_eta0;
    math::MatrixMxP cov_eta_z = ROOT::Math::Transpose(cov_z_eta);

    // CV: compute chi^2
    double chi2 = ROOT::Math::Dot(lambda, Vinv_D*lambda);
    if ( verbosity_ >= 1 )
    {
      std::cout << "chi^2 = " << chi2 << "\n";
    }

    // CV: compute residuals
    math::VectorC residuals = D*deta + E*dz + d;
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
      math::VectorM pulls;
      for ( int idx = 0; idx < numMeasurements; ++idx )
      {      
        pulls(idx) = deta(idx)/std::sqrt(V_eta0(idx,idx));
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
      reco::Candidate::Point pv(pvX + deta(0), pvY + deta(1), pvZ + deta(2));
      kineEvt_kinfit.pv_ = pv;
      kineEvt_kinfit.pvCov_ = V_eta.Sub<math::Matrix3x3>(0,0);
      reco::Candidate::LorentzVector recoilP4(recoilPx + deta(9), recoilPy + deta(10), recoilPz + deta(11), recoilE + deta(12));
      kineEvt_kinfit.recoilP4_ = recoilP4;
      kineEvt_kinfit.recoilCov_ = V_eta.Sub<math::Matrix4x4>(9,9);
      reco::Candidate::LorentzVector tauPlusP4(tauPlusPx + dz(0), tauPlusPy + dz(1), tauPlusPz + dz(2), tauPlusE + dz(3));
      kineEvt_kinfit.tauPlusP4_ = tauPlusP4;
      kineEvt_kinfit.tauPlusP4_isValid_ = true;
      kineEvt_kinfit.tauPlusCov_ = V_z.Sub<math::Matrix4x4>(0,0);
      reco::Candidate::Point svTauPlus(svTauPlusX + deta(3), svTauPlusY + deta(4), svTauPlusZ + deta(5));
      kineEvt_kinfit.svTauPlus_ = svTauPlus;
      kineEvt_kinfit.svTauPlusCov_ = V_eta.Sub<math::Matrix3x3>(3,3);
      reco::Candidate::LorentzVector tauMinusP4(tauMinusPx + dz(4), tauMinusPy + dz(5), tauMinusPz + dz(6), tauMinusE + dz(7));
      kineEvt_kinfit.tauMinusP4_ = tauMinusP4;
      kineEvt_kinfit.tauMinusP4_isValid_ = true;
      kineEvt_kinfit.tauMinusCov_ = V_z.Sub<math::Matrix4x4>(4,4);
      reco::Candidate::Point svTauMinus(svTauMinusX + deta(6), svTauMinusY + deta(7), svTauMinusZ + deta(8));
      kineEvt_kinfit.svTauMinus_ = svTauMinus;
      kineEvt_kinfit.svTauMinusCov_ = V_eta.Sub<math::Matrix3x3>(6,6);
      math::MatrixMpPxMpP cov;
      cov.Place_at(V_eta,                   0,               0);
      cov.Place_at(cov_z_eta, numMeasurements,               0);
      cov.Place_at(cov_eta_z,               0, numMeasurements);
      cov.Place_at(V_z,       numMeasurements, numMeasurements);
      kineEvt_kinfit.kinFitCov_ = cov;
      kineEvt_kinfit.kinFitChi2_ = chi2;
      kineEvt_kinfit.kinFit_isValid_ = true;
      if ( verbosity_ >= 1 )
      {
        std::cout << "constraint equations:\n";
        reco::Candidate::LorentzVector higgsP4 = tauPlusP4 + tauMinusP4;
        std::cout << "Higgs mass = " << higgsP4.mass() << "\n";
        std::cout << "tau+ mass = " << tauPlusP4.mass() << "\n";
        double nuTauPlusPx   = tauPlusPx - visTauPlusPx;
        double nuTauPlusPy   = tauPlusPy - visTauPlusPy;
        double nuTauPlusPz   = tauPlusPz - visTauPlusPz;
        double nuTauPlusE    = tauPlusE  - visTauPlusE;
        std::cout << "neutrino from tau+ decay:" 
                  << " E = " << nuTauPlusE << ", Px = " << nuTauPlusPx << ", Py = " << nuTauPlusPy << ", Pz = " << nuTauPlusPz << ","
                  << " mass = " << (tauPlusP4 - visTauPlusP4).mass() << "\n";
        std::cout << "tau- mass = " << tauMinusP4.mass() << "\n";
        double nuTauMinusPx   = tauMinusPx - visTauMinusPx;
        double nuTauMinusPy   = tauMinusPy - visTauMinusPy;
        double nuTauMinusPz   = tauMinusPz - visTauMinusPz;
        double nuTauMinusE    = tauMinusE  - visTauMinusE;
        std::cout << "neutrino from tau- decay:" 
                  << " E = " << nuTauMinusE << ", Px = " << nuTauMinusPx << ", Py = " << nuTauMinusPy << ", Pz = " << nuTauMinusPz << ","
                  << " mass = " << (tauMinusP4 - visTauMinusP4).mass() << "\n";
        auto tauPlusD3 = svTauPlus - pv;
        std::cout << "phi of tau+: four-vector = " << tauPlusP4.phi() << ", decay vertex = " << tauPlusD3.phi() << "\n";
        std::cout << "theta of tau+: four-vector = " << tauPlusP4.theta() << ", decay vertex = " << tauPlusD3.theta() << "\n";
        std::cout << "phi of tau-: four-vector = " << tauMinusP4.phi() << ", decay vertex = " << tauMinusD3.phi() << "\n";
        std::cout << "theta of tau-: four-vector = " << tauMinusP4.theta() << ", decay vertex = " << tauMinusD3.theta() << "\n";
        std::cout << "Higgs - recoil:"
                  " E = "  << higgsP4.energy() - recoilP4.energy() << ","
                  " Px = " << higgsP4.px()     - recoilP4.px()     << ","
                  " Py = " << higgsP4.py()     - recoilP4.py()     << ","
                  " Pz = " << higgsP4.pz()     - recoilP4.pz()     << "\n";
      }
    }

    double diff = std::sqrt(ROOT::Math::Dot(deta, deta) + ROOT::Math::Dot(dz, dz));
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

