#include "TauAnalysis/Entanglement/interface/KinematicFit.h"

#include "DataFormats/Math/interface/deltaR.h"                    // deltaR2()
#include "DataFormats/Math/interface/Matrix.h"                    // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                    // math::Vector
#include "DataFormats/TauReco/interface/PFTau.h"                  // reco::PFTau::kOneProng0PiZero

#include "TauAnalysis/Entanglement/interface/cmsException.h"      // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"         // mHiggs, mTau
#include "TauAnalysis/Entanglement/interface/KinematicParticle.h" // math::Matrix7x7
#include "TauAnalysis/Entanglement/interface/square.h"            // square()

#include "Math/Functions.h"                                       // ROOT::Math::Dot(), ROOT::Math::Similarity(), ROOT::Math::Transpose() 

#include <cmath>                                                  // std::abs(), std::sqrt()
#include <iostream>                                               // std::cout
#include <string>                                                 // std::string

const int numParameters = 32; 
// CV: the parameters are defined in the following order:
//    (primaryVtx       (3);
//     tauPlusP4        (4)                    ;
       tauPlusDecayVtx  (3);
//     tauMinusP4       (4)                    ;
//     tauMinusDecayVtx (3);
//     recoilP4)        (4))
// where:
//     tauPlusVtx = tauMinusVtx = primary event vertex (PV)
//     visTauPlusVtx            = tau+ decay vertex    (SV+)
//     visTauMinusVtx           = tau- decay vertex    (SV-)
// and energy and momentum components of four-vectors are given in the order:
//    (px, py, pz, E)
// while the position of vertices are given in the order:
//    (x, y, z)
// cf. Section II of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf

const int numConstraints = 16;

namespace math
{
  typedef Matrix<3,3>::type                           Matrix3x3;
  typedef Vector<3>::type                             Vector3;
  
  typedef Matrix<4,4>::type                           Matrix4x4;
  typedef Vector<4>::type                             Vector4;

  typedef Matrix<numParameters,numParameters>::type   MatrixPxP;
  typedef Matrix<numParameters,numConstraints>::type  MatrixPxC;
  typedef Matrix<numConstraints,numParameters>::type  MatrixCxP;
  typedef Matrix<numConstraints,numConstraints>::type MatrixCxC;
  typedef Vector<numParameters>::type                 VectorP;
  typedef Vector<numConstraints>::type                VectorC;
}

KinematicFit::KinematicFit(const edm::ParameterSet& cfg)
  : resolutions_(nullptr)
  , spinAnalyzer_(cfg)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  edm::ParameterSet cfg_resolutions = cfg.getParameterSet("resolutions");
  resolutions_ = new Resolutions(cfg_resolutions);
}

KinematicFit::~KinematicFit()
{
  delete resolutions_;
}

namespace
{
  template <typename T>
  void
  printCovMatrix(const std::string& label, const T& cov)
  {
    std::cout << label << ":\n";
    std::cout << cov << "\n";
    double det = -1.;
    cov.Det2(det);
    std::cout << " det = " << det << "\n";
  }
}

KinematicEvent
KinematicFit::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<KinematicFit::operator()>:\n"; 
  }

  reco::Candidate::LorentzVector tauPlusP4 = kineEvt.get_tauPlusP4();
  double tauPlusPx = tauPlusP4.px();
  double tauPlusPy = tauPlusP4.py();
  double tauPlusPz = tauPlusP4.pz();
  double tauPlusPt = tauPlusP4.pt();
  double tauPlusP  = tauPlusP4.P();
  double tauPlusE  = tauPlusP4.energy();
  // CV: the four-vector of the tau+ is not measured;
  //     we acccount for this by setting the diagonal elements of the covariance matrix to large values,
  //     following the "huge error method" described in Section 6 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
  math::Matrix4x4 tauPlus_cov4x4;
  tauPlus_cov4x4(0,0) = square(1.e+3);
  tauPlus_cov4x4(1,1) = square(1.e+3);
  tauPlus_cov4x4(2,2) = square(1.e+3);
  tauPlus_cov4x4(3,3) = square(1.e+3);
  assert(kineEvt.get_svTauPlus_isValid());
  auto tauPlusD3 = kineEvt.get_svTauPlus() - kineEvt.get_pv();
  double tauPlusDx = tauPlusD3.X();
  double tauPlusDy = tauPlusD3.Y();
  double tauPlusDz = tauPlusD3.Z();
  double tauPlusDt = std::sqrt(tauPlusD3.Perp2());
  double tauPlusD  = std::sqrt(tauPlusD3.Mag2());
  if ( verbosity_ >= 1 )
  {
    printCovMatrix("tauPlus_cov4x4", tauPlus_cov4x4);
  }

  reco::Candidate::LorentzVector tauMinusP4 = kineEvt.get_tauMinusP4();
  double tauMinusPx = tauMinusP4.px();
  double tauMinusPy = tauMinusP4.py();
  double tauMinusPz = tauMinusP4.pz();
  double tauMinusPt = tauMinusP4.pt();
  double tauMinusP  = tauMinusP4.P();
  double tauMinusE  = tauMinusP4.energy();
  // CV: the four-vector of the tau- is not measured;
  //     we acccount for this by setting the diagonal elements of the covariance matrix to large values,
  //     following the "huge error method" described in Section 6 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
  //    (same as for tau+)
  math::Matrix4x4 tauMinus_cov4x4 = tauPlus_cov4x4;
  assert(kineEvt.get_svTauMinus_isValid());
  auto tauMinusD3 = kineEvt.get_svTauMinus() - kineEvt.get_pv();
  double tauMinusDx = tauMinusD3.X();
  double tauMinusDy = tauMinusD3.Y();
  double tauMinusDz = tauMinusD3.Z();
  double tauMinusDt = std::sqrt(tauMinusD3.Perp2());
  double tauMinusD  = std::sqrt(tauMinusD3.Mag2());
  if ( verbosity_ >= 1 )
  {
    printCovMatrix("tauMinus_cov4x4", tauMinus_cov4x4);
  }

  math::Matrix3x3 pv_cov;
  pv_cov(0,0) = square(resolutions_->get_pvResolution_xy());
  pv_cov(1,1) = square(resolutions_->get_pvResolution_xy());
  pv_cov(2,2) = square(resolutions_->get_pvResolution_z());
  if ( verbosity_ >= 1 )
  {
    printCovMatrix("pv_cov", pv_cov);
  }

  reco::Candidate::LorentzVector visTauPlusP4 = kineEvt.get_visTauPlusP4();
  double visTauPlusPx = visTauPlusP4.px();
  double visTauPlusPy = visTauPlusP4.py();
  double visTauPlusPz = visTauPlusP4.pz();
  double visTauPlusE  = visTauPlusP4.energy();
  const std::vector<KinematicParticle>& daughtersTauPlus = kineEvt.get_daughtersTauPlus();
  math::Matrix4x4 visTauPlus_cov4x4;
  math::Matrix3x3 visTauPlus_cov3x3;
  int idxDaughterTauPlus = 0;
  bool visTauPlus_isFirst = true;
  for ( const KinematicParticle& daughter : daughtersTauPlus )
  {
    int daughter_absPdgId = std::abs(daughter.get_pdgId());
    if ( daughter_absPdgId == 12 || daughter_absPdgId == 14 || daughter_absPdgId == 16 ) continue;
    const math::Matrix7x7& cov7x7 = daughter.get_cov7x7();
    if ( verbosity_ >= 1 )
    {
      std::cout << "daughterTauPlus #" << idxDaughterTauPlus << "\n";
      printCovMatrix("cov7x7", cov7x7);
    }
    // CV: compute 4x4 covariance matrix of the tauh four-vector by summing the covariance matrices of all visible tau decay products;
    //     for the 3x3 covariance matrix of the tauh vertex position, simply take covariance matrix of first tau decay product
    //    (assuming the vertex positions of all tau decay products to be equal)
    //     For the syntax of retrieving a small covariance matrix from a larger one, 
    //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
    visTauPlus_cov4x4 += cov7x7.Sub<math::Matrix4x4>(0,0);
    if ( visTauPlus_isFirst )
    {
      visTauPlus_cov3x3 = cov7x7.Sub<math::Matrix3x3>(4,4);
      visTauPlus_isFirst = false;
    }
    ++idxDaughterTauPlus;
  }
  // CV: add uncertainty on mass of visible decay products.
  //     This uncertainty needs to be added even if tau is reconstructed in OneProng0Pi0 decay mode,
  //     to avoid that covariance matrix has rank zero
  //    (matrices of rank zero cannot be inverted !!)
  const double visTauPlus_sigma2_mass = square(1.e-1);
  visTauPlus_cov4x4(3,3) += square(visTauPlusP4.mass()/visTauPlusP4.energy())*visTauPlus_sigma2_mass;
  if ( verbosity_ >= 1 )
  {
    printCovMatrix("visTauPlus_cov4x4", visTauPlus_cov4x4);
    printCovMatrix("visTauPlus_cov3x3", visTauPlus_cov3x3);
  }

  reco::Candidate::LorentzVector visTauMinusP4 = kineEvt.get_visTauMinusP4();
  double visTauMinusPx = visTauMinusP4.px();
  double visTauMinusPy = visTauMinusP4.py();
  double visTauMinusPz = visTauMinusP4.pz();
  double visTauMinusE  = visTauMinusP4.energy();
  const std::vector<KinematicParticle>& daughtersTauMinus = kineEvt.get_daughtersTauMinus();
  math::Matrix4x4 visTauMinus_cov4x4;
  math::Matrix3x3 visTauMinus_cov3x3;
  int idxDaughterTauMinus = 0;
  bool visTauMinus_isFirst = true;
  for ( const KinematicParticle& daughter : daughtersTauMinus )
  {
    int daughter_absPdgId = std::abs(daughter.get_pdgId());
    if ( daughter_absPdgId == 12 || daughter_absPdgId == 14 || daughter_absPdgId == 16 ) continue;
    const math::Matrix7x7& cov7x7 = daughter.get_cov7x7();
    if ( verbosity_ >= 1 )
    {
      std::cout << "daughterTauMinus #" << idxDaughterTauMinus << "\n";
      printCovMatrix("cov7x7", cov7x7);
    }
    // CV: compute 4x4 covariance matrix of the tauh four-vector by summing the covariance matrices of all visible tau decay products;
    //     for the 3x3 covariance matrix of the tauh vertex position, simply take covariance matrix of first tau decay product
    //    (assuming the vertex positions of all tau decay products to be equal)
    //     For the syntax of retrieving a small covariance matrix from a larger one, 
    //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
    visTauMinus_cov4x4 += cov7x7.Sub<math::Matrix4x4>(0,0);
    if ( visTauMinus_isFirst )
    {
      visTauMinus_cov3x3 = cov7x7.Sub<math::Matrix3x3>(4,4);
      visTauMinus_isFirst = false;
    }
    ++idxDaughterTauMinus;
  }
  // CV: add uncertainty on mass of visible decay products.
  //     This uncertainty needs to be added even if tau is reconstructed in OneProng0Pi0 decay mode,
  //     to avoid that covariance matrix has rank zero
  //    (matrices of rank zero cannot be inverted !!)
  const double visTauMinus_sigma2_mass = square(1.e-1);
  visTauMinus_cov4x4(3,3) += square(visTauMinusP4.mass()/visTauMinusP4.energy())*visTauMinus_sigma2_mass;
  if ( verbosity_ >= 1 )
  {
    printCovMatrix("visTauMinus_cov4x4", visTauMinus_cov4x4);
    printCovMatrix("visTauMinus_cov3x3", visTauMinus_cov3x3);
  }

  math::Matrix4x4 recoil_cov;
  recoil_cov(0,0) = square(resolutions_->get_recoilResolution_px());
  recoil_cov(1,1) = square(resolutions_->get_recoilResolution_py());
  recoil_cov(2,2) = square(resolutions_->get_recoilResolution_pz());
  recoil_cov(3,3) = square(resolutions_->get_recoilResolution_energy());
  if ( verbosity_ >= 1 )
  {
    printCovMatrix("recoil_cov", recoil_cov);
  }

  math::MatrixPxP V_alpha0;
  // CV: add resolution on tau+ four-vector
  V_alpha0.Place_at(tauPlus_cov4x4    , 0*7 + 0, 0*7 + 0);
  // CV: add resolution on tau- four-vector
  V_alpha0.Place_at(tauMinus_cov4x4   , 1*7 + 0, 1*7 + 0);
  // CV: add resolution on production vertex of tau+ and tau-
  //     Note that the production vertex of tau+ and tau- equal the PV and are thus 100% correlated;
  //     in order to avoid that the resolution on the primary event vertex constrains the tau+ and tau- production vertex twice,
  //     we scale the covariance matrix of the PV by a factor two.
  //     For the syntax of "embedding" a small covariance matrix into a larger one, 
  //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
  V_alpha0.Place_at(2.*pv_cov         , 0*7 + 4, 0*7 + 4);
  // CV: ignore correlation between position of tau+ and tau- production vertex,
  //     as covariance matrix attains rank zero and cannot be inverted otherwise !!
  //V_alpha0.Place_at(2.*pv_cov         , 0*7 + 4, 1*7 + 4);
  //V_alpha0.Place_at(2.*pv_cov         , 1*7 + 4, 0*7 + 4);
  V_alpha0.Place_at(2.*pv_cov         , 1*7 + 4, 1*7 + 4);
  // CV: add resolutions on four-vector of visible decay products and on decay vertex of tau+
  V_alpha0.Place_at(visTauPlus_cov4x4 , 2*7 + 0, 2*7 + 0);
  V_alpha0.Place_at(visTauPlus_cov3x3 , 2*7 + 4, 2*7 + 4);
  // CV: add resolutions on four-vector of visible decay products and on decay vertex of tau-
  V_alpha0.Place_at(visTauMinus_cov4x4, 3*7 + 0, 3*7 + 0);
  V_alpha0.Place_at(visTauMinus_cov3x3, 3*7 + 4, 3*7 + 4);
  // CV: add resolutions on four-vector of recoil
  V_alpha0.Place_at(recoil_cov        , 4*7 + 0, 4*7 + 0);
  
  math::MatrixCxP D;
  math::VectorC   d;
  // CV: add Higgs mass constraint
  D( 0,0*7 + 0) = -2.*(tauPlusPx + tauMinusPx);
  D( 0,0*7 + 1) = -2.*(tauPlusPy + tauMinusPy);
  D( 0,0*7 + 2) = -2.*(tauPlusPz + tauMinusPz);
  D( 0,0*7 + 3) = +2.*(tauPlusE  + tauMinusE );
  D( 0,1*7 + 0) = -2.*(tauPlusPx + tauMinusPx);
  D( 0,1*7 + 1) = -2.*(tauPlusPy + tauMinusPy);
  D( 0,1*7 + 2) = -2.*(tauPlusPz + tauMinusPz);
  D( 0,1*7 + 3) = +2.*(tauPlusE  + tauMinusE );
  d( 0) = (tauPlusP4 + tauMinusP4).mass() - mHiggs;
  // CV: add tau+ mass constraint
  D( 1,0*7 + 0) = -2.*tauPlusPx;
  D( 1,0*7 + 1) = -2.*tauPlusPy;
  D( 1,0*7 + 2) = -2.*tauPlusPz;
  D( 1,0*7 + 3) = +2.*tauPlusE;
  d( 1) = tauPlusP4.mass() - mTau;
  // CV: add tau- mass constraint
  D( 2,1*7 + 0) = -2.*tauMinusPx;
  D( 2,1*7 + 1) = -2.*tauMinusPy;
  D( 2,1*7 + 2) = -2.*tauMinusPz;
  D( 2,1*7 + 3) = +2.*tauMinusE;
  d( 2) = tauMinusP4.mass() - mTau;
  // CV: add neutrino mass constraint for tau+
  D( 3,1*7 + 0) = -2.*(tauPlusPx  - visTauPlusPx );
  D( 3,1*7 + 1) = -2.*(tauPlusPy  - visTauPlusPy );
  D( 3,1*7 + 2) = -2.*(tauPlusPz  - visTauPlusPz );
  D( 3,1*7 + 3) = +2.*(tauPlusE   - visTauPlusE  );
  D( 3,3*7 + 0) = +2.*(tauPlusPx  - visTauPlusPx );
  D( 3,3*7 + 1) = +2.*(tauPlusPy  - visTauPlusPy );
  D( 3,3*7 + 2) = +2.*(tauPlusPz  - visTauPlusPz );
  D( 3,3*7 + 3) = -2.*(tauPlusE   - visTauPlusE  );
  d( 3) = (tauPlusP4 - visTauPlusP4).mass();
  // CV: add neutrino mass constraint for tau-
  D( 4,2*7 + 0) = -2.*(tauMinusPx - visTauMinusPx);
  D( 4,2*7 + 1) = -2.*(tauMinusPy - visTauMinusPy);
  D( 4,2*7 + 2) = -2.*(tauMinusPz - visTauMinusPz);
  D( 4,2*7 + 3) = +2.*(tauMinusE  - visTauMinusE );
  D( 4,4*7 + 0) = +2.*(tauMinusPx - visTauMinusPx);
  D( 4,4*7 + 1) = +2.*(tauMinusPy - visTauMinusPy);
  D( 4,4*7 + 2) = +2.*(tauMinusPz - visTauMinusPz);
  D( 4,4*7 + 3) = -2.*(tauMinusE  - visTauMinusE );
  d( 4) = (tauMinusP4 - visTauMinusP4).mass();
  // CV: add constraints that tau+ and tau- originate from the same vertex
  D( 5,0*7 + 4) = +1.;
  D( 5,1*7 + 4) = -1.;
  d( 5) = 0.;
  D( 6,0*7 + 5) = +1.;
  D( 6,1*7 + 5) = -1.;
  d( 6) = 0.;
  D( 7,0*7 + 6) = +1.;
  D( 7,1*7 + 6) = -1.;
  d( 7) = 0.;
  // CV: add "parallelism" constraints;
  //     cf. Section 4.1.3.3 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
  D( 8,0*7 + 0) =  (1. - tauPlusPx/tauPlusPt)*tauPlusPy;
  D( 8,0*7 + 1) = -(tauPlusPx/tauPlusPy)*D(8,0*7 + 0);
  D( 8,0*7 + 4) =  (1. - tauPlusDx/tauPlusDt)/tauPlusDy;
  D( 8,0*7 + 5) = -(tauPlusDx/tauPlusDy)*D(8,0*7 + 4);
  D( 8,2*7 + 4) = -D( 8,0*7 + 4);
  D( 8,2*7 + 5) = -D( 8,0*7 + 5);
  d( 8)         =  (tauPlusDt - tauPlusDx)/tauPlusDy + (tauPlusPx - tauPlusPt)/tauPlusPy;
  D( 9,0*7 + 0) =  (tauPlusPx/tauPlusPz)*(1./tauPlusPt - 1./tauPlusP);
  D( 9,0*7 + 1) =  (tauPlusPy/tauPlusPz)*(1./tauPlusPt - 1./tauPlusP);
  D( 9,0*7 + 2) = -1./tauPlusP + (tauPlusP - tauPlusPt)/square(tauPlusPz);
  D( 9,0*7 + 4) =  (tauPlusDx/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);
  D( 9,0*7 + 5) =  (tauPlusDy/tauPlusDz)*(1./tauPlusDt - 1./tauPlusD);
  D( 9,0*7 + 6) = -1./tauPlusD + (tauPlusD - tauPlusDt)/square(tauPlusDz);
  D( 9,2*7 + 4) = -D( 9,0*7 + 4);
  D( 9,2*7 + 5) = -D( 9,0*7 + 5);
  D( 9,2*7 + 6) = -D( 9,0*7 + 6);
  d( 9)         =  (tauPlusD - tauPlusDt)/tauPlusDz + (tauPlusPt - tauPlusP)/tauPlusPz;
  D(10,1*7 + 0) =  (1. - tauMinusPx/tauMinusPt)*tauMinusPy;
  D(10,1*7 + 1) = -(tauMinusPx/tauMinusPy)*D(8,0*7 + 0);
  D(10,1*7 + 4) =  (1. - tauMinusDx/tauMinusDt)/tauMinusDy;
  D(10,1*7 + 5) = -(tauMinusDx/tauMinusDy)*D(8,0*7 + 4);
  D(10,3*7 + 4) = -D(10,1*7 + 4);
  D(10,3*7 + 5) = -D(10,1*7 + 5);
  d(10)         =  (tauMinusDt - tauMinusDx)/tauMinusDy + (tauMinusPx - tauMinusPt)/tauMinusPy;
  D(11,1*7 + 0) =  (tauMinusPx/tauMinusPz)*(1./tauMinusPt - 1./tauMinusP);
  D(11,1*7 + 1) =  (tauMinusPy/tauMinusPz)*(1./tauMinusPt - 1./tauMinusP);
  D(11,1*7 + 2) = -1./tauMinusP + (tauMinusP - tauMinusPt)/square(tauMinusPz);
  D(11,1*7 + 4) =  (tauMinusDx/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);
  D(11,1*7 + 5) =  (tauMinusDy/tauMinusDz)*(1./tauMinusDt - 1./tauMinusD);
  D(11,1*7 + 6) = -1./tauMinusD + (tauMinusD - tauMinusDt)/square(tauMinusDz);
  D(11,3*7 + 4) = -D(11,1*7 + 4);
  D(11,3*7 + 5) = -D(11,1*7 + 5);
  D(11,3*7 + 6) = -D(11,1*7 + 6);
  d(11)         =  (tauMinusD - tauMinusDt)/tauMinusDz + (tauMinusPt - tauMinusP)/tauMinusPz;
  // CV: add constraint on four-vector of recoil 
  //    (idx=0: px, idx=1: py, idx=2: pz, idx=3: energy)
  for ( int idx = 0; idx < 4; ++idx )
  {
    D(12 + idx,      idx) = +1.;
    D(12 + idx,1*7 + idx) = +1.;
    D(12 + idx,4*7 + idx) = -1.;
  }

  if ( verbosity_ >= 1 )
  {
    std::cout << "V_alpha0:\n";
    std::cout << V_alpha0 << "\n";
    double det = -1.;
    V_alpha0.Det2(det);
    std::cout << " det = " << det << "\n";
    std::cout << "D:\n";
    std::cout << D << "\n";
    std::cout << "d:\n";
    std::cout << d << "\n";
  }

  //-------------------------------------------------------------------------------------------------
  // CV: compute solution to minimization problem
  //     using formulas given in Section 3 of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
  //     see Section "Matrix and vector functions" of the ROOT documentation https://root.cern.ch/doc/v608/MatVecFunctions.html 
  //     and Section "Linear algebra functions" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html 
  //     for the syntax of matrix operations

  math::MatrixPxC DT = ROOT::Math::Transpose(D);

  math::MatrixCxC Vinv_D = D*V_alpha0*DT;
  if ( verbosity_ >= 1 )
  {
    std::cout << "Vinv_D:\n";
    std::cout << Vinv_D << "\n";
    double det = -1.;
    Vinv_D.Det2(det);
    std::cout << " det = " << det << "\n";
  }

  int errorFlag = 0;
  math::MatrixCxC V_D = Vinv_D.Inverse(errorFlag);
  if ( errorFlag != 0 )
    throw cmsException("KinematicFit::operator()", __LINE__)
      << "Failed to invert matrix Vinv_D !!\n";

  // CV: compute matrix lambda.
  //     Note that deltaalpha0 is zero, since the initial parameter vector and the expansion point are the same, i.e. alpha0 = alphaA
  math::VectorC lambda = V_D*d;

  // CV: compute change in parameters alpha,
  //     where dalpha := alpha - alpha0
  math::VectorP dalpha = -V_alpha0*DT*lambda;
  if ( verbosity_ >= 1 )
  {
    std::cout << "dalpha:\n";
    std::cout << dalpha << "\n";
  }

  // CV: compute covariance matrix on parameters alpha
  math::MatrixPxP V_alpha = V_alpha0 - V_alpha0*DT*V_D*D*V_alpha0;
  if ( verbosity_ >= 1 )
  {
    std::cout << "V_alpha:\n";
    std::cout << V_alpha << "\n";
  }

  // CV: compute chi^2
  double chi2 = ROOT::Math::Dot(lambda, d);
  if ( verbosity_ >= 1 )
  {
    std::cout << "chi^2 = " << chi2 << "\n";

    math::VectorP pulls;
    for ( int idx = 0; idx < numParameters; ++idx )
    {     
      pulls(idx) = dalpha(idx)/std::sqrt(V_alpha0(idx,idx));
    }
    std::cout << "pulls:\n";
    std::cout << pulls << "\n";

    math::VectorC residuals = D*dalpha + d;
    std::cout << "residuals of constraint equations:\n";
    std::cout << residuals << "\n";
    reco::Candidate::LorentzVector tauPlusP4(tauPlusPx + dalpha(0), tauPlusPy + dalpha(1), tauPlusPz + dalpha(2), tauPlusE + dalpha(3));
    reco::Candidate::LorentzVector tauMinusP4(tauMinusPx + dalpha(7), tauMinusPy + dalpha(8), tauMinusPz + dalpha(9), tauMinusE + dalpha(10));
    std::cout << "tau-pair mass = " << (tauPlusP4 + tauMinusP4).mass() << "\n";
  }
  //-------------------------------------------------------------------------------------------------

  assert(0);

  // TO BE IMPLEMENTED:
  //  1) run the code above iteratively until convergence
  //  2) store result of kinematic fit in KinematicEvent class

  KinematicEvent kineEvt_fitted = kineEvt;

  return kineEvt_fitted;
}

