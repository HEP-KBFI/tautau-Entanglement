#include "TauAnalysis/Entanglement/interface/KinFitConstraint.h"

#include "DataFormats/Math/interface/deltaPhi.h"                                 // deltaPhi()
#include "DataFormats/Math/interface/Vector.h"                                   // math::Vector

#include "TauAnalysis/Entanglement/interface/comp_mT.h"                          // comp_mT()
#include "TauAnalysis/Entanglement/interface/comp_nuPz.h"                        // build_nuP4(), comp_nuPz()
#include "TauAnalysis/Entanglement/interface/constants.h"                        // kLHC, kSuperKEKB, mHiggs, mTau
#include "TauAnalysis/Entanglement/interface/KinFitParameters_and_Constraints.h" // kinFit::numParameters
#include "TauAnalysis/Entanglement/interface/printVector.h"                      // printVector()
#include "TauAnalysis/Entanglement/interface/rotateVector.h"                     // rotateVector()
#include "TauAnalysis/Entanglement/interface/square.h"                           // square()

#include <cmath>                                                                 // std::sqrt()

KinFitConstraint::KinFitConstraint(int collider, const KinematicEvent& kineEvt, double signTauPlus, double signTauMinus, int verbosity)
  : KinFitConstraintBase(kinFit::numParameters, collider == kSuperKEKB ? kinFit::numConstraints_SuperKEKB : kinFit::numConstraints_LHC, 2, verbosity)
  , collider_(collider)
  , kineEvt_(kineEvt)
  , signTauPlus_(signTauPlus)
  , signTauMinus_(signTauMinus)
{
  d_eq_metric_(0,0) = 1.e+2;                                    // H_Pphi(tau+)
  d_eq_metric_(1,1) = 1.e+2;                                    // H_Ptheta(tau+)
  d_eq_metric_(2,2) = 1.e+2;                                    // H_Pphi(tau-)
  d_eq_metric_(3,3) = 1.e+2;                                    // H_Ptheta(tau-)
  d_eq_metric_(4,4) = ( collider_ == kSuperKEKB ) ? 1. : 1.e-2; // H_px(recoil)
  d_eq_metric_(5,5) = ( collider_ == kSuperKEKB ) ? 1. : 1.e-2; // H_py(recoil)
  d_eq_metric_(6,6) = ( collider_ == kSuperKEKB ) ? 1. : 1.e-2; // H_pz(recoil)
  d_eq_metric_(7,7) = ( collider_ == kSuperKEKB ) ? 1. : 1.e-2; // H_energy(recoil)
  if ( collider_ != kSuperKEKB )
  {
    d_eq_metric_(8,8) = 1.e-4;                                  // H(Higgs mass)
  }

  d_ineq_metric_(0,0) = 1.;                                     // H_mT(tau+)
  d_ineq_metric_(1,1) = 1.;                                     // H_mT(tau-)
}

KinFitConstraint::~KinFitConstraint()
{}

void
KinFitConstraint::set_alphaA(const TVectorD& alphaA)
{
  // CV: The parameters alphaA at the expansion point are defined in the following order:
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

  pv0_ = kineEvt_.pv();
  pv_ = reco::Candidate::Point(alphaA(0), alphaA(1), alphaA(2));

  // CV: Reconstruct event kinematics
  visTauPlusP4_ = kineEvt_.visTauPlusP4();
  double nuTauPlusPx = alphaA(3);
  double nuTauPlusPy = alphaA(4);
  bool nuTauPlus_errorFlag = false;
  double nuTauPlusPz = comp_nuPz(visTauPlusP4_, nuTauPlusPx, nuTauPlusPy, signTauPlus_, nuTauPlus_dPzdPx_, nuTauPlus_dPzdPy_, nuTauPlus_errorFlag, verbosity_);
  nuTauPlusP4_ = build_nuP4(nuTauPlusPx, nuTauPlusPy, nuTauPlusPz);
  tauPlusP4_ = visTauPlusP4_ + nuTauPlusP4_;
  tauPlusD3_rnk_ = reco::Candidate::Vector(alphaA(5), alphaA(6), alphaA(7));

  visTauMinusP4_ = kineEvt_.visTauMinusP4();
  double nuTauMinusPx = alphaA(8);
  double nuTauMinusPy = alphaA(9);
  bool nuTauMinus_errorFlag = false;
  double nuTauMinusPz = comp_nuPz(visTauMinusP4_, nuTauMinusPx, nuTauMinusPy, signTauMinus_, nuTauMinus_dPzdPx_, nuTauMinus_dPzdPy_, nuTauMinus_errorFlag, verbosity_);
  nuTauMinusP4_ = build_nuP4(nuTauMinusPx, nuTauMinusPy, nuTauMinusPz);
  tauMinusP4_ = visTauMinusP4_ + nuTauMinusP4_;
  tauMinusD3_rnk_ = reco::Candidate::Vector(alphaA(10), alphaA(11), alphaA(12));

  recoilP4_ = reco::Candidate::LorentzVector(alphaA(13), alphaA(14), alphaA(15), alphaA(16));

  // CV: Initialize matrix of derrivatives dH/dalpha of constraint equations
  //     and values H(alphaA) of constraint equations at expansion point alphaA.
  // 
  //     The constraints are defined in the following order:
  //       "parallelism" constraint for tau+ [1] (2)
  //       "parallelism" constraint for tau- [1] (2) 
  //       constraint that recoil = tau+ + tau-  (4)
  //       Higgs mass constraint                 (1)
  //     Note that the Higgs mass constraintis applied at the LHC,
  //     but not at the SuperKEKB collider (Belle)
  //  [1] cf. Section 4.1.3.3 of https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf

  // CV: add "parallelism" constraint for tau+ and tau-
  set_parlConstraint("tau+", 0, 3, 
    tauPlusP4_, tauPlusD3_rnk_,
    kineEvt_.tauPlus_r(), kineEvt_.tauPlus_n(), kineEvt_.tauPlus_k(), kineEvt_.tauPlus_rotMatrix_rnk2xyz(),
    nuTauPlusP4_, nuTauPlus_dPzdPx_, nuTauPlus_dPzdPy_, signTauPlus_);
  set_parlConstraint("tau-", 2, 8, 
    tauMinusP4_, tauMinusD3_rnk_,
    kineEvt_.tauMinus_r(), kineEvt_.tauMinus_n(), kineEvt_.tauMinus_k(), kineEvt_.tauMinus_rotMatrix_rnk2xyz(),
    nuTauMinusP4_, nuTauMinus_dPzdPx_, nuTauMinus_dPzdPy_, signTauMinus_);

  // CV: add constraint that recoil = tau+ + tau-
  set_recoilConstraint();

  // CV: add Higgs mass constraint (not for SuperKEKB!)
  if ( collider_ != kSuperKEKB )
  {
    set_higgsMassConstraint();
  }

  // CV: add inequality constraints to ensure that computation of neutrino Pz yields a physical solution
  //     The inequality constraints are expected to be of the form f(parameters) <= 0. 
  //     The case f(parameters) >= 0 needs to be transformed into the case f(parameters) <= 0 
  //     by multiplication of the constraint equation f(parameters) by -1.
  set_mTConstraint("tau+", 0, 3, 
    visTauPlusP4_,
    nuTauPlusP4_);
  set_mTConstraint("tau-", 1, 8,
    visTauMinusP4_,
    nuTauMinusP4_);
  
  errorFlag_ = nuTauPlus_errorFlag || nuTauMinus_errorFlag;    
  if ( verbosity_ >= 2 )
  {
    std::cout << "errorFlag = " << errorFlag_ << "\n";
  }
}

reco::Candidate::Vector
KinFitConstraint::get_flightlength_up(unsigned int idxPar, unsigned int idxOffsetP, double epsilon,
                                      const reco::Candidate::Vector& tauD3_rnk,
                                      const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k, 
                                      const math::Matrix3x3& rotMatrix_rnk2xyz)
{
  reco::Candidate::Vector flightlength_up;
  if ( idxPar == 0 || idxPar == 1 || idxPar == 2 )
  {
    reco::Candidate::Vector tauD3_xyz = rotateVector(tauD3_rnk, rotMatrix_rnk2xyz);
    auto sv = tauD3_xyz + pv0_;
    double pvX = pv_.x();
    if ( idxPar == 0 ) pvX += epsilon;
    double pvY = pv_.y();
    if ( idxPar == 1 ) pvY += epsilon;
    double pvZ = pv_.z();
    if ( idxPar == 2 ) pvZ += epsilon;
    flightlength_up = sv - reco::Candidate::Vector(pvX, pvY, pvZ);
  }
  else if ( idxPar == (idxOffsetP+2) || idxPar == (idxOffsetP+3) || idxPar == (idxOffsetP+4) )
  {
    double tauD3_r = tauD3_rnk.x();
    if ( idxPar == (idxOffsetP+2) ) tauD3_r += epsilon;
    double tauD3_n = tauD3_rnk.y();
    if ( idxPar == (idxOffsetP+3) ) tauD3_n += epsilon;
    double tauD3_k = tauD3_rnk.z();
    if ( idxPar == (idxOffsetP+4) ) tauD3_k += epsilon;
    reco::Candidate::Vector tauD3_xyz = rotateVector(reco::Candidate::Vector(tauD3_r, tauD3_n, tauD3_k), rotMatrix_rnk2xyz);
    auto sv = tauD3_xyz + pv0_;
    flightlength_up = sv - reco::Candidate::Vector(pv_.x(), pv_.y(), pv_.z());
  }
  return flightlength_up;
}

reco::Candidate::LorentzVector
KinFitConstraint::get_tauP4_up(unsigned int idxPar, unsigned int idxOffsetP, double epsilon,
                               const reco::Candidate::LorentzVector& tauP4,
                               const reco::Candidate::LorentzVector& nuP4, double sign)
{
  reco::Candidate::LorentzVector tauP4_up;
  if ( idxPar == idxOffsetP || idxPar == (idxOffsetP+1) )
  {
    double nuPx = nuP4.px();
    if ( idxPar ==  idxOffsetP    ) nuPx += epsilon;
    double nuPy = nuP4.py();
    if ( idxPar == (idxOffsetP+1) ) nuPy += epsilon;
    reco::Candidate::LorentzVector visTauP4 = tauP4 - nuP4;
    double dummy;
    bool errorFlag;
    double nuPz = comp_nuPz(visTauP4, nuPx, nuPy, sign, dummy, dummy, errorFlag);      
    tauP4_up = visTauP4 + build_nuP4(nuPx, nuPy, nuPz);
  }
  return tauP4_up;
}

double
KinFitConstraint::get_Dnum_H_Pphi(unsigned int idxPar, unsigned int idxOffsetC, unsigned int idxOffsetP, 
                                  const reco::Candidate::LorentzVector& tauP4,
                                  const reco::Candidate::Vector& tauD3_rnk,
                                  const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k, 
                                  const math::Matrix3x3& rotMatrix_rnk2xyz,
                                  const reco::Candidate::LorentzVector& nuP4, double sign)
{
  double H_nom = d_eq_(idxOffsetC);
  double H_up  = H_nom;
  double epsilon = 1.;
  if ( idxPar == 0 || idxPar == 1 || idxPar == (idxOffsetP+2) || idxPar == (idxOffsetP+3) || idxPar == (idxOffsetP+4) )
  {
    epsilon = 1.e-4;
    reco::Candidate::Vector flightlength_up = get_flightlength_up(idxPar, idxOffsetP, epsilon, tauD3_rnk, r, n, k, rotMatrix_rnk2xyz);
    H_up = flightlength_up.phi() - tauP4.phi();
  }
  else if ( idxPar == idxOffsetP || idxPar == (idxOffsetP+1) )
  {
    epsilon = 1.e-2;
    reco::Candidate::LorentzVector tauP4_up = get_tauP4_up(idxPar, idxOffsetP, epsilon, tauP4, nuP4, sign);
    H_up -= (tauP4_up.phi() - tauP4.phi());
  }
  return (H_up - H_nom)/epsilon;
}

double
KinFitConstraint::get_Dnum_H_Ptheta(unsigned int idxPar, unsigned int idxOffsetC, unsigned int idxOffsetP, 
                                    const reco::Candidate::LorentzVector& tauP4,
                                    const reco::Candidate::Vector& tauD3_rnk,
                                    const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k, 
                                    const math::Matrix3x3& rotMatrix_rnk2xyz,
                                    const reco::Candidate::LorentzVector& nuP4, double sign)
{
  double H_nom = d_eq_(idxOffsetC);
  double H_up  = H_nom;
  double epsilon = 1.;
  if ( idxPar == 0 || idxPar == 1 || idxPar == 2 || idxPar == (idxOffsetP+2) || idxPar == (idxOffsetP+3) || idxPar == (idxOffsetP+4) )
  {
    epsilon = 1.e-4;
    reco::Candidate::Vector flightlength_up = get_flightlength_up(idxPar, idxOffsetP, epsilon, tauD3_rnk, r, n, k, rotMatrix_rnk2xyz);
    H_up = flightlength_up.theta() - tauP4.theta();
  }
  else if ( idxPar == idxOffsetP || idxPar == (idxOffsetP+1) )
  {
    epsilon = 1.e-2;
    reco::Candidate::LorentzVector tauP4_up = get_tauP4_up(idxPar, idxOffsetP, epsilon, tauP4, nuP4, sign);
    H_up -= (tauP4_up.theta() - tauP4.theta());
  }
  return (H_up - H_nom)/epsilon;
}

void
KinFitConstraint::set_parlConstraint(const std::string& label, unsigned int idxOffsetC, unsigned int idxOffsetP,
                                     const reco::Candidate::LorentzVector& tauP4, 
                                     const reco::Candidate::Vector& tauD3_rnk,
                                     const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k, 
                                     const math::Matrix3x3& rotMatrix_rnk2xyz,
                                     const reco::Candidate::LorentzVector& nuP4, double nu_dPzdPx, double nu_dPzdPy,
                                     double sign)
{
  if ( verbosity_ >= 3 )
  {
    std::cout << "<set_parlConstraint (" << label << ")>:\n";
  }

  double tauPx = tauP4.px();
  double tauPy = tauP4.py();
  double tauPz = tauP4.pz();
  double tauPt = tauP4.pt();
  double tauPt2 = square(tauPt);
  double tauP2  = square(tauP4.P());

  reco::Candidate::Vector tauD3_xyz = rotateVector(tauD3_rnk, rotMatrix_rnk2xyz);
  auto sv = tauD3_xyz + pv0_;
  auto flightlength = sv - pv_;
  double flightlengthX   = flightlength.x();
  double flightlengthY   = flightlength.y();
  double flightlengthZ   = flightlength.z();    
  double flightlengthDt2 = flightlength.perp2();
  double flightlengthDt  = std::sqrt(flightlengthDt2);
  double flightlengthD2  = flightlength.mag2();

  // CV: add "parallelism" constraint for tau+ and tau-
  //     The constraint equation H_Pphi (H_Ptheta) demands that the azimuthal (polar) angle of the tau momentum 
  //     equals the azimuthal (polar) angle of the flightlength = sv - pv.
  //     The small correction for the helix bend between tau lepton production and decay in Eq. (4.13) in Section 4.1.3.3 of [1] is not used, 
  //     as this correction is expected to be on the level of 50 microrad, which we consider negligible.
  //
  //   [1] https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
  D_eq_(idxOffsetC  ,0)            =   flightlengthY/flightlengthDt2;                               // dH_Pphi/dpvX
  D_eq_(idxOffsetC  ,1)            =  -flightlengthX/flightlengthDt2;                               // dH_Pphi/dpvY
  D_eq_(idxOffsetC  ,idxOffsetP)   =   tauPy/square(tauPt);                                         // dH_Pphi/dnuPx
  D_eq_(idxOffsetC  ,idxOffsetP+1) =  -tauPx/square(tauPt);                                         // dH_Pphi/dnuPy     
  D_eq_(idxOffsetC  ,idxOffsetP+2) = (-r.x()*flightlengthY + r.y()*flightlengthX)/flightlengthDt2;  // dH_Pphi/dsvR
  D_eq_(idxOffsetC  ,idxOffsetP+3) = (-n.x()*flightlengthY + n.y()*flightlengthX)/flightlengthDt2;  // dH_Pphi/dsvN
  D_eq_(idxOffsetC  ,idxOffsetP+4) = (-k.x()*flightlengthY + k.y()*flightlengthX)/flightlengthDt2;  // dH_Pphi/dsvK
  d_eq_(idxOffsetC  ) = deltaPhi(flightlength.phi(), tauP4.phi());
  if ( verbosity_ >= 3 )
  {
    TVectorD Dana(Np_);
    TVectorD Dnum(Np_);
    for ( unsigned int idxPar = 0; idxPar < Np_; ++idxPar )
    {
      Dana(idxPar) = D_eq_(idxOffsetC,idxPar);
      Dnum(idxPar) = get_Dnum_H_Pphi(idxPar, idxOffsetC, idxOffsetP, tauP4, tauD3_rnk, r, n, k, rotMatrix_rnk2xyz, nuP4, sign);
    }
    std::cout << "H_Pphi:\n";
    std::cout << "flightlength: phi = " << flightlength.phi() << "\n";
    std::cout << "tauP4: phi = " << tauP4.phi() << "\n";
    printVector("analytic derrivative", Dana);
    printVector("numeric derrivative", Dnum);
  }
  D_eq_(idxOffsetC+1,0)            =  -flightlengthX*flightlengthZ/(flightlengthDt*flightlengthD2); // dH_Ptheta/dpvX
  D_eq_(idxOffsetC+1,1)            =  -flightlengthY*flightlengthZ/(flightlengthDt*flightlengthD2); // dH_Ptheta/dpvY
  D_eq_(idxOffsetC+1,2)            =   flightlengthDt/flightlengthD2;                               // dH_Ptheta/dpvZ
  D_eq_(idxOffsetC+1,idxOffsetP)   = (-tauPx*tauPz + tauPt2*nu_dPzdPx)/(tauPt*tauP2);               // dH_Ptheta/dnuPx
  D_eq_(idxOffsetC+1,idxOffsetP+1) = (-tauPy*tauPz + tauPt2*nu_dPzdPy)/(tauPt*tauP2);               // dH_Ptheta/dnuPy
  D_eq_(idxOffsetC+1,idxOffsetP+2) = (-r.z()*flightlengthDt2 + (r.x()*flightlengthX 
    + r.y()*flightlengthY)*flightlengthZ)/(flightlengthDt*flightlengthD2);                          // dH_Ptheta/dsvR
  D_eq_(idxOffsetC+1,idxOffsetP+3) = (-n.z()*flightlengthDt2 + (n.x()*flightlengthX 
    + n.y()*flightlengthY)*flightlengthZ)/(flightlengthDt*flightlengthD2);                          // dH_Ptheta/dsvN
  D_eq_(idxOffsetC+1,idxOffsetP+4) = (-k.z()*flightlengthDt2 + (k.x()*flightlengthX 
    + k.y()*flightlengthY)*flightlengthZ)/(flightlengthDt*flightlengthD2);                          // dH_Ptheta/dsvK
  d_eq_(idxOffsetC+1) = flightlength.theta() - tauP4.theta();
  if ( verbosity_ >= 3 )
  {
    TVectorD Dana(Np_);
    TVectorD Dnum(Np_);
    for ( unsigned int idxPar = 0; idxPar < Np_; ++idxPar )
    {
      Dana(idxPar) = D_eq_(idxOffsetC+1,idxPar);
      Dnum(idxPar) = get_Dnum_H_Ptheta(idxPar, idxOffsetC+1, idxOffsetP, tauP4, tauD3_rnk, r, n, k, rotMatrix_rnk2xyz, nuP4, sign);
    }
    std::cout << "H_Ptheta:\n";
    std::cout << "flightlength: theta = " << flightlength.theta() << "\n";
    std::cout << "tauP4: theta = " << tauP4.theta() << "\n";
    printVector("analytic derrivative", Dana);
    printVector("numeric derrivative", Dnum);
  }
}

double
KinFitConstraint::get_Dnum_H_recoilPx(unsigned int idxPar)
{
  double epsilon = 1.e-2;
  double H_nom = d_eq_(4);
  double H_up  = H_nom;
  if ( idxPar ==  3 || idxPar == 8 ) H_up += epsilon;
  if ( idxPar == 13                ) H_up -= epsilon;
  return (H_up - H_nom)/epsilon;
}

double
KinFitConstraint::get_Dnum_H_recoilPy(unsigned int idxPar)
{
  double epsilon = 1.e-2;
  double H_nom = d_eq_(5);
  double H_up  = H_nom;
  if ( idxPar ==  4 || idxPar == 9 ) H_up += epsilon;
  if ( idxPar == 14                ) H_up -= epsilon;
  return (H_up - H_nom)/epsilon;
}

double
KinFitConstraint::get_Dnum_H_recoilPz(unsigned int idxPar)
{
  double epsilon = 1.e-2;
  double H_nom = d_eq_(6);
  double H_up  = H_nom;
  double dummy;
  bool errorFlag;
  if ( idxPar == 3 || idxPar == 4 )
  {
    double nuTauPlusPx = nuTauPlusP4_.px();
    if ( idxPar == 3 ) nuTauPlusPx += epsilon;
    double nuTauPlusPy = nuTauPlusP4_.py();
    if ( idxPar == 4 ) nuTauPlusPy += epsilon;
    double nuTauPlusPz = comp_nuPz(visTauPlusP4_, nuTauPlusPx, nuTauPlusPy, signTauPlus_, dummy, dummy, errorFlag);
    H_up += (nuTauPlusPz - nuTauPlusP4_.pz());
  }
  else if ( idxPar == 8 || idxPar == 9 )
  {
    double nuTauMinusPx = nuTauMinusP4_.px();
    if ( idxPar == 8 ) nuTauMinusPx += epsilon;
    double nuTauMinusPy = nuTauMinusP4_.py();
    if ( idxPar == 9 ) nuTauMinusPy += epsilon;
    double nuTauMinusPz = comp_nuPz(visTauMinusP4_, nuTauMinusPx, nuTauMinusPy, signTauMinus_, dummy, dummy, errorFlag);
    H_up += (nuTauMinusPz - nuTauMinusP4_.pz());
  }
  else if ( idxPar == 15 )
  {
    H_up -= epsilon;
  }
  return (H_up - H_nom)/epsilon;
}

double
KinFitConstraint::get_Dnum_H_recoilEn(unsigned int idxPar)
{
  double epsilon = 1.e-2;
  double H_nom = d_eq_(7);
  double H_up  = H_nom;
  double dummy;
  bool errorFlag;
  if ( idxPar == 3 || idxPar == 4 )
  {
    double nuTauPlusPx = nuTauPlusP4_.px();
    if ( idxPar == 3 ) nuTauPlusPx += epsilon;
    double nuTauPlusPy = nuTauPlusP4_.py();
    if ( idxPar == 4 ) nuTauPlusPy += epsilon;
    double nuTauPlusPz = comp_nuPz(visTauPlusP4_, nuTauPlusPx, nuTauPlusPy, signTauPlus_, dummy, dummy, errorFlag);
    reco::Candidate::LorentzVector nuTauPlusP4 = build_nuP4(nuTauPlusPx, nuTauPlusPy, nuTauPlusPz);
    H_up += (nuTauPlusP4.energy() - nuTauPlusP4_.energy());
  }
  else if ( idxPar == 8 || idxPar == 9 )
  {
    double nuTauMinusPx = nuTauMinusP4_.px();
    if ( idxPar == 8 ) nuTauMinusPx += epsilon;
    double nuTauMinusPy = nuTauMinusP4_.py();
    if ( idxPar == 9 ) nuTauMinusPy += epsilon;
    double nuTauMinusPz = comp_nuPz(visTauMinusP4_, nuTauMinusPx, nuTauMinusPy, signTauMinus_, dummy, dummy, errorFlag);
    reco::Candidate::LorentzVector nuTauMinusP4 = build_nuP4(nuTauMinusPx, nuTauMinusPy, nuTauMinusPz);
    H_up += (nuTauMinusP4.energy() - nuTauMinusP4_.energy());
  }
  else if ( idxPar == 16 )
  {
    H_up -= epsilon;
  }
  return (H_up - H_nom)/epsilon;
}

void
KinFitConstraint::set_recoilConstraint()
{
  if ( verbosity_ >= 3 )
  {
    std::cout << "<set_recoilConstraint>:\n";
  }

  // CV: add constraint that recoil = tau+ + tau-
  //       H_px     = nuTauPlusPx + visTauPlusPx + nuTauMinusPx + visTauMinusPx - recoilPx = 0
  //       H_py     = nuTauPlusPy + visTauPlusPy + nuTauMinusPy + visTauMinusPy - recoilPy = 0
  //       H_pz     = nuTauPlusPz + visTauPlusPz + nuTauMinusPz + visTauMinusPz - recoilPz = 0
  //       H_energy = nuTauPlusE  + visTauPlusE  + nuTauMinusE  + visTauMinusE  - recoilE  = 0
  //     where:
  //       nuTauPlusE  = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2 )
  //       nuTauMinusE = sqrt(nuTauMisusPx^2 + nuTauMinusPy^2 + nuTauMinusPz^2)
  D_eq_(4, 3) = +1.;                                                                                // dH_recoilPx/dnuTauPlusPx
  D_eq_(4, 8) = +1.;                                                                                // dH_recoilPx/dnuTauMinusPx
  D_eq_(4,13) = -1.;                                                                                // dH_recoilPx/drecoilPx
  d_eq_(4) = tauPlusP4_.px() + tauMinusP4_.px() - recoilP4_.px();
  if ( verbosity_ >= 3 )
  {
    TVectorD Dana(Np_);
    TVectorD Dnum(Np_);
    for ( unsigned int idxPar = 0; idxPar < Np_; ++idxPar )
    {
      Dana(idxPar) = D_eq_(4,idxPar);
      Dnum(idxPar) = get_Dnum_H_recoilPx(idxPar);
    }
    std::cout << "H_recoilPx:\n";
    printVector("analytic derrivative", Dana);
    printVector("numeric derrivative", Dnum);
  }
  D_eq_(5, 4) = +1.;                                                                                // dH_recoilPy/dnuTauPlusPy
  D_eq_(5, 9) = +1.;                                                                                // dH_recoilPy/dnuTauMinusPy
  D_eq_(5,14) = -1.;                                                                                // dH_recoilPy/drecoilPy
  d_eq_(5) = tauPlusP4_.py() + tauMinusP4_.py() - recoilP4_.py();
  if ( verbosity_ >= 3 )
  {
    TVectorD Dana(Np_);
    TVectorD Dnum(Np_);
    for ( unsigned int idxPar = 0; idxPar < Np_; ++idxPar )
    {
      Dana(idxPar) = D_eq_(5,idxPar);
      Dnum(idxPar) = get_Dnum_H_recoilPy(idxPar);
    }
    std::cout << "H_recoilPy:\n";
    printVector("analytic derrivative", Dana);
    printVector("numeric derrivative", Dnum);
  }
  D_eq_(6, 3) = nuTauPlus_dPzdPx_;                                                                  // dH_recoilPz/dnuTauPlusPx
  D_eq_(6, 4) = nuTauPlus_dPzdPy_;                                                                  // dH_recoilPz/dnuTauPlusPy
  D_eq_(6, 8) = nuTauMinus_dPzdPx_;                                                                 // dH_recoilPz/dnuTauMinusPx
  D_eq_(6, 9) = nuTauMinus_dPzdPy_;                                                                 // dH_recoilPz/dnuTauMinusPy
  D_eq_(6,15) = -1.;                                                                                // dH_recoilPz/drecoilPz
  d_eq_(6) = tauPlusP4_.pz() + tauMinusP4_.pz() - recoilP4_.pz();
  if ( verbosity_ >= 3 )
  {
    TVectorD Dana(Np_);
    TVectorD Dnum(Np_);
    for ( unsigned int idxPar = 0; idxPar < Np_; ++idxPar )
    {
      Dana(idxPar) = D_eq_(6,idxPar);
      Dnum(idxPar) = get_Dnum_H_recoilPz(idxPar);
    }
    std::cout << "H_recoilPz:\n";
    printVector("analytic derrivative", Dana);
    printVector("numeric derrivative", Dnum);
  }
  D_eq_(7, 3) = (tauPlusP4_.px()  + nuTauPlus_dPzdPx_*tauPlusP4_.pz())/tauPlusP4_.energy();         // dH_recoilEn/dnuTauPlusPx
  D_eq_(7, 4) = (tauPlusP4_.py()  + nuTauPlus_dPzdPy_*tauPlusP4_.pz())/tauPlusP4_.energy();         // dH_recoilEn/dnuTauPlusPy
  D_eq_(7, 8) = (tauMinusP4_.px() + nuTauMinus_dPzdPx_*tauMinusP4_.pz())/tauMinusP4_.energy();      // dH_recoilEn/dnuTauPlusPx
  D_eq_(7, 9) = (tauMinusP4_.py() + nuTauMinus_dPzdPy_*tauMinusP4_.pz())/tauMinusP4_.energy();      // dH_recoilEn/dnuTauPlusPy
  D_eq_(7,16) = -1.;                                                                                // dH_recoilEn/drecoilE
  d_eq_(7) = tauPlusP4_.energy() + tauMinusP4_.energy() - recoilP4_.energy();
  if ( verbosity_ >= 3 )
  {
    TVectorD Dana(Np_);
    TVectorD Dnum(Np_);
    for ( unsigned int idxPar = 0; idxPar < Np_; ++idxPar )
    {
      Dana(idxPar) = D_eq_(7,idxPar);
      Dnum(idxPar) = get_Dnum_H_recoilEn(idxPar);
    }
    std::cout << "H_recoilEn:\n";
    printVector("analytic derrivative", Dana);
    printVector("numeric derrivative", Dnum);
  }
}

double
KinFitConstraint::get_Dnum_H_mHiggs(unsigned int idxPar)
{
  double epsilon = 1.e-2;
  double H_nom = d_eq_(8);
  double H_up  = H_nom;
  double dummy;
  bool errorFlag;
  if ( idxPar == 3 || idxPar == 4 )
  {
    double nuTauPlusPx = nuTauPlusP4_.px();
    if ( idxPar == 3 ) nuTauPlusPx += epsilon;
    double nuTauPlusPy = nuTauPlusP4_.py();
    if ( idxPar == 4 ) nuTauPlusPy += epsilon;
    double nuTauPlusPz = comp_nuPz(visTauPlusP4_, nuTauPlusPx, nuTauPlusPy, signTauPlus_, dummy, dummy, errorFlag);
    reco::Candidate::LorentzVector nuTauPlusP4 = build_nuP4(nuTauPlusPx, nuTauPlusPy, nuTauPlusPz);
    H_up = (visTauPlusP4_ + nuTauPlusP4 + tauMinusP4_).mass();
  }
  else if ( idxPar == 8 || idxPar == 9 )
  {
    double nuTauMinusPx = nuTauMinusP4_.px();
    if ( idxPar == 8 ) nuTauMinusPx += epsilon;
    double nuTauMinusPy = nuTauMinusP4_.py();
    if ( idxPar == 9 ) nuTauMinusPy += epsilon;
    double nuTauMinusPz = comp_nuPz(visTauMinusP4_, nuTauMinusPx, nuTauMinusPy, signTauMinus_, dummy, dummy, errorFlag);
    reco::Candidate::LorentzVector nuTauMinusP4 = build_nuP4(nuTauMinusPx, nuTauMinusPy, nuTauMinusPz);
    H_up = (tauPlusP4_ + visTauMinusP4_ + nuTauMinusP4).mass(); 
  }
  return (H_up - H_nom)/epsilon;
}
  
void
KinFitConstraint::set_higgsMassConstraint()
{
  if ( verbosity_ >= 3 )
  {
    std::cout << "<set_higgsMassConstraint>:\n";
  }

  double tauPlusPx  = tauPlusP4_.px();
  double tauPlusPy  = tauPlusP4_.py();
  double tauPlusPz  = tauPlusP4_.pz();
  double tauPlusE   = tauPlusP4_.energy();

  double tauMinusPx = tauMinusP4_.px();
  double tauMinusPy = tauMinusP4_.py();
  double tauMinusPz = tauMinusP4_.pz();
  double tauMinusE  = tauMinusP4_.energy();

  // CV: add Higgs mass constraint
  //       H = sqrt(2*mTau^2 
  //              + 2*(visTauPlusE  + nuTauPlusE )*(visTauMinusE  + nuTauMinusE ) 
  //              - 2*(visTauPlusPx + nuTauPlusPx)*(visTauMinusPx + nuTauMinusPx)
  //              - 2*(visTauPlusPy + nuTauPlusPy)*(visTauMinusPy + nuTauMinusPy)
  //              - 2*(visTauPlusPz + nuTauPlusPz)*(visTauMinusPz + nuTauMinusPz)) - mHiggs = 0
  //     where
  //       nuTauPlusE  = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2 )
  //       nuTauMinusE = sqrt(nuTauMisusPx^2 + nuTauMinusPy^2 + nuTauMinusPz^2)
  double denominator0 = (tauPlusP4_ + tauMinusP4_).mass();
  D_eq_(8,3) = ((tauMinusE/tauPlusE)*(tauPlusPx  + nuTauPlus_dPzdPx_*tauPlusPz)   
    - tauMinusPx - nuTauPlus_dPzdPx_*tauMinusPz)/denominator0;                                   // dH_mHiggs/dnuTauPlusPx
  D_eq_(8,4) = ((tauMinusE/tauPlusE)*(tauPlusPy  + nuTauPlus_dPzdPy_*tauPlusPz)   
    - tauMinusPy - nuTauPlus_dPzdPy_*tauMinusPz)/denominator0;                                   // dH_mHiggs/dnuTauPlusPy
  D_eq_(8,8) = ((tauPlusE/tauMinusE)*(tauMinusPx + nuTauMinus_dPzdPx_*tauMinusPz) 
    - tauPlusPx  - nuTauMinus_dPzdPx_*tauPlusPz)/denominator0;                                   // dH_mHiggs/dnuTauMinusPx
  D_eq_(8,9) = ((tauPlusE/tauMinusE)*(tauMinusPy + nuTauMinus_dPzdPy_*tauMinusPz) 
    - tauPlusPy  - nuTauMinus_dPzdPy_*tauPlusPz)/denominator0;                                   // dH_mHiggs/dnuTauMinusPy
  d_eq_(8) = (tauPlusP4_ + tauMinusP4_).mass() - mHiggs;
  if ( verbosity_ >= 3 )
  {
    TVectorD Dana(Np_);
    TVectorD Dnum(Np_);
    for ( unsigned int idxPar = 0; idxPar < Np_; ++idxPar )
    {
      Dana(idxPar) = D_eq_(8,idxPar);
      Dnum(idxPar) = get_Dnum_H_mHiggs(idxPar);
    }
    std::cout << "H_mHiggs:\n";
    printVector("analytic derrivative", Dana);
    printVector("numeric derrivative", Dnum);
  }
}

double
KinFitConstraint::get_Dnum_H_mT(unsigned int idxPar, unsigned int idxOffsetC, unsigned int idxOffsetP,
                                const reco::Candidate::LorentzVector& visP4,
                                const reco::Candidate::LorentzVector& nuP4)
{
  double epsilon = 1.e-2;
  double H_nom = d_ineq_(idxOffsetC);
  double H_up = H_nom;
  double nuPx = nuP4.px();
  if ( idxPar == 3 || idxPar == 8 ) nuPx += epsilon;
  double nuPy = nuP4.py();
  if ( idxPar == 4 || idxPar == 9 ) nuPy += epsilon;
  H_up = comp_mT(visP4.mass(), visP4.px(), visP4.py(), 0., nuPx, nuPy) - mTau;
  return (H_up - H_nom)/epsilon;
}

void
KinFitConstraint::set_mTConstraint(const std::string& label, unsigned int idxOffsetC, unsigned int idxOffsetP,
                                   const reco::Candidate::LorentzVector& visP4,
                                   const reco::Candidate::LorentzVector& nuP4)
{
  if ( verbosity_ >= 3 )
  {
    std::cout << "<set_mTConstraint (" << label << ")>:\n";
  }

  double visPx   = visP4.px();
  double visPy   = visP4.py();
  double visMass = visP4.mass();
  double visEt   = std::sqrt(square(visMass) + square(visPx) + square(visPy));

  double nuPx    = nuP4.px();
  double nuPy    = nuP4.py();
  double nuEt    = std::sqrt(square(nuPx) + square(nuPy));

  double mT      = comp_mT(visMass, visPx, visPy, 0., nuPx, nuPy);

  // CV: add transverse mass (mT) constraint for tau+ and tau-
  //     The purpose of this constraint is to ensure that computation of neutrino Pz yields a physical solution
  //     Note that this constraint is an inequiality constraint: mT(visP4,nuP4) <= mTau
  D_ineq_(idxOffsetC  ,idxOffsetP)   =  -(1./mT)*(visPx - (visEt/nuEt)*nuPx);                       // dH_mT/dnuPx
  D_ineq_(idxOffsetC  ,idxOffsetP+1) =  -(1./mT)*(visPy - (visEt/nuEt)*nuPy);                       // dH_mT/dnuPy
  d_ineq_(idxOffsetC  ) = mT - mTau;
  if ( verbosity_ >= 3 )
  {
    TVectorD Dana(Np_);
    TVectorD Dnum(Np_);
    for ( unsigned int idxPar = 0; idxPar < Np_; ++idxPar )
    {
      Dana(idxPar) = D_ineq_(idxOffsetC,idxPar);
      Dnum(idxPar) = get_Dnum_H_mT(idxPar, idxOffsetC, idxOffsetP, visP4, nuP4);
    }
    std::cout << "H_mT:\n";
    printVector("analytic derrivative", Dana);
    printVector("numeric derrivative", Dnum);
  }
}

