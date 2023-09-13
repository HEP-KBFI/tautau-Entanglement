#ifndef TauAnalysis_Entanglement_KinFitConstraint_h
#define TauAnalysis_Entanglement_KinFitConstraint_h

#include "DataFormats/Math/interface/Matrix.h"                       // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                       // math::Vector

#include "TauAnalysis/Entanglement/interface/comp_nuPz.h"            // build_nuP4(), comp_nuPz()
#include "TauAnalysis/Entanglement/interface/constants.h"            // kLHC, kSuperKEKB, mHiggs
#include "TauAnalysis/Entanglement/interface/KinFitConstraintBase.h" // KinFitConstraintBase<>
#include "TauAnalysis/Entanglement/interface/square.h"               // square()

// CV: The template parameters P and C represent 
//     the number of parameters and the number of constraints, respectively.
//     It is neccessary to make the number of parameters and constraints template parameters,
//     because the ROOT::Math::SVector and ROOT::Math::SMatrix classes
//     on which the kinematic fitting code is based require these template parameters.

template <unsigned int P, unsigned int C>
class KinFitConstraint : public KinFitConstraintBase<P,C>
{
 public:
  KinFitConstraint(int collider, const KinematicEvent& kineEvt, double signTauPlus, double signTauMinus, int verbosity = -1)
    : KinFitConstraintBase<P,C>(verbosity)
    , collider_(collider)
    , kineEvt_(kineEvt)
    , signTauPlus_(signTauPlus)
    , signTauMinus_(signTauMinus)
  {
    d_metric_(  0,  0) = 1.e+2;                                    // H_Pphi(tau+)
    d_metric_(  1,  1) = 1.e+2;                                    // H_Ptheta(tau+)
    d_metric_(  2,  2) = 1.e+2;                                    // H_Pphi(tau-)
    d_metric_(  3,  3) = 1.e+2;                                    // H_Ptheta(tau-)
    d_metric_(  4,  4) = ( collider_ == kSuperKEKB ) ? 1. : 1.e-2; // H_px(recoil)
    d_metric_(  5,  5) = ( collider_ == kSuperKEKB ) ? 1. : 1.e-2; // H_py(recoil)
    d_metric_(  6,  6) = ( collider_ == kSuperKEKB ) ? 1. : 1.e-2; // H_pz(recoil)
    d_metric_(  7,  7) = ( collider_ == kSuperKEKB ) ? 1. : 1.e-2; // H_energy(recoil)
    d_metric_(C-1,C-1) = 1.e-4;                                    // H(Higgs mass)
  }
  ~KinFitConstraint()
  {}

  void
  set_alphaA(const typename KinFitConstraintBase<P,C>::VectorP& alphaA)
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

    // CV: Reconstruct event kinematics
    reco::Candidate::Point pv = reco::Candidate::Point(alphaA(0), alphaA(1), alphaA(2));

    visTauPlusP4_ = kineEvt_.visTauPlusP4();
    double nuTauPlusPx = alphaA(3);
    double nuTauPlusPy = alphaA(4);
    bool nuTauPlus_errorFlag = false;
    double nuTauPlusPz = comp_nuPz(visTauPlusP4_, nuTauPlusPx, nuTauPlusPy, signTauPlus_, nuTauPlus_dPzdPx_, nuTauPlus_dPzdPy_, nuTauPlus_errorFlag, verbosity_);
    nuTauPlusP4_ = build_nuP4(nuTauPlusPx, nuTauPlusPy, nuTauPlusPz);
    tauPlusP4_ = visTauPlusP4_ + nuTauPlusP4_;
    reco::Candidate::Point svTauPlus = reco::Candidate::Point(alphaA(5), alphaA(6), alphaA(7));
    tauPlusD3_ = svTauPlus - pv;

    visTauMinusP4_ = kineEvt_.visTauMinusP4();
    double nuTauMinusPx = alphaA(8);
    double nuTauMinusPy = alphaA(9);
    bool nuTauMinus_errorFlag = false;
    double nuTauMinusPz = comp_nuPz(visTauMinusP4_, nuTauMinusPx, nuTauMinusPy, signTauMinus_, nuTauMinus_dPzdPx_, nuTauMinus_dPzdPy_, nuTauMinus_errorFlag, verbosity_);
    nuTauMinusP4_ = build_nuP4(nuTauMinusPx, nuTauMinusPy, nuTauMinusPz);
    tauMinusP4_ = visTauMinusP4_ + nuTauMinusP4_;
    reco::Candidate::Point svTauMinus = reco::Candidate::Point(alphaA(10), alphaA(11), alphaA(12));
    tauMinusD3_ = svTauMinus - pv;

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
    set_parlConstraint(0, 3, tauPlusP4_,  tauPlusD3_,  nuTauPlus_dPzdPx_,  nuTauPlus_dPzdPy_);
    set_parlConstraint(2, 8, tauMinusP4_, tauMinusD3_, nuTauMinus_dPzdPx_, nuTauMinus_dPzdPy_);

    // CV: add constraint that recoil = tau+ + tau-
    set_recoilConstraint();

    // CV: add Higgs mass constraint (not for SuperKEKB!)
    if ( collider_ != kSuperKEKB )
    {
      set_higgsMassConstraint();
    }

    errorFlag_ = nuTauPlus_errorFlag || nuTauMinus_errorFlag;    
    if ( verbosity_ >= 2 )
    {
      std::cout << "errorFlag = " << errorFlag_ << "\n";
    }
  }

 protected:
  void
  set_parlConstraint(unsigned int idxOffsetC, unsigned int idxOffsetP,
                     const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::Vector& tauD3,
                     double nu_dPzdPx, double nu_dPzdPy)
  {
    double tauPx = tauP4.px();
    double tauPy = tauP4.py();
    double tauPz = tauP4.pz();
    double tauPt = tauP4.pt();
    double tauP  = tauP4.P();

    double tauDx = tauD3.X();
    double tauDy = tauD3.Y();
    double tauDz = tauD3.Z();
    double tauDt = std::sqrt(tauD3.Perp2());
    double tauD  = std::sqrt(tauD3.Mag2());

    // CV: add "parallelism" constraint for tau+ and tau-
    //     The constraint equations H_Pphi and H_Ptheta are given by Eqs. (4.10) and (4.11)
    //     in Section 4.1.3.3 of [1].
    //     The small correction for the helix bend between tau lepton production and decay in Eq. (4.13) is not used, 
    //     as this correction is expected to be on the level of 50 microrad, which we consider negligible.
    //     Following [1], we perform a transformation of variables from phi and theta to 1/2*pi - phi and 1/2*pi - theta,
    //     in order to reduce the magnitude of derivatives (by avoiding the division by small numbers)
    //
    //   [1] https://cds.cern.ch/record/1358627/files/CERN-THESIS-2011-028.pdf
    D_(idxOffsetC  ,0)            =  (1. - tauDx/tauDt)/tauDy;                                                                                // dH_Pphi/dpvX
    D_(idxOffsetC  ,1)            = -(tauDx/tauDy)*D_(idxOffsetC,0);                                                                          // dH_Pphi/dpvY
    D_(idxOffsetC  ,idxOffsetP)   =  (1. - tauPx/tauPt)/tauPy;                                                                                // dH_Pphi/dnuPx
    D_(idxOffsetC  ,idxOffsetP+1) = -(tauPx/tauPy)*D_(idxOffsetC,idxOffsetP);                                                                 // dH_Pphi/dnuPy     
    D_(idxOffsetC  ,idxOffsetP+2) = -D_(idxOffsetC,0);                                                                                        // dH_Pphi/dsvX
    D_(idxOffsetC  ,idxOffsetP+3) = -D_(idxOffsetC,1);                                                                                        // dH_Pphi/dsvY
    d_(idxOffsetC  ) = (tauDt - tauDx)/tauDy + (tauPx - tauPt)/tauPy;
    D_(idxOffsetC+1,0)            =  (tauDx/tauDz)*(1./tauDt - 1./tauD);                                                                      // dH_Ptheta/dpvX
    D_(idxOffsetC+1,1)            =  (tauDy/tauDz)*(1./tauDt - 1./tauD);                                                                      // dH_Ptheta/dpvY
    D_(idxOffsetC+1,2)            = -1./tauD + (tauD - tauDt)/square(tauDz);                                                                  // dH_Ptheta/dpvZ
    D_(idxOffsetC+1,idxOffsetP)   =  (1./square(tauPz))*nu_dPzdPx*(tauP - tauPt) + (1./tauPz)*(tauPx/tauPt - (tauPx + nu_dPzdPx*tauPz)/tauP); // dH_Ptheta/dnuPx
    D_(idxOffsetC+1,idxOffsetP+1) =  (1./square(tauPz))*nu_dPzdPy*(tauP - tauPt) + (1./tauPz)*(tauPy/tauPt - (tauPy + nu_dPzdPy*tauPz)/tauP); // dH_Ptheta/dnuPy
    D_(idxOffsetC+1,idxOffsetP+2) = -D_(idxOffsetC+1,0);                                                                                      // dH_Ptheta/dsvX
    D_(idxOffsetC+1,idxOffsetP+3) = -D_(idxOffsetC+1,1);                                                                                      // dH_Ptheta/dsvY
    D_(idxOffsetC+1,idxOffsetP+4) = -D_(idxOffsetC+1,2);                                                                                      // dH_Ptheta/dsvZ
    d_(idxOffsetC+1) = (tauD - tauDt)/tauDz + (tauPt - tauP)/tauPz;
  }

  void
  set_recoilConstraint()
  {
    // CV: add constraint that recoil = tau+ + tau-
    //       H_px     = nuTauPlusPx + visTauPlusPx + nuTauMinusPx + visTauMinusPx - recoilPx = 0
    //       H_py     = nuTauPlusPy + visTauPlusPy + nuTauMinusPy + visTauMinusPy - recoilPy = 0
    //       H_pz     = nuTauPlusPz + visTauPlusPz + nuTauMinusPz + visTauMinusPz - recoilPz = 0
    //       H_energy = nuTauPlusE  + visTauPlusE  + nuTauMinusE  + visTauMinusE  - recoilE  = 0
    //     where:
    //       nuTauPlusE  = sqrt(nuTauPlusPx^2  + nuTauPlusPy^2  + nuTauPlusPz^2 )
    //       nuTauMinusE = sqrt(nuTauMisusPx^2 + nuTauMinusPy^2 + nuTauMinusPz^2)
    D_(4, 3) = +1.;                                                                           // dH_px/dnuTauPlusPx
    D_(4, 8) = +1.;                                                                           // dH_px/dnuTauMinusPx
    D_(4,13) = -1.;                                                                           // dH_px/drecoilPx
    d_(4) = tauPlusP4_.px() + tauMinusP4_.px() - recoilP4_.px();
    D_(5, 4) = +1.;                                                                           // dH_py/dnuTauPlusPy
    D_(5, 9) = +1.;                                                                           // dH_py/dnuTauMinusPy
    D_(5,14) = -1.;                                                                           // dH_py/drecoilPy
    d_(5) = tauPlusP4_.py() + tauMinusP4_.py() - recoilP4_.py();
    D_(6, 3) = nuTauPlus_dPzdPx_;                                                             // dH_pz/dnuTauPlusPx
    D_(6, 4) = nuTauPlus_dPzdPy_;                                                             // dH_pz/dnuTauPlusPy
    D_(6, 8) = nuTauMinus_dPzdPx_;                                                            // dH_pz/dnuTauMinusPx
    D_(6, 9) = nuTauMinus_dPzdPy_;                                                            // dH_pz/dnuTauMinusPy
    D_(6,15) = -1.;                                                                           // dH_pz/drecoilPz
    d_(6) = tauPlusP4_.pz() + tauMinusP4_.pz() - recoilP4_.pz();
    D_(7, 3) = (tauPlusP4_.px()  + nuTauPlus_dPzdPx_*tauPlusP4_.pz())/tauPlusP4_.energy();    // dH_energy/dnuTauPlusPx
    D_(7, 4) = (tauPlusP4_.py()  + nuTauPlus_dPzdPy_*tauPlusP4_.pz())/tauPlusP4_.energy();    // dH_energy/dnuTauPlusPy
    D_(7, 8) = (tauMinusP4_.px() + nuTauMinus_dPzdPx_*tauMinusP4_.pz())/tauMinusP4_.energy(); // dH_energy/dnuTauPlusPx
    D_(7, 9) = (tauMinusP4_.py() + nuTauMinus_dPzdPy_*tauMinusP4_.pz())/tauMinusP4_.energy(); // dH_energy/dnuTauPlusPy
    D_(7,16) = -1.;                                                                           // dH_energy/drecoilE
    d_(7) = tauPlusP4_.energy() + tauMinusP4_.energy() - recoilP4_.energy();
  }
  
  void
  set_higgsMassConstraint()
  {
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
    // CV: Setting idx to 8 generates compile-time errors that array indices are out-of-bounds
    //     for ROOT::Math::SMatrix and ROOT::Math::SVector objects.
    //     The error disappears when setting idx to numConstraints - 1 instead.
    //     Note that the index will be wrong when collider = SuperKEKB,
    //     but this does not matter as the code will actually be executed for collider = LHC only.
    const int idxC = C - 1;
    D_(idxC,3) = ((tauMinusE/tauPlusE)*(tauPlusPx  + nuTauPlus_dPzdPx_*tauPlusPz)   - tauMinusPx - nuTauPlus_dPzdPx_*tauMinusPz)/denominator0; // dH/dnuTauPlusPx
    D_(idxC,4) = ((tauMinusE/tauPlusE)*(tauPlusPy  + nuTauPlus_dPzdPy_*tauPlusPz)   - tauMinusPy - nuTauPlus_dPzdPy_*tauMinusPz)/denominator0; // dH/dnuTauPlusPy
    D_(idxC,8) = ((tauPlusE/tauMinusE)*(tauMinusPx + nuTauMinus_dPzdPx_*tauMinusPz) - tauPlusPx  - nuTauMinus_dPzdPx_*tauPlusPz)/denominator0; // dH/dnuTauMinusPx
    D_(idxC,9) = ((tauPlusE/tauMinusE)*(tauMinusPy + nuTauMinus_dPzdPy_*tauMinusPz) - tauPlusPy  - nuTauMinus_dPzdPy_*tauPlusPz)/denominator0; // dH/dnuTauMinusPy
    d_(idxC) = (tauPlusP4_ + tauMinusP4_).mass() - mHiggs;
  }

  using KinFitConstraintBase<P,C>::D_;
  using KinFitConstraintBase<P,C>::d_;
  using KinFitConstraintBase<P,C>::d_metric_;
  using KinFitConstraintBase<P,C>::errorFlag_;
  using KinFitConstraintBase<P,C>::verbosity_;

  int collider_;

  const KinematicEvent& kineEvt_;

  double signTauPlus_;
  double signTauMinus_;

  reco::Candidate::LorentzVector tauPlusP4_;
  reco::Candidate::Vector tauPlusD3_;
  reco::Candidate::LorentzVector visTauPlusP4_;
  reco::Candidate::LorentzVector nuTauPlusP4_;
  double nuTauPlus_dPzdPx_;
  double nuTauPlus_dPzdPy_;

  reco::Candidate::LorentzVector tauMinusP4_;
  reco::Candidate::Vector tauMinusD3_;
  reco::Candidate::LorentzVector visTauMinusP4_;
  reco::Candidate::LorentzVector nuTauMinusP4_;
  double nuTauMinus_dPzdPx_;
  double nuTauMinus_dPzdPy_;

  reco::Candidate::LorentzVector recoilP4_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraintBase_h
