#ifndef TauAnalysis_Entanglement_KinFitConstraint_h
#define TauAnalysis_Entanglement_KinFitConstraint_h

#include "DataFormats/Candidate/interface/Candidate.h"               // reco::Candidate::LorentzVector, reco::Candidate::Vector

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"       // KinematicEvent
#include "TauAnalysis/Entanglement/interface/KinFitConstraintBase.h" // KinFitConstraintBase

#include <TMatrixD.h>                                                // TMatrixD
#include <TVectorD.h>                                                // TVectorD

class KinFitConstraint : public KinFitConstraintBase
{
 public:
  KinFitConstraint(int collider, const KinematicEvent& kineEvt, double signTauPlus, double signTauMinus, int verbosity = -1);
  ~KinFitConstraint();

  void
  set_alphaA(const TVectorD& alphaA);

 protected:
  reco::Candidate::Vector
  get_flightlength_up(unsigned int idxPar, unsigned int idxOffsetP, double epsilon,
                      const reco::Candidate::Vector& tauD3_rnk,
                      const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k, 
                      const math::Matrix3x3& rotMatrix_rnk2xyz);

  reco::Candidate::LorentzVector
  get_tauP4_up(unsigned int idxPar, unsigned int idxOffsetP, double epsilon,
               const reco::Candidate::LorentzVector& tauP4,
               const reco::Candidate::LorentzVector& nuP4, double sign);

  double
  get_Dnum_H_Pphi(unsigned int idxPar, unsigned int idxOffsetC, unsigned int idxOffsetP, 
                  const reco::Candidate::LorentzVector& tauP4,
                  const reco::Candidate::Vector& tauD3_rnk,
                  const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k, 
                  const math::Matrix3x3& rotMatrix_rnk2xyz,
                  const reco::Candidate::LorentzVector& nuP4, double sign);

  double
  get_Dnum_H_Ptheta(unsigned int idxPar, unsigned int idxOffsetC, unsigned int idxOffsetP, 
                    const reco::Candidate::LorentzVector& tauP4,
                    const reco::Candidate::Vector& tauD3_rnk,
                    const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k, 
                    const math::Matrix3x3& rotMatrix_rnk2xyz,
                    const reco::Candidate::LorentzVector& nuP4, double sign);

  void
  set_parlConstraint(const std::string& label, unsigned int idxOffsetC, unsigned int idxOffsetP,
                     const reco::Candidate::LorentzVector& tauP4, 
                     const reco::Candidate::Vector& tauD3_rnk,
                     const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k, 
                     const math::Matrix3x3& rotMatrix_rnk2xyz,
                     const reco::Candidate::LorentzVector& nuP4, double nu_dPzdPx, double nu_dPzdPy,
                     double sign);

  double
  get_Dnum_H_recoilPx(unsigned int idxPar);

  double
  get_Dnum_H_recoilPy(unsigned int idxPar);

  double
  get_Dnum_H_recoilPz(unsigned int idxPar);

  double
  get_Dnum_H_recoilEn(unsigned int idxPar);

  void
  set_recoilConstraint();

  double
  get_Dnum_H_mHiggs(unsigned int idxPar);
  
  void
  set_higgsMassConstraint();

  double
  get_Dnum_H_mT(unsigned int idxPar, unsigned int idxOffsetC, unsigned int idxOffsetP,
                const reco::Candidate::LorentzVector& visP4,
                const reco::Candidate::LorentzVector& nuP4);

  void
  set_mTConstraint(const std::string& label, unsigned int idxOffsetC, unsigned int idxOffsetP,
                   const reco::Candidate::LorentzVector& visP4,
                   const reco::Candidate::LorentzVector& nuP4);

  double
  get_Dnum_H_flightlength_k(unsigned int idxPar, unsigned int idxOffsetC, unsigned int idxOffsetP,
                            const reco::Candidate::Vector& tauD3_rnk,
                            const reco::Candidate::Vector& k,
                            const math::Matrix3x3& rotMatrix_rnk2xyz);

  void
  set_flightlength_kConstraint(const std::string& label, unsigned int idxOffsetC, unsigned int idxOffsetP,
                               const reco::Candidate::Vector& tauD3_rnk,
                               const reco::Candidate::Vector& k,
                               const math::Matrix3x3& rotMatrix_rnk2xyz);

  int collider_;

  const KinematicEvent& kineEvt_;

  double signTauPlus_;
  double signTauMinus_;

  reco::Candidate::Point pv0_;                   // original (unfitted) position of primary event vertex
  reco::Candidate::Point pv_;                    // fitted position of primary event vertex

  reco::Candidate::LorentzVector tauPlusP4_;
  reco::Candidate::Vector tauPlusD3_rnk_;        // flightlength of tau+ = sv(tau+) - pv, in helicity-frame coordinates { r, n, k }
  reco::Candidate::LorentzVector visTauPlusP4_;
  reco::Candidate::LorentzVector nuTauPlusP4_;
  double nuTauPlus_dPzdPx_;
  double nuTauPlus_dPzdPy_;

  reco::Candidate::LorentzVector tauMinusP4_;
  reco::Candidate::Vector tauMinusD3_rnk_;       // flightlength of tau- = sv(tau-) - pv, in helicity-frame coordinates { r, n, k }
  reco::Candidate::LorentzVector visTauMinusP4_;
  reco::Candidate::LorentzVector nuTauMinusP4_;
  double nuTauMinus_dPzdPx_;
  double nuTauMinus_dPzdPy_;

  reco::Candidate::LorentzVector recoilP4_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraint_h
