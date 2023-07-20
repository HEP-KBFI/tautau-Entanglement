#ifndef TauAnalysis_Entanglement_KinematicEvent_h
#define TauAnalysis_Entanglement_KinematicEvent_h

/** KinematicEvent
 *
 * Representation of an event, used as input and output data format in kinematic fits.
 *
 * \authors Christian Veelken, Tallinn
 *
 */

#include "DataFormats/Candidate/interface/Candidate.h"            // reco::Candidate::LorentzVector, reco::Candidate::Point

#include "TauAnalysis/Entanglement/interface/KinematicParticle.h" // KinematicParticle
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3, math::Matrix4x4

class KinematicEvent
{
 public:
  KinematicEvent();
  ~KinematicEvent();

  const reco::Candidate::Point&
  pv() const;

  const math::Matrix3x3&
  pvCov() const;

  const reco::Candidate::LorentzVector&
  recoilP4() const;

  const math::Matrix4x4&
  recoilCov() const;

  const reco::Candidate::LorentzVector&
  tauPlusP4() const;

  bool
  tauPlusP4_isValid() const;

  const math::Matrix4x4&
  tauPlusCov() const;

  bool
  tauPlusCov_isValid() const;

  const reco::Candidate::LorentzVector&
  visTauPlusP4() const;

  const math::Matrix4x4&
  visTauPlusCov() const;

  int
  tauPlus_decayMode() const;

  const std::vector<KinematicParticle>&
  daughtersTauPlus() const;

  const reco::Candidate::Point&
  tipPCATauPlus() const;

  const reco::Candidate::Point&
  svTauPlus() const;

  const math::Matrix3x3&
  svTauPlusCov() const;

  bool
  svTauPlus_isValid() const;

  const reco::Candidate::Vector&
  hPlus() const;

  bool
  hPlus_isValid() const;

  const reco::Candidate::LorentzVector&
  tauMinusP4() const;

  bool
  tauMinusP4_isValid() const;

  const math::Matrix4x4&
  tauMinusCov() const;

  bool
  tauMinusCov_isValid() const;

  const reco::Candidate::LorentzVector&
  visTauMinusP4() const;

  const math::Matrix4x4&
  visTauMinusCov() const;

  int
  tauMinus_decayMode() const;

  const std::vector<KinematicParticle>&
  daughtersTauMinus() const;

  const reco::Candidate::Point&
  tipPCATauMinus() const;

  const reco::Candidate::Point&
  svTauMinus() const;

  const math::Matrix3x3&
  svTauMinusCov() const;

  bool
  svTauMinus_isValid() const;

  const reco::Candidate::Vector&
  hMinus() const;

  bool
  hMinus_isValid() const;

  double
  kinFitChi2() const;

  const math::MatrixPxP&
  kinFitCov() const;

  bool
  kinFit_isValid() const;

  friend class GenKinematicEventBuilder;
  friend class KinematicFit;
  friend class KinematicFitStartPosFinder;
  friend class Smearing;

 private:
  reco::Candidate::Point pv_;
  math::Matrix3x3 pvCov_;

  reco::Candidate::LorentzVector recoilP4_;
  math::Matrix4x4 recoilCov_;

  reco::Candidate::LorentzVector tauPlusP4_;
  bool tauPlusP4_isValid_;
  math::Matrix4x4 tauPlusCov_;
  bool tauPlusCov_isValid_;
  reco::Candidate::LorentzVector visTauPlusP4_;
  math::Matrix4x4 visTauPlusCov_;
  int tauPlus_decayMode_;
  std::vector<KinematicParticle> daughtersTauPlus_;
  reco::Candidate::Point tipPCATauPlus_;
  reco::Candidate::Point svTauPlus_;
  math::Matrix3x3 svTauPlusCov_;
  bool svTauPlus_isValid_;
  reco::Candidate::Vector hPlus_;
  bool hPlus_isValid_;

  reco::Candidate::LorentzVector tauMinusP4_;
  bool tauMinusP4_isValid_;
  math::Matrix4x4 tauMinusCov_;
  bool tauMinusCov_isValid_;
  reco::Candidate::LorentzVector visTauMinusP4_;
  math::Matrix4x4 visTauMinusCov_;
  int tauMinus_decayMode_;
  std::vector<KinematicParticle> daughtersTauMinus_;
  reco::Candidate::Point tipPCATauMinus_;
  reco::Candidate::Point svTauMinus_;
  math::Matrix3x3 svTauMinusCov_;
  bool svTauMinus_isValid_;
  reco::Candidate::Vector hMinus_;
  bool hMinus_isValid_;

  math::MatrixPxP kinFitCov_;
  double kinFitChi2_;
  bool kinFit_isValid_;
};

void
printKinematicEvent(const std::string& label,
                    const KinematicEvent& kineEvt,
                    bool cartesian = true);

#endif // TauAnalysis_Entanglement_KinematicEvent_h
