#ifndef TauAnalysis_Entanglement_KinematicEvent_h
#define TauAnalysis_Entanglement_KinematicEvent_h

/** KinematicEvent
 *
 * Representation of an event, used as input and output data format in kinematic fits.
 *
 * \authors Christian Veelken, Tallinn
 *
 */

#include "DataFormats/Candidate/interface/Candidate.h"            // Candidate::LorentzVector, Candidate::Point

#include "TauAnalysis/Entanglement/interface/KinematicParticle.h" // KinematicParticle

class KinematicEvent
{
 public:
  KinematicEvent(const reco::Candidate::Point& pv,
                 const reco::Candidate::LorentzVector& recoilP4,
                 const std::vector<KinematicParticle>& daughtersTauPlus,
                 const std::vector<KinematicParticle>& daughtersTauMinus);
  ~KinematicEvent();

  void
  set_tauPlusP4(const reco::Candidate::LorentzVector& tauPlusP4);

  void
  set_svTauPlus(const reco::Candidate::Point& svTauPlus);

  void
  set_tauMinusP4(const reco::Candidate::LorentzVector& tauMinusP4);

  void
  set_svTauMinus(const reco::Candidate::Point& svTauMinus);

  const reco::Candidate::Point&
  get_pv() const;

  const reco::Candidate::LorentzVector&
  get_recoilP4() const;

  const std::vector<KinematicParticle>&
  get_daughtersTauPlus() const;

  const reco::Candidate::Point&
  get_svTauPlus() const;

  bool
  get_svTauPlus_isValid() const;

  const std::vector<KinematicParticle>&
  get_daughtersTauMinus() const;

  const reco::Candidate::Point&
  get_svTauMinus() const;

  bool
  get_svTauMinus_isValid() const;

  friend class Smearing;

 private:
  reco::Candidate::Point pv_;

  reco::Candidate::LorentzVector recoilP4_;

  reco::Candidate::LorentzVector tauPlusP4_;
  bool tauPlusP4_isValid_;
  std::vector<KinematicParticle> daughtersTauPlus_;
  reco::Candidate::Point svTauPlus_;
  bool svTauPlus_isValid_;

  reco::Candidate::LorentzVector tauMinusP4_;
  bool tauMinusP4_isValid_;
  std::vector<KinematicParticle> daughtersTauMinus_;
  reco::Candidate::Point svTauMinus_;
  bool svTauMinus_isValid_;
};

#endif // TauAnalysis_Entanglement_KinematicEvent_h
