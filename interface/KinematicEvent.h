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

class KinematicEvent
{
 public:
  KinematicEvent();
  ~KinematicEvent();

  const reco::Candidate::Point&
  get_pv() const;

  const reco::Candidate::LorentzVector&
  get_recoilP4() const;

  reco::Candidate::LorentzVector
  get_tauPlusP4() const;

  bool
  get_tauPlusP4_isValid() const;

  const reco::Candidate::LorentzVector&
  get_visTauPlusP4() const;

  int
  get_tauPlus_decayMode() const;

  const std::vector<KinematicParticle>&
  get_daughtersTauPlus() const;

  const reco::Candidate::Point&
  get_tipPCATauPlus() const;

  const reco::Candidate::Point&
  get_svTauPlus() const;

  bool
  get_svTauPlus_isValid() const;

  const reco::Candidate::Vector&
  get_hPlus() const;

  bool
  get_hPlus_isValid() const;

  reco::Candidate::LorentzVector
  get_tauMinusP4() const;

  bool
  get_tauMinusP4_isValid() const;

  const reco::Candidate::LorentzVector&
  get_visTauMinusP4() const;

  int
  get_tauMinus_decayMode() const;

  const std::vector<KinematicParticle>&
  get_daughtersTauMinus() const;

  const reco::Candidate::Point&
  get_tipPCATauMinus() const;

  const reco::Candidate::Point&
  get_svTauMinus() const;

  bool
  get_svTauMinus_isValid() const;

  const reco::Candidate::Vector&
  get_hMinus() const;

  bool
  get_hMinus_isValid() const;

  friend class GenKinematicEventBuilder;
  friend class KinematicFitStartPosFinder;
  friend class Smearing;

 private:
  reco::Candidate::Point pv_;

  reco::Candidate::LorentzVector recoilP4_;

  reco::Candidate::LorentzVector tauPlusP4_;
  bool tauPlusP4_isValid_;
  reco::Candidate::LorentzVector visTauPlusP4_;
  int tauPlus_decayMode_;
  std::vector<KinematicParticle> daughtersTauPlus_;
  reco::Candidate::Point tipPCATauPlus_;
  reco::Candidate::Point svTauPlus_;
  bool svTauPlus_isValid_;
  reco::Candidate::Vector hPlus_;
  bool hPlus_isValid_;

  reco::Candidate::LorentzVector tauMinusP4_;
  bool tauMinusP4_isValid_;
  reco::Candidate::LorentzVector visTauMinusP4_;
  int tauMinus_decayMode_;
  std::vector<KinematicParticle> daughtersTauMinus_;
  reco::Candidate::Point tipPCATauMinus_;
  reco::Candidate::Point svTauMinus_;
  bool svTauMinus_isValid_;
  reco::Candidate::Vector hMinus_;
  bool hMinus_isValid_;
};

void
printKinematicEvent(const std::string& label,
                    const KinematicEvent& kineEvt,
                    bool cartesian = true);

#endif // TauAnalysis_Entanglement_KinematicEvent_h
