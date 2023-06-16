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
  KinematicEvent(const reco::Candidate::Point& pv, const reco::Candidate::LorentzVector& recoilP4);
  ~KinematicEvent();

  void
  set_tauPlusP4(const reco::Candidate::LorentzVector& tauPlusP4);

  void
  set_visTauPlus(const reco::Candidate::LorentzVector& visTauPlusP4, int tauPlus_decaymode, 
                 const std::vector<KinematicParticle>& daughtersTauPlus,
                 const reco::Candidate::Point& tipPCATauPlus);

  void
  set_svTauPlus(const reco::Candidate::Point& svTauPlus);

  void
  set_tauMinusP4(const reco::Candidate::LorentzVector& tauMinusP4);

  void
  set_visTauMinus(const reco::Candidate::LorentzVector& visTauMinusP4, int tauMinus_decaymode, 
                  const std::vector<KinematicParticle>& daughtersTauMinus,
                  const reco::Candidate::Point& tipPCATauMinus);

  void
  set_svTauMinus(const reco::Candidate::Point& svTauMinus);

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
  get_tauPlus_decaymode() const;

  const std::vector<KinematicParticle>&
  get_daughtersTauPlus() const;

  const reco::Candidate::Point&
  get_tipPCATauPlus() const;

  const reco::Candidate::Point&
  get_svTauPlus() const;

  bool
  get_svTauPlus_isValid() const;

  reco::Candidate::LorentzVector
  get_tauMinusP4() const;

  bool
  get_tauMinusP4_isValid() const;

  const reco::Candidate::LorentzVector&
  get_visTauMinusP4() const;

  int
  get_tauMinus_decaymode() const;

  const std::vector<KinematicParticle>&
  get_daughtersTauMinus() const;

  const reco::Candidate::Point&
  get_tipPCATauMinus() const;

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
  reco::Candidate::LorentzVector visTauPlusP4_;
  bool visTauPlus_isValid_;
  int tauPlus_decaymode_;
  std::vector<KinematicParticle> daughtersTauPlus_;
  reco::Candidate::Point tipPCATauPlus_;
  reco::Candidate::Point svTauPlus_;
  bool svTauPlus_isValid_;

  reco::Candidate::LorentzVector tauMinusP4_;
  bool tauMinusP4_isValid_;
  reco::Candidate::LorentzVector visTauMinusP4_;
  bool visTauMinus_isValid_;
  int tauMinus_decaymode_;
  std::vector<KinematicParticle> daughtersTauMinus_;
  reco::Candidate::Point tipPCATauMinus_;
  reco::Candidate::Point svTauMinus_;
  bool svTauMinus_isValid_;
};

void
printKinematicEvent(const std::string& label,
                    const KinematicEvent& evt,
                    bool cartesian = true);

#endif // TauAnalysis_Entanglement_KinematicEvent_h
