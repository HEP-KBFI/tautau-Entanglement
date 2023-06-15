#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException

KinematicEvent::KinematicEvent(const reco::Candidate::Point& pv,
                               const reco::Candidate::LorentzVector& recoilP4,
                               const reco::Candidate::LorentzVector& visTauPlusP4, int tauPlus_decaymode, 
                               const std::vector<KinematicParticle>& daughtersTauPlus,
                               const reco::Candidate::Point& tipPCATauPlus,
                               const reco::Candidate::LorentzVector& visTauMinusP4, int tauMinus_decaymode, 
                               const std::vector<KinematicParticle>& daughtersTauMinus,
                               const reco::Candidate::Point& tipPCATauMinus)
  : pv_(pv)
  , recoilP4_(recoilP4)
  , tauPlusP4_isValid_(false)
  , visTauPlusP4_(visTauPlusP4)
  , tauPlus_decaymode_(tauPlus_decaymode)
  , daughtersTauPlus_(daughtersTauPlus)
  , tipPCATauPlus_(tipPCATauPlus)
  , svTauPlus_isValid_(false)
  , tauMinusP4_isValid_(false)
  , visTauMinusP4_(visTauMinusP4)
  , tauMinus_decaymode_(tauMinus_decaymode)
  , daughtersTauMinus_(daughtersTauMinus)
  , tipPCATauMinus_(tipPCATauMinus)
  , svTauMinus_isValid_(false)
{}

KinematicEvent::~KinematicEvent()
{}

void
KinematicEvent::set_tauPlusP4(const reco::Candidate::LorentzVector& tauPlusP4)
{
  tauPlusP4_ = tauPlusP4;
  tauPlusP4_isValid_ = true;
}

void
KinematicEvent::set_svTauPlus(const reco::Candidate::Point& svTauPlus)
{
  svTauPlus_ = svTauPlus;
  svTauPlus_isValid_ = true;
}

void
KinematicEvent::set_tauMinusP4(const reco::Candidate::LorentzVector& tauMinusP4)
{
  tauMinusP4_ = tauMinusP4;
  tauMinusP4_isValid_ = true;
}

void
KinematicEvent::set_svTauMinus(const reco::Candidate::Point& svTauMinus)
{
  svTauMinus_ = svTauMinus;
  svTauMinus_isValid_ = true;
}

const reco::Candidate::Point&
KinematicEvent::get_pv() const
{
  return pv_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::get_recoilP4() const
{
  return recoilP4_;
}

reco::Candidate::LorentzVector
KinematicEvent::get_tauPlusP4() const
{
  return tauPlusP4_;
}

bool
KinematicEvent::get_tauPlusP4_isValid() const
{
  return tauPlusP4_isValid_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::get_visTauPlusP4() const
{
  return visTauPlusP4_;
}

int
KinematicEvent::get_tauPlus_decaymode() const
{
  return tauPlus_decaymode_;
}

const std::vector<KinematicParticle>&
KinematicEvent::get_daughtersTauPlus() const
{
  return daughtersTauPlus_;
}

const reco::Candidate::Point&
KinematicEvent::get_tipPCATauPlus() const
{
  return tipPCATauPlus_;
}

const reco::Candidate::Point&
KinematicEvent::get_svTauPlus() const
{
  if ( !svTauPlus_isValid_ )
    throw cmsException("KinematicParticle", __LINE__)
      << "Tau+ decay vertex not initialized !!\n";
  return svTauPlus_;
}

bool
KinematicEvent::get_svTauPlus_isValid() const
{
  return svTauPlus_isValid_;
}

reco::Candidate::LorentzVector
KinematicEvent::get_tauMinusP4() const
{
  return tauMinusP4_;
}

bool
KinematicEvent::get_tauMinusP4_isValid() const
{
  return tauMinusP4_isValid_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::get_visTauMinusP4() const
{
  return visTauMinusP4_;
}

int
KinematicEvent::get_tauMinus_decaymode() const
{
  return tauMinus_decaymode_;
}

const std::vector<KinematicParticle>&
KinematicEvent::get_daughtersTauMinus() const
{
  return daughtersTauMinus_;
}

const reco::Candidate::Point&
KinematicEvent::get_tipPCATauMinus() const
{
  return tipPCATauMinus_;
}

const reco::Candidate::Point&
KinematicEvent::get_svTauMinus() const
{
  if ( !svTauMinus_isValid_ )
    throw cmsException("KinematicParticle", __LINE__)
      << "Tau- decay vertex not initialized !!\n";
  return svTauMinus_;
}

bool
KinematicEvent::get_svTauMinus_isValid() const
{
  return svTauMinus_isValid_;
}
