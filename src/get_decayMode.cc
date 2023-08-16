#include "TauAnalysis/Entanglement/interface/get_decayMode.h"

#include "DataFormats/TauReco/interface/PFTau.h" // reco::PFTau::hadronicDecayMode

int
get_decayMode(const std::vector<const reco::GenParticle*>& tau_ch,
              const std::vector<const reco::GenParticle*>& tau_pi0,
              const std::vector<const reco::GenParticle*>& tau_nu)
{
  size_t nch = tau_ch.size();
  size_t npi0 = tau_pi0.size();
  size_t nnu = tau_nu.size();
  // CV: there are rare (on the level of 10^-4) events in the official CMS Monte Carlo samples
  //     in which hadronic tau decays involve more than one neutrino;
  //     let's for now simply discard those
  if      ( nch == 1 && npi0 == 0 && nnu == 1 ) return reco::PFTau::kOneProng0PiZero;
  else if ( nch == 1 && npi0 == 1 && nnu == 1 ) return reco::PFTau::kOneProng1PiZero;
  else if ( nch == 1 && npi0 == 2 && nnu == 1 ) return reco::PFTau::kOneProng2PiZero;
  else if ( nch == 3 && npi0 == 0 && nnu == 1 ) return reco::PFTau::kThreeProng0PiZero;
  else                                          return reco::PFTau::kNull;
}

bool
is1Prong(int decayMode)
{
  if ( decayMode == reco::PFTau::kOneProng0PiZero || 
       decayMode == reco::PFTau::kOneProng1PiZero || 
       decayMode == reco::PFTau::kOneProng2PiZero ) return true;
  else                                              return false;
}

bool
is3Prong(int decayMode)
{
  if ( decayMode == reco::PFTau::kThreeProng0PiZero ) return true;
  else                                                return false;
}
