#include "TauAnalysis/Entanglement/interface/compVisP4.h"

reco::Candidate::LorentzVector
compVisP4(const std::vector<const reco::GenParticle*>& decayProducts)
{
  reco::Candidate::LorentzVector visP4;
  for ( const reco::GenParticle* decayProduct : decayProducts )
  {
    int absPdgId = std::abs(decayProduct->pdgId());
    if ( absPdgId == 12 || absPdgId == 14 || absPdgId == 16) continue;
    visP4 += decayProduct->p4();
  }
  return visP4;
}
