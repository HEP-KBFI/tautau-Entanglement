#include "TauAnalysis/Entanglement/interface/comp_nuP4.h"

reco::Candidate::LorentzVector
comp_nuP4(const std::vector<const reco::GenParticle*>& decayProducts)
{
  reco::Candidate::LorentzVector nuP4;
  for ( const reco::GenParticle* decayProduct : decayProducts )
  {
    int absPdgId = std::abs(decayProduct->pdgId());
    if ( absPdgId == 12 || absPdgId == 14 || absPdgId == 16 ) nuP4 += decayProduct->p4();
  }
  return nuP4;
}
