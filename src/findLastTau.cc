#include "TauAnalysis/Entanglement/interface/findLastTau.h"

const reco::GenParticle*
findLastTau(const reco::GenParticle* mother)
{
  size_t numDaughters = mother->numberOfDaughters();
  for ( size_t idxDaughter = 0; idxDaughter < numDaughters; ++idxDaughter )
  {
    const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(mother->daughter(idxDaughter));
    assert(daughter);
    if ( daughter->pdgId() == mother->pdgId() )
    {
      return findLastTau(daughter);
    }
  }
  return mother;
}
