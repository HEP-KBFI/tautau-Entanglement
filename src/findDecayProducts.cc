#include "TauAnalysis/Entanglement/interface/findDecayProducts.h"

#include <algorithm> // std::sort()

namespace
{
  bool
  isHigherPt(const reco::GenParticle* particle1, const reco::GenParticle* particle2)
  {
    return particle1->pt() > particle2->pt();
  }
}

void
findDecayProducts(const reco::GenParticle* mother, 
                  std::vector<const reco::GenParticle*>& daughters)
{
  size_t numDaughters = mother->numberOfDaughters();
  for ( size_t idxDaughter = 0; idxDaughter < numDaughters; ++idxDaughter )
  {
    const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(mother->daughter(idxDaughter));
    assert(daughter);
    // CV: treat neutral pions as stable particles
    if ( daughter->status() == 1 || daughter->pdgId() == 111 )
    {
      daughters.push_back(daughter);
    }
    else 
    {
      findDecayProducts(daughter, daughters);
    }
  }
  // CV: sort tau decay products by decreasing pT
  std::sort(daughters.begin(), daughters.end(), isHigherPt);
}
