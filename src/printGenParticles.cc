#include "TauAnalysis/Entanglement/interface/printGenParticles.h"

void
printGenParticles(const std::vector<const reco::GenParticle*>& decayProducts,
                  bool cartesian)
{
  size_t numDecayProducts = decayProducts.size();
  for ( size_t idxDecayProduct = 0; idxDecayProduct < numDecayProducts; ++idxDecayProduct )
  {
    const reco::GenParticle* decayProduct = decayProducts[idxDecayProduct];
    std::cout << "GenParticle #" << idxDecayProduct << ":";
    reco::Candidate::LorentzVector p4 = decayProduct->p4();
    if ( cartesian )
    {
      std::cout << " E = " << p4.energy() << ", Px = " << p4.px() << ", Py = " << p4.py() << ", Pz = " << p4.pz();
    }
    else
    {
      std::cout << " pT = " << p4.pt() << ", eta = " << p4.eta() << ", phi = " << p4.phi() << ", mass = " << p4.mass();
    }
    std::cout << ", pdgId = " << decayProduct->pdgId() << "\n";
  } 
}
