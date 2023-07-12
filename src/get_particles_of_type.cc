#include "TauAnalysis/Entanglement/interface/get_particles_of_type.h"

#include <cmath> // std::abs(), std::fabs()

std::vector<const reco::GenParticle*>
get_chargedHadrons(const std::vector<const reco::GenParticle*>& decayProducts)
{
  std::vector<const reco::GenParticle*> chargedHadrons;
  for ( const reco::GenParticle* decayProduct : decayProducts )
  {
    int absPdgId = std::abs(decayProduct->pdgId());
    if ( absPdgId == 11 || absPdgId == 13 ) continue;
    if ( std::fabs(decayProduct->charge()) > 0.5 )
    {
      chargedHadrons.push_back(decayProduct);
    }
  }
  return chargedHadrons;
}

std::vector<const reco::GenParticle*>
get_chargedKaons(const std::vector<const reco::GenParticle*>& decayProducts)
{
  return get_particles_of_type(decayProducts, { 321 });
}

std::vector<const reco::GenParticle*>
get_neutralKaons(const std::vector<const reco::GenParticle*>& decayProducts)
{
  return get_particles_of_type(decayProducts, { 130, 310, 311 });
}

std::vector<const reco::GenParticle*>
get_neutralPions(const std::vector<const reco::GenParticle*>& decayProducts)
{
  return get_particles_of_type(decayProducts, { 111 });
}

std::vector<const reco::GenParticle*>
get_neutrinos(const std::vector<const reco::GenParticle*>& decayProducts)
{
  return get_particles_of_type(decayProducts, { 12, 14, 16 });
}

std::vector<const reco::GenParticle*>
get_particles_of_type(const std::vector<const reco::GenParticle*>& decayProducts, const std::vector<int>& selPdgIds)
{
  std::vector<const reco::GenParticle*> selDecayProducts;
  for ( const reco::GenParticle* decayProduct : decayProducts )
  {
    int absPdgId = std::abs(decayProduct->pdgId());
    bool isSelPdgId = false;
    for ( int selPdgId : selPdgIds )
    {
      if ( absPdgId == selPdgId )
      {
        isSelPdgId = true;
        break;
      }
    }
    if ( isSelPdgId )
    {
      selDecayProducts.push_back(decayProduct);
    }
  }
  return selDecayProducts;
}

std::vector<const reco::GenParticle*>
get_photons(const std::vector<const reco::GenParticle*>& decayProducts)
{
  return get_particles_of_type(decayProducts, { 22 });
}
