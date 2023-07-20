#include "TauAnalysis/Entanglement/interface/get_leadTrack.h"

const reco::GenParticle*
get_leadTrack(const std::vector<const reco::GenParticle*>& decayProducts)
{
  const reco::GenParticle* leadTrack = nullptr;
  double max_pt = -1.;
  for ( const reco::GenParticle* decayProduct : decayProducts )
  {
    if ( std::fabs(decayProduct->charge()) > 0.5 && decayProduct->pt() > max_pt )
    {
      leadTrack = decayProduct;
      max_pt = decayProduct->pt();
    }
  }
  return leadTrack;
}

const KinematicParticle*
get_leadTrack(const std::vector<KinematicParticle>& daughters)
{
  const KinematicParticle* leadTrack = nullptr;
  double max_pt = -1.;
  for ( const KinematicParticle& daughter : daughters )
  {
    if ( std::fabs(daughter.charge()) > 0.5 && daughter.p4().pt() > max_pt )
    {
      leadTrack = &daughter;
      max_pt = daughter.p4().pt();
    }
  }
  return leadTrack;
}
