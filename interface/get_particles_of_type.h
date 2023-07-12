#ifndef TauAnalysis_Entanglement_get_particles_of_type_h
#define TauAnalysis_Entanglement_get_particles_of_type_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle

#include <vector>                                             // std::vector<>

std::vector<const reco::GenParticle*>
get_chargedHadrons(const std::vector<const reco::GenParticle*>& decayProducts);

std::vector<const reco::GenParticle*>
get_chargedKaons(const std::vector<const reco::GenParticle*>& decayProducts);

std::vector<const reco::GenParticle*>
get_neutralKaons(const std::vector<const reco::GenParticle*>& decayProducts);

std::vector<const reco::GenParticle*>
get_neutralPions(const std::vector<const reco::GenParticle*>& decayProducts);

std::vector<const reco::GenParticle*>
get_neutrinos(const std::vector<const reco::GenParticle*>& decayProducts);

std::vector<const reco::GenParticle*>
get_particles_of_type(const std::vector<const reco::GenParticle*>& decayProducts, const std::vector<int>& selPdgIds);

std::vector<const reco::GenParticle*>
get_photons(const std::vector<const reco::GenParticle*>& decayProducts);

#endif // TauAnalysis_Entanglement_get_particles_of_type_h
