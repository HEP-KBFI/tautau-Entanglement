#include "TauAnalysis/Entanglement/plugins/HepMCGenParticleProducer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::*

#include "HepMC/IO_GenEvent.h" // HepMC::*

namespace
{
std::unique_ptr<std::vector<reco::GenParticle>> convertHepMCtoGenParticle(const HepMC::GenEvent* hepmc_event)
  {
    auto genParticles = std::make_unique<std::vector<reco::GenParticle>>();

    for(HepMC::GenEvent::particle_const_iterator p = hepmc_event->particles_begin(); p != hepmc_event->particles_end(); ++p)
    {
      const HepMC::GenParticle* hepmc_particle = *p;

      int status = hepmc_particle->status();
      int pdgId = hepmc_particle->pdg_id();
      const HepMC::FourVector& p4 = hepmc_particle->momentum();

      reco::Candidate::Point vertex;
      if (hepmc_particle->production_vertex())
      {
        vertex.SetXYZ(hepmc_particle->production_vertex()->point3d().x(),
                      hepmc_particle->production_vertex()->point3d().y(),
                      hepmc_particle->production_vertex()->point3d().z());
      }

      reco::GenParticle genParticle(status, reco::Candidate::LorentzVector(p4.px(), p4.py(), p4.pz(), p4.e()), vertex, pdgId, hepmc_particle->barcode(), true);

      genParticles->push_back(genParticle);
    }

    return genParticles;
  }
}

HepMCGenParticleProducer::HepMCGenParticleProducer(const edm::ParameterSet& params)
    : fileName_(params.getParameter<std::string>("inputFile"))
    , asciiInput_(new HepMC::IO_GenEvent(fileName_.c_str(), std::ios::in))
{
  produces<std::vector<reco::GenParticle>>("genParticles");
}

HepMCGenParticleProducer::~HepMCGenParticleProducer()
{
  delete asciiInput_;
}

void HepMCGenParticleProducer::produce(edm::Event& evt, const edm::EventSetup&)
{
  HepMC::GenEvent* evt_hepmc = new HepMC::GenEvent();

  if(! asciiInput_->fill_next_event(evt_hepmc))
  {
    delete evt_hepmc;
    return;
  }
  
  auto genParticles = ::convertHepMCtoGenParticle(evt_hepmc);
  evt.put(std::move(genParticles));
}

#include "FWCore/Framework/interface/MakerMacros.h" // DEFINE_FWK_MODULE

DEFINE_FWK_MODULE(HepMCGenParticleProducer);
