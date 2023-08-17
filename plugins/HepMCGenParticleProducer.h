#ifndef TauAnalysis_Entanglement_HepMCGenParticleProducer_h
#define TauAnalysis_Entanglement_HepMCGenParticleProducer_h

#include "FWCore/Framework/interface/EDProducer.h" // edm::EDProducer
#include "FWCore/Framework/interface/Event.h" // edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet

// forward declarations
namespace HepMC {
  class IO_GenEvent;
}

class HepMCGenParticleProducer
  : public edm::EDProducer
{
public:
  explicit HepMCGenParticleProducer(const edm::ParameterSet&);
  ~HepMCGenParticleProducer() override;
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  std::string fileName_;
  HepMC::IO_GenEvent* asciiInput_;
};

#endif // TauAnalysis_Entanglement_HepMCGenParticleProducer_h
