#ifndef TauAnalysis_Entanglement_KinFitConstraintAnalyzer_h
#define TauAnalysis_Entanglement_KinFitConstraintAnalyzer_h

#include "FWCore/Framework/interface/one/EDAnalyzer.h"                   // edm::one::EDAnalyzer<>
#include "FWCore/Framework/interface/Event.h"                            // edm::Event
#include "FWCore/Framework/interface/EventSetup.h"                       // edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"                  // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h"                         // edm::InputTag<>

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"            // reco::GenParticle
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"         // reco::GenParticleCollection

#include "TauAnalysis/Entanglement/interface/EntanglementNtuple.h"       // EntanglementNtuple
#include "TauAnalysis/Entanglement/interface/GenKinematicEventBuilder.h" // GenKinematicEventBuilder

#include <vector>                                                        // std::vector<>
#include <string>                                                        // std::string

class KinFitConstraintAnalyzer : public edm::one::EDAnalyzer<>
{
 public:
  explicit KinFitConstraintAnalyzer(const edm::ParameterSet&);
  ~KinFitConstraintAnalyzer();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::GenParticleCollection> token_;

  GenKinematicEventBuilder* genKineEvtBuilder_;

  bool applySmearing_;

  int collider_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraintAnalyzer_h
