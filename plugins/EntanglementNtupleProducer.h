#ifndef TauAnalysis_Entanglement_EntanglementNtupleProducer_h
#define TauAnalysis_Entanglement_EntanglementNtupleProducer_h

#include "FWCore/Framework/interface/one/EDAnalyzer.h"                   // edm::one::EDAnalyzer<>
#include "FWCore/Framework/interface/Event.h"                            // edm::Event
#include "FWCore/Framework/interface/EventSetup.h"                       // edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"                  // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h"                         // edm::InputTag<>

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"            // reco::GenParticle
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"         // reco::GenParticleCollection

#include "TauAnalysis/Entanglement/interface/EntanglementNtuple.h"       // EntanglementNtuple
#include "TauAnalysis/Entanglement/interface/GenKinematicEventBuilder.h" // GenKinematicEventBuilder
#include "TauAnalysis/Entanglement/interface/KinematicFit.h"             // KinematicFit
#include "TauAnalysis/Entanglement/interface/StartPosFinder.h"           // StartPosFinder

#include <TTree.h>                                                       // TTree

#include <vector>                                                        // std::vector<>
#include <string>                                                        // std::string

class EntanglementNtupleProducer : public edm::one::EDAnalyzer<>
{
 public:
  explicit EntanglementNtupleProducer(const edm::ParameterSet&);
  ~EntanglementNtupleProducer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::GenParticleCollection> token_;

  GenKinematicEventBuilder* genKineEvtBuilder_woSmearing_;
  GenKinematicEventBuilder* genKineEvtBuilder_wSmearing_;

  StartPosFinder* startPosFinder_;

  KinematicFit* kinematicFit_;

  int collider_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights_;
  std::vector<edm::EDGetTokenT<double>> tokenWeights_;

  TTree* ntuple_piPlus_piMinus_;
  EntanglementNtuple* ntupleFiller_piPlus_piMinus_;

  TTree* ntuple_piPlus_rhoMinus_;
  EntanglementNtuple* ntupleFiller_piPlus_rhoMinus_;

  TTree* ntuple_rhoPlus_piMinus_;
  EntanglementNtuple* ntupleFiller_rhoPlus_piMinus_;

  TTree* ntuple_rhoPlus_rhoMinus_;
  EntanglementNtuple* ntupleFiller_rhoPlus_rhoMinus_;

  TTree* ntuple_had_had_;
  EntanglementNtuple* ntupleFiller_had_had_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_EntanglementNtupleProducer_h
