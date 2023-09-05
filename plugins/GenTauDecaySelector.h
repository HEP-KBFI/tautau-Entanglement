#ifndef TauAnalysis_Entanglement_GenTauDecaySelector_h
#define TauAnalysis_Entanglement_GenTauDecaySelector_h

#include "FWCore/Framework/interface/stream/EDProducer.h"            // edm::stream::EDProducer
#include "FWCore/Framework/interface/Event.h"                        // edm::Event
#include "FWCore/Framework/interface/EventSetup.h"                   // edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"              // edm::ParameterSet
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h" // edm::ConfigurationDescriptions
#include "FWCore/Utilities/interface/InputTag.h"                     // edm::InputTag

#include "DataFormats/JetReco/interface/GenJet.h"                    // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h"          // reco::GenJetCollection

class GenTauDecaySelector : public edm::stream::EDProducer<>
{
 public:
  explicit GenTauDecaySelector(const edm::ParameterSet& cfg);
  ~GenTauDecaySelector();

  static 
  void 
  fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  void 
  produce(edm::Event& evt, const edm::EventSetup& es);

  edm::InputTag src_;
  edm::EDGetTokenT<reco::GenJetCollection> token_;

  int maxNumChargedKaons_;
  int maxNumNeutralKaons_;
  int maxNumPhotons_;
  double maxSumPhotonEn_;
};

#endif // TauAnalysis_Entanglement_GenTauDecaySelector_h
