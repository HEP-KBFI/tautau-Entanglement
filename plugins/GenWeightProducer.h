#ifndef TauAnalysis_Entanglement_GenWeightProducer_h
#define TauAnalysis_Entanglement_GenWeightProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" // GenEventInfoProduct

class GenWeightProducer : public edm::EDProducer 
{
 public:
  explicit GenWeightProducer(const edm::ParameterSet& cfg);
  ~GenWeightProducer();

 private:
  void produce(edm::Event& evt, const edm::EventSetup& es);

  edm::InputTag src_;
  edm::EDGetTokenT<GenEventInfoProduct> token_;
};

#endif
