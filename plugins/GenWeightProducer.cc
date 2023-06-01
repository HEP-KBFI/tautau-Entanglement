#include "TauAnalysis/Entanglement/plugins/GenWeightProducer.h"

#include "FWCore/Utilities/interface/InputTag.h"

GenWeightProducer::GenWeightProducer(const edm::ParameterSet& cfg) 
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<GenEventInfoProduct>(src_);

  produces<double>();
}

GenWeightProducer::~GenWeightProducer()
{}

void GenWeightProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<GenEventInfoProduct> genEventInfoProduct;
  evt.getByLabel(src_, genEventInfoProduct);

  std::unique_ptr<double> genWeight(new double(genEventInfoProduct->weight()));
  evt.put(std::move(genWeight));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenWeightProducer);
