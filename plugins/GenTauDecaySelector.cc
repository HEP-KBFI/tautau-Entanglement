#include "TauAnalysis/Entanglement/plugins/GenTauDecaySelector.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"             // edm::ConsumesCollector

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"         // reco::GenParticle
#pragma GCC diagnostic pop

#include "TauAnalysis/Entanglement/interface/comp_visP4.h"            // comp_visP4()
#include "TauAnalysis/Entanglement/interface/get_particles_of_type.h" // get_chargedHadrons(), get_neutralPions(), get_neutrinos()

#include <memory>                                                     // std::unique_ptr

GenTauDecaySelector::GenTauDecaySelector(const edm::ParameterSet& cfg)
  : maxNumChargedKaons_(cfg.getParameter<int>("maxNumChargedKaons"))
  , maxNumNeutralKaons_(cfg.getParameter<int>("maxNumNeutralKaons"))
  , maxNumPhotons_(cfg.getParameter<int>("maxNumPhotons"))
  , maxSumPhotonEn_(cfg.getParameter<double>("maxSumPhotonEn"))
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenJetCollection>(src_);

  produces<reco::GenJetCollection>();
}

GenTauDecaySelector::~GenTauDecaySelector()
{}

void
GenTauDecaySelector::produce(edm::Event& evt, const edm::EventSetup& es)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<GenTauDecaySelector::produce>:" << "\n";
  }

  edm::Handle<reco::GenJetCollection> genTaus;
  evt.getByToken(token_, genTaus);
  if ( verbosity_ >= 1 )
  {
    std::cout << "#genTaus = " << genTaus->size() << "\n";
  }

  std::unique_ptr<reco::GenJetCollection> selectedGenTaus(new reco::GenJetCollection());

  size_t numGenTaus = genTaus->size();
  for ( size_t idxGenTau = 0; idxGenTau < numGenTaus; ++idxGenTau )
  {
    const reco::GenJet& genTau = genTaus->at(idxGenTau);
    if ( verbosity_ >= 2 )
    {
      std::cout << "genTau #" << idxGenTau << ": pT = " << genTau.pt() << ", theta = " << genTau.theta() << ", phi = " << genTau.phi() << "\n";
    }

    std::vector<const reco::GenParticle*> daughters = genTau.getGenConstituents();
    int numChargedKaons = get_chargedKaons(daughters).size();
    int numNeutralKaons = get_neutralKaons(daughters).size();
    std::vector<const reco::GenParticle*> photons = get_photons(daughters);
    int numPhotons = photons.size();
    double sumPhotonEn = comp_visP4(photons).energy();
    if ( verbosity_ >= 2 )
    {
      std::cout << "numChargedKaons = " << numChargedKaons << " (maxNumChargedKaons = " << maxNumChargedKaons_ << ")\n";
      std::cout << "numNeutralKaons = " << numNeutralKaons << " (maxNumNeutralKaons = " << maxNumNeutralKaons_ << ")\n";
      std::cout << "numPhotons = " << numPhotons << " (maxNumPhotons = " << maxNumPhotons_ << ")\n";
      std::cout << "sumPhotonEn = " << sumPhotonEn << " (maxSumPhotonEn = " << maxSumPhotonEn_ << ")\n";
    }

    if ( (maxNumChargedKaons_ == -1  || numChargedKaons <= maxNumChargedKaons_) &&
         (maxNumNeutralKaons_ == -1  || numNeutralKaons <= maxNumNeutralKaons_) &&
         (maxNumPhotons_      == -1  || numPhotons      <= maxNumPhotons_     ) &&
         (maxSumPhotonEn_     <   0. || sumPhotonEn     <= maxSumPhotonEn_    ) ) {
      selectedGenTaus->push_back(genTau);
    }
  }
  if ( verbosity_ >= 1 )
  {
    std::cout << "#selectedGenTaus = " << selectedGenTaus->size() << "\n";
  }

  evt.put(std::move(selectedGenTaus));
}

void
GenTauDecaySelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag(""));
  desc.add<int>("maxNumChargedKaons", 0);
  desc.add<int>("maxNumNeutralKaons", 0);
  desc.add<int>("maxNumPhotons", -1);
  desc.add<double>("maxSumPhotonEn", 0.5);
  desc.addUntracked<int>("verbosity", 0);
  descriptions.add("genTauDecaySelector", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
 
DEFINE_FWK_MODULE(GenTauDecaySelector);
