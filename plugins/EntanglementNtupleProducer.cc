#include "TauAnalysis/Entanglement/plugins/EntanglementNtupleProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"                 // edm::Service<>
#include "CommonTools/UtilAlgos/interface/TFileService.h"             // TFileService

#include "DataFormats/Common/interface/Handle.h"                      // edm::Handle<>

#include "DataFormats/TauReco/interface/PFTau.h"                      // reco::PFTau::hadronicDecayMode

#include "TauAnalysis/Entanglement/interface/compVisP4.h"             // compVisP4()
#include "TauAnalysis/Entanglement/interface/findDecayProducts.h"     // findDecayProducts()
#include "TauAnalysis/Entanglement/interface/findLastTau.h"           // findLastTau()
#include "TauAnalysis/Entanglement/interface/get_decayMode.h"         // get_decayMode()
#include "TauAnalysis/Entanglement/interface/get_particles_of_type.h" // get_chargedHadrons(), get_neutralPions(), get_neutrinos()

#include <iostream>                                                   // std::cout

EntanglementNtupleProducer::EntanglementNtupleProducer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , genKineEvtBuilder_woSmearing_(nullptr)
  , genKineEvtBuilder_wSmearing_(nullptr)
  , kineFitStartPosFinder_(cfg)
  , kineFit_(cfg)
  , ntuple_piPlus_piMinus_(nullptr)
  , ntupleFiller_piPlus_piMinus_(nullptr)
  , ntuple_piPlus_rhoMinus_(nullptr)
  , ntupleFiller_piPlus_rhoMinus_(nullptr)
  , ntuple_rhoPlus_piMinus_(nullptr)
  , ntupleFiller_rhoPlus_piMinus_(nullptr)
  , ntuple_rhoPlus_rhoMinus_(nullptr)
  , ntupleFiller_rhoPlus_rhoMinus_(nullptr)
  , ntuple_had_had_(nullptr)
  , ntupleFiller_had_had_(nullptr)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenParticleCollection>(src_);

  bool applySmearing = cfg.getParameter<bool>("applySmearing");

  edm::ParameterSet cfg_woSmearing = cfg;
  cfg_woSmearing.addParameter<bool>("applySmearing", false);
  genKineEvtBuilder_woSmearing_ = new GenKinematicEventBuilder(cfg_woSmearing);

  edm::ParameterSet cfg_wSmearing = cfg;
  cfg_wSmearing.addParameter<bool>("applySmearing", applySmearing);
  genKineEvtBuilder_wSmearing_ = new GenKinematicEventBuilder(cfg_wSmearing);

  srcWeights_ = cfg.getParameter<vInputTag>("srcEvtWeights");
  for ( const edm::InputTag& srcWeight : srcWeights_ )
  {
    tokenWeights_.push_back(consumes<double>(srcWeight));
  }
}

EntanglementNtupleProducer::~EntanglementNtupleProducer()
{
  delete genKineEvtBuilder_woSmearing_;
  delete genKineEvtBuilder_wSmearing_;

  // CV: don't delete TTree objects, as these are handled by TFileService

  delete ntupleFiller_piPlus_piMinus_;
  delete ntupleFiller_piPlus_rhoMinus_;
  delete ntupleFiller_rhoPlus_piMinus_;
  delete ntupleFiller_rhoPlus_rhoMinus_;
  delete ntupleFiller_had_had_;
}

void EntanglementNtupleProducer::beginJob()
{
  edm::Service<TFileService> fs;

  ntuple_piPlus_piMinus_ = fs->make<TTree>("piPlus_piMinus", "piPlus_piMinus");
  ntupleFiller_piPlus_piMinus_ = new EntanglementNtuple(ntuple_piPlus_piMinus_);
  ntuple_piPlus_rhoMinus_ = fs->make<TTree>("piPlus_rhoMinus", "piPlus_rhoMinus");
  ntupleFiller_piPlus_rhoMinus_ = new EntanglementNtuple(ntuple_piPlus_rhoMinus_);
  ntuple_rhoPlus_piMinus_ = fs->make<TTree>("rhoPlus_piMinus", "rhoPlus_piMinus");
  ntupleFiller_rhoPlus_piMinus_ = new EntanglementNtuple(ntuple_rhoPlus_piMinus_);
  ntuple_rhoPlus_rhoMinus_ = fs->make<TTree>("rhoPlus_rhoMinus", "rhoPlus_rhoMinus");
  ntupleFiller_rhoPlus_rhoMinus_ = new EntanglementNtuple(ntuple_rhoPlus_rhoMinus_);

  ntuple_had_had_ = fs->make<TTree>("had_had", "had_had");
  ntupleFiller_had_had_ = new EntanglementNtuple(ntuple_had_had_);
}

void EntanglementNtupleProducer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<EntanglementNtupleProducer::analyze>:\n";
    std::cout << " moduleLabel = " << moduleLabel_ << "\n";
  }

  double evtWeight = 1.;
  for ( const edm::EDGetTokenT<double>& tokenWeight : tokenWeights_ )
  {
    edm::Handle<double> weight;
    evt.getByToken(tokenWeight, weight);
    evtWeight *= (*weight);
  }
  if ( verbosity_ >= 1 )
  {
    std::cout << "evtWeight = " << evtWeight << "\n";
  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(token_, genParticles);

  KinematicEvent kineEvt_gen = (*genKineEvtBuilder_woSmearing_)(*genParticles);
  if ( verbosity_ >= 1 )
  {
    printKinematicEvent("kineEvt_gen", kineEvt_gen, cartesian_);
  }

  const reco::GenParticle* tauPlus  = nullptr;
  const reco::GenParticle* tauMinus = nullptr;
  for ( const reco::GenParticle& genParticle : *genParticles )
  {
    if ( genParticle.pdgId() == -15 && !tauPlus  ) tauPlus  = findLastTau(&genParticle);
    if ( genParticle.pdgId() == +15 && !tauMinus ) tauMinus = findLastTau(&genParticle);
  }
  if ( !(tauPlus && tauMinus) ) 
  {
    std::cerr << "WARNING: Failed to find tau+ tau- pair --> skipping the event !!\n";
    return;
  }

  std::vector<const reco::GenParticle*> tauPlus_daughters;
  findDecayProducts(tauPlus, tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_ch = get_chargedHadrons(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_pi0 = get_neutralPions(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_nu = get_neutrinos(tauPlus_daughters);
  int tauPlus_decaymode = get_decayMode(tauPlus_ch, tauPlus_pi0, tauPlus_nu);
  int tauPlus_nChargedKaons = get_chargedKaons(tauPlus_daughters).size();
  int tauPlus_nNeutralKaons = get_neutralKaons(tauPlus_daughters).size();
  std::vector<const reco::GenParticle*> tauPlus_y = get_photons(tauPlus_daughters);
  int tauPlus_nPhotons = tauPlus_y.size();
  double tauPlus_sumPhotonEn = compVisP4(tauPlus_y).energy();

  std::vector<const reco::GenParticle*> tauMinus_daughters;
  findDecayProducts(tauMinus, tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_ch = get_chargedHadrons(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_pi0 = get_neutralPions(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_nu = get_neutrinos(tauMinus_daughters);
  int tauMinus_decaymode = get_decayMode(tauMinus_ch, tauMinus_pi0, tauMinus_nu);
  int tauMinus_nChargedKaons = get_chargedKaons(tauMinus_daughters).size();
  int tauMinus_nNeutralKaons = get_neutralKaons(tauMinus_daughters).size();
  std::vector<const reco::GenParticle*> tauMinus_y = get_photons(tauMinus_daughters);
  int tauMinus_nPhotons = tauMinus_y.size();
  double tauMinus_sumPhotonEn = compVisP4(tauMinus_y).energy();

  KinematicEvent kineEvt_gen_smeared = (*genKineEvtBuilder_wSmearing_)(*genParticles);
  if ( verbosity_ >= 1 )
  {
    printKinematicEvent("kineEvt_gen_smeared", kineEvt_gen_smeared, cartesian_);
  }

  KinematicEvent kineEvt_startPos = kineFitStartPosFinder_(kineEvt_gen_smeared);
  if ( verbosity_ >= 1 )
  {
    printKinematicEvent("kineEvt_startPos", kineEvt_startPos, cartesian_);
  }

  KinematicEvent kineEvt_kinFit = kineFit_(kineEvt_startPos);
  if ( verbosity_ >= 1 )
  {
    printKinematicEvent("kineEvt_kinFit", kineEvt_kinFit, cartesian_);
  }

  EntanglementNtuple* ntupleFiller = nullptr;
  if ( tauPlus_decaymode == reco::PFTau::kOneProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    ntupleFiller = ntupleFiller_piPlus_piMinus_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero )
  {
    ntupleFiller = ntupleFiller_piPlus_rhoMinus_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng1PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    ntupleFiller = ntupleFiller_rhoPlus_piMinus_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng1PiZero && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero )
  {
    ntupleFiller = ntupleFiller_rhoPlus_rhoMinus_;
  }
  if ( ntupleFiller )
  {
    ntupleFiller->fillBranches(
      evt, 
      &kineEvt_gen,
      tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons, tauPlus_sumPhotonEn,
      tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons, tauMinus_sumPhotonEn,
      &kineEvt_gen_smeared, &kineEvt_startPos, &kineEvt_kinFit,
      evtWeight);
  }
  ntupleFiller_had_had_->fillBranches(
    evt, 
    &kineEvt_gen,
    tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons, tauPlus_sumPhotonEn,
    tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons, tauMinus_sumPhotonEn,
    &kineEvt_gen_smeared, &kineEvt_startPos, &kineEvt_kinFit,
    evtWeight);
}

void EntanglementNtupleProducer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(EntanglementNtupleProducer);
