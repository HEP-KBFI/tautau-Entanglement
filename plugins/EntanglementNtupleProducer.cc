#include "TauAnalysis/Entanglement/plugins/EntanglementNtupleProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"                 // edm::Service<>
#include "CommonTools/UtilAlgos/interface/TFileService.h"             // TFileService

#include "DataFormats/Common/interface/Handle.h"                      // edm::Handle<>

#include "DataFormats/TauReco/interface/PFTau.h"                      // reco::PFTau::hadronicDecayMode

#include "TauAnalysis/Entanglement/interface/cmsException.h"          // cmsException
#include "TauAnalysis/Entanglement/interface/comp_visP4.h"            // comp_visP4()
#include "TauAnalysis/Entanglement/interface/constants.h"             // kLHC, kSuperKEKB
#include "TauAnalysis/Entanglement/interface/findDecayProducts.h"     // findDecayProducts()
#include "TauAnalysis/Entanglement/interface/findLastTau.h"           // findLastTau()
#include "TauAnalysis/Entanglement/interface/get_decayMode.h"         // get_decayMode()
#include "TauAnalysis/Entanglement/interface/get_particles_of_type.h" // get_chargedKaons(), get_neutralKaons(), get_photons()

#include <iostream>                                                   // std::cout

EntanglementNtupleProducer::EntanglementNtupleProducer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , genKineEvtBuilder_woSmearing_(nullptr)
  , genKineEvtBuilder_wSmearing_(nullptr)
  , kinematicFit_(nullptr)
  , ntuple_pi_pi_(nullptr)
  , ntupleFiller_pi_pi_(nullptr)
  , ntuple_pi_rho_(nullptr)
  , ntupleFiller_pi_rho_(nullptr)
  , ntuple_pi_a1_(nullptr)
  , ntupleFiller_pi_a1_(nullptr)
  , ntuple_rho_rho_(nullptr)
  , ntupleFiller_rho_rho_(nullptr)
  , ntuple_rho_a1_(nullptr)
  , ntupleFiller_rho_a1_(nullptr)
  , ntuple_a1_a1_(nullptr)
  , ntupleFiller_a1_a1_(nullptr)
  , ntuple_had_had_(nullptr)
  , ntupleFiller_had_had_(nullptr)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenParticleCollection>(src_);

  edm::ParameterSet cfg_resolutions = cfg.getParameter<edm::ParameterSet>("resolutions");

  bool applySmearing = cfg.getParameter<bool>("applySmearing");

  edm::ParameterSet cfg_woSmearing = cfg;
  cfg_woSmearing.addParameter<edm::ParameterSet>("resolutions", cfg_resolutions);
  cfg_woSmearing.addParameter<bool>("applySmearing", false);
  genKineEvtBuilder_woSmearing_ = new GenKinematicEventBuilder(cfg_woSmearing);

  edm::ParameterSet cfg_wSmearing = cfg;
  cfg_wSmearing.addParameter<edm::ParameterSet>("resolutions", cfg_resolutions);
  cfg_wSmearing.addParameter<bool>("applySmearing", applySmearing);
  genKineEvtBuilder_wSmearing_ = new GenKinematicEventBuilder(cfg_wSmearing);
  
  std::string hAxis = cfg.getParameter<std::string>("hAxis");

  std::string collider = cfg.getParameter<std::string>("collider");
  if      ( collider == "LHC"       ) collider_ = kLHC;
  else if ( collider == "SuperKEKB" ) collider_ = kSuperKEKB;
  else throw cmsException("EntanglementNtupleProducer", __LINE__)
    << "Invalid Configuration parameter 'collider' = " << collider << " !!\n";

  edm::ParameterSet cfg_startPosFinder = cfg.getParameter<edm::ParameterSet>("startPosFinder");
  cfg_startPosFinder.addParameter<edm::ParameterSet>("resolutions", cfg_resolutions);
  cfg_startPosFinder.addParameter<std::string>("hAxis", hAxis);
  cfg_startPosFinder.addParameter<std::string>("collider", collider);
  cfg_startPosFinder.addUntrackedParameter<int>("verbosity", verbosity_);
  cfg_startPosFinder.addUntrackedParameter<bool>("cartesian", cartesian_);
  std::vector<int> startPos_algos = cfg_startPosFinder.getParameter<std::vector<int>>("algos");
  for ( int startPos_algo : startPos_algos )
  {
    edm::ParameterSet cfg_startPosFinder_algo = cfg_startPosFinder;
    cfg_startPosFinder_algo.addParameter<int>("algo", startPos_algo);
    startPosFinders_.push_back(new StartPosFinder(cfg_startPosFinder_algo));
  }

  edm::ParameterSet cfg_kinematicFit = cfg.getParameter<edm::ParameterSet>("kinematicFit");
  cfg_kinematicFit.addParameter<edm::ParameterSet>("resolutions", cfg_resolutions);
  cfg_kinematicFit.addParameter<std::string>("hAxis", hAxis);
  cfg_kinematicFit.addParameter<std::string>("collider", collider);
  cfg_kinematicFit.addUntrackedParameter<int>("verbosity", verbosity_);
  cfg_kinematicFit.addUntrackedParameter<bool>("cartesian", cartesian_);
  kinematicFit_ = new KinematicFit(cfg_kinematicFit);

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

  for ( StartPosFinder* startPosFinder : startPosFinders_ )
  {
    delete startPosFinder;
  }

  delete kinematicFit_;

  // CV: don't delete TTree objects, as these are handled by TFileService

  delete ntupleFiller_pi_pi_;
  delete ntupleFiller_pi_rho_;
  delete ntupleFiller_pi_a1_;
  delete ntupleFiller_rho_rho_;
  delete ntupleFiller_rho_a1_;
  delete ntupleFiller_a1_a1_;
  delete ntupleFiller_had_had_;
}

void EntanglementNtupleProducer::beginJob()
{
  edm::Service<TFileService> fs;

  ntuple_pi_pi_ = fs->make<TTree>("pi_pi", "pi_pi");
  ntupleFiller_pi_pi_ = new EntanglementNtuple(ntuple_pi_pi_);
  ntuple_pi_rho_ = fs->make<TTree>("pi_rho", "pi_rho");
  ntupleFiller_pi_rho_ = new EntanglementNtuple(ntuple_pi_rho_);
  ntuple_pi_a1_ = fs->make<TTree>("pi_a1", "pi_a1");
  ntupleFiller_pi_a1_ = new EntanglementNtuple(ntuple_pi_a1_);
  ntuple_rho_rho_ = fs->make<TTree>("rho_rho", "rho_rho");
  ntupleFiller_rho_rho_ = new EntanglementNtuple(ntuple_rho_rho_);
  ntuple_rho_a1_ = fs->make<TTree>("rho_a1", "rho_a1");
  ntupleFiller_rho_a1_ = new EntanglementNtuple(ntuple_rho_a1_);
  ntuple_a1_a1_ = fs->make<TTree>("a1_a1", "a1_a1");
  ntupleFiller_a1_a1_ = new EntanglementNtuple(ntuple_a1_a1_);
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
    printKinematicEvent("kineEvt_gen", kineEvt_gen, verbosity_, cartesian_);
  }

  KinematicEvent kineEvt_gen_smeared = (*genKineEvtBuilder_wSmearing_)(*genParticles);
  if ( verbosity_ >= 1 )
  {
    printKinematicEvent("kineEvt_gen_smeared", kineEvt_gen_smeared, verbosity_, cartesian_);
  }

  KinematicEvent kineEvt_startPos_bestfit;
  KinematicEvent kineEvt_kinFit_bestfit;
  double kinFitChi2_bestfit = -1.;
  int kinFitStatus_bestfit = -1;
  bool isFirst = true;
  for ( StartPosFinder* startPosFinder : startPosFinders_ )
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "executing startPosFinder algo #" << startPosFinder->get_algo() << "...";
    }
    std::vector<KinematicEvent> kineEvts_startPos = (*startPosFinder)(kineEvt_gen_smeared);

    size_t numSolutions = kineEvts_startPos.size();
    if ( verbosity_ >= 1 )
    {
      std::cout << "#solutions = " << numSolutions<< "\n";
    }
    for ( size_t idxSolution = 0; idxSolution < numSolutions; ++idxSolution )
    {
      const KinematicEvent& kineEvt_startPos = kineEvts_startPos.at(idxSolution);

      if ( verbosity_ >= 1 )
      {
        std::cout << "executing KinematicFit for startPosFinder algo #" << startPosFinder->get_algo() << ", solution #" << idxSolution << "...\n";
        printKinematicEvent("kineEvt_startPos", kineEvt_startPos, verbosity_, cartesian_);
      }
      KinematicEvent kineEvt_kinFit = (*kinematicFit_)(kineEvt_startPos);

      if ( verbosity_ >= 1 )
      {
        printKinematicEvent("kineEvt_kinFit", kineEvt_kinFit, verbosity_, cartesian_);
      }
 
      bool isBetterFit = kineEvt_kinFit.kinFitStatus() >  kinFitStatus_bestfit || 
                        (kineEvt_kinFit.kinFitStatus() == kinFitStatus_bestfit && kineEvt_kinFit.kinFitChi2() < kinFitChi2_bestfit);
      if ( verbosity_ >= 1 )
      {
        std::cout << "isBetterFit = " << isBetterFit << "\n";
      }
      if ( isFirst || isBetterFit )
      {
        kineEvt_startPos_bestfit = kineEvt_startPos;
        kineEvt_kinFit_bestfit = kineEvt_kinFit;
        kinFitChi2_bestfit = kineEvt_kinFit.kinFitChi2();
        kinFitStatus_bestfit = kineEvt_kinFit.kinFitStatus();
        isFirst = false;
      }
    }
  }
  if ( !(kinFitStatus_bestfit == 1 || (kinFitStatus_bestfit == 0 && kinFitChi2_bestfit < 1.e+2)) )
  {
    std::cerr << "WARNING: KinematicFit failed to converge !!" << std::endl;
  }

  if ( verbosity_ >= 1 )
  {
    std::cout << "best fit:\n";
    printKinematicEvent("kineEvt_startPos", kineEvt_startPos_bestfit, verbosity_, cartesian_);
    printKinematicEvent("kineEvt_kinFit", kineEvt_kinFit_bestfit, verbosity_, cartesian_);
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
  double tauPlus_sumPhotonEn = comp_visP4(tauPlus_y).energy();

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
  double tauMinus_sumPhotonEn = comp_visP4(tauMinus_y).energy();

  EntanglementNtuple* ntupleFiller = nullptr;
  if ( tauPlus_decaymode == reco::PFTau::kOneProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    ntupleFiller = ntupleFiller_pi_pi_;
  }
  else if ( (tauPlus_decaymode == reco::PFTau::kOneProng0PiZero   && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero  ) ||
            (tauPlus_decaymode == reco::PFTau::kOneProng1PiZero   && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero  ) )
  {
    ntupleFiller = ntupleFiller_pi_rho_;
  }
  else if ( (tauPlus_decaymode == reco::PFTau::kOneProng0PiZero   && tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero) ||
            (tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero  ) )
  {
    ntupleFiller = ntupleFiller_pi_a1_;
  }
  else if (  tauPlus_decaymode == reco::PFTau::kOneProng1PiZero   && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero  )
  {
    ntupleFiller = ntupleFiller_rho_rho_;
  }
  else if ( (tauPlus_decaymode == reco::PFTau::kOneProng1PiZero   && tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero) ||
            (tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero  ) )
  {
    ntupleFiller = ntupleFiller_rho_a1_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero && tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero )
  {
    ntupleFiller = ntupleFiller_a1_a1_;
  }
  
  if ( ntupleFiller )
  {
    ntupleFiller->fillBranches(
      evt, 
      &kineEvt_gen,
      tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons, tauPlus_sumPhotonEn,
      tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons, tauMinus_sumPhotonEn,
      &kineEvt_gen_smeared, &kineEvt_startPos_bestfit, &kineEvt_kinFit_bestfit,
      evtWeight);
  }
  ntupleFiller_had_had_->fillBranches(
    evt, 
    &kineEvt_gen,
    tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons, tauPlus_sumPhotonEn,
    tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons, tauMinus_sumPhotonEn,
    &kineEvt_gen_smeared, &kineEvt_startPos_bestfit, &kineEvt_kinFit_bestfit,
    evtWeight);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(EntanglementNtupleProducer);
