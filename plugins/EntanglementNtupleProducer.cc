#include "TauAnalysis/Entanglement/plugins/EntanglementNtupleProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // reco::Candidate::LorentzVector, reco::Candidate::Vector
#include "DataFormats/TauReco/interface/PFTau.h"                          // reco::PFTau::hadronicDecayMode

#include "TauAnalysis/Entanglement/interface/auxFunctions.h"              // printLorentzVector(), printVector(), square()
#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/constants.h"                 // mProton, mTau, gamma_va
#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"              // SpinAnalyzer::kTauPlus, SpinAnalyzer::kTauMinus

#include <Math/Boost.h>                                                   // Boost

#include <iostream>
#include <iomanip>
#include <cmath>

EntanglementNtupleProducer::EntanglementNtupleProducer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , mode_(-1)
  , resolutions_(nullptr)
  , smearing_(cfg)
  , applySmearing_(cfg.getParameter<bool>("applySmearing"))
  , kineFit_(cfg)
  , spinAnalyzerOneProng0Pi0_(cfg)
  , spinAnalyzerOneProng1Pi0_(cfg)
  , spinAnalyzerThreeProng0Pi0_(cfg)
  , ntuple_piPlus_piMinus_(nullptr)
  , branches_piPlus_piMinus_(reco::PFTau::kOneProng0PiZero, reco::PFTau::kOneProng0PiZero)
  , ntuple_piPlus_rhoMinus_(nullptr)
  , branches_piPlus_rhoMinus_(reco::PFTau::kOneProng0PiZero, reco::PFTau::kOneProng1PiZero)
  , ntuple_rhoPlus_piMinus_(nullptr)
  , branches_rhoPlus_piMinus_(reco::PFTau::kOneProng1PiZero, reco::PFTau::kOneProng0PiZero)
  , ntuple_rhoPlus_rhoMinus_(nullptr)
  , branches_rhoPlus_rhoMinus_(reco::PFTau::kOneProng1PiZero, reco::PFTau::kOneProng1PiZero)
  , ntuple_had_had_(nullptr)
  , branches_had_had_(reco::PFTau::kNull, reco::PFTau::kNull)
  , verbosity_(-1)
  , cartesian_(false)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenParticleCollection>(src_);

  std::string mode = cfg.getParameter<std::string>("mode");
  if      ( mode == "gen" ) mode_ = kGen;
  else if ( mode == "rec" ) mode_ = kRec;
  else throw cmsException("EntanglementNtupleProducer", __LINE__)
    << "Invalid Configuration parameter 'mode' = " << mode << " !!\n";

  edm::ParameterSet cfg_resolutions = cfg.getParameterSet("resolutions");
  resolutions_ = new Resolutions(cfg_resolutions);

  srcWeights_ = cfg.getParameter<vInputTag>("srcEvtWeights");
  for ( const edm::InputTag& srcWeight : srcWeights_ )
  {
    tokenWeights_.push_back(consumes<double>(srcWeight));
  }

  verbosity_ = cfg.getUntrackedParameter<int>("verbosity");
  cartesian_ = cfg.getUntrackedParameter<bool>("cartesian");
}

EntanglementNtupleProducer::~EntanglementNtupleProducer()
{
  delete resolutions_;
}

void EntanglementNtupleProducer::beginJob()
{
  edm::Service<TFileService> fs;

  ntuple_piPlus_piMinus_ = fs->make<TTree>("piPlus_piMinus", "piPlus_piMinus");
  branches_piPlus_piMinus_.initBranches(ntuple_piPlus_piMinus_);
  ntuple_piPlus_rhoMinus_ = fs->make<TTree>("piPlus_rhoMinus", "piPlus_rhoMinus");
  branches_piPlus_rhoMinus_.initBranches(ntuple_piPlus_rhoMinus_);
  ntuple_rhoPlus_piMinus_ = fs->make<TTree>("rhoPlus_piMinus", "rhoPlus_piMinus");
  branches_rhoPlus_piMinus_.initBranches(ntuple_rhoPlus_piMinus_);
  ntuple_rhoPlus_rhoMinus_ = fs->make<TTree>("rhoPlus_rhoMinus", "rhoPlus_rhoMinus");
  branches_rhoPlus_rhoMinus_.initBranches(ntuple_rhoPlus_rhoMinus_);

  ntuple_had_had_ = fs->make<TTree>("had_had", "had_had");
  branches_had_had_.initBranches(ntuple_had_had_);
}

namespace
{
  const reco::GenParticle*
  findLastTau(const reco::GenParticle* mother)
  {
    size_t numDaughters = mother->numberOfDaughters();
    for ( size_t idxDaughter = 0; idxDaughter < numDaughters; ++idxDaughter )
    {
      const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(mother->daughter(idxDaughter));
      assert(daughter);
      if ( daughter->pdgId() == mother->pdgId() )
      {
        return findLastTau(daughter);
      }
    }
    return mother;
  }

  void
  findDecayProducts(const reco::GenParticle* mother, 
                    std::vector<const reco::GenParticle*>& daughters)
  {
    size_t numDaughters = mother->numberOfDaughters();
    for ( size_t idxDaughter = 0; idxDaughter < numDaughters; ++idxDaughter )
    {
      const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(mother->daughter(idxDaughter));
      assert(daughter);
      // CV: treat neutral pions as stable particles
      if ( daughter->status() == 1 || daughter->pdgId() == 111 )
      {
        daughters.push_back(daughter);
      }
      else 
      {
        findDecayProducts(daughter, daughters);
      }
    }
  }

  std::vector<const reco::GenParticle*>
  get_chargedHadrons(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    std::vector<const reco::GenParticle*> chargedHadrons;
    for ( const reco::GenParticle* decayProduct : decayProducts )
    {
      int absPdgId = std::abs(decayProduct->pdgId());
      if ( absPdgId == 11 || absPdgId == 13 ) continue;
      if ( std::fabs(decayProduct->charge()) > 0.5 )
      {
        chargedHadrons.push_back(decayProduct);
      }
    }
    return chargedHadrons;
  }
  
  std::vector<const reco::GenParticle*>
  get_particles_of_type(const std::vector<const reco::GenParticle*>& decayProducts, const std::vector<int>& selPdgIds)
  {
    std::vector<const reco::GenParticle*> selDecayProducts;
    for ( const reco::GenParticle* decayProduct : decayProducts )
    {
      int absPdgId = std::abs(decayProduct->pdgId());
      bool isSelPdgId = false;
      for ( int selPdgId : selPdgIds )
      {
        if ( absPdgId == selPdgId )
        {
          isSelPdgId = true;
          break;
        }
      }
      if ( isSelPdgId )
      {
        selDecayProducts.push_back(decayProduct);
      }
    }
    return selDecayProducts;
  }

  std::vector<const reco::GenParticle*>
  get_neutralPions(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    return get_particles_of_type(decayProducts, { 111 });
  }

  std::vector<const reco::GenParticle*>
  get_neutrinos(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    return get_particles_of_type(decayProducts, { 12, 14, 16 });
  }

  std::vector<const reco::GenParticle*>
  get_chargedKaons(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    return get_particles_of_type(decayProducts, { 321 });
  }

  std::vector<const reco::GenParticle*>
  get_neutralKaons(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    return get_particles_of_type(decayProducts, { 130, 310, 311 });
  }

  std::vector<const reco::GenParticle*>
  get_photons(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    return get_particles_of_type(decayProducts, { 22 });
  }

  void
  printGenParticles(const std::vector<const reco::GenParticle*>& decayProducts,
                    bool cartesian = true)
  {
    size_t numDecayProducts = decayProducts.size();
    for ( size_t idxDecayProduct = 0; idxDecayProduct < numDecayProducts; ++idxDecayProduct )
    {
      const reco::GenParticle* decayProduct = decayProducts[idxDecayProduct];
      std::cout << "GenParticle #" << idxDecayProduct << ":";
      reco::Candidate::LorentzVector p4 = decayProduct->p4();
      if ( cartesian )
      {
        std::cout << " E = " << p4.energy() << ", Px = " << p4.px() << ", Py = " << p4.py() << ", Pz = " << p4.pz();
      }
      else
      {
        std::cout << " pT = " << p4.pt() << ", eta = " << p4.eta() << ", phi = " << p4.phi() << ", mass = " << p4.mass();
      }
      std::cout << ", pdgId = " << decayProduct->pdgId() << "\n";
    } 
  }

  reco::Candidate::LorentzVector
  compVisP4(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    reco::Candidate::LorentzVector visP4;
    for ( const reco::GenParticle* decayProduct : decayProducts )
    {
      int absPdgId = std::abs(decayProduct->pdgId());
      if ( !(absPdgId == 12 || absPdgId == 14 || absPdgId == 16) )
      {
        visP4 += decayProduct->p4();
      }
    }
    return visP4;
  }

  int
  get_decayMode(const std::vector<const reco::GenParticle*>& tau_ch,
                const std::vector<const reco::GenParticle*>& tau_pi0,
                const std::vector<const reco::GenParticle*>& tau_nu)
  {
    size_t nch = tau_ch.size();
    size_t npi0 = tau_pi0.size();
    size_t nnu = tau_nu.size();
    // CV: there are rare (on the level of 10^-4) events in the official CMS Monte Carlo samples
    //     in which hadronic tau decays involve more than one neutrino;
    //     let's for now simply discard those
    if      ( nch == 1 && npi0 == 0 && nnu == 1 ) return reco::PFTau::kOneProng0PiZero;
    else if ( nch == 1 && npi0 == 1 && nnu == 1 ) return reco::PFTau::kOneProng1PiZero;
    else if ( nch == 1 && npi0 == 2 && nnu == 1 ) return reco::PFTau::kOneProng2PiZero;
    else if ( nch == 3 && npi0 == 0 && nnu == 1 ) return reco::PFTau::kThreeProng0PiZero;
    else                                          return reco::PFTau::kNull;
  }

  const reco::GenParticle*
  get_leadTrack(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    const reco::GenParticle* leadTrack = nullptr;
    double max_pt = -1.;
    for ( const reco::GenParticle* decayProduct : decayProducts )
    {
      if ( std::fabs(decayProduct->charge()) > 0.5 && decayProduct->pt() > max_pt )
      {
        leadTrack = decayProduct;
        max_pt = decayProduct->pt();
      }
    }
    return leadTrack;
  }

  std::vector<KinematicParticle>
  get_kineDaughters(const std::vector<const reco::GenParticle*>& decayProducts, const Resolutions& resolutions)
  {
    std::vector<KinematicParticle> kineDaughters;
    for ( const reco::GenParticle* decayProduct : decayProducts )
    {
      KinematicParticle kineDaughter(decayProduct->pdgId());
      TVectorD params7(7);
      params7(0) = decayProduct->px();
      params7(1) = decayProduct->py();
      params7(2) = decayProduct->pz();
      params7(3) = decayProduct->energy();
      params7(4) = decayProduct->vertex().x();
      params7(5) = decayProduct->vertex().y();
      params7(6) = decayProduct->vertex().z();
      // CV: the covariance matrix is only a crude approximation;
      //     to be implemented in more detail later !!
      TMatrixD cov7x7(7,7);
      reco::Candidate::Vector r, n, k; 
      get_localCoordinateSystem(decayProduct->p4(), nullptr, nullptr, kBeam, r, n, k);
      double dk = 0.;
      double dr = resolutions.get_tipResolutionPerp();
      double dn = resolutions.get_tipResolutionPerp();
      reco::Candidate::Vector x(1., 0., 0.);
      reco::Candidate::Vector y(0., 1., 0.);
      reco::Candidate::Vector z(0., 0., 1.);
      double dx = dr*r.Dot(x) + dn*n.Dot(x) + dk*k.Dot(x);
      double dy = dr*r.Dot(y) + dn*n.Dot(y) + dk*k.Dot(y);
      double dz = dr*r.Dot(z) + dn*n.Dot(z) + dk*k.Dot(z);
      cov7x7(4,4) = square(dr*r.Dot(x) + dn*n.Dot(x) + dk*k.Dot(x));
      cov7x7(4,5) = dx*dy;
      cov7x7(4,6) = dx*dz;
      cov7x7(5,4) = dy*dx;
      cov7x7(5,5) = square(dy);
      cov7x7(5,6) = dy*dz;
      cov7x7(6,4) = dz*dx;
      cov7x7(6,5) = dz*dy;
      cov7x7(6,6) = square(dz);
      kineDaughter.set_params7(params7, cov7x7);
      kineDaughters.push_back(kineDaughter);
    }
    return kineDaughters;
  }

  reco::Candidate::Point
  get_tipPCA(const reco::Candidate::Point& pv, const reco::GenParticle* leadTrack)
  {
    reco::Candidate::Point sv = leadTrack->vertex();
    auto flightlength = sv - pv;
    auto e_trk = leadTrack->momentum().unit();
    reco::Candidate::Point tipPCA = pv + flightlength - flightlength.Dot(e_trk)*e_trk;
    return tipPCA;
  }

  reco::Candidate::Point
  get_linearPCA(const reco::Candidate::Point& pv,
                const reco::Candidate::Vector& tauP3,
                const reco::Candidate::Point& sv,
                const reco::Candidate::Vector& visTauP3,
                int verbosity)
  {
    // CV: code based on https://math.stackexchange.com/questions/1993953/closest-points-between-two-lines
    reco::Candidate::Vector d = tauP3.Cross(visTauP3).unit();
    TMatrixD v(3,3);
    v(0,0) =  tauP3.x();
    v(0,1) = -visTauP3.x();
    v(0,2) =  d.x();
    v(1,0) =  tauP3.y();
    v(1,1) = -visTauP3.y();
    v(1,2) =  d.y();
    v(2,0) =  tauP3.z();
    v(2,1) = -visTauP3.z();
    v(2,2) =  d.z();
    TVectorD r(3);
    r(0) = -pv.x() + sv.x();
    r(1) = -pv.y() + sv.y();
    r(2) = -pv.z() + sv.z();
    TVectorD lambda = v.Invert()*r;
    reco::Candidate::Point pca1 = pv + lambda(0)*tauP3;
    reco::Candidate::Point pca2 = sv + lambda(1)*visTauP3;
    if ( verbosity >= 1 )
    {
      printPoint("pca1", pca1);
      printPoint("pca2", pca2);
      std::cout << "|d| = " << lambda(2)*std::sqrt(d.mag2()) << std::endl;
    }
    return pca1;
  }

  double
  get_z(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visTauP4, const ROOT::Math::Boost& boost_ttrf)
  {
    reco::Candidate::LorentzVector tauP4_ttrf = getP4_rf(tauP4, boost_ttrf);
    reco::Candidate::LorentzVector visTauP4_ttrf = getP4_rf(visTauP4, boost_ttrf);
    double z = visTauP4_ttrf.energy()/tauP4_ttrf.energy();
    return z;
  }
}

void EntanglementNtupleProducer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<EntanglementNtupleProducer::analyze>:" << std::endl;
    std::cout << " moduleLabel = " << moduleLabel_ << std::endl;
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
    std::cout << "evtWeight = " << evtWeight << std::endl;
  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(token_, genParticles);

  const reco::GenParticle* tauPlus  = nullptr;
  const reco::GenParticle* tauMinus = nullptr;
  for ( const reco::GenParticle& genParticle : *genParticles )
  {
    if ( genParticle.pdgId() == -15 && !tauPlus  ) tauPlus  = findLastTau(&genParticle);
    if ( genParticle.pdgId() == +15 && !tauMinus ) tauMinus = findLastTau(&genParticle);
  }

  if ( !(tauPlus && tauMinus) ) 
  {
    std::cerr << "WARNING: Failed to find tau+ tau- pair --> skipping the event !!" << std::endl;
    return;
  }
  reco::Candidate::LorentzVector tauPlusP4 = tauPlus->p4();
  reco::Candidate::LorentzVector tauMinusP4 = tauMinus->p4();

  std::vector<const reco::GenParticle*> tauPlus_daughters;
  findDecayProducts(tauPlus, tauPlus_daughters);
  reco::Candidate::LorentzVector visTauPlusP4 = compVisP4(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_ch = get_chargedHadrons(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_pi0 = get_neutralPions(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_nu = get_neutrinos(tauPlus_daughters);
  int tauPlus_decaymode = get_decayMode(tauPlus_ch, tauPlus_pi0, tauPlus_nu);
  if ( verbosity_ >= 1 )
  {
    printLorentzVector("tauPlusP4", tauPlusP4, cartesian_);
    if ( cartesian_ )
    {
      std::cout << " mass = " << tauPlusP4.mass() << "\n";
    }
    std::cout << "tau+ decay products:" << "\n";
    printGenParticles(tauPlus_daughters, cartesian_);
    printLorentzVector("visTauPlusP4", visTauPlusP4, cartesian_);
    if ( cartesian_ )
    {
      std::cout << " mass = " << visTauPlusP4.mass() << "\n";
    }
    std::cout << "tauPlus_decaymode = " << tauPlus_decaymode << "\n";
  }
  int tauPlus_nChargedKaons = get_chargedKaons(tauPlus_daughters).size();
  int tauPlus_nNeutralKaons = get_neutralKaons(tauPlus_daughters).size();
  std::vector<const reco::GenParticle*> tauPlus_y = get_photons(tauPlus_daughters);
  int tauPlus_nPhotons = tauPlus_y.size();
  double tauPlus_sumPhotonEn = compVisP4(tauPlus_y).energy();

  std::vector<const reco::GenParticle*> tauMinus_daughters;
  findDecayProducts(tauMinus, tauMinus_daughters);
  reco::Candidate::LorentzVector visTauMinusP4 = compVisP4(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_ch = get_chargedHadrons(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_pi0 = get_neutralPions(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_nu = get_neutrinos(tauMinus_daughters);
  int tauMinus_decaymode = get_decayMode(tauMinus_ch, tauMinus_pi0, tauMinus_nu);
  if ( verbosity_ >= 1 )
  {
    printLorentzVector("tauMinusP4", tauMinusP4, cartesian_);
    if ( cartesian_ )
    {
      std::cout << " mass = " << tauPlusP4.mass() << "\n";
    }
    std::cout << "tau- decay products:" << "\n";
    printGenParticles(tauMinus_daughters, cartesian_);
    printLorentzVector("visTauMinusP4", visTauMinusP4, cartesian_);
    if ( cartesian_ )
    {
      std::cout << " mass = " << visTauMinusP4.mass() << "\n";
    }
    std::cout << "tauMinus_decaymode = " << tauMinus_decaymode << "\n";
  }
  int tauMinus_nChargedKaons = get_chargedKaons(tauMinus_daughters).size();
  int tauMinus_nNeutralKaons = get_neutralKaons(tauMinus_daughters).size();
  std::vector<const reco::GenParticle*> tauMinus_y = get_photons(tauMinus_daughters);
  int tauMinus_nPhotons = tauMinus_y.size();
  double tauMinus_sumPhotonEn = compVisP4(tauMinus_y).energy();
 
  if ( !(std::fabs(tauPlus->vertex().x() - tauMinus->vertex().x()) < 1.e-3 &&
         std::fabs(tauPlus->vertex().y() - tauMinus->vertex().y()) < 1.e-3 &&
         std::fabs(tauPlus->vertex().z() - tauMinus->vertex().z()) < 1.e-3) )
  {
    std::cerr << "WARNING: Failed to find vertex of tau+ tau- pair --> skipping the event !!" << std::endl;
    return;
  }
  // CV: set primary event vertex (PV),
  //     using generator-level information
  reco::Candidate::Point pv = tauMinus->vertex();
  if ( verbosity_ >= 1 )
  {
    printPoint("GenPV", pv);
  }

  reco::Candidate::LorentzVector recoilP4 = tauPlusP4 + tauMinusP4;

  const reco::GenParticle* tauPlus_leadTrack = get_leadTrack(tauPlus_daughters);
  const reco::GenParticle* tauMinus_leadTrack = get_leadTrack(tauMinus_daughters);
  if ( !(tauPlus_leadTrack && tauMinus_leadTrack) )
  {
    std::cerr << "WARNING: Failed to find leading tracks of tau+ tau- pair --> skipping the event !!" << std::endl;
    return;
  }
  if ( verbosity_ >= 1 )
  {
    printPoint("GenSV(tau+)", tauPlus_ch[0]->vertex());
    printPoint("GenSV(tau-)", tauMinus_ch[0]->vertex());
  }

  std::vector<KinematicParticle> daughtersTauPlus = get_kineDaughters(tauPlus_daughters, *resolutions_);
  reco::Candidate::Point tipPCATauPlus = get_tipPCA(pv, tauPlus_leadTrack);
  std::vector<KinematicParticle> daughtersTauMinus = get_kineDaughters(tauMinus_daughters, *resolutions_);
  reco::Candidate::Point tipPCATauMinus = get_tipPCA(pv, tauMinus_leadTrack);
  KinematicEvent kineEvt(pv, recoilP4);
  kineEvt.set_visTauPlus(visTauPlusP4,  tauPlus_decaymode,  daughtersTauPlus,  tipPCATauPlus);
  kineEvt.set_visTauMinus(visTauMinusP4, tauMinus_decaymode, daughtersTauMinus, tipPCATauMinus);
  if ( mode_ == kGen )
  {
    // CV: set tau decay vertex for all taus,
    //     using generator-level information
    kineEvt.set_tauPlusP4(tauPlusP4);
    if ( tauPlus_ch.size() >= 1 )
    {  
      kineEvt.set_svTauPlus(tauPlus_ch[0]->vertex());
    }
    kineEvt.set_tauMinusP4(tauMinusP4);
    if ( tauMinus_ch.size() >= 1 )
    { 
      kineEvt.set_svTauMinus(tauMinus_ch[0]->vertex());
    }
  }
  else if ( mode_ == kRec )
  {
    // CV: use generator-level information to set tau decay vertex (SV) for three-prongs only,
    //     smear position of tau decay vertex to simulate experimental resolution
    //     using generator-level information,
    //     which subsequently gets smeared to simulate experimental resolution
    if ( tauPlus_ch.size() >= 3 )
    {
      kineEvt.set_svTauPlus(tauPlus_ch[0]->vertex());
    }
    if ( tauMinus_ch.size() >= 3 )
    {
      kineEvt.set_svTauMinus(tauMinus_ch[0]->vertex());
    }
    if ( applySmearing_ )
    {
      kineEvt = smearing_(kineEvt);
    }
    // CV: set tau decay vertex for one-prongs
    //     to point of closest approach between tau lepton momentum vector and track of charged pion
    if ( tauPlus_ch.size() < 3 )
    {
      reco::Candidate::Point svTauPlus = get_linearPCA(pv, tauPlus->momentum(), tauPlus_ch[0]->vertex(), tauPlus_ch[0]->momentum(), verbosity_);
      if ( verbosity_ >= 1 )
      {
        printPoint("RecSV(tau+)", svTauPlus);
      }
      kineEvt.set_svTauPlus(svTauPlus);
    }
    if ( tauMinus_ch.size() < 3 )
    {
      reco::Candidate::Point svTauMinus = get_linearPCA(pv, tauMinus->momentum(), tauMinus_ch[0]->vertex(), tauMinus_ch[0]->momentum(), verbosity_);
      if ( verbosity_ >= 1 )
      {
        printPoint("RecSV(tau-)", svTauMinus);
      }
      kineEvt.set_svTauMinus(svTauMinus);
    }
    KinematicEvent fitted_kineEvt = kineFit_(kineEvt);
    kineEvt.set_tauPlusP4(fitted_kineEvt.get_tauPlusP4());
    kineEvt.set_tauMinusP4(fitted_kineEvt.get_tauMinusP4());
  } else assert(0);
  if ( verbosity_ >= 1 )
  {
    printKinematicEvent("kineEvt", kineEvt, cartesian_);
  }

  reco::Candidate::Vector hPlus;
  if ( tauPlus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    hPlus = spinAnalyzerOneProng0Pi0_(kineEvt, SpinAnalyzer::kTauPlus);
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng1PiZero )
  {
    hPlus = spinAnalyzerOneProng1Pi0_(kineEvt, SpinAnalyzer::kTauPlus);
  }
  //else if ( tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero )
  //{
  //  hPlus = spinAnalyzerThreeProng0Pi0_(kineEvt, SpinAnalyzer::kTauPlus);
  //}

  reco::Candidate::Vector hMinus;
  if ( tauMinus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    hMinus = spinAnalyzerOneProng0Pi0_(kineEvt, SpinAnalyzer::kTauMinus);
  }
  else if ( tauMinus_decaymode == reco::PFTau::kOneProng1PiZero )
  {
    hMinus = spinAnalyzerOneProng1Pi0_(kineEvt, SpinAnalyzer::kTauMinus);
  }
  //else if ( tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero )
  //{
  //  hMinus = spinAnalyzerThreeProng0Pi0_(kineEvt, SpinAnalyzer::kTauMinus);
  //}

  branchType* branches = nullptr;
  if ( tauPlus_decaymode == reco::PFTau::kOneProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    branches = &branches_piPlus_piMinus_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero )
  {
    branches = &branches_piPlus_rhoMinus_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng1PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    branches = &branches_rhoPlus_piMinus_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng1PiZero && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero )
  {
    branches = &branches_rhoPlus_rhoMinus_;
  }
  if ( branches )
  {
    ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(recoilP4.BoostToCM());
    double zPlus = 0.;
    bool zPlus_isValid = false;
    double zMinus = 0.;
    bool zMinus_isValid = false;
    if ( kineEvt.get_tauPlusP4_isValid() )
    {
      zPlus = get_z(tauPlusP4, visTauPlusP4, boost_ttrf);
      zPlus_isValid = true;
    }
    if ( kineEvt.get_tauMinusP4_isValid() )
    {
      zMinus = get_z(tauMinusP4, visTauMinusP4, boost_ttrf);
      zMinus_isValid = true;
    }
    reco::Candidate::LorentzVector tauMinusP4_ttrf = getP4_rf(kineEvt.get_tauMinusP4(), boost_ttrf);
    double cosTheta = cos(tauMinusP4_ttrf.theta());
    if ( verbosity_ >= 1 )
    {
      printVector("hPlus", hPlus, cartesian_);
      if ( zPlus_isValid )
      {
        std::cout << "zPlus = " << zPlus << std::endl;
      }
      else
      {
        std::cout << "zPlus: N/A\n";
      }
      printVector("hMinus", hMinus, cartesian_);
      if ( zMinus_isValid )
      {
        std::cout << "zMinus = " << zMinus << std::endl;
      }
      else
      {
        std::cout << "zMinus: N/A\n";
      }
      std::cout << "C_rr = " << hPlus.x()*hMinus.x() << "\n";
      std::cout << "C_rn = " << hPlus.x()*hMinus.y() << "\n";
      std::cout << "C_rk = " << hPlus.x()*hMinus.z() << "\n";
      std::cout << "C_nr = " << hPlus.y()*hMinus.x() << "\n";
      std::cout << "C_nn = " << hPlus.y()*hMinus.y() << "\n";
      std::cout << "C_nk = " << hPlus.y()*hMinus.z() << "\n";
      std::cout << "C_kr = " << hPlus.z()*hMinus.x() << "\n";
      std::cout << "C_kn = " << hPlus.z()*hMinus.y() << "\n";
      std::cout << "C_kk = " << hPlus.z()*hMinus.z() << "\n";
      double mTauTau = (tauPlusP4 + tauMinusP4).mass();
      double mVis = (visTauPlusP4 + visTauMinusP4).mass();
      std::cout << "mTauTau = " << mTauTau << ", mVis = " << mVis << "\n";
      std::cout << "cosTheta = " << cosTheta << "\n";
    }
    branches->fillBranches(evt,
      hPlus, tauPlusP4, tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons, tauPlus_sumPhotonEn, visTauPlusP4, zPlus, 
      hMinus, tauMinusP4, tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons, tauMinus_sumPhotonEn, visTauMinusP4, zMinus, 
      cosTheta, evtWeight);

    branches_had_had_.tauPlus_decaymode_ = tauPlus_decaymode;
    branches_had_had_.tauMinus_decaymode_ = tauMinus_decaymode;
    branches_had_had_.fillBranches(evt,
      hPlus, tauPlusP4, tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons, tauPlus_sumPhotonEn, visTauPlusP4, zPlus, 
      hMinus, tauMinusP4, tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons, tauMinus_sumPhotonEn, visTauMinusP4, zMinus, 
      cosTheta, evtWeight);
  }
}

void EntanglementNtupleProducer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(EntanglementNtupleProducer);
