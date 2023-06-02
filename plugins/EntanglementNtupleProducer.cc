#include "TauAnalysis/Entanglement/plugins/EntanglementNtupleProducer.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::LorentzVector, reco::Candidate::Vector
#include "DataFormats/TauReco/interface/PFTau.h"       // reco::PFTau::hadronicDecayMode

#include <Math/Boost.h>                                // Boost

#include <iostream>
#include <iomanip>
#include <cmath>

const double mProton = 0.938272;
const double mTau = 1.77686;

const double gamma_va = 1.;

EntanglementNtupleProducer::EntanglementNtupleProducer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , hAxis_(-1)
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
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenParticleCollection>(src_);

  std::string hAxis = cfg.getParameter<std::string>("hAxis");
  if      ( hAxis == "beam"  ) hAxis_ = kBeam;
  else if ( hAxis == "higgs" ) hAxis_ = kHiggs;
  else throw cms::Exception("EntanglementNtupleProducer") 
    << "Invalid Configuration parameter 'hAxis' = " << hAxis << " !!\n";

  srcWeights_ = cfg.getParameter<vInputTag>("srcEvtWeights");
  for ( const edm::InputTag& srcWeight : srcWeights_ )
  {
    tokenWeights_.push_back(consumes<double>(srcWeight));
  }

  verbosity_ = cfg.getUntrackedParameter<int>("verbosity");
}

EntanglementNtupleProducer::~EntanglementNtupleProducer()
{}

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
  getChargedHadrons(const std::vector<const reco::GenParticle*>& decayProducts)
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
  getNeutralPions(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    std::vector<const reco::GenParticle*> neutralPions;
    for ( const reco::GenParticle* decayProduct : decayProducts )
    {
      if ( decayProduct->pdgId() == 111 )
      {
        neutralPions.push_back(decayProduct);
      }
    }
    return neutralPions;
  }

  std::vector<const reco::GenParticle*>
  getNeutrinos(const std::vector<const reco::GenParticle*>& decayProducts)
  {
    std::vector<const reco::GenParticle*> neutrinos;
    for ( const reco::GenParticle* decayProduct : decayProducts )
    {
      int absPdgId = std::abs(decayProduct->pdgId());
      if ( absPdgId == 12 || absPdgId == 14 || absPdgId == 16 )
      {
        neutrinos.push_back(decayProduct);
      }
    }
    return neutrinos;
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

  void
  printLorentzVector(const std::string& label,
                     const reco::Candidate::LorentzVector& p4,
                     bool cartesian = true)
  {
    std::cout << label << ":";
    if ( cartesian )
    {
      std::cout << " E = " << p4.energy() << ", Px = " << p4.px() << ", Py = " << p4.py() << ", Pz = " << p4.pz() << "\n";
    }
    else
    {
      std::cout << " pT = " << p4.pt() << ", eta = " << p4.eta() << ", phi = " << p4.phi() << ", mass = " << p4.mass() << "\n";
    }
  }

  void
  printVector(const std::string& label,
              const reco::Candidate::Vector& p3,
              bool cartesian = true)
  {
    std::cout << label << ":";
    if ( cartesian )
    {
      std::cout << " Px = " << p3.x() << ", Py = " << p3.y() << ", Pz = " << p3.z() << "\n";
    }
    else
    {
      std::cout << " pT = " << p3.r() << ", eta = " << p3.eta() << ", phi = " << p3.phi() << "\n";
    }
  }

  reco::Candidate::LorentzVector
  getP4_rf(const reco::Candidate::LorentzVector& p4,
           const ROOT::Math::Boost& boost)
  {
    // CV: boost given four-vector to restframe
    reco::Candidate::LorentzVector p4_rf = boost(p4);
    return p4_rf;
  }

  reco::Candidate::LorentzVector
  getP4_hf(const reco::Candidate::LorentzVector& p4,
           const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k)
  {
    // CV: rotate given four-vector to helicity frame
    reco::Candidate::Vector p3 = p4.Vect();
    double Pr = p3.Dot(r);
    double Pn = p3.Dot(n);
    double Pk = p3.Dot(k);
    reco::Candidate::LorentzVector p4_hf(Pr, Pn, Pk, p4.energy());
    return p4_hf;
  }

  int
  getDecayMode(const std::vector<const reco::GenParticle*>& tau_ch,
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

  reco::Candidate::LorentzVector
  getP4_ttrf_hf_trf(const reco::Candidate::LorentzVector& p4,
                    const ROOT::Math::Boost& boost_ttrf,
                    const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                    const ROOT::Math::Boost& boost_trf)
  {
    // CV: boost given four-vector to restframe of tau pair,
    //     rotate to helicity frame,
    //     and finally boost to tau restframe
    reco::Candidate::LorentzVector p4_ttrf = getP4_rf(p4, boost_ttrf);
    reco::Candidate::LorentzVector p4_hf = getP4_hf(p4_ttrf, r, n, k);
    reco::Candidate::LorentzVector p4_trf = getP4_rf(p4_hf, boost_trf);
    return p4_trf;
  }

  double
  square(double x)
  {
    return x*x;
  }
  
  double
  cube(double x)
  {
    return x*x*x;
  }

  reco::Candidate::Vector
  getPolarimetricVec_OneProng0PiZero(const reco::GenParticle* tau,
                                     const std::vector<const reco::GenParticle*> tau_ch,
                                     const std::vector<const reco::GenParticle*> tau_pi0,
                                     const std::vector<const reco::GenParticle*> tau_nu,
                                     const ROOT::Math::Boost& boost_ttrf,
                                     const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                                     const ROOT::Math::Boost& boost_trf)
  {
    //std::cout << "<getPolarimetricVec_OneProng0PiZero>:" << std::endl;
    //std::cout << " #tau_ch = " << tau_ch.size() << std::endl;
    //std::cout << " #tau_pi0 = " << tau_pi0.size() << std::endl;
    //std::cout << " #tau_nu = " << tau_nu.size() << std::endl;

    assert(tau_ch.size() == 1 && tau_pi0.size() == 0 && tau_nu.size() == 1);

    // CV: notation of four-vectors chosen according to Section 3.3 of the paper
    //       Comput.Phys.Commun. 64 (1990) 275
    const reco::GenParticle* ch = tau_ch[0];
    reco::Candidate::LorentzVector chP4 = ch->p4();

    reco::Candidate::LorentzVector Q = getP4_ttrf_hf_trf(chP4, boost_ttrf, r, n, k, boost_trf);
    const double f1 = 0.1284;
    double omega = (square(mTau) - chP4.mass2())*square(mTau);
    reco::Candidate::Vector h = -(2.*gamma_va*square(f1)*cube(mTau)/omega)*Q.Vect();
    return h.unit();
  }

  reco::Candidate::Vector
  getPolarimetricVec_OneProng1PiZero(const reco::GenParticle* tau,
                                     const std::vector<const reco::GenParticle*> tau_ch,
                                     const std::vector<const reco::GenParticle*> tau_pi0,
                                     const std::vector<const reco::GenParticle*> tau_nu,
                                     const ROOT::Math::Boost& boost_ttrf,
                                     const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                                     const ROOT::Math::Boost& boost_trf)
  {
    //std::cout << "<getPolarimetricVec_OneProng1PiZero>:" << std::endl;
    //std::cout << " #tau_ch = " << tau_ch.size() << std::endl;
    //std::cout << " #tau_pi0 = " << tau_pi0.size() << std::endl;
    //std::cout << " #tau_nu = " << tau_nu.size() << std::endl;

    assert(tau_ch.size() == 1 && tau_pi0.size() == 1 && tau_nu.size() == 1);

    // CV: notation of four-vectors chosen according to Section 3.4 of the paper
    //       Comput.Phys.Commun. 64 (1990) 275
    const reco::GenParticle* ch = tau_ch[0];
    reco::Candidate::LorentzVector chP4 = ch->p4();
    reco::Candidate::LorentzVector q1 = getP4_ttrf_hf_trf(chP4, boost_ttrf, r, n, k, boost_trf);

    const reco::GenParticle* pi0 = tau_pi0[0];
    reco::Candidate::LorentzVector pi0P4 = pi0->p4();
    reco::Candidate::LorentzVector q2 = getP4_ttrf_hf_trf(pi0P4, boost_ttrf, r, n, k, boost_trf);

    const reco::GenParticle* nu = tau_nu[0];
    reco::Candidate::LorentzVector nuP4 = nu->p4();
    reco::Candidate::LorentzVector N = getP4_ttrf_hf_trf(nuP4, boost_ttrf, r, n, k, boost_trf);

    reco::Candidate::LorentzVector tauP4 = tau->p4();
    reco::Candidate::LorentzVector P = getP4_ttrf_hf_trf(tauP4, boost_ttrf, r, n, k, boost_trf);

    reco::Candidate::LorentzVector q = q1 - q2;
    double omega = 2.*(q.Dot(N))*(q.Dot(P)) - q.mass2()*(N.Dot(P));
    // CV: term 2.*|f2|^2 appears in expression for h as well as in expression for omega
    //     and drops out
    reco::Candidate::Vector h = -(gamma_va*mTau/omega)*(2.*(q.Dot(N))*q.Vect() - q.mass2()*N.Vect());
    return h.unit();
  }

  reco::Candidate::Vector
  getPolarimetricVec_ThreeProng0PiZero(const reco::GenParticle* tau,
                                       const std::vector<const reco::GenParticle*> tau_ch,
                                       const std::vector<const reco::GenParticle*> tau_pi0,
                                       const std::vector<const reco::GenParticle*> tau_nu,
                                       const ROOT::Math::Boost& boost_ttrf,
                                       const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                                       const ROOT::Math::Boost& boost_trf)
  {
    //std::cout << "<getPolarimetricVec_ThreeProng0PiZero>:" << std::endl;
    //std::cout << " #tau_ch = " << tau_ch.size() << std::endl;
    //std::cout << " #tau_pi0 = " << tau_pi0.size() << std::endl;
    //std::cout << " #tau_nu = " << tau_nu.size() << std::endl;

    assert(tau_ch.size() == 3 && tau_pi0.size() == 0 && tau_nu.size() == 1);

    throw cms::Exception("EntanglementNtupleProducer") 
      << "Function 'getPolarimetricVec_ThreeProng0PiZero' not implemented yet !!\n";
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
    if ( genParticle.pdgId() == -15 && !tauPlus  ) tauPlus  = &genParticle;
    if ( genParticle.pdgId() == +15 && !tauMinus ) tauMinus = &genParticle;
  }

  if ( !(tauPlus && tauMinus) ) 
  {
    std::cerr << "Failed to find tau+ tau- pair --> skipping the event !!" << std::endl;
    return;
  }

  reco::Candidate::LorentzVector tauPlusP4 = tauPlus->p4();
  reco::Candidate::LorentzVector tauMinusP4 = tauMinus->p4();
  reco::Candidate::LorentzVector taupairP4 = tauPlusP4 + tauMinusP4;
  ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(taupairP4.BoostToCM());
  reco::Candidate::LorentzVector tauPlusP4_ttrf = getP4_rf(tauPlusP4, boost_ttrf);
  reco::Candidate::LorentzVector tauMinusP4_ttrf = getP4_rf(tauMinusP4, boost_ttrf);
  if ( verbosity_ >= 1 )
  {
    reco::Candidate::LorentzVector taupairP4_ttrf = getP4_rf(taupairP4, boost_ttrf);
    printLorentzVector("taupairP4_ttrf", taupairP4_ttrf);
    printLorentzVector("tauPlusP4_ttrf", tauPlusP4_ttrf);
    printLorentzVector("tauMinusP4_ttrf", tauMinusP4_ttrf);
  }

  auto k = tauMinusP4_ttrf.Vect().unit();
  if ( verbosity_ >= 1 )
  {
    printVector("k", k);
  }

  reco::Candidate::Vector h;
  if ( hAxis_ == kBeam )
  {
    double beamE = 7.e+3;
    double beamPx = 0.;
    double beamPy = 0.;
    double beamPz = std::sqrt(square(beamE) - square(mProton));
    reco::Candidate::LorentzVector beamP4(beamPx, beamPy, beamPz, beamE);
    reco::Candidate::LorentzVector beamP4_ttrf = getP4_rf(beamP4, boost_ttrf);
    h = beamP4_ttrf.Vect().unit();
  }
  else if ( hAxis_ == kHiggs )
  {
    // CV: this code does not depend on the assumption that the tau pair originates from a Higgs boson decay;
    //     it also works for tau pairs originating from Z/gamma* -> tau+ tau- decays
    const double sf = 1.01;
    double higgsPx = sf*taupairP4.px();
    double higgsPy = sf*taupairP4.py();
    double higgsPz = sf*taupairP4.pz();
    double higgsE = std::sqrt(square(higgsPx) + square(higgsPy) + square(higgsPz) + square(taupairP4.mass()));
    reco::Candidate::LorentzVector higgsP4(higgsPx, higgsPy, higgsPz, higgsE);
    reco::Candidate::LorentzVector higgsP4_ttrf = getP4_rf(higgsP4, boost_ttrf);
    h = higgsP4_ttrf.Vect().unit();
  }
  else assert(0);
  if ( verbosity_ >= 1 )
  {
    printVector("h", h);
  }

  double cosTheta = k.Dot(h);
  assert(cosTheta >= -1. && cosTheta <= +1.);
  double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  auto r = (h - k*cosTheta)*(1./sinTheta);
  if ( verbosity_ >= 1 )
  {
    printVector("r", r);
  }

  auto n = k.Cross(r);
  if ( verbosity_ >= 1 )
  {
    printVector("n", n);
  }

  if ( verbosity_ >= 1 )
  {
    std::cout << "r*n = " << r.Dot(n) << "\n";
    std::cout << "r*k = " << r.Dot(k) << "\n";
    std::cout << "n*k = " << n.Dot(k) << "\n";
  }

  reco::Candidate::LorentzVector tauPlusP4_hf = getP4_hf(tauPlusP4_ttrf, r, n, k);
  ROOT::Math::Boost boost_tprf = ROOT::Math::Boost(tauPlusP4_hf.BoostToCM());
  reco::Candidate::LorentzVector tauMinusP4_hf = getP4_hf(tauMinusP4_ttrf, r, n, k);
  ROOT::Math::Boost boost_tmrf = ROOT::Math::Boost(tauMinusP4_hf.BoostToCM());

  std::vector<const reco::GenParticle*> tauPlus_daughters;
  findDecayProducts(tauPlus, tauPlus_daughters);
  if ( verbosity_ >= 1 )
  {
    std::cout << "tau+ decay products:" << "\n";
    printGenParticles(tauPlus_daughters, false);
  }
  reco::Candidate::LorentzVector visTauPlusP4 = compVisP4(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_ch = getChargedHadrons(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_pi0 = getNeutralPions(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_nu = getNeutrinos(tauPlus_daughters);
  int tauPlus_decaymode = getDecayMode(tauPlus_ch, tauPlus_pi0, tauPlus_nu);
  reco::Candidate::LorentzVector visTauPlusP4_ttrf = getP4_rf(visTauPlusP4, boost_ttrf);
  double zPlus = visTauPlusP4_ttrf.energy()/tauPlusP4_ttrf.energy();
  if ( verbosity_ >= 1 )
  {
    std::cout << "#tauPlus_ch = " << tauPlus_ch.size() << "\n";
    std::cout << "#tauPlus_pi0 = " << tauPlus_pi0.size() << "\n";
    std::cout << "#tauPlus_nu = " << tauPlus_nu.size() << "\n";
    std::cout << "tauPlus_decaymode = " << tauPlus_decaymode << "\n";
    std::cout << "zPlus = " << zPlus << "\n";
  }

  std::vector<const reco::GenParticle*> tauMinus_daughters;
  findDecayProducts(tauMinus, tauMinus_daughters);
  if ( verbosity_ >= 1 )
  {
    std::cout << "tau- decay products:" << "\n";
    printGenParticles(tauMinus_daughters, false);
  }
  reco::Candidate::LorentzVector visTauMinusP4 = compVisP4(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_ch = getChargedHadrons(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_pi0 = getNeutralPions(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_nu = getNeutrinos(tauMinus_daughters);
  int tauMinus_decaymode = getDecayMode(tauMinus_ch, tauMinus_pi0, tauMinus_nu);
  reco::Candidate::LorentzVector visTauMinusP4_ttrf = getP4_rf(visTauMinusP4, boost_ttrf);
  double zMinus = visTauMinusP4_ttrf.energy()/tauMinusP4_ttrf.energy();
  if ( verbosity_ >= 1 )
  {
    std::cout << "#tauMinus_ch = " << tauMinus_ch.size() << "\n";
    std::cout << "#tauMinus_pi0 = " << tauMinus_pi0.size() << "\n";
    std::cout << "#tauMinus_nu = " << tauMinus_nu.size() << "\n";
    std::cout << "tauMinus_decaymode = " << tauMinus_decaymode << "\n";
    std::cout << "zMinus = " << zMinus << "\n";
  }
  
  reco::Candidate::Vector hPlus;
  reco::Candidate::Vector hMinus;
  branchType* branches = nullptr;
  if ( tauPlus_decaymode == reco::PFTau::kOneProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    if ( verbosity_ >= 1 )
    {
      reco::Candidate::LorentzVector tauPlusP4_tprf = getP4_ttrf_hf_trf(tauPlusP4, boost_ttrf, r, n, k, boost_tprf);
      printLorentzVector("tauPlusP4_tprf", tauPlusP4_tprf);
      reco::Candidate::LorentzVector chPlusP4_tprf = getP4_ttrf_hf_trf(tauPlus_ch[0]->p4(), boost_ttrf, r, n, k, boost_tprf);
      printLorentzVector("chPlusP4_tprf", chPlusP4_tprf);
      std::cout << "tauPlusP4*chPlusP4: laboratory frame = " << tauPlusP4.Dot(tauPlus_ch[0]->p4()) << "," 
                << " tau+ restframe = " <<  tauPlusP4_tprf.Dot(chPlusP4_tprf) << "\n";
      reco::Candidate::LorentzVector nubarP4_tprf = getP4_ttrf_hf_trf(tauPlus_nu[0]->p4(), boost_ttrf, r, n, k, boost_tprf);
      printLorentzVector("nubarP4_tprf", nubarP4_tprf);
      reco::Candidate::LorentzVector tauMinusP4_tmrf = getP4_ttrf_hf_trf(tauMinusP4, boost_ttrf, r, n, k, boost_tmrf);
      printLorentzVector("tauMinusP4_tmrf", tauMinusP4_tmrf);
      reco::Candidate::LorentzVector chMinusP4_tmrf = getP4_ttrf_hf_trf(tauMinus_ch[0]->p4(), boost_ttrf, r, n, k, boost_tmrf);
      printLorentzVector("chMinusP4_tmrf", chMinusP4_tmrf);
      std::cout << "tauMinusP4*chMinusP4: laboratory frame = " << tauMinusP4.Dot(tauMinus_ch[0]->p4()) << "," 
                << " tau- restframe = " <<  tauMinusP4_tmrf.Dot(chMinusP4_tmrf) << "\n";
      reco::Candidate::LorentzVector nuP4_tmrf = getP4_ttrf_hf_trf(tauMinus_nu[0]->p4(), boost_ttrf, r, n, k, boost_tmrf);
      printLorentzVector("nuP4_tmrf", nuP4_tmrf);
    }
    hPlus = getPolarimetricVec_OneProng0PiZero(tauPlus, tauPlus_ch, tauPlus_pi0, tauPlus_nu, boost_ttrf, r, n, k, boost_tprf);
    hMinus = getPolarimetricVec_OneProng0PiZero(tauMinus, tauMinus_ch, tauMinus_pi0, tauMinus_nu, boost_ttrf, r, n, k, boost_tmrf);
    branches = &branches_piPlus_piMinus_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero )
  {
    hPlus = getPolarimetricVec_OneProng0PiZero(tauPlus, tauPlus_ch, tauPlus_pi0, tauPlus_nu, boost_ttrf, r, n, k, boost_tprf);
    hMinus = getPolarimetricVec_OneProng1PiZero(tauMinus, tauMinus_ch, tauMinus_pi0, tauMinus_nu, boost_ttrf, r, n, k, boost_tmrf);
    branches = &branches_piPlus_rhoMinus_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng1PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    hPlus = getPolarimetricVec_OneProng1PiZero(tauPlus, tauPlus_ch, tauPlus_pi0, tauPlus_nu, boost_ttrf, r, n, k, boost_tprf);
    hMinus = getPolarimetricVec_OneProng0PiZero(tauMinus, tauMinus_ch, tauMinus_pi0, tauMinus_nu, boost_ttrf, r, n, k, boost_tmrf);
    branches = &branches_rhoPlus_piMinus_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kOneProng1PiZero && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero )
  {
    hPlus = getPolarimetricVec_OneProng1PiZero(tauPlus, tauPlus_ch, tauPlus_pi0, tauPlus_nu, boost_ttrf, r, n, k, boost_tprf);
    hMinus = getPolarimetricVec_OneProng1PiZero(tauMinus, tauMinus_ch, tauMinus_pi0, tauMinus_nu, boost_ttrf, r, n, k, boost_tmrf);
    branches = &branches_rhoPlus_rhoMinus_;
  }
  if ( branches )
  {
    double cosTheta = cos(tauMinusP4_ttrf.theta());
    if ( verbosity_ >= 1 )
    {
      printLorentzVector("tauPlusP4", tauPlusP4, false);
      printLorentzVector("visTauPlusP4", visTauPlusP4, false);
      printVector("hPlus", hPlus);
      printLorentzVector("tauMinusP4", tauMinusP4, false);
      printLorentzVector("visTauMinusP4", visTauMinusP4, false);
      printVector("hMinus", hMinus);
      double mTauTau = (tauPlusP4 + tauMinusP4).mass();
      double mVis = (visTauPlusP4 + visTauMinusP4).mass();
      std::cout << "mTauTau = " << mTauTau << ", mVis = " << mVis << "\n";
      std::cout << "cosTheta = " << cosTheta << "\n";
    }
    branches->fillBranches(evt, hPlus, tauPlusP4, visTauPlusP4, zPlus, hMinus, tauMinusP4, visTauMinusP4, zMinus, cosTheta, evtWeight);

    branches_had_had_.tauPlus_decaymode_ = tauPlus_decaymode;
    branches_had_had_.tauMinus_decaymode_ = tauMinus_decaymode;
    branches_had_had_.fillBranches(evt, hPlus, tauPlusP4, visTauPlusP4, zPlus, hMinus, tauMinusP4, visTauMinusP4, zMinus, cosTheta, evtWeight);
  }
}

void EntanglementNtupleProducer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(EntanglementNtupleProducer);
