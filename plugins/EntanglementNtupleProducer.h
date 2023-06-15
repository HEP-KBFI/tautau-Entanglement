#ifndef TauAnalysis_Entanglement_EntanglementNtupleProducer_h
#define TauAnalysis_Entanglement_EntanglementNtupleProducer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"                         // edm::EDAnalyzer
#include "FWCore/Framework/interface/Event.h"                              // edm::Event
#include "FWCore/Framework/interface/EventSetup.h"                         // edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"                    // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h"                           // edm::InputTag

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"              // reco::GenParticle
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"           // reco::GenParticleCollection

#include "TauAnalysis/Entanglement/interface/KinematicFit.h"               // KinematicFitOneProng0Pi0
#include "TauAnalysis/Entanglement/interface/Smearing.h"                   // Smearing
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerOneProng0Pi0.h"   // SpinAnalyzerOneProng0Pi0
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerOneProng1Pi0.h"   // SpinAnalyzerOneProng1Pi0
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerThreeProng0Pi0.h" // SpinAnalyzerThreeProng0Pi0

#include <TTree.h>                                                         // TTree
#include <Rtypes.h>                                                        // Float_t, Int_t

#include <vector>                                                          // std::vector
#include <string>                                                          // std::string

class EntanglementNtupleProducer : public edm::EDAnalyzer 
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

  enum { kGen, kRec };
  int mode_;

  Smearing smearing_;

  KinematicFit kineFit_;

  SpinAnalyzerOneProng0Pi0 spinAnalyzerOneProng0Pi0_;
  SpinAnalyzerOneProng1Pi0 spinAnalyzerOneProng1Pi0_;
  SpinAnalyzerThreeProng0Pi0 spinAnalyzerThreeProng0Pi0_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights_;
  std::vector<edm::EDGetTokenT<double>> tokenWeights_;

  class branchType
  {
   public:
    branchType(int tauPlus_decaymode, int tauMinus_decaymode)
      : ntuple_(nullptr)
      , run_(0)
      , lumi_(0)
      , event_(0)
      , hPlus_r_(0.)
      , hPlus_n_(0.)
      , hPlus_k_(0.)
      , tauPlus_pt_(0.)
      , tauPlus_eta_(0.)
      , tauPlus_phi_(0.)
      , tauPlus_mass_(0.)
      , tauPlus_decaymode_(tauPlus_decaymode)
      , tauPlus_nChargedKaons_(0)
      , tauPlus_nNeutralKaons_(0)
      , tauPlus_nPhotons_(0)
      , tauPlus_sumPhotonEn_(0.)
      , visTauPlus_pt_(0.)
      , visTauPlus_eta_(0.)
      , visTauPlus_phi_(0.)
      , visTauPlus_mass_(0.)
      , zPlus_(0.)
      , hMinus_r_(0.)
      , hMinus_n_(0.)
      , hMinus_k_(0.)
      , tauMinus_pt_(0.)
      , tauMinus_eta_(0.)
      , tauMinus_phi_(0.)
      , tauMinus_mass_(0.)
      , tauMinus_decaymode_(tauMinus_decaymode)
      , tauMinus_nChargedKaons_(0)
      , tauMinus_nNeutralKaons_(0)
      , tauMinus_nPhotons_(0)
      , tauMinus_sumPhotonEn_(0.)
      , visTauMinus_pt_(0.)
      , visTauMinus_eta_(0.)
      , visTauMinus_phi_(0.)
      , visTauMinus_mass_(0.)
      , zMinus_(0.)
      , mTauTau_(0.)
      , mVis_(0.)
      , cosTheta_(0.)
      , evtWeight_(1.)
    {}
    ~branchType()
    {}

    void
    initBranches(TTree* ntuple)
    {
      ntuple_ = ntuple;

      ntuple->Branch("run", &run_, "run/i");
      ntuple->Branch("lumi", &lumi_, "lumi/i");
      ntuple->Branch("event", &event_, "event/l");

      ntuple->Branch("hPlus_r", &hPlus_r_, "hPlus_r/F");
      ntuple->Branch("hPlus_n", &hPlus_n_, "hPlus_n/F");
      ntuple->Branch("hPlus_k", &hPlus_k_, "hPlus_k/F");
      ntuple->Branch("tauPlus_pt", &tauPlus_pt_, "tauPlus_pt/F");
      ntuple->Branch("tauPlus_eta", &tauPlus_eta_, "tauPlus_eta/F");
      ntuple->Branch("tauPlus_phi", &tauPlus_phi_, "tauPlus_phi/F");
      ntuple->Branch("tauPlus_mass", &tauPlus_mass_, "tauPlus_mass/F");
      ntuple->Branch("tauPlus_decaymode", &tauPlus_decaymode_, "tauPlus_decaymode/I");
      ntuple->Branch("tauPlus_nChargedKaons", &tauPlus_nChargedKaons_, "tauPlus_nChargedKaons/I");
      ntuple->Branch("tauPlus_nNeutralKaons", &tauPlus_nNeutralKaons_, "tauPlus_nNeutralKaons/I");
      ntuple->Branch("tauPlus_nPhotons", &tauPlus_nPhotons_, "tauPlus_nPhotons/I");
      ntuple->Branch("tauPlus_sumPhotonEn", &tauPlus_sumPhotonEn_, "tauPlus_sumPhotonEn/F");
      ntuple->Branch("visTauPlus_pt", &visTauPlus_pt_, "visTauPlus_pt/F");
      ntuple->Branch("visTauPlus_eta", &visTauPlus_eta_, "visTauPlus_eta/F");
      ntuple->Branch("visTauPlus_phi", &visTauPlus_phi_, "visTauPlus_phi/F");
      ntuple->Branch("visTauPlus_mass", &visTauPlus_mass_, "visTauPlus_mass/F");
      ntuple->Branch("zPlus", &zPlus_, "zPlus/F");

      ntuple->Branch("hMinus_r", &hMinus_r_, "hMinus_r/F");
      ntuple->Branch("hMinus_n", &hMinus_n_, "hMinus_n/F");
      ntuple->Branch("hMinus_k", &hMinus_k_, "hMinus_k/F");
      ntuple->Branch("tauMinus_pt", &tauMinus_pt_, "tauMinus_pt/F");
      ntuple->Branch("tauMinus_eta", &tauMinus_eta_, "tauMinus_eta/F");
      ntuple->Branch("tauMinus_phi", &tauMinus_phi_, "tauMinus_phi/F");
      ntuple->Branch("tauMinus_mass", &tauMinus_mass_, "tauMinus_mass/F");
      ntuple->Branch("tauMinus_decaymode", &tauMinus_decaymode_, "tauMinus_decaymode/I");
      ntuple->Branch("tauMinus_nChargedKaons", &tauMinus_nChargedKaons_, "tauminus_nChargedKaons/I");
      ntuple->Branch("tauMinus_nNeutralKaons", &tauMinus_nNeutralKaons_, "tauMinus_nNeutralKaons/I");
      ntuple->Branch("tauMinus_nPhotons", &tauMinus_nPhotons_, "tauMinus_nPhotons/I");
      ntuple->Branch("tauMinus_sumPhotonEn", &tauMinus_sumPhotonEn_, "tauMinus_sumPhotonEn/F");
      ntuple->Branch("visTauMinus_pt", &visTauMinus_pt_, "visTauMinus_pt/F");
      ntuple->Branch("visTauMinus_eta", &visTauMinus_eta_, "visTauMinus_eta/F");
      ntuple->Branch("visTauMinus_phi", &visTauMinus_phi_, "visTauMinus_phi/F");
      ntuple->Branch("visTauMinus_mass", &visTauMinus_mass_, "visTauMinus_mass/F");
      ntuple->Branch("zMinus", &zMinus_, "zMinus/F");

      ntuple->Branch("mTauTau", &mTauTau_, "mTauTau/F");
      ntuple->Branch("mVis", &mVis_, "mVis/F");
      ntuple->Branch("cosTheta", &cosTheta_, "cosTheta/F");

      ntuple->Branch("evtWeight", &evtWeight_, "evtWeight/F");
    }

    void
    fillBranches(const edm::Event & event,
                 const reco::Candidate::Vector& hPlus,
                 const reco::Candidate::LorentzVector& tauPlusP4,
                 Int_t tauPlus_nChargedKaons, Int_t tauPlus_nNeutralKaons, Int_t tauPlus_nPhotons, double tauPlus_sumPhotonEn,
                 const reco::Candidate::LorentzVector& visTauPlusP4,
                 double zPlus,
                 const reco::Candidate::Vector& hMinus,
                 const reco::Candidate::LorentzVector& tauMinusP4,
                 Int_t tauMinus_nChargedKaons, Int_t tauMinus_nNeutralKaons, Int_t tauMinus_nPhotons, double tauMinus_sumPhotonEn,
                 const reco::Candidate::LorentzVector& visTauMinusP4,
                 double zMinus,
                 double cosTheta,
                 double evtWeight)
    {
      assert(ntuple_);

      run_ = event.id().run();
      lumi_ = event.id().luminosityBlock();
      event_ = event.id().event();

      hPlus_r_ = hPlus.x();
      hPlus_n_ = hPlus.y();
      hPlus_k_ = hPlus.z();
      tauPlus_pt_ = tauPlusP4.pt();
      tauPlus_eta_ = tauPlusP4.eta();
      tauPlus_phi_ = tauPlusP4.phi();
      tauPlus_mass_ = tauPlusP4.mass();
      tauPlus_nChargedKaons_ = tauPlus_nChargedKaons;
      tauPlus_nNeutralKaons_ = tauPlus_nNeutralKaons;
      tauPlus_nPhotons_ = tauPlus_nPhotons;
      tauPlus_sumPhotonEn_ = tauPlus_sumPhotonEn;
      visTauPlus_pt_ = visTauPlusP4.pt();
      visTauPlus_eta_ = visTauPlusP4.eta();
      visTauPlus_phi_ = visTauPlusP4.phi();
      visTauPlus_mass_ = visTauPlusP4.mass();
      zPlus_ = zPlus;

      hMinus_r_ = hMinus.x();
      hMinus_n_ = hMinus.y();
      hMinus_k_ = hMinus.z();
      tauMinus_pt_ = tauMinusP4.pt();
      tauMinus_eta_ = tauMinusP4.eta();
      tauMinus_phi_ = tauMinusP4.phi();
      tauMinus_mass_ = tauMinusP4.mass();
      tauMinus_nChargedKaons_ = tauMinus_nChargedKaons;
      tauMinus_nNeutralKaons_ = tauMinus_nNeutralKaons;
      tauMinus_nPhotons_ = tauMinus_nPhotons;
      tauMinus_sumPhotonEn_ = tauMinus_sumPhotonEn;
      visTauMinus_pt_ = visTauMinusP4.pt();
      visTauMinus_eta_ = visTauMinusP4.eta();
      visTauMinus_phi_ = visTauMinusP4.phi();
      visTauMinus_mass_ = visTauMinusP4.mass();
      zMinus_ = zMinus;

      mTauTau_ = (tauPlusP4 + tauMinusP4).mass();
      mVis_ = (visTauPlusP4 + visTauMinusP4).mass();
      cosTheta_ = cosTheta;

      evtWeight_ = evtWeight;

      ntuple_->Fill();
    }
    
    TTree* ntuple_;

    UInt_t run_;                   // run number
    UInt_t lumi_;                  // luminosity-section number
    ULong64_t event_;              // event number    

    Float_t hPlus_r_;              // r component (in helicity frame) of polarimetric vector of tau+
    Float_t hPlus_n_;              // n component (in helicity frame) of polarimetric vector of tau+
    Float_t hPlus_k_;              // k component (in helicity frame) of polarimetric vector of tau+
    Float_t tauPlus_pt_;           // transverse momentum (in laboratory frame) of tau+
    Float_t tauPlus_eta_;          // pseudo-rapidity (in laboratory frame) of tau+
    Float_t tauPlus_phi_;          // azimuthal angle (in laboratory frame) of tau+
    Float_t tauPlus_mass_;         // mass component of tau+ four-vector
    Int_t tauPlus_decaymode_;      // tau+ decay mode
    Int_t tauPlus_nChargedKaons_;  // number of charged Kaons produced in decay of tau-
    Int_t tauPlus_nNeutralKaons_;  // number of neutral Kaons produced in decay of tau-
    Int_t tauPlus_nPhotons_;       // number of photons radiated from tau- or tau- decay products
    Float_t tauPlus_sumPhotonEn_;  // energy sum of photons radiated from tau- or tau- decay products
    Float_t visTauPlus_pt_;        // transverse momentum (in laboratory frame) of visible decay products of tau+
    Float_t visTauPlus_eta_;       // pseudo-rapidity (in laboratory frame) of visible decay products of tau+
    Float_t visTauPlus_phi_;       // azimuthal angle (in laboratory frame) of visible decay products of tau+
    Float_t visTauPlus_mass_;      // mass of visible decay products of tau+
    Float_t zPlus_;                // fraction of tau+ energy (in tau-pair restframe) carried by visible decay products of tau+

    Float_t hMinus_r_;             // r component (in helicity frame) of polarimetric vector of tau-
    Float_t hMinus_n_;             // n component (in helicity frame) of polarimetric vector of tau-
    Float_t hMinus_k_;             // k component (in helicity frame) of polarimetric vector of tau-
    Float_t tauMinus_pt_;          // transverse momentum (in laboratory frame) of tau-
    Float_t tauMinus_eta_;         // pseudo-rapidity (in laboratory frame) of tau-
    Float_t tauMinus_phi_;         // azimuthal angle (in laboratory frame) of tau-
    Float_t tauMinus_mass_;        // mass component of tau- four-vector
    Int_t tauMinus_decaymode_;     // tau- decay mode
    Int_t tauMinus_nChargedKaons_; // number of charged Kaons produced in decay of tau-
    Int_t tauMinus_nNeutralKaons_; // number of neutral Kaons produced in decay of tau-
    Int_t tauMinus_nPhotons_;      // number of photons radiated from tau- or tau- decay products
    Float_t tauMinus_sumPhotonEn_; // energy sum of photons radiated from tau- or tau- decay products
    Float_t visTauMinus_pt_;       // transverse momentum (in laboratory frame) of visible decay products of tau-
    Float_t visTauMinus_eta_;      // pseudo-rapidity (in laboratory frame) of visible decay products of tau-
    Float_t visTauMinus_phi_;      // azimuthal angle (in laboratory frame) of visible decay products of tau-
    Float_t visTauMinus_mass_;     // mass of visible decay products of tau-
    Float_t zMinus_;               // fraction of tau- energy (in tau-pair restframe) carried by visible decay products of tau-

    Float_t mTauTau_;              // mass of tau pair
    Float_t mVis_;                 // mass of visible decay products of tau pair
    Float_t cosTheta_;             // polar angle of tau- in tau-pair restframe

    Float_t evtWeight_;            // event weight of Monte Carlo generator
  };

  TTree* ntuple_piPlus_piMinus_;
  branchType branches_piPlus_piMinus_;

  TTree* ntuple_piPlus_rhoMinus_;
  branchType branches_piPlus_rhoMinus_;

  TTree* ntuple_rhoPlus_piMinus_;
  branchType branches_rhoPlus_piMinus_;

  TTree* ntuple_rhoPlus_rhoMinus_;
  branchType branches_rhoPlus_rhoMinus_;

  TTree* ntuple_had_had_;
  branchType branches_had_had_;

  int verbosity_;
};

#endif // TauAnalysis_Entanglement_EntanglementNtupleProducer_h
