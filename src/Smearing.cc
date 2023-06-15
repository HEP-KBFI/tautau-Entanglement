#include "TauAnalysis/Entanglement/interface/Smearing.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // Candidate::LorentzVector, Candidate::Point, Candidate::Vector

#include "TauAnalysis/Entanglement/interface/auxFunctions.h"              // square()
#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"                 // mChargedPion
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/KinematicParticle.h"         // KinematicParticle

#include <cmath>                                                          // std::sqrt()

Smearing::Smearing(const edm::ParameterSet& cfg)
  : resolutions_(nullptr)
{
  edm::ParameterSet cfg_smearing = cfg.getParameterSet("smearing");
  recoilSmearPx_ = cfg.getParameter<bool>("recoilSmearPx");
  recoilSmearPy_ = cfg.getParameter<bool>("recoilSmearPy");
  recoilSmearPz_ = cfg.getParameter<bool>("recoilSmearPz");
  recoilSmearE_  = cfg.getParameter<bool>("recoilSmearE");
  pvSmearXY_     = cfg.getParameter<bool>("pvSmearXY");
  pvSmearZ_      = cfg.getParameter<bool>("pvSmearZ");
  svSmearPerp_   = cfg.getParameter<bool>("svSmearPerp");
  svSmearParl_   = cfg.getParameter<bool>("svSmearParl");
  tipSmearPerp_  = cfg.getParameter<bool>("tipSmearPerp");

  edm::ParameterSet cfg_resolutions = cfg.getParameterSet("resolutions");
  resolutions_   = new Resolutions(cfg_resolutions);
}

Smearing::~Smearing()
{
  delete resolutions_;
}

KinematicEvent
Smearing::operator()(const KinematicEvent& evt)
{
  KinematicEvent smeared_evt(evt);

  const reco::Candidate::Point& pv = evt.get_pv();
  double smeared_pvX = pv.x() + rnd_.Gaus(0., resolutions_->get_pvResolutionXY());
  double smeared_pvY = pv.y() + rnd_.Gaus(0., resolutions_->get_pvResolutionXY());
  double smeared_pvZ = pv.z() + rnd_.Gaus(0., resolutions_->get_pvResolutionZ());
  smeared_evt.pv_ = reco::Candidate::Point(smeared_pvX, smeared_pvY, smeared_pvZ);

  const reco::Candidate::LorentzVector& recoilP4 = evt.get_recoilP4();
  double smeared_recoilPx = recoilP4.px()     + rnd_.Gaus(0., resolutions_->get_recoilResolutionPx());
  double smeared_recoilPy = recoilP4.py()     + rnd_.Gaus(0., resolutions_->get_recoilResolutionPy());
  double smeared_recoilPz = recoilP4.pz()     + rnd_.Gaus(0., resolutions_->get_recoilResolutionPz());
  double smeared_recoilE  = recoilP4.energy() + rnd_.Gaus(0., resolutions_->get_recoilResolutionE());
  smeared_evt.recoilP4_ = reco::Candidate::LorentzVector(smeared_recoilPx, smeared_recoilPy, smeared_recoilPz, smeared_recoilE);

  smeared_evt.tauPlusP4_ = reco::Candidate::LorentzVector(0., 0., 0., 0.);
  smeared_evt.tauPlusP4_isValid_ = false;
  const std::vector<KinematicParticle>& daughtersTauPlus = evt.get_daughtersTauPlus();
  smeared_evt.daughtersTauPlus_.clear();
  for ( const KinematicParticle& daughterTauPlus : daughtersTauPlus )
  {
    smeared_evt.daughtersTauPlus_.push_back(smear_daughter(daughterTauPlus));
  }
  smeared_evt.tipPCATauPlus_ = smear_tipPCA(daughtersTauPlus, evt.get_tipPCATauPlus());
  if ( evt.get_svTauPlus_isValid() )
  {
    //const reco::Candidate::Point& svTauPlus = evt.get_svTauPlus();
    //double Px = svTauPlus.x() - pv.x();
    //double Py = svTauPlus.y() - pv.y();
    //double Pz = svTauPlus.z() - pv.z();
    //double E  = std::sqrt(square(Px) + square(Py) + square(Pz) + square(mChargedPion));
    //reco::Candidate::LorentzVector p4(Px, Py, Pz, E);
    //smeared_evt.svTauPlus_ = smear_sv(p4, svTauPlus);
    smeared_evt.svTauPlus_ = smear_sv(evt.get_visTauPlusP4(), evt.get_svTauPlus());
  }

  smeared_evt.tauMinusP4_ = reco::Candidate::LorentzVector(0., 0., 0., 0.);
  smeared_evt.tauMinusP4_isValid_ = false;
  const std::vector<KinematicParticle>& daughtersTauMinus = evt.get_daughtersTauMinus();
  smeared_evt.daughtersTauMinus_.clear();
  for ( const KinematicParticle& daughterTauMinus : daughtersTauMinus )
  {
    smeared_evt.daughtersTauMinus_.push_back(smear_daughter(daughterTauMinus));
  }
  smeared_evt.tipPCATauMinus_ = smear_tipPCA(daughtersTauMinus, evt.get_tipPCATauMinus());
  if ( evt.get_svTauMinus_isValid() )
  {
    //const reco::Candidate::Point& svTauMinus = evt.get_svTauMinus();
    //double Px = svTauMinus.x() - pv.x();
    //double Py = svTauMinus.y() - pv.y();
    //double Pz = svTauMinus.z() - pv.z();
    //double E  = std::sqrt(square(Px) + square(Py) + square(Pz) + square(mChargedPion));
    //reco::Candidate::LorentzVector p4(Px, Py, Pz, E);
    //smeared_evt.svTauMinus_ = smear_sv(p4, svTauMinus);
    smeared_evt.svTauMinus_ = smear_sv(evt.get_visTauMinusP4(), evt.get_svTauMinus());
  }

  return smeared_evt;
}

KinematicParticle
Smearing::smear_daughter(const KinematicParticle& daughter)
{
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(daughter.get_p4(), nullptr, nullptr, kBeam, r, n, k);
  double dr = rnd_.Gaus(0., resolutions_->get_svResolutionPerp());
  double dn = rnd_.Gaus(0., resolutions_->get_svResolutionPerp());
  double dk = rnd_.Gaus(0., resolutions_->get_svResolutionParl());
  const reco::Candidate::Point& vertex = daughter.get_vertex();
  double smeared_vertexX = vertex.x() + dr*r.x() + dn*n.x() + dk*k.x();
  double smeared_vertexY = vertex.y() + dr*r.y() + dn*n.y() + dk*k.y();
  double smeared_vertexZ = vertex.z() + dr*r.z() + dn*n.z() + dk*k.z();
  KinematicParticle smeared_daughter(daughter);
  smeared_daughter.vertex_ = reco::Candidate::Point(smeared_vertexX, smeared_vertexY, smeared_vertexZ);
  return smeared_daughter;
}

namespace
{
  const KinematicParticle*
  get_leadTrack(const std::vector<KinematicParticle>& daughters)
  {
    const KinematicParticle* leadTrack = nullptr;
    double max_pt = -1.;
    for ( const KinematicParticle& daughter : daughters )
    {
      if ( std::fabs(daughter.get_charge()) > 0.5 && daughter.get_p4().pt() > max_pt )
      {
        leadTrack = &daughter;
        max_pt = daughter.get_p4().pt();
      }
    }
    return leadTrack;
  }
}

reco::Candidate::Point
Smearing::smear_tipPCA(const std::vector<KinematicParticle>& daughters, const reco::Candidate::Point& tipPCA)
{
  const KinematicParticle* leadTrack = get_leadTrack(daughters);
  if ( !leadTrack )
  {
    std::cerr << "WARNING: Failed to find leading track of tau --> returning null vector !!" << std::endl;
    return reco::Candidate::Point(0., 0., 0.);
  }
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(leadTrack->get_p4(), nullptr, nullptr, kBeam, r, n, k);
  double dr = rnd_.Gaus(0., resolutions_->get_tipResolutionPerp());
  double dn = rnd_.Gaus(0., resolutions_->get_tipResolutionPerp());
  double smeared_pcaX = tipPCA.x() + dr*r.x() + dn*n.x();
  double smeared_pcaY = tipPCA.y() + dr*r.y() + dn*n.y();
  double smeared_pcaZ = tipPCA.z() + dr*r.z() + dn*n.z();
  reco::Candidate::Point smeared_tipPCA(smeared_pcaX, smeared_pcaY, smeared_pcaZ);
  return smeared_tipPCA;
}

reco::Candidate::Point
Smearing::smear_sv(const reco::Candidate::LorentzVector& p4, const reco::Candidate::Point& sv)
{
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(p4, nullptr, nullptr, kBeam, r, n, k);
  double dr = rnd_.Gaus(0., resolutions_->get_svResolutionPerp());
  double dn = rnd_.Gaus(0., resolutions_->get_svResolutionPerp());
  double dk = rnd_.Gaus(0., resolutions_->get_svResolutionParl());
  double smeared_svX = sv.x() + dr*r.x() + dn*n.x() + dk*k.x();
  double smeared_svY = sv.y() + dr*r.y() + dn*n.y() + dk*k.y();
  double smeared_svZ = sv.z() + dr*r.z() + dn*n.z() + dk*k.z();
  reco::Candidate::Point smeared_sv(smeared_svX, smeared_svY, smeared_svZ);
  return smeared_sv;
}
