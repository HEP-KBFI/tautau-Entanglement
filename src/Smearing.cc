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
  applySmearing_recoilPx_ = cfg_smearing.getParameter<bool>("applySmearing_recoilPx");
  applySmearing_recoilPy_ = cfg_smearing.getParameter<bool>("applySmearing_recoilPy");
  applySmearing_recoilPz_ = cfg_smearing.getParameter<bool>("applySmearing_recoilPz");
  applySmearing_recoilE_  = cfg_smearing.getParameter<bool>("applySmearing_recoilE");
  applySmearing_pvXY_     = cfg_smearing.getParameter<bool>("applySmearing_pvXY");
  applySmearing_pvZ_      = cfg_smearing.getParameter<bool>("applySmearing_pvZ");
  applySmearing_svPerp_   = cfg_smearing.getParameter<bool>("applySmearing_svPerp");
  applySmearing_svParl_   = cfg_smearing.getParameter<bool>("applySmearing_svParl");
  applySmearing_tipPerp_  = cfg_smearing.getParameter<bool>("applySmearing_tipPerp");

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
  double smeared_pvX = pv.x();
  double smeared_pvY = pv.y();
  double smeared_pvZ = pv.z();
  if ( applySmearing_pvXY_ )
  {
    smeared_pvX += rnd_.Gaus(0., resolutions_->get_pvResolutionXY());
    smeared_pvY += rnd_.Gaus(0., resolutions_->get_pvResolutionXY());
  }
  if ( applySmearing_pvZ_ )
  {
    smeared_pvZ += rnd_.Gaus(0., resolutions_->get_pvResolutionZ());
  }
  smeared_evt.pv_ = reco::Candidate::Point(smeared_pvX, smeared_pvY, smeared_pvZ);

  const reco::Candidate::LorentzVector& recoilP4 = evt.get_recoilP4();
  double smeared_recoilPx = recoilP4.px();
  double smeared_recoilPy = recoilP4.py();
  double smeared_recoilPz = recoilP4.pz();
  double smeared_recoilE  = recoilP4.energy();
  if ( applySmearing_recoilPx_ )
  {
    smeared_recoilPx += rnd_.Gaus(0., resolutions_->get_recoilResolutionPx());
  }
  if ( applySmearing_recoilPy_ )
  {
    smeared_recoilPy += rnd_.Gaus(0., resolutions_->get_recoilResolutionPy());
  }
  if ( applySmearing_recoilPz_ )
  {
    smeared_recoilPz += rnd_.Gaus(0., resolutions_->get_recoilResolutionPz());
  }
  if ( applySmearing_recoilE_ )
  {
    smeared_recoilE  += rnd_.Gaus(0., resolutions_->get_recoilResolutionE());
  }
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
  const reco::Candidate::Point& sv = daughter.get_vertex();
  double smeared_svX = sv.x();
  double smeared_svY = sv.y();
  double smeared_svZ = sv.z();
  if ( applySmearing_svPerp_ )
  {
    smeared_svX += dr*r.x() + dn*n.x();
    smeared_svY += dr*r.y() + dn*n.y();
    smeared_svZ += dr*r.z() + dn*n.z();
  }
  if ( applySmearing_svParl_ )
  {
    smeared_svX += dk*k.x();
    smeared_svY += dk*k.y();
    smeared_svZ += dk*k.z();
  }
  KinematicParticle smeared_daughter(daughter);
  smeared_daughter.vertex_ = reco::Candidate::Point(smeared_svX, smeared_svY, smeared_svZ);
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
  double smeared_pcaX = tipPCA.x();
  double smeared_pcaY = tipPCA.y();
  double smeared_pcaZ = tipPCA.z();
  if ( applySmearing_tipPerp_ )
  {
    smeared_pcaX += dr*r.x() + dn*n.x();
    smeared_pcaY += dr*r.y() + dn*n.y();
    smeared_pcaZ += dr*r.z() + dn*n.z();
  }
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
  double smeared_svX = sv.x();
  double smeared_svY = sv.y();
  double smeared_svZ = sv.z();
  if ( applySmearing_svPerp_ )
  {
    smeared_svX += dr*r.x() + dn*n.x();
    smeared_svY += dr*r.y() + dn*n.y();
    smeared_svZ += dr*r.z() + dn*n.z();
  }
  if ( applySmearing_svParl_ )
  {
    smeared_svX += dk*k.x();
    smeared_svY += dk*k.y();
    smeared_svZ += dk*k.z();
  }
  reco::Candidate::Point smeared_sv(smeared_svX, smeared_svY, smeared_svZ);
  return smeared_sv;
}
