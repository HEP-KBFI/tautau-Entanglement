#include "TauAnalysis/Entanglement/interface/Smearing.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // Candidate::LorentzVector, Candidate::Point, Candidate::Vector

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"                 // mChargedPion
#include "TauAnalysis/Entanglement/interface/get_leadTrack.h"             // get_leadTrack()
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/KinematicParticle.h"         // KinematicParticle
#include "TauAnalysis/Entanglement/interface/square.h"                    // square()

#include <TMath.h>                                                        // TMath::Pi()

#include <cmath>                                                          // std::sqrt()

Smearing::Smearing(const edm::ParameterSet& cfg)
  : resolutions_(nullptr)
{
  edm::ParameterSet cfg_smearing = cfg.getParameterSet("smearing");
  applySmearing_recoil_px_     = cfg_smearing.getParameter<bool>("applySmearing_recoil_px");
  applySmearing_recoil_py_     = cfg_smearing.getParameter<bool>("applySmearing_recoil_py");
  applySmearing_recoil_pz_     = cfg_smearing.getParameter<bool>("applySmearing_recoil_pz");
  applySmearing_recoil_energy_ = cfg_smearing.getParameter<bool>("applySmearing_recoil_energy");

  applySmearing_pv_xy_         = cfg_smearing.getParameter<bool>("applySmearing_pv_xy");
  applySmearing_pv_z_          = cfg_smearing.getParameter<bool>("applySmearing_pv_z");

  applySmearing_track_pt_      = cfg_smearing.getParameter<bool>("applySmearing_track_pt");
  applySmearing_track_theta_   = cfg_smearing.getParameter<bool>("applySmearing_track_theta");
  applySmearing_track_phi_     = cfg_smearing.getParameter<bool>("applySmearing_track_phi");

  applySmearing_ecal_energy_   = cfg_smearing.getParameter<bool>("applySmearing_ecal_energy");
  applySmearing_ecal_theta_    = cfg_smearing.getParameter<bool>("applySmearing_ecal_theta");
  applySmearing_ecal_phi_      = cfg_smearing.getParameter<bool>("applySmearing_ecal_phi");

  applySmearing_sv_perp_       = cfg_smearing.getParameter<bool>("applySmearing_sv_perp");
  applySmearing_sv_parl_       = cfg_smearing.getParameter<bool>("applySmearing_sv_parl");

  applySmearing_tip_perp_      = cfg_smearing.getParameter<bool>("applySmearing_tip_perp");

  edm::ParameterSet cfg_resolutions = cfg.getParameterSet("resolutions");
  resolutions_ = new Resolutions(cfg_resolutions);
}

Smearing::~Smearing()
{
  delete resolutions_;
}

KinematicEvent
Smearing::operator()(const KinematicEvent& kineEvt)
{
  KinematicEvent kineEvt_smeared(kineEvt);

  kineEvt_smeared.pv_ = smear_pv(kineEvt.get_pv());

  kineEvt_smeared.recoilP4_ = smear_recoil_p4(kineEvt.get_recoilP4());

  kineEvt_smeared.tauPlusP4_ = reco::Candidate::LorentzVector(0.,0.,0.,0.);
  kineEvt_smeared.tauPlusP4_isValid_ = false;
  if ( kineEvt.get_svTauPlus_isValid() )
  {
    kineEvt_smeared.svTauPlus_ = smear_sv(kineEvt.get_visTauPlusP4(), kineEvt.get_svTauPlus());
  }
  const std::vector<KinematicParticle>& daughtersTauPlus = kineEvt.get_daughtersTauPlus();
  kineEvt_smeared.daughtersTauPlus_.clear();
  for ( const KinematicParticle& daughter : daughtersTauPlus )
  {
    KinematicParticle daughter_smeared = daughter;
    daughter_smeared.p4_ = smear_daughter_p4(daughter);
    daughter_smeared.vertex_ = kineEvt_smeared.svTauPlus_;
    kineEvt_smeared.daughtersTauPlus_.push_back(daughter_smeared);
  }
  kineEvt_smeared.tipPCATauPlus_ = smear_tipPCA(daughtersTauPlus, kineEvt.get_tipPCATauPlus());

  kineEvt_smeared.tauMinusP4_ = reco::Candidate::LorentzVector(0.,0.,0.,0.);
  kineEvt_smeared.tauMinusP4_isValid_ = false;
  if ( kineEvt.get_svTauMinus_isValid() )
  {
    kineEvt_smeared.svTauMinus_ = smear_sv(kineEvt.get_visTauMinusP4(), kineEvt.get_svTauMinus());
  }
  const std::vector<KinematicParticle>& daughtersTauMinus = kineEvt.get_daughtersTauMinus();
  kineEvt_smeared.daughtersTauMinus_.clear();
  for ( const KinematicParticle& daughter : daughtersTauMinus )
  {
    KinematicParticle daughter_smeared = daughter;
    daughter_smeared.p4_ = smear_daughter_p4(daughter);
    daughter_smeared.vertex_ = kineEvt_smeared.svTauPlus_;
    kineEvt_smeared.daughtersTauMinus_.push_back(daughter_smeared);
  }
  kineEvt_smeared.tipPCATauMinus_ = smear_tipPCA(daughtersTauMinus, kineEvt.get_tipPCATauMinus());

  return kineEvt_smeared;
}

reco::Candidate::Point
Smearing::smear_pv(const reco::Candidate::Point& pv)
{
  double smeared_pvX = pv.x();
  double smeared_pvY = pv.y();
  double smeared_pvZ = pv.z();
  if ( applySmearing_pv_xy_ )
  {
    smeared_pvX += rnd_.Gaus(0., resolutions_->get_pvResolution_xy());
    smeared_pvY += rnd_.Gaus(0., resolutions_->get_pvResolution_xy());
  }
  if ( applySmearing_pv_z_ )
  {
    smeared_pvZ += rnd_.Gaus(0., resolutions_->get_pvResolution_z());
  }
  return reco::Candidate::Point(smeared_pvX, smeared_pvY, smeared_pvZ);
}

reco::Candidate::LorentzVector
Smearing::smear_recoil_p4(const reco::Candidate::LorentzVector& recoilP4)
{
  double smeared_recoilPx = recoilP4.px();
  double smeared_recoilPy = recoilP4.py();
  double smeared_recoilPz = recoilP4.pz();
  double smeared_recoilE  = recoilP4.energy();
  if ( applySmearing_recoil_px_ )
  {
    smeared_recoilPx += rnd_.Gaus(0., resolutions_->get_recoilResolution_px());
  }
  if ( applySmearing_recoil_py_ )
  {
    smeared_recoilPy += rnd_.Gaus(0., resolutions_->get_recoilResolution_py());
  }
  if ( applySmearing_recoil_pz_ )
  {
    smeared_recoilPz += rnd_.Gaus(0., resolutions_->get_recoilResolution_pz());
  }
  if ( applySmearing_recoil_energy_ )
  {
    smeared_recoilE  += rnd_.Gaus(0., resolutions_->get_recoilResolution_energy());
  }
  return reco::Candidate::LorentzVector(smeared_recoilPx, smeared_recoilPy, smeared_recoilPz, smeared_recoilE);
}

reco::Candidate::Point
Smearing::smear_sv(const reco::Candidate::LorentzVector& p4, const reco::Candidate::Point& sv)
{
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(p4, nullptr, nullptr, kBeam, r, n, k);
  double dr = rnd_.Gaus(0., resolutions_->get_svResolution_perp());
  double dn = rnd_.Gaus(0., resolutions_->get_svResolution_perp());
  double dk = rnd_.Gaus(0., resolutions_->get_svResolution_parl());
  double smeared_svX = sv.x();
  double smeared_svY = sv.y();
  double smeared_svZ = sv.z();
  if ( applySmearing_sv_perp_ )
  {
    smeared_svX += dr*r.x() + dn*n.x();
    smeared_svY += dr*r.y() + dn*n.y();
    smeared_svZ += dr*r.z() + dn*n.z();
  }
  if ( applySmearing_sv_parl_ )
  {
    smeared_svX += dk*k.x();
    smeared_svY += dk*k.y();
    smeared_svZ += dk*k.z();
  }
  reco::Candidate::Point smeared_sv(smeared_svX, smeared_svY, smeared_svZ);
  return smeared_sv;
}

reco::Candidate::LorentzVector
Smearing::smear_daughter_p4(const KinematicParticle& daughter)
{
  double sigma_pt    = 0.;
  double sigma_theta = 0.;
  double sigma_phi   = 0.;      
  if ( std::fabs(daughter.get_charge()) > 0.5 )
  {
    sigma_pt         = get_trackResolution_pt(daughter.get_p4(), *resolutions_);
    sigma_theta      = TMath::Pi()*resolutions_->get_trackResolution_theta();
    sigma_phi        = TMath::Pi()*resolutions_->get_trackResolution_phi();
  }
  else if ( daughter.get_pdgId() == 111 )
  {
    sigma_pt         = get_ecalResolution_pt(daughter.get_p4(), *resolutions_);
    sigma_theta      = TMath::Pi()*resolutions_->get_ecalResolution_theta();
    sigma_phi        = TMath::Pi()*resolutions_->get_ecalResolution_phi();
  }
  const reco::Candidate::LorentzVector& daughterP4 = daughter.get_p4();
  double smeared_daughterPt    = rnd_.Gaus(daughterP4.pt(), sigma_pt);
  double smeared_daughterTheta = rnd_.Gaus(daughterP4.theta(), sigma_theta);
  double smeared_daughterPhi   = rnd_.Gaus(daughterP4.phi(), sigma_phi);
  double smeared_daughterPx    = smeared_daughterPt*cos(smeared_daughterPhi);
  double smeared_daughterPy    = smeared_daughterPt*sin(smeared_daughterPhi);
  double smeared_daughterPz    = smeared_daughterPt/tan(smeared_daughterTheta);
  double smeared_daughterE     = std::sqrt(square(smeared_daughterPx) + square(smeared_daughterPy) + square(smeared_daughterPz) + square(daughterP4.mass()));
  return reco::Candidate::LorentzVector(smeared_daughterPx, smeared_daughterPy, smeared_daughterPz, smeared_daughterE);
}

reco::Candidate::Point
Smearing::smear_tipPCA(const std::vector<KinematicParticle>& daughters, const reco::Candidate::Point& tipPCA)
{
  const KinematicParticle* leadTrack = get_leadTrack(daughters);
  if ( !leadTrack )
  {
    std::cerr << "WARNING: Failed to find leading track of tau --> returning null vector !!" << std::endl;
    return reco::Candidate::Point(0.,0.,0.);
  }
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(leadTrack->get_p4(), nullptr, nullptr, kBeam, r, n, k);
  double dr = rnd_.Gaus(0., resolutions_->get_tipResolution_perp());
  double dn = rnd_.Gaus(0., resolutions_->get_tipResolution_perp());
  double smeared_pcaX = tipPCA.x();
  double smeared_pcaY = tipPCA.y();
  double smeared_pcaZ = tipPCA.z();
  if ( applySmearing_tip_perp_ )
  {
    smeared_pcaX += dr*r.x() + dn*n.x();
    smeared_pcaY += dr*r.y() + dn*n.y();
    smeared_pcaZ += dr*r.z() + dn*n.z();
  }
  reco::Candidate::Point smeared_tipPCA(smeared_pcaX, smeared_pcaY, smeared_pcaZ);
  return smeared_tipPCA;
}

