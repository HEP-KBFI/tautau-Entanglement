#include "TauAnalysis/Entanglement/interface/Smearing.h"

#include "DataFormats/Candidate/interface/Candidate.h"            // Candidate::LorentzVector, Candidate::Point, Candidate::Vector

#include "TauAnalysis/Entanglement/interface/KinematicParticle.h" // KinematicParticle

Smearing::Smearing(const edm::ParameterSet& cfg, const Resolutions& resolutions)
  : recoilSmearPx_(cfg.getParameter<bool>("recoilSmearPx"))
  , recoilSmearPy_(cfg.getParameter<bool>("recoilSmearPy"))
  , recoilSmearPz_(cfg.getParameter<bool>("recoilSmearPz"))
  , recoilSmearE_(cfg.getParameter<bool>("recoilSmearE"))
  , pvSmearXY_(cfg.getParameter<bool>("pvSmearXY"))
  , pvSmearZ_(cfg.getParameter<bool>("pvSmearZ"))
  , svSmearPerp_(cfg.getParameter<bool>("svSmearPerp"))
  , svSmearParl_(cfg.getParameter<bool>("svSmearParl"))
  , tipSmearPerp_(cfg.getParameter<bool>("tipSmearPerp"))
  , resolutions_(resolutions)
{}

Smearing::~Smearing()
{}

namespace
{
  void
  get_localCoordinateSystem(const reco::Candidate::Vector& p,
                            reco::Candidate::Vector& r, reco::Candidate::Vector& n, reco::Candidate::Vector& k)
  {
    k = p.unit();

    const double mProton = 0.938272;

    const double beamE = 7.e+3;
    const double beamPx = 0.;
    const double beamPy = 0.;
    const double beamPz = std::sqrt(beamE*beamE - mProton*mProton);
    reco::Candidate::LorentzVector beamP4(beamPx, beamPy, beamPz, beamE);
    reco::Candidate::Vector h = beamP4.Vect().unit();

    double cosTheta = k.Dot(h);
    assert(cosTheta >= -1. && cosTheta <= +1.);
    double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
    r = (h - k*cosTheta)*(1./sinTheta);

    n = k.Cross(r);
  }
}

KinematicEvent
Smearing::operator()(const KinematicEvent& evt)
{
  KinematicEvent smeared_evt(evt);

  const reco::Candidate::Point& pv = evt.get_pv();
  double smeared_pvX = pv.x() + rnd_.Gaus(0., resolutions_.get_pvResolutionXY());
  double smeared_pvY = pv.y() + rnd_.Gaus(0., resolutions_.get_pvResolutionXY());
  double smeared_pvZ = pv.z() + rnd_.Gaus(0., resolutions_.get_pvResolutionZ());
  smeared_evt.pv_ = reco::Candidate::Point(smeared_pvX, smeared_pvY, smeared_pvZ);

  const reco::Candidate::LorentzVector& recoilP4 = evt.get_recoilP4();
  double smeared_recoilPx = recoilP4.px()     + rnd_.Gaus(0., resolutions_.get_recoilResolutionPx());
  double smeared_recoilPy = recoilP4.py()     + rnd_.Gaus(0., resolutions_.get_recoilResolutionPy());
  double smeared_recoilPz = recoilP4.pz()     + rnd_.Gaus(0., resolutions_.get_recoilResolutionPz());
  double smeared_recoilE  = recoilP4.energy() + rnd_.Gaus(0., resolutions_.get_recoilResolutionE());
  smeared_evt.recoilP4_ = reco::Candidate::LorentzVector(smeared_recoilPx, smeared_recoilPy, smeared_recoilPz, smeared_recoilE);

  smeared_evt.tauPlusP4_ = reco::Candidate::LorentzVector(0., 0., 0., 0.);
  smeared_evt.tauPlusP4_isValid_ = false;
  const std::vector<KinematicParticle>& daughtersTauPlus = evt.get_daughtersTauPlus();
  smeared_evt.daughtersTauPlus_.clear();
  for ( const KinematicParticle& daughterTauPlus : daughtersTauPlus )
  {
    smeared_evt.daughtersTauPlus_.push_back(smear_daughter(daughterTauPlus));
  }
  if ( evt.get_svTauPlus_isValid() )
  {
    const reco::Candidate::Point& svTauPlus = evt.get_svTauPlus();
    reco::Candidate::Vector p(svTauPlus.x() - pv.x(), svTauPlus.y() - pv.y(), svTauPlus.z() - pv.z());
    smeared_evt.svTauPlus_ = smear_sv(p, evt.get_svTauPlus());
  }

  smeared_evt.tauMinusP4_ = reco::Candidate::LorentzVector(0., 0., 0., 0.);
  smeared_evt.tauMinusP4_isValid_ = false;
  const std::vector<KinematicParticle>& daughtersTauMinus = evt.get_daughtersTauMinus();
  smeared_evt.daughtersTauMinus_.clear();
  for ( const KinematicParticle& daughterTauMinus : daughtersTauMinus )
  {
    smeared_evt.daughtersTauMinus_.push_back(smear_daughter(daughterTauMinus));
  }
  if ( evt.get_svTauMinus_isValid() )
  {
    const reco::Candidate::Point& svTauMinus = evt.get_svTauMinus();
    reco::Candidate::Vector p(svTauMinus.x() - pv.x(), svTauMinus.y() - pv.y(), svTauMinus.z() - pv.z());
    smeared_evt.svTauMinus_ = smear_sv(p, svTauMinus);
  }

  return smeared_evt;
}

KinematicParticle
Smearing::smear_daughter(const KinematicParticle& daughter)
{
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(daughter.get_p4().Vect(), r, n, k);
  double dr = rnd_.Gaus(0., resolutions_.get_svResolutionPerp());
  double dn = rnd_.Gaus(0., resolutions_.get_svResolutionPerp());
  double dk = rnd_.Gaus(0., resolutions_.get_svResolutionParl());
  const reco::Candidate::Point& vertex = daughter.get_vertex();
  double smeared_vertexX = vertex.x() + dr*r.x() + dn*n.x() + dk*k.x();
  double smeared_vertexY = vertex.y() + dr*r.y() + dn*n.y() + dk*k.y();
  double smeared_vertexZ = vertex.z() + dr*r.z() + dn*n.z() + dk*k.z();
  KinematicParticle smeared_daughter(daughter);
  smeared_daughter.vertex_ = reco::Candidate::Point(smeared_vertexX, smeared_vertexY, smeared_vertexZ);
  return smeared_daughter;
}

reco::Candidate::Point
Smearing::smear_sv(const reco::Candidate::Vector& p, const reco::Candidate::Point& sv)
{
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(p, r, n, k);
  double dr = rnd_.Gaus(0., resolutions_.get_svResolutionPerp());
  double dn = rnd_.Gaus(0., resolutions_.get_svResolutionPerp());
  double dk = rnd_.Gaus(0., resolutions_.get_svResolutionParl());
  double smeared_svX = sv.x() + dr*r.x() + dn*n.x() + dk*k.x();
  double smeared_svY = sv.y() + dr*r.y() + dn*n.y() + dk*k.y();
  double smeared_svZ = sv.z() + dr*r.z() + dn*n.z() + dk*k.z();
  return reco::Candidate::Point(smeared_svX, smeared_svY, smeared_svZ);
}
