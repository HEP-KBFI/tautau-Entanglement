#include "TauAnalysis/Entanglement/interface/SpinAnalyzerThreeProng0Pi0.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/getP4_hf.h"                  // getP4_hf()
#include "TauAnalysis/Entanglement/interface/getP4_rf.h"                  // getP4_rf()

#include <Math/Boost.h>                                                   // Boost

SpinAnalyzerThreeProng0Pi0::SpinAnalyzerThreeProng0Pi0(const edm::ParameterSet& cfg)
  : SpinAnalyzerBase(cfg)
{}

SpinAnalyzerThreeProng0Pi0::~SpinAnalyzerThreeProng0Pi0()
{}

namespace
{
  reco::Candidate::Vector
  getPolarimetricVec_ThreeProng0PiZero(const reco::Candidate::LorentzVector& tauP4,
                                     const std::vector<KinematicParticle>& daughters,
                                     const ROOT::Math::Boost& boost_ttrf,
                                     const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                                     const ROOT::Math::Boost& boost_trf,
                                     int verbosity = 0, bool cartesian = true)
  {
    if ( verbosity >= 2 )
    {
      std::cout << "<getPolarimetricVec_ThreeProng0PiZero>:\n";
    }

    throw cmsException("SpinAnalyzerThreeProng0Pi0", __LINE__)
      << "Function 'getPolarimetricVec_ThreeProng0PiZero' not implemented yet !!\n";

    return reco::Candidate::Vector(0., 0., 0.);
  }
}

reco::Candidate::Vector
SpinAnalyzerThreeProng0Pi0::operator()(const KinematicEvent& evt, int tau)
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<SpinAnalyzerThreeProng0Pi0::operator()>:\n";
  }

  reco::Candidate::LorentzVector tauP4;
  const std::vector<KinematicParticle>* daughters = nullptr;
  if ( tau == SpinAnalyzerBase::kTauPlus )
  {
    tauP4 = evt.tauPlusP4();
    daughters = &evt.daughtersTauPlus();
  }
  else if ( tau == SpinAnalyzerBase::kTauMinus )
  {
    tauP4 = evt.tauMinusP4();
    daughters = &evt.daughtersTauMinus();
  }
  else assert(0);
  reco::Candidate::LorentzVector higgsP4 = evt.tauPlusP4() + evt.tauMinusP4();
  ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(higgsP4.BoostToCM());
  reco::Candidate::LorentzVector tauP4_ttrf = getP4_rf(tauP4, boost_ttrf);
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(evt.tauMinusP4(), &higgsP4, &boost_ttrf, hAxis_, r, n, k, verbosity_, cartesian_);
  reco::Candidate::LorentzVector tauP4_hf = getP4_hf(tauP4_ttrf, r, n, k);
  ROOT::Math::Boost boost_trf = ROOT::Math::Boost(tauP4_hf.BoostToCM());

  reco::Candidate::Vector h = getPolarimetricVec_ThreeProng0PiZero(tauP4, *daughters, boost_ttrf, r, n, k, boost_trf, verbosity_, cartesian_);
  return h;
}
