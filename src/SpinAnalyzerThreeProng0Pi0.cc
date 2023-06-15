#include "TauAnalysis/Entanglement/interface/SpinAnalyzerThreeProng0Pi0.h"

#include "TauAnalysis/Entanglement/interface/auxFunctions.h"              // getP4_rf(), getP4_hf(), getP4_ttrf_hf_trf()
#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()

#include <Math/Boost.h>                                                   // Boost

SpinAnalyzerThreeProng0Pi0::SpinAnalyzerThreeProng0Pi0(const edm::ParameterSet& cfg)
  : SpinAnalyzer(cfg)
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
                                     const ROOT::Math::Boost& boost_trf)
  {
    //std::cout << "<getPolarimetricVec_ThreeProng0PiZero>:" << std::endl;

    throw cmsException("SpinAnalyzerThreeProng0Pi0", __LINE__)
      << "Function 'getPolarimetricVec_ThreeProng0PiZero' not implemented yet !!\n";

    return reco::Candidate::Vector(0., 0., 0.);
  }
}

reco::Candidate::Vector
SpinAnalyzerThreeProng0Pi0::operator()(const KinematicEvent& evt, int tau)
{
  reco::Candidate::LorentzVector tauP4;
  const std::vector<KinematicParticle>* daughters = nullptr;
  if ( tau == kTauPlus )
  {
    tauP4 = evt.get_tauPlusP4();
    daughters = &evt.get_daughtersTauPlus();
  }
  else if ( tau == kTauMinus )
  {
    tauP4 = evt.get_tauMinusP4();
    daughters = &evt.get_daughtersTauMinus();
  }
  else assert(0);
  reco::Candidate::LorentzVector recoilP4 = evt.get_recoilP4();
  ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(recoilP4.BoostToCM());
  reco::Candidate::LorentzVector tauP4_ttrf = getP4_rf(tauP4, boost_ttrf);
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(evt.get_tauMinusP4(), &recoilP4, &boost_ttrf, hAxis_, r, n, k, verbosity_);
  reco::Candidate::LorentzVector tauP4_hf = getP4_hf(tauP4_ttrf, r, n, k);
  ROOT::Math::Boost boost_trf = ROOT::Math::Boost(tauP4_hf.BoostToCM());

  reco::Candidate::Vector h = getPolarimetricVec_ThreeProng0PiZero(tauP4, *daughters, boost_ttrf, r, n, k, boost_trf);
  return h;
}
