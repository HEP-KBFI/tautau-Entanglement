#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"

#include "DataFormats/TauReco/interface/PFTau.h"                          // reco::PFTau::hadronicDecayMode

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // kBeam, kHiggs
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerBase.h"          // SpinAnalyzerBase::kTauPlus, SpinAnalyzerBase::kTauMinus

SpinAnalyzer::SpinAnalyzer(const edm::ParameterSet& cfg)
  : spinAnalyzerOneProng0Pi0_(cfg)
  , spinAnalyzerOneProng1Pi0_(cfg)
  , spinAnalyzerThreeProng0Pi0_(cfg)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{}

SpinAnalyzer::~SpinAnalyzer()
{}

reco::Candidate::Vector
SpinAnalyzer::operator()(const KinematicEvent& kineEvt, int tau)
{
  int tau_decayMode = -1;
  if      ( tau == SpinAnalyzerBase::kTauPlus  ) tau_decayMode = kineEvt.tauPlus_decayMode();
  else if ( tau == SpinAnalyzerBase::kTauMinus ) tau_decayMode = kineEvt.tauMinus_decayMode();
  else assert(0);

  reco::Candidate::Vector h;
  if ( tau_decayMode == reco::PFTau::kOneProng0PiZero )
  {
    h = spinAnalyzerOneProng0Pi0_(kineEvt, tau);
  }
  else if ( tau_decayMode == reco::PFTau::kOneProng1PiZero )
  {
    h = spinAnalyzerOneProng1Pi0_(kineEvt, tau);
  }
  //else if ( tau_decayMode == reco::PFTau::kThreeProng0PiZero )
  //{
  //  h = spinAnalyzerThreeProng0Pi0_(kineEvt, tau);
  //}
  return h;
}
