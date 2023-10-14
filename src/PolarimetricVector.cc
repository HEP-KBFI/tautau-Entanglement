#include "TauAnalysis/Entanglement/interface/PolarimetricVector.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#include "DataFormats/TauReco/interface/PFTau.h"                           // reco::PFTau::hadronicDecayMode
#pragma GCC diagnostic pop

#include "TauAnalysis/Entanglement/interface/cmsException.h"               // cmsException
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h"  // kBeam, kHiggs
#include "TauAnalysis/Entanglement/interface/printVector.h"                // printVector()
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoBase.h" // PolarimetricVector::kTauPlus, PolarimetricVector::kTauMinus

PolarimetricVector::PolarimetricVector(const edm::ParameterSet& cfg)
  : algoOneProng0Pi0_(cfg)
  , algoOneProng1Pi0_(cfg)
  , algoThreeProng0Pi0_(cfg)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{}

PolarimetricVector::~PolarimetricVector()
{}

reco::Candidate::Vector
PolarimetricVector::operator()(const KinematicEvent& kineEvt, int tau)
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVector::operator()>:\n";
  }

  reco::Candidate::Vector h;
  if ( kineEvt.tauPlusP4_isValid() && kineEvt.tauMinusP4_isValid() )
  {
    int tau_decayMode = -1;
    if      ( tau == pol::kTauPlus  ) tau_decayMode = kineEvt.tauPlus_decayMode();
    else if ( tau == pol::kTauMinus ) tau_decayMode = kineEvt.tauMinus_decayMode();
    else assert(0);
    if ( verbosity_ >= 2 )
    {
      std::cout << " tau_decayMode = " << tau_decayMode << "\n";
    }
 
    if ( tau_decayMode == reco::PFTau::kOneProng0PiZero )
    {
      h = algoOneProng0Pi0_(kineEvt, tau);
    }
    else if ( tau_decayMode == reco::PFTau::kOneProng1PiZero )
    {
      h = algoOneProng1Pi0_(kineEvt, tau);
    }
    //else if ( tau_decayMode == reco::PFTau::kThreeProng0PiZero )
    //{
    //  h = algoThreeProng0Pi0_(kineEvt, tau);
    //}
  }
  if ( verbosity_ >= 1 )
  {
    printVector("h", h);
  }

  return h;
}
