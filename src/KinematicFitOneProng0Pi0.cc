#include "TauAnalysis/Entanglement/interface/KinematicFitOneProng0Pi0.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException

KinematicFitOneProng0Pi0::KinematicFitOneProng0Pi0()
{}

KinematicFitOneProng0Pi0::~KinematicFitOneProng0Pi0()
{}

KinematicEvent
KinematicFitOneProng0Pi0::operator()(const KinematicEvent& evt)
{
  const std::vector<KinematicParticle>& daughtersTauPlus = evt.get_daughtersTauPlus();
  for ( const KinematicParticle& daughterTauPlus : daughtersTauPlus )
  {
    if ( daughterTauPlus.get_pdgId() == +211 )
    {
      piPlus_ = &daughterTauPlus;
    }
  }
  
  const std::vector<KinematicParticle>& daughtersTauMinus = evt.get_daughtersTauMinus();
  for ( const KinematicParticle& daughterTauMinus : daughtersTauMinus )
  {
    if ( daughterTauMinus.get_pdgId() == +211 )
    {
      piMinus_ = &daughterTauMinus;
    }
  }

  KinematicEvent fitted_evt(evt);
  
  return fitted_evt;
}

void
KinematicFitOneProng0Pi0::findStartPosition()
{

}
