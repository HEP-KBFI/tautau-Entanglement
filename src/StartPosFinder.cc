#include "TauAnalysis/Entanglement/interface/StartPosFinder.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"               // cmsException
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoBase.h" // PolarimetricVector::kTauPlus, PolarimetricVector::kTauMinus
#include "TauAnalysis/Entanglement/interface/StartPosAlgo1.h"              // StartPosAlgo1
#include "TauAnalysis/Entanglement/interface/StartPosAlgo2.h"              // StartPosAlgo2

#include <iostream>                                                        // std::cout

StartPosFinder::StartPosFinder(const edm::ParameterSet& cfg)
  : algo_(nullptr)
  , polarimetricVector_(cfg)
  , skip_(cfg.getParameter<bool>("skip"))
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  int algo = cfg.getParameter<int>("algo");
  if      ( algo == 1 ) algo_ = new StartPosAlgo1(cfg);
  else if ( algo == 2 ) algo_ = new StartPosAlgo2(cfg);
  else throw cmsException("StartPosFinder::StartPosFinder", __LINE__) 
    << "Invalid Configuration parameter 'algo' = " << algo << " !!\n";
}

StartPosFinder::~StartPosFinder()
{
  delete algo_;
}

int
StartPosFinder::get_algo() const
{
  return algo_->get_algo();
}

std::vector<KinematicEvent>
StartPosFinder::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<StartPosFinder::operator()>:\n";
  }

  if ( skip_ )
  {
    return { kineEvt };
  }

  std::vector<KinematicEvent> kineEvts_startPos = (*algo_)(kineEvt);
  for ( KinematicEvent& kineEvt_startPos : kineEvts_startPos )
  {
    reco::Candidate::Vector hPlus = polarimetricVector_(kineEvt_startPos, pol::kTauPlus);
    kineEvt_startPos.hPlus_ = hPlus;
    kineEvt_startPos.hPlus_isValid_ = true;
    reco::Candidate::Vector hMinus = polarimetricVector_(kineEvt_startPos, pol::kTauMinus);
    kineEvt_startPos.hMinus_ = hMinus;
    kineEvt_startPos.hMinus_isValid_ = true;
  }

  return kineEvts_startPos;
}
