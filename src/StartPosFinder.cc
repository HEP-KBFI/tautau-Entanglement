#include "TauAnalysis/Entanglement/interface/StartPosFinder.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"               // cmsException
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoBase.h" // PolarimetricVector::kTauPlus, PolarimetricVector::kTauMinus
#include "TauAnalysis/Entanglement/interface/StartPosAlgo1.h"              // StartPosAlgo1
#include "TauAnalysis/Entanglement/interface/StartPosAlgo2.h"              // StartPosAlgo2

#include <iostream>                                                        // std::cout

StartPosFinder::StartPosFinder(const edm::ParameterSet& cfg)
  : algo_(nullptr)
  , mode_(cfg.getParameter<int>("startPosMode"))
  , polarimetricVector_(cfg)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  if      ( mode_ == 1 ) algo_ = new StartPosAlgo1(cfg);
  else if ( mode_ == 2 ) algo_ = new StartPosAlgo2(cfg);
  else throw cmsException("StartPosFinder::StartPosFinder", __LINE__) 
    << "Invalid Configuration parameter 'mode' = " << mode_ << " !!\n";
}

StartPosFinder::~StartPosFinder()
{
  delete algo_;
}

KinematicEvent
StartPosFinder::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<StartPosFinder::operator()>:\n";
  }

  KinematicEvent kineEvt_startpos = (*algo_)(kineEvt);

  reco::Candidate::Vector hPlus = polarimetricVector_(kineEvt_startpos, pol::kTauPlus);
  kineEvt_startpos.hPlus_ = hPlus;
  kineEvt_startpos.hPlus_isValid_ = true;
  reco::Candidate::Vector hMinus = polarimetricVector_(kineEvt_startpos, pol::kTauMinus);
  kineEvt_startpos.hMinus_ = hMinus;
  kineEvt_startpos.hMinus_isValid_ = true;

  return kineEvt_startpos;
}
