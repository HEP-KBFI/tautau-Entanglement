#include "TauAnalysis/Entanglement/interface/StartPosFinder.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"     // cmsException
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerBase.h" // SpinAnalyzerBase::kTauPlus, SpinAnalyzerBase::kTauMinus
#include "TauAnalysis/Entanglement/interface/StartPosFinder1.h"  // StartPosFinder1
#include "TauAnalysis/Entanglement/interface/StartPosFinder2.h"  // StartPosFinder2

#include <iostream>                                              // std::cout

StartPosFinder::StartPosFinder(const edm::ParameterSet& cfg)
  : algo_(nullptr)
  , mode_(cfg.getParameter<int>("startPosMode"))
  , spinAnalyzer_(cfg)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  if      ( mode_ == 1 ) algo_ = new StartPosFinder1(cfg);
  else if ( mode_ == 2 ) algo_ = new StartPosFinder2(cfg);
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

  reco::Candidate::Vector hPlus = spinAnalyzer_(kineEvt_startpos, SpinAnalyzerBase::kTauPlus);
  kineEvt_startpos.hPlus_ = hPlus;
  kineEvt_startpos.hPlus_isValid_ = true;
  reco::Candidate::Vector hMinus = spinAnalyzer_(kineEvt_startpos, SpinAnalyzerBase::kTauMinus);
  kineEvt_startpos.hMinus_ = hMinus;
  kineEvt_startpos.hMinus_isValid_ = true;

  return kineEvt_startpos;
}
