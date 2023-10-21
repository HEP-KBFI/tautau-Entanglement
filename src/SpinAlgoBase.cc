#include "TauAnalysis/Entanglement/interface/SpinAlgoBase.h"

#include "TauAnalysis/Entanglement/interface/comp_concurrence.h"  // comp_concurrence()
#include "TauAnalysis/Entanglement/interface/comp_Ek.h"           // comp_Ek()
#include "TauAnalysis/Entanglement/interface/comp_Rchsh.h"        // comp_Rchsh()
#include "TauAnalysis/Entanglement/interface/comp_steerability.h" // comp_steerability()

using namespace spin;

SpinAlgoBase::SpinAlgoBase(const edm::ParameterSet& cfg)
  : verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
{}

SpinAlgoBase::~SpinAlgoBase()
{}

void
SpinAlgoBase::set_verbosity(int verbosity)
{
  verbosity_ = verbosity;
}

void
SpinAlgoBase::addEntanglementVariables(spin::Measurement& measurement)
{
  measurement.concurrence_ = comp_concurrence(measurement.Bp_, measurement.Bm_, measurement.C_, verbosity_);
  measurement.Ek_ = comp_Ek(measurement.C_);
  measurement.Rchsh_ = comp_Rchsh(measurement.C_, verbosity_);
  measurement.steerability_ = comp_steerability(measurement.C_, 360, 360);
}
