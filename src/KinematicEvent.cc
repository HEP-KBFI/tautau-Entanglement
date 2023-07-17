#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"       // cmsException
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h" // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printPoint.h"         // printPoint()
#include "TauAnalysis/Entanglement/interface/printVector.h"        // printVector()

#include <TString.h>                                               // Form()

KinematicEvent::KinematicEvent()
  : tauPlusP4_isValid_(false)
  , svTauPlus_isValid_(false)
  , tauMinusP4_isValid_(false)
  , svTauMinus_isValid_(false)
{}

KinematicEvent::~KinematicEvent()
{}

const reco::Candidate::Point&
KinematicEvent::get_pv() const
{
  return pv_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::get_recoilP4() const
{
  return recoilP4_;
}

reco::Candidate::LorentzVector
KinematicEvent::get_tauPlusP4() const
{
  return tauPlusP4_;
}

bool
KinematicEvent::get_tauPlusP4_isValid() const
{
  return tauPlusP4_isValid_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::get_visTauPlusP4() const
{
  return visTauPlusP4_;
}

int
KinematicEvent::get_tauPlus_decayMode() const
{
  return tauPlus_decayMode_;
}

const std::vector<KinematicParticle>&
KinematicEvent::get_daughtersTauPlus() const
{
  return daughtersTauPlus_;
}

const reco::Candidate::Point&
KinematicEvent::get_tipPCATauPlus() const
{
  return tipPCATauPlus_;
}

const reco::Candidate::Point&
KinematicEvent::get_svTauPlus() const
{
  return svTauPlus_;
}

bool
KinematicEvent::get_svTauPlus_isValid() const
{
  return svTauPlus_isValid_;
}

const reco::Candidate::Vector&
KinematicEvent::get_hPlus() const
{
  return hPlus_;
}

bool
KinematicEvent::get_hPlus_isValid() const
{
  return hPlus_isValid_;
}

reco::Candidate::LorentzVector
KinematicEvent::get_tauMinusP4() const
{
  return tauMinusP4_;
}

bool
KinematicEvent::get_tauMinusP4_isValid() const
{
  return tauMinusP4_isValid_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::get_visTauMinusP4() const
{
  return visTauMinusP4_;
}

int
KinematicEvent::get_tauMinus_decayMode() const
{
  return tauMinus_decayMode_;
}

const std::vector<KinematicParticle>&
KinematicEvent::get_daughtersTauMinus() const
{
  return daughtersTauMinus_;
}

const reco::Candidate::Point&
KinematicEvent::get_tipPCATauMinus() const
{
  return tipPCATauMinus_;
}

const reco::Candidate::Point&
KinematicEvent::get_svTauMinus() const
{
  return svTauMinus_;
}

bool
KinematicEvent::get_svTauMinus_isValid() const
{
  return svTauMinus_isValid_;
}

const reco::Candidate::Vector&
KinematicEvent::get_hMinus() const
{
  return hMinus_;
}

bool
KinematicEvent::get_hMinus_isValid() const
{
  return hMinus_isValid_;
}

void
printKinematicEvent(const std::string& label,
                    const KinematicEvent& kineEvt,
                    bool cartesian)
{
  std::cout << label << ":\n";

  printPoint("pv", kineEvt.get_pv());

  printLorentzVector("recoilP4", kineEvt.get_recoilP4(), cartesian);

  if ( kineEvt.get_tauPlusP4_isValid() )
  {
    printLorentzVector("tauPlusP4", kineEvt.get_tauPlusP4(), cartesian);
  }
  else
  {
    std::cout << "tauPlusP4: N/A\n";
  }
  printLorentzVector("visTauPlusP4", kineEvt.get_visTauPlusP4(), cartesian);
  std::cout << "tauPlus_decayMode = " << kineEvt.get_tauPlus_decayMode() << "\n";
  std::cout << "daughtersTauPlus:\n";
  const std::vector<KinematicParticle>& daughtersTauPlus = kineEvt.get_daughtersTauPlus();
  for ( size_t idx = 0; idx < daughtersTauPlus.size(); ++idx )
  {
    const KinematicParticle& daughter = daughtersTauPlus[idx];
    std::string label = Form("#%i", (int)idx);
    printKinematicParticle(label, daughter, cartesian);
  }
  printPoint("tipPCATauPlus", kineEvt.get_tipPCATauPlus());
  if ( kineEvt.get_svTauPlus_isValid() )
  {
    printPoint("svTauPlus", kineEvt.get_svTauPlus());
  }
  else
  {
    std::cout << "svTauPlus: N/A\n";
  }
  if ( kineEvt.get_hPlus_isValid() )
  {
    printVector("hPlus", kineEvt.get_hPlus());
  }
  else
  {
    std::cout << "hPlus: N/A\n";
  }

  if ( kineEvt.get_tauMinusP4_isValid() )
  {
    printLorentzVector("tauMinusP4", kineEvt.get_tauMinusP4(), cartesian);
  }
  else
  {
    std::cout << "tauMinusP4: N/A\n";
  }
  printLorentzVector("visTauMinusP4", kineEvt.get_visTauMinusP4(), cartesian);
  std::cout << "tauMinus_decayMode = " << kineEvt.get_tauMinus_decayMode() << "\n";
  std::cout << "daughtersTauMinus:\n";
  const std::vector<KinematicParticle>& daughtersTauMinus = kineEvt.get_daughtersTauMinus();
  for ( size_t idx = 0; idx < daughtersTauMinus.size(); ++idx )
  {
    const KinematicParticle& daughter = daughtersTauMinus[idx];
    std::string label = Form("#%i", (int)idx);
    printKinematicParticle(label, daughter, cartesian);
  }
  printPoint("tipPCATauMinus", kineEvt.get_tipPCATauMinus());
  if ( kineEvt.get_svTauMinus_isValid() )
  {
    printPoint("svTauMinus", kineEvt.get_svTauMinus());
  }
  else
  {
    std::cout << "svTauMinus: N/A\n";
  }
  if ( kineEvt.get_hMinus_isValid() )
  {
    printVector("hMinus", kineEvt.get_hMinus());
  }
  else
  {
    std::cout << "hMinus: N/A\n";
  }
}
