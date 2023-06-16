#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"

#include "TauAnalysis/Entanglement/interface/auxFunctions.h" // printLorentzVector(), printPoint()
#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException

#include <TString.h>                                         // Form()

KinematicEvent::KinematicEvent(const reco::Candidate::Point& pv, const reco::Candidate::LorentzVector& recoilP4)
  : pv_(pv)
  , recoilP4_(recoilP4)
  , tauPlusP4_isValid_(false)
  , visTauPlus_isValid_(false)
  , svTauPlus_isValid_(false)
  , tauMinusP4_isValid_(false)
  , visTauMinus_isValid_(false)
  , svTauMinus_isValid_(false)
{}

KinematicEvent::~KinematicEvent()
{}

void
KinematicEvent::set_tauPlusP4(const reco::Candidate::LorentzVector& tauPlusP4)
{
  tauPlusP4_ = tauPlusP4;
  tauPlusP4_isValid_ = true;
}

void
KinematicEvent::set_visTauPlus(const reco::Candidate::LorentzVector& visTauPlusP4, int tauPlus_decaymode, 
                               const std::vector<KinematicParticle>& daughtersTauPlus,
                               const reco::Candidate::Point& tipPCATauPlus)
{
  visTauPlusP4_ = visTauPlusP4;
  tauPlus_decaymode_ = tauPlus_decaymode;
  daughtersTauPlus_ = daughtersTauPlus;
  tipPCATauPlus_ = tipPCATauPlus;
  visTauPlus_isValid_ = true;
}

void
KinematicEvent::set_svTauPlus(const reco::Candidate::Point& svTauPlus)
{
  svTauPlus_ = svTauPlus;
  svTauPlus_isValid_ = true;
}

void
KinematicEvent::set_tauMinusP4(const reco::Candidate::LorentzVector& tauMinusP4)
{
  tauMinusP4_ = tauMinusP4;
  tauMinusP4_isValid_ = true;
}

void
KinematicEvent::set_visTauMinus(const reco::Candidate::LorentzVector& visTauMinusP4, int tauMinus_decaymode, 
                                const std::vector<KinematicParticle>& daughtersTauMinus,
                                const reco::Candidate::Point& tipPCATauMinus)
{
  visTauMinusP4_ = visTauMinusP4;
  tauMinus_decaymode_ = tauMinus_decaymode;
  daughtersTauMinus_ = daughtersTauMinus;
  tipPCATauMinus_ = tipPCATauMinus;
  visTauMinus_isValid_ = true;
}

void
KinematicEvent::set_svTauMinus(const reco::Candidate::Point& svTauMinus)
{
  svTauMinus_ = svTauMinus;
  svTauMinus_isValid_ = true;
}

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
  if ( !tauPlusP4_isValid_ )
    throw cmsException("KinematicParticle", __LINE__)
      << "Tau+ four-vector not initialized !!\n";
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
KinematicEvent::get_tauPlus_decaymode() const
{
  return tauPlus_decaymode_;
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
  if ( !svTauPlus_isValid_ )
    throw cmsException("KinematicParticle", __LINE__)
      << "Tau+ decay vertex not initialized !!\n";
  return svTauPlus_;
}

bool
KinematicEvent::get_svTauPlus_isValid() const
{
  return svTauPlus_isValid_;
}

reco::Candidate::LorentzVector
KinematicEvent::get_tauMinusP4() const
{
  if ( !tauMinusP4_isValid_ )
    throw cmsException("KinematicParticle", __LINE__)
      << "Tau- four-vector not initialized !!\n";
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
KinematicEvent::get_tauMinus_decaymode() const
{
  return tauMinus_decaymode_;
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
  if ( !svTauMinus_isValid_ )
    throw cmsException("KinematicParticle", __LINE__)
      << "Tau- decay vertex not initialized !!\n";
  return svTauMinus_;
}

bool
KinematicEvent::get_svTauMinus_isValid() const
{
  return svTauMinus_isValid_;
}

void
printKinematicEvent(const std::string& label,
                    const KinematicEvent& evt,
                    bool cartesian)
{
  std::cout << label << ":\n";

  printPoint("pv", evt.get_pv());

  printLorentzVector("recoilP4", evt.get_recoilP4(), cartesian);

  if ( evt.get_tauPlusP4_isValid() )
  {
    printLorentzVector("tauPlusP4", evt.get_tauPlusP4(), cartesian);
  }
  else
  {
    std::cout << "tauPlusP4: N/A\n";
  }
  printLorentzVector("visTauPlusP4", evt.get_visTauPlusP4(), cartesian);
  std::cout << "tauPlus_decaymode = " << evt.get_tauPlus_decaymode() << "\n";
  std::cout << "daughtersTauPlus:\n";
  const std::vector<KinematicParticle>& daughtersTauPlus = evt.get_daughtersTauPlus();
  for ( size_t idx = 0; idx < daughtersTauPlus.size(); ++idx )
  {
    const KinematicParticle& daughter = daughtersTauPlus[idx];
    std::string label = Form("#%i", (int)idx);
    printKinematicParticle(label, daughter, cartesian);
  }
  printPoint("tipPCATauPlus", evt.get_tipPCATauPlus());
  if ( evt.get_svTauPlus_isValid() )
  {
    printPoint("svTauPlus", evt.get_svTauPlus());
  }
  else
  {
    std::cout << "svTauPlus: N/A\n";
  }
  
  if ( evt.get_tauMinusP4_isValid() )
  {
    printLorentzVector("tauMinusP4", evt.get_tauMinusP4(), cartesian);
  }
  else
  {
    std::cout << "tauMinusP4: N/A\n";
  }
  printLorentzVector("visTauMinusP4", evt.get_visTauMinusP4(), cartesian);
  std::cout << "tauMinus_decaymode = " << evt.get_tauMinus_decaymode() << "\n";
  std::cout << "daughtersTauMinus:\n";
  const std::vector<KinematicParticle>& daughtersTauMinus = evt.get_daughtersTauMinus();
  for ( size_t idx = 0; idx < daughtersTauMinus.size(); ++idx )
  {
    const KinematicParticle& daughter = daughtersTauMinus[idx];
    std::string label = Form("#%i", (int)idx);
    printKinematicParticle(label, daughter, cartesian);
  }
  printPoint("tipPCATauMinus", evt.get_tipPCATauMinus());
  if ( evt.get_svTauMinus_isValid() )
  {
    printPoint("svTauMinus", evt.get_svTauMinus());
  }
  else
  {
    std::cout << "svTauMinus: N/A\n";
  }
}
