#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"       // cmsException
#include "TauAnalysis/Entanglement/interface/comp_cosThetaGJ.h"    // comp_cosThetaGJ(), comp_cosThetaGJ_solution()
#include "TauAnalysis/Entanglement/interface/printDistance.h"      // printDistance()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h" // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printPoint.h"         // printPoint()
#include "TauAnalysis/Entanglement/interface/printVector.h"        // printVector()

#include <TString.h>                                               // Form()

#include <cmath>                                                   // std::acos()

KinematicEvent::KinematicEvent()
  : tauPlusP4_isValid_(false)
  , nuTauPlusCov_isValid_(false)
  , svTauPlus_isValid_(false)
  , tauMinusP4_isValid_(false)
  , nuTauMinusCov_isValid_(false)
  , svTauMinus_isValid_(false)
  , kinFitStatus_(-1)
  , kinFitChi2_(-1.)
  , kinFit_isValid_(false)
{}

KinematicEvent::~KinematicEvent()
{}

const reco::Candidate::Point&
KinematicEvent::pv() const
{
  return pv_;
}

const math::Matrix3x3&
KinematicEvent::pvCov() const
{
  return pvCov_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::recoilP4() const
{
  return recoilP4_;
}

const math::Matrix4x4&
KinematicEvent::recoilCov() const
{
  return recoilCov_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::tauPlusP4() const
{
  return tauPlusP4_;
}

bool
KinematicEvent::tauPlusP4_isValid() const
{
  return tauPlusP4_isValid_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::visTauPlusP4() const
{
  return visTauPlusP4_;
}

const math::Matrix4x4&
KinematicEvent::visTauPlusCov() const
{
  return visTauPlusCov_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::nuTauPlusP4() const
{
  return nuTauPlusP4_;
}

bool
KinematicEvent::nuTauPlusP4_isValid() const
{
  return nuTauPlusP4_isValid_;
}

const math::Matrix3x3&
KinematicEvent::nuTauPlusCov() const
{
  return nuTauPlusCov_;
}

bool
KinematicEvent::nuTauPlusCov_isValid() const
{
  return nuTauPlusCov_isValid_;
}

int
KinematicEvent::tauPlus_decayMode() const
{
  return tauPlus_decayMode_;
}

const std::vector<KinematicParticle>&
KinematicEvent::daughtersTauPlus() const
{
  return daughtersTauPlus_;
}

const reco::Candidate::Point&
KinematicEvent::tipPCATauPlus() const
{
  return tipPCATauPlus_;
}

const reco::Candidate::Point&
KinematicEvent::svTauPlus() const
{
  return svTauPlus_;
}

const math::Matrix3x3&
KinematicEvent::svTauPlusCov() const
{
  return svTauPlusCov_;
}

bool
KinematicEvent::svTauPlus_isValid() const
{
  return svTauPlus_isValid_;
}

const reco::Candidate::Vector&
KinematicEvent::hPlus() const
{
  return hPlus_;
}

bool
KinematicEvent::hPlus_isValid() const
{
  return hPlus_isValid_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::tauMinusP4() const
{
  return tauMinusP4_;
}

bool
KinematicEvent::tauMinusP4_isValid() const
{
  return tauMinusP4_isValid_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::visTauMinusP4() const
{
  return visTauMinusP4_;
}

const reco::Candidate::LorentzVector&
KinematicEvent::nuTauMinusP4() const
{
  return nuTauMinusP4_;
}

bool
KinematicEvent::nuTauMinusP4_isValid() const
{
  return nuTauMinusP4_isValid_;
}

const math::Matrix3x3&
KinematicEvent::nuTauMinusCov() const
{
  return nuTauMinusCov_;
}

bool
KinematicEvent::nuTauMinusCov_isValid() const
{
  return nuTauMinusCov_isValid_;
}

int
KinematicEvent::tauMinus_decayMode() const
{
  return tauMinus_decayMode_;
}

const std::vector<KinematicParticle>&
KinematicEvent::daughtersTauMinus() const
{
  return daughtersTauMinus_;
}

const reco::Candidate::Point&
KinematicEvent::tipPCATauMinus() const
{
  return tipPCATauMinus_;
}

const reco::Candidate::Point&
KinematicEvent::svTauMinus() const
{
  return svTauMinus_;
}

const math::Matrix3x3&
KinematicEvent::svTauMinusCov() const
{
  return svTauMinusCov_;
}

bool
KinematicEvent::svTauMinus_isValid() const
{
  return svTauMinus_isValid_;
}

const reco::Candidate::Vector&
KinematicEvent::hMinus() const
{
  return hMinus_;
}

bool
KinematicEvent::hMinus_isValid() const
{
  return hMinus_isValid_;
}

int
KinematicEvent::kinFitStatus() const
{
  return kinFitStatus_;
}

double
KinematicEvent::kinFitChi2() const
{
  return kinFitChi2_;
}

const math::MatrixPxP&
KinematicEvent::kinFitCov() const
{
  return kinFitCov_;
}

bool
KinematicEvent::kinFit_isValid() const
{
  return kinFit_isValid_;
}

void
printKinematicEvent(const std::string& label,
                    const KinematicEvent& kineEvt,
                    bool cartesian)
{
  std::cout << label << ":\n";

  printPoint("pv", kineEvt.pv());

  printLorentzVector("recoilP4", kineEvt.recoilP4(), cartesian);
  std::cout << " mass = " << kineEvt.recoilP4().mass() << "\n";

  if ( kineEvt.tauPlusP4_isValid() )
  {
    printLorentzVector("tauPlusP4", kineEvt.tauPlusP4(), cartesian);
    std::cout << " mass = " << kineEvt.tauPlusP4().mass() << "\n";
    printLorentzVector("tauPlusP4", kineEvt.tauPlusP4(), false);
  }
  else
  {
    std::cout << "tauPlusP4: N/A\n";
  }
  printLorentzVector("visTauPlusP4", kineEvt.visTauPlusP4(), cartesian);
  std::cout << " mass = " << kineEvt.visTauPlusP4().mass() << "\n";
  printLorentzVector("nuTauPlusP4", kineEvt.nuTauPlusP4(), cartesian);
  std::cout << " mass = " << kineEvt.nuTauPlusP4().mass() << "\n";
  std::cout << "tauPlus_decayMode = " << kineEvt.tauPlus_decayMode() << "\n";
  std::cout << "daughtersTauPlus:\n";
  const std::vector<KinematicParticle>& daughtersTauPlus = kineEvt.daughtersTauPlus();
  for ( size_t idx = 0; idx < daughtersTauPlus.size(); ++idx )
  {
    const KinematicParticle& daughter = daughtersTauPlus[idx];
    std::string label = Form("#%i", (int)idx);
    printKinematicParticle(label, daughter, cartesian);
  }
  double tauPlus_cosThetaGJ = comp_cosThetaGJ(kineEvt.tauPlusP4(), kineEvt.visTauPlusP4());
  std::cout << "Gottfied-Jackson angle = " << std::acos(tauPlus_cosThetaGJ) << " "
            << "(expected = " << std::acos(comp_cosThetaGJ_solution(kineEvt.tauPlusP4(), kineEvt.visTauPlusP4())) << ")\n";
  printPoint("tipPCATauPlus", kineEvt.tipPCATauPlus());
  if ( kineEvt.svTauPlus_isValid() )
  {
    printPoint("svTauPlus", kineEvt.svTauPlus());
    printDistance("svTauPlus - pv", kineEvt.svTauPlus() - kineEvt.pv(), cartesian);
    printDistance("svTauPlus - pv", kineEvt.svTauPlus() - kineEvt.pv(), false);
  }
  else
  {
    std::cout << "svTauPlus: N/A\n";
  }
  if ( kineEvt.hPlus_isValid() )
  {
    printVector("hPlus", kineEvt.hPlus());
  }
  else
  {
    std::cout << "hPlus: N/A\n";
  }

  if ( kineEvt.tauMinusP4_isValid() )
  {
    printLorentzVector("tauMinusP4", kineEvt.tauMinusP4(), cartesian);
    std::cout << " mass = " << kineEvt.tauMinusP4().mass() << "\n";
    printLorentzVector("tauMinusP4", kineEvt.tauMinusP4(), false);
  }
  else
  {
    std::cout << "tauMinusP4: N/A\n";
  }
  printLorentzVector("visTauMinusP4", kineEvt.visTauMinusP4(), cartesian);
  std::cout << " mass = " << kineEvt.visTauMinusP4().mass() << "\n";
  printLorentzVector("nuTauMinusP4", kineEvt.nuTauMinusP4(), cartesian);
  std::cout << " mass = " << kineEvt.nuTauMinusP4().mass() << "\n";
  std::cout << "tauMinus_decayMode = " << kineEvt.tauMinus_decayMode() << "\n";
  std::cout << "daughtersTauMinus:\n";
  const std::vector<KinematicParticle>& daughtersTauMinus = kineEvt.daughtersTauMinus();
  for ( size_t idx = 0; idx < daughtersTauMinus.size(); ++idx )
  {
    const KinematicParticle& daughter = daughtersTauMinus[idx];
    std::string label = Form("#%i", (int)idx);
    printKinematicParticle(label, daughter, cartesian);
  }
  double tauMinus_cosThetaGJ = comp_cosThetaGJ(kineEvt.tauMinusP4(), kineEvt.visTauMinusP4());
  std::cout << "Gottfied-Jackson angle = " << std::acos(tauMinus_cosThetaGJ) << " "
            << "(expected = " << std::acos(comp_cosThetaGJ_solution(kineEvt.tauMinusP4(), kineEvt.visTauMinusP4())) << ")\n";
  printPoint("tipPCATauMinus", kineEvt.tipPCATauMinus());
  if ( kineEvt.svTauMinus_isValid() )
  {
    printPoint("svTauMinus", kineEvt.svTauMinus());
    printDistance("svTauMinus - pv", kineEvt.svTauMinus() - kineEvt.pv(), cartesian);
    printDistance("svTauMinus - pv", kineEvt.svTauMinus() - kineEvt.pv(), false);
  }
  else
  {
    std::cout << "svTauMinus: N/A\n";
  }
  if ( kineEvt.hMinus_isValid() )
  {
    printVector("hMinus", kineEvt.hMinus());
  }
  else
  {
    std::cout << "hMinus: N/A\n";
  }
}
