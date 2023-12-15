#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoThreeProng0Pi0.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/getP4_hf.h"                  // getP4_hf()
#include "TauAnalysis/Entanglement/interface/getP4_rf.h"                  // getP4_rf()
#include "TauAnalysis/Entanglement/interface/getP4_ttrf_hf_trf.h"         // getP4_ttrf_hf_trf()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"        // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printVector.h"               // printVector()

#include <Math/Boost.h>                                                   // Boost
#include <TLorentzVector.h>                                               // TLorentzVector
#include <TVector3.h>                                                     // TVector3

#include <algorithm>                                                      // std::sort()
#include <cmath>                                                          // std::abs()

PolarimetricVectorAlgoThreeProng0Pi0::PolarimetricVectorAlgoThreeProng0Pi0(const edm::ParameterSet& cfg)
  : PolarimetricVectorAlgoBase(cfg)
{}

PolarimetricVectorAlgoThreeProng0Pi0::~PolarimetricVectorAlgoThreeProng0Pi0()
{}

namespace
{
  bool
  isHigherPt(const KinematicParticle* particle1, const KinematicParticle* particle2)
  {
    return particle1->p4().pt() > particle2->p4().pt();
  }

  reco::Candidate::Vector
  getPolarimetricVec_ThreeProng0PiZero(const PolarimetricVectorTau2a1& a1pol,
                                       const reco::Candidate::LorentzVector& tauP4,
                                       const std::vector<KinematicParticle>& daughters,
                                       const ROOT::Math::Boost& boost_ttrf,
                                       const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                                       const ROOT::Math::Boost& boost_trf,
                                       int verbosity = 0, bool cartesian = true)
  {
    if ( verbosity >= 4 )
    {
      std::cout << "<getPolarimetricVec_ThreeProng0PiZero>:\n";
    }

    std::vector<const KinematicParticle*> chsPlus;
    std::vector<const KinematicParticle*> chsMinus;
    int charge_sum = 0;
    for ( const KinematicParticle& daughter : daughters )
    {
      if ( abs(daughter.pdgId()) == 211 || abs(daughter.pdgId()) == 321 )
      {
        if ( daughter.charge() > 0. ) chsPlus.push_back(&daughter);
        if ( daughter.charge() < 0. ) chsMinus.push_back(&daughter);
        charge_sum += daughter.charge();
      }
    }
    if ( !((chsPlus.size() + chsMinus.size()) == 3 && std::abs(charge_sum) == 1) )
      throw cmsException("getPolarimetricVec_ThreeProng0PiZero", __LINE__)
        << "Failed to find three charged pions with |sum(charge)| = 1 !!\n";

    std::vector<const KinematicParticle*> chsOS;
    std::vector<const KinematicParticle*> chsSS;
    if ( charge_sum > 0. )
    {
      chsOS = chsMinus;
      chsSS = chsPlus;
    }
    else
    {
      chsOS = chsPlus;
      chsSS = chsMinus;
    }
    assert(chsOS.size() == 1 && chsSS.size() == 2);
    const reco::Candidate::LorentzVector& chOSP4 = chsOS[0]->p4();
    std::sort(chsSS.begin(), chsSS.end(), isHigherPt);
    const reco::Candidate::LorentzVector& chSS1P4 = chsSS[0]->p4();
    const reco::Candidate::LorentzVector& chSS2P4 = chsSS[1]->p4();
    if ( verbosity >= 4 )
    { 
      std::cout << "in laboratory frame:\n";
      printLorentzVector("tauP4",   tauP4,       cartesian);
      printLorentzVector("tauP4",   tauP4,       false);
      printLorentzVector("chOSP4",  chOSP4,      cartesian);
      printLorentzVector("chOSP4",  chOSP4,      false);
      printLorentzVector("chSS1P4", chSS1P4,     cartesian);
      printLorentzVector("chSS1P4", chSS1P4,     false);
      printLorentzVector("chSS2P4", chSS2P4,     cartesian);
      printLorentzVector("chSS2P4", chSS2P4,     false);
    }

    reco::Candidate::LorentzVector tauP4_trf   = getP4_ttrf_hf_trf(tauP4,   boost_ttrf, r, n, k, boost_trf);
    reco::Candidate::LorentzVector chOSP4_trf  = getP4_ttrf_hf_trf(chOSP4,  boost_ttrf, r, n, k, boost_trf);
    reco::Candidate::LorentzVector chSS1P4_trf = getP4_ttrf_hf_trf(chSS1P4, boost_ttrf, r, n, k, boost_trf);
    reco::Candidate::LorentzVector chSS2P4_trf = getP4_ttrf_hf_trf(chSS2P4, boost_ttrf, r, n, k, boost_trf);
    if ( verbosity >= 4 )
    {
      std::cout << "in tau restframe:\n";
      printLorentzVector("tauP4",   tauP4_trf,   cartesian);
      printLorentzVector("tauP4",   tauP4_trf,   false);
      printLorentzVector("chOSP4",  chOSP4_trf,  cartesian);
      printLorentzVector("chOSP4",  chOSP4_trf,  false);
      printLorentzVector("chSS1P4", chSS1P4_trf, cartesian);
      printLorentzVector("chSS1P4", chSS1P4_trf, false);
      printLorentzVector("chSS2P4", chSS2P4_trf, cartesian);
      printLorentzVector("chSS2P4", chSS2P4_trf, false);
    }

    reco::Candidate::Vector h = a1pol(chSS1P4_trf, chSS2P4_trf, chOSP4_trf, charge_sum, PolarimetricVectorTau2a1::k3ChargedPi);
    if ( verbosity >= 4 )
    { 
      printVector("h", h, cartesian);
      printVector("h", h, false);
    }
    return h.unit();
  }
}

reco::Candidate::Vector
PolarimetricVectorAlgoThreeProng0Pi0::operator()(const KinematicEvent& evt, int tau) const
{
  if ( verbosity_ >= 4 )
  {
    std::cout << "<PolarimetricVectorAlgoThreeProng0Pi0::operator()>:\n";
  }

  reco::Candidate::LorentzVector tauP4;
  const std::vector<KinematicParticle>* daughters = nullptr;
  if ( tau == pol::kTauPlus )
  {
    tauP4 = evt.tauPlusP4();
    daughters = &evt.daughtersTauPlus();
  }
  else if ( tau == pol::kTauMinus )
  {
    tauP4 = evt.tauMinusP4();
    daughters = &evt.daughtersTauMinus();
  }
  else assert(0);
  reco::Candidate::LorentzVector higgsP4 = evt.tauPlusP4() + evt.tauMinusP4();
  ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(higgsP4.BoostToCM());
  reco::Candidate::LorentzVector tauP4_ttrf = getP4_rf(tauP4, boost_ttrf);
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(evt.tauMinusP4(), &higgsP4, &boost_ttrf, hAxis_, collider_, r, n, k, verbosity_, cartesian_);
  reco::Candidate::LorentzVector tauP4_hf = getP4_hf(tauP4_ttrf, r, n, k);
  ROOT::Math::Boost boost_trf = ROOT::Math::Boost(tauP4_hf.BoostToCM());

  reco::Candidate::Vector h = getPolarimetricVec_ThreeProng0PiZero(a1pol_, tauP4, *daughters, boost_ttrf, r, n, k, boost_trf, verbosity_, cartesian_);
  return h;
}
