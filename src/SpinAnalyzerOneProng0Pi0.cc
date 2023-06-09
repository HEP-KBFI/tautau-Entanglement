#include "TauAnalysis/Entanglement/interface/SpinAnalyzerOneProng0Pi0.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // Candidate::LorentzVector, Candidate::Point

#include "TauAnalysis/Entanglement/interface/auxFunctions.h"              // getP4_rf(), getP4_hf(), getP4_ttrf_hf_trf(), printVector()
#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"                 // mTau, gamma_va
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()

#include <Math/Boost.h>                                                   // Boost

SpinAnalyzerOneProng0Pi0::SpinAnalyzerOneProng0Pi0(const edm::ParameterSet& cfg)
  : SpinAnalyzer(cfg)
{}

SpinAnalyzerOneProng0Pi0::~SpinAnalyzerOneProng0Pi0()
{}

namespace
{
  reco::Candidate::Vector
  getPolarimetricVec_OneProng0PiZero(const reco::Candidate::LorentzVector& tauP4,
                                     const std::vector<KinematicParticle>& daughters,
                                     const ROOT::Math::Boost& boost_ttrf,
                                     const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                                     const ROOT::Math::Boost& boost_trf,
                                     int verbosity = 0, bool cartesian = true)
  {
    if ( verbosity >= 1 )
    { 
      std::cout << "<getPolarimetricVec_OneProng0PiZero>:\n";
    }

    const KinematicParticle* ch = nullptr;
    for ( const KinematicParticle& daughter : daughters )
    {
      if ( abs(daughter.get_pdgId()) == 211 || abs(daughter.get_pdgId()) == 321 )
      {
        ch = &daughter;
      }
    }
    if ( !ch )
      throw cmsException("getPolarimetricVec_OneProng0PiZero", __LINE__)
        << "Failed to find charged pion !!\n";

    // CV: notation of four-vectors chosen according to Section 3.3 of the paper
    //       Comput.Phys.Commun. 64 (1990) 275
    reco::Candidate::LorentzVector chP4 = ch->get_p4();
    if ( verbosity >= 1 )
    { 
      printLorentzVector("chP4", chP4, cartesian);
    }
    
    reco::Candidate::LorentzVector Q = getP4_ttrf_hf_trf(chP4, boost_ttrf, r, n, k, boost_trf);
    if ( verbosity >= 1 )
    { 
      printLorentzVector("Q", Q, cartesian);
    }
    const double f1 = 0.1284;
    double omega = (square(mTau) - chP4.mass2())*square(mTau);
    reco::Candidate::Vector h = -(2.*gamma_va*square(f1)*cube(mTau)/omega)*Q.Vect();
    if ( verbosity >= 1 )
    { 
      printVector("h", h, cartesian);
    }
    return h.unit();
  }
}

reco::Candidate::Vector
SpinAnalyzerOneProng0Pi0::operator()(const KinematicEvent& evt, int tau)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<SpinAnalyzerOneProng0Pi0::operator()>:\n";
  }

  reco::Candidate::LorentzVector tauP4;
  const std::vector<KinematicParticle>* daughters = nullptr;
  if ( tau == kTauPlus )
  {
    tauP4 = evt.get_tauPlusP4();
    daughters = &evt.get_daughtersTauPlus();
  }
  else if ( tau == kTauMinus )
  {
    tauP4 = evt.get_tauMinusP4();
    daughters = &evt.get_daughtersTauMinus();
  }
  else assert(0);
  reco::Candidate::LorentzVector recoilP4 = evt.get_recoilP4();
  ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(recoilP4.BoostToCM());
  reco::Candidate::LorentzVector tauP4_ttrf = getP4_rf(tauP4, boost_ttrf);
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(evt.get_tauMinusP4(), &recoilP4, &boost_ttrf, hAxis_, r, n, k, verbosity_, cartesian_);
  reco::Candidate::LorentzVector tauP4_hf = getP4_hf(tauP4_ttrf, r, n, k);
  ROOT::Math::Boost boost_trf = ROOT::Math::Boost(tauP4_hf.BoostToCM());

  reco::Candidate::Vector h = getPolarimetricVec_OneProng0PiZero(tauP4, *daughters, boost_ttrf, r, n, k, boost_trf, verbosity_, cartesian_);
  return h;
}
