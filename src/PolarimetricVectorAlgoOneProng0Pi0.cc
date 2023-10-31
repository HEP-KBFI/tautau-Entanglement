#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoOneProng0Pi0.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // Candidate::LorentzVector, Candidate::Point

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/cube.h"                      // cube()
#include "TauAnalysis/Entanglement/interface/constants.h"                 // mTau, gamma_va
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/getP4_hf.h"                  // getP4_hf()
#include "TauAnalysis/Entanglement/interface/getP4_rf.h"                  // getP4_rf()
#include "TauAnalysis/Entanglement/interface/getP4_ttrf_hf_trf.h"         // getP4_ttrf_hf_trf()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"        // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printVector.h"               // printVector()
#include "TauAnalysis/Entanglement/interface/square.h"                    // square()

#include <Math/Boost.h>                                                   // Boost

PolarimetricVectorAlgoOneProng0Pi0::PolarimetricVectorAlgoOneProng0Pi0(const edm::ParameterSet& cfg)
  : PolarimetricVectorAlgoBase(cfg)
{}

PolarimetricVectorAlgoOneProng0Pi0::~PolarimetricVectorAlgoOneProng0Pi0()
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
    if ( verbosity >= 4 )
    { 
      std::cout << "<getPolarimetricVec_OneProng0PiZero>:\n";
    }

    const KinematicParticle* ch = nullptr;
    for ( const KinematicParticle& daughter : daughters )
    {
      if ( abs(daughter.pdgId()) == 211 || abs(daughter.pdgId()) == 321 )
      {
        ch = &daughter;
      }
    }
    if ( !ch )
      throw cmsException("getPolarimetricVec_OneProng0PiZero", __LINE__)
        << "Failed to find charged pion !!\n";

    // CV: notation of four-vectors chosen according to Section 3.3 of the paper
    //       Comput.Phys.Commun. 64 (1990) 275
    reco::Candidate::LorentzVector chP4 = ch->p4();
    if ( verbosity >= 4 )
    { 
      printLorentzVector("chP4", chP4, cartesian);
    }
    
    reco::Candidate::LorentzVector Q = getP4_ttrf_hf_trf(chP4, boost_ttrf, r, n, k, boost_trf);
    if ( verbosity >= 4 )
    { 
      printLorentzVector("Q", Q, cartesian);
    }
    const double f1 = 0.1284;
    double omega = (square(mTau) - chP4.mass2())*square(mTau);
    reco::Candidate::Vector h = -(2.*gamma_va*square(f1)*cube(mTau)/omega)*Q.Vect();
    if ( verbosity >= 4 )
    { 
      printVector("h", h, cartesian);
    }
    return h.unit();
  }
}

reco::Candidate::Vector
PolarimetricVectorAlgoOneProng0Pi0::operator()(const KinematicEvent& evt, int tau) const
{
  if ( verbosity_ >= 4 )
  {
    std::cout << "<PolarimetricVectorAlgoOneProng0Pi0::operator()>:\n";
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

  reco::Candidate::Vector h = getPolarimetricVec_OneProng0PiZero(tauP4, *daughters, boost_ttrf, r, n, k, boost_trf, verbosity_, cartesian_);
  return h;
}
