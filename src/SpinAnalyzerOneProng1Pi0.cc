#include "TauAnalysis/Entanglement/interface/SpinAnalyzerOneProng1Pi0.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // Candidate::LorentzVector, Candidate::Point

#include "TauAnalysis/Entanglement/interface/auxFunctions.h"              // getP4_rf(), getP4_hf(), getP4_ttrf_hf_trf(), printVector()
#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"                 // mTau, gamma_va
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()

#include <Math/Boost.h>                                                   // Boost

SpinAnalyzerOneProng1Pi0::SpinAnalyzerOneProng1Pi0(const edm::ParameterSet& cfg)
  : SpinAnalyzer(cfg)
{}

SpinAnalyzerOneProng1Pi0::~SpinAnalyzerOneProng1Pi0()
{}

namespace
{
  reco::Candidate::LorentzVector
  fixNeutrinoMass(const reco::Candidate::LorentzVector& nuP4)
  {
    double nuPx = nuP4.px();
    double nuPy = nuP4.py();
    double nuPz = nuP4.pz();
    double nuE  = std::sqrt(nuPx*nuPx + nuPy*nuPy + nuPz*nuPz);
    reco::Candidate::LorentzVector nuP4_fixed(nuPx, nuPy, nuPz, nuE);
    return nuP4_fixed;
  }

  reco::Candidate::Vector
  getPolarimetricVec_OneProng1PiZero(const reco::Candidate::LorentzVector& tauP4,
                                     const std::vector<KinematicParticle>& daughters,
                                     const reco::Candidate::LorentzVector& visTauP4,
                                     const ROOT::Math::Boost& boost_ttrf,
                                     const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                                     const ROOT::Math::Boost& boost_trf,
                                     int verbosity = 0, bool cartesian = true)
  {
    if ( verbosity >= 1 )
    {
      std::cout << "<getPolarimetricVec_OneProng1PiZero>:\n";
    }

    const KinematicParticle* ch = nullptr;
    const KinematicParticle* pi0 = nullptr;
    for ( const KinematicParticle& daughter : daughters )
    {
      if ( abs(daughter.get_pdgId()) == 211 || abs(daughter.get_pdgId()) == 321 )
      {
        ch = &daughter;
      }
      if ( daughter.get_pdgId() == 111 )
      {
        pi0 = &daughter;
      }
    }
    if ( !ch )
      throw cmsException("getPolarimetricVec_OneProng0PiZero", __LINE__)
        << "Failed to find charged pion !!\n";
    if ( !pi0 )
      throw cmsException("getPolarimetricVec_OneProng1PiZero", __LINE__)
        << "Failed to find neutral pion !!\n";

    // CV: notation of four-vectors chosen according to Section 3.4 of the paper
    //       Comput.Phys.Commun. 64 (1990) 275
    reco::Candidate::LorentzVector chP4 = ch->get_p4();
    reco::Candidate::LorentzVector q1 = getP4_ttrf_hf_trf(chP4, boost_ttrf, r, n, k, boost_trf);
    if ( verbosity >= 1 )
    { 
      printLorentzVector("chP4", chP4, cartesian);
      printLorentzVector("q1", q1, cartesian);
    }

    reco::Candidate::LorentzVector pi0P4 = pi0->get_p4();
    reco::Candidate::LorentzVector q2 = getP4_ttrf_hf_trf(pi0P4, boost_ttrf, r, n, k, boost_trf);
    if ( verbosity >= 1 )
    {
      printLorentzVector("pi0P4", pi0P4, cartesian);
      printLorentzVector("q2", q2, cartesian);
    }

    // CV: the neutrino four-vector computed by taking the difference 
    //     between the tau four-vector and the four-vector of the visible tau decay products
    //     yields mass values that a few GeV off, presumably due to rounding errors.
    //     Adjust the energy component of the neutrino four-vector such that its mass value equals zero,
    //     while keeping the Px, Py, Pz momentum components fixed
    reco::Candidate::LorentzVector nuP4 = fixNeutrinoMass(tauP4 - visTauP4);
    reco::Candidate::LorentzVector N = fixNeutrinoMass(getP4_ttrf_hf_trf(nuP4, boost_ttrf, r, n, k, boost_trf));
    if ( verbosity >= 1 )
    {
      printLorentzVector("nuP4", nuP4, cartesian);
      std::cout << " mass = " << nuP4.mass() << "\n";
      printLorentzVector("N", N, cartesian);
      std::cout << " mass = " << N.mass() << "\n";
    }
    assert(nuP4.energy() >= 0. && N.energy() >= 0.);

    reco::Candidate::LorentzVector P = getP4_ttrf_hf_trf(tauP4, boost_ttrf, r, n, k, boost_trf);
    if ( verbosity >= 1 )
    {
      printLorentzVector("P", P, cartesian);
    }

    reco::Candidate::LorentzVector q = q1 - q2;
    double omega = 2.*(q.Dot(N))*(q.Dot(P)) - q.mass2()*(N.Dot(P));
    // CV: term 2.*|f2|^2 appears in expression for h as well as in expression for omega
    //     and drops out
    reco::Candidate::Vector h = -(gamma_va*mTau/omega)*(2.*(q.Dot(N))*q.Vect() - q.mass2()*N.Vect());
    if ( verbosity >= 1 )
    { 
      printVector("h", h, cartesian);
    }
    return h.unit();
  }
}

reco::Candidate::Vector
SpinAnalyzerOneProng1Pi0::operator()(const KinematicEvent& evt, int tau)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<SpinAnalyzerOneProng1Pi0::operator()>:\n";
  }

  reco::Candidate::LorentzVector tauP4;
  const std::vector<KinematicParticle>* daughters = nullptr;
  reco::Candidate::LorentzVector visTauP4;
  if ( tau == kTauPlus )
  {
    tauP4 = evt.get_tauPlusP4();
    daughters = &evt.get_daughtersTauPlus();
    visTauP4 = evt.get_visTauPlusP4();
  }
  else if ( tau == kTauMinus )
  {
    tauP4 = evt.get_tauMinusP4();
    daughters = &evt.get_daughtersTauMinus();
    visTauP4 = evt.get_visTauMinusP4();
  }
  else assert(0);
  reco::Candidate::LorentzVector recoilP4 = evt.get_recoilP4();
  ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(recoilP4.BoostToCM());
  reco::Candidate::LorentzVector tauP4_ttrf = getP4_rf(tauP4, boost_ttrf);
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(evt.get_tauMinusP4(), &recoilP4, &boost_ttrf, hAxis_, r, n, k, verbosity_, cartesian_);
  reco::Candidate::LorentzVector tauP4_hf = getP4_hf(tauP4_ttrf, r, n, k);
  ROOT::Math::Boost boost_trf = ROOT::Math::Boost(tauP4_hf.BoostToCM());

  reco::Candidate::Vector h = getPolarimetricVec_OneProng1PiZero(tauP4, *daughters, visTauP4, boost_ttrf, r, n, k, boost_trf, verbosity_, cartesian_);
  return h;
}
