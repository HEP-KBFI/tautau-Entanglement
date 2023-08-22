#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoOneProng1Pi0.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // Candidate::LorentzVector, Candidate::Point

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"                 // mTau, gamma_va
#include "TauAnalysis/Entanglement/interface/fixMass.h"                   // fixTauMass()
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/getP4_hf.h"                  // getP4_hf()
#include "TauAnalysis/Entanglement/interface/getP4_rf.h"                  // getP4_rf()
#include "TauAnalysis/Entanglement/interface/getP4_ttrf_hf_trf.h"         // getP4_ttrf_hf_trf()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"        // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printVector.h"               // printVector()

#include <Math/Boost.h>                                                   // Boost

PolarimetricVectorAlgoOneProng1Pi0::PolarimetricVectorAlgoOneProng1Pi0(const edm::ParameterSet& cfg)
  : PolarimetricVectorAlgoBase(cfg)
{}

PolarimetricVectorAlgoOneProng1Pi0::~PolarimetricVectorAlgoOneProng1Pi0()
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
    if ( verbosity >= 2 )
    {
      std::cout << "<getPolarimetricVec_OneProng1PiZero>:\n";
    }

    const KinematicParticle* ch = nullptr;
    const KinematicParticle* pi0 = nullptr;
    for ( const KinematicParticle& daughter : daughters )
    {
      if ( abs(daughter.pdgId()) == 211 || abs(daughter.pdgId()) == 321 )
      {
        ch = &daughter;
      }
      if ( daughter.pdgId() == 111 )
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
    reco::Candidate::LorentzVector chP4 = ch->p4();
    reco::Candidate::LorentzVector q1 = getP4_ttrf_hf_trf(chP4, boost_ttrf, r, n, k, boost_trf);
    if ( verbosity >= 2 )
    { 
      printLorentzVector("chP4", chP4, cartesian);
      printLorentzVector("q1", q1, cartesian);
    }

    reco::Candidate::LorentzVector pi0P4 = pi0->p4();
    reco::Candidate::LorentzVector q2 = getP4_ttrf_hf_trf(pi0P4, boost_ttrf, r, n, k, boost_trf);
    if ( verbosity >= 2 )
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
    if ( verbosity >= 2 )
    {
      printLorentzVector("nuP4", nuP4, cartesian);
      std::cout << " mass = " << nuP4.mass() << "\n";
      printLorentzVector("N", N, cartesian);
      std::cout << " mass = " << N.mass() << "\n";
    }
    assert(nuP4.energy() >= 0. && N.energy() >= 0.);

    reco::Candidate::LorentzVector P = getP4_ttrf_hf_trf(tauP4, boost_ttrf, r, n, k, boost_trf);
    if ( verbosity >= 2 )
    {
      printLorentzVector("P", P, cartesian);
    }

    reco::Candidate::LorentzVector q = q1 - q2;
    double omega = 2.*(q.Dot(N))*(q.Dot(P)) - q.mass2()*(N.Dot(P));
    // CV: term 2.*|f2|^2 appears in expression for h as well as in expression for omega
    //     and drops out
    reco::Candidate::Vector h = -(gamma_va*mTau/omega)*(2.*(q.Dot(N))*q.Vect() - q.mass2()*N.Vect());
    if ( verbosity >= 2 )
    { 
      printVector("h", h, cartesian);
    }
    return h.unit();
  }
}

reco::Candidate::Vector
PolarimetricVectorAlgoOneProng1Pi0::operator()(const KinematicEvent& evt, int tau)
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorAlgoOneProng1Pi0::operator()>:\n";
  }

  reco::Candidate::LorentzVector tauP4;
  const std::vector<KinematicParticle>* daughters = nullptr;
  reco::Candidate::LorentzVector visTauP4;
  if ( tau == pol::kTauPlus )
  {
    tauP4 = evt.tauPlusP4();
    daughters = &evt.daughtersTauPlus();
    visTauP4 = evt.visTauPlusP4();
  }
  else if ( tau == pol::kTauMinus )
  {
    tauP4 = evt.tauMinusP4();
    daughters = &evt.daughtersTauMinus();
    visTauP4 = evt.visTauMinusP4();
  }
  else assert(0);
  if ( verbosity_ >= 2 )
  {
    printLorentzVector("tauP4", tauP4);
    std::cout << " mass = " << tauP4.mass() << "\n";
    printLorentzVector("visTauP4", visTauP4);
    std::cout << " mass = " << visTauP4.mass() << "\n";
  } 
  reco::Candidate::LorentzVector higgsP4 = evt.tauPlusP4() + evt.tauMinusP4();
  if ( verbosity_ >= 2 )
  {
    printLorentzVector("higgsP4", higgsP4);
    std::cout << " mass = " << higgsP4.mass() << "\n";
  }
  ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(higgsP4.BoostToCM());
  reco::Candidate::LorentzVector tauP4_ttrf = fixTauMass(getP4_rf(tauP4, boost_ttrf));
  if ( verbosity_ >= 2 )
  {
    printLorentzVector("tauP4_ttrf", tauP4_ttrf);
    std::cout << " mass = " << tauP4_ttrf.mass() << "\n";
  }
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(evt.tauMinusP4(), &higgsP4, &boost_ttrf, hAxis_, collider_, r, n, k, verbosity_, cartesian_);
  reco::Candidate::LorentzVector tauP4_hf = fixTauMass(getP4_hf(tauP4_ttrf, r, n, k));
  if ( verbosity_ >= 2 )
  {
    printLorentzVector("tauP4_hf", tauP4_hf);
    std::cout << " mass = " << tauP4_hf.mass() << "\n";
  }
  ROOT::Math::Boost boost_trf = ROOT::Math::Boost(tauP4_hf.BoostToCM());

  reco::Candidate::Vector h = getPolarimetricVec_OneProng1PiZero(tauP4, *daughters, visTauP4, boost_ttrf, r, n, k, boost_trf, verbosity_, cartesian_);
assert(getP4_hf(tauP4_ttrf, r, n, k).mass() >= 0.);
  return h;
}
