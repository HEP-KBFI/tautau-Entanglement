#include "TauAnalysis/Entanglement/interface/StartPosFinder2.h"

#include "TauAnalysis/Entanglement/interface/comp_cosThetaGJ.h"           // comp_cosThetaGJ_solution()
#include "TauAnalysis/Entanglement/interface/comp_PCA_line2line.h"        // comp_PCA_line2line()
#include "TauAnalysis/Entanglement/interface/constants.h"                 // mHiggs, mTau
#include "TauAnalysis/Entanglement/interface/fixMass.h"                   // fixNuMass(), fixTauMass()
#include "TauAnalysis/Entanglement/interface/get_leadTrack.h"             // get_leadTrack()
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"        // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printVector.h"               // printVector()
#include "TauAnalysis/Entanglement/interface/square.h"                    // square()

#include <cmath>                                                          // std::atan(), std::cos(), std::sin(), std::sqrt()
#include <iostream>                                                       // std::cout

StartPosFinder2::StartPosFinder2(const edm::ParameterSet& cfg)
  : StartPosFinderBase(cfg)
  , resolutions_(nullptr)
  , applyHiggsMassConstraint_(cfg.getParameter<bool>("applyHiggsMassConstraint"))
{
  edm::ParameterSet cfg_resolutions = cfg.getParameterSet("resolutions");
  resolutions_ = new Resolutions(cfg_resolutions);
}

StartPosFinder2::~StartPosFinder2()
{
  delete resolutions_;
}

namespace
{
  reco::Candidate::LorentzVector
  comp_tauP4(double tauP, 
             const reco::Candidate::LorentzVector& visTauP4, 
             const reco::Candidate::Point& pv, const reco::Candidate::Point& tipPCA, 
             int verbosity, bool cartesian = true)
  {
    if ( verbosity >= 1 )
    {
      std::cout << "<comp_tauP4>:\n";
    }

    reco::Candidate::Vector r, n, k; 
    get_localCoordinateSystem(visTauP4, nullptr, nullptr, kBeam, r, n, k);
    if ( verbosity >= 3 )
    {
      printVector("r", r, false);
      printVector("n", n, false);
      printVector("k", k, false);
    }

    double tauPx_CA = tauP*std::cos(visTauP4.phi())*std::sin(visTauP4.theta());
    double tauPy_CA = tauP*std::sin(visTauP4.phi())*std::sin(visTauP4.theta());
    double tauPz_CA = tauP*std::cos(visTauP4.theta());
    double tauE     = std::sqrt(square(tauP) + square(mTau));
    double cosThetaGJ = comp_cosThetaGJ_solution(reco::Candidate::LorentzVector(tauPx_CA, tauPy_CA, tauPz_CA, tauE), visTauP4);

    auto tmp = tipPCA - pv;
    reco::Candidate::Vector flightlength(tmp.x(), tmp.y(), tmp.z()); 
    double flightlength_r = flightlength.Dot(r);
    double flightlength_n = flightlength.Dot(n);
   
    double alpha = std::atan(flightlength_n/flightlength_r);
    double beta  = std::acos(cosThetaGJ);
    if ( verbosity >= 1 )
    {
      std::cout << "alpha = " << alpha << ": cos(alpha) = " << cos(alpha) << ", sin(alpha) = " << sin(alpha) << "\n";
      std::cout << "beta = " << beta << ": cos(beta) = " << cos(beta) << ", sin(beta) = " << sin(beta) << "\n";
    }

    reco::Candidate::Vector tauP3 = tauP*(std::cos(alpha)*std::sin(beta)*r + std::sin(alpha)*std::sin(beta)*n + std::cos(beta)*k);
    reco::Candidate::LorentzVector tauP4(tauP3.x(), tauP3.y(), tauP3.z(), tauE);
    if ( verbosity >= 1 )
    {
      printLorentzVector("tauP4", tauP4, cartesian);
      if ( cartesian ) std::cout << " mass = " << tauP4.mass() << "\n";
    }

    return tauP4;
  }
}

KinematicEvent
StartPosFinder2::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<StartPosFinder2::operator()>:\n";
  }

  const reco::Candidate::LorentzVector& visTauPlusP4 = kineEvt.visTauPlusP4();
  double visTauPlusPx = visTauPlusP4.px();
  double visTauPlusPy = visTauPlusP4.py();
  const reco::Candidate::LorentzVector& visTauMinusP4 = kineEvt.visTauMinusP4();
  double visTauMinusPx = visTauMinusP4.px();
  double visTauMinusPy = visTauMinusP4.py();
  double mVis = (visTauPlusP4 + visTauMinusP4).mass();

  const reco::Candidate::LorentzVector& recoilP4 = kineEvt.recoilP4();
  double recoilPx = recoilP4.px();
  double recoilPy = recoilP4.py();
  double sigma_px = resolutions_->recoilResolution_px();
  double sigma_py = resolutions_->recoilResolution_px();

  double xPlus_bestfit = -1.;
  double xMinus_bestfit = -1.;
  double min_chi2 = -1.;
  const int numSteps = 10000;
  for ( int idxStep = 0; idxStep < numSteps; ++idxStep )
  {
    double xPlus  = idxStep*1./numSteps;
    double tauPlusPx = visTauPlusPx/xPlus;
    double tauPlusPy = visTauPlusPy/xPlus;
    double xMinus = square(mVis/mHiggs)/xPlus;
    double tauMinusPx = visTauMinusPx/xMinus;
    double tauMinusPy = visTauMinusPy/xMinus;
    double chi2 = square((tauPlusPx + tauMinusPx - recoilPx)/sigma_px) + square((tauPlusPy + tauMinusPy - recoilPy)/sigma_py); 
    if ( min_chi2 == -1. || chi2 < min_chi2 )
    {
      xPlus_bestfit = xPlus;
      xMinus_bestfit = xMinus;
      min_chi2 = chi2;
    }
  }

  reco::Candidate::LorentzVector tauPlusP4  = comp_tauP4(visTauPlusP4.P()/xPlus_bestfit,   visTauPlusP4,  kineEvt.pv(), kineEvt.tipPCATauPlus(), verbosity_, cartesian_);
  reco::Candidate::LorentzVector tauMinusP4 = comp_tauP4(visTauMinusP4.P()/xMinus_bestfit, visTauMinusP4, kineEvt.pv(), kineEvt.tipPCATauMinus(), verbosity_, cartesian_);

  if ( applyHiggsMassConstraint_ )
  {
    reco::Candidate::LorentzVector higgsP4 = tauPlusP4 + tauMinusP4;
    double sf = mHiggs/higgsP4.mass();
    double tauPlusPx = sf*tauPlusP4.px();
    double tauPlusPy = sf*tauPlusP4.py();
    double tauPlusPz = sf*tauPlusP4.pz();
    double tauPlusE  = std::sqrt(square(tauPlusPx) + square(tauPlusPy) + square(tauPlusPz) + square(mTau));
    tauPlusP4 = reco::Candidate::LorentzVector(tauPlusPx, tauPlusPy, tauPlusPz, tauPlusE);
    double tauMinusPx = sf*tauMinusP4.px();
    double tauMinusPy = sf*tauMinusP4.py();
    double tauMinusPz = sf*tauMinusP4.pz();
    double tauMinusE  = std::sqrt(square(tauMinusPx) + square(tauMinusPy) + square(tauMinusPz) + square(mTau));
    tauMinusP4 = reco::Candidate::LorentzVector(tauMinusPx, tauMinusPy, tauMinusPz, tauMinusE);
  }

  KinematicEvent kineEvt_startpos = kineEvt;
  
  kineEvt_startpos.tauPlusP4_ = tauPlusP4;
  kineEvt_startpos.tauPlusP4_isValid_ = true;
  //kineEvt_startpos.nuTauPlusP4_ = fixNuMass(tauPlusP4 - visTauPlusP4);
  kineEvt_startpos.nuTauPlusP4_ = tauPlusP4 - visTauPlusP4;
  kineEvt_startpos.nuTauPlusP4_isValid_ = true;
  if ( !kineEvt_startpos.svTauPlus_isValid() )
  {
    const KinematicParticle* leadTrack = get_leadTrack(kineEvt_startpos.daughtersTauPlus());
    assert(leadTrack);
    const reco::Candidate::Point& tipPCA = kineEvt_startpos.tipPCATauPlus();
    kineEvt_startpos.svTauPlus_ = comp_PCA_line2line(kineEvt_startpos.pv(), tauPlusP4, tipPCA, leadTrack->p4(), verbosity_);
    kineEvt_startpos.svTauPlus_isValid_ = true;
  }

  kineEvt_startpos.tauMinusP4_ = tauMinusP4;
  kineEvt_startpos.tauMinusP4_isValid_ = true;
  //kineEvt_startpos.nuTauMinusP4_ = fixNuMass(tauMinusP4 - visTauMinusP4);
  kineEvt_startpos.nuTauMinusP4_ = tauMinusP4 - visTauMinusP4;
  kineEvt_startpos.nuTauMinusP4_isValid_ = true;
  if ( !kineEvt_startpos.svTauMinus_isValid() )
  {
    const KinematicParticle* leadTrack = get_leadTrack(kineEvt_startpos.daughtersTauMinus());
    assert(leadTrack);
    const reco::Candidate::Point& tipPCA = kineEvt_startpos.tipPCATauMinus();
    kineEvt_startpos.svTauMinus_ = comp_PCA_line2line(kineEvt_startpos.pv(), tauMinusP4, tipPCA, leadTrack->p4(), verbosity_);
    kineEvt_startpos.svTauMinus_isValid_ = true;
  }

  return kineEvt_startpos;
}
