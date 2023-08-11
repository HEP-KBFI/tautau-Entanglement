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

#include <algorithm>                                                      // std::min()
#include <cmath>                                                          // std::atan(), std::cos(), std::sin(), std::sqrt()
#include <iostream>                                                       // std::cout

StartPosFinder2::StartPosFinder2(const edm::ParameterSet& cfg)
  : StartPosFinderBase(cfg)
  , resolutions_(nullptr)
  , applyRecoilEnergy_and_PzConstraint_(cfg.getParameter<bool>("applyRecoilEnergy_and_PzConstraint"))
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
             bool& errorFlag,
             int verbosity, bool cartesian = true)
  {
    if ( verbosity >= 3 )
    {
      std::cout << "<comp_tauP4>:\n";
    }

    reco::Candidate::Vector r, n, k; 
    get_localCoordinateSystem(visTauP4, nullptr, nullptr, kBeam, r, n, k);
    if ( verbosity >= 4 )
    {
      printVector("r", r, false);
      printVector("n", n, false);
      printVector("k", k, false);
    }

    // CV: assume angle between tau lepton and visible decay products to be negligible, 
    //     to obtain estimate for tau lepton four-vector, which is used to compute Gottfried-Jackson angle
    //    (this assumption is referred to as collinear approximation in the literature, cf. Section 5.2 of CMS AN-2011/165)
    double tauPx_CA = tauP*std::cos(visTauP4.phi())*std::sin(visTauP4.theta());
    double tauPy_CA = tauP*std::sin(visTauP4.phi())*std::sin(visTauP4.theta());
    double tauPz_CA = tauP*std::cos(visTauP4.theta());
    double tauE     = std::sqrt(square(tauP) + square(mTau));
    reco::Candidate::LorentzVector tauP4_CA(tauPx_CA, tauPy_CA, tauPz_CA, tauE);
    double cosThetaGJ = comp_cosThetaGJ_solution(tauP4_CA, visTauP4, &errorFlag);
    if ( errorFlag )
    {
      return reco::Candidate::LorentzVector();
    }

    auto tmp = tipPCA - pv;
    reco::Candidate::Vector flightlength(tmp.x(), tmp.y(), tmp.z()); 
    double flightlength_r = flightlength.Dot(r);
    double flightlength_n = flightlength.Dot(n);
   
    double alpha = std::atan(flightlength_n/flightlength_r);
    double beta  = std::acos(cosThetaGJ);
    if ( verbosity >= 3 )
    {
      std::cout << "alpha = " << alpha << ": cos(alpha) = " << cos(alpha) << ", sin(alpha) = " << sin(alpha) << "\n";
      std::cout << "beta = " << beta << ": cos(beta) = " << cos(beta) << ", sin(beta) = " << sin(beta) << "\n";
    }

    reco::Candidate::Vector tauP3 = tauP*(std::cos(alpha)*std::sin(beta)*r + std::sin(alpha)*std::sin(beta)*n + std::cos(beta)*k);
    reco::Candidate::LorentzVector tauP4(tauP3.x(), tauP3.y(), tauP3.z(), tauE);
    if ( verbosity >= 3 )
    {
      printLorentzVector("tauP4", tauP4, cartesian);
      if ( cartesian ) std::cout << " mass = " << tauP4.mass() << "\n";
    }

    return tauP4;
  }

  double
  comp_higgsMass(const reco::Candidate::LorentzVector& tauPlusP4,
                 const reco::Candidate::LorentzVector& visTauMinusP4, const reco::Candidate::Point& pv, const reco::Candidate::Point& tipPCATauMinus, double xMinus,
                 bool& errorFlag,
                 int verbosity = 0, bool cartesian = true)
  {
    double tauMinusP = visTauMinusP4.P()/xMinus;
    reco::Candidate::LorentzVector tauMinusP4 = comp_tauP4(tauMinusP, visTauMinusP4, pv, tipPCATauMinus, errorFlag, verbosity, cartesian);
    if ( errorFlag )
    {
      return 0.;
    }
    double higgsMass = (tauPlusP4 + tauMinusP4).mass();
    return higgsMass;
  }

  double
  comp_higgsMass_CA(const reco::Candidate::LorentzVector& tauPlusP4,
                    const reco::Candidate::LorentzVector& visTauMinusP4, double xMinus,
                    int verbosity = 0, bool cartesian = true)
  {
    double tauMinusPx = visTauMinusP4.px()/xMinus;
    double tauMinusPy = visTauMinusP4.py()/xMinus;
    double tauMinusPz = visTauMinusP4.pz()/xMinus;
    double tauMinusE  = std::sqrt(square(tauMinusPx) + square(tauMinusPy) + square(tauMinusPz) + square(mTau));
    reco::Candidate::LorentzVector tauMinusP4(tauMinusPx, tauMinusPy, tauMinusPz, tauMinusE);
    double higgsMass = (tauPlusP4 + tauMinusP4).mass();
    return higgsMass;
  }
  
  

  double
  comp_xMinus(const reco::Candidate::LorentzVector& tauPlusP4, const reco::Candidate::LorentzVector& visTauPlusP4, double xPlus,
              const reco::Candidate::LorentzVector& visTauMinusP4, const reco::Candidate::Point& pv, const reco::Candidate::Point& tipPCATauMinus,
              bool& errorFlag,
              int verbosity, bool cartesian = true)
  {
    if ( verbosity >= 3 )
    {
      std::cout << "<comp_xMinus>:\n";
      std::cout << " xPlus = " << xPlus << "\n";
    }

    double mVis = (visTauPlusP4 + visTauMinusP4).mass();

    double xMinus = square(mVis/mHiggs)/xPlus;
    double xMinus_lo = 1.e-6;
    double xMinus_hi = 1.;

    int iteration = 0;
    const int max_iterations = 1000;
    double xMinus_bestfit = -1.;
    bool hasConverged = false;    
    while ( !hasConverged && iteration < max_iterations )
    {
      if ( verbosity >= 3 )
      {
        std::cout << "iteration #" << iteration << ":\n";
      }

      bool higgsMass_errorFlag = false;
      double higgsMass = comp_higgsMass(tauPlusP4, visTauMinusP4, pv, tipPCATauMinus, xMinus, higgsMass_errorFlag, verbosity, cartesian);
      if ( higgsMass_errorFlag )
      {
        higgsMass = comp_higgsMass_CA(tauPlusP4, visTauMinusP4, xMinus);
      }
      // CV: high visible energy fraction means low Higgs mass and vice versa
      bool higgsMass_lo_errorFlag = false;
      double higgsMass_lo = comp_higgsMass(tauPlusP4, visTauMinusP4, pv, tipPCATauMinus, xMinus_hi, higgsMass_lo_errorFlag, verbosity, cartesian);
      if ( higgsMass_lo_errorFlag )
      {
        higgsMass_lo = comp_higgsMass_CA(tauPlusP4, visTauMinusP4, xMinus_hi);
      }
      bool higgsMass_hi_errorFlag = false;
      double higgsMass_hi = comp_higgsMass(tauPlusP4, visTauMinusP4, pv, tipPCATauMinus, xMinus_lo, higgsMass_hi_errorFlag, verbosity, cartesian);
      if ( higgsMass_hi_errorFlag )
      {
        higgsMass_hi = comp_higgsMass_CA(tauPlusP4, visTauMinusP4, xMinus_lo);
      }
      if ( verbosity >= 3 )
      {
        std::cout << "xMinus = " << xMinus << ": higgsMass = " << higgsMass << " (errorFlag = " << higgsMass_errorFlag << ")\n";
        std::cout << " lo = " << xMinus_lo << ": higgsMass = " << higgsMass_hi << "\n";
        std::cout << " hi = " << xMinus_hi << ": higgsMass = " << higgsMass_lo << "\n";
      }
      if ( !(higgsMass_lo <= mHiggs && higgsMass_hi >= mHiggs) )
      {
        if ( verbosity >= 3 )
        {
          std::cout << "Higgs mass not within low and high bounds --> returning error !!\n";
        }
        errorFlag = true;
        return 0.;
      }
    
      if ( std::fabs(higgsMass - mHiggs) < 1.e-2 )
      {
        if ( !higgsMass_errorFlag )
        {
          xMinus_bestfit = xMinus;
          hasConverged = true;
        }
        else
        {
          if ( verbosity >= 3 )
          {
            std::cout << "Higgs mass has converged, but is not a physical solution --> returning error !!\n";
          }
          errorFlag = true;
          return 0.;
        }
      }

      // CV: "Intervallhalbierungsverfahren";
      //     note again that high visible energy fraction means low Higgs mass and vice versa
      if ( higgsMass > mHiggs )
      {
        xMinus_lo = xMinus;
        xMinus = 0.5*(xMinus + xMinus_hi);
      }
      else
      {
        xMinus_hi = xMinus;
        xMinus = 0.5*(xMinus + xMinus_lo);
      }

      ++iteration;
    }
    assert(hasConverged);
    
    errorFlag = false;
    return xMinus_bestfit;
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
  const reco::Candidate::LorentzVector& visTauMinusP4 = kineEvt.visTauMinusP4();
  double mVis = (visTauPlusP4 + visTauMinusP4).mass();
  if ( verbosity_ >= 1 )
  {
    std::cout << "mVis = " << mVis << "\n";
  }

  const reco::Candidate::LorentzVector& recoilP4 = kineEvt.recoilP4();
  double recoilPx   = recoilP4.px();
  double recoilPy   = recoilP4.py();
  double recoilPz   = recoilP4.pz();
  double recoilM    = recoilP4.mass();
  double sigma_px   = resolutions_->recoilResolution_px();
  double sigma_py   = resolutions_->recoilResolution_py();
  double sigma_pz   = resolutions_->recoilResolution_pz();
  double sigma_mass = resolutions_->recoilResolution_mass();
  if ( verbosity_ >= 1 )
  {
    std::cout << "recoil:\n";
    std::cout << " Px = " << recoilPx << ", sigma = " << sigma_px << "\n";
    std::cout << " Py = " << recoilPy << ", sigma = " << sigma_py << "\n";
    if ( applyRecoilEnergy_and_PzConstraint_ )
    {
      std::cout << " Pz = " << recoilPz << ", sigma = " << sigma_pz << "\n";
      std::cout << " mass = " << recoilM << ", sigma = " << sigma_mass << "\n";
    }
  }

  double xPlus_bestfit = -1.;
  reco::Candidate::LorentzVector tauPlusP4_bestfit;
  double xMinus_bestfit = -1.;
  reco::Candidate::LorentzVector tauMinusP4_bestfit;
  double min_chi2 = -1.;
  const int numSteps = ( verbosity_ >= 2 ) ? 100 : 10000;
  double min_xPlus = square(mVis/mHiggs);
  double max_xPlus = 1.;
  for ( int idxStep = 0; idxStep < numSteps; ++idxStep )
  {
    if ( verbosity_ >= 3 )
    {
      std::cout << "step #" << idxStep << ":\n";
    }

    double xPlus = min_xPlus + idxStep*(max_xPlus - min_xPlus)/numSteps;    
    assert(xPlus >= 0. && xPlus <= 1.);
    double tauPlusP = visTauPlusP4.P()/xPlus;
    bool tauPlus_errorFlag = false;
    reco::Candidate::LorentzVector tauPlusP4 = comp_tauP4(tauPlusP, visTauPlusP4, kineEvt.pv(), kineEvt.tipPCATauPlus(), tauPlus_errorFlag, verbosity_, cartesian_);
    if ( tauPlus_errorFlag )
    {
      if ( verbosity_ >= 1 )
      {
        std::cout << "No physical solution exists for Gottfried-Jackson angle of tau+ -> skipping !!\n";
      }
      continue;
    }
    
    bool tauMinus_errorFlag = false;
    double xMinus = comp_xMinus(tauPlusP4, visTauPlusP4, xPlus, visTauMinusP4, kineEvt.pv(), kineEvt.tipPCATauMinus(), tauMinus_errorFlag, verbosity_, cartesian_);
    reco::Candidate::LorentzVector tauMinusP4;
    if ( !tauMinus_errorFlag )
    {
      assert(xMinus >= 0. && xMinus <= 1.);
      double tauMinusP = visTauMinusP4.P()/xMinus;
      tauMinusP4 = comp_tauP4(tauMinusP, visTauMinusP4, kineEvt.pv(), kineEvt.tipPCATauMinus(), tauMinus_errorFlag, verbosity_, cartesian_);
    }
    if ( tauMinus_errorFlag )
    {
      if ( verbosity_ >= 1 )
      {
        std::cout << "No physical solution exists for Gottfried-Jackson angle of tau- -> skipping !!\n";
      }
      continue;
    }
    
    double chi2 = square((tauPlusP4.px() + tauMinusP4.px() - recoilPx)/sigma_px) + square((tauPlusP4.py() + tauMinusP4.py() - recoilPy)/sigma_py);
    if ( applyRecoilEnergy_and_PzConstraint_ )
    {
      chi2 += square((tauPlusP4.pz() + tauMinusP4.pz() - recoilPz)/sigma_pz) + square(((tauPlusP4 + tauMinusP4).mass() - recoilM)/sigma_mass);
    }
    if ( verbosity_ >= 3 )
    {
      std::cout << "xPlus = " << xPlus << ", xMinus = " << xMinus << "\n";
      std::cout << "tauPlus + tauMinus:\n";
      std::cout << " Px = " << tauPlusP4.px() + tauMinusP4.px() << "\n";
      std::cout << " Py = " << tauPlusP4.py() + tauMinusP4.py() << "\n";
      if ( applyRecoilEnergy_and_PzConstraint_ )
      {
        std::cout << " Pz = " << tauPlusP4.pz() + tauMinusP4.pz() << "\n";
        std::cout << " mass = " << (tauPlusP4 + tauMinusP4).mass() << "\n";
      }
      std::cout << "chi^2 = " << chi2 << "\n";
      std::cout << "higgsMass = " << (tauPlusP4 + tauMinusP4).mass() << "\n";
    }
   
    if ( min_chi2 == -1. || chi2 < min_chi2 )
    {
      xPlus_bestfit = xPlus;
      tauPlusP4_bestfit = tauPlusP4;
      xMinus_bestfit = xMinus;
      tauMinusP4_bestfit = tauMinusP4;
      min_chi2 = chi2;
    }
  }
  if ( verbosity_ >= 1 )
  {
    std::cout << "best fit:\n";
    std::cout << "xPlus = " << xPlus_bestfit << ", xMinus = " << xMinus_bestfit << "\n";
    std::cout << "tauPlus + tauMinus:\n";
    std::cout << " Px = " << tauPlusP4_bestfit.px() + tauMinusP4_bestfit.px() << "\n";
    std::cout << " Py = " << tauPlusP4_bestfit.py() + tauMinusP4_bestfit.py() << "\n";
    if ( applyRecoilEnergy_and_PzConstraint_ )
    {
      std::cout << " Pz = " << tauPlusP4_bestfit.pz() + tauMinusP4_bestfit.pz() << "\n";
      std::cout << " mass = " << (tauPlusP4_bestfit + tauMinusP4_bestfit).mass() << "\n";
    }
    std::cout << "chi^2 = " << min_chi2 << "\n";
    printLorentzVector("tauPlusP4", tauPlusP4_bestfit, cartesian_);
    printLorentzVector("tauMinusP4", tauMinusP4_bestfit, cartesian_);
    std::cout << "higgsMass = " << (tauPlusP4_bestfit + tauMinusP4_bestfit).mass() << "\n";
  }

  KinematicEvent kineEvt_startpos = kineEvt;
  
  kineEvt_startpos.tauPlusP4_ = tauPlusP4_bestfit;
  kineEvt_startpos.tauPlusP4_isValid_ = true;
  //kineEvt_startpos.nuTauPlusP4_ = fixNuMass(tauPlusP4_bestfit - visTauPlusP4);
  kineEvt_startpos.nuTauPlusP4_ = tauPlusP4_bestfit - visTauPlusP4;
  kineEvt_startpos.nuTauPlusP4_isValid_ = true;
  const KinematicParticle* tauPlus_leadTrack = get_leadTrack(kineEvt_startpos.daughtersTauPlus());
  assert(tauPlus_leadTrack);
  const reco::Candidate::Point& tauPlus_tipPCA = kineEvt_startpos.tipPCATauPlus();
  kineEvt_startpos.svTauPlus_ = comp_PCA_line2line(kineEvt_startpos.pv(), tauPlusP4_bestfit, tauPlus_tipPCA, tauPlus_leadTrack->p4(), verbosity_);
  kineEvt_startpos.svTauPlus_isValid_ = true;

  kineEvt_startpos.tauMinusP4_ = tauMinusP4_bestfit;
  kineEvt_startpos.tauMinusP4_isValid_ = true;
  //kineEvt_startpos.nuTauMinusP4_ = fixNuMass(tauMinusP4_bestfit - visTauMinusP4);
  kineEvt_startpos.nuTauMinusP4_ = tauMinusP4_bestfit - visTauMinusP4;
  kineEvt_startpos.nuTauMinusP4_isValid_ = true;
  const KinematicParticle* tauMinus_leadTrack = get_leadTrack(kineEvt_startpos.daughtersTauMinus());
  assert(tauMinus_leadTrack);
  const reco::Candidate::Point& tauMinus_tipPCA = kineEvt_startpos.tipPCATauMinus();
  kineEvt_startpos.svTauMinus_ = comp_PCA_line2line(kineEvt_startpos.pv(), tauMinusP4_bestfit, tauMinus_tipPCA, tauMinus_leadTrack->p4(), verbosity_);
  kineEvt_startpos.svTauMinus_isValid_ = true;

  return kineEvt_startpos;
}
