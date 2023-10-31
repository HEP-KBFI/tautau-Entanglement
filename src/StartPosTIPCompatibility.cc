
#include "TauAnalysis/Entanglement/interface/StartPosTIPCompatibility.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#include "DataFormats/TauReco/interface/PFTau.h"                          // reco::PFTau::hadronicDecayMode
#pragma GCC diagnostic pop

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/comp_PCA_line2line.h"        // comp_PCA_line2line()
#include "TauAnalysis/Entanglement/interface/get_leadTrack.h"             // get_leadTrack()
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_n(), get_r()
#include "TauAnalysis/Entanglement/interface/getP4_rf.h"                  // getP4_rf()
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h"         // math::Matrix3x3, math::Vector3
#include "TauAnalysis/Entanglement/interface/printCovMatrix.h"            // printCovMatrix()
#include "TauAnalysis/Entanglement/interface/printDistance.h"             // printDistance()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"        // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printVector.h"               // printVector()

#include <cmath>                                                          // std::fabs(), std::sqrt()

StartPosTIPCompatibility::StartPosTIPCompatibility(const edm::ParameterSet& cfg)
  : verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{}

StartPosTIPCompatibility::~StartPosTIPCompatibility()
{}
   
namespace
{
  double
  comp_tipCompatibilityThreeProng(const reco::Candidate::LorentzVector& tauP4,
                                  const reco::Candidate::Point& pv, const math::Matrix3x3& pvCov,
                                  const reco::Candidate::Point& sv, const math::Matrix3x3& svCov,
                                  int verbosity = -1, bool cartesian = true)
  {
    if ( verbosity >= 3 )
    {
      std::cout << "<comp_tipCompatibilityThreeProng>:\n";
      printLorentzVector("tauP4", tauP4, cartesian);
    }

    double flightlength = std::sqrt((sv - pv).mag2());
    auto sv_exp = pv + flightlength*tauP4.Vect().unit();
    auto residual = sv - sv_exp;
    math::Matrix3x3 cov = pvCov + svCov;
    int errorFlag = 0;
    math::Matrix3x3 covInv = cov.Inverse(errorFlag);
    if ( errorFlag != 0 )
    {
      if ( verbosity >= 0 )
      {
        printCovMatrix("cov", cov);
      }
      throw cmsException("comp_tipCompatibilityThreeProng", __LINE__) 
        << "Failed to invert matrix cov !!\n";
    }
    // CV: large negative values of tipCompatibility indicate tension between the tau momentum vector and the vector flightlength = sv - pv,
    //     while a tipCompatibility of zero indicates perfect agreement
    math::Vector3 d(residual.x(), residual.y(), residual.z());
    double tipCompatibility = -ROOT::Math::Dot(d, covInv*d);
    return tipCompatibility;
  }

  double
  comp_tipCompatibilityOneProng(const reco::Candidate::LorentzVector& tauP4,
                                const reco::Candidate::Point& pv, const math::Matrix3x3& pvCov,
                                const reco::Candidate::LorentzVector& leadTrackP4,
                                const reco::Candidate::Point& tipPCA, const math::Matrix3x3& svCov,
                                int verbosity = -1, bool cartesian = true)
  {
    if ( verbosity >= 3 )
    {
      std::cout << "<comp_tipCompatibilityOneProng>:\n";
    }

    std::pair<reco::Candidate::Point, reco::Candidate::Point> Qs = comp_PCA_line2line(pv, tauP4.Vect(), tipPCA, leadTrackP4.Vect(), verbosity);
    const reco::Candidate::Point& Q1 = Qs.first;
    const reco::Candidate::Point& Q2 = Qs.second;
    auto residual = Q2 - Q1; 
    math::Matrix3x3 cov = pvCov + svCov;
    int errorFlag = 0;
    math::Matrix3x3 covInv = cov.Inverse(errorFlag);
    if ( errorFlag != 0 )
    {
      if ( verbosity >= 0 )
      {
        printCovMatrix("cov", cov);
      }
      throw cmsException("comp_tipCompatibilityOneProng", __LINE__) 
        << "Failed to invert matrix cov !!\n";
    }
    // CV: large negative values of tipCompatibility indicate sizeable distances between the tau momentum vector and the trajectory of the track,
    //     while a tipCompatibility of zero indicates that tau momentum vector and trajectory intersect
    math::Vector3 d(residual.x(), residual.y(), residual.z());
    double tipCompatibility = -ROOT::Math::Dot(d, covInv*d);
    return tipCompatibility;
  }
}

double
StartPosTIPCompatibility::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 3 )
  {
    std::cout << "<StartPosTIPCompatibility::operator()>:\n"; 
  }

  // CV: compute compatibility of StartPosFinder solutions with transverse impact parameters,
  //     using the procedure described in the paper arXiv:hep-ph/9307269

  int tauPlus_decaymode = kineEvt.tauPlus_decayMode();
  int tauMinus_decaymode = kineEvt.tauMinus_decayMode();

  double tipCompatibility = 0.;
  if ( tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero )
  {
    tipCompatibility += comp_tipCompatibilityThreeProng(kineEvt.tauPlusP4(), kineEvt.pv(), kineEvt.pvCov(), kineEvt.svTauPlus(), kineEvt.svTauPlusCov(),
                                                        verbosity_, cartesian_);
  }
  else
  {
    const KinematicParticle* tauPlus_leadTrack = get_leadTrack(kineEvt.daughtersTauPlus());
    assert(tauPlus_leadTrack);
    tipCompatibility += comp_tipCompatibilityOneProng(kineEvt.tauPlusP4(),
                                                      kineEvt.pv(), kineEvt.pvCov(),
                                                      tauPlus_leadTrack->p4(),
                                                      kineEvt.tipPCATauPlus(), kineEvt.svTauPlusCov(),
                                                      verbosity_, cartesian_);
  }
  if ( tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero )
  {
    tipCompatibility += comp_tipCompatibilityThreeProng(kineEvt.tauMinusP4(), kineEvt.pv(), kineEvt.pvCov(), kineEvt.svTauMinus(), kineEvt.svTauMinusCov(),
                                                        verbosity_, cartesian_);
  }
  else
  {
    const KinematicParticle* tauMinus_leadTrack = get_leadTrack(kineEvt.daughtersTauMinus());
    assert(tauMinus_leadTrack);
    tipCompatibility += comp_tipCompatibilityOneProng(kineEvt.tauMinusP4(),
                                                      kineEvt.pv(), kineEvt.pvCov(),
                                                      tauMinus_leadTrack->p4(),
                                                      kineEvt.tipPCATauMinus(), kineEvt.svTauMinusCov(),
                                                      verbosity_, cartesian_);
  }
  if ( verbosity_ >= 3 )
  {
    std::cout << "tipCompatibility = " << tipCompatibility << "\n";
  }
  return tipCompatibility;
}
