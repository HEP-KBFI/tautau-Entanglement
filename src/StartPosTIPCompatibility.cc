
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
/*
  reco::Candidate::Point 
  getPoint_rf(const reco::Candidate::Point& p,
              const ROOT::Math::Boost& boost)
  {
    reco::Candidate::LorentzVector p4(p.x(), p.y(), p4.z(), 0.);
    reco::Candidate::LorentzVector p4_rf = getP4_rf(p4, boost);
    reco::Candidate::Vector p_rf = p4_rf.Vect();
    return reco::Candidate::Point(p_rf.x(), p_rf.y(), p_rf.z());
  }

  double
  comp_tipCompatibilityOneProng(const KinematicEvent& kineEvt,
                                int verbosity = -1, bool cartesian = true)
  {
    // CV: compute { x, y, z } coordinate system as described in the text before Eq. (4) and in Fig. 1 of the paper arXiv:hep-ph/9307269
    if ( verbosity >= 3 )
    {
      std::cout << "<comp_tipCompatibilityOneProng>:\n";
    }

    const reco::Candidate::LorentzVector& tauPlusP4 = kineEvt.tauPlusP4();
    const reco::Candidate::LorentzVector& tauMinusP4 = kineEvt.tauMinusP4();
    if ( verbosity >= 3 )
    {
      printLorentzVector("tauPlusP4", tauPlusP4, cartesian);
      printLorentzVector("tauMinusP4", tauMinusP4, cartesian);
    }
    reco::Candidate::LorentzVector higgsP4 = tauPlusP4 + tauMinusP4;
    ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(higgsP4.BoostToCM());
    reco::Candidate::LorentzVector tauPlusP4_ttrf = getP4_rf(tauMinusP4, boost_ttrf);
    reco::Candidate::LorentzVector tauMinusP4_ttrf = getP4_rf(tauMinusP4, boost_ttrf);
    if ( verbosity >= 3 )
    {
      printLorentzVector("tauPlusP4_ttrf", tauPlusP4_ttrf, cartesian);
      printLorentzVector("tauMinusP4_ttrf", tauMinusP4_ttrf, cartesian);
    }

    reco::Candidate::Vector e_tauMinus = tauMinusP4_ttrf.Vect().unit();
    const reco::Candidate::LorentzVector& visTauMinusP4 = kineEvt.visTauMinusP4();
    reco::Candidate::LorentzVector visTauMinusP4_ttrf = getP4_rf(visTauMinusP4, boost_ttrf);
    if ( verbosity >= 4 )
    {
      printLorentzVector("visTauMinusP4_ttrf", visTauMinusP4_ttrf, cartesian);
    }
    reco::Candidate::Vector e_visTauMinus = visTauMinusP4_ttrf.Vect().unit();
    reco::Candidate::Vector z = e_tauMinus.unit();
    reco::Candidate::Vector x = get_r(z, e_visTauMinus.unit(), verbosity, cartesian);
    reco::Candidate::Vector y = get_n(x, z, verbosity, cartesian);
    if ( verbosity >= 4 )
    {
      printVector("xXy", x.Cross(y), cartesian);
      printVector("z", z, cartesian);
    }

    double nMinus_x = e_visTauMinus.Dot(x);
    double nMinus_y = e_visTauMinus.Dot(y);
    double nMinus_z = e_visTauMinus.Dot(z);
    reco::Candidate::Vector nMinus(nMinus_x, nMinus_y, nMinus_z);
    if ( verbosity >= 3 )
    {
      printVector("nMinus", nMinus, true);
      printVector("nMinus", nMinus, false);
    }

    const reco::Candidate::LorentzVector& visTauPlusP4 = kineEvt.visTauPlusP4();
    reco::Candidate::LorentzVector visTauPlusP4_ttrf = getP4_rf(visTauPlusP4, boost_ttrf);
    if ( verbosity >= 4 )
    {
      printLorentzVector("visTauPlusP4_ttrf", visTauPlusP4_ttrf, cartesian);
    }
    reco::Candidate::Vector e_visTauPlus = visTauPlusP4_ttrf.Vect().unit();
    double nPlus_x = e_visTauPlus.Dot(x);
    double nPlus_y = e_visTauPlus.Dot(y);
    double nPlus_z = e_visTauPlus.Dot(z);
    reco::Candidate::Vector nPlus(nPlus_x, nPlus_y, nPlus_z);
    if ( verbosity >= 3 )
    {
      printVector("nPlus", nPlus, true);
      printVector("nPlus", nPlus, false);
    }

    reco::Candidate::Vector e_perp = nPlus.Cross(nMinus);

    reco::Candidate::Vector d(0., 0., -1.);

    double term1 = 1. - pow(nMinus.Dot(nPlus), 2);
    double term2 = d.Dot(nPlus)*nPlus.Dot(nMinus) - d.Dot(nMinus);
    double term3 = d.Dot(nMinus)*nPlus.Dot(nMinus) - d.Dot(nPlus);
    reco::Candidate::Vector dmin_kine = d + (1./term1)*(term2*nMinus + term3*nPlus);
    double perp_kine = e_perp.Dot(dmin_kine);
    if ( verbosity >= 3 )
    {
      printDistance("dmin_kine", dmin_kine);
      std::cout << "perp_kine = " << perp_kine << "\n";
    }

    const reco::Candidate::Point& tipPCATauPlus = kineEvt.tipPCATauPlus();
    reco::Candidate::Point tipPCATauPlus_ttrf = getPoint_rf(tipPCATauPlus, boost_ttrf);
    const KinematicParticle* tauPlus_leadTrack = get_leadTrack(kineEvt.daughtersTauPlus());
    assert(tauPlus_leadTrack);
    const reco::Candidate::LorentzVector& tauPlus_leadTrackP4 = tauPlus_leadTrack->p4();
    reco::Candidate::LorentzVector tauPlus_leadTrackP4_ttrf = getP4_rf(tauPlus_leadTrackP4, boost_ttrf);
    
    const reco::Candidate::Point& tipPCATauMinus = kineEvt.tipPCATauMinus();
    reco::Candidate::Point tipPCATauMinus_ttrf = getPoint_rf(tipPCATauMinus, boost_ttrf);
    const KinematicParticle* tauMinus_leadTrack = get_leadTrack(kineEvt.daughtersTauMinus());
    assert(tauMinus_leadTrack);
    const reco::Candidate::LorentzVector& tauMinus_leadTrackP4 = tauMinus_leadTrack->p4();
    reco::Candidate::LorentzVector tauMinus_leadTrackP4_ttrf = getP4_rf(tauMinus_leadTrackP4, boost_ttrf);

    std::pair<reco::Candidate::Point, reco::Candidate::Point> Qs = comp_PCA_line2line(
      tipPCATauPlus_ttrf, tauPlus_leadTrackP4_ttrf.Vect(),
      tipPCATauMinus_ttrf, tauMinus_leadTrackP4_ttrf.Vect(), 
      verbosity);
    const reco::Candidate::Point& Q1 = Qs.first;
    const reco::Candidate::Point& Q2 = Qs.second;
    reco::Candidate::Vector dmin_tip = Q1 - Q2;
    double perp_tip = e_perp.Dot(dmin_tip);
    if ( verbosity >= 3 )
    {
      printDistance("dmin_tip", dmin_tip);
      std::cout << "perp_tip = " << perp_tip << "\n";
      std::cout << "dmin_kine*dmin_tip = " << dmin_kine.Dot(dmin_tip) << "\n";
    }

    // CV: tipCompatibility will be positive if the sign of sin(phi) is the same in Eqs. (4) and (8),
    //     and negative otherwise
    double tipCompatibility = perp_kine*perp_tip;
    return tipCompatibility;
  }
 */
  double
  get_tipPerp(const reco::Candidate::LorentzVector& tauP4, 
              const reco::Candidate::LorentzVector& leadTrackP4, 
              const reco::Candidate::Point& pv, const reco::Candidate::Point& tipPCA,
              int verbosity = -1, bool cartesian = true)
  {
    if ( verbosity >= 3 )
    {
      std::cout << "<get_tipPerp>:\n";
      printLorentzVector("tauP4", tauP4, cartesian);
      //printLorentzVector("leadTrackP4", leadTrackP4, cartesian);
      //printPoint("pv", pv, true);
      //printPoint("tipPCA", tipPCA, true);
    }
    reco::Candidate::Vector e_tau  = tauP4.Vect().unit();
    reco::Candidate::Vector e_vis  = leadTrackP4.Vect().unit();
    reco::Candidate::Vector e_perp = e_tau.Cross(e_vis).unit();
    reco::Candidate::Vector flightlength = tipPCA - pv;
    double tipPerp = std::fabs(e_perp.Dot(flightlength));
    if ( verbosity >= 3 )
    {
      std::cout << "tipPerp = " << tipPerp << "\n";
    }
    return tipPerp;
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
  if ( tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero || tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero )
  {
    if ( tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero )
    {
      tipCompatibility += comp_tipCompatibilityThreeProng(kineEvt.tauPlusP4(), kineEvt.pv(), kineEvt.pvCov(), kineEvt.svTauPlus(), kineEvt.svTauPlusCov(),
                                                          verbosity_, cartesian_);
    }
    if ( tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero )
    {
      tipCompatibility = comp_tipCompatibilityThreeProng(kineEvt.tauMinusP4(), kineEvt.pv(), kineEvt.pvCov(), kineEvt.svTauMinus(), kineEvt.svTauMinusCov(),
                                                         verbosity_, cartesian_);
    }
  }
  else
  {
    // CV: The procedure described in the paper arXiv:hep-ph/9307269 does not seem to work (yet);
    //     at least my implementation of the procedure does not seem to work :|
    //     For the time being, take solution with minimal components of transverse impact parameters 
    //     perpendicular to plane spanned by tau and leadTrack momentum vectors;
    //     the negative sign "penalizes" solutions for which the perpendicular components are large
    //tipCompatibility = comp_tipCompatibilityOneProng(kineEvt, verbosity_, cartesian_);
    const KinematicParticle* tauPlus_leadTrack = get_leadTrack(kineEvt.daughtersTauPlus());
    assert(tauPlus_leadTrack);
    tipCompatibility -= get_tipPerp(kineEvt.tauPlusP4(), tauPlus_leadTrack->p4(), kineEvt.pv(), kineEvt.tipPCATauPlus(),
                                    verbosity_, cartesian_);
    const KinematicParticle* tauMinus_leadTrack = get_leadTrack(kineEvt.daughtersTauMinus());
    assert(tauMinus_leadTrack);
    tipCompatibility -= get_tipPerp(kineEvt.tauMinusP4(), tauMinus_leadTrack->p4(), kineEvt.pv(), kineEvt.tipPCATauMinus(),
                                    verbosity_, cartesian_);
  }
  if ( verbosity_ >= 3 )
  {
    std::cout << "tipCompatibility = " << tipCompatibility << "\n";
  }
  return tipCompatibility;
}
