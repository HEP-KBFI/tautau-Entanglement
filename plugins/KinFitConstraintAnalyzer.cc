#include "TauAnalysis/Entanglement/plugins/KinFitConstraintAnalyzer.h"

#include "DataFormats/Common/interface/Handle.h"                                 // edm::Handle<>

#include "TauAnalysis/Entanglement/interface/cmsException.h"                     // cmsException
#include "TauAnalysis/Entanglement/interface/comp_nuPz.h"                        // comp_nuPz()
#include "TauAnalysis/Entanglement/interface/constants.h"                        // kLHC, kSuperKEKB
#include "TauAnalysis/Entanglement/interface/KinFitConstraint.h"                 // KinFitConstraint
#include "TauAnalysis/Entanglement/interface/KinFitConstraintTester.h"           // KinFitConstraintTester
#include "TauAnalysis/Entanglement/interface/KinFitParameters_and_Constraints.h" // kinFit::numParameters

#include <TString.h>                                                             // Form()
#include <TVectorD.h>                                                            // TVectorD

#include <iostream>                                                              // std::cout

using namespace math;

const int numConstraints_LHC = 9;
const int numConstraints_SuperKEKB = 8;

KinFitConstraintAnalyzer::KinFitConstraintAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , genKineEvtBuilder_(nullptr)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenParticleCollection>(src_);

  genKineEvtBuilder_ = new GenKinematicEventBuilder(cfg);

  applySmearing_ = cfg.getParameter<bool>("applySmearing"); 
  if ( applySmearing_ )
  {
    std::cerr << "WARNING: Smearing enabled in Configuration file - this will neutrino Px, Py not to be initialized by GenKinematicEventBuilder !!\n";
  }

  std::string collider = cfg.getParameter<std::string>("collider");
  if      ( collider == "LHC"       ) collider_ = kLHC;
  else if ( collider == "SuperKEKB" ) collider_ = kSuperKEKB;
  else throw cmsException("EntanglementNtupleProducer", __LINE__)
    << "Invalid Configuration parameter 'collider' = " << collider << " !!\n";
}

KinFitConstraintAnalyzer::~KinFitConstraintAnalyzer()
{
  delete genKineEvtBuilder_;
}

void KinFitConstraintAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<KinFitConstraintAnalyzer::analyze>:\n";
    std::cout << " moduleLabel = " << moduleLabel_ << "\n";
  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(token_, genParticles);

  KinematicEvent kineEvt_gen = (*genKineEvtBuilder_)(*genParticles);
  if ( verbosity_ >= 1 )
  {
    std::string label = ( applySmearing_ ) ? "kineEvt_gen(smeared)" : "kineEvt_gen";
    printKinematicEvent(label, kineEvt_gen, verbosity_, cartesian_);
  }

  for ( int idxSignTauPlus = 0; idxSignTauPlus <= 1; ++idxSignTauPlus )
  {
    double signTauPlus = 0.;
    if      ( idxSignTauPlus == 0 ) signTauPlus = -1.;
    else if ( idxSignTauPlus == 1 ) signTauPlus = +1.;
    else assert(0);
    for ( int idxSignTauMinus = 0; idxSignTauMinus <= 1; ++idxSignTauMinus )
    {
      double signTauMinus = 0.;
      if      ( idxSignTauMinus == 0 ) signTauMinus = -1.;
      else if ( idxSignTauMinus == 1 ) signTauMinus = +1.;
      else assert(0);

      if ( verbosity_ >= 1 )
      {
        std::cout << "sign(tau+) = " << signTauPlus << ", sign(tau-) = " << signTauMinus << "\n";
      }

      const reco::Candidate::Point& pv_gen = kineEvt_gen.pv();
      const reco::Candidate::LorentzVector& visTauPlusP4_gen = kineEvt_gen.visTauPlusP4();
      const reco::Candidate::LorentzVector& nuTauPlusP4_gen = kineEvt_gen.nuTauPlusP4();
      const reco::Candidate::Point& svTauPlus_gen = kineEvt_gen.svTauPlus();
      const reco::Candidate::LorentzVector& visTauMinusP4_gen = kineEvt_gen.visTauMinusP4();
      const reco::Candidate::LorentzVector& nuTauMinusP4_gen = kineEvt_gen.nuTauMinusP4();
      const reco::Candidate::Point& svTauMinus_gen = kineEvt_gen.svTauMinus();
      const reco::Candidate::LorentzVector& recoilP4_gen = kineEvt_gen.recoilP4();

      // CV: Check that chosen signs for tau+ and tau- reproduce neutrino Pz of start position; skip sign combination if not.
      //     Skipping these sign combinations saves computing time and avoids running into unphysical solutions.
      double nuTauPlusPx = nuTauPlusP4_gen.px();
      double nuTauPlusPy = nuTauPlusP4_gen.py();
      double nuTauPlus_dPzdPx, nuTauPlus_dPzdPy;
      bool nuTauPlus_errorFlag = false; 
      double nuTauPlusPz = comp_nuPz(visTauPlusP4_gen, nuTauPlusPx, nuTauPlusPy, signTauPlus, nuTauPlus_dPzdPx, nuTauPlus_dPzdPy, nuTauPlus_errorFlag, verbosity_);
      double nuTauMinusPx = nuTauMinusP4_gen.px();
      double nuTauMinusPy = nuTauMinusP4_gen.py();
      double nuTauMinus_dPzdPx, nuTauMinus_dPzdPy;
      bool nuTauMinus_errorFlag = false;
      double nuTauMinusPz = comp_nuPz(visTauMinusP4_gen, nuTauMinusPx, nuTauMinusPy, signTauMinus, nuTauMinus_dPzdPx, nuTauMinus_dPzdPy, nuTauMinus_errorFlag, verbosity_);
      if ( nuTauPlus_errorFlag || nuTauMinus_errorFlag )
      {
	if ( verbosity_ >= 1 )
        {
          std::cout << "--> skipping this sign combination, because computation of neutrino Pz returned an error !!\n";
        }
        continue;
      }  
      if ( std::fabs(nuTauPlusPz  - nuTauPlusP4_gen.pz())  > std::max(1., 0.10*nuTauPlusP4_gen.pz())  ||
           std::fabs(nuTauMinusPz - nuTauMinusP4_gen.pz()) > std::max(1., 0.10*nuTauMinusP4_gen.pz()) )
      {
        if ( verbosity_ >= 1 )
        {
          std::cout << "--> skipping this sign combination, because it does not reproduce the neutrino Pz of the start position !!\n";
        }
        continue;
      }

      TVectorD alpha0(kinFit::numParameters);
      alpha0( 0) = pv_gen.x();
      alpha0( 1) = pv_gen.y();
      alpha0( 2) = pv_gen.z();
      alpha0( 3) = nuTauPlusP4_gen.px();
      alpha0( 4) = nuTauPlusP4_gen.py();
      alpha0( 5) = svTauPlus_gen.x();
      alpha0( 6) = svTauPlus_gen.y();
      alpha0( 7) = svTauPlus_gen.z();
      alpha0( 8) = nuTauMinusP4_gen.px();
      alpha0( 9) = nuTauMinusP4_gen.py();
      alpha0(10) = svTauMinus_gen.x();
      alpha0(11) = svTauMinus_gen.y();
      alpha0(12) = svTauMinus_gen.z();
      alpha0(13) = recoilP4_gen.px();
      alpha0(14) = recoilP4_gen.py();
      alpha0(15) = recoilP4_gen.pz();
      alpha0(16) = recoilP4_gen.energy();
      if ( verbosity_ >= 2 )
      {
        std::cout << "alpha0:\n";
        alpha0.Print();
      }

      KinFitConstraint constraint(collider_, kineEvt_gen, signTauPlus, signTauMinus, verbosity_);
      KinFitConstraintTester constraintTester(alpha0, verbosity_);
      std::string outputFileName = Form("testConstraint_r%uls%uev%llu.png", evt.id().run(), evt.id().luminosityBlock(), evt.id().event());
      constraintTester(constraint, outputFileName);
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(KinFitConstraintAnalyzer);
