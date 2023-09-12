#include "TauAnalysis/Entanglement/plugins/KinFitConstraintAnalyzer.h"

#include "DataFormats/Common/interface/Handle.h"                       // edm::Handle<>

#include "TauAnalysis/Entanglement/interface/cmsException.h"           // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"              // kLHC, kSuperKEKB
#include "TauAnalysis/Entanglement/interface/KinFitConstraint.h"       // KinFitConstraint
#include "TauAnalysis/Entanglement/interface/KinFitConstraintTester.h" // KinFitConstraintTester

#include <TString.h>                                                   // Form()

#include <iostream>                                                    // std::cout

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
    printKinematicEvent("kineEvt_gen(smeared)", kineEvt_gen, verbosity_, cartesian_);
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

      const reco::Candidate::LorentzVector& visTauPlusP4_gen = kineEvt_gen.visTauPlusP4();
      const reco::Candidate::LorentzVector& nuTauPlusP4_gen = kineEvt_gen.nuTauPlusP4();
      const reco::Candidate::LorentzVector& visTauMinusP4_gen = kineEvt_gen.visTauMinusP4();
      const reco::Candidate::LorentzVector& nuTauMinusP4_gen = kineEvt_gen.nuTauMinusP4();

      // CV: Check that chosen signs for tau+ and tau- reproduce neutrino Pz of start position; skip sign combination if not.
      //     Skipping these sign combinations saves computing time and avoids running into unphysical solutions.
      double nuTauPlusPx = nuTauPlusP4_gen.px();
      double nuTauPlusPy = nuTauPlusP4_gen.py();
      double nuTauPlus_dPzdPx, nuTauPlus_dPzdPy;
      double nuTauPlusPz = comp_nuPz(visTauPlusP4_gen, nuTauPlusPx, nuTauPlusPy, signTauPlus, nuTauPlus_dPzdPx, nuTauPlus_dPzdPy, verbosity_);
      double nuTauMinusPx = nuTauMinusP4_gen.px();
      double nuTauMinusPy = nuTauMinusP4_gen.py();
      double nuTauMinus_dPzdPx, nuTauMinus_dPzdPy;
      double nuTauMinusPz = comp_nuPz(visTauMinusP4_gen, nuTauMinusPx, nuTauMinusPy, signTauMinus, nuTauMinus_dPzdPx, nuTauMinus_dPzdPy, verbosity_);
      if ( std::fabs(nuTauPlusPz  - nuTauPlusP4_gen.pz())  > std::max(1., 0.10*nuTauPlusP4_gen.pz())  ||
           std::fabs(nuTauMinusPz - nuTauMinusP4_gen.pz()) > std::max(1., 0.10*nuTauMinusP4_gen.pz()) )
      {
        if ( verbosity_ >= 1 )
        {
          std::cout << "--> skipping this sign combination, because it does not reproduce the neutrino Pz of the start position !!\n";
        }
        continue;
      }

      if ( collider_ == kLHC ) 
      {
        const int C = numConstraints_LHC;
        KinFitConstraint<P,C> constraint(collider_, kineEvt_gen, signTauPlus, signTauMinus, verbosity_);
        KinFitConstraintTester<P,C> constraintTester(verbosity_);
        std::string outputFileName = Form("testConstraint_r%uls%uev%llu.png", evt.id().run(), evt.id().luminosityBlock(), evt.id().event());
        constraintTester(constraint, outputFileName);
      }
      else if ( collider_ == kSuperKEKB )
      {
        const int C = numConstraints_SuperKEKB;
        KinFitConstraint<P,C> constraint(collider_, kineEvt_gen, signTauPlus, signTauMinus, verbosity_);
        KinFitConstraintTester<P,C> constraintTester(verbosity_);
        std::string outputFileName = Form("testConstraint_r%uls%uev%llu.png", evt.id().run(), evt.id().luminosityBlock(), evt.id().event());
        constraintTester(constraint, outputFileName);
      }
      else assert(0);
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(KinFitConstraintAnalyzer);
