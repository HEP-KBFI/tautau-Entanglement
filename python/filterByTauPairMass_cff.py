import FWCore.ParameterSet.Config as cms

# CV: Veto events in which the tau leptons radiate high pT photons,
#     causing the mass of the tau pair to be significantly lower than the Higgs/Z boson mass.
#
#     The KinematicFit will likely fail for these events and may cause the neutrino momentum to become "not-a-number" (NaN),
#     triggering the following assert statement:
#       cmsRun: /home/veelken/Entanglement/CMSSW_12_4_8/src/TauAnalysis/Entanglement/src/SpinAnalyzerOneProng1Pi0.cc:103: reco::Candidate::Vector {anonymous}::getPolarimetricVec_OneProng1PiZero(const LorentzVector&, const std::vector<KinematicParticle>&, const LorentzVector&, const ROOT::Math::Boost&, const Vector&, const Vector&, const Vector&, const ROOT::Math::Boost&, int, bool): Assertion `nuP4.energy() >= 0. && N.energy() >= 0.' failed.
#
genTaus = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag('genParticles'),
    cut = cms.string("abs(pdgId) = 15 & status = 2"),
    filter = cms.bool(False)
)

genTauPair = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("genTaus@+ genTaus@-"),
    cut = cms.string("mass > 0."),
    filter = cms.bool(False)
)

selectedTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('genTauPair'),
    minNumber = cms.uint32(1)
)

filterByTauPairMass = cms.Sequence(genTaus + genTauPair + selectedTauPairFilter)

