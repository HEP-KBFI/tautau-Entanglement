import FWCore.ParameterSet.Config as cms

# CV: Require both generator-level tau leptons to decay into "oneProng0Pi0", "oneProng1Pi0", or "threeProng0Pi0".
#     The tau decay mode strings are defined in:
#       https://cmssdt.cern.ch/lxr/source/PhysicsTools/JetMCUtils/src/JetMCTag.cc
# 
#     Note that the KinematicFit will likely fail in case the event contains >= leptonic tau decay 
#     and the spin analyzer vectors will be zero,
#     so it is recommended to ALWAYS apply a tau decay mode selection 
#     and require that BOTH taus decay hadronically and to the selected decay modes
#
#     Similarly, it is recommended to ALWAYS apply the pT > 20 GeV and |eta| < 2.3 cuts,
#     in order to remove pathological events in which one of the tau leptons have very high longitudinal momentum
#     or the visible decay products have very low transverse momentum
#
from PhysicsTools.JetMCAlgos.TauGenJets_cfi import tauGenJets
#tauGenJets.GenParticles = cms.InputTag('genTaus')
tauGenJets.verbose = cms.untracked.bool(True)

from PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi import tauGenJetsSelectorAllHadrons
tauGenJetsSelectorAllHadrons.select = cms.vstring(
    "oneProng0Pi0", 
    "oneProng1Pi0", 
    "oneProng2Pi0", 
    "oneProngOther",
    "threeProng0Pi0", 
    "threeProng1Pi0", 
    "threeProngOther", 
    "rare"
)

# CV: Reject tau decays with charged and neutral Kaons 
#     or with high pT photon radiation
selectedGenHadTaus = cms.EDProducer("GenTauDecaySelector",
    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    maxNumChargedKaons = cms.int32(0),
    maxNumNeutralKaons = cms.int32(0),
    maxNumPhotons = cms.int32(-1),
    #maxSumPhotonEn = cms.double(0.5),
    maxSumPhotonEn = cms.double(-1.),
    verbosity = cms.untracked.int32(2)
)

selectedGenHadTauFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedGenHadTaus'),
    minNumber = cms.uint32(2)
)

filterByTauDecayMode = cms.Sequence(tauGenJets + tauGenJetsSelectorAllHadrons + selectedGenHadTaus + selectedGenHadTauFilter)
