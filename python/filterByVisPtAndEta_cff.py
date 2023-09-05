import FWCore.ParameterSet.Config as cms

selectedGenHadTausByPtAndEta = cms.EDFilter("GenJetSelector",
    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    cut = cms.string('pt > 20. & abs(eta) < 2.3'),
    filter = cms.bool(False)
)

selectedGenHadTauFilterByPtAndEta = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedGenHadTausByPtAndEta'),
    minNumber = cms.uint32(2)
)

filterByVisPtAndEta = cms.Sequence(selectedGenHadTausByPtAndEta + selectedGenHadTauFilterByPtAndEta)
