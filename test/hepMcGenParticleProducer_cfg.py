import FWCore.ParameterSet.Config as cms

process = cms.Process("READHEPMC")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

process.source = cms.Source("EmptySource")

process.hepMcGenParticleProducer = cms.EDProducer(
    "HepMCGenParticleProducer",
    inputFile = cms.string("file:///home/karl/sandbox/tag_1_pythia8_events.hepmc")
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('hepmc.root'),
    outputCommands = cms.untracked.vstring(
        'keep *',
    )
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 999

process.p = cms.Path(process.hepMcGenParticleProducer)
process.e = cms.EndPath(process.out)
