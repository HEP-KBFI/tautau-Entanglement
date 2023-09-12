import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeKinFitConstraint")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2018_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(-1)
    input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)
#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1:10' 
#)

inputFileNames = [ 'file:/local/karl/ee2tt_aod_unwgt/aodsim_1.root' ]

collider = "SuperKEKB"

applyTauPairMassSelection  = True
applyTauDecayModeSelection = True
applyVisPtAndEtaSelection  = False

srcGenParticles = None
tauPairMassCut = None
from TauAnalysis.Entanglement.resolutions_cfi import resolutions_LHC, resolutions_SuperKEKB
resolutions = None
if collider == "LHC":
    srcGenParticles = 'prunedGenParticles'
    tauPairMassCut = "mass > 120. & mass < 130."
    resolutions = resolutions_LHC
elif collider == "SuperKEKB":
    srcGenParticles = 'genParticles'
    tauPairMassCut = "mass > 0."
    resolutions = resolutions_SuperKEKB
else:
    raise ValueError("Invalid Configuration parameter 'collider' = '%s' !!" % collider)

process.source.fileNames = cms.untracked.vstring(inputFileNames)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun2_v2', '')

process.analysisSequence = cms.Sequence()

process.dumpGenParticles = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag(srcGenParticles),
    maxEventsToPrint = cms.untracked.int32(10) 
)
process.analysisSequence += process.dumpGenParticles

#--------------------------------------------------------------------------------
# apply event selection

if applyTauPairMassSelection:
    process.load("TauAnalysis.Entanglement.filterByTauPairMass_cff")
    process.genTaus.src = cms.InputTag(srcGenParticles)
    process.genTauPair.cut = cms.string(tauPairMassCut)
    process.analysisSequence += process.filterByTauPairMass

if applyTauDecayModeSelection:
    process.load("TauAnalysis.Entanglement.filterByTauDecayMode_cff")
    process.tauGenJetsSelectorAllHadrons.select = cms.vstring("oneProng0Pi0")
    process.analysisSequence += process.filterByTauDecayMode

if applyVisPtAndEtaSelection:
    if not applyTauDecayModeSelection:
        process.load("TauAnalysis.Entanglement.filterByTauDecayMode_cff")
        process.analysisSequence += process.tauGenJets
        process.analysisSequence += process.tauGenJetsSelectorAllHadrons
    process.load("TauAnalysis.Entanglement.filterByVisPtAndEta_cff")
    process.analysisSequence += process.filterByVisPtAndEta
#--------------------------------------------------------------------------------

from TauAnalysis.Entanglement.smearing_cfi import smearing
process.analyzeKinFitConstraint = cms.EDAnalyzer("KinFitConstraintAnalyzer",
    src = cms.InputTag(srcGenParticles),
    collider = cms.string(collider),
    resolutions = resolutions,
    smearing = smearing.clone(
        rndSeed = cms.uint64(1)
    ),
    #applySmearing = cms.bool(False),
    applySmearing = cms.bool(True),
    verbosity = cms.untracked.int32(-1),
    #verbosity = cms.untracked.int32(3),
    cartesian = cms.untracked.bool(True)
)
process.analysisSequence += process.analyzeKinFitConstraint

process.p = cms.Path(process.analysisSequence)
