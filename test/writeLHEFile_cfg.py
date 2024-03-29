
import FWCore.ParameterSet.Config as cms

process = cms.Process("produceEntanglementNtuple")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2018_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    #input = cms.untracked.int32(10)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)
#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1:10' 
#)

#inputFilePath = '/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/100000/'
#inputFileNames = None
inputFilePath = None
inputFileNames = [ 'file:/local/karl/ee2tt_aod_unwgt/aodsim_1.root' ]

applyTauPairMassSelection  = True
applyTauDecayModeSelection = True
applyVisPtAndEtaSelection  = False

##inputFilePath = None
##inputFileNames = $inputFileNames

inputFile_regex = r"[a-zA-Z0-9-_]+.root"

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
    print("Found %i input files." % len(inputFileNames))
    #process.source.fileNames = cms.untracked.vstring(inputFileNames)
else:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun2_v2', '')

process.analysisSequence = cms.Sequence()

process.dumpGenParticles = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag('genParticles'),
    maxEventsToPrint = cms.untracked.int32(10) 
)
process.analysisSequence += process.dumpGenParticles

#--------------------------------------------------------------------------------
# apply event selection

if applyTauPairMassSelection:
    process.load("TauAnalysis.Entanglement.filterByTauPairMass_cff")
    process.genTaus.src = cms.InputTag('genParticles')
    process.genTauPair.cut = cms.string("mass > 0.")
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

#process.lheWriter = cms.EDAnalyzer("LHEWriter",
#    moduleLabel = cms.untracked.InputTag('externalLHEProducer')
#)
process.lheWriter = cms.EDAnalyzer("MyLHEWriter",
    srcLHERunInfo = cms.InputTag('externalLHEProducer'),
    srcLHEEvent = cms.InputTag('externalLHEProducer'),
    srcGenParticles = cms.InputTag('genParticles'),
    outputFileName = cms.string("aodsim_1_sel_piPlus_piMinus.lhe")
)
process.analysisSequence += process.lheWriter

process.p = cms.Path(process.analysisSequence)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
