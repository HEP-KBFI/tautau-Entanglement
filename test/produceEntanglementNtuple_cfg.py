
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
    ##input = cms.untracked.int32(10000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    ##eventsToProcess = cms.untracked.VEventRange(
    ##    '1:97:96087'
    ##) 
)

inputFilePath = '/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/100000/'
inputFileNames = None
processName = "qqH_htt"
hAxis = "higgs"
outputFileName = "entanglementNtuple_%s_DEBUG.root" % processName

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##hAxis = "$hAxis"
##outputFileName = "$outputFileName"

inputFile_regex = r"[a-zA-Z0-9-_]+.root"

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
    print("Found %i input files." % len(inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
else:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v24', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# CV: generator-level tau decay mode selection
##process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
##process.tauGenJets.GenParticles = cms.InputTag('prunedGenParticles')
##process.analysisSequence += process.tauGenJets
##
##process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
##process.tauGenJetsSelectorAllHadrons.select = cms.vstring(
##    'oneProng0Pi0', 
##    #'oneProng1Pi0', 
##    #'oneProng2Pi0', 
##    #'oneProngOther',
##    #'threeProng0Pi0', 
##    #'threeProng1Pi0', 
##    #'threeProngOther', 
##    #'rare'
##)
##process.analysisSequence += process.tauGenJetsSelectorAllHadrons
##
##process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
##    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
##    cut = cms.string('pt > 20. & abs(eta) < 2.3'),
##    filter = cms.bool(False)
##)
##process.analysisSequence += process.selectedGenHadTaus
##
##process.selectedGenHadTauFilter = cms.EDFilter("CandViewCountFilter",
##    src = cms.InputTag('selectedGenHadTaus'),
##    minNumber = cms.uint32(2)
##)
##process.analysisSequence += process.selectedGenHadTauFilter
#--------------------------------------------------------------------------------

process.dumpGenParticles = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag('prunedGenParticles'),
    maxEventsToPrint = cms.untracked.int32(10) 
)
process.analysisSequence += process.dumpGenParticles

process.genWeight = cms.EDProducer("GenWeightProducer",
    src = cms.InputTag('generator')
)
process.analysisSequence += process.genWeight

process.ntupleProducer = cms.EDAnalyzer("EntanglementNtupleProducer",
    src = cms.InputTag('prunedGenParticles'),
    hAxis = cms.string(hAxis),
    srcEvtWeights = cms.VInputTag('genWeight'),
    verbosity = cms.untracked.int32(-1)
    ##verbosity = cms.untracked.int32(1)
)
process.analysisSequence += process.ntupleProducer

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(outputFileName),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.analysisSequence)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)
