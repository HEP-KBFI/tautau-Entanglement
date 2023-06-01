import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    maxEvents_beforeCuts = cms.int32(-1),
    maxEvents_afterCuts = cms.int32(10000),
    outputEvery = cms.uint32(10000)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('')
)

process.analyzeEntanglementNtuple = cms.PSet(
    treeName = cms.string('ntupleProducer/piPlus_piMinus'),

    minVisTauPt = cms.double(20.),
    maxAbsVisTauEta = cms.double(2.3),
    ##minVisTauPt = cms.double(0.),
    ##maxAbsVisTauEta = cms.double(1.e+3),

    branchName_evtWeight = cms.string('evtWeight'),

    isDEBUG = cms.bool(False)
)

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames

processName = "ggH_htt"
hAxis = "higgs"

inputFilePath = '/scratch/persistent/veelken/Entanglement/ntuples/2023Jun01/'
inputFile_regex = r"entanglementNtuple_%s_%sAxis_[0-9]+.root" % (processName, hAxis)
inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
print("Found %i input files." % len(inputFileNames))

outputFileName = 'analyzeEntanglementNtuple_%s_%s.root' % (processName, hAxis)

process.fwliteInput.fileNames = cms.vstring(inputFileNames)
process.fwliteOutput.fileName = cms.string(outputFileName)

