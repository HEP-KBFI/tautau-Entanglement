import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    maxEvents_beforeCuts = cms.int32(-1),
    ##maxEvents_afterCuts = cms.int32(10000),
    maxEvents_afterCuts = cms.int32(-1),
    outputEvery = cms.uint32(10000)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('')
)

process.analyzeEntanglementNtuple = cms.PSet(
    treeName = cms.string('ntupleProducer/piPlus_piMinus'),

    mode = cms.string('gen'),

    #minVisTauPt = cms.double(20.),
    #maxAbsVisTauEta = cms.double(2.3),
    minVisTauPt = cms.double(0.),
    maxAbsVisTauEta = cms.double(1.e+3),
    maxNumChargedKaons = cms.int32(0),
    maxNumNeutralKaons = cms.int32(0),
    maxNumPhotons = cms.int32(-1),
    maxSumPhotonEn = cms.double(0.5),

    branchName_evtWeight = cms.string('evtWeight'),

    isDEBUG = cms.bool(False)
)

inputFilePath = '/scratch/persistent/veelken/Entanglement/ntuples/2023Aug02/'
inputFileNames = None
processName = "ggH_htt_pythia8"
mode = 'gen'
hAxis = "higgs"
outputFileName = 'analyzeEntanglementNtuple_%s_%sMode_%sAxis.root' % (processName, mode, hAxis)

treeName = 'ntupleProducer/piPlus_piMinus'
#minVisTauPt = 20.
#maxAbsVisTauEta = 2.3
minVisTauPt = 0.
maxAbsVisTauEta = 1.e+3
maxNumChargedKaons = 0
maxNumNeutralKaons = 0
maxNumPhotons = -1
maxSumPhotonEn = 5.

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##mode = "$mode"
##hAxis = "$hAxis"
##outputFileName = "$outputFileName"

##treeName = "$treeName"
##minVisTauPt = $minVisTauPt
##maxAbsVisTauEta = $maxAbsVisTauEta
##maxNumChargedKaons = $maxNumChargedKaons
##maxNumNeutralKaons = $maxNumNeutralKaons
##maxNumPhotons = $maxNumPhotons
##maxSumPhotonEn = $maxSumPhotonEn

inputFile_regex = r"entanglementNtuple_%s_%sAxis_[0-9]+.root" % (processName, hAxis)

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
    print("Found %i input files." % len(inputFileNames))
else:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
#--------------------------------------------------------------------------------

process.fwliteInput.fileNames = cms.vstring(inputFileNames)

process.fwliteOutput.fileName = cms.string(outputFileName)

process.analyzeEntanglementNtuple.treeName = treeName
process.analyzeEntanglementNtuple.mode = mode
process.analyzeEntanglementNtuple.minVisTauPt = minVisTauPt
process.analyzeEntanglementNtuple.maxAbsVisTauEta = maxAbsVisTauEta
process.analyzeEntanglementNtuple.maxNumChargedKaons = maxNumChargedKaons
process.analyzeEntanglementNtuple.maxNumNeutralKaons = maxNumNeutralKaons
process.analyzeEntanglementNtuple.maxNumPhotons = maxNumPhotons
process.analyzeEntanglementNtuple.maxSumPhotonEn = maxSumPhotonEn
