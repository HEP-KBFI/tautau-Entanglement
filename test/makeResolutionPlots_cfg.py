import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    maxEvents_beforeCuts = cms.int32(-1),
    maxEvents_afterCuts = cms.int32(-1),
    outputEvery = cms.uint32(10000)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('')
)

process.makeResolutionPlots = cms.PSet(
    treeName = cms.string('ntupleProducer/piPlus_piMinus'),

    mode = cms.string('gen'),

    collider = cms.string(''),

    #minVisTauPt = cms.double(20.),
    #maxAbsVisTauEta = cms.double(2.3),
    minVisTauPt = cms.double(-1.),
    maxAbsVisTauEta = cms.double(1.e+3),
    #minTauTIP = cms.double(0.0030),
    minTauTIP = cms.double(-1.),
    maxNumChargedKaons = cms.int32(0),
    maxNumNeutralKaons = cms.int32(0),
    maxNumPhotons = cms.int32(-1),
    maxSumPhotonEn = cms.double(0.5),

    maxChi2 = cms.double(1.e+2),
    statusSelection = cms.vint32(0,1),

    branchName_evtWeight = cms.string('evtWeight'),

    isDEBUG = cms.bool(False)
)

inputFilePath = '/scratch/persistent/veelken/Entanglement/ntuples/2023Aug02/'
inputFileNames = None
processName = "ggH_htt_pythia8"
mode = 'gen'
collider = "LHC"
hAxis = "higgs"
decayMode = "piPlus_piMinus"
outputFileName = 'makeResolutionPlots_%s_%sMode_%sAxis.root' % (processName, mode, hAxis)

#minVisTauPt = 20.
#maxAbsVisTauEta = 2.3
minVisTauPt = -1.
maxAbsVisTauEta = 1.e+3
#minTauTIP = 0.0030
minTauTIP = -1.
maxNumChargedKaons = 0
maxNumNeutralKaons = 0
maxNumPhotons = -1
maxSumPhotonEn = 5.

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##mode = "$mode"
##collider = "$collider"
##hAxis = "$hAxis"
##decayMode = "$decayMode"
##outputFileName = "$outputFileName"

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

process.makeResolutionPlots.treeName = 'ntupleProducer/%s' % decayMode
process.makeResolutionPlots.mode = mode
process.makeResolutionPlots.collider = collider
process.makeResolutionPlots.minVisTauPt = minVisTauPt
process.makeResolutionPlots.maxAbsVisTauEta = maxAbsVisTauEta
process.makeResolutionPlots.minTauTIP = minTauTIP
process.makeResolutionPlots.maxNumChargedKaons = maxNumChargedKaons
process.makeResolutionPlots.maxNumNeutralKaons = maxNumNeutralKaons
process.makeResolutionPlots.maxNumPhotons = maxNumPhotons
process.makeResolutionPlots.maxSumPhotonEn = maxSumPhotonEn

