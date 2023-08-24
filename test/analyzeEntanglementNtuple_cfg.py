import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    maxEvents_beforeCuts = cms.int32(-1),
    maxEvents_afterCuts = cms.int32(2000),
    ##maxEvents_afterCuts = cms.int32(-1),
    outputEvery = cms.uint32(10000)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('')
)

process.analyzeEntanglementNtuple = cms.PSet(
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
    apply_evtWeight = cms.bool(True),

    par_gen = cms.vdouble(),

    spinAnalyzer = cms.string('by_summation'), # CV: either 'by_summation' or 'by_mlfit'
    numBootstrapSamples = cms.uint32(1000),

    mlfit_outputFileName = cms.string(""),
    mlfit_scan_likelihood = cms.bool(False),

    verbosity = cms.untracked.int32(1)
)

inputFilePath = '/scratch/persistent/veelken/Entanglement/ntuples/2023Aug02/'
inputFileNames = None
processName = "ggH_htt_pythia8"
mode = 'gen'
collider = "LHC"
hAxis = "higgs"
decayMode = "piPlus_piMinus"
apply_evtWeight = True
# CV: define Standard Model expectation for matrix C, given by Eq. (69) of arXiv:2208:11723, 
#     for comparison with measured matrix elements
par_gen = [ 0., 0., 0., 0., 0., 0., +1., 0., 0., 0., +1., 0., 0., 0., -1. ]
spinAnalyzer = "by_summation"
outputFileName = 'analyzeEntanglementNtuple_%s_%sMode_%sAxis.root' % (processName, mode, hAxis)

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
##apply_evtWeight = $apply_evtWeight
##par_gen = $par_gen
##spinAnalyzer = "$spinAnalyzer"
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

process.analyzeEntanglementNtuple.treeName = 'ntupleProducer/%s' % decayMode
process.analyzeEntanglementNtuple.mode = mode
process.analyzeEntanglementNtuple.collider = collider
process.analyzeEntanglementNtuple.minVisTauPt = minVisTauPt
process.analyzeEntanglementNtuple.maxAbsVisTauEta = maxAbsVisTauEta
process.analyzeEntanglementNtuple.minTauTIP = minTauTIP
process.analyzeEntanglementNtuple.maxNumChargedKaons = maxNumChargedKaons
process.analyzeEntanglementNtuple.maxNumNeutralKaons = maxNumNeutralKaons
process.analyzeEntanglementNtuple.maxNumPhotons = maxNumPhotons
process.analyzeEntanglementNtuple.maxSumPhotonEn = maxSumPhotonEn
process.analyzeEntanglementNtuple.apply_evtWeight = apply_evtWeight
process.analyzeEntanglementNtuple.par_gen = par_gen
process.analyzeEntanglementNtuple.spinAnalyzer = spinAnalyzer
process.analyzeEntanglementNtuple.mlfit_outputFileName = outputFileName
