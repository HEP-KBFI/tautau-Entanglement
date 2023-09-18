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
    minVisTauZ = cms.double(-1.),
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

    # CV: Disabled correction for selection bias by default, as current implementation of TMVA 
    #     does not allow to run multiple analyzeEntanglementNtuple jobs in parallel.
    #     The issue is that the XML files containing the trained KNN configuration clash,
    #     as each analyzeEntanglementNtuple job tries to write to and read from the same XML file !!
    read_selBiasCorrection = cms.bool(False),
    write_selBiasCorrection = cms.bool(False),
    selBiasCorrection_outputFileName = cms.string("dataset/weights/TMVAClassification_KNN_weights.xml"),

    mlfit_outputFileName = cms.string(""),
    mlfit_scan_likelihood = cms.bool(False),

    verbosity = cms.untracked.int32(1)
)

inputFilePath = '/scratch/persistent/veelken/Entanglement/ntuples/SuperKEKB/2023Sep14_wSmearing/'
inputFileNames = None
processName = "dy_lo_pythia8"
# CV: define Standard Model expectation for matrix C, given by Eq. (69) of arXiv:2208:11723, 
#     for comparison with measured matrix elements
par_gen = [ 0., 0., 0., 0., 0., 0., -0.249686, 0.00508408, 0.00391348, 0.00653464, 0.138683, -0.276996, -0.00186421, 0.243658, 0.885224 ]
mode = 'gen'
collider = "SuperKEKB"
hAxis = "beam"
decayMode = "piPlus_piMinus"
apply_evtWeight = True
spinAnalyzer = "by_summation"
read_selBiasCorrection = False
write_selBiasCorrection = False
selBiasCorrection_outputFileName = 'dataset/weights/TMVAClassification_KNN_weights.xml'
outputFileName = 'analyzeEntanglementNtuple_%s_%sMode_%sAxis_%sDecayMode_%s.root' % (processName, mode, hAxis, decayMode, spinAnalyzer)

#minVisTauPt = 2.5
#maxAbsVisTauEta = 2.3
minVisTauPt = -1.
maxAbsVisTauEta = 1.e+3
minVisTauZ = -1.
#minTauTIP = 0.0030
minTauTIP = -1.
maxNumChargedKaons = 0
maxNumNeutralKaons = 0
maxNumPhotons = -1
maxSumPhotonEn = None
if collider == "LHC":
    maxSumPhotonEn = 5.
elif collider == "SuperKEKB":
    maxSumPhotonEn = 0.2
else:
    raise ValueError("Invalid Configuration parameter 'collider' = '%s' !!" % collider)

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##par_gen = $par_gen
##mode = "$mode"
##collider = "$collider"
##hAxis = "$hAxis"
##decayMode = "$decayMode"
##apply_evtWeight = $apply_evtWeight
##spinAnalyzer = "$spinAnalyzer"
##read_selBiasCorrection = $read_selBiasCorrection
##write_selBiasCorrection = $write_selBiasCorrection
##selBiasCorrection_outputFileName = "$selBiasCorrection_outputFileName"
##outputFileName = "$outputFileName"

##minVisTauPt = $minVisTauPt
##maxAbsVisTauEta = $maxAbsVisTauEta
##minVisTauZ = $minVisTauZ
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
process.analyzeEntanglementNtuple.read_selBiasCorrection = read_selBiasCorrection
process.analyzeEntanglementNtuple.write_selBiasCorrection = write_selBiasCorrection
process.analyzeEntanglementNtuple.selBiasCorrection_outputFileName = selBiasCorrection_outputFileName
