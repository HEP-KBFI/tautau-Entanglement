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

process.dumpEntanglementNtuple = cms.PSet(
    treeName = cms.string('ntupleProducer/piPlus_piMinus'),

    mode = cms.string('gen'),

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

    branchName_evtWeight = cms.string('evtWeight'),
    apply_evtWeight = cms.bool(True),

    verbosity = cms.untracked.int32(1)
)

#inputFilePath = '/scratch/persistent/veelken/Entanglement/ntuples/2023Aug02/'
#inputFileNames = None
inputFilePath = None
inputFileNames = [ '/home/veelken/Entanglement/CMSSW_12_4_8/src/TauAnalysis/Entanglement/test/DEBUG/entanglementNtuple_ee2tt_unwgt_DEBUG_forLuca.root' ]
processName = "ggH_htt_pythia8"
mode = 'gen'
hAxis = "beam"
decayMode = "piPlus_piMinus"
apply_evtWeight = False

minVisTauPt = -1.
maxAbsVisTauEta = 1.e+3
minTauTIP = -1.
maxNumChargedKaons = 0
maxNumNeutralKaons = 0
maxNumPhotons = -1
maxSumPhotonEn = 5.

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##mode = "$mode"
##hAxis = "$hAxis"
##decayMode = "$decayMode"
##apply_evtWeight = $apply_evtWeight

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

process.fwliteOutput.fileName = cms.string("")

process.dumpEntanglementNtuple.treeName = 'ntupleProducer/%s' % decayMode
process.dumpEntanglementNtuple.mode = mode
process.dumpEntanglementNtuple.minVisTauPt = minVisTauPt
process.dumpEntanglementNtuple.maxAbsVisTauEta = maxAbsVisTauEta
process.dumpEntanglementNtuple.minTauTIP = minTauTIP
process.dumpEntanglementNtuple.maxNumChargedKaons = maxNumChargedKaons
process.dumpEntanglementNtuple.maxNumNeutralKaons = maxNumNeutralKaons
process.dumpEntanglementNtuple.maxNumPhotons = maxNumPhotons
process.dumpEntanglementNtuple.maxSumPhotonEn = maxSumPhotonEn
process.dumpEntanglementNtuple.apply_evtWeight = apply_evtWeight
