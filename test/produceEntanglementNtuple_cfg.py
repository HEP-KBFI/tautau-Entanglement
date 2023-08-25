
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
    #input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)
#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1:10' 
#)

inputFilePath = '/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/100000/'
inputFileNames = None
processName = "qqH_htt_pythia8"
hAxis = "beam"
rndSeed = 1
outputFileName = "entanglementNtuple_%s_DEBUG.root" % processName
#collider = "LHC"
collider = "SuperKEKB"

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##collider = "$collider"
##hAxis = "$hAxis"
##rndSeed = $rndSeed
##outputFileName = "$outputFileName"

inputFile_regex = r"[a-zA-Z0-9-_]+.root"

srcGenParticles = None
tauPairMassCut = None
from TauAnalysis.Entanglement.resolutions_cfi import resolutions_LHC, resolutions_SuperKEKB
resolutions = None
applyHiggsMassConstraint = None
if collider == "LHC":
    process.source.fileNames = cms.untracked.vstring(
        'file:/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/100000/4FC3731A-E9C4-DD47-B222-83083ECF5684.root'
    )
    srcGenParticles = 'prunedGenParticles'
    tauPairMassCut = 'mass > 120. & mass < 130.'
    resolutions = resolutions_LHC
    applyHiggsMassConstraint = True
elif collider == "SuperKEKB":
    process.source.fileNames = cms.untracked.vstring(
        'file:/local/karl/ee2tt_aod/aodsim_1.root'
    )
    srcGenParticles = 'genParticles'
    tauPairMassCut = 'mass > 0.'
    resolutions = resolutions_SuperKEKB
    applyHiggsMassConstraint = False
else:
    raise ValueError("Invalid Configuration parameter 'collider' = '%s' !!" % collider)

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
    src = cms.InputTag(srcGenParticles),
    maxEventsToPrint = cms.untracked.int32(10) 
)
process.analysisSequence += process.dumpGenParticles

#--------------------------------------------------------------------------------
# CV: Veto events in which the tau leptons radiate high pT photons,
#     causing the mass of the tau pair to be significantly lower than the Higgs/Z boson mass.
#
#     The KinematicFit will likely fail for these events and may cause the neutrino momentum to become "not-a-number" (NaN),
#     triggering the following assert statement:
#       cmsRun: /home/veelken/Entanglement/CMSSW_12_4_8/src/TauAnalysis/Entanglement/src/SpinAnalyzerOneProng1Pi0.cc:103: reco::Candidate::Vector {anonymous}::getPolarimetricVec_OneProng1PiZero(const LorentzVector&, const std::vector<KinematicParticle>&, const LorentzVector&, const ROOT::Math::Boost&, const Vector&, const Vector&, const Vector&, const ROOT::Math::Boost&, int, bool): Assertion `nuP4.energy() >= 0. && N.energy() >= 0.' failed.
#
process.genTaus = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag(srcGenParticles),
    cut = cms.string("abs(pdgId) = 15 & status = 2"),
    filter = cms.bool(False)
)
process.analysisSequence += process.genTaus

process.genTauPair = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("genTaus@+ genTaus@-"),
    cut = cms.string(tauPairMassCut),
    filter = cms.bool(False)
)
process.analysisSequence += process.genTauPair

process.selectedTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('genTauPair'),
    minNumber = cms.uint32(1)
)
process.analysisSequence += process.selectedTauPairFilter
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: Require both generator-level tau leptons to decay into "oneProng0Pi0", "oneProng1Pi0", or "threeProng0Pi0".
#     The tau decay mode strings are defined in:
#       https://cmssdt.cern.ch/lxr/source/PhysicsTools/JetMCUtils/src/JetMCTag.cc
# 
#     Note that the KinematicFit will likely fail in case the event contains >= leptonic tau decay 
#     and the spin analyzer vectors will be zero,
#     so it is recommended to ALWAYS apply a tau decay mode selection 
#     and require that BOTH taus decay hadronically and to the selected decay modes
#
#     Similarly, it is recommended to ALWAYS apply the pT > 20 GeV and |eta| < 2.3 cuts,
#     in order to remove pathological events in which one of the tau leptons have very high longitudinal momentum
#     or the visible decay products have very low transverse momentum
#
process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.tauGenJets.GenParticles = cms.InputTag('genTaus')
process.analysisSequence += process.tauGenJets

process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.tauGenJetsSelectorAllHadrons.select = cms.vstring(
    'oneProng0Pi0', 
    'oneProng1Pi0', 
##    #'oneProng2Pi0', 
##    #'oneProngOther',
    'threeProng0Pi0', 
##    #'threeProng1Pi0', 
##    #'threeProngOther', 
##    #'rare'
)
process.analysisSequence += process.tauGenJetsSelectorAllHadrons

## CV: disabled pT and eta cuts on visible tau decay products 
##     to test analyzeEntanglementNtuples.cc code (17/08/2023)
##process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
##    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
##    cut = cms.string('pt > 20. & abs(eta) < 2.3'),
##    filter = cms.bool(False)
##)
##process.analysisSequence += process.selectedGenHadTaus

process.selectedGenHadTauFilter = cms.EDFilter("CandViewCountFilter",
    ##src = cms.InputTag('selectedGenHadTaus'),
    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    minNumber = cms.uint32(2)
)
process.analysisSequence += process.selectedGenHadTauFilter
#--------------------------------------------------------------------------------

process.genWeight = cms.EDProducer("GenWeightProducer",
    src = cms.InputTag('generator')
)
process.analysisSequence += process.genWeight
 
from TauAnalysis.Entanglement.smearing_cfi import smearing
process.ntupleProducer = cms.EDAnalyzer("EntanglementNtupleProducer",
    src = cms.InputTag(srcGenParticles),
    collider = cms.string(collider),
    hAxis = cms.string(hAxis),
    resolutions = resolutions,
    smearing = smearing.clone(
        rndSeed = cms.uint64(rndSeed)
    ),
    #applySmearing = cms.bool(False),
    applySmearing = cms.bool(True),
    startPosFinder = cms.PSet(
        algos = cms.vint32(1),
        applyHiggsMassConstraint = cms.bool(applyHiggsMassConstraint),
        applyRecoilEnergy_and_PzConstraint = cms.bool(True)
    ),
    kinematicFit = cms.PSet(
        applyLifetimeConstraint = cms.bool(False)
    ),
    srcEvtWeights = cms.VInputTag('genWeight'),
    verbosity = cms.untracked.int32(-1),
    #verbosity = cms.untracked.int32(3),
    cartesian = cms.untracked.bool(True)
)
process.analysisSequence += process.ntupleProducer

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(outputFileName),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.analysisSequence)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
