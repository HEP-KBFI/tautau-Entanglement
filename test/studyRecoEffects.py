#!/usr/bin/env python3

import getpass
import os

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames, build_Makefile
from TauAnalysis.Entanglement.samples import samples_LHC, samples_SuperKEKB

modes = [ "gen", "gen_smeared", "startPos", "kinFit" ]
hAxes = [ "beam", "higgs" ]
#hAxes = [ "higgs" ]
#collider = "LHC"
collider = "SuperKEKB"
#decayModes = [ "piPlus_piMinus", "piPlus_rhoMinus", "rhoPlus_piMinus", "rhoPlus_rhoMinus" ]
decayModes = [ "piPlus_piMinus" ]
spinAnalyzers = [ "by_summation", "by_mlfit", "by_differentialXsec1d", "by_differentialXsec2d", "by_asymmetry" ]
#spinAnalyzers = [ "by_summation" ]

version = "2023Sep14_wSmearing"

samples = None
if collider == "LHC":
    samples = samples_LHC
elif collider == "SuperKEKB":
    samples = samples_SuperKEKB
else:
    raise ValueError("Invalid Configuration parameter 'collider' = '%s' !!" % collider)

recoEffects = {
  'noVisTauPtCut' : {
    'minVisTauPt' : -1. 
  },
  'minVisTauPtGt0p5' : {
    'minVisTauPt' : 0.5 
  },
  'minVisTauPtGt1p0' : {
    'minVisTauPt' : 1.0 
  },
  'minVisTauPtGt1p5' : {
    'minVisTauPt' : 1.5 
  },
  'minVisTauPtGt2p0' : {
    'minVisTauPt' : 2.0 
  },
  'minVisTauPtGt2p5' : {
    'minVisTauPt' : 2.5 
  },
  'minVisTauPtGt3p0' : {
    'minVisTauPt' : 3.0 
  },
  'noVisTauEtaCut' : {
    'maxAbsVisTauEta' : -1.
  },
  'maxAbsVisTauEtaLt1p5' : {
    'maxAbsVisTauEta' : 1.5
  },
  'maxAbsVisTauEtaLt2p5' : {
    'maxAbsVisTauEta' : 2.5
  },
  'minVisTauZgt0p4' : {
    'minVisTauZ' : 0.4 
  },
  'minVisTauZgt0p5' : {
    'minVisTauZ' : 0.5 
  },
  'minVisTauZgt0p6' : {
    'minVisTauZ' : 0.6 
  },
  'noChargedKaonCut' : {
    'maxNumChargedKaons' : -1
  },
  'maxNumChargedKaonsEq0' : {
    'maxNumChargedKaons' : 0
  },
  'noNeutralKaonCut' : {
    'maxNumNeutralKaons' : -1
  },
  'maxNumNeutralKaonsEq0' : {
    'maxNumNeutralKaons' : 0
  },
  'noPhotonEnCut' : {
    'maxSumPhotonEn' : -1. 
  },
  'maxSumPhotonEnEq0' : {
    'maxSumPhotonEn' : 0. 
  },
  'maxSumPhotonEnLt0p1' : {
    'maxSumPhotonEn' : 0.1 
  },
  'maxSumPhotonEnLt0p5' : {
    'maxSumPhotonEn' : 0.5 
  }
}

inputFilePath = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/ntuples/", collider, version)

configDir  = os.path.join("/home",               getpass.getuser(), "Entanglement/studyRecoEffects/", collider, version)
outputDir  = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/studyRecoEffects/", collider, version)
workingDir = os.getcwd()
cmsswDir   = os.getenv('CMSSW_BASE')

def run_command(command):
  #print("executing command = '%s'" % command)
  os.system(command)

run_command('mkdir -p %s' % configDir)
run_command('mkdir -p %s' % os.path.join(configDir, "plots"))
run_command('mkdir -p %s' % outputDir)

def build_cfgFile(cfgFile_original, cfgFile_modified, 
                  inputFileNames, process, par_gen,
                  mode, collider, hAxis, decayMode, apply_evtWeight, spinAnalyzer,
                  replacements,
                  outputFileName):
  print("Building configFile = '%s'" % cfgFile_modified)

  rmCommand   = 'rm -f %s' % cfgFile_modified
  run_command(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##par_gen/par_gen/; s/\$par_gen/%s/;' % par_gen
  sedCommand += '  s/##mode/mode/; s/\$mode/%s/;' % mode
  sedCommand += '  s/##collider/collider/; s/\$collider/%s/;' % collider
  sedCommand += '  s/##hAxis/hAxis/; s/\$hAxis/%s/;' % hAxis
  sedCommand += '  s/##decayMode/decayMode/; s/\$decayMode/%s/;' % decayMode
  sedCommand += '  s/##apply_evtWeight/apply_evtWeight/; s/\$apply_evtWeight/%s/;' % apply_evtWeight
  sedCommand += '  s/##spinAnalyzer/spinAnalyzer/; s/\$spinAnalyzer/%s/;' % spinAnalyzer
  for replacement_key, replacement_value in replacements.items():
      sedCommand += '  s/##%s/%s/; s/\$%s/%s/;' % (replacement_key, replacement_key, replacement_key, replacement_value)
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName.replace("/", "\/")
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  run_command(sedCommand)
 
jobOptions_analysis  = {} # key = sample, mode, hAxis, recoEffect_key
jobOptions_ctrlPlots = {} # key = sample, mode, hAxis, recoEffect_key
for sampleName, sample in samples.items():
  for hAxis in hAxes:
    print("processing sample = '%s', hAxis = '%s'" % (sampleName, hAxis))
    print(" inputFilePath = '%s'" % inputFilePath)
    inputFile_regex = r"entanglementNtuple_%s_%sAxis_[0-9]+.root" % (sampleName, hAxis)
    inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
    numInputFiles = len(inputFileNames)
    print("Found %i input files." % numInputFiles)
    for mode in modes:
      for decayMode in decayModes:
        for recoEffect_key, recoEffect_value in recoEffects.items():
          if mode != "gen_smeared":
            for spinAnalyzer in spinAnalyzers:
              cfgFileName_analysis_modified = os.path.join(configDir, "analyzeEntanglementNtuple_%s_%sMode_%sAxis_%sDecayMode_%s_%s_cfg.py" % \
                (sampleName, mode, hAxis, decayMode, spinAnalyzer, recoEffect_key))
              outputFileName_analysis = "analyzeEntanglementNtuple_%s_%sMode_%sAxis_%sDecayMode_%s_%s.root" % \
                (sampleName, mode, hAxis, decayMode, spinAnalyzer, recoEffect_key)
              build_cfgFile(
                "analyzeEntanglementNtuple_cfg.py", cfgFileName_analysis_modified, 
                inputFileNames, sample['process'], sample['par_gen'],
                mode, collider, hAxis, decayMode, sample['apply_evtWeight'], spinAnalyzer, 
                recoEffect_value,
                outputFileName_analysis)
              logFileName_analysis = cfgFileName_analysis_modified.replace("_cfg.py", ".log")
              job_key_analysis = '%s_%s_%s_%s_%s_%s_analysis' % (sampleName, mode, hAxis, decayMode, spinAnalyzer, recoEffect_key)
              dependencies_analysis = [ cfgFileName_analysis_modified ]
              dependencies_analysis.extend(inputFileNames)
              jobOptions_analysis[job_key_analysis] = {
                'inputFileNames' : dependencies_analysis,
                'cfgFileName'    : cfgFileName_analysis_modified,
                'outputFilePath' : outputDir,
                'outputFileName' : outputFileName_analysis,
                'logFileName'    : logFileName_analysis,
              }
          cfgFileName_ctrlPlots_modified = os.path.join(configDir, "makeControlPlots_%s_%sMode_%sAxis_%sDecayMode_%s_cfg.py" % \
            (sampleName, mode, hAxis, decayMode, recoEffect_key))
          outputFileName_ctrlPlots = os.path.join(configDir, "plots", "makeControlPlots_%s_%sMode_%sAxis_%sDecayMode_%s.root" % \
            (sampleName, mode, hAxis, decayMode, recoEffect_key))
          build_cfgFile(
            "makeControlPlots_cfg.py", cfgFileName_ctrlPlots_modified, 
            inputFileNames, sample['process'], sample['par_gen'],
            mode, collider, hAxis, decayMode, sample['apply_evtWeight'], "",
            recoEffect_value,
            outputFileName_ctrlPlots)
          logFileName_ctrlPlots = cfgFileName_ctrlPlots_modified.replace("_cfg.py", ".log")
          job_key_ctrlPlots = '%s_%s_%s_%s_%s_ctrlPlots' % (sampleName, mode, hAxis, decayMode, recoEffect_key)
          dependencies_ctrlPlots = [ cfgFileName_ctrlPlots_modified ]
          dependencies_ctrlPlots.extend(inputFileNames) 
          jobOptions_ctrlPlots[job_key_ctrlPlots] = {
            'inputFileNames' : dependencies_ctrlPlots,
            'cfgFileName'    : cfgFileName_ctrlPlots_modified,
            'outputFilePath' : outputDir,
            'outputFileName' : outputFileName_ctrlPlots,
            'logFileName'    : logFileName_ctrlPlots,
          }

jobOptions_Makefile = []
for job_key, job in jobOptions_analysis.items():
  commands = []
  commands.append('rm -f %s' % job['outputFileName'])
  commands.append('rm -f %s' % job['logFileName'])
  commands.append('analyzeEntanglementNtuple %s >& %s' % (job['cfgFileName'], job['logFileName']))
  commands.append('cp %s %s' % (job['outputFileName'], os.path.join(outputDir, job['outputFileName'])))
  commands.append('rm -f %s' % job['outputFileName'])
  jobOptions_Makefile.append({
    'target'          : os.path.join(outputDir, job['outputFileName']),
    'dependencies'    : [ inputFileName.replace("file:", "") for inputFileName in job['inputFileNames'] ],
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, job['outputFileName']) ],
  })
for job_key, job in jobOptions_ctrlPlots.items():
  commands = []
  commands.append('rm -f %s' % job['outputFileName'])
  commands.append('rm -f %s' % job['logFileName'])
  commands.append('makeControlPlots %s >& %s' % (job['cfgFileName'], job['logFileName']))
  jobOptions_Makefile.append({
    'target'          : os.path.join(outputDir, job['outputFileName']),
    'dependencies'    : [ inputFileName.replace("file:", "") for inputFileName in job['inputFileNames'] ],
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, job['outputFileName']) ],
  })
makeFileName = os.path.join(configDir, "Makefile")
build_Makefile(makeFileName, jobOptions_Makefile)

message  = "Finished building config files."
message += " Now execute 'make -j 12 -f %s' to start the jobs." % makeFileName
print(message)
