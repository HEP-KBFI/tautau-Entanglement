#!/usr/bin/env python

import getpass
import os

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames, build_Makefile

samples = [ 'ggH_htt_pythia8' ]

hAxes = [ "beam" ]

version = "2023Jun02"

recoEffects = {
  'noVisTauPtCut' : {
    'minVisTauPt' : -1. 
  },
  'minVisTauPtGt10' : {
    'minVisTauPt' : 10. 
  },
  'minVisTauPtGt20' : {
    'minVisTauPt' : 20. 
  },
  'minVisTauPtGt30' : {
    'minVisTauPt' : 30. 
  },
  'minVisTauPtGt40' : {
    'minVisTauPt' : 40. 
  },
  'noVisTauEtaCut' : {
    'maxAbsVisTauEta' : -1.
  },
  'maxAbsVisTauEtaLt2p0' : {
    'maxAbsVisTauEta' : 2.
  },
  'maxAbsVisTauEtaLt3p0' : {
    'maxAbsVisTauEta' : 3.
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
  'maxSumPhotonEnLt2' : {
    'maxSumPhotonEn' : 2. 
  },
  'maxSumPhotonEnLt5' : {
    'maxSumPhotonEn' : 5. 
  }
}

inputFilePath = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/ntuples/", version)

configDir  = os.path.join("/home",               getpass.getuser(), "Entanglement/studyRecoEffects/", version)
outputDir  = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/studyRecoEffects/", version)
workingDir = os.getcwd()
cmsswDir   = os.getenv('CMSSW_BASE')

def run_command(command):
  #print("executing command = '%s'" % command)
  os.system(command)

run_command('mkdir -p %s' % configDir)
run_command('mkdir -p %s' % outputDir)

def build_cfgFile(cfgFile_original, cfgFile_modified, 
                  inputFileNames, process,
                  hAxis,
                  replacements,
                  outputFileName):
  print("Building configFile = '%s'" % cfgFile_modified)

  rmCommand   = 'rm -f %s' % cfgFile_modified
  run_command(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##hAxis/hAxis/; s/\$hAxis/%s/;' % hAxis
  for replacement_key, replacement_value in replacements.items():
      sedCommand += '  s/##%s/%s/; s/\$%s/%s/;' % (replacement_key, replacement_key, replacement_key, replacement_value)
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  run_command(sedCommand)
 
jobOptions_analysis  = {} # key = sample, hAxis, recoEffect_key
jobOptions_ctrlPlots = {} # key = sample, hAxis, recoEffect_key
for sample in samples:  
  for hAxis in hAxes:
    print("processing sample = '%s', hAxis = '%s'" % (sample, hAxis))
    print(" inputFilePath = '%s'" % inputFilePath)
    inputFile_regex = r"entanglementNtuple_%s_%sAxis_[0-9]+.root" % (sample, hAxis)
    inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
    numInputFiles = len(inputFileNames)
    print("Found %i input files." % numInputFiles)
    for recoEffect_key, recoEffect_value in recoEffects.items():  
      cfgFileName_analysis_modified = os.path.join(configDir, "analyzeEntanglementNtuple_%s_%sAxis_%s_cfg.py" % \
        (sample, hAxis, recoEffect_key))
      outputFileName_analysis = "analyzeEntanglementNtuple_%s_%sAxis_%s.root" % \
        (sample, hAxis, recoEffect_key)
      build_cfgFile(
        "analyzeEntanglementNtuple_cfg.py", cfgFileName_analysis_modified, 
        inputFileNames, sample,
        hAxis,
        recoEffect_value,
        outputFileName_analysis)
      logFileName_analysis = cfgFileName_analysis_modified.replace("_cfg.py", ".log")
      job_key_analysis = '%s_%s_%s_analysis' % (sample, hAxis, recoEffect_key)
      jobOptions_analysis[job_key_analysis] = {
        'inputFileNames' : inputFileNames,
        'cfgFileName'    : cfgFileName_analysis_modified,
        'outputFilePath' : outputDir,
        'outputFileName' : outputFileName_analysis,
        'logFileName'    : logFileName_analysis,
      }
      cfgFileName_ctrlPlots_modified = os.path.join(configDir, "makeControlPlots_%s_%sAxis_%s_cfg.py" % \
        (sample, hAxis, recoEffect_key))
      outputFileName_ctrlPlots = "makeControlPlots_%s_%sAxis_%s.root" % \
        (sample, hAxis, recoEffect_key)
      build_cfgFile(
        "makeControlPlots_cfg.py", cfgFileName_ctrlPlots_modified, 
        inputFileNames, sample,
        hAxis,
        recoEffect_value,
        outputFileName_ctrlPlots)
      logFileName_ctrlPlots = cfgFileName_ctrlPlots_modified.replace("_cfg.py", ".log")
      job_key_ctrlPlots = '%s_%s_%s_ctrlPlots' % (sample, hAxis, recoEffect_key)
      jobOptions_ctrlPlots[job_key_ctrlPlots] = {
        'inputFileNames' : inputFileNames,
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
  commands.append('cp %s %s' % (job['outputFileName'], os.path.join(outputDir, job['outputFileName'])))
  commands.append('rm -f %s' % job['outputFileName'])
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
