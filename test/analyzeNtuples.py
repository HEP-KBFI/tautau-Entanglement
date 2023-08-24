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
decayModes = [ "piPlus_piMinus", "rhoPlus_rhoMinus" ]
spinAnalyzers = [ "by_summation", "by_mlfit" ]
#spinAnalyzers = [ "by_summation" ]

version = "2023Aug24_woSmearing"

samples = None
if collider == "LHC":
    samples = samples_LHC
elif collider == "SuperKEKB":
    samples = samples_SuperKEKB
else:
    raise ValueError("Invalid Configuration parameter 'collider' = '%s' !!" % collider)

inputFilePath = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/ntuples/", collider, version)

configDir  = os.path.join("/home",               getpass.getuser(), "Entanglement/analysis/", collider, version)
outputDir  = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/analysis/", collider, version)
workingDir = os.getcwd()
cmsswDir   = os.getenv('CMSSW_BASE')

def run_command(command):
  #print("executing command = '%s'" % command)
  os.system(command)

run_command('mkdir -p %s' % configDir)
run_command('mkdir -p %s' % outputDir)

def build_cfgFile(cfgFile_original, cfgFile_modified, 
                  inputFileNames, process,
                  mode, collider, hAxis, decayMode, apply_evtWeight, spinAnalyzer,
                  outputFileName):
  print("Building configFile = '%s'" % cfgFile_modified)

  rmCommand   = 'rm -f %s' % cfgFile_modified
  run_command(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##mode/mode/; s/\$mode/%s/;' % mode
  sedCommand += '  s/##collider/collider/; s/\$collider/%s/;' % collider
  sedCommand += '  s/##hAxis/hAxis/; s/\$hAxis/%s/;' % hAxis
  sedCommand += '  s/##decayMode/decayMode/; s/\$decayMode/%s/;' % decayMode
  sedCommand += '  s/##apply_evtWeight/apply_evtWeight/; s/\$apply_evtWeight/%s/;' % apply_evtWeight
  sedCommand += '  s/##spinAnalyzer/spinAnalyzer/; s/\$spinAnalyzer/%s/;' % spinAnalyzer
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  run_command(sedCommand)
 
jobOptions_analysis  = {} # key = sample, mode, hAxis
jobOptions_ctrlPlots = {} # key = sample, mode, hAxis
jobOptions_resPlots  = {} # key = sample, mode, hAxis
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
        for spinAnalyzer in spinAnalyzers:
          cfgFileName_analysis_modified = os.path.join(configDir, "analyzeEntanglementNtuple_%s_%sMode_%sAxis_%sDecayMode_%s_cfg.py" % \
            (sampleName, mode, hAxis, decayMode, spinAnalyzer))
          outputFileName_analysis = "analyzeEntanglementNtuple_%s_%sMode_%sAxis_%sDecayMode_%s.root" % \
            (sampleName, mode, hAxis, decayMode, spinAnalyzer)
          build_cfgFile(
            "analyzeEntanglementNtuple_cfg.py", cfgFileName_analysis_modified, 
            inputFileNames, sample['process'],
            mode, collider, hAxis, decayMode, sample['apply_evtWeight'], spinAnalyzer, 
            outputFileName_analysis)
          logFileName_analysis = cfgFileName_analysis_modified.replace("_cfg.py", ".log")
          job_key_analysis = '%s_%s_%s_%s_%s_analysis' % (sampleName, mode, hAxis, decayMode, spinAnalyzer)
          dependencies_analysis = [ cfgFileName_analysis_modified ]
          dependencies_analysis.extend(inputFileNames)
          jobOptions_analysis[job_key_analysis] = {
            'inputFileNames' : dependencies_analysis,
            'cfgFileName'    : cfgFileName_analysis_modified,
            'outputFilePath' : outputDir,
            'outputFileName' : outputFileName_analysis,
            'logFileName'    : logFileName_analysis,
          }
        cfgFileName_ctrlPlots_modified = os.path.join(configDir, "makeControlPlots_%s_%sMode_%sAxis_%sDecayMode_cfg.py" % \
          (sampleName, mode, hAxis, decayMode))
        outputFileName_ctrlPlots = "makeControlPlots_%s_%sMode_%sAxis_%sDecayMode.root" % \
          (sampleName, mode, hAxis, decayMode)
        build_cfgFile(
          "makeControlPlots_cfg.py", cfgFileName_ctrlPlots_modified, 
          inputFileNames, sample['process'],
          mode, collider, hAxis, decayMode, sample['apply_evtWeight'], "",
          outputFileName_ctrlPlots)
        logFileName_ctrlPlots = cfgFileName_ctrlPlots_modified.replace("_cfg.py", ".log")
        job_key_ctrlPlots = '%s_%s_%s_%s_ctrlPlots' % (sampleName, mode, hAxis, decayMode)
        dependencies_ctrlPlots = [ cfgFileName_ctrlPlots_modified ]
        dependencies_ctrlPlots.extend(inputFileNames) 
        jobOptions_ctrlPlots[job_key_ctrlPlots] = {
          'inputFileNames' : dependencies_ctrlPlots,
          'cfgFileName'    : cfgFileName_ctrlPlots_modified,
          'outputFilePath' : outputDir,
          'outputFileName' : outputFileName_ctrlPlots,
          'logFileName'    : logFileName_ctrlPlots,
        }
        if mode != "gen":
          cfgFileName_resPlots_modified = os.path.join(configDir, "makeResolutionPlots_%s_%sMode_%sAxis_%sDecayMode_cfg.py" % \
            (sampleName, mode, hAxis, decayMode))
          outputFileName_resPlots = "makeResolutionPlots_%s_%sMode_%sAxis_%sDecayMode.root" % \
            (sampleName, mode, hAxis, decayMode)
          build_cfgFile(
            "makeResolutionPlots_cfg.py", cfgFileName_resPlots_modified, 
            inputFileNames, sample['process'],
            mode, collider, hAxis, decayMode, sample['apply_evtWeight'], "", 
            outputFileName_resPlots)
          logFileName_resPlots = cfgFileName_resPlots_modified.replace("_cfg.py", ".log")
          job_key_resPlots = '%s_%s_%s_%s_resPlots' % (sampleName, mode, hAxis, decayMode)
          dependencies_resPlots = [ cfgFileName_resPlots_modified ]
          dependencies_resPlots.extend(inputFileNames)
          jobOptions_resPlots[job_key_resPlots] = {
            'inputFileNames' : dependencies_resPlots,
            'cfgFileName'    : cfgFileName_resPlots_modified,
            'outputFilePath' : outputDir,
            'outputFileName' : outputFileName_resPlots,
            'logFileName'    : logFileName_resPlots,
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
for job_key, job in jobOptions_resPlots.items():
  commands = []
  commands.append('rm -f %s' % job['outputFileName'])
  commands.append('rm -f %s' % job['logFileName'])
  commands.append('makeResolutionPlots %s >& %s' % (job['cfgFileName'], job['logFileName']))
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
