#!/usr/bin/env python

import getpass
import os

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames, build_Makefile

samples = [ 'ggH_htt_pythia8', 'ggH_htt_tauola', 'dy_lo_pythia8', 'dy_nlo_pythia8' ]

hAxes = [ "beam", "higgs" ]

version = "2023Jun02"

inputFilePath = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/ntuples/", version)

configDir  = os.path.join("/home",               getpass.getuser(), "Entanglement/analysis/", version)
outputDir  = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/analysis/", version)
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
                  outputFileName):
  print("Building configFile = '%s'" % cfgFileName_modified)

  rmCommand   = 'rm -f %s' % cfgFileName_modified
  run_command(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##hAxis/hAxis/; s/\$hAxis/%s/;' % hAxis
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  run_command(sedCommand)
 
jobOptions_analysis  = {} # key = sample, hAxis
jobOptions_ctrlPlots = {} # key = sample, hAxis
for sample in samples:  
  for hAxis in hAxes:
    print("processing sample = '%s', hAxis = '%s'" % (sample, hAxis))
    print(" inputFilePath = '%s'" % inputFilePath)
    inputFile_regex = r"entanglementNtuple_%s_%sAxis_[0-9]+.root" % (sample, hAxis)
    inputFileNames = getInputFileNames(inputFilePath)
    numInputFiles = len(inputFileNames)
    print("Found %i input files." % numInputFiles)
    cfgFileName_analysis_modified = os.path.join(configDir, "analyzeEntanglementNtuple_%s_%sAxis_cfg.py" % \
      (sample, hAxis))
    outputFileName_analysis = "analyzeEntanglementNtuple_%s_%s.root" % \
      (sample, hAxis)
    build_cfgFile(
      "analyzeEntanglementNtuple_cfg.py", cfgFileName_analysis_modified, 
      inputFileNames_job, sample,
      hAxis, 
      outputFileName)
    logFileName_analysis = cfgFileName_analysis_modified.replace("_cfg.py", ".log")
    job_key_analysis = '%s_%s_analysis' % (sample, hAxis)
    jobOptions_analysis[job_key_analysis] = {
      'inputFileNames' : inputFileNames_job,
      'cfgFileName'    : cfgFileName_modified,
      'outputFilePath' : outputDir,
      'outputFileName' : outputFileName_analysis,
      'logFileName'    : logFileName_analysis,
    }
    cfgFileName_ctrlPlots_modified = os.path.join(configDir, "makeControlPlots_%s_%sAxis_cfg.py" % \
      (sample, hAxis))
    outputFileName_ctrlPlots = "makeControlPlots_%s_%s.root" % \
      (sample, hAxis)
    build_cfgFile(
      "makeControlPlots_cfg.py", cfgFileName_ctrlPlots_modified, 
      inputFileNames_job, sample,
      hAxis, 
      outputFileName)
    logFileName_analysis = cfgFileName_ctrlPlots_modified.replace("_cfg.py", ".log")
    job_key_ctrlPlots = '%s_%s_ctrlPlots' % (sample, hAxis)
    jobOptions_ctrlPlots[job_key_ctrlPlots] = {
      'inputFileNames' : inputFileNames_job,
      'cfgFileName'    : cfgFileName_modified,
      'outputFilePath' : outputDir,
      'outputFileName' : outputFileName_analysis,
      'logFileName'    : logFileName_analysis,
    }

jobOptions_Makefile = []
for job_key, job in jobOptions_analysis.items():
  commands = []
  commands.append('rm -f %s' % job['outputFileName'])
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
