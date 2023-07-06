#!/usr/bin/env python

import getpass
import os

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames, build_Makefile

#samples = [ 'ggH_htt_pythia8', 'ggH_htt_tauola', 'dy_lo_pythia8', 'dy_nlo_pythia8' ]
samples = [ 'ggH_htt_pythia8' ]

modes = [ "gen", "rec" ]
#hAxes = [ "beam", "higgs" ]
hAxes = [ "higgs" ]

version = "2023Jun16"

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
                  mode, hAxis,
                  outputFileName):
  print("Building configFile = '%s'" % cfgFile_modified)

  rmCommand   = 'rm -f %s' % cfgFile_modified
  run_command(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##mode/mode/; s/\$mode/%s/;' % mode
  sedCommand += '  s/##hAxis/hAxis/; s/\$hAxis/%s/;' % hAxis
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  run_command(sedCommand)
 
jobOptions_analysis  = {} # key = sample, mode, hAxis
jobOptions_ctrlPlots = {} # key = sample, mode, hAxis
for sample in samples:
  for mode in modes:
    for hAxis in hAxes:
      print("processing sample = '%s', mode = '%s', hAxis = '%s'" % (sample, mode, hAxis))
      print(" inputFilePath = '%s'" % inputFilePath)
      inputFile_regex = r"entanglementNtuple_%s_%sMode_%sAxis_[0-9]+.root" % (sample, mode, hAxis)
      inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
      numInputFiles = len(inputFileNames)
      print("Found %i input files." % numInputFiles)
      cfgFileName_analysis_modified = os.path.join(configDir, "analyzeEntanglementNtuple_%s_%sMode_%sAxis_cfg.py" % \
        (sample, mode, hAxis))
      outputFileName_analysis = "analyzeEntanglementNtuple_%s_%sMode_%sAxis.root" % \
        (sample, mode, hAxis)
      build_cfgFile(
        "analyzeEntanglementNtuple_cfg.py", cfgFileName_analysis_modified, 
        inputFileNames, sample,
        mode, hAxis, 
        outputFileName_analysis)
      logFileName_analysis = cfgFileName_analysis_modified.replace("_cfg.py", ".log")
      job_key_analysis = '%s_%s_%s_analysis' % (sample, mode, hAxis)
      jobOptions_analysis[job_key_analysis] = {
        'inputFileNames' : inputFileNames,
        'cfgFileName'    : cfgFileName_analysis_modified,
        'outputFilePath' : outputDir,
        'outputFileName' : outputFileName_analysis,
        'logFileName'    : logFileName_analysis,
      }
      cfgFileName_ctrlPlots_modified = os.path.join(configDir, "makeControlPlots_%s_%sMode_%sAxis_cfg.py" % \
        (sample, mode, hAxis))
      outputFileName_ctrlPlots = "makeControlPlots_%s_%sMode_%sAxis.root" % \
        (sample, mode, hAxis)
      build_cfgFile(
        "makeControlPlots_cfg.py", cfgFileName_ctrlPlots_modified, 
        inputFileNames, sample,
        mode, hAxis, 
        outputFileName_ctrlPlots)
      logFileName_ctrlPlots = cfgFileName_ctrlPlots_modified.replace("_cfg.py", ".log")
      job_key_ctrlPlots = '%s_%s_%s_ctrlPlots' % (sample, mode, hAxis)
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
