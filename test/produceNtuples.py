#!/usr/bin/env python

import getpass
import os

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames, build_Makefile

samples = {
  'ggH_htt_pythia8' : {
    'inputFilePath' : '/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/100000/',
    'numJobs' : 10,
    'process' : "ggH_htt_pythia8"
  },
  ##'ggH_htt_tauola' : {
  ##  'inputFilePath' : '/store/mc/RunIISummer16MiniAODv2/GluGluToHToTauTau_M125_13TeV_amcatnloFXFX_pythia8_tauola/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/',
  ##  'numJobs' : 10,
  ##  'process' : "ggH_htt_tauola"
  ##},
  #'dy_lo_pythia8' : {
  #  'inputFilePath' : '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/270000/',
  #  'numJobs' : 10,
  #  'process' : "dy_lo_pythia8"
  #},
  #'dy_nlo_pythia8' : {
  #  'inputFilePath' : '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/230000/',
  #  'numJobs' : 10,
  #  'process' : "dy_nlo_pythia8"
  #},
}

hAxes = [ "beam", "higgs" ]

version = "2023Aug17_startPosMode1_woSmearing"

configDir  = os.path.join("/home",               getpass.getuser(), "Entanglement/ntuples/", version)
outputDir  = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/ntuples/", version)
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
                  rndSeed,
                  outputFileName):
  print("Building configFile = '%s'" % cfgFile_modified)
  #print(" rndSeed = %i" % rndSeed)

  rmCommand   = 'rm -f %s' % cfgFile_modified
  run_command(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##hAxis/hAxis/; s/\$hAxis/%s/;' % hAxis
  sedCommand += '  s/##rndSeed/rndSeed/; s/\$rndSeed/%i/;' % rndSeed
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  run_command(sedCommand)
 
jobOptions = {} # key = process, hAxis, jobId
for sampleName, sample in samples.items():
  print("processing sample = '%s'" % sampleName)
  process = sample['process']
  inputFilePath = sample['inputFilePath']
  print(" inputFilePath = '%s'" % inputFilePath)
  inputFileNames = getInputFileNames(inputFilePath)
  numInputFiles = len(inputFileNames)
  print("Found %i input files." % numInputFiles)
  numJobs = sample['numJobs']
  for hAxis in hAxes:
    for jobId in range(numJobs):
      idxFirstFile = jobId*numInputFiles/numJobs
      idxLastFile = (jobId + 1)*numInputFiles/numJobs - 1
      inputFileNames_job = inputFileNames[idxFirstFile:idxLastFile + 1]
      cfgFileName_modified = os.path.join(configDir, "produceEntanglementNtuple_%s_%sAxis_%i_cfg.py" % \
        (sampleName, hAxis, jobId))
      rndSeed = jobId + 1
      outputFileName = "entanglementNtuple_%s_%sAxis_%i.root" % \
        (sampleName, hAxis, jobId)
      build_cfgFile(
        "produceEntanglementNtuple_cfg.py", cfgFileName_modified, 
        inputFileNames_job, sample['process'],
        hAxis,
        rndSeed,
        outputFileName)
      logFileName = cfgFileName_modified.replace("_cfg.py", ".log")
      job_key = '%s_%s_%i' % (process, hAxis, jobId)
      dependencies = [ cfgFileName_modified ]
      dependencies.extend(inputFileNames_job)
      jobOptions[job_key] = {
        'inputFileNames' : dependencies,
        'cfgFileName'    : cfgFileName_modified,
        'outputFilePath' : outputDir,
        'outputFileName' : outputFileName,
        'logFileName'    : logFileName,
      }

jobOptions_Makefile = []
for job_key, job in jobOptions.items():
  commands = []
  commands.append('rm -f %s' % job['outputFileName'])
  commands.append('rm -f %s' % job['logFileName'])
  commands.append('cmsRun %s >& %s' % (job['cfgFileName'], job['logFileName']))
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
