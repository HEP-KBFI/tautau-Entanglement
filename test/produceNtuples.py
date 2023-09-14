#!/usr/bin/env python3

import getpass
import os

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames, build_Makefile
from TauAnalysis.Entanglement.samples import samples_LHC, samples_SuperKEKB

hAxes = [ "beam", "higgs" ]
#hAxes = [ "beam" ]
#collider = "LHC"
collider = "SuperKEKB"

version = "2023Sep14_wSmearing"

samples = None
if collider == "LHC":
    samples = samples_LHC
elif collider == "SuperKEKB":
    samples = samples_SuperKEKB
else:
    raise ValueError("Invalid Configuration parameter 'collider' = '%s' !!" % collider)

configDir  = os.path.join("/home",               getpass.getuser(), "Entanglement/ntuples/", collider, version)
outputDir  = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/ntuples/", collider, version)
workingDir = os.getcwd()
cmsswDir   = os.getenv('CMSSW_BASE')

def run_command(command):
  #print("executing command = '%s'" % command)
  os.system(command)

run_command('mkdir -p %s' % configDir)
run_command('mkdir -p %s' % outputDir)

def build_cfgFile(cfgFile_original, cfgFile_modified, 
                  inputFileNames, process,
                  collider, hAxis,
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
  sedCommand += '  s/##collider/collider/; s/\$collider/%s/;' % collider
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
      idxFirstFile = int(jobId*numInputFiles/numJobs)
      idxLastFile = int((jobId + 1)*numInputFiles/numJobs - 1)
      inputFileNames_job = inputFileNames[idxFirstFile:idxLastFile + 1]
      cfgFileName_modified = os.path.join(configDir, "produceEntanglementNtuple_%s_%sAxis_%i_cfg.py" % \
        (sampleName, hAxis, jobId))
      rndSeed = jobId + 1
      outputFileName = "entanglementNtuple_%s_%sAxis_%i.root" % \
        (sampleName, hAxis, jobId)
      build_cfgFile(
        "produceEntanglementNtuple_cfg.py", cfgFileName_modified, 
        inputFileNames_job, sample['process'],
        collider, hAxis,
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
