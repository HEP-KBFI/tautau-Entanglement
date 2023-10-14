#!/usr/bin/env python3

import datetime
import argparse
import getpass
import sys
import os

# Example usage:
# ./test/produceNtuples.py -v 2023Oct06 -s dy_lo_pythia8

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames, build_Makefile, query_yes_no, \
  build_cfg, mkdir, read_contents, save_cmd, positive_int_type, build_sbatchSubmission

mode_choices = [ "gen", "gen_smeared", "startPos", "kinFit" ]
hAxes_choices = [ "beam", "higgs" ]
collider_choices = [ "LHC", "SuperKEKB" ]
decayMode_choices = [ "piPlus_piMinus", "piPlus_rhoMinus", "rhoPlus_piMinus", "rhoPlus_rhoMinus" ] #TODO to be changed
spinAnalyzer_choices = [ "by_summation", "by_mlfit", "by_differentialXsec1d", "by_differentialXsec2d", "by_asymmetry" ]

parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--version', type = str, required = True, help = f'Version, e.g. {datetime.datetime.now().strftime("%Y%b%d")}')
parser.add_argument('-V', '--ntuple-version', type = str, default = '', help = 'Version of input Ntuples (not needed if the same as --version)')
parser.add_argument('-c', '--collider', type = str, choices = collider_choices, default = 'SuperKEKB', help = 'Collider')
parser.add_argument('-a', '--axes', nargs = '*', type = str, choices = collider_choices, default = [ 'beam' ], help = 'Axes')
parser.add_argument('-S', '--spin-analyzers', nargs = '*', type = str, choices = spinAnalyzer_choices, default = spinAnalyzer_choices, help = 'Spin analyzers')
parser.add_argument('-d', '--decay-modes', nargs = '*', type = str, choices = decayMode_choices, default = decayMode_choices, help = 'Tau decay modes')
parser.add_argument('-m', '--modes', nargs = '*', type = str, choices = mode_choices, default = mode_choices, help = 'Input data')
parser.add_argument('-s', '--samples', nargs = '*', default = [], help = 'Whitelisted samples')
parser.add_argument('-b', '--bootstrap-size', type = int, default = 2000, help = 'Size of bootstrap dataset (use -1 to consider all events from the input sample)')
parser.add_argument('-B', '--bootstrap-count', type = positive_int_type, default = 1000, help = 'Number of bootstrap datasets')
parser.add_argument('-j', '--job-type', type = str, choices = ['local', 'cluster'], required = True, help = 'Job type')
args = parser.parse_args()

version = args.version
ntuple_version = args.ntuple_version
if not ntuple_version:
  ntuple_version = version

hAxes = args.axes
collider = args.collider
modes = args.modes
decayModes = args.decay_modes
spinAnalyzers = args.spin_analyzers
whitelist = args.samples
bootstrap_size = args.bootstrap_size
bootstrap_count = args.bootstrap_count
run_makefile = args.job_type == 'local'

if collider == "LHC":
    from TauAnalysis.Entanglement.samples import samples_LHC as samples, PAR_GEN_LHC as par_gen, \
      MAX_SUM_PHOTON_EN_LHC as maxSumPhotonEn
elif collider == "SuperKEKB":
    from TauAnalysis.Entanglement.samples import samples_SuperKEKB as samples, PAR_GEN_SUPERKEKB as par_gen, \
      MAX_SUM_PHOTON_EN_SUPERKEKB as maxSumPhotonEn
else:
    assert(False)
# Karl: ignore FSR cut, see commit d5e08c353cb1d91639f804198e349667ec312e16
maxSumPhotonEn = -1

if not whitelist:
  run_all_samples = query_yes_no(
    "Do you really want to process all those samples: {}?".format(', '.join(samples.keys())),
    default = "no",
  )
  if not run_all_samples:
    sys.exit(0)

if not whitelist:
  run_all_samples = query_yes_no(
    "Do you really want to process all those samples: {}?".format(', '.join(samples.keys())),
    default = "no",
  )
  if not run_all_samples:
    sys.exit(0)

inputFilePath = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/ntuples", collider, ntuple_version)

configDir  = os.path.join("/home",               getpass.getuser(), "Entanglement/analysis", collider, version)
outputDir  = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/analysis", collider, version)
testDir    = os.path.dirname(os.path.abspath(__file__))
cmsswDir   = os.getenv('CMSSW_BASE')

outputDir_ctrlPlots = os.path.join(configDir, "plots")
outputDir_resPlots = os.path.join(configDir, "plots")

analyzeNtuple_template = read_contents(os.path.join(testDir, "analyzeEntanglementNtuple_cfg.py"))
makeControlPlot_template = read_contents(os.path.join(testDir, "makeControlPlots_cfg.py"))
makeResolutionPlot_template = read_contents(os.path.join(testDir, "makeResolutionPlots_cfg.py"))

analyzeEntanglementNtuple_cmd = 'analyzeEntanglementNtuple'
makeControlPlots_cmd = 'makeControlPlots'
makeResolutionPlots_cmd = 'makeResolutionPlots'

mkdir(configDir)
mkdir(os.path.join(configDir, "plots"))
mkdir(outputDir)

if not save_cmd(os.path.join(configDir, "cmd.txt")):
  sys.exit(0)

def init_dict(dictionary, keys):
    """Auxiliary function to initialize dictionary for access with multiple levels of keys
    """
    dictionary_at_keylevel = dictionary
    numKeys = len(keys)
    for idxKey in range(numKeys - 1):
        key = keys[idxKey]
        if key not in dictionary_at_keylevel.keys():
            dictionary_at_keylevel[key] = {}
        dictionary_at_keylevel = dictionary_at_keylevel[key]

def get_all_values(dictionary):
    """Auxiliary function to get all values stored in dictionary with multiple levels of keys
    """
    all_values = []
    for key in dictionary.keys():
        if isinstance(dictionary[key], dict):
            values = get_all_values(dictionary[key])
            all_values.extend(values)
        else:
            all_values.append(dictionary[key])
    return all_values

jobOptions_analysis              = {} # key = job_key_analysis
jobOptions_ctrlPlots             = {} # key = job_key_ctrlPlots
jobOptions_resPlots              = {} # key = job_key_resPlots

# use list of all previous writeSelBiasCorrection jobs as dependencies for next writeSelBiasCorrection job, to ensure that
# only one writeSelBiasCorrection job runs at a time; also use all writeSelBiasCorrection jobs as dependencies for all analysis jobs
outputFiles_writeSelBiasCorrection = []
for sampleName, sample in samples.items():
  if whitelist and sampleName not in whitelist:
    continue
  for hAxis in hAxes:
    print(f"processing sample = '{sampleName}', hAxis = '{hAxis}'")
    print(f" inputFilePath = '{inputFilePath}'")
    inputFile_regex = rf"entanglementNtuple_{sampleName}_{hAxis}Axis_[0-9]+.root"
    inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
    numInputFiles = len(inputFileNames)
    print("Found %i input files." % numInputFiles)
    for mode in modes:
      for decayMode in decayModes:
        if mode != "gen_smeared":
          for spinAnalyzer in spinAnalyzers:
            cfgFileName_analysis_modified = os.path.join(
              configDir, f"analyzeEntanglementNtuple_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode_{spinAnalyzer}_cfg.py"
            )
            outputFileName_analysis = f"analyzeEntanglementNtuple_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode_{spinAnalyzer}.root"
            args_analysis = {
              'inputFileNames'      : inputFileNames,
              'par_gen'             : par_gen,
              'mode'                : mode,
              'collider'            : collider,
              'decayMode'           : decayMode,
              'apply_evtWeight'     : sample['apply_evtWeight'],
              'spinAnalyzer'        : spinAnalyzer,
              'maxSumPhotonEn'      : maxSumPhotonEn,
              'bootstrapSize'       : bootstrap_size,
              'numBootstrapSamples' : bootstrap_count,
              'outputFileName'      : outputFileName_analysis,
            }
            build_cfg(analyzeNtuple_template, cfgFileName_analysis_modified, args_analysis)

            logFileName_analysis = cfgFileName_analysis_modified.replace("_cfg.py", ".log")
            job_key_analysis = f'{sampleName}_{mode}_{hAxis}_{decayMode}_{spinAnalyzer}_analysis'
            dependencies_analysis = [ cfgFileName_analysis_modified ]
            dependencies_analysis.extend(inputFileNames)
            jobOptions_analysis[job_key_analysis] = {
              'inputFileNames' : dependencies_analysis,
              'cfgFileName'    : cfgFileName_analysis_modified,
              'outputFilePath' : outputDir,
              'outputFileName' : outputFileName_analysis,
              'logFileName'    : logFileName_analysis,
              'cmd'            : analyzeEntanglementNtuple_cmd,
            }
        cfgFileName_ctrlPlots_modified = os.path.join(
          configDir, f"makeControlPlots_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode_cfg.py"
        )
        outputFileName_ctrlPlots = f"makeControlPlots_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode.root"
        if run_makefile:
          outputFileName_ctrlPlots = os.path.join(outputDir_ctrlPlots, outputFileName_ctrlPlots)
        args_ctrlPlots = {
          'inputFileNames'  : inputFileNames,
          'mode'            : mode,
          'collider'        : collider,
          'decayMode'       : decayMode,
          'apply_evtWeight' : sample['apply_evtWeight'],
          'maxSumPhotonEn'  : maxSumPhotonEn,
          'outputFileName'  : outputFileName_ctrlPlots,
        }
        build_cfg(makeControlPlot_template, cfgFileName_ctrlPlots_modified, args_ctrlPlots)

        logFileName_ctrlPlots = cfgFileName_ctrlPlots_modified.replace("_cfg.py", ".log")
        job_key_ctrlPlots = f'{sampleName}_{mode}_{hAxis}_{decayMode}_ctrlPlots'
        dependencies_ctrlPlots = [ cfgFileName_ctrlPlots_modified ]
        dependencies_ctrlPlots.extend(inputFileNames) 
        jobOptions_ctrlPlots[job_key_ctrlPlots] = {
          'inputFileNames' : dependencies_ctrlPlots,
          'cfgFileName'    : cfgFileName_ctrlPlots_modified,
          'outputFilePath' : outputDir_ctrlPlots,
          'outputFileName' : outputFileName_ctrlPlots,
          'logFileName'    : logFileName_ctrlPlots,
          'cmd'            : makeControlPlots_cmd,
        }
        if mode != "gen":
          cfgFileName_resPlots_modified = os.path.join(
            configDir, f"makeResolutionPlots_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode_cfg.py"
          )

          outputFileName_resPlots = f"makeResolutionPlots_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode.root"
          if run_makefile:
            outputFileName_resPlots = os.path.join(outputDir_resPlots, outputFileName_resPlots)
          args_resPlots = {
            'inputFileNames'  : inputFileNames,
            'mode'            : mode,
            'collider'        : collider,
            'decayMode'       : decayMode,
            'apply_evtWeight' : sample['apply_evtWeight'],
            'maxSumPhotonEn'  : maxSumPhotonEn,
            'outputFileName'  : outputFileName_resPlots,
          }
          build_cfg(makeResolutionPlot_template, cfgFileName_resPlots_modified, args_resPlots)

          logFileName_resPlots = cfgFileName_resPlots_modified.replace("_cfg.py", ".log")
          job_key_resPlots = f'{sampleName}_{mode}_{hAxis}_{decayMode}_resPlots'
          dependencies_resPlots = [ cfgFileName_resPlots_modified ]
          dependencies_resPlots.extend(inputFileNames)
          jobOptions_resPlots[job_key_resPlots] = {
            'inputFileNames' : dependencies_resPlots,
            'cfgFileName'    : cfgFileName_resPlots_modified,
            'outputFilePath' : outputDir,
            'outputFileName' : outputFileName_resPlots,
            'logFileName'    : logFileName_resPlots,
            'cmd'            : makeResolutionPlots_cmd,
          }

jobOptions_Makefile = []
for job_key, job in jobOptions_analysis.items():
  commands = []
  commands.append('rm -f {}'.format(job['outputFileName']))
  commands.append('/usr/bin/time --verbose {} {} &> {}'.format(analyzeEntanglementNtuple_cmd, job['cfgFileName'], job['logFileName']))
  commands.append('cp -v {} {}'.format(job['outputFileName'], os.path.join(outputDir, job['outputFileName'])))
  commands.append('rm -f {}'.format(job['outputFileName']))
  jobOptions_Makefile.append({
    'target'          : os.path.join(outputDir, job['outputFileName']),
    'dependencies'    : [ inputFileName.replace("file:", "") for inputFileName in job['inputFileNames'] ],
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, job['outputFileName']) ],
  })
for job_key, job in jobOptions_ctrlPlots.items():
  commands = []
  commands.append('rm -f {}'.format(job['outputFileName']))
  commands.append('/usr/bin/time --verbose {} {} &> {}'.format(makeControlPlots_cmd, job['cfgFileName'], job['logFileName']))
  jobOptions_Makefile.append({
    'target'          : os.path.join(outputDir, job['outputFileName']),
    'dependencies'    : [ inputFileName.replace("file:", "") for inputFileName in job['inputFileNames'] ],
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, job['outputFileName']) ],
  })
for job_key, job in jobOptions_resPlots.items():
  commands = []
  commands.append('rm -f {}'.format(job['outputFileName']))
  commands.append('/usr/bin/time --verbose {} {} &> {}'.format(makeResolutionPlots_cmd, job['cfgFileName'], job['logFileName']))
  jobOptions_Makefile.append({
    'target'          : os.path.join(outputDir, job['outputFileName']),
    'dependencies'    : [ inputFileName.replace("file:", "") for inputFileName in job['inputFileNames'] ],
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, job['outputFileName']) ],
  })

jobOptions = { **jobOptions_analysis, **jobOptions_ctrlPlots, **jobOptions_resPlots }
message = f"Finished building config files for {len(jobOptions)} job(s)."
if run_makefile:
  makeFileName = os.path.join(configDir, "Makefile")
  build_Makefile(makeFileName, jobOptions_Makefile)
  message += f" Now execute 'make -j 12 -f {makeFileName}' to start the jobs."
else:
  sbatchSubmissionFileName = os.path.join(configDir, "sbatch_submission.sh")
  build_sbatchSubmission(sbatchSubmissionFileName, jobOptions, 'produceEntanglementNtuple')
  os.chmod(sbatchSubmissionFileName, 0o755)
  message += f" Now execute '{sbatchSubmissionFileName}' to submit the jobs to SLURM."

print(message)
